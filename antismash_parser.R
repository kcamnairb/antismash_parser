#!/usr/bin/env Rscript 

library(rvest, quietly=TRUE)
library(janitor, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(argparser, quietly=TRUE)

parser = arg_parser('Parse Antismash results into a csv or optionally a bed file. \
                    An antismash results directory or url can be provided')
parser = add_argument(parser, 'antismash_results', help='antismash directory or url to parse')
parser = add_argument(parser, 'output_csv', help='output csv file')
parser = add_argument(parser, '--bed_file', help='output bed file')
args =  parse_args(parser)
antismash_results = args$antismash_results
results_type = if_else(str_ends(antismash_results, 'index.html'), 'html', 'directory')
output_file = args$output_csv
parse_webpage = function(antismash_results){
  if (results_type == 'directory'){
    html = read_html(paste0(antismash_results, '/index.html'))
  } else {
    html = read_html(antismash_results)
  }
  text =  html %>% html_text2()
  chroms = str_match_all(text, '\\n(?<chrom>.+?)\\nRegion Type From To Most')[[1]] %>%
    as_tibble(.name_repair='unique') %>% pull(chrom) %>% str_remove(' .*|\\(.*')
  #chroms = chroms[1:length(chroms)-1]
  tables = html %>% html_table()
  clusters = tables[1:length(chroms)]
  names(clusters) = chroms
  clusters = clusters %>%
    imap(~janitor::clean_names(.x) %>%
          rename(backbone_type = most_similar_known_cluster_2, start=from, end=to) %>%
          mutate(region = str_replace(region, 'Region&nbsp', 'region_')) %>%
          mutate(across(c(start, end), ~str_remove_all(.x, ',') %>% as.numeric())) %>%
          mutate(similarity = str_remove(similarity, '%') %>% as.numeric()) %>%
          mutate(chrom = .y)) %>%
    bind_rows() %>%
    distinct() %>%
    arrange(chrom, start)
  return(clusters)
}
get_genes_in_clusters = function(antismash_directory){
  cluster_blasts = list.files(paste0(antismash_directory, '/clusterblast'), full.names = T) %>%
    map(~read_lines(.x, skip=3) %>%
          tibble(text = .) %>%
          filter(cumsum(str_detect(text, 'Significant hits:')) < 1) %>%
          mutate(filename = .x)) %>%
    bind_rows()
  genes_in_clusters = cluster_blasts %>%
    filter(text != '') %>%
    mutate(gene_id = str_remove(text, '\t.*'),
           cluster_num = str_replace(filename, '.*/.*?_c(\\d+).txt', '\\1'),
           cluster_blast_key = str_remove(filename, '.*/')) %>%
    group_by(cluster_blast_key) %>%
    summarize(gene_id = str_c(gene_id, collapse=';'))
  return(genes_in_clusters)
}

webpage_table = parse_webpage(antismash_results)
## If a directory is provided, fetch the gene ids from the antismash clusterblast results directory
if (results_type == 'directory'){
  cluster_genes = get_genes_in_clusters(antismash_results)
  res = webpage_table %>%
    group_by(chrom) %>%
    arrange(start) %>%
    mutate(chrom_clust_num = row_number()) %>%
    ungroup() %>% 
    mutate(cluster_blast_key = str_glue('{chrom}_c{chrom_clust_num}.txt')) %>%
    left_join(cluster_genes, by='cluster_blast_key') %>% 
    select(-cluster_blast_key, -chrom_clust_num)
} else {
  res = webpage_table
}
res %>% write_csv(output_file, na='')
if (! is.na(args$bed_file)){
  res %>%
    select(chrom, start, end, region) %>%
    write_tsv(args$bed_file, col_names = FALSE)
}

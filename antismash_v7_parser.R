#!/usr/bin/env Rscript 
library(rvest)
library(janitor)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
antismash_directory = args[1]
parse_antismash_7 = function(antismash_directory){
  tables = read_html(paste0(antismash_directory, '/index.html')) %>% 
    html_table()
  clusters = tables %>% 
    keep(~'Region' %in% colnames(.x)) %>%
    map(janitor::clean_names) %>%
    bind_rows() %>%
    distinct() %>%
    rename(backbone_type = most_similar_known_cluster_2, start=from, end=to) %>%
    mutate(region = str_replace(region, 'Region&nbsp', 'region_')) %>%
    mutate(across(c(start, end), ~str_remove_all(.x, ',') %>% as.numeric())) %>%
    mutate(similarity = str_remove(similarity, '%') %>% as.numeric())
  cluster_blasts = list.files(paste0(antismash_directory, '/clusterblast'), full.names = T) %>%
    map(~read_lines(.x, skip=3) %>%
          tibble(text = .) %>%
          filter(cumsum(str_detect(text, 'Significant hits:')) < 1) %>%
          mutate(filename = .x)) %>%
    bind_rows()
  region_chrom_and_gene = cluster_blasts %>%
    filter(text != '') %>%
    mutate(gene_id = str_remove(text, '\t.*'),
           chrom = str_replace(filename, '.*/(.*?)_c\\d+.txt', '\\1'),
           cluster_num = str_replace(filename, '.*/.*?_c(\\d+).txt', '\\1')) %>%
    select(-text, -filename) %>%
    group_by(chrom) %>%
    mutate(chrom_num = if_else(row_number() == 1, 1, 0)) %>%
    ungroup() %>%
    mutate(chrom_num = cumsum(chrom_num),
           region = str_glue('region_{chrom_num}.{cluster_num}')) %>%
    group_by(region) %>%
    summarize(gene_id = str_c(gene_id, collapse=';'), chrom=first(chrom))
  clusters %>% left_join(region_chrom_and_gene, by='region')
}
antismash_results = parse_antismash_7(antismash_directory)
antismash_results %>% write_csv(paste0(basename(antismash_directory), '_BGC_summary.csv', na=''))
#!/usr/bin/env python3
import argparse
import json
import os
import re
import pandas as pd
from glob import glob
parser = argparse.ArgumentParser(description="""Extracts regions from antismash version 6 output and saves to a csv file, 
and to 4 column bed file. Csv output will be similar to this:


+------------+-------------+------------------+---------+---------+-------------+----------------------------+------------+-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|   contig   | contig_edge |     product      |  start  |   end   | region_name | Most similar known cluster | Similarity |     known_cluster_type      |                                                                                                                                 genes_in_cluster                                                                                                                                 |
+------------+-------------+------------------+---------+---------+-------------+----------------------------+------------+-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| EQ963472.1 | FALSE       | NRPS-like        | 1875147 | 1917746 |         1.1 |                            |            |                             | A_flavus_3357_000682;A_flavus_3357_000683;A_flavus_3357_000684;A_flavus_3357_000685;A_flavus_3357_000686;A_flavus_3357_000687;A_flavus_3357_000688                                                                                                                               |
| EQ963472.1 | FALSE       | T1PKS            | 2207321 | 2247495 |         1.2 |                            |            |                             | A_flavus_3357_000787;A_flavus_3357_000788;A_flavus_3357_000789;A_flavus_3357_000790;A_flavus_3357_000791;A_flavus_3357_000792;A_flavus_3357_000793;A_flavus_3357_000794;A_flavus_3357_000795;A_flavus_3357_000796;A_flavus_3357_000797;A_flavus_3357_000798;A_flavus_3357_000799 |
| EQ963472.1 | FALSE       | NRPS-like, T1PKS | 2594700 | 2672837 |         1.3 | asparasone A               | 75%        | Polyketide:Iterative type I | A_flavus_3357_000936;A_flavus_3357_000937;A_flavus_3357_000938;A_flavus_3357_000939;A_flavus_3357_000940;A_flavus_3357_000941;                                                                                                                                                   |
| EQ963472.1 | FALSE       | NRPS             | 2708378 | 2756022 |         1.4 |                            |            |                             | A_flavus_3357_000978;A_flavus_3357_000979;A_flavus_3357_000980;A_flavus_3357_000981;A_flavus_3357_000982;A_flavus_3357_00098357_000989;A_flavus_3357_000990;A_flavus_3357_000991;A_flavus_3357_000992;A_flavus_3357_000993                                                       |
| EQ963472.1 | FALSE       | indole           | 2926771 | 2948131 |         1.5 |                            |            |                             | A_flavus_3357_001057;A_flavus_3357_001058;A_flavus_3357_001059;A_flavus_3357_001060;A_flavus_3357_001061                                                                                                                                                                         |
+------------+-------------+------------------+---------+---------+-------------+----------------------------+------------+-----------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

""", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('antismash_v6_directory')
args = parser.parse_args()
json_file = glob(args.antismash_v6_directory + '/*.json')[0]
html_index = os.path.join(args.antismash_v6_directory, 'index.html')
knownclusterblasts = glob(os.path.join(args.antismash_v6_directory, 'knownclusterblast') + '/*.txt')
as5 = json.load(open(json_file))
regions = []
for contig in as5['records']:
    for f in contig['features']:
        if 'region_number' in f['qualifiers']:
            regions.append(dict(contig=contig['name'], location=f['location'], **f['qualifiers']))
regions_df = pd.DataFrame(regions)
contig_to_num = {record['name']:num for num, record in enumerate(as5['records'], 1)}
regions_df['product'] = regions_df['product'].apply(', '.join)
regions_df['contig_edge'] = regions_df['contig_edge'].apply(', '.join)
regions_df['region_number'] = regions_df['region_number'].apply(', '.join)
regions_df[['start','end']] = regions_df['location'].str.replace('[\[\]]', '').str.split(':', expand=True)
regions_df['start'] = regions_df['start'].str.replace('<|>','')
regions_df['end'] = regions_df['end'].str.replace('<|>','')
def contig_and_region_num_to_name(row):
    contig_num = contig_to_num[row['contig']]
    return str(contig_num) + '.' + row['region_number']
regions_df['region_name'] = regions_df.apply(contig_and_region_num_to_name, axis=1)
regions_df = regions_df.drop(['candidate_cluster_numbers', 'subregion_numbers', 'tool', 'rules', 'region_number', 'location'], axis=1)
base = os.path.basename(args.antismash_v6_directory).replace('json', '')
def get_known_clusterblast_results(html_index):
    dfs = []
    for df in pd.read_html(html_index)[:-1]:
        if df.shape[1] == 7:
            df.columns = ['region_name','Type','From','To','Most similar known cluster','known_cluster_type','Similarity']
            dfs.append(df)
        else:
            dfs.append(df)
    df = pd.concat(dfs)
    df.region_name = df.region_name.str.replace('Region&nbsp','')
    df = df[['region_name','Most similar known cluster','Similarity','known_cluster_type']]         
    return df
known_clusters = get_known_clusterblast_results(html_index)
def get_genes_in_cluster(clusterblast):
    region_name = re.search('(.*)_c(\d+?).txt', os.path.basename(clusterblast))
    contig, region_number = region_name.group(1), region_name.group(2)
    region_name = str(contig_to_num[contig]) + '.' + region_number
    genes_in_cluster = []
    for idx, line in enumerate(open(clusterblast)):
        if idx > 2:
            if line.startswith('\n'): 
                break
            genes_in_cluster.append(line.split('\t')[0])
    return (region_name, ';'.join(genes_in_cluster))
genes_in_cluster = pd.DataFrame.from_records([get_genes_in_cluster(clusterblast) for clusterblast in knownclusterblasts], columns=['region_name', 'genes_in_cluster']).drop_duplicates()
regions_df = regions_df.merge(known_clusters.drop_duplicates(), on='region_name', how='left')
regions_df = regions_df.merge(genes_in_cluster.drop_duplicates(), on='region_name', how='left')
regions_df = regions_df
regions_df.to_csv(os.path.join(args.antismash_v6_directory,base+'.csv'), index=False)
regions_df[['contig','start', 'end', 'region_name']].to_csv(os.path.join(args.antismash_v6_directory,base+'.bed'), index=False, header=False, sep='\t')

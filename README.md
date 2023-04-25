# Antismash Results Parser

This is an R script that parses an antismash results directory or URL into a summary csv file. Antismash is a tool for identifying biosynthetic gene clusters in bacterial and fungal genomes. The summary csv file contains information on the cluster number, type, coordinates, the genes in the cluster, and the most similar known cluster for each gene cluster. If a directory is give as input, the genes in each cluster will also be reported in the summary file. An optional bed file can also be created. This is able to parse results from antismash v5 - v7.

## Installation

To use this script, you need to have R installed on your system. You also need to install the following R packages: tidyverse, rvest, janitor, and argparser. You can install them using the command `install.packages(c("tidyverse", "rvest", "janitor", "argparser"))`.

## Usage

To run the script, you need to provide two arguments: the path to the antismash results directory or URL and the name of the output csv file.

For example:

`Rscript antismash_parser.R /path/to/antismash/results output.csv`

Include a bed file as output:

`Rscript antismash_parser.R --bed_file output.bed /path/to/antismash/results output.csv`

It can also parse directly from the antismash URL:

`Rscript antismash_parser.R http://antismash_results/index.html output.csv`

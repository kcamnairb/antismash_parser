# Antismash Results Parser

This is an R script that parses an antismash results directory into a summary csv file. Antismash is a tool for identifying biosynthetic gene clusters in bacterial and fungal genomes. The summary csv file contains information on the cluster number, type, coordinates, the genes in the cluster, and the most similar known cluster for each gene cluster. This is able to parse results from antismash v5 - v7.

## Installation

To use this script, you need to have R installed on your system. You also need to install the following R packages: tidyverse, rvest, and janitor. You can install them using the command `install.packages(c("tidyverse", "rvest", "janitor"))`.

## Usage

To run the script, you need to provide two arguments: the path to the antismash results directory and the name of the output csv file.

For example:

`Rscript antismash_parser.R /path/to/antismash/results output.csv`

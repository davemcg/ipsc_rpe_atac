#!/usr/bin/env Rscript
# Read in fimo motif file, print out counts of each motif
args = commandArgs(trailingOnly=TRUE)

file_name = args[1]

library(tidyverse)

input = read_tsv(file_name)

input <- input %>% mutate(sequence_name = as.character(sequence_name)) %>% mutate(sample = sample, bootstrap = bootstrap_num)
# collapse to just counts per motif

output <- input %>% group_by(sample, motif_alt_id) %>% summarise(Count=n()) %>% ungroup()

write_tsv(output, path = args[2])

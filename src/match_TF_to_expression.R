#Rscript

library(rvest)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
homer_html <- args[1]
gene_col <- args[2]
output_file <- args[3]

lsTPM <- read_tsv('/home/mcgaugheyd/git/ipsc_rpe_RNA-seq/data/lsTPM_by_Line.tsv')
delta <- read_csv('/home/mcgaugheyd/git/ipsc_rpe_RNA-seq/data/RPE_vs_iPSC.results.csv')
homer <- read_html(homer_html)
processed <- (homer %>% html_table(header = T))[[1]]
processed$Gene = sapply(processed[,gene_col], function(x) str_split(x, '\\(|\\/', n = 2)[[1]][1])
processed$HomerName = processed[,gene_col]
processed <- as.tibble(processed) 

# left join's homer output with ipsc-RNA RNA-seq data for lsTPM values by GFP/RFP/iPSC
# and RNA-seq differential expression of RPE(GFP+RFP) vs iPSC
merged <- processed %>% 
			mutate(Gene = toupper(Gene)) %>% 
			left_join(., lsTPM %>% dplyr::select(Gene, Line, `log2(lsTPM)`) %>%  
			spread(Line, `log2(lsTPM)`), by='Gene') %>% 
			left_join(., delta, by = 'Gene')

write_tsv(merged, path = args[3])


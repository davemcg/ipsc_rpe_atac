#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

###########################
# Integrates common ATAC peaks with closest TSS and homer ID'ed TFBS motifs
###########################

library(tidyverse)
library(data.table)
library(annotables)

common_peaks_file <- args[1] # common_peaks_file <- '/Volumes/data/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/macs_peak/all_common_peaks.blackListed.narrowPeak.closestTSS.bed'
motif_bed_file <- args[2] # motif_bed_file <- '/Volumes/data/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/homer/motif5.motif.bed'
gtf_file <- args[3]
closest_num <- args[4] #3 closest?
output_file <- args[5] 

common_peaks_raw <- fread(cmd = paste('gzcat', common_peaks_file))
# 1. filter to closest n (args[4]) genes within 500,000bp
# 2. add gene name and description
# 3. keep most useful columns
common_peaks <- common_peaks_raw %>% 
  # part 1
  filter(V22 < 500000) %>% 
  group_by(V1, V2, V3) %>% 
  arrange(V1, V2, V3, V22) %>% 
  top_n(closest_num) %>% 
  # part 2
  mutate(enstxp = gsub('\\..*','', V13)) %>% 
  left_join(grch37_tx2gene) %>% 
  left_join(grch37) %>% 
  # part 3
  select(chrom = V1, start = V2, end = V3, id = V4, peak_score = V5, strand = V6, thickStart = V7, thickEnd = V8, rgb = V9, distance = V22, symbol, description) %>% 
  unique() %>% 
  mutate(id = case_when(rgb == '255,0,0' ~ 'RFP',
                        rgb == '0,255,0' ~ 'GFP',
                        TRUE ~ 'iPSC'))
  
# load in motif bed
motif_bed <- fread(motif_bed_file, skip = 1)
motif <- motif_bed$V4[1]
motif_bed[,motif] <- motif_bed$V5
motif_bed <- motif_bed %>% as.tibble() %>% select(chrom = V1, start = V2, end = V3, 7)

# merge!
output <- left_join(common_peaks, motif_bed)
#!/usr/bin/env Rscript
# reads in the closest TSS file and parses and creates
# a slimmed down Rdata data frame
args = commandArgs(trailingOnly=TRUE)

file_name = args[1]
gtf_path = args[2]
output_name = args[3]

library(tidyverse)

# load gtf to get gene name 
gtf <- read_tsv(gtf_path, col_names = c('Transcript','Gene'))
print(head(gtf))
# load in nearest tss motif info
sample <- str_split(str_split(file_name, '/')[[1]][2], '\\.')[[1]][1]
input <- read_tsv(file_name, col_names = c('sequence_name', 'start', 'end', 'motif', 'fimo_pvalue', 'strand', 'tss_seq','tss_start','tss_end', 'transcript',  'blank','tss_strand','coord1','coord2','blank2','exon_num','size','exon_pos','distance'))
input <- input %>% mutate(sequence_name = as.character(sequence_name),
                          tss_seq = as.character(tss_seq),
                          coord1=as.numeric(coord1),
                          coord2=as.numeric(coord2),
                          blank2 = as.character(blank2),
                          exon_num=as.integer(exon_num),
                          size=as.numeric(size)) %>% 
  mutate(sample = as.factor(sample)) %>% 
  rowwise() %>% 
  mutate(Transcript = str_split(transcript, '_')[[1]][1])
print(head(input))
both <- left_join(input, gtf) %>% 
  filter(!is.na(Gene)) %>% 
  mutate(motif_loc = paste(sequence_name, start, end, sep='_'))
print(head(both))
out <- both %>%
  # remove any TSS over 500kb away
  filter(distance < 500000) %>%
  #  one gene per motif
  group_by(motif_loc, Gene) %>%
  top_n(1, distance) %>%
  ungroup() %>%
  # add up to two genes (total) per motif
  group_by(motif_loc) %>%
  top_n(2, distance) %>%
  ungroup() %>%
  # arrange by genes most linked to motif
  group_by(Gene) %>% 
  select(sample, sequence_name, start, end, motif, fimo_pvalue, strand, Transcript, Gene, distance, motif_loc)

# %>%
#   summarise(Count=n(), Motif_Regions = paste(motif_loc, collapse=', ')) %>%
#   arrange(-Count)
print(head(out))

print(dim(out))
write_tsv(out, 
          path = output_name)



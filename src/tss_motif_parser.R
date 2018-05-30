#!/usr/bin/env Rscript
# reads in the closest TSS file and parses and creates
# a slimmed down Rdata data frame
args = commandArgs(trailingOnly=TRUE)

file_name = args[1]
gtf = args[2]
output_name = args[3]

library(tidyverse)

# load gtf to get gene name 
gtf <- read_tsv(gtf, col_names = c('Transcript','Gene')) #"/data/mcgaugheyd/genomes/1000G_phase2_GRCh37/gencode.v28lift37.metadata.HGNC.gz"

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

both <- left_join(input, gtf) %>% 
  filter(!is.na(Gene)) %>% 
  mutate(motif_loc = paste(sequence_name, start, end, sep='_'))

write_tsv(both %>% 
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
            summarise(Count=n(), paste(motif_loc, collapse=', ')) %>% 
            arrange(-Count), 
          path = output_name)



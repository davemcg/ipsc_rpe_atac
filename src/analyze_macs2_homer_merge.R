library(tidyverse)
peak <- read_tsv('/Volumes/data/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/peak_full/all_common_peaks.blackListed.narrowPeak.closestTSS__interesting_homer_motif.bed.gz', col_names = F)
expression <- read_csv('~/git/ipsc_rpe_RNA-seq/data/lsTPM_by_Line_with_Diff_Exp.tsv')

colnames(peak) <- c('chrom','start','stop','type','macs2_score','strand','thickStart','thickStop','RGB','Distance_to_Gene','Gene','Description','chrom_motif','start_motif','stop_motif','motif','motif_score','strand_motif')

motifs <- peak$motif %>% unique()

# create gene centric info on peak/motif and differential gene expression iPSC <-> RPE
grouped_stats <- list()

for (i in motifs){
  out <- peak %>% 
    filter(!is.na(Gene)) %>% # remove NA genes where the peak didn't get a match (no genes within 500kb of a peak)
    filter(motif == i) %>% # only retain the motif TFBS homer hits
    group_by(type, chrom, start, stop, Gene) %>% # group by GFP/RFP,iPSC, the macs2 peak, and the gene
    summarise(Num_TFBS_motifs=n()) %>%  # count total number of motifs
    group_by(type, Gene) %>% # group now by GFP/RFP/IPSC and gene
    summarise(Total_TFBS_motifs = sum(Num_TFBS_motifs), Num_Peaks = n()) %>% # now add num of peaks near gene
    select(type, Gene, Total_TFBS_motifs, Num_Peaks)
  grouped_stats[[i]] <- out
  print(i)
}

# further cut down to most interesting genes by differential peak activity and 
# expression by TF

TF_gene_stats_list <- list()

for (i in motifs){
  
  # ID differentially peak activity around genes for a TF
  out <- grouped_stats[[i]] %>% select(-Total_TFBS_motifs) %>% 
    spread(type, Num_Peaks) %>% 
    mutate(ratio = GFP/iPSC, delta = GFP - iPSC) %>% 
    arrange(-delta) %>% 
    left_join(expression %>% 
                select(Gene, log2FoldChange__RPE_vs_iPSC, padj__RPE_vs_iPSC) %>% 
                unique()) %>% 
    filter(!is.na(padj__RPE_vs_iPSC), 
           abs(log2FoldChange__RPE_vs_iPSC) > 1, 
           padj__RPE_vs_iPSC < 0.01) %>% 
    filter(delta > 4, ratio > 1.5) # only keep 
  TF_gene_stats_list[[i]] <- out
  print(i)
}

TF_gene_stats <- bind_rows(TF_gene_stats_list, .id = 'motif') %>% rowwise() %>% mutate(motif_core = toupper(gsub('\\(.*|11-GGGTGTGT,BestGuess:','',motif))) %>% ungroup()

# count num of genes associated with TF
TF_gene_stats %>% group_by(motif_core) %>% summarise(count = n(), Gene = paste(Gene, collapse = ',') ) %>% arrange(-count) %>% data.frame()

# most common genes associated with TF
TF_gene_stats %>% group_by(Gene) %>% summarise(count = n(), motif = paste(motif_core, collapse=',')) %>% arrange(-count) %>% data.frame()

# cut down for play
library(visNetwork)
TF_gene_stats_n <- TF_gene_stats %>% 
  group_by(motif_core) %>% 
  top_n(10, delta) %>% # no more than n by tF
  filter(!Gene %in% c('FZD8','OPCML','PFKP','NRIP1','TSHZ3', 'SOX6')) %>% # trimming out most common genes between TF
  ungroup()

TF <- TF_gene_stats_n %>%
  distinct(motif_core) %>%
  rename(label = motif_core)

genes <- TF_gene_stats_n %>%
  distinct(Gene) %>%
  rename(label = Gene)

nodes <- full_join(TF, genes, by = "label") %>% rowid_to_column("id") %>% mutate(group = case_when(label %in% TF$label ~ 'TF', 
                                                                                                   TRUE ~ 'Gene'),
                                                                                 shape = 'ellipse') %>% 
  left_join(., expression %>% 
              select(Gene, value = log2FoldChange__RPE_vs_iPSC) %>% 
              unique(), by = c("label" = "Gene"))
  


edges <- TF_gene_stats_n %>%  
  select(TF = motif_core, Gene, weight = log2FoldChange__RPE_vs_iPSC) %>% 
  unique() %>% 
  left_join(nodes, by = c("TF" = "label")) %>% 
  rename(from = id) %>% 
  left_join(nodes, by = c("Gene" = "label")) %>% 
  rename(to = id) %>% 
  ungroup() %>% 
  select(from, to, weight)

visNetwork(nodes, edges) %>% 
  visGroups(groupname = 'TF', color='red') 

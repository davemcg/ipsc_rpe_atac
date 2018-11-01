peak <- read_tsv('/Volumes/data/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/peak_full/all_common_peaks.blackListed.narrowPeak.closestTSS__interesting_homer_motif.bed.gz', col_names = F)
expression <- read_csv('~/git/ipsc_rpe_RNA-seq/data/lsTPM_by_Line_with_Diff_Exp.tsv')

otx2 <- peak %>% 
  filter(!is.na(X11)) %>% # remove NA genes where the peak didn't get a match (no genes within 500kb of a peak)
  filter(grepl('OTX2', X16, ignore.case = T)) %>% # only retain OTX2 TFBS homer hits
  group_by(X4, X1, X2, X3, X11) %>% # group by GFP/RFP,iPSC, the macs2 peak, and the gene
  summarise(Num_OTX2_motifs=n()) %>%  # count total number of OTX2 motifs
  group_by(X4, X11) %>% # group now by GFP/RFP/IPSC and gene
  summarise(Total_OTX2_motifs = sum(Num_OTX2_motifs), Num_Peaks = n()) %>% # now add num of peaks near gene
  select(Tissue = X4, Gene = X11, Total_OTX2_motifs, Num_Peaks)



mitf <- peak %>% 
  filter(!is.na(X11)) %>% # remove NA genes where the peak didn't get a match (no genes within 500kb of a peak)
  filter(grepl('MITF', X16, ignore.case = T)) %>% # only retain OTX2 TFBS homer hits
  group_by(X4, X1, X2, X3, X11) %>% # group by GFP/RFP,iPSC, the macs2 peak, and the gene
  summarise(Num_TFBS_motifs=n()) %>%  # count total number of OTX2 motifs
  group_by(X4, X11) %>% # group now by GFP/RFP/IPSC and gene
  summarise(Total_TFBS_motifs = sum(Num_TFBS_motifs), Num_Peaks = n()) %>% # now add num of peaks near gene
  select(Tissue = X4, Gene = X11, Total_TFBS_motifs, Num_Peaks)

mitf %>% spread(Tissue, Num_Peaks) %>% arrange(-GFP)

mitf %>% select(-Total_TFBS_motifs) %>% spread(Tissue, Num_Peaks) %>% mutate(GFP_iPSC = GFP/iPSC) %>% arrange(-GFP_iPSC)

mitf %>% select(-Total_TFBS_motifs) %>% spread(Tissue, Num_Peaks) %>% mutate(GFP_iPSC = GFP/iPSC) %>% arrange(-GFP_iPSC) %>% left_join(expression %>% select(-Line, -lsTPM, -`log2(lsTPM)`, -Rank) %>% unique())

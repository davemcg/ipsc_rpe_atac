# Rscript to do the heavy analysis


library(tidyverse)
library(ggridges)
library(cowplot)
library(DESeq2)
library(ggrepel)

##########################
# Load in real data
###############################
samples_full = list.files('/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/fimo_motifs/', pattern = '*dat.gz', full.names = T)

#samples_full = list.files('/Volumes/McGaughey_S/fimo_motifs/', pattern = '*dat.gz', full.names = T)

fimo_data <- list()
for (i in samples_full){
  sample <- str_split(i, '/')[[1]] %>% tail(1) %>% str_split(., '\\.') %>% unlist() %>% head(1)
  input <- read_tsv(i)
  input <- input %>% mutate(sequence_name = as.character(sequence_name)) %>% mutate(sample = sample)
  fimo_data[[sample]] <- input
}
sample_motifs = bind_rows(fimo_data) %>% mutate(sample=as.factor(sample))
rm(fimo_data)

#tf_motif <- sample_motifs %>% select('TF' = `# motif_id`, motif_alt_id) %>% unique()
#save(tf_motif, '~/git/ipsc_rpe_atac/data/tf_motif.Rdata')
load('~/git/ipsc_rpe_atac/data/tf_motif.Rdata')

#########################
# Load in bootstraps
#########################
bootstrap_stats = list.files('/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/fimo_motifs/bootstrapping_stats/', pattern = '*dat.gz', full.names = T)
#bootstrap_stats = list.files('/Volumes/McGaughey_S/fimo_motifs/bootstrapping_stats/', pattern = '*dat.gz', full.names = T)

# big difference here vs above
# not saving full file, as each bootstrap is ~200mb COMPRESSED
# just saving the stats
# right now, counts by motif/sample
fimo_data <- list()
for (i in bootstrap_stats){
  sample <-  str_split(i, '/')[[1]] %>% tail(1) %>% str_split(., '\\.') %>% unlist() %>% head(1)
  bootstrap_num <-  str_split(i, '/')[[1]] %>% 
    tail(1) %>% 
    str_split(., '\\.') %>%
    unlist() %>% 
    head(2) %>% tail(1) %>% 
    str_split(., '_') %>% 
    unlist() %>% tail(1)
  sample_boot <- paste0(sample,'_',bootstrap_num)
  stats <- read_tsv(i) %>% mutate(sample = as.factor(sample), bootstrap = bootstrap_num)
  fimo_data[[sample_boot]] <- stats
}

bootstrap_counts = bind_rows(fimo_data) %>% mutate(sample=as.factor(sample)) %>% left_join(tf_motif)
rm(fimo_data)

# merge sample counts with bootstrap counts
all <- bind_rows(sample_motifs %>% group_by(sample, motif_alt_id) %>% summarise(Count=n()) %>% ungroup() %>% mutate(bootstrap = 'real'),
                 bootstrap_counts)

#############################
# DESeq2-based analysis
############################
cts <- all %>% filter(bootstrap=='real') %>% spread(sample, Count) %>% select(-bootstrap, -TF) %>% data.frame()
row.names(cts) <- cts$motif_alt_id
cts %>% head()
cts <- cts[,-1]

coldata <- all %>% select(sample) %>% unique() %>% filter(sample!='RFP_IIE_1') %>% mutate(Type =  case_when(grepl('GFP', sample) ~ 'GFP',
                                                                                                            grepl('RFP', sample) ~ 'RFP',
                                                                                                            TRUE ~ 'iPSC'))

dds <- DESeqDataSetFromMatrix(countData = cts %>% select(one_of(as.character(coldata$sample)) ),
                              colData = coldata,
                              design = ~ Type)
dds$Type <- relevel(dds$Type, ref='RFP')
dds <- DESeq(dds)

res <- results(DESeq(dds))
resLFC <- lfcShrink(dds, coef='Type_GFP_vs_RFP', type='apeglm')


###########################################
# Load in RNA-seq data
##############################################
gfp_rfp <- read_csv('~/git/ipsc_rpe_RNA-seq/data/GFP_vs_RFP.results.csv')
rpe_ipsc <- read_csv('~/git/ipsc_rpe_RNA-seq/data/RPE_vs_iPSC.results.csv')
#gfp_rfp %>% head()


##########################
# Gene <-> motif matching
# load data
###########################
samples_full = list.files('/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/closest_TSS_motifs/processed/',
                          pattern = '*rdata', full.names = T)

closestTSS_data <- list()
for (i in samples_full){
  sample <- str_split(i, '/')[[1]] %>% tail(1) %>% str_split(., '\\.') %>% unlist() %>% head(1)
  #input <- read_tsv(i, col_names = c('sample','sequence_name','start','end','motif','fimo_pvalue','strand','Transcript','Gene','distance','motif_loc'))
  load(i)
  out_some <- out #%>% filter(motif == 'M5700_1.02')
  closestTSS_data[[sample]] <- out_some
}
rm(out)
rm(out_some)

sample_closestTSS = bind_rows(closestTSS_data) %>% mutate(sample=as.factor(sample))


#############################
# Save data
#############################
save(sample_closestTSS, file='/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/sample_closestTSS.Rdata')
save(resLFC, file='/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/resLFC.Rdata')
sample_bootstrap_counts <- all
save(sample_bootstrap_counts, file='/data/mcgaugheyd/projects/nei/hufnagel/iPSC_RPE_ATAC_Seq/sample_bootstrap_counts.Rdata')
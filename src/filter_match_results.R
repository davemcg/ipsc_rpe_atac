library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

known_file <- args[1]
novel_file <- args[2]

known <- read_tsv(known_file)
novel <- read_tsv(novel_file)

both <- bind_rows(known %>% select(HomerName, `q-value (Benjamini)`,  `P-value`, padj, log2FoldChange, GFP, RFP, iPSC))

#output <- both %>% filter(`P-value` < 1e-20, padj < 0.01, abs(log2FoldChange) > 1) %>% pull(HomerName)
output <- both %>% filter(GFP > 1, `q-value (Benjamini)` < 0.001) %>% pull(HomerName)

write(output, file = args[3])


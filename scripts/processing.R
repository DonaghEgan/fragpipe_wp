library(readr)
library(tidyverse)
library(stringr)
library(grid)
library(PhosR)
library(sva)
library(pheatmap)
library(limma)
library(plyr)

# read input file
raw <- read_delim("/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/inputs/combined_protein_wp.tsv", delim = "\t", escape_double = FALSE, 
                  col_names = TRUE, trim_ws = TRUE) %>%
  select("Protein ID", "Gene", contains("MaxLFQ"))

# delete rows with all NA
raw <- raw[rowSums(is.na(raw)) != ncol(raw), ]

# gene label meta data
gene_meta <-  raw %>% select("Protein ID", "Gene")
saveRDS(gene_meta, "/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/outputs/gene_meta.Rds")

# change rownames
raw$Gene <- NULL
raw <- raw %>% column_to_rownames("Protein ID")

# Setting up experimental design data frame
experimental_design <- data.frame(colnames(raw))
colnames(experimental_design) <- "label"
experimental_design$ID <- str_extract(experimental_design$label, "\\d+(?=_[a-zA-Z]+.+$)")

# formatting sample conditions using Martina's e-mail (21 Jan) 
conditions <- data.frame(c(1,11,21,31,2,12,22,32,3,13,23,33,4,14,24,34,5,15,25,35,6,16,26,36,7,17,27,37,8,18,28,38,9,19,29,39,10,20,30,40), 
                         c(rep("CNTR_0min",4),
                           rep("CNTR_5min",4),
                           rep("CNTR_20min",4),
                           rep("CNTR_24h_unstimulated",4),
                           rep("CNTR_24h_stimulated",4),
                           rep("PD1_0min",4),
                           rep("PD1_5min",4),
                           rep("PD1_20min",4),
                           rep("PD1_24h_unstimulated",4),
                           rep("PD1_24h_stimulated",4)),rep(1:4, 10))
colnames(conditions) <- c("ID","condition","batch")

# Join conditions with experimental design 
experimental_design <- plyr::join(experimental_design, conditions, by = "ID", match = "all")

# assigning each condition to a column 
experimental_design$type <- sub("_.*", "", experimental_design$condition) 
experimental_design$time <- sapply(strsplit(experimental_design$condition, "_"), "[", 2)
experimental_design$stim = ifelse(grepl("_stimulated|_5min|_20",experimental_design$condition),"stim","unstim")
rownames(experimental_design) <- experimental_design$label
experimental_design$label <- NULL

#setting 0 values to NA 
raw[raw == "0"] <- NA

# Log normalize
log2_norm <- log2(raw)

# remove proteins not found in 80% of samples per condition
protein_filter <- selectGrps(as.matrix(log2_norm), experimental_design$condition, 0.8, n=1)

## median centering and imputation ####
impute_protein <- PhosR::medianScaling(protein_filter, scale = FALSE, grps = experimental_design$condition, reorder = FALSE, assay = NULL)
impute_protein <- PhosR::scImpute(impute_protein, 0.5, experimental_design$condition)
impute_protein <- PhosR::tImpute(impute_protein)

## Dealing with batch effects ####
modcombat <- model.matrix(~1, data = experimental_design)
batch <- experimental_design$batch
protein_normalized <- ComBat(dat = impute_protein, batch = batch,
                             mod = modcombat, par.prior = TRUE, prior.plots=F)

saveRDS(protein_normalized, "/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/outputs/protein_normalized.Rds")
saveRDS(experimental_design, "/mnt/8TB/Projects/POIAZ/Donagh/fragpipe_wp/outputs/experimental_design.Rds")

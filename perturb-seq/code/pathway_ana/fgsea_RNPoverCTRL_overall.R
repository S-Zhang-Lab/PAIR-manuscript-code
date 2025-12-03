rm(list = ls())
set.seed(7)
library(data.table) 
library(Seurat)
library(dplyr)
library(fgsea)
library(qs)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")
meta_data <- scRNA_seq@meta.data

# Convert loge to log2
scRNA_seq <- log2(exp(scRNA_seq[["RNA"]]$data)) %>% as.matrix()

all(colnames(scRNA_seq) == rownames(meta_data)) # TRUE

pwy <- c("c2")

RNP_exp <- scRNA_seq[, meta_data$treatment == "RNP"]
CTRL_exp <- scRNA_seq[, meta_data$treatment == "CTRL"]

# Calculate log2 fold changes
RNP_over_CTRL <- sort(rowMeans(RNP_exp) - rowMeans(CTRL_exp), decreasing = TRUE)

signature <- readRDS(paste0("./data/pathways/", pwy, ".rds"))

# Run fGSEA----
fgseaRes_RNP_over_CTRL <- fgsea(signature, RNP_over_CTRL, minSize=15, maxSize=500)
fgseaRes_RNP_over_CTRL <- as.data.frame(fgseaRes_RNP_over_CTRL)
rownames(fgseaRes_RNP_over_CTRL) <- fgseaRes_RNP_over_CTRL$pathway

# only keep significant pathways
fwrite(fgseaRes_RNP_over_CTRL, paste0("./res/2025_1203/RNP_over_CTRL_overall.csv"))

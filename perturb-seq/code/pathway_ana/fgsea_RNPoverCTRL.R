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

treatments <- c("CTRL", "RNP")
PAIR_tag <- unique(meta_data$assigned_tag)

# Load the pathway signatures
pwy_collections <- c("c2")
# pwy_collections <- c("h", "c2", "c6", "c7", "c8")

for (i in 1:length(PAIR_tag)) {
  celltype <- PAIR_tag[i]
  
  # Chosen cell types
  celltype_idx <- meta_data$assigned_tag == celltype
  meta_sub <- meta_data[celltype_idx, ]
  scRNA_seq_sub <- scRNA_seq[, celltype_idx]
  all(colnames(scRNA_seq_sub) == rownames(meta_sub)) # TRUE
  
  RNP_exp <- scRNA_seq_sub[, meta_sub$treatment == "RNP"]
  CTRL_exp <- scRNA_seq_sub[, meta_sub$treatment == "CTRL"]
  
  # Calculate log2 fold changes
  RNP_over_CTRL <- sort(rowMeans(RNP_exp) - rowMeans(CTRL_exp), decreasing = TRUE)

  for (j in 1:length(pwy_collections)) {
    pwy <- pwy_collections[j]
    signature <- readRDS(paste0("./data/pathways/", pwy, ".rds"))
    
    # Run fGSEA----
    fgseaRes_RNP_over_CTRL <- fgsea(signature, RNP_over_CTRL, minSize=15, maxSize=500)
    
    # only keep significant pathways
    fgseaRes_RNP_over_CTRL <- as.data.frame(fgseaRes_RNP_over_CTRL)
    rownames(fgseaRes_RNP_over_CTRL) <- fgseaRes_RNP_over_CTRL$pathway
    
    fwrite(fgseaRes_RNP_over_CTRL, paste0("./res/2025_1203/pwy_results_RNPoverCTRL/RNP_over_CTRL_", celltype, "_", pwy, ".csv"))

  }
}

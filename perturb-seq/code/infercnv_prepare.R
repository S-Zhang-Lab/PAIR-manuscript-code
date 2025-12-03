rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(Matrix)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")
meta_data <- scRNA_seq@meta.data

# seu_sub <- subset(
#   scRNA_seq,
#   subset = (
#     treatment == "RNP" |
#       (treatment == "CTRL" & assigned_tag == "NT_CRISPRa_NT_CasRx")
#   )
# )

seu_sub <- scRNA_seq

meta <- seu_sub@meta.data
meta$treat_PAIR <- paste0(meta$treatment, "_", meta$assigned_tag)

seu_sub@meta.data <- meta

expr_mat <- GetAssayData(seu_sub, slot = "data")
meta <- seu_sub@meta.data

writeMM(expr_mat, "./data/infercnv/expr_normalized.mtx")

write.csv(data.frame(gene = rownames(expr_mat)), "./data/infercnv/genes.csv", row.names = FALSE)
write.csv(data.frame(cell = colnames(expr_mat)), "./data/infercnv/cells.csv", row.names = FALSE)
write.csv(meta, "./data/infercnv/metadata.csv")

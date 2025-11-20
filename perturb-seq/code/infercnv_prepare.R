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
Idents(scRNA_seq) <- "treatment"

seu_RNP <- subset(scRNA_seq, subset = treatment == "RNP")

expr_mat <- GetAssayData(seu_RNP, slot = "data")
meta <- seu_RNP@meta.data

writeMM(expr_mat, "./data/infercnv/expr_normalized.mtx")

write.csv(data.frame(gene = rownames(expr_mat)), "./data/infercnv/genes.csv", row.names = FALSE)
write.csv(data.frame(cell = colnames(expr_mat)), "./data/infercnv/cells.csv", row.names = FALSE)
write.csv(meta, "./data/infercnv/metadata.csv")

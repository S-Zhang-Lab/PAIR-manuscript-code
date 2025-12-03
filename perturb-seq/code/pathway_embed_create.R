rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(tidyr)

pwy_ls <- readRDS("./data/pathways/c2.rds")
gene_embed <- readRDS("./data/embedding/hs_ProtTrans_embed_All.rds")

emb_mat <- as.matrix(gene_embed)

# Function to compute one pathway embedding (mean of its gene embeddings)
get_pathway_embedding <- function(genes) {
  genes_in_data <- intersect(genes, rownames(emb_mat))
  
  if (length(genes_in_data) == 0L) {
    # no genes found -> return NA vector
    return(rep(NA_real_, ncol(emb_mat)))
  }
  
  # subset and take column means
  colMeans(emb_mat[genes_in_data, , drop = FALSE])
}

# Apply to each pathway and bind into a matrix (pathways x 1024)
pathway_emb_mat <- t(vapply(
  pwy_ls,
  FUN = get_pathway_embedding,
  FUN.VALUE = numeric(ncol(emb_mat))
))


# Add pathway names as rownames if not already present
if (!is.null(names(pwy_ls))) {
  rownames(pathway_emb_mat) <- names(pwy_ls)
}

pathway_emb_mat_clean <- pathway_emb_mat[!apply(pathway_emb_mat, 1, anyNA), ]

qsave(pathway_emb_mat_clean, file = "./data/embedding/c2_pathway_embeddings.qs")

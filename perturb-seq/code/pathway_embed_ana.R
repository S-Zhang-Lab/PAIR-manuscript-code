rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(tidyr)
library(tibble)

# ---- 1. Load embedding matrix ----
pathway_emb <- qread("./data/embedding/c2_pathway_embeddings.qs")
# Assume: rows = pathways, cols = embedding dimensions
# We want rows = features (dims), cols = cells (pathways)
mat <- t(as.matrix(pathway_emb))

# Optional: give feature names if missing
if (is.null(rownames(mat))) {
  rownames(mat) <- paste0("dim_", seq_len(nrow(mat)))
}

# ---- 2. Create Seurat object using embeddings as data ----
seu <- CreateSeuratObject(
  counts = mat,
  min.cells = 0,
  min.features = 0
)

DefaultAssay(seu) <- "RNA"

# Store pathway names as metadata (one per "cell")
seu$pathway <- colnames(mat)

# IMPORTANT: in Seurat v5, use layers instead of @counts/@data/@scale.data
# Here we force all three layers to be exactly the embedding matrix
seu@assays[["RNA"]]@layers[["counts"]]     <- mat
seu@assays[["RNA"]]@layers[["data"]]       <- mat
seu@assays[["RNA"]]@layers[["scale.data"]] <- mat

seu <- RunPCA(
  seu,
  features = rownames(seu),
  layer = "counts",
  npcs = 30,
  verbose = FALSE
)

pcs_to_use <- 1:20
seu <- FindNeighbors(seu, dims = pcs_to_use)
seu <- FindClusters(seu, resolution = 1.0)

# ---- UMAP ----
seu <- RunUMAP(seu, dims = pcs_to_use, verbose = FALSE)

p <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("./res/2025_1203/pathway_embedding_umap_clusters.png", plot = p, width = 6, height = 5)

umap_df <- Embeddings(seu, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_id")   # here each "cell" = pathway

# Make sure the UMAP columns have nice names
colnames(umap_df)[colnames(umap_df) == "UMAP_1"] <- "UMAP_1"
colnames(umap_df)[colnames(umap_df) == "UMAP_2"] <- "UMAP_2"

# ---- 2. Extract metadata ----
meta_df <- seu@meta.data %>%
  rownames_to_column(var = "cell_id")

# ---- 3. Merge metadata + UMAP ----
merged_df <- left_join(meta_df, umap_df, by = "cell_id")

qsave(merged_df, "./data/pathway_embedding_metadata_umap.qs")

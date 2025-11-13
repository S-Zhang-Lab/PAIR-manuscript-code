rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(biomaRt)
library(Matrix)

seu <- qread("./data/raw_seu_with_hto_pair_tags.qs")
# ==== convert Ensembl to gene symbols ====
# Use Ensembl human gene database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract your Ensembl IDs
ens_ids <- rownames(seu)

# Query symbols
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ens_ids,
  mart = mart
)

# Keep unique and clean map
gene_map <- gene_map %>% distinct(ensembl_gene_id, .keep_all = TRUE)
rownames(gene_map) <- gene_map$ensembl_gene_id

# Replace rownames in Seurat object
new_symbols <- gene_map[ens_ids, "hgnc_symbol"]

# If symbol is missing, keep original ID
new_symbols[is.na(new_symbols) | new_symbols == ""] <- ens_ids[is.na(new_symbols) | new_symbols == ""]

# Assign
rownames(seu) <- new_symbols

# ==== merge duplicated gene names ====
mat1 <- GetAssayData(seu, layer = "counts.1")
mat2 <- GetAssayData(seu, layer = "counts.2")
mat3 <- GetAssayData(seu, layer = "counts.3")
mat_all <- cbind(mat1, mat2, mat3)

genes <- rownames(mat_all)

mat_merged <- rowsum(as.matrix(mat_all), group = genes)
mat_merged <- Matrix(mat_merged, sparse = TRUE)

seu_merged <- CreateSeuratObject(counts = mat_merged)
seu_merged@meta.data <- seu@meta.data[colnames(seu_merged), , drop = FALSE]

seu <- seu_merged
# ==== seu_merged preprocessing ====
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RPL|^RPS")

q <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
ggsave("./res/2025_1113/quality_metrics_before_qc.png", plot = q, width = 12, height = 4)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

seu <- ScaleData(seu, features = rownames(seu))

# PCA
seu <- RunPCA(seu, features = VariableFeatures(seu))
q <- ElbowPlot(seu)
ggsave("./res/2025_1113/ElbowPlot.png", plot = q, width = 6, height = 4)

dims_use <- 1:20

seu <- FindNeighbors(seu, dims = dims_use)
seu <- FindClusters(seu, resolution = 0.4)

seu <- RunUMAP(seu, dims = dims_use)

qsave(seu, "./data/seu_prep.qs")

# ==== umap plots ====
meta <- seu@meta.data
q <- DimPlot(seu, reduction = "umap", group.by = "treatment", pt.size = 0.5)
ggsave("./res/2025_1113/umap_treatment.png", plot = q, width = 5, height = 4)

q <- DimPlot(seu, reduction = "umap", group.by = "assigned_tag", pt.size = 0.5)
ggsave("./res/2025_1113/umap_assigned_tag.png", plot = q, width = 7, height = 4)

q <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5)
ggsave("./res/2025_1113/umap_seurat_clusters.png", plot = q, width = 5, height = 4)

# ==== umap plots for chosen PAIR====
tags_to_keep <- unique(meta$assigned_tag)  # replace with actual tags of interest
outdir <- "./res/2025_1113/"
for (tag in tags_to_keep) {
  
  q <- DimPlot(
    seu,
    group.by = "assigned_tag",
    cells.highlight = WhichCells(seu, expression = assigned_tag == tag),
    cols.highlight = "red",
    cols = "gray80"
  ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(tag)
  
  ggsave(
    filename = file.path(outdir, paste0("umap_", tag, ".png")),
    plot = q,
    width = 4,
    height = 4,
    dpi = 300
  )
}

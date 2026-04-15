rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Matrix)

# --- Source config.R (sets DATA_ROOT, SEU_RAW, SEU_PREP; see code/config.R) -
if (file.exists("code/config.R")) {
  source("code/config.R")                 # run from repo/
} else if (file.exists("config.R")) {
  source("config.R")                      # run from repo/code/
} else {
  stop("config.R not found. Run from repo/ or repo/code/.")
}
# ---------------------------------------------------------------------------

prep_fig_dir <- file.path(FIGURES_DIR, "legacy_prep")
prep_meta_dir <- file.path(OUTPUT_DIR, "legacy_prep")
reference_dir <- file.path(DATA_DIR, "reference")
gene_map_file <- file.path(reference_dir, "ensembl_to_hgnc_hsapiens.csv")

dir.create(prep_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prep_meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(SEU_PREP), recursive = TRUE, showWarnings = FALSE)

seu <- qread(SEU_RAW)
seu <- JoinLayers(seu)

# ==== convert Ensembl to gene symbols ====
ens_ids <- rownames(seu)
if (file.exists(gene_map_file)) {
  cat("Loading cached gene annotation:", gene_map_file, "\n")
  gene_map <- read.csv(gene_map_file, stringsAsFactors = FALSE)
} else {
  cat("Cached gene annotation not found. Building local cache...\n")
  local_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ens_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  gene_map <- data.frame(
    ensembl_gene_id = names(local_symbols),
    hgnc_symbol = unname(local_symbols),
    stringsAsFactors = FALSE
  )

  missing_local <- mean(is.na(gene_map$hgnc_symbol) | gene_map$hgnc_symbol == "")
  if (missing_local > 0.2) {
    cat("Local annotation missing", round(100 * missing_local, 1),
        "% of symbols; keeping Ensembl IDs for unmapped genes in offline mode.\n")
  }

  gene_map <- gene_map %>% distinct(ensembl_gene_id, .keep_all = TRUE)
  write.csv(gene_map, gene_map_file, row.names = FALSE)
}

gene_map <- gene_map %>% distinct(ensembl_gene_id, .keep_all = TRUE)
rownames(gene_map) <- gene_map$ensembl_gene_id

# Replace rownames in Seurat object
new_symbols <- gene_map[ens_ids, "hgnc_symbol"]

# If symbol is missing, keep original ID
new_symbols[is.na(new_symbols) | new_symbols == ""] <- ens_ids[is.na(new_symbols) | new_symbols == ""]

# Assign
rownames(seu) <- new_symbols

# ==== merge duplicated gene names after layer join ====
mat_all <- GetAssayData(seu, layer = "counts")
genes <- rownames(mat_all)

mat_merged <- rowsum(as.matrix(mat_all), group = genes)
mat_merged <- Matrix(mat_merged, sparse = TRUE)

meta_raw <- seu@meta.data
seu_merged <- CreateSeuratObject(counts = mat_merged, meta.data = meta_raw)

seu <- seu_merged
# ==== seu_merged preprocessing ====
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RPL|^RPS")

q <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
ggsave(file.path(prep_fig_dir, "quality_metrics_before_qc.png"), plot = q, width = 12, height = 4)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

seu <- ScaleData(seu, features = rownames(seu))

# PCA
seu <- RunPCA(seu, features = VariableFeatures(seu))
q <- ElbowPlot(seu)
ggsave(file.path(prep_fig_dir, "ElbowPlot.png"), plot = q, width = 6, height = 4)

dims_use <- 1:20

seu <- FindNeighbors(seu, dims = dims_use)
seu <- FindClusters(seu, resolution = 0.4)

seu <- RunUMAP(seu, dims = dims_use)

qsave(seu, SEU_PREP)

# ==== umap plots ====
meta <- seu@meta.data
q <- DimPlot(seu, reduction = "umap", group.by = "treatment", pt.size = 0.5)
ggsave(file.path(prep_fig_dir, "umap_treatment.png"), plot = q, width = 5, height = 4)

q <- DimPlot(seu, reduction = "umap", group.by = "assigned_tag", pt.size = 0.5)
ggsave(file.path(prep_fig_dir, "umap_assigned_tag.png"), plot = q, width = 7, height = 4)

q <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5)
ggsave(file.path(prep_fig_dir, "umap_seurat_clusters.png"), plot = q, width = 5, height = 4)

# ==== umap plots for chosen PAIR====
tags_to_keep <- sort(unique(na.omit(meta$assigned_tag)))
outdir <- prep_fig_dir
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

writeLines(capture.output(sessionInfo()), file.path(prep_meta_dir, "sessionInfo.txt"))

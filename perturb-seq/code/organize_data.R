rm(list = ls())
set.seed(7)

library(dplyr)
library(Matrix)
library(Seurat)
library(qs)

load_sample <- function(dir, sample_name) {
  
  # Read raw 10x files
  mat <- Matrix::readMM(file.path(dir, "matrix.mtx.gz"))
  barcodes <- read.table(file.path(dir, "barcodes.tsv.gz"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  features <- read.table(file.path(dir, "features.tsv.gz"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Assign row/column names
  rownames(mat) <- features$V1
  colnames(mat) <- barcodes$V1
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = mat)
  
  # Add metadata
  seu$sample <- sample_name
  seu$barcode_original <- colnames(seu)
  
  return(seu)
}

ctl1_path <- "./data/mRNA_HTO_Raw_Counts/Ctl1/"
rnp2_path <- "./data/mRNA_HTO_Raw_Counts/RNP2/"
rnp3_path <- "./data/mRNA_HTO_Raw_Counts/RNP3/"

ctl1_seu <- load_sample(ctl1_path, "Ctl1")
rnp2_seu <- load_sample(rnp2_path, "RNP2")
rnp3_seu <- load_sample(rnp3_path, "RNP3")

seu <- merge(ctl1_seu, y = c(rnp2_seu, rnp3_seu))
meta <- seu@meta.data
meta$treatment <- "RNP"
meta$treatment[meta$sample == "Ctl1"] <- "CTRL"

seu@meta.data <- meta

# ==== load PAIR results ====
# Specify file paths
barcodes_file <- "./data/PAIR_output/PAIR_matched_sparse_UMI_matrix_rows.txt"
features_file <- "./data/PAIR_output/PAIR_matched_sparse_UMI_matrix_columns.txt"
matrix_file   <- "./data/PAIR_output/PAIR_matched_sparse_UMI_matrix.mtx"

# Read the matrix (sparse format)
mat <- readMM(matrix_file)

# Read features and barcodes
features <- read.delim(features_file, header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim(barcodes_file, header = FALSE, stringsAsFactors = FALSE)

# Assign row and column names
rownames(mat) <- barcodes[, 1]   # gene IDs
colnames(mat) <- features[, 1]   # cell barcodes

# Convert to a dense matrix if desired (warning: may use a lot of memory)
# mat_dense <- as.matrix(mat)

# Check dimensions
dim(mat)

mat <- as.matrix(mat)

assigned_tag <- apply(mat, 1, function(x) {
  if (all(x == 0)) {
    return(NA)  # no tag detected
  } else {
    return(colnames(mat)[which.max(x)])
  }
})

# Combine with cell barcodes (rownames)
tag_assignment <- data.frame(
  cell_barcode = rownames(mat),
  assigned_tag = assigned_tag,
  stringsAsFactors = FALSE
)

tag_assignment$cell_barcode <- paste0(tag_assignment$cell_barcode, "-1")

# ==== subset seurat object ====
barcodes_to_keep <- tag_assignment$cell_barcode

seu_subset <- subset(seu, cells = barcodes_to_keep)

rownames(tag_assignment) <- tag_assignment$cell_barcode

meta_to_add <- tag_assignment[, "assigned_tag", drop = FALSE]

seu_subset <- AddMetaData(seu_subset, metadata = meta_to_add)

meta <- seu_subset@meta.data

qsave(seu_subset, file = "./data/raw_seu_with_hto_pair_tags.qs")

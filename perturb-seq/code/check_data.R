rm(list = ls())
set.seed(7)

library(dplyr)
library(Matrix)

# Specify file paths
barcodes_file <- "./data/HTO_raw_Run2/barcodes.tsv.gz"
features_file <- "./data/HTO_raw_Run2/features.tsv.gz"
matrix_file   <- "./data/HTO_raw_Run2/matrix.mtx.gz"

# Read the matrix (sparse format)
mat <- readMM(matrix_file)

# Read features and barcodes
features <- read.delim(features_file, header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.delim(barcodes_file, header = FALSE, stringsAsFactors = FALSE)

# Assign row and column names
rownames(mat) <- features[, 1]   # gene IDs
colnames(mat) <- barcodes[, 1]   # cell barcodes

# Convert to a dense matrix if desired (warning: may use a lot of memory)
# mat_dense <- as.matrix(mat)

# Check dimensions
dim(mat)

rownames_interested <- c("C0254_cmo", "C0255_cmo", "C0256_cmo", 
                         "C0257_cmo", "C0258_cmo", "C0259_cmo",
                         "C0254_hto", "C0255_hto", "C0256_hto",
                         "C0257_hto", "C0258_hto", "C0259_hto")

mat_sub <- mat[rownames_interested, ]





barcode <- read.csv("./code/PAIR_extraction/mRNA_barcodes.csv", header = FALSE)
barcode <- cbind("GRCh38", barcode)
write.table(
  barcode,
  file = "./code/PAIR_extraction/mRNA_barcodes_with_prefix.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

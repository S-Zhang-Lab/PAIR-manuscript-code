rm(list = ls())
set.seed(7)

library(dplyr)
library(Matrix)

# Specify file paths
barcodes_file <- "./data/HTO_raw/barcodes.tsv.gz"
features_file <- "./data/HTO_raw/features.tsv.gz"
matrix_file   <- "./data/HTO_raw/matrix.mtx.gz"

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

mat_sub <- mat[, 1:100] %>% as.matrix()

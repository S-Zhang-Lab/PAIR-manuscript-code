rm(list = ls())
set.seed(7)

library(dplyr)
library(Matrix)

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

hto_tab <- read.csv("./data/mRNA_HTO_Raw_Counts/tag_calls_per_cell.csv", 
                    header = TRUE)
rownames(hto_tab) <- sub("-1$", "", hto_tab$cell_barcode)

overlap_bc <- intersect(rownames(tag_assignment),
                        rownames(hto_tab))

tag_assignment_sub <- tag_assignment[overlap_bc, ]
hto_tab_sub <- hto_tab[overlap_bc, ]

hto_pair_tab <- cbind(tag_assignment_sub,
                     hto_tab_sub)

# Initialize treatment column as NA
hto_pair_tab$treatment <- NA

# Assign treatments by tag
hto_pair_tab$treatment[hto_pair_tab$feature_call %in% c("C0255_cmo", "C0256_cmo")] <- "RNP"
hto_pair_tab$treatment[hto_pair_tab$feature_call == "C0257_cmo"] <- "CTRL"




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

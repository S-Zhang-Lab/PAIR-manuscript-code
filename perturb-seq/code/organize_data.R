rm(list = ls())
set.seed(7)

library(dplyr)
library(Matrix)
library(Seurat)
library(qs)

# --- Source config.R (sets DATA_ROOT, SEU_RAW; see code/config.R) ----------
if (file.exists("code/config.R")) {
  source("code/config.R")                 # run from repo/
} else if (file.exists("config.R")) {
  source("config.R")                      # run from repo/code/
} else {
  stop("config.R not found. Run from repo/ or repo/code/.")
}
# ---------------------------------------------------------------------------

if (!exists("PAIR_ASSIGNMENT_MODE")) PAIR_ASSIGNMENT_MODE <- "legacy"
if (!exists("PAIR_MIN_TOP_UMI")) PAIR_MIN_TOP_UMI <- 2
if (!exists("PAIR_MIN_MARGIN_UMI")) PAIR_MIN_MARGIN_UMI <- 1
if (!exists("PAIR_MAX_SECOND_TO_TOP_RATIO")) PAIR_MAX_SECOND_TO_TOP_RATIO <- 0.5

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

ctl1_path <- file.path(DATA_DIR, "mRNA_HTO_Raw_Counts", "Ctl1")
rnp2_path <- file.path(DATA_DIR, "mRNA_HTO_Raw_Counts", "RNP2")
rnp3_path <- file.path(DATA_DIR, "mRNA_HTO_Raw_Counts", "RNP3")

legacy_sample_dirs_present <- all(dir.exists(c(ctl1_path, rnp2_path, rnp3_path)))

if (legacy_sample_dirs_present) {
  ctl1_seu <- load_sample(ctl1_path, "Ctl1")
  rnp2_seu <- load_sample(rnp2_path, "RNP2")
  rnp3_seu <- load_sample(rnp3_path, "RNP3")

  seu <- merge(ctl1_seu, y = c(rnp2_seu, rnp3_seu))
  meta <- seu@meta.data
  meta$treatment <- "RNP"
  meta$treatment[meta$sample == "Ctl1"] <- "CTRL"
  seu@meta.data <- meta
} else if (file.exists(SEU_RAW)) {
  cat("Legacy per-sample raw directories not found; using existing SEU_RAW as the HTO-demultiplexed base object.\n")
  seu <- qread(SEU_RAW)
} else {
  stop("Could not find either legacy sample directories or an existing SEU_RAW object to use as the base Seurat object.")
}

# ==== load PAIR results ====
# Specify file paths
# NOTE: Use the full PAIR matrix (6,227 barcodes), NOT the middle10 subset
# (6,211 barcodes).  The original script used middle10_rows/columns/matrix
# but seu_qc.qs was built from the full matrix; using middle10 here would
# drop 324 cells and cause a barcode mismatch.  See docs/REVIEW_SUMMARY.md
# §1 for the full explanation.
barcodes_file <- file.path(DATA_DIR, "PAIR_output", "PAIR_matched_sparse_UMI_matrix_rows.txt")
features_file <- file.path(DATA_DIR, "PAIR_output", "PAIR_matched_sparse_UMI_matrix_columns.txt")
matrix_file   <- file.path(DATA_DIR, "PAIR_output", "PAIR_matched_sparse_UMI_matrix.mtx")

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

assign_pair_tag <- function(x, tags) {
  total_pair_umi <- sum(x)
  if (total_pair_umi == 0) {
    return(data.frame(
      assigned_tag_legacy = NA_character_,
      top_tag = NA_character_,
      top_umi = 0,
      second_tag = NA_character_,
      second_umi = 0,
      assignment_margin = 0,
      second_to_top_ratio = NA_real_,
      total_pair_umi = 0,
      assignment_status = "zero_signal_unassigned",
      assigned_tag_confident = NA_character_,
      assigned_tag = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  ord <- order(x, decreasing = TRUE)
  top_idx <- ord[1]
  second_idx <- if (length(ord) >= 2) ord[2] else NA_integer_

  top_umi <- unname(x[top_idx])
  second_umi <- if (!is.na(second_idx)) unname(x[second_idx]) else 0
  top_tag <- tags[top_idx]
  second_tag <- if (!is.na(second_idx)) tags[second_idx] else NA_character_
  assignment_margin <- top_umi - second_umi
  second_to_top_ratio <- if (top_umi > 0) second_umi / top_umi else NA_real_
  assigned_tag_legacy <- top_tag

  assignment_status <- "assigned_confident"
  if (top_umi < PAIR_MIN_TOP_UMI) {
    assignment_status <- "low_umi_unassigned"
  } else if (assignment_margin < PAIR_MIN_MARGIN_UMI ||
             (!is.na(second_to_top_ratio) &&
              second_to_top_ratio > PAIR_MAX_SECOND_TO_TOP_RATIO)) {
    assignment_status <- "ambiguous_unassigned"
  }

  assigned_tag_confident <- if (assignment_status == "assigned_confident") top_tag else NA_character_
  assigned_tag <- if (identical(PAIR_ASSIGNMENT_MODE, "confident")) {
    assigned_tag_confident
  } else {
    assigned_tag_legacy
  }

  data.frame(
    assigned_tag_legacy = assigned_tag_legacy,
    top_tag = top_tag,
    top_umi = top_umi,
    second_tag = second_tag,
    second_umi = second_umi,
    assignment_margin = assignment_margin,
    second_to_top_ratio = second_to_top_ratio,
    total_pair_umi = total_pair_umi,
    assignment_status = assignment_status,
    assigned_tag_confident = assigned_tag_confident,
    assigned_tag = assigned_tag,
    stringsAsFactors = FALSE
  )
}

tag_assignment_list <- lapply(seq_len(nrow(mat)), function(i) {
  assign_pair_tag(mat[i, ], colnames(mat))
})

# Combine with cell barcodes (rownames)
tag_assignment <- bind_rows(tag_assignment_list)
tag_assignment$cell_barcode <- rownames(mat)

tag_assignment$cell_barcode <- paste0(tag_assignment$cell_barcode, "-1")

pair_output_dir <- file.path(DATA_DIR, "PAIR_output")
dir.create(pair_output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  tag_assignment,
  file.path(pair_output_dir, "PAIR_assignment_diagnostics.csv"),
  row.names = FALSE
)

assignment_qc_summary <- tag_assignment %>%
  count(assignment_status, name = "n_cells") %>%
  mutate(
    pct_cells = round(100 * n_cells / sum(n_cells), 2),
    assignment_mode = PAIR_ASSIGNMENT_MODE,
    min_top_umi = PAIR_MIN_TOP_UMI,
    min_margin_umi = PAIR_MIN_MARGIN_UMI,
    max_second_to_top_ratio = PAIR_MAX_SECOND_TO_TOP_RATIO
  )

write.csv(
  assignment_qc_summary,
  file.path(pair_output_dir, "PAIR_assignment_qc_summary.csv"),
  row.names = FALSE
)

# ==== subset seurat object ====
barcodes_to_keep <- tag_assignment$cell_barcode

seu_subset <- subset(seu, cells = barcodes_to_keep)

rownames(tag_assignment) <- tag_assignment$cell_barcode

meta_to_add <- tag_assignment[, c(
  "assigned_tag",
  "assigned_tag_legacy",
  "assigned_tag_confident",
  "top_tag",
  "top_umi",
  "second_tag",
  "second_umi",
  "assignment_margin",
  "second_to_top_ratio",
  "total_pair_umi",
  "assignment_status"
), drop = FALSE]

seu_subset <- AddMetaData(seu_subset, metadata = meta_to_add)

meta <- seu_subset@meta.data

dir.create(dirname(SEU_RAW), recursive = TRUE, showWarnings = FALSE)
qsave(seu_subset, file = SEU_RAW)

## =============================================================================
## 02_count_matrix_analysis.R
## Exploratory analysis of PAIR-DDR count matrix
## =============================================================================
##
## Purpose:
##   Perform initial exploratory analysis of the PAIR count matrix:
##   (1) Normalize raw counts by library size (CPM-like scaling)
##   (2) Compute log fold changes between sorted and unsorted cells
##   (3) Identify PAIRs with extreme fold changes in both replicates
##   (4) Generate a heatmap of concordant hits
##   (5) Reformat the count matrix for downstream MAGeCK input
##
## Input:
##   - ./Misc/kmer20_HD1_count_matrix.csv
##     (Raw count matrix with k-mer error correction, Hamming distance 1)
##
## Experimental design:
##   - DDR1 (unsorted, replicate 1) vs DDR2 (BFP sorted, replicate 1)
##   - DDR3 (unsorted, replicate 2) vs DDR4 (BFP sorted, replicate 2)
##   - DDR5 (post-puromycin selection baseline; excluded from DE analysis)
##
## Output:
##   - ./Results/PAIR_DDR_count_matrix.csv (normalized counts and fold changes)
##   - ./Results/Differential_Genes_DDR1_DDR2.csv
##   - ./Results/Differential_Genes_DDR3_DDR4.csv
##   - ./Results/Common_Positive_Genes.csv
##   - ./Results/Common_Negative_Genes.csv
##   - ./Figures/Heatmap_of_Common_Genes_kmer20.pdf
##   - ./Results/kmer20_HD1_count_matrix_for_mageck.txt (MAGeCK input)
##   - ./Results/02_count_matrix_analysis_sessionInfo.txt
##
## Note:
##   This script uses simple CPM normalization and arbitrary logFC thresholds
##   for initial exploration. Formal statistical testing is performed in
##   scripts 03 (DESeq2) and 04 (edgeR). MAGeCK is the primary analysis.
## =============================================================================

# --- Load required libraries --------------------------------------------------
library(dplyr)
library(pheatmap)

# --- 1. Load and clean count matrix -------------------------------------------

# Load error-corrected count matrix (k-mer mapping with Hamming distance 1)
count_matrix <- read.csv("./Misc/kmer20_HD1_count_matrix.csv")

# Remove DDR5 (post-puromycin baseline; not used in sorted vs unsorted comparison)
count_matrix$DDR5 <- NULL

# Create unique row identifiers from PAIR components
row.names(count_matrix) <- paste(count_matrix$CRISPRa_name,
                                  count_matrix$CasRx_name, sep = "_")
count_matrix$combined_name <- paste(count_matrix$CRISPRa_name,
                                     count_matrix$CasRx_name, sep = "_")

# --- 2. Normalize counts by library size -------------------------------------
# Simple CPM-like normalization: divide each count by column total, then
# scale to a common factor (100,000) for interpretability.
# This is a rough normalization for exploratory purposes only.

columns_to_normalize <- c("DDR1", "DDR2", "DDR3", "DDR4")
scaling_factor <- 100000

normalized_counts <- scale(
  count_matrix[columns_to_normalize],
  center = FALSE,
  scale  = colSums(count_matrix[columns_to_normalize], na.rm = TRUE)
) * scaling_factor

colnames(normalized_counts) <- paste0(columns_to_normalize, "_normalized")
count_matrix <- cbind(count_matrix, normalized_counts)

# --- 3. Compute log fold changes ----------------------------------------------
# Log2 fold change: sorted (treatment) / unsorted (control) within each replicate
# Pseudocount of 1 added to handle zeros

count_matrix$logFC_DDR1_DDR2 <- log2(count_matrix$DDR2_normalized + 1) -
                                 log2(count_matrix$DDR1_normalized + 1)
count_matrix$logFC_DDR3_DDR4 <- log2(count_matrix$DDR4_normalized + 1) -
                                 log2(count_matrix$DDR3_normalized + 1)

cat("Range of DDR2 normalized counts:", range(count_matrix$DDR2_normalized), "\n")
cat("Range of logFC (rep1):", range(count_matrix$logFC_DDR1_DDR2), "\n")
cat("Range of logFC (rep2):", range(count_matrix$logFC_DDR3_DDR4), "\n")

# Save full matrix with normalized counts and fold changes
write.csv(count_matrix, file = "./Results/PAIR_DDR_count_matrix.csv", row.names = FALSE)

# --- 4. Identify PAIRs with extreme fold changes -----------------------------
# NOTE: These thresholds are arbitrary and used for exploratory visualization.
# Formal statistical testing is done in scripts 03 and 04.

TH12 <- 6   # logFC threshold for replicate 1
TH34 <- 6   # logFC threshold for replicate 2

diff_genes_DDR1_DDR2 <- subset(count_matrix, abs(logFC_DDR1_DDR2) > TH12)
diff_genes_DDR3_DDR4 <- subset(count_matrix, abs(logFC_DDR3_DDR4) > TH34)

# --- 5. Find concordant PAIRs (same direction in both replicates) -------------

# Enriched in both replicates
common_positive <- intersect(
  rownames(subset(diff_genes_DDR1_DDR2, logFC_DDR1_DDR2 > 0)),
  rownames(subset(diff_genes_DDR3_DDR4, logFC_DDR3_DDR4 > 0))
)

# Depleted in both replicates
common_negative <- intersect(
  rownames(subset(diff_genes_DDR1_DDR2, logFC_DDR1_DDR2 < 0)),
  rownames(subset(diff_genes_DDR3_DDR4, logFC_DDR3_DDR4 < 0))
)

common_genes <- c(common_positive, common_negative)
cat("Concordant enriched PAIRs:", length(common_positive), "\n")
cat("Concordant depleted PAIRs:", length(common_negative), "\n")

# --- 6. Heatmap of concordant PAIRs -------------------------------------------

heatmap_data <- count_matrix[common_genes, paste0(columns_to_normalize, "_normalized")]

# Reorder columns: controls first (DDR1, DDR3), then treatments (DDR2, DDR4)
heatmap_data <- heatmap_data[, c("DDR1_normalized", "DDR3_normalized",
                                  "DDR2_normalized", "DDR4_normalized")]

# Row annotation: enriched vs depleted
annotation <- data.frame(
  Direction = c(rep("Enriched", length(common_positive)),
                rep("Depleted", length(common_negative)))
)
rownames(annotation) <- common_genes

# Generate heatmap (row-scaled)
p <- pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,    # Preserve manual column order
  scale = "row",
  annotation_row = annotation,
  main = "Heatmap of Concordant PAIRs (Exploratory)"
)

# Save heatmap
pdf("./Figures/Heatmap_of_Common_Genes_kmer20.pdf", width = 10, height = 10)
grid::grid.draw(p$gtable)
dev.off()

# --- 7. Export results ---------------------------------------------------------

write.csv(data.frame(Gene = row.names(diff_genes_DDR1_DDR2), diff_genes_DDR1_DDR2),
          "./Results/Differential_Genes_DDR1_DDR2.csv", row.names = FALSE)
write.csv(data.frame(Gene = row.names(diff_genes_DDR3_DDR4), diff_genes_DDR3_DDR4),
          "./Results/Differential_Genes_DDR3_DDR4.csv", row.names = FALSE)
write.csv(data.frame(Gene = common_positive),
          "./Results/Common_Positive_Genes.csv", row.names = FALSE)
write.csv(data.frame(Gene = common_negative),
          "./Results/Common_Negative_Genes.csv", row.names = FALSE)

# Display summary
list(
  positive_common = common_positive,
  negative_common = common_negative
)

# --- 8. Reformat count matrix for MAGeCK input --------------------------------
# MAGeCK expects: sgRNA | Gene | sample1 | sample2 | ...
# sgRNA = unique PAIR identifier
# Gene  = gene-level grouping (GeneA_GeneB)

data <- read.csv("./Misc/kmer20_HD1_count_matrix.csv")

data <- data %>%
  mutate(
    sgRNA = paste(CRISPRa_name, CasRx_name, sep = "_"),
    Gene  = paste(sub("_.*", "", CRISPRa_name),
                  sub("_.*", "", CasRx_name),
                  sep = "_")
  )

# Select columns by name for MAGeCK format
data <- data[, c("sgRNA", "Gene", "DDR1", "DDR2", "DDR3", "DDR4")]

# Remove rows where all count columns are zero (uninformative)
data1 <- data %>%
  filter(!if_all(where(is.numeric), ~ . == 0))

cat("PAIRs in MAGeCK input (after removing all-zero rows):", nrow(data1), "\n")

# Save as tab-delimited text for MAGeCK
write.table(data1, "./Results/kmer20_HD1_count_matrix_for_mageck.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- 9. Save session info for reproducibility ---------------------------------
writeLines(capture.output(sessionInfo()), "./Results/02_count_matrix_analysis_sessionInfo.txt")

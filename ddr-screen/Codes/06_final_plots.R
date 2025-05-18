# 06 final plot of the common Up and Down PAIRs using DE-seq2 DEG
# Load required libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# Load the data
data <- read.delim("./Results/kmer20_HD1_count_matrix_for_mageck.txt", header = TRUE, sep = "\t")

# Prepare count matrix and metadata
count_data <- as.matrix(data[, c("DDR1", "DDR2", "DDR3", "DDR4")])
rownames(count_data) <- data$sgRNA

# Filter out rows with three or more zero counts
count_data <- count_data[rowSums(count_data == 0) < 3, ]

# Define metadata with batch information
col_data <- data.frame(
  sample = colnames(count_data),
  condition = c("Control", "Treatment", "Control", "Treatment"),  # Define conditions
  batch = c("batch1", "batch1", "batch2", "batch2")               # Add batch information
)
rownames(col_data) <- col_data$sample

# Create DESeq2 dataset with batch in design
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ batch + condition
)

# Run DESeq2
dds <- DESeq(dds)

# Extract DESeq2 results
res <- results(dds)

# Filter significant DE genes (adjusted p-value < X, |log2FoldChange| > Y)
de_genes <- res[!is.na(res$padj) & res$padj < 0.001 & abs(res$log2FoldChange) > 15, ]
de_genes_df <- as.data.frame(de_genes)

# Save significant DE genes
write.csv(as.data.frame(de_genes), "./Results/DEseq2_Significant_DE_genes_with_batch.csv")

# Extract normalized counts for heatmap
normalized_counts <- counts(dds, normalized = TRUE)

# Log-transform normalized counts (add pseudocount to avoid log(0))
log_counts <- log2(normalized_counts + 1)

# Filter DE genes for heatmap
heatmap_data <- log_counts[rownames(de_genes), ]

# load the common PAIR list between DE-seq2 and EdgeR 
DE.seq_EdgeR_common_genes <- read.csv("./Results/DE-seq_EdgeR_common_genes.csv")

# Combine 'Up' and 'Down' columns into a single vector of row names to keep
rows_to_keep <- unique(c(DE.seq_EdgeR_common_genes$Up, DE.seq_EdgeR_common_genes$Down))

# Filter the heatmap_data by row names
heatmap_data <- heatmap_data[rownames(heatmap_data) %in% rows_to_keep, ]


# Prepare annotation data for columns
annotation_col <- data.frame(
  Condition = col_data$condition,
  Batch = col_data$batch
)
rownames(annotation_col) <- colnames(heatmap_data)  # Ensure row names match

# Specify colors for annotation
ann_colors = list(
  Condition = c(Control = "grey", Treatment = "darkgreen"),
  Batch = c(batch1 = "lightblue", batch2 = "lightcoral")
)

# Set PDF output file
pdf("./Figures/Final_heatmap_common_PAIR_DE_genes.pdf", width = 8, height = 6) 
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap of Significant DE Genes (Batch Corrected)",
  annotation_col = annotation_col,
  annotation_colors = ann_colors
)
dev.off()

pdf("./Figures/heatmap_common_PAIR_DE_genes_large.pdf", width = 10, height = 16) 
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap of Significant DE Genes (Batch Corrected)",
  annotation_col = annotation_col,
  annotation_colors = ann_colors
)
dev.off()

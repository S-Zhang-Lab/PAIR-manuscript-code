#### using DEseq2 for DE gene analysis ####

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
  condition = c("Control", "Treatment", "Control", "Treatment"),  # Define conditions Control = BFP-, Treatment = BFP+ HDR
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

# Order by decreasing log2FoldChange (upregulated)
top_up <- de_genes[order(-de_genes$log2FoldChange), ][1:20, ]

# Order by increasing log2FoldChange (downregulated)
top_down <- de_genes[order(de_genes$log2FoldChange), ][1:20, ]

# Combine top 20 in each direction
top_genes <- rbind(top_up, top_down)

# Ensure gene names match
heatmap_subset <- heatmap_data[rownames(heatmap_data) %in% rownames(top_genes), ]

# Extract normalized counts for heatmap
normalized_counts <- counts(dds, normalized = TRUE)

# Log-transform normalized counts (add pseudocount to avoid log(0))
log_counts <- log2(normalized_counts + 1)

# Filter DE genes for heatmap
heatmap_data <- log_counts[rownames(heatmap_subset), ]

# Prepare annotation data for columns
annotation_col <- data.frame(
  Condition = col_data$condition,
  Batch = col_data$batch
)
rownames(annotation_col) <- colnames(heatmap_data)  # Ensure row names match

# Specify colors for annotation
ann_colors = list(
  Condition = c(Control = "grey", Treatment = "blue"),
  Batch = c(batch1 = "lightgreen", batch2 = "lightcoral")
)

# Reorder columns by condition: Control first, then Treatment
control_first_order <- rownames(annotation_col)[annotation_col$Condition == "Control"]
treatment_order <- rownames(annotation_col)[annotation_col$Condition == "Treatment"]
ordered_cols <- c(control_first_order, treatment_order)

# Reorder heatmap data and annotation accordingly
heatmap_data <- heatmap_data[, ordered_cols]
annotation_col <- annotation_col[ordered_cols, , drop = FALSE]

# Generate PDF with pheatmap (disable column clustering to preserve order)
pdf("./Figures/DEseq2_heatmap_PAIR_DE_genes_large.pdf", width = 10, height = 10) 
pheatmap(
  heatmap_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # Disable to preserve your manual column order
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap of Significant DE Genes (Batch Corrected)",
  annotation_col = annotation_col,
  annotation_colors = ann_colors
)
dev.off()


# Generate volcano plot with top 20 gene labels
volcano_data <- as.data.frame(res)

# Add a column to identify significant genes
volcano_data$significant <- !is.na(volcano_data$padj) & 
  volcano_data$padj < 0.001 & 
  abs(volcano_data$log2FoldChange) > 17

# Sort the data by significance (-log10(padj)) and select the top 10 DE genes
top_genes <- volcano_data %>%
  arrange(padj) %>%
  head(20)

# Create the volcano plot
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of DE Genes (Batch Corrected)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-10, 10), linetype = "dashed", color = "blue") +
  # Add labels for the top 10 genes
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)), # Adjust if gene names are stored in another column
    size = 3,
    max.overlaps = 10
  )

# Save volcano plot
ggsave("./Figures/DEseq2_Volcano_Plot_with_Batch.pdf", width = 5, height = 5)
ggsave("./Figures/DEseq2_Volcano_Plot_with_Batch_large.pdf", width = 20, height = 20)



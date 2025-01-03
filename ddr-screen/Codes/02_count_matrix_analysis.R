# Load count matrix table
count_matrix <- read.csv("./Misc/count_matrix.csv")

# Step 1: Clean up dataframe
# Remove DDR5 as it is not required
count_matrix$DDR5 <- NULL

# Set row names using CRISPRa_name and CasRx_name
row.names(count_matrix) <- paste(count_matrix$CRISPRa_name, count_matrix$CasRx_name, sep = "_")
count_matrix$combined_name <- paste(count_matrix$CRISPRa_name, count_matrix$CasRx_name, sep = "_")
  
# Step 2: Normalize DDR columns by dividing by the total counts of each column
columns_to_normalize <- c("DDR1", "DDR2", "DDR3", "DDR4")

# Define scaling factor
scaling_factor <- 100000

# Create a new data frame for normalized values
normalized_counts <- scale(
  count_matrix[columns_to_normalize], 
  center = FALSE, 
  scale = colSums(count_matrix[columns_to_normalize], na.rm = TRUE)
) * scaling_factor

# Add scaled normalized counts as new columns
colnames(normalized_counts) <- paste0(columns_to_normalize, "_normalized")
count_matrix <- cbind(count_matrix, normalized_counts)

# Calculate log fold changes using scaled normalized values
count_matrix$logFC_DDR1_DDR2 <- log2(count_matrix$DDR2_normalized + 1) - log2(count_matrix$DDR1_normalized + 1)
count_matrix$logFC_DDR3_DDR4 <- log2(count_matrix$DDR4_normalized + 1) - log2(count_matrix$DDR3_normalized + 1)

range(count_matrix$DDR2_normalized)

# Step 3: Compare Control vs Treatment (DDR1 vs DDR2 and DDR3 vs DDR4) using normalized values
# Calculate log fold changes using normalized values
count_matrix$logFC_DDR1_DDR2 <- log2(count_matrix$DDR2_normalized + 1) - log2(count_matrix$DDR1_normalized + 1)
count_matrix$logFC_DDR3_DDR4 <- log2(count_matrix$DDR4_normalized + 1) - log2(count_matrix$DDR3_normalized + 1)

range(count_matrix$logFC_DDR1_DDR2)
range(count_matrix$logFC_DDR3_DDR4)

# write matrix
write.csv(count_matrix, file = "./Results/PAIR_DDR_count_matrix.csv", row.names = FALSE)


# Set thresholds for differential expression
TH12 <- 6
TH34 <- 6

# Identify differentially expressed rows
diff_genes_DDR1_DDR2 <- subset(count_matrix, abs(logFC_DDR1_DDR2) > TH12)
diff_genes_DDR3_DDR4 <- subset(count_matrix, abs(logFC_DDR3_DDR4) > TH34)

# Step 4: Identify common genes with consistent directional changes
# Positively regulated genes in both comparisons
common_positive <- intersect(
  rownames(subset(diff_genes_DDR1_DDR2, logFC_DDR1_DDR2 > 0)),
  rownames(subset(diff_genes_DDR3_DDR4, logFC_DDR3_DDR4 > 0))
)

# Negatively regulated (depleted) genes in both comparisons
common_negative <- intersect(
  rownames(subset(diff_genes_DDR1_DDR2, logFC_DDR1_DDR2 < 0)),
  rownames(subset(diff_genes_DDR3_DDR4, logFC_DDR3_DDR4 < 0))
)

# Step 5: Extract normalized counts for visualization
# Combine common genes for positive and negative into one list
common_genes <- c(common_positive, common_negative)

# Extract normalized counts for common genes
heatmap_data <- count_matrix[common_genes, paste0(columns_to_normalize, "_normalized")]

# Step 6: Plot heatmap using pheatmap
library(pheatmap)

# Define annotation for positive and negative genes
annotation <- data.frame(
  Direction = c(rep("Positive", length(common_positive)), rep("Negative", length(common_negative)))
)
rownames(annotation) <- common_genes

# Plot heatmap
p <- pheatmap(
  heatmap_data, 
  cluster_rows = TRUE,  # Cluster rows
  cluster_cols = TRUE,  # Cluster columns
  scale = "row",        # Scale data by row for better visualization
  annotation_row = annotation,
  main = "Heatmap of Common Genes"
)

p

# Save the heatmap to a PDF file
pdf("./Figures/Heatmap_of_Common_Genes.pdf", width = 10, height = 10)  # Set appropriate dimensions
grid::grid.draw(p$gtable)  # Draw the pheatmap object
dev.off()  # Close the PDF device

# Step 7: Export results (Optional)
write.csv(data.frame(Gene=row.names(diff_genes_DDR1_DDR2), diff_genes_DDR1_DDR2),
          "./Results/Differential_Genes_DDR1_DDR2.csv", row.names = FALSE)
write.csv(data.frame(Gene=row.names(diff_genes_DDR3_DDR4), diff_genes_DDR3_DDR4),
          "./Results/Differential_Genes_DDR3_DDR4.csv", row.names = FALSE)

write.csv(data.frame(Gene=common_positive), "./Results/Common_Positive_Genes.csv", row.names = FALSE)
write.csv(data.frame(Gene=common_negative), "./Results/Common_Negative_Genes.csv", row.names = FALSE)

# Step 8: Display results
list(
  positive_common = common_positive,
  negative_common = common_negative
)



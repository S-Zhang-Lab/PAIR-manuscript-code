## =============================================================================
## 03_DE-seq2.R
## PRIMARY: DESeq2 differential abundance analysis of PAIR-DDR screen
## =============================================================================
##
## Purpose:
##   Identify differentially enriched/depleted PAIRs using DESeq2 as the
##   PRIMARY parametric method for this screen. DESeq2's shrinkage-based
##   dispersion estimation is the most robust approach available for n=2
##   designs, partially compensating for low replication through empirical
##   Bayes moderation of per-gene variance. With a sparse count matrix,
##   fold changes may still be inflated for zero-to-nonzero transitions,
##   but the shrinkage prior mitigates this better than alternatives.
##   Results are used in 05_integration_analysis.R as the primary hit set.
##
## Input:
##   - ./Results/kmer20_HD1_count_matrix_for_mageck.txt
##
## Experimental design:
##   - DDR1 (unsorted, batch 1) vs DDR2 (BFP sorted, batch 1)
##   - DDR3 (unsorted, batch 2) vs DDR4 (BFP sorted, batch 2)
##   - Design formula: ~ batch + condition
##
## Output:
##   - ./Results/DEseq2_all_results.csv         (full result table)
##   - ./Results/DEseq2_Significant_DE_genes_with_batch.csv (filtered hits)
##   - ./Results/Normalized_counts.csv
##   - ./Figures/DEseq2_DE-gene_Heatmap.pdf
##   - ./Figures/DEseq2_Volcano_Plot_with_Batch.pdf
##   - ./Results/03_DEseq2_sessionInfo.txt
## =============================================================================

# --- Load required libraries --------------------------------------------------
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- Official HUGO gene symbol mapping ----------------------------------------
# Maps common/colloquial names (as used in sgRNA identifiers) to official HUGO
# symbols. Applied when extracting gene names for plot labels and annotations.
# Common name → HUGO: DNA-PKcs → PRKDC, 53BP1 → TP53BP1,
# KU70 → XRCC6, CtIP → RBBP8, XLF → NHEJ1.
gene_symbol_map <- c(
  "DNA-PKcs" = "PRKDC",
  "53BP1"    = "TP53BP1",
  "KU70"     = "XRCC6",
  "CtIP"     = "RBBP8",
  "XLF"      = "NHEJ1"
)

to_hugo <- function(x) {
  mapped <- gene_symbol_map[x]
  ifelse(is.na(mapped), x, mapped)
}

# --- Shared thresholds --------------------------------------------------------
# These thresholds are harmonized with edgeR (script 04) so the two methods
# are compared on equal footing in the integration analysis (script 05).
padj_threshold <- 0.01
lfc_threshold  <- 5

# --- 1. Load and prepare count data -------------------------------------------

data <- read.delim("./Results/kmer20_HD1_count_matrix_for_mageck.txt",
                    header = TRUE, sep = "\t")

count_data <- as.matrix(data[, c("DDR1", "DDR2", "DDR3", "DDR4")])
rownames(count_data) <- data$sgRNA

# Pre-filter: remove PAIRs with 3 or more zero counts across samples
count_data <- count_data[rowSums(count_data == 0) < 3, ]
cat("PAIRs retained after zero-count filter:", nrow(count_data), "\n")

# --- 2. Define experimental metadata ------------------------------------------

col_data <- data.frame(
  sample    = colnames(count_data),
  condition = factor(c("Control", "Treatment", "Control", "Treatment"),
                     levels = c("Control", "Treatment")),
  batch     = factor(c("batch1", "batch1", "batch2", "batch2"))
)
rownames(col_data) <- col_data$sample

# --- 3. Run DESeq2 with batch-corrected design --------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = col_data,
  design    = ~ batch + condition
)

dds <- DESeq(dds)
res <- results(dds)

# --- 4. Save full results (used by script 05 for integration) -----------------
res_df <- as.data.frame(res)
res_df$sgRNA <- rownames(res_df)
write.csv(res_df, "./Results/DEseq2_all_results.csv", row.names = FALSE)

# --- 5. Filter significant DE PAIRs -------------------------------------------

de_genes <- res[!is.na(res$padj) &
                  res$padj < padj_threshold &
                  abs(res$log2FoldChange) > lfc_threshold, ]
de_genes_df <- as.data.frame(de_genes)

cat("DESeq2 significant PAIRs (padj <", padj_threshold,
    ", |log2FC| >", lfc_threshold, "):", nrow(de_genes_df), "\n")

write.csv(de_genes_df, "./Results/DEseq2_Significant_DE_genes_with_batch.csv")

# --- 6. Normalized counts -----------------------------------------------------

normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "./Results/Normalized_counts.csv")
log_counts <- log2(normalized_counts + 1)

# --- 7. Heatmap ---------------------------------------------------------------
# Column order: controls (DDR1, DDR3) on left, treatments (DDR2, DDR4) on right

if (nrow(de_genes_df) > 1) {
  heatmap_data <- log_counts[rownames(de_genes), ]
  col_order <- c("DDR1", "DDR3", "DDR2", "DDR4")
  heatmap_data <- heatmap_data[, col_order]

  annotation_col <- data.frame(
    Condition = col_data[col_order, "condition"],
    Batch     = col_data[col_order, "batch"]
  )
  rownames(annotation_col) <- col_order

  ann_colors <- list(
    Condition = c(Control = "grey70", Treatment = "darkgreen"),
    Batch     = c(batch1 = "lightblue", batch2 = "lightcoral")
  )

  heatmap_obj <- pheatmap(
    heatmap_data,
    scale          = "row",
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    main           = "DESeq2: Significant DE PAIRs (Supplementary)",
    annotation_col = annotation_col,
    annotation_colors = ann_colors
  )

  pdf("./Figures/DEseq2_DE-gene_Heatmap.pdf", width = 10, height = 10)
  grid::grid.draw(heatmap_obj$gtable)
  dev.off()
} else {
  cat("Fewer than 2 significant PAIRs; skipping heatmap.\n")
}

# --- 8. Volcano plot -----------------------------------------------------------

volcano_data <- as.data.frame(res)
volcano_data$gene_name <- rownames(volcano_data)

volcano_data$significant <- !is.na(volcano_data$padj) &
  volcano_data$padj < padj_threshold &
  abs(volcano_data$log2FoldChange) > lfc_threshold

# Clean gene pair labels: GeneA_CRISPRa_N_GeneB_CasRx_N -> HUGO1_HUGO2
volcano_data$gene_label <- sapply(volcano_data$gene_name, function(x) {
  gene1 <- to_hugo(sub("_CRISPRa.*", "", x))
  gene2 <- sub(".*_CRISPRa_\\d+_", "", x)
  gene2 <- to_hugo(sub("_CasRx_.*", "", gene2))
  paste(gene1, gene2, sep = "_")
})

top_genes <- volcano_data %>%
  filter(significant == TRUE) %>%
  arrange(padj) %>%
  head(20)

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(volcano_data, !significant),
             color = "gray70", alpha = 0.5, size = 1) +
  geom_point(data = subset(volcano_data, significant),
             color = "red", alpha = 0.7, size = 1.5) +
  geom_hline(yintercept = -log10(padj_threshold),
             linetype = "dashed", color = "blue", linewidth = 0.5) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
             linetype = "dashed", color = "blue", linewidth = 0.5) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_label),
    size = 3, max.overlaps = 15,
    segment.color = "gray40", segment.size = 0.3, fontface = "italic"
  ) +
  labs(
    title = "DESeq2: Volcano Plot of Significantly Changed PAIR",
    x     = "Log2 Fold Change (Treatment vs Control)",
    y     = expression(-Log[10] ~ "Adjusted P-value")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(color = "#C71585", face = "bold", size = 14, hjust = 0.5),
    legend.position  = "none",
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = min(volcano_data$log2FoldChange, na.rm = TRUE) * 0.6,
           y = max(-log10(volcano_data$padj), na.rm = TRUE) * 0.95,
           label = expression("" %<-% " BFP"^"-"),
           size = 4, color = "black", hjust = 0) +
  annotate("text", x = max(volcano_data$log2FoldChange, na.rm = TRUE) * 0.6,
           y = max(-log10(volcano_data$padj), na.rm = TRUE) * 0.95,
           label = expression("BFP"^"+" ~ "" %->% ""),
           size = 4, color = "black", hjust = 1)

print(volcano_plot)
ggsave("./Figures/DEseq2_Volcano_Plot_with_Batch.pdf",
       plot = volcano_plot, width = 8, height = 8)

# --- 9. Session info -----------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./Results/03_DEseq2_sessionInfo.txt")

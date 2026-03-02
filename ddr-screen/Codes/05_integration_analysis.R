## =============================================================================
## 05_integration_analysis.R
## PRIMARY: DESeq2-based hit integration with edgeR/MAGeCK cross-validation
## =============================================================================
##
## Purpose:
##   Integrate results from the three differential abundance methods and
##   perform gene-level frequency analysis with a formal null model.
##
##   Analysis hierarchy:
##     PRIMARY   – DESeq2 (script 03): Best-available parametric method for
##                 n=2. Dispersion shrinkage partially compensates for low
##                 replication. padj < 0.01, |log2FC| > 5.
##     SECONDARY – edgeR (script 04): GLM with batch correction + replicate
##                 concordance filter. FDR < 0.05, |logFC| > 5.
##     REFERENCE – MAGeCK (HPC): Included for completeness but its FDR
##                 calibration is unreliable for this focused library
##                 (10,525 PAIRs, limited sgRNA diversity per gene pair).
##                 MAGeCK calls 83% of tested PAIRs as significant at
##                 FDR < 0.05, indicating miscalibrated p-values.
##
##   Gene-frequency analysis uses DESeq2-significant PAIRs (the primary
##   hit set) and applies a hypergeometric test against library composition
##   as the null model.
##
## Input:
##   - ./Results/DEseq2_all_results.csv                 (from script 03)
##   - ./Results/edgeR_all_results.csv                  (from script 04)
##   - ./Results/MAGeCK/DDR_default.sgrna_summary.txt   (MAGeCK, optional)
##   - ./Results/MAGeCK/DDR_default.gene_summary.txt    (MAGeCK, optional)
##   - ./Results/kmer20_HD1_count_matrix_for_mageck.txt  (count matrix)
##
## Output:
##   - ./Figures/DESeq2_Primary_Heatmap.pdf
##   - ./Figures/Cross_Validation_Summary.pdf
##   - ./Figures/Gene_Frequency_Enriched.pdf
##   - ./Figures/Gene_Frequency_Depleted.pdf
##   - ./Results/DESeq2_primary_significant_PAIRs.csv
##   - ./Results/cross_validation_summary.csv
##   - ./Results/gene_frequency_enriched.csv
##   - ./Results/gene_frequency_depleted.csv
##   - ./Results/05_integration_sessionInfo.txt
## =============================================================================

# --- Load required libraries --------------------------------------------------
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)

# --- Official HUGO gene symbol mapping ----------------------------------------
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

# --- Configuration ------------------------------------------------------------
# DESeq2 thresholds (primary)
deseq2_padj_threshold <- 0.01
deseq2_lfc_threshold  <- 5

# edgeR thresholds (secondary; relaxed FDR due to n=2 + batch correction)
edger_fdr_threshold <- 0.05
edger_lfc_threshold <- 5

# MAGeCK thresholds (reference only)
mageck_fdr_threshold <- 0.05

# =============================================================================
# SECTION 1: Load results from all three methods
# =============================================================================

# --- 1. DESeq2 results (PRIMARY) ---------------------------------------------

deseq2_all <- read.csv("./Results/DEseq2_all_results.csv")
cat("DESeq2 results loaded:", nrow(deseq2_all), "PAIRs\n")

deseq2_sig <- deseq2_all %>%
  filter(!is.na(padj) & padj < deseq2_padj_threshold &
         abs(log2FoldChange) > deseq2_lfc_threshold)

# Assign direction based on fold change sign
deseq2_sig$direction <- ifelse(deseq2_sig$log2FoldChange > 0,
                                "enriched", "depleted")

cat("DESeq2 significant PAIRs (padj <", deseq2_padj_threshold,
    ", |log2FC| >", deseq2_lfc_threshold, "):", nrow(deseq2_sig), "\n")
cat("  Enriched (BFP+):", sum(deseq2_sig$direction == "enriched"), "\n")
cat("  Depleted (BFP-):", sum(deseq2_sig$direction == "depleted"), "\n")

write.csv(deseq2_sig, "./Results/DESeq2_primary_significant_PAIRs.csv",
          row.names = FALSE)

# --- 2. edgeR results (SECONDARY) --------------------------------------------

edger_all <- read.csv("./Results/edgeR_all_results.csv")
cat("\nedgeR results loaded:", nrow(edger_all), "PAIRs\n")

edger_sig <- edger_all %>%
  filter(FDR < edger_fdr_threshold & abs(logFC) > edger_lfc_threshold)

cat("edgeR significant PAIRs (FDR <", edger_fdr_threshold,
    ", |logFC| >", edger_lfc_threshold, "):", nrow(edger_sig), "\n")

# --- 3. MAGeCK results (REFERENCE, optional) ----------------------------------

mageck_sgrna_file <- "./Results/MAGeCK/DDR_default.sgrna_summary.txt"
mageck_gene_file  <- "./Results/MAGeCK/DDR_default.gene_summary.txt"
mageck_available  <- file.exists(mageck_sgrna_file)

if (mageck_available) {
  mageck_sgrna <- read.delim(mageck_sgrna_file, header = TRUE, sep = "\t")
  cat("\nMAGeCK results loaded:", nrow(mageck_sgrna), "PAIRs\n")

  mageck_sgrna$direction <- ifelse(
    mageck_sgrna$high_in_treatment == "True", "enriched", "depleted"
  )
  mageck_sgrna$significant <- mageck_sgrna$FDR < mageck_fdr_threshold

  cat("MAGeCK significant PAIRs (FDR <", mageck_fdr_threshold, "):",
      sum(mageck_sgrna$significant), "\n")
  cat("  NOTE: MAGeCK calls", round(mean(mageck_sgrna$significant) * 100),
      "% of PAIRs as significant — likely miscalibrated for this library.\n")
} else {
  cat("\nMAGeCK results not found; cross-validation will use DESeq2 and edgeR only.\n")
}

# =============================================================================
# SECTION 2: DESeq2 primary heatmap
# =============================================================================

# --- 4. Heatmap of DESeq2-significant PAIRs -----------------------------------

count_data <- read.delim("./Results/kmer20_HD1_count_matrix_for_mageck.txt",
                          header = TRUE, sep = "\t")
count_mat <- as.matrix(count_data[, c("DDR1", "DDR2", "DDR3", "DDR4")])
rownames(count_mat) <- count_data$sgRNA

if (nrow(deseq2_sig) > 1) {
  # Log2(CPM + 1) for visualization
  lib_sizes <- colSums(count_mat)
  cpm_mat <- t(t(count_mat) / lib_sizes * 1e6)
  log_cpm <- log2(cpm_mat + 1)

  sig_pairs_in_mat <- intersect(deseq2_sig$sgRNA, rownames(log_cpm))

  if (length(sig_pairs_in_mat) > 1) {
    heatmap_data <- log_cpm[sig_pairs_in_mat, ]

    # Column order: controls left, treatments right
    col_order <- c("DDR1", "DDR3", "DDR2", "DDR4")
    heatmap_data <- heatmap_data[, col_order]

    annotation_col <- data.frame(
      Condition = factor(c("Control", "Control", "Treatment", "Treatment"),
                         levels = c("Control", "Treatment")),
      Batch = factor(c(1, 2, 1, 2)),
      row.names = col_order
    )

    ann_colors <- list(
      Condition = c(Control = "grey70", Treatment = "darkgreen"),
      Batch     = c("1" = "lightblue", "2" = "lightcoral")
    )

    heatmap_obj <- pheatmap(
      heatmap_data,
      scale          = "row",
      cluster_rows   = TRUE,
      cluster_cols   = FALSE,
      show_rownames  = (nrow(heatmap_data) <= 80),
      show_colnames  = TRUE,
      main           = "DESeq2: Significant DE PAIRs (Primary Analysis)",
      annotation_col = annotation_col,
      annotation_colors = ann_colors
    )

    pdf("./Figures/DESeq2_Primary_Heatmap.pdf", width = 10, height = 10)
    grid::grid.draw(heatmap_obj$gtable)
    dev.off()
  } else {
    cat("Fewer than 2 DESeq2-significant PAIRs in count matrix; skipping heatmap.\n")
  }
} else {
  cat("Fewer than 2 DESeq2-significant PAIRs; skipping heatmap.\n")
}

# =============================================================================
# SECTION 3: Cross-validation
# =============================================================================

# --- 5. Build cross-validation summary ----------------------------------------

deseq2_set <- deseq2_sig$sgRNA
edger_set  <- edger_sig$sgRNA

if (mageck_available) {
  mageck_set <- mageck_sgrna$sgrna[mageck_sgrna$significant]

  cat("\n--- Cross-validation summary ---\n")
  cat("DESeq2 significant (primary):", length(deseq2_set), "\n")
  cat("edgeR significant (secondary):", length(edger_set), "\n")
  cat("MAGeCK significant (reference):", length(mageck_set), "\n")

  deseq2_edger  <- intersect(deseq2_set, edger_set)
  deseq2_mageck <- intersect(deseq2_set, mageck_set)
  edger_mageck  <- intersect(edger_set, mageck_set)
  all_three     <- Reduce(intersect, list(deseq2_set, edger_set, mageck_set))

  cat("DESeq2 & edgeR:", length(deseq2_edger), "\n")
  cat("DESeq2 & MAGeCK:", length(deseq2_mageck), "\n")
  cat("edgeR & MAGeCK:", length(edger_mageck), "\n")
  cat("All three:", length(all_three), "\n")

  cross_val <- data.frame(
    Category = c("DESeq2 only", "edgeR only", "MAGeCK only",
                  "DESeq2 & edgeR", "DESeq2 & MAGeCK",
                  "edgeR & MAGeCK", "All three"),
    Count = c(
      length(setdiff(deseq2_set, union(edger_set, mageck_set))),
      length(setdiff(edger_set, union(deseq2_set, mageck_set))),
      length(setdiff(mageck_set, union(deseq2_set, edger_set))),
      length(setdiff(deseq2_edger, mageck_set)),
      length(setdiff(deseq2_mageck, edger_set)),
      length(setdiff(edger_mageck, deseq2_set)),
      length(all_three)
    )
  )
} else {
  cat("\n--- Cross-validation summary (DESeq2 vs edgeR) ---\n")
  cat("DESeq2 significant (primary):", length(deseq2_set), "\n")
  cat("edgeR significant (secondary):", length(edger_set), "\n")

  deseq2_edger <- intersect(deseq2_set, edger_set)
  cat("DESeq2 & edgeR:", length(deseq2_edger), "\n")

  cross_val <- data.frame(
    Category = c("DESeq2 only", "edgeR only", "DESeq2 & edgeR"),
    Count = c(
      length(setdiff(deseq2_set, edger_set)),
      length(setdiff(edger_set, deseq2_set)),
      length(deseq2_edger)
    )
  )
}

write.csv(cross_val, "./Results/cross_validation_summary.csv", row.names = FALSE)

# --- 6. Cross-validation bar chart --------------------------------------------

cross_val$Category <- factor(cross_val$Category, levels = cross_val$Category)

# Color palette
if (mageck_available) {
  bar_colors <- c("#377EB8", "#4DAF4A", "#E41A1C",
                   "#984EA3", "#FF7F00", "#A65628", "#F781BF")
} else {
  bar_colors <- c("#377EB8", "#4DAF4A", "#984EA3")
}

cross_plot <- ggplot(cross_val, aes(x = Category, y = Count)) +
  geom_bar(stat = "identity", fill = bar_colors[seq_len(nrow(cross_val))]) +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  labs(
    title = "Cross-validation of Significant PAIRs",
    subtitle = paste0("DESeq2 (padj<", deseq2_padj_threshold,
                      ", |LFC|>", deseq2_lfc_threshold,
                      ") | edgeR (FDR<", edger_fdr_threshold,
                      ", |LFC|>", edger_lfc_threshold, ")",
                      if (mageck_available) paste0(" | MAGeCK (FDR<",
                        mageck_fdr_threshold, ")") else ""),
    x = NULL, y = "Number of significant PAIRs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40")
  )

print(cross_plot)
ggsave("./Figures/Cross_Validation_Summary.pdf",
       plot = cross_plot, width = 8, height = 6)

# =============================================================================
# SECTION 4: Gene-level frequency analysis (based on DESeq2 primary hits)
# =============================================================================

# --- 7. Extract gene identities from PAIR names ------------------------------
# PAIR names: GeneA_CRISPRa_N_GeneB_CasRx_N
# Position 1 (CRISPRa target, transcriptionally activated) = GeneA
# Position 2 (CasRx target, transcriptionally knocked down) = GeneB

extract_genes <- function(pair_names) {
  p1 <- to_hugo(sub("_CRISPRa.*", "", pair_names))
  p2 <- sub(".*_CRISPRa_\\d+_", "", pair_names)
  p2 <- to_hugo(sub("_CasRx_.*", "", p2))
  data.frame(pair = pair_names, position1 = p1, position2 = p2,
             stringsAsFactors = FALSE)
}

# --- 8. Build library composition (background) --------------------------------
# The null expectation depends on how many PAIRs in the full library contain
# each gene. A gene represented by more PAIRs has a higher chance of appearing
# in any random hit list.
#
# We use the full count matrix (all 10,525 PAIRs) as the background, not the
# filtered subset, because the library composition is a property of the
# experimental design, not of the statistical analysis.

all_pairs <- count_data$sgRNA
all_genes <- extract_genes(all_pairs)

lib_p1_counts <- table(all_genes$position1)
lib_p2_counts <- table(all_genes$position2)
n_total <- length(all_pairs)

cat("\nLibrary composition:", n_total, "total PAIRs,",
    length(lib_p1_counts), "unique CRISPRa genes,",
    length(lib_p2_counts), "unique CasRx genes\n")

# --- 9. Hypergeometric gene frequency test ------------------------------------
# For each gene at each position, test whether its frequency among DESeq2-
# significant PAIRs exceeds what is expected given its library representation.
#
# phyper(q-1, m, n, k, lower.tail=FALSE) = P(X >= q)
#   q = observed count for this gene among hits
#   m = total PAIRs containing this gene in the library
#   n = total PAIRs NOT containing this gene
#   k = total significant PAIRs (draws without replacement)

compute_gene_freq <- function(sig_pair_names, direction_label) {
  if (length(sig_pair_names) == 0) {
    return(data.frame(Gene = character(0), Position = character(0),
                      Observed = integer(0), Expected = numeric(0),
                      Library_total = integer(0), Fold_enrichment = numeric(0),
                      pvalue = numeric(0), padj = numeric(0),
                      Direction = character(0)))
  }

  sig_genes <- extract_genes(sig_pair_names)
  k <- length(sig_pair_names)

  # Position 1 (CRISPRa) frequencies
  p1_freq <- as.data.frame(table(sig_genes$position1), stringsAsFactors = FALSE)
  colnames(p1_freq) <- c("Gene", "Observed")
  p1_freq$Library_total <- as.integer(lib_p1_counts[p1_freq$Gene])
  p1_freq$Expected <- p1_freq$Library_total * k / n_total
  p1_freq$Fold_enrichment <- p1_freq$Observed / pmax(p1_freq$Expected, 1e-10)
  p1_freq$pvalue <- phyper(
    q = p1_freq$Observed - 1,
    m = p1_freq$Library_total,
    n = n_total - p1_freq$Library_total,
    k = k,
    lower.tail = FALSE
  )
  p1_freq$Position <- "Position 1 (CRISPRa)"

  # Position 2 (CasRx) frequencies
  p2_freq <- as.data.frame(table(sig_genes$position2), stringsAsFactors = FALSE)
  colnames(p2_freq) <- c("Gene", "Observed")
  p2_freq$Library_total <- as.integer(lib_p2_counts[p2_freq$Gene])
  p2_freq$Expected <- p2_freq$Library_total * k / n_total
  p2_freq$Fold_enrichment <- p2_freq$Observed / pmax(p2_freq$Expected, 1e-10)
  p2_freq$pvalue <- phyper(
    q = p2_freq$Observed - 1,
    m = p2_freq$Library_total,
    n = n_total - p2_freq$Library_total,
    k = k,
    lower.tail = FALSE
  )
  p2_freq$Position <- "Position 2 (CasRx)"

  combined <- rbind(p1_freq, p2_freq)
  combined$padj <- p.adjust(combined$pvalue, method = "BH")
  combined$Direction <- direction_label
  combined[order(combined$padj), ]
}

# Separate enriched and depleted DESeq2 hits
sig_enriched <- deseq2_sig %>%
  filter(direction == "enriched") %>%
  pull(sgRNA)

sig_depleted <- deseq2_sig %>%
  filter(direction == "depleted") %>%
  pull(sgRNA)

freq_enriched <- compute_gene_freq(sig_enriched, "Enriched (BFP+)")
freq_depleted <- compute_gene_freq(sig_depleted, "Depleted (BFP-)")

write.csv(freq_enriched, "./Results/gene_frequency_enriched.csv", row.names = FALSE)
write.csv(freq_depleted, "./Results/gene_frequency_depleted.csv", row.names = FALSE)

cat("\n--- Gene frequency analysis (DESeq2 enriched PAIRs) ---\n")
cat("Genes significantly over-represented (hypergeometric padj < 0.05):\n")
sig_genes_enriched <- freq_enriched %>% filter(padj < 0.05)
if (nrow(sig_genes_enriched) > 0) {
  print(sig_genes_enriched[, c("Gene", "Position", "Observed", "Expected",
                                 "Fold_enrichment", "padj")])
} else {
  cat("  None\n")
}

cat("\n--- Gene frequency analysis (DESeq2 depleted PAIRs) ---\n")
cat("Genes significantly over-represented (hypergeometric padj < 0.05):\n")
sig_genes_depleted <- freq_depleted %>% filter(padj < 0.05)
if (nrow(sig_genes_depleted) > 0) {
  print(sig_genes_depleted[, c("Gene", "Position", "Observed", "Expected",
                                 "Fold_enrichment", "padj")])
} else {
  cat("  None\n")
}

# --- 10. Gene frequency bar plots with expected baseline ----------------------

plot_gene_freq <- function(freq_df, title_label) {
  if (nrow(freq_df) == 0) {
    cat("No genes to plot for:", title_label, "\n")
    return(NULL)
  }

  freq_df$sig_label <- ifelse(freq_df$padj < 0.05, "*", "")
  freq_df$sig_label[freq_df$padj < 0.01] <- "**"
  freq_df$sig_label[freq_df$padj < 0.001] <- "***"

  p <- ggplot(freq_df, aes(x = reorder(Gene, -Observed), y = Observed)) +
    geom_bar(stat = "identity", aes(fill = padj < 0.05), show.legend = TRUE) +
    geom_point(aes(y = Expected), color = "red", size = 2, shape = 4) +
    geom_text(aes(label = sig_label), vjust = -0.3, size = 4) +
    scale_fill_manual(
      values = c("FALSE" = "gray70", "TRUE" = "steelblue"),
      labels = c("Not significant", "padj < 0.05"),
      name = "Hypergeometric test"
    ) +
    facet_wrap(~Position, scales = "free", ncol = 1) +
    labs(
      title = title_label,
      subtitle = "Red × = expected count given library composition",
      x = "Gene", y = "Number of DESeq2-significant PAIRs containing gene"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
    )

  p
}

p_enriched <- plot_gene_freq(freq_enriched,
  "Gene Frequency Among DESeq2-Enriched PAIRs (BFP+)")
if (!is.null(p_enriched)) {
  print(p_enriched)
  ggsave("./Figures/Gene_Frequency_Enriched.pdf",
         plot = p_enriched, width = 10, height = 10)
}

p_depleted <- plot_gene_freq(freq_depleted,
  "Gene Frequency Among DESeq2-Depleted PAIRs (BFP-)")
if (!is.null(p_depleted)) {
  print(p_depleted)
  ggsave("./Figures/Gene_Frequency_Depleted.pdf",
         plot = p_depleted, width = 10, height = 10)
}

# =============================================================================
# SECTION 5: MAGeCK gene-level summary (reference, if available)
# =============================================================================

if (mageck_available && file.exists(mageck_gene_file)) {
  mageck_gene <- read.delim(mageck_gene_file, header = TRUE, sep = "\t")
  cat("\n--- MAGeCK gene-level summary (reference) ---\n")
  cat("Gene pairs tested:", nrow(mageck_gene), "\n")

  gene_pos_sig <- mageck_gene %>% filter(pos.fdr < mageck_fdr_threshold)
  gene_neg_sig <- mageck_gene %>% filter(neg.fdr < mageck_fdr_threshold)

  cat("Positively selected gene pairs (FDR <", mageck_fdr_threshold, "):",
      nrow(gene_pos_sig), "\n")
  cat("Negatively selected gene pairs (FDR <", mageck_fdr_threshold, "):",
      nrow(gene_neg_sig), "\n")

  if (nrow(gene_pos_sig) > 0) {
    cat("\nTop positively selected gene pairs:\n")
    print(head(gene_pos_sig[order(gene_pos_sig$pos.fdr),
                             c("id", "pos.score", "pos.fdr")], 20))
  }
  if (nrow(gene_neg_sig) > 0) {
    cat("\nTop negatively selected gene pairs:\n")
    print(head(gene_neg_sig[order(gene_neg_sig$neg.fdr),
                             c("id", "neg.score", "neg.fdr")], 20))
  }
}

# --- 11. Session info ---------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           "./Results/05_integration_sessionInfo.txt")

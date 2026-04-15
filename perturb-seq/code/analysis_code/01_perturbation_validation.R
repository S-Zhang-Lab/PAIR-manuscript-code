##############################################################################
# Step 01: Perturbation Validation (Dual-Directional)
# PAIR-Perturb-Seq Analysis Pipeline
#
# Input:  data/seu_qc.qs
# Output: res/01_validation/ (DotPlots, VlnPlots, cell counts, UMAPs)
#
# Validates that CRISPRa activates NBN and CasRx suppresses partners
# at the transcriptional level.
##############################################################################

library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

set.seed(42)

# --- Source config.R (sets DATA_ROOT, RES_DIR; see code/config.R) -----------
if (file.exists("code/config.R")) {
  source("code/config.R")                 # run from repo/
} else if (file.exists("../config.R")) {
  source("../config.R")                   # run from repo/code/analysis_code/
} else {
  stop("config.R not found. Run from repo/ or repo/code/analysis_code/.")
}
# ---------------------------------------------------------------------------

base_dir <- DATA_ROOT
out_dir <- file.path(OUTPUT_DIR, "01_validation")
fig_dir <- file.path(FIGURES_DIR, "01_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load QC-filtered object ---
cat("Loading seu_qc.qs...\n")
seu <- qread(SEU_QC)
cat("Cells:", ncol(seu), "| Genes:", nrow(seu), "\n")

# --- Define target genes and conditions ---
target_genes <- c("NBN", "XRCC6", "TP53BP1", "POLQ", "PRKDC")

# Create simplified condition label
seu$condition <- paste(seu$treatment, seu$assigned_tag, sep = "_")

# Master tag → HUGO display-name map
tag_hugo_display <- c(
  "NBN_CRISPRa_53BP1_CasRx"  = "NBN_CRISPRa_TP53BP1_CasRx",
  "NBN_CRISPRa_KU70_CasRx"   = "NBN_CRISPRa_XRCC6_CasRx",
  "NBN_CRISPRa_PQ_CasRx"     = "NBN_CRISPRa_POLQ_CasRx",
  "NT_CRISPRa_53BP1_CasRx"   = "NT_CRISPRa_TP53BP1_CasRx",
  "NT_CRISPRa_KU70_CasRx"    = "NT_CRISPRa_XRCC6_CasRx",
  "NT_CRISPRa_PQ_CasRx"      = "NT_CRISPRa_POLQ_CasRx",
  "NT_CRISPRa_DNPKcs_CasRx"  = "NT_CRISPRa_PRKDC_CasRx",
  "NBN_CRISPRa_DNPKcs_CasRx" = "NBN_CRISPRa_PRKDC_CasRx"
)
seu$display_tag <- as.character(seu$assigned_tag)
for (old in names(tag_hugo_display)) {
  seu$display_tag[seu$display_tag == old] <- tag_hugo_display[old]
}

# --- 1. Cell Count Summary ---
cat("\n=== Cell Counts by Tag x Treatment ===\n")
count_tab <- as.data.frame(table(seu$display_tag, seu$treatment))
colnames(count_tab) <- c("Tag", "Treatment", "Count")
count_wide <- count_tab %>%
  pivot_wider(names_from = Treatment, values_from = Count, values_fill = 0) %>%
  mutate(Total = CTRL + RNP) %>%
  arrange(desc(Total))
print(count_wide)
write.csv(count_wide, file.path(out_dir, "cell_counts_summary.csv"), row.names = FALSE)

# Bar chart
p_bar <- ggplot(count_tab, aes(x = reorder(Tag, -Count), y = Count, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("CTRL" = "#4DBBD5", "RNP" = "#E64B35")) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(x = "PAIR Tag", y = "Cell Count", title = "Cell Counts by Tag and Treatment")
ggsave(file.path(fig_dir, "cell_counts_barplot.pdf"), p_bar, width = 10, height = 6)

# --- 2. UMAP Plots ---

# UMAP by tag
p_umap1 <- DimPlot(seu, group.by = "display_tag", pt.size = 1.5, label = FALSE) +
  ggtitle("UMAP: PAIR Tag Assignment") +
  theme(legend.text = element_text(size = 10))

# UMAP by treatment
p_umap2 <- DimPlot(seu, group.by = "treatment", pt.size = 1.5) +
  scale_color_manual(values = c("CTRL" = "#4DBBD5", "RNP" = "#E64B35")) +
  ggtitle("UMAP: Treatment")

p_umap <- p_umap1 | p_umap2
ggsave(file.path(fig_dir, "umap_overview.pdf"), p_umap, width = 16, height = 7)

# --- 3. DotPlot: Target Gene Expression ---
# Subset to key NBN-axis conditions for clarity
nbn_conditions <- c(
  "NBN_CRISPRa_NT_CasRx",
  "NBN_CRISPRa_53BP1_CasRx",
  "NBN_CRISPRa_KU70_CasRx",
  "NBN_CRISPRa_PQ_CasRx",
  "NT_CRISPRa_NT_CasRx"
)

seu_sub <- subset(seu, assigned_tag %in% nbn_conditions)

# Create readable condition labels
condition_labels <- c(
  "NBN_CRISPRa_NT_CasRx"    = "NBN + NT",
  "NBN_CRISPRa_53BP1_CasRx" = "NBN + TP53BP1 KD",
  "NBN_CRISPRa_KU70_CasRx"  = "NBN + XRCC6 KD",
  "NBN_CRISPRa_PQ_CasRx"    = "NBN + POLQ KD",
  "NT_CRISPRa_NT_CasRx"     = "NT + NT (Control)"
)
seu_sub$label <- unname(condition_labels[as.character(seu_sub$assigned_tag)])

# Condition order: NT+NT control first, then NBN+NT, then dual combinations
# In DotPlot + coord_flip: level 1 appears at LEFT
cond_levels <- c("NT + NT (Control)", "NBN + NT", "NBN + TP53BP1 KD", "NBN + XRCC6 KD", "NBN + POLQ KD")
seu_sub$label <- factor(seu_sub$label, levels = cond_levels)

# DotPlot
Idents(seu_sub) <- "label"
p_dot <- DotPlot(seu_sub, features = target_genes, dot.scale = 8) +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        aspect.ratio = 1) +
  ggtitle("Target Gene Expression by Condition")
ggsave(file.path(fig_dir, "dotplot_target_genes.pdf"), p_dot, width = 10, height = 7)

# Split by treatment (CTRL before RNP within each condition)
label_treat_levels <- as.vector(outer(cond_levels, c("CTRL", "RNP"), paste, sep = " | "))
seu_sub$label_treat <- paste(as.character(seu_sub$label), as.character(seu_sub$treatment), sep = " | ")
seu_sub$label_treat <- factor(seu_sub$label_treat, levels = label_treat_levels)
Idents(seu_sub) <- "label_treat"
p_dot_split <- DotPlot(seu_sub, features = target_genes, dot.scale = 8) +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 12),
        aspect.ratio = 1) +
  ggtitle("Target Gene Expression by Condition x Treatment")
ggsave(file.path(fig_dir, "dotplot_target_genes_by_treatment.pdf"), p_dot_split,
       width = 14, height = 8)

# --- 4. VlnPlot: Target Gene Expression ---
Idents(seu_sub) <- "label"
p_vln <- VlnPlot(seu_sub, features = target_genes, split.by = "treatment",
                  pt.size = 0, ncol = 3) &
  theme(axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 13))
p_vln <- p_vln + plot_annotation(title = "Target Gene Expression: CTRL vs RNP")
ggsave(file.path(fig_dir, "vlnplot_target_genes.pdf"), p_vln, width = 16, height = 12)

# --- 5. Quantitative Summary Table ---
# Compute percent expressing and mean expression per condition per gene
meta <- seu_sub@meta.data
expr_mat <- GetAssayData(seu_sub, layer = "data")

quant_list <- list()
for (gene in target_genes) {
  if (!gene %in% rownames(expr_mat)) next
  gene_expr <- as.numeric(expr_mat[gene, ])
  for (cond in unique(meta$label)) {
    for (treat in c("CTRL", "RNP")) {
      idx <- which(meta$label == cond & meta$treatment == treat)
      if (length(idx) == 0) next
      vals <- gene_expr[idx]
      quant_list[[length(quant_list) + 1]] <- data.frame(
        gene = gene,
        condition = cond,
        treatment = treat,
        n_cells = length(vals),
        pct_expressing = round(100 * mean(vals > 0), 1),
        mean_expression = round(mean(vals), 4)
      )
    }
  }
}
quant_df <- do.call(rbind, quant_list)
write.csv(quant_df, file.path(out_dir, "target_gene_quantification.csv"), row.names = FALSE)

cat("\n=== Target Gene Expression Summary ===\n")
print(quant_df %>% filter(gene == "NBN") %>% arrange(condition, treatment))

# --- 6. PAIR Assignment Diagnostics -------------------------------------------
# Critique §D: document the argmax rule, zero-UMI fate, tie handling, margin
# distribution, and double-positive fraction so assignment robustness is
# transparent.  Reads the raw PAIR sparse matrix directly.
cat("\n=== PAIR Assignment Diagnostics ===\n")

pair_dir  <- file.path(DATA_DIR, "PAIR_output")
# NOTE: Code discovery — the current organize_data.R points at the "middle10"
# PAIR subset (6211 barcodes), but the existing seu_qc.qs was built when the
# FULL matrix (6227 barcodes) was used: all 2271 QC-passing cells match the
# full matrix, while 324 do not match middle10. Diagnostics therefore use the
# full matrix for complete, accurate coverage. If organize_data.R is re-run
# with the current middle10 code, ~324 cells would be lost from seu_qc.
# This code/data inconsistency is documented in DATA_ANALYSIS_GUIDE.md Part II.
rows_file <- file.path(pair_dir, "PAIR_matched_sparse_UMI_matrix_rows.txt")
cols_file <- file.path(pair_dir, "PAIR_matched_sparse_UMI_matrix_columns.txt")
mtx_file  <- file.path(pair_dir, "PAIR_matched_sparse_UMI_matrix.mtx")

if (!file.exists(mtx_file)) {
  cat("PAIR matrix not found at", mtx_file, "— skipping diagnostics.\n")
} else {
  library(Matrix)

  pair_barcodes <- read.delim(rows_file, header = FALSE, stringsAsFactors = FALSE)$V1
  pair_tags     <- read.delim(cols_file, header = FALSE, stringsAsFactors = FALSE)$V1
  pair_mat      <- readMM(mtx_file)          # rows = cells, columns = tags
  rownames(pair_mat) <- pair_barcodes
  colnames(pair_mat) <- pair_tags
  pair_mat_dense <- as.matrix(pair_mat)      # 6211 x 11

  cat("Raw PAIR matrix dimensions:", nrow(pair_mat_dense), "cells x",
      ncol(pair_mat_dense), "tags\n")
  cat("Tags:", paste(pair_tags, collapse = ", "), "\n")

  # Per-cell summary statistics
  total_umi   <- rowSums(pair_mat_dense)
  max_umi     <- apply(pair_mat_dense, 1, max)
  runner_up   <- apply(pair_mat_dense, 1, function(x) sort(x, decreasing = TRUE)[2])
  margin      <- max_umi - runner_up
  n_nonzero   <- rowSums(pair_mat_dense > 0)

  diag_df <- data.frame(
    cell_barcode   = paste0(pair_barcodes, "-1"),   # match seu_qc format
    total_pair_umi = total_umi,
    max_tag_umi    = max_umi,
    runner_up_umi  = runner_up,
    margin         = margin,
    n_nonzero_tags = n_nonzero,
    zero_umi_cell  = total_umi == 0,
    stringsAsFactors = FALSE
  )

  # Annotate with QC status + assigned tag from seu_qc
  qc_meta <- seu@meta.data[, c("sample", "treatment", "assigned_tag"), drop = FALSE]
  qc_meta$cell_barcode <- rownames(qc_meta)
  qc_meta$in_seu_qc    <- TRUE
  diag_df <- merge(diag_df, qc_meta, by = "cell_barcode", all.x = TRUE)
  diag_df$in_seu_qc[is.na(diag_df$in_seu_qc)] <- FALSE

  write.csv(diag_df, file.path(out_dir, "pair_assignment_diagnostics.csv"),
            row.names = FALSE)

  # ── Print summary ──────────────────────────────────────────────────────────
  n_total   <- nrow(diag_df)
  n_zero    <- sum(diag_df$zero_umi_cell)
  n_qc      <- sum(diag_df$in_seu_qc)
  n_single  <- sum(!diag_df$zero_umi_cell & diag_df$n_nonzero_tags == 1)
  n_multi   <- sum(!diag_df$zero_umi_cell & diag_df$n_nonzero_tags > 1)

  diag_qc <- diag_df[diag_df$in_seu_qc, ]
  med_mar  <- median(diag_qc$margin, na.rm = TRUE)
  n_low_mg <- sum(diag_qc$margin < 3, na.rm = TRUE)

  cat(sprintf("Total cells in raw PAIR matrix : %d\n",   n_total))
  cat(sprintf("  Zero-UMI cells (-> NA -> dropped): %d (%.1f%%)\n",
              n_zero, 100 * n_zero / n_total))
  cat(sprintf("  Single-tag cells (clean assign): %d (%.1f%%)\n",
              n_single, 100 * n_single / n_total))
  cat(sprintf("  Multi-tag cells (>=2 non-zero): %d (%.1f%%)\n",
              n_multi, 100 * n_multi / n_total))
  cat(sprintf("Cells passing QC (in seu_qc)    : %d\n",   n_qc))
  cat(sprintf("  Median top/runner-up margin   : %.0f UMIs\n", med_mar))
  cat(sprintf("  Cells with margin < 3 UMIs    : %d (%.1f%% of QC cells)\n",
              n_low_mg, 100 * n_low_mg / n_qc))

  # ── Figure 1: Total PAIR UMI per cell (all cells, capped at 200) ──────────
  library(ggplot2)
  p_hist_total <- ggplot(
    diag_df %>% filter(!zero_umi_cell),
    aes(x = pmin(total_pair_umi, 200))
  ) +
    geom_histogram(bins = 60, fill = "#4DBBD5", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = median(diag_df$total_pair_umi[!diag_df$zero_umi_cell]),
               linetype = "dashed", color = "#E64B35", linewidth = 0.8) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 100, 150, 200),
      labels = c("0", "25", "50", "100", "150", "≥200")
    ) +
    labs(
      title = "Total PAIR UMI per cell (non-zero cells)",
      subtitle = sprintf("n = %d cells; %d zero-UMI cells excluded",
                         n_total - n_zero, n_zero),
      x = "Total PAIR UMI (capped at 200)",
      y = "Number of cells",
      caption = "Dashed red line = median"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "pair_diag_total_umi_hist.pdf"),
         p_hist_total, width = 7, height = 4.5)

  # ── Figure 2: Top-vs-runner-up margin distribution (non-zero cells) ───────
  p_hist_margin <- ggplot(
    diag_df %>% filter(!zero_umi_cell),
    aes(x = pmin(margin, 100))
  ) +
    geom_histogram(bins = 50, fill = "#91D1C2", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 3, linetype = "dashed", color = "#E64B35",
               linewidth = 0.8) +
    scale_x_continuous(
      breaks = c(0, 3, 10, 25, 50, 100),
      labels = c("0", "3", "10", "25", "50", "≥100")
    ) +
    labs(
      title = "PAIR assignment margin (top tag − runner-up UMI)",
      subtitle = sprintf("n = %d non-zero cells; margin = 0 means single non-zero tag",
                         n_total - n_zero),
      x = "Margin (capped at 100)",
      y = "Number of cells",
      caption = "Dashed red line = margin 3 (ambiguity threshold)"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "pair_diag_margin_hist.pdf"),
         p_hist_margin, width = 7, height = 4.5)

  # ── Figure 3: Number of non-zero tags per cell ────────────────────────────
  nz_tab <- as.data.frame(table(n_nonzero_tags = diag_df$n_nonzero_tags))
  nz_tab$n_nonzero_tags <- as.integer(as.character(nz_tab$n_nonzero_tags))
  nz_tab$pct <- 100 * nz_tab$Freq / n_total
  nz_tab$label_col <- ifelse(nz_tab$n_nonzero_tags == 0, "#cccccc",
                      ifelse(nz_tab$n_nonzero_tags == 1, "#4DBBD5", "#E64B35"))

  p_bar_nz <- ggplot(nz_tab, aes(x = factor(n_nonzero_tags), y = pct,
                                  fill = label_col)) +
    geom_col(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 3.2) +
    scale_fill_identity() +
    labs(
      title = "Number of non-zero PAIR tags per cell",
      subtitle = "Grey=zero-UMI (dropped), blue=single tag, red=multi-tag",
      x = "# Non-zero PAIR tags", y = "% of all cells"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "pair_diag_nonzero_tags_bar.pdf"),
         p_bar_nz, width = 6, height = 4)

  # ── Figure 4: QC-passing cells — margin stratified by partner arm ─────────
  if (n_qc > 0) {
    p_vio_margin <- ggplot(
      diag_qc %>% filter(!is.na(assigned_tag)),
      aes(x = reorder(assigned_tag, margin, median), y = pmin(margin, 100),
          fill = assigned_tag)
    ) +
      geom_violin(alpha = 0.7, trim = TRUE) +
      geom_boxplot(width = 0.08, outlier.shape = NA, color = "grey30") +
      coord_flip() +
      labs(
        title = "PAIR assignment margin per tag arm (QC-passing cells)",
        subtitle = "margin = max UMI − runner-up UMI; values capped at 100",
        x = NULL, y = "Margin (capped at 100)"
      ) +
      theme_classic(base_size = 11) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold"))
    ggsave(file.path(fig_dir, "pair_diag_margin_by_arm.pdf"),
           p_vio_margin, width = 8, height = 6)
    cat("Saved: pair_diag_margin_by_arm.pdf\n")
  }

  cat("Saved: pair_assignment_diagnostics.csv\n")
  cat("Saved: pair_diag_total_umi_hist.pdf\n")
  cat("Saved: pair_diag_margin_hist.pdf\n")
  cat("Saved: pair_diag_nonzero_tags_bar.pdf\n")
}

cat("\nStep 01 complete.\n")
cat("Tables written to:", out_dir, "\n")
cat("Figures written to:", fig_dir, "\n")

##############################################################################
# Step 02: Tier 1 - NBN Activation Baseline DE
# PAIR-Perturb-Seq Analysis Pipeline
#
# Input:  data/seu_qc.qs
# Output: res/02_tier1/ (DE table, volcano plot)
#
# Question: What transcriptomic changes does NBN CRISPRa activation
# cause at baseline (without exogenous stress)?
# Comparison: CTRL_NBN_NT vs CTRL_NT_NT
##############################################################################

library(Seurat)
library(qs)
library(ggplot2)
library(dplyr)
library(ggrepel)

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
out_dir <- file.path(OUTPUT_DIR, "02_tier1")
fig_dir <- file.path(FIGURES_DIR, "02_tier1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading seu_qc.qs...\n")
seu <- qread(SEU_QC)

# --- Define groups ---
# NBN activation group: CTRL cells with NBN_CRISPRa_NT_CasRx
# Control group: CTRL cells with NT_CRISPRa_NT_CasRx
seu_ctrl <- subset(seu, treatment == "CTRL")

group1_tag <- "NBN_CRISPRa_NT_CasRx"
group2_tag <- "NT_CRISPRa_NT_CasRx"

n_group1 <- sum(seu_ctrl$assigned_tag == group1_tag)
n_group2 <- sum(seu_ctrl$assigned_tag == group2_tag)
cat("Group 1 (NBN_NT, CTRL):", n_group1, "cells\n")
cat("Group 2 (NT_NT, CTRL):", n_group2, "cells\n")

if (n_group2 < 10) {
  cat("WARNING: Control group has fewer than 10 cells. Results may be unreliable.\n")
}

# --- Run DE ---
Idents(seu_ctrl) <- "assigned_tag"

cat("Running FindMarkers (Wilcoxon)...\n")
de_results <- FindMarkers(
  seu_ctrl,
  ident.1 = group1_tag,
  ident.2 = group2_tag,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0, # Keep all genes (don't pre-filter)
  verbose = FALSE
)

de_results$gene <- rownames(de_results)
de_results <- de_results %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

cat("Total genes tested:", nrow(de_results), "\n")
cat("Significant (adj.P < 0.05):", sum(de_results$p_val_adj < 0.05), "\n")
cat(
  "Significant (adj.P < 0.05, |FC| > 0.5):",
  sum(de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 0.5), "\n"
)

# Save full results
write.csv(de_results, file.path(out_dir, "Tier1_NBN_baseline_DE.csv"), row.names = FALSE)


de_results <- read.csv(file.path(out_dir, "Tier1_NBN_baseline_DE.csv"))


# --- Volcano Plot ---
# Classify genes
de_results$class <- "Not Sig"
de_results$class[de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 0.5] <- "Up in NBN"
de_results$class[de_results$p_val_adj < 0.05 & de_results$avg_log2FC < -0.5] <- "Down in NBN"

# Handle zero p-values for -log10 plotting
min_nonzero_p <- min(de_results$p_val_adj[de_results$p_val_adj > 0], na.rm = TRUE) # Use p_val_adj
de_results$p_plot <- ifelse(de_results$p_val_adj == 0, min_nonzero_p * 0.1, de_results$p_val_adj)

# Top genes to label
top_up <- de_results %>%
  filter(class == "Up in NBN") %>%
  top_n(10, wt = -p_val_adj)
top_down <- de_results %>%
  filter(class == "Down in NBN") %>%
  top_n(10, wt = -p_val_adj)
# If no significant genes, label top by raw p-value
if (nrow(top_up) + nrow(top_down) == 0) {
  top_up <- de_results %>%
    filter(avg_log2FC > 0) %>%
    top_n(5, wt = -p_val)
  top_down <- de_results %>%
    filter(avg_log2FC < 0) %>%
    top_n(5, wt = -p_val)
}
label_genes <- rbind(top_up, top_down)

p_volcano <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_plot), color = class)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c(
    "Up in NBN" = "#E64B35", "Down in NBN" = "#4DBBD5",
    "Not Sig" = "grey70"
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = label_genes, aes(label = gene),
    size = 3.5, max.overlaps = 20, color = "black"
  ) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 1) +
  labs(
    title = "Tier 1: NBN Activation Baseline (CTRL only)",
    subtitle = paste0("NBN_NT vs NT_NT | 370 vs 41 cells"),
    x = "Log2 Fold Change", y = "-Log10(P-value)", color = ""
  )

ggsave(file.path(fig_dir, "Tier1_Volcano.pdf"), p_volcano, width = 7, height = 7)

# --- Summary ---
cat("\n=== Tier 1 Summary ===\n")
cat("Comparison: CTRL_NBN_NT (n=", n_group1, ") vs CTRL_NT_NT (n=", n_group2, ")\n")
cat("Total genes tested:", nrow(de_results), "\n")
cat(
  "Significant (adj.P < 0.05 & |FC| > 0.5):",
  sum(de_results$class != "Not Sig"), "\n"
)
cat("  Up in NBN:", sum(de_results$class == "Up in NBN"), "\n")
cat("  Down in NBN:", sum(de_results$class == "Down in NBN"), "\n")

# Check NBN itself
if ("NBN" %in% de_results$gene) {
  nbn_row <- de_results[de_results$gene == "NBN", ]
  cat(
    "\nNBN gene: log2FC =", round(nbn_row$avg_log2FC, 3),
    "| adj.P =", format(nbn_row$p_val_adj, digits = 3), "\n"
  )
}

cat("\nStep 02 complete.\n")
cat("Tables written to:", out_dir, "\n")
cat("Figures written to:", fig_dir, "\n")

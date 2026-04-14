##############################################################################
# Step 07: Tier 5 - AUCell Pathway Activity Mapping
# PAIR-Perturb-Seq Analysis Pipeline
#
# Input:  data/seu_qc.qs, data/pathways/hallmark_pathways.rds
# Output: output/07_tier5/ (AUCell activity matrices)
#         figures/07_tier5/ (AUCell dot plots)
#
# Uses rank-based AUCell scoring to compute per-cell pathway activity,
# then summarizes as per-condition dot plots showing:
#   - Color: z-scored mean AUC
#   - Size: percent cells active (above global median)
##############################################################################

library(Seurat)
library(qs)
library(AUCell)
library(ggplot2)
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
out_dir <- file.path(OUTPUT_DIR, "07_tier5")
fig_dir <- file.path(FIGURES_DIR, "07_tier5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading seu_qc.qs...\n")
seu <- qread(file.path(base_dir, "data", "seu_qc.qs"))

# --- Load gene sets ---
# Hallmark pathways
hallmark_path <- file.path(base_dir, "data", "pathways", "hallmark_pathways.rds")
if (file.exists(hallmark_path)) {
  pathways <- readRDS(hallmark_path)
} else {
  library(msigdbr)
  msig_df <- msigdbr(species = "Homo sapiens", collection = "H")
  pathways <- split(msig_df$gene_symbol, msig_df$gs_name)
  saveRDS(pathways, hallmark_path)
}
cat("Loaded", length(pathways), "Hallmark pathways\n")

# Custom repair pathways
repair_path <- file.path(OUTPUT_DIR, "06_tier4", "repair_module_gene_sets.rds")
if (file.exists(repair_path)) {
  repair_sets <- readRDS(repair_path)
  # Prefix to distinguish from Hallmark
  names(repair_sets) <- paste0("REPAIR_", names(repair_sets))
  pathways <- c(pathways, repair_sets)
  cat("Added", length(repair_sets), "custom repair pathway sets\n")
}

# --- Build AUCell Rankings ---
cat("Building AUCell rankings...\n")
expr_mat <- GetAssayData(seu, layer = "data")
cells_rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE, verbose = FALSE)

# --- Compute AUC Scores ---
cat("Computing AUC scores for", length(pathways), "gene sets...\n")
cells_AUC <- AUCell_calcAUC(pathways, cells_rankings, verbose = FALSE)
auc_mat <- getAUC(cells_AUC)

cat("AUC matrix:", nrow(auc_mat), "pathways x", ncol(auc_mat), "cells\n")

# --- Define conditions ---
nbn_conditions <- c(
  "NBN_CRISPRa_NT_CasRx", "NBN_CRISPRa_53BP1_CasRx",
  "NBN_CRISPRa_KU70_CasRx", "NBN_CRISPRa_PQ_CasRx",
  "NT_CRISPRa_NT_CasRx"
)

condition_labels <- c(
  "NBN_CRISPRa_NT_CasRx"    = "NBN+NT",
  "NBN_CRISPRa_53BP1_CasRx" = "NBN+TP53BP1 KD",
  "NBN_CRISPRa_KU70_CasRx"  = "NBN+XRCC6 KD",
  "NBN_CRISPRa_PQ_CasRx"    = "NBN+POLQ KD",
  "NT_CRISPRa_NT_CasRx"     = "NT+NT (Ctrl)"
)

# Condition order: NT+NT control first, NBN+NT second, then dual combinations
# Treatment order: CTRL before RNP within each condition
cond_order <- c("NT+NT (Ctrl)", "NBN+NT", "NBN+TP53BP1 KD", "NBN+XRCC6 KD", "NBN+POLQ KD")
cond_treat_levels <- as.vector(outer(cond_order, c("CTRL", "RNP"), paste, sep = " | "))

# Subset to NBN-axis conditions
keep_cells <- colnames(seu)[seu$assigned_tag %in% nbn_conditions]
auc_sub <- auc_mat[, keep_cells]

# Get metadata for these cells
meta <- seu@meta.data[keep_cells, ]
meta$condition <- unname(condition_labels[as.character(meta$assigned_tag)])
meta$cond_treat <- paste(meta$condition, meta$treatment, sep = " | ")

# --- Per-condition summary statistics ---
cat("Computing per-condition AUC summaries...\n")

conditions <- unique(meta$cond_treat)
summary_list <- list()

for (cond in conditions) {
  cells <- rownames(meta)[meta$cond_treat == cond]
  if (length(cells) < 5) next

  auc_cond <- auc_sub[, cells, drop = FALSE]
  mean_auc <- rowMeans(auc_cond)
  n_cells <- length(cells)

  summary_list[[cond]] <- data.frame(
    condition = cond,
    pathway = names(mean_auc),
    mean_auc = as.numeric(mean_auc),
    n_cells = n_cells
  )
}

summary_df <- do.call(rbind, summary_list)

# Compute global median AUC per pathway (for percent active)
global_median <- apply(auc_sub, 1, median)

# Compute percent active per condition
pct_active_list <- list()
for (cond in conditions) {
  cells <- rownames(meta)[meta$cond_treat == cond]
  if (length(cells) < 5) next

  auc_cond <- auc_sub[, cells, drop = FALSE]
  pct <- rowMeans(auc_cond > global_median[rownames(auc_cond)]) * 100

  pct_active_list[[cond]] <- data.frame(
    condition = cond,
    pathway = names(pct),
    pct_active = as.numeric(pct)
  )
}
pct_df <- do.call(rbind, pct_active_list)

# Merge
dot_df <- merge(summary_df, pct_df, by = c("condition", "pathway"))

# Z-score mean AUC per pathway across conditions
dot_df <- dot_df %>%
  group_by(pathway) %>%
  mutate(z_auc = (mean_auc - mean(mean_auc)) / sd(mean_auc)) %>%
  ungroup()

# Replace NaN z-scores (pathways with 0 variance)
dot_df$z_auc[is.nan(dot_df$z_auc)] <- 0

write.csv(dot_df, file.path(out_dir, "aucell_summary.csv"), row.names = FALSE)

# --- Dot Plot: Hallmark pathways ---
cat("Generating dot plots...\n")

# Filter to Hallmark pathways only
hallmark_df <- dot_df %>% filter(grepl("^HALLMARK_", pathway))

# Select pathways with most variance across conditions
pathway_var <- hallmark_df %>%
  group_by(pathway) %>%
  summarise(var_z = var(z_auc, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(var_z))

top_pathways <- head(pathway_var$pathway, 25)
hallmark_plot <- hallmark_df %>% filter(pathway %in% top_pathways)

# Clean pathway names
hallmark_plot$pathway_clean <- gsub("HALLMARK_", "", hallmark_plot$pathway)
hallmark_plot$pathway_clean <- gsub("_", " ", hallmark_plot$pathway_clean)

# Apply correct condition order (NT+NT first, CTRL before RNP)
existing_levels <- cond_treat_levels[cond_treat_levels %in% unique(hallmark_plot$condition)]
hallmark_plot$condition <- factor(hallmark_plot$condition, levels = existing_levels)

p_dot <- ggplot(hallmark_plot, aes(x = condition, y = pathway_clean)) +
  geom_point(aes(size = pct_active, color = z_auc)) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#E64B35",
    midpoint = 0, name = "Z-scored\nMean AUC"
  ) +
  scale_size_continuous(range = c(3, 12), name = "% Active") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11)
  ) +
  labs(
    title = "AUCell: Top Variable Hallmark Pathways",
    subtitle = "Size = % cells active (above global median), Color = z-scored mean AUC",
    x = "", y = ""
  )

ggsave(file.path(fig_dir, "aucell_dotplot_hallmark.pdf"), p_dot,
  width = 16, height = max(7, length(top_pathways) * 0.4 + 3)
)

# --- Dot Plot: RNP only (stress conditions) ---
rnp_df <- hallmark_df %>%
  filter(grepl("RNP", condition)) %>%
  filter(pathway %in% top_pathways)

rnp_df$pathway_clean <- gsub("HALLMARK_", "", rnp_df$pathway)
rnp_df$pathway_clean <- gsub("_", " ", rnp_df$pathway_clean)

# Apply correct condition order
rnp_levels <- cond_treat_levels[grepl("RNP", cond_treat_levels)]
rnp_levels <- rnp_levels[rnp_levels %in% unique(rnp_df$condition)]
rnp_df$condition <- factor(rnp_df$condition, levels = rnp_levels)

p_dot_rnp <- ggplot(rnp_df, aes(x = condition, y = pathway_clean)) +
  geom_point(aes(size = pct_active, color = z_auc)) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#E64B35",
    midpoint = 0, name = "Z-scored\nMean AUC"
  ) +
  scale_size_continuous(range = c(3, 12), name = "% Active") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11)
  ) +
  labs(
    title = "AUCell: Hallmark Pathways (RNP Stress Only)",
    x = "", y = ""
  )

ggsave(file.path(fig_dir, "aucell_dotplot_hallmark_RNP.pdf"), p_dot_rnp,
  width = 11, height = max(7, length(top_pathways) * 0.4 + 3)
)

# --- Dot Plot: Custom repair pathways ---
repair_df <- dot_df %>% filter(grepl("^REPAIR_", pathway))

if (nrow(repair_df) > 0) {
  repair_df$pathway_clean <- gsub("REPAIR_", "", repair_df$pathway)

  existing_repair_levels <- cond_treat_levels[cond_treat_levels %in% unique(repair_df$condition)]
  repair_df$condition <- factor(repair_df$condition, levels = existing_repair_levels)

  p_dot_repair <- ggplot(repair_df, aes(x = condition, y = pathway_clean)) +
    geom_point(aes(size = pct_active, color = z_auc)) +
    scale_color_gradient2(
      low = "#2166AC", mid = "white", high = "#E64B35",
      midpoint = 0, name = "Z-scored\nMean AUC"
    ) +
    scale_size_continuous(range = c(3, 12), name = "% Active") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 12)
    ) +
    labs(
      title = "AUCell: DNA Repair Module Activity",
      x = "", y = ""
    )

  ggsave(file.path(fig_dir, "aucell_dotplot_repair.pdf"), p_dot_repair, width = 14, height = 6)
}

cat("\nStep 07 complete. Tables in:", out_dir, "| Figures in:", fig_dir, "\n")

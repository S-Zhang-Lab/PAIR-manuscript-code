##############################################################################
# Step 08: Phase 2 - Epistasis Analysis (Modules A-D)
# PAIR-Perturb-Seq Analysis Pipeline
#
# Input:  data/seu_qc.qs
# Output: output/08_phase2/ (epistasis tables)
#         figures/08_phase2/ (heatmaps, PCA manifold, scatter plots)
#
# Module A: Gene-level additive model (FC residuals)
# Module B: Pathway-level tau scores (AUCell-based)
# Module C: Transcriptomic PCA manifold (pseudobulk)
# Module D: Gatekeeper hierarchy (rescue percentages)
##############################################################################

library(Seurat)
library(qs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(AUCell)

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
out_dir <- file.path(OUTPUT_DIR, "08_phase2")
fig_dir <- file.path(FIGURES_DIR, "08_phase2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading seu_qc.qs...\n")
seu <- qread(file.path(base_dir, "data", "seu_qc.qs"))

# --- Define conditions ---
# Key conditions for epistasis:
# Reference:    RNP_NT_NT          (stressed WT)
# CRISPRa only: RNP_NBN_NT         (stressed + NBN activation)
# Dual:         RNP_NBN_TP53BP1    (stressed + NBN + TP53BP1 KD)
# Dual:         RNP_NBN_XRCC6      (stressed + NBN + XRCC6 KD)
# Dual:         RNP_NBN_POLQ       (stressed + NBN + POLQ KD)

conditions <- c(
  "NT_CRISPRa_NT_CasRx",
  "NBN_CRISPRa_NT_CasRx",
  "NBN_CRISPRa_53BP1_CasRx",
  "NBN_CRISPRa_KU70_CasRx",
  "NBN_CRISPRa_PQ_CasRx"
)

condition_short <- c(
  "NT_CRISPRa_NT_CasRx"     = "NT_NT",
  "NBN_CRISPRa_NT_CasRx"    = "NBN_NT",
  "NBN_CRISPRa_53BP1_CasRx" = "NBN_TP53BP1",
  "NBN_CRISPRa_KU70_CasRx"  = "NBN_XRCC6",
  "NBN_CRISPRa_PQ_CasRx"    = "NBN_POLQ"
)

partners <- c("TP53BP1", "XRCC6", "POLQ")

##############################################################################
# MODULE A: Gene-level Additive Model
##############################################################################
cat("\n=== MODULE A: Gene-level Additive Model ===\n")

# Compute pseudobulk using log-of-means (correct method)
# AverageExpression returns mean of raw counts when using "RNA" counts layer
# Then compute log2FC = log2(mean_test / mean_ref + pseudocount)

# Helper: compute pseudobulk log2FC correctly
compute_log2fc <- function(seu, cells_test, cells_ref, pseudocount = 1) {
  expr <- GetAssayData(seu, layer = "counts")
  mean_test <- rowMeans(expr[, cells_test, drop = FALSE])
  mean_ref <- rowMeans(expr[, cells_ref, drop = FALSE])
  log2((mean_test + pseudocount) / (mean_ref + pseudocount))
}

# Get cells per condition x treatment
get_cells <- function(seu, tag, treatment) {
  colnames(seu)[seu$assigned_tag == tag & seu$treatment == treatment]
}

# Focus on RNP-treated cells for stress-response epistasis
cells_rnp_nt_nt <- get_cells(seu, "NT_CRISPRa_NT_CasRx", "RNP")
cells_rnp_nbn_nt <- get_cells(seu, "NBN_CRISPRa_NT_CasRx", "RNP")
cells_rnp_nbn_53bp1 <- get_cells(seu, "NBN_CRISPRa_53BP1_CasRx", "RNP")
cells_rnp_nbn_ku70 <- get_cells(seu, "NBN_CRISPRa_KU70_CasRx", "RNP")
cells_rnp_nbn_pq <- get_cells(seu, "NBN_CRISPRa_PQ_CasRx", "RNP")

# Also need CTRL cells for CasRx-effect proxy
cells_ctrl_nbn_nt <- get_cells(seu, "NBN_CRISPRa_NT_CasRx", "CTRL")
cells_ctrl_nbn_53bp1 <- get_cells(seu, "NBN_CRISPRa_53BP1_CasRx", "CTRL")
cells_ctrl_nbn_ku70 <- get_cells(seu, "NBN_CRISPRa_KU70_CasRx", "CTRL")
cells_ctrl_nbn_pq <- get_cells(seu, "NBN_CRISPRa_PQ_CasRx", "CTRL")

cat("RNP cell counts:\n")
cat("  NT_NT:", length(cells_rnp_nt_nt), "\n")
cat("  NBN_NT:", length(cells_rnp_nbn_nt), "\n")
cat("  NBN_TP53BP1:", length(cells_rnp_nbn_53bp1), "\n")
cat("  NBN_XRCC6:", length(cells_rnp_nbn_ku70), "\n")
cat("  NBN_PQ:", length(cells_rnp_nbn_pq), "\n")

# For each partner, compute:
# FC_A = log2FC(RNP_NBN_NT / RNP_NT_NT) — CRISPRa effect under stress
# FC_B = log2FC(CTRL_NBN_X / CTRL_NBN_NT) — CasRx effect (proxy, within NBN context)
# FC_AB = log2FC(RNP_NBN_X / RNP_NT_NT) — Combined dual perturbation
# Residual = FC_AB - (FC_A + FC_B)

FC_A <- compute_log2fc(seu, cells_rnp_nbn_nt, cells_rnp_nt_nt)

partner_cells_rnp <- list(
  "TP53BP1" = cells_rnp_nbn_53bp1,
  "XRCC6"   = cells_rnp_nbn_ku70,
  "POLQ"    = cells_rnp_nbn_pq
)
partner_cells_ctrl <- list(
  "TP53BP1" = cells_ctrl_nbn_53bp1,
  "XRCC6"   = cells_ctrl_nbn_ku70,
  "POLQ"    = cells_ctrl_nbn_pq
)

# Filter to variable features for meaningful analysis
var_features <- VariableFeatures(seu)
if (length(var_features) == 0) {
  var_features <- rownames(seu)
}

epistasis_results <- list()

for (partner in partners) {
  cat("\n--- Partner:", partner, "---\n")

  FC_B <- compute_log2fc(seu, partner_cells_ctrl[[partner]], cells_ctrl_nbn_nt)
  FC_AB <- compute_log2fc(seu, partner_cells_rnp[[partner]], cells_rnp_nt_nt)

  # Compute residual
  residual <- FC_AB - (FC_A + FC_B)

  # Combine into data frame (variable features only)
  genes <- intersect(var_features, names(residual))
  epi_df <- data.frame(
    gene = genes,
    FC_A = FC_A[genes],
    FC_B = FC_B[genes],
    FC_AB = FC_AB[genes],
    expected = FC_A[genes] + FC_B[genes],
    residual = residual[genes],
    partner = partner
  )

  # Classify
  threshold <- 0.3
  epi_df$class <- "Additive"
  epi_df$class[epi_df$residual > threshold] <- "Synergistic"
  epi_df$class[epi_df$residual < -threshold] <- "Buffering"

  cat("  Synergistic:", sum(epi_df$class == "Synergistic"), "\n")
  cat("  Buffering:", sum(epi_df$class == "Buffering"), "\n")
  cat("  Additive:", sum(epi_df$class == "Additive"), "\n")

  epistasis_results[[partner]] <- epi_df
}

# Combine and save
epi_all <- do.call(rbind, epistasis_results)
write.csv(epi_all, file.path(out_dir, "ModuleA_epistasis_gene_level.csv"), row.names = FALSE)

# Summary table
epi_summary <- epi_all %>%
  group_by(partner, class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(partner) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

write.csv(epi_summary, file.path(out_dir, "ModuleA_epistasis_summary.csv"), row.names = FALSE)
cat("\nModule A Summary:\n")
print(epi_summary)

# Scatter: Expected vs Observed
for (partner in partners) {
  epi_df <- epistasis_results[[partner]]

  p_scatter <- ggplot(epi_df, aes(x = expected, y = FC_AB, color = class)) +
    geom_point(size = 1.5, alpha = 0.4) +
    scale_color_manual(values = c(
      "Synergistic" = "#E64B35",
      "Buffering" = "#2166AC",
      "Additive" = "grey70"
    )) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_abline(slope = 1, intercept = threshold, linetype = "dotted", color = "#E64B35", alpha = 0.5) +
    geom_abline(slope = 1, intercept = -threshold, linetype = "dotted", color = "#2166AC", alpha = 0.5) +
    coord_equal() +
    theme_classic(base_size = 14) +
    theme(aspect.ratio = 1) +
    labs(
      title = paste0("Module A: Epistasis - ", partner),
      x = "Expected (FC_A + FC_B)", y = "Observed (FC_AB)", color = ""
    )

  ggsave(file.path(fig_dir, paste0("ModuleA_scatter_", partner, ".pdf")),
    p_scatter,
    width = 7, height = 7
  )
}

##############################################################################
# MODULE B: Pathway-level Tau Scores
##############################################################################
cat("\n=== MODULE B: Pathway-level Tau Scores ===\n")

# Load AUCell results if available
aucell_file <- file.path(OUTPUT_DIR, "07_tier5", "aucell_summary.csv")
if (file.exists(aucell_file)) {
  auc_df <- read.csv(aucell_file)

  # Compute tau = Observed_dual - Expected_additive
  # Need mean AUC per condition for RNP cells
  auc_rnp <- auc_df %>% filter(grepl("RNP", condition))

  # Extract condition components
  auc_rnp$tag <- gsub(" \\| RNP", "", auc_rnp$condition)

  # Get reference values
  ref_nt <- auc_rnp %>%
    filter(tag == "NT+NT (Ctrl)") %>%
    select(pathway, mean_auc_ref = mean_auc)
  ref_nbn <- auc_rnp %>%
    filter(tag == "NBN+NT") %>%
    select(pathway, mean_auc_nbn = mean_auc)

  tau_list <- list()
  for (partner in c("NBN+TP53BP1 KD", "NBN+XRCC6 KD", "NBN+POLQ KD")) {
    dual <- auc_rnp %>%
      filter(tag == partner) %>%
      select(pathway, mean_auc_dual = mean_auc)

    tau_df <- merge(dual, ref_nt, by = "pathway")
    tau_df <- merge(tau_df, ref_nbn, by = "pathway")

    # Expected = AUC_NBN (effect of CRISPRa) + (AUC_dual_proxy - AUC_NBN) estimated
    # Simplified: tau = AUC_dual - AUC_NBN (deviation from single perturbation)
    tau_df$tau <- tau_df$mean_auc_dual - tau_df$mean_auc_nbn
    tau_df$partner <- partner

    tau_list[[partner]] <- tau_df
  }

  tau_all <- do.call(rbind, tau_list)
  write.csv(tau_all, file.path(out_dir, "ModuleB_tau_scores.csv"), row.names = FALSE)

  # Heatmap of tau scores
  tau_wide <- tau_all %>%
    select(pathway, partner, tau) %>%
    pivot_wider(names_from = partner, values_from = tau) %>%
    as.data.frame()
  rownames(tau_wide) <- tau_wide$pathway
  tau_wide$pathway <- NULL

  # Clean names
  rownames(tau_wide) <- gsub("HALLMARK_", "", rownames(tau_wide))
  rownames(tau_wide) <- gsub("_", " ", rownames(tau_wide))

  # Filter to pathways with substantial tau
  max_tau <- apply(abs(as.matrix(tau_wide)), 1, max)
  sig_rows <- max_tau > quantile(max_tau, 0.75)

  if (sum(sig_rows) >= 3) {
    pdf(file.path(fig_dir, "ModuleB_tau_heatmap.pdf"),
      width = 8, height = max(6, sum(sig_rows) * 0.4 + 2)
    )
    pheatmap(as.matrix(tau_wide[sig_rows, ]),
      color = colorRampPalette(c("#2166AC", "white", "#E64B35"))(100),
      breaks = seq(-max(abs(tau_wide[sig_rows, ])),
        max(abs(tau_wide[sig_rows, ])),
        length.out = 101
      ),
      cluster_cols = FALSE,
      fontsize_row = 11,
      fontsize_col = 12,
      main = "Module B: Pathway Tau Scores (Dual - NBN_NT)"
    )
    dev.off()
  }

  cat("Module B complete:", nrow(tau_all), "pathway-partner pairs\n")
} else {
  cat("WARNING: AUCell results not found. Run Step 07 first.\n")
}

##############################################################################
# MODULE C: Transcriptomic PCA Manifold
##############################################################################
cat("\n=== MODULE C: Transcriptomic PCA Manifold ===\n")

# Compute pseudobulk profiles for each condition x treatment
all_conditions <- unique(paste(seu$assigned_tag, seu$treatment, sep = "_"))

# Filter to NBN-axis conditions
keep_conds <- paste(rep(conditions, each = 2),
  rep(c("CTRL", "RNP"), length(conditions)),
  sep = "_"
)

expr_counts <- GetAssayData(seu, layer = "counts")

pb_list <- list()
for (cond in keep_conds) {
  parts <- strsplit(cond, "_(?=(CTRL|RNP)$)", perl = TRUE)[[1]]
  if (length(parts) != 2) next
  tag <- parts[1]
  treat <- parts[2]

  cells <- colnames(seu)[seu$assigned_tag == tag & seu$treatment == treat]
  if (length(cells) < 10) {
    cat("  Skipping", cond, "- only", length(cells), "cells\n")
    next
  }

  pb_list[[cond]] <- rowMeans(expr_counts[, cells])
}

if (length(pb_list) >= 3) {
  pb_mat <- do.call(cbind, pb_list)

  # Log-transform pseudobulk
  pb_log <- log1p(pb_mat)

  # PCA on top variable genes
  gene_var <- apply(pb_log, 1, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(2000, sum(gene_var > 0))]
  pb_scaled <- t(scale(t(pb_log[top_genes, ])))

  # Handle NaN from zero-variance genes
  pb_scaled[is.nan(pb_scaled)] <- 0

  pca_res <- prcomp(t(pb_scaled), center = TRUE, scale. = FALSE)

  # Extract PC coordinates
  pc_df <- as.data.frame(pca_res$x[, 1:min(3, ncol(pca_res$x))])
  pc_df$sample <- rownames(pc_df)

  # Parse condition and treatment
  pc_df$treatment <- ifelse(grepl("_RNP$", pc_df$sample), "RNP", "CTRL")
  pc_df$tag <- gsub("_(CTRL|RNP)$", "", pc_df$sample)
  pc_df$tag_short <- unname(condition_short[pc_df$tag])

  # Variance explained
  var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  p_pca <- ggplot(pc_df, aes(x = PC1, y = PC2, color = tag_short, shape = treatment)) +
    geom_point(size = 5) +
    scale_color_brewer(palette = "Set2") +
    theme_classic(base_size = 14) +
    theme(aspect.ratio = 1) +
    labs(
      title = "Module C: Pseudobulk PCA Manifold",
      subtitle = paste0(ncol(pb_mat), " pseudobulk profiles, ", length(top_genes), " genes"),
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)"),
      color = "Condition", shape = "Treatment"
    )

  ggsave(file.path(fig_dir, "ModuleC_pca_manifold.pdf"), p_pca, width = 7, height = 7)

  # Correlation matrix
  cor_mat <- cor(pb_log[top_genes, ])
  colnames(cor_mat) <- gsub("_CRISPRa|_CasRx", "", colnames(cor_mat))
  rownames(cor_mat) <- gsub("_CRISPRa|_CasRx", "", rownames(cor_mat))

  pdf(file.path(fig_dir, "ModuleC_correlation_heatmap.pdf"), width = 8, height = 8)
  pheatmap(cor_mat,
    color = colorRampPalette(c("#2166AC", "white", "#E64B35"))(100),
    main = "Module C: Pseudobulk Correlation Matrix",
    fontsize = 11
  )
  dev.off()

  cat("Module C complete:", ncol(pb_mat), "pseudobulk profiles\n")
} else {
  cat("WARNING: Fewer than 3 pseudobulk profiles. Skipping PCA.\n")
}

##############################################################################
# MODULE D: Gatekeeper Hierarchy
##############################################################################
cat("\n=== MODULE D: Gatekeeper Hierarchy ===\n")

# Compute rescue percentages using Tier 3 DE results
# Rescue = partner KD reverses genes that are changed by NBN activation
# Use DE from Tier 3 (partner vs NT within RNP_NBN)

tier3_dir <- file.path(OUTPUT_DIR, "04_tier3")
tier2_file <- file.path(OUTPUT_DIR, "03_tier2", "Tier2_interaction_table.csv")

if (file.exists(tier2_file)) {
  tier2 <- read.csv(tier2_file)

  # Genes with substantial interaction effect
  affected_genes <- tier2$gene[abs(tier2$interaction_score) > 0.25]
  cat("Genes with |interaction_score| > 0.25:", length(affected_genes), "\n")

  rescue_results <- list()
  for (partner in partners) {
    de_file <- file.path(tier3_dir, paste0("Tier3_", partner, "_DE.csv"))
    if (!file.exists(de_file)) next

    de <- read.csv(de_file)

    # Merge with interaction data
    merged <- merge(
      tier2[, c("gene", "interaction_score")],
      de[, c("gene", "avg_log2FC", "p_val_adj")],
      by = "gene"
    )

    # A gene is "rescued" if:
    # 1. It had a strong interaction effect (|interaction_score| > 0.25)
    # 2. Partner KD moves it in the opposite direction
    # (i.e., interaction_score * avg_log2FC < 0)
    merged$interaction_affected <- abs(merged$interaction_score) > 0.25
    merged$opposite_direction <- merged$interaction_score * merged$avg_log2FC < 0
    merged$rescued <- merged$interaction_affected & merged$opposite_direction

    n_affected <- sum(merged$interaction_affected)
    n_rescued <- sum(merged$rescued, na.rm = TRUE)
    pct_rescued <- round(100 * n_rescued / max(n_affected, 1), 1)

    rescue_results[[partner]] <- data.frame(
      partner = partner,
      n_total_genes = nrow(merged),
      n_interaction_affected = n_affected,
      n_rescued = n_rescued,
      pct_rescued = pct_rescued
    )

    cat(partner, ": ", n_rescued, "/", n_affected,
      " genes rescued (", pct_rescued, "%)\n",
      sep = ""
    )
  }

  rescue_df <- do.call(rbind, rescue_results)
  write.csv(rescue_df, file.path(out_dir, "ModuleD_rescue_summary.csv"), row.names = FALSE)

  # Bar chart of rescue percentages
  p_rescue <- ggplot(rescue_df, aes(x = partner, y = pct_rescued, fill = partner)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = paste0(pct_rescued, "%")), vjust = -0.5, size = 5) +
    scale_fill_manual(values = c("TP53BP1" = "#E64B35", "XRCC6" = "#2166AC", "POLQ" = "#00A087")) +
    theme_classic(base_size = 14) +
    labs(
      title = "Module D: Gene-Level Rescue by Partner KD",
      subtitle = paste0("Genes with |interaction_score| > 0.25: ", length(affected_genes)),
      x = "Partner Knocked Down", y = "% Genes Rescued", fill = ""
    ) +
    ylim(0, max(rescue_df$pct_rescued) * 1.2)

  ggsave(file.path(fig_dir, "ModuleD_rescue_barplot.pdf"), p_rescue, width = 6, height = 6)
} else {
  cat("WARNING: Tier 2 results not found. Run Step 03 first.\n")
}

cat("\nStep 08 complete. Tables in:", out_dir, "| Figures in:", fig_dir, "\n")

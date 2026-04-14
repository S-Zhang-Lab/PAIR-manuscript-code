##############################################################################
# Step 03: Tier 2 - Stress Capacity Interaction
# PAIR-Perturb-Seq Analysis Pipeline
#
# Input:  data/seu_qc.qs
# Output: res/03_tier2/ (DE tables, interaction table, scatter plot)
#
# Question: Does NBN activation alter the cell's transcriptional response
# to exogenous DNA damage (RNP stress)?
#
# Approach:
#   Vector A (WT stress response):  RNP_NT_NT   vs CTRL_NT_NT
#   Vector B (NBN stress response): RNP_NBN_NT  vs CTRL_NBN_NT
#   Interaction = log2FC_NBN - log2FC_WT
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
out_dir <- file.path(OUTPUT_DIR, "03_tier2")
fig_dir <- file.path(FIGURES_DIR, "03_tier2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading seu_qc.qs...\n")
seu <- qread(file.path(base_dir, "data", "seu_qc.qs"))

# --- Define subsets ---
Idents(seu) <- "assigned_tag"

# Subset to NT_NT and NBN_NT cells
seu_nt <- subset(seu, assigned_tag == "NT_CRISPRa_NT_CasRx")
seu_nbn <- subset(seu, assigned_tag == "NBN_CRISPRa_NT_CasRx")

# Cell counts
n_nt_ctrl <- sum(seu_nt$treatment == "CTRL")
n_nt_rnp <- sum(seu_nt$treatment == "RNP")
n_nbn_ctrl <- sum(seu_nbn$treatment == "CTRL")
n_nbn_rnp <- sum(seu_nbn$treatment == "RNP")

cat(
  "Vector A (WT stress): RNP_NT_NT (n=", n_nt_rnp,
  ") vs CTRL_NT_NT (n=", n_nt_ctrl, ")\n"
)
cat(
  "Vector B (NBN stress): RNP_NBN_NT (n=", n_nbn_rnp,
  ") vs CTRL_NBN_NT (n=", n_nbn_ctrl, ")\n"
)

# --- Run DE for Vector A: WT stress response ---
Idents(seu_nt) <- "treatment"
cat("\nRunning FindMarkers: Vector A (WT stress response)...\n")
de_wt <- FindMarkers(
  seu_nt,
  ident.1 = "RNP",
  ident.2 = "CTRL",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0,
  verbose = FALSE
)
de_wt$gene <- rownames(de_wt)
cat("Vector A: ", nrow(de_wt), " genes tested\n")

# --- Run DE for Vector B: NBN stress response ---
Idents(seu_nbn) <- "treatment"
cat("Running FindMarkers: Vector B (NBN stress response)...\n")
de_nbn <- FindMarkers(
  seu_nbn,
  ident.1 = "RNP",
  ident.2 = "CTRL",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0,
  verbose = FALSE
)
de_nbn$gene <- rownames(de_nbn)
cat("Vector B: ", nrow(de_nbn), " genes tested\n")

# Save individual DE results
write.csv(de_wt, file.path(out_dir, "Tier2_VectorA_WT_stress_DE.csv"), row.names = FALSE)
write.csv(de_nbn, file.path(out_dir, "Tier2_VectorB_NBN_stress_DE.csv"), row.names = FALSE)

# --- Merge with full outer join ---
# IMPORTANT: Use full outer join to keep genes detected in only one comparison
merged <- merge(
  de_wt[, c("gene", "avg_log2FC", "p_val", "p_val_adj")],
  de_nbn[, c("gene", "avg_log2FC", "p_val", "p_val_adj")],
  by = "gene",
  all = TRUE, # Full outer join
  suffixes = c("_WT", "_NBN")
)

# Replace NA log2FC with 0 (gene not detected = no change)
merged$avg_log2FC_WT[is.na(merged$avg_log2FC_WT)] <- 0
merged$avg_log2FC_NBN[is.na(merged$avg_log2FC_NBN)] <- 0

# Compute interaction score
merged$interaction_score <- merged$avg_log2FC_NBN - merged$avg_log2FC_WT

# Classify interaction
merged$class <- "Similar Response"
merged$class[merged$interaction_score > 0.5] <- "Hyper-responsive in NBN"
merged$class[merged$interaction_score < -0.5] <- "Fails to respond in NBN"

merged <- merged %>% arrange(desc(abs(interaction_score)))

cat("\nTotal genes in merged table:", nrow(merged), "\n")
cat("Hyper-responsive:", sum(merged$class == "Hyper-responsive in NBN"), "\n")
cat("Fails to respond:", sum(merged$class == "Fails to respond in NBN"), "\n")
cat("Similar:", sum(merged$class == "Similar Response"), "\n")

# Compute p_plot column before any plotting
merged$min_p <- pmin(merged$p_val_WT, merged$p_val_NBN, na.rm = TRUE)
min_nonzero_p <- min(merged$min_p[merged$min_p > 0], na.rm = TRUE)
merged$p_plot <- ifelse(is.na(merged$min_p) | merged$min_p == 0,
  min_nonzero_p * 0.1, merged$min_p
)

write.csv(merged, file.path(out_dir, "Tier2_interaction_table.csv"), row.names = FALSE)

# --- Scatter Plot: WT vs NBN stress response ---
# Label top interaction genes
top_hyper <- merged %>%
  filter(class == "Hyper-responsive in NBN") %>%
  top_n(10, wt = abs(interaction_score))
top_fail <- merged %>%
  filter(class == "Fails to respond in NBN") %>%
  top_n(10, wt = abs(interaction_score))
label_genes <- rbind(top_hyper, top_fail)

p_scatter <- ggplot(merged, aes(x = avg_log2FC_WT, y = avg_log2FC_NBN, color = class)) +
  geom_point(size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c(
    "Hyper-responsive in NBN" = "#E64B35",
    "Fails to respond in NBN" = "#4DBBD5",
    "Similar Response" = "grey70"
  )) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_abline(slope = 1, intercept = 0.5, linetype = "dotted", color = "#E64B35", alpha = 0.5) +
  geom_abline(slope = 1, intercept = -0.5, linetype = "dotted", color = "#4DBBD5", alpha = 0.5) +
  geom_text_repel(
    data = label_genes, aes(label = gene),
    size = 3.5, max.overlaps = 20, color = "black"
  ) +
  theme_classic(base_size = 14) +
  coord_equal() +
  labs(
    title = "Tier 2: Stress Response Interaction",
    subtitle = paste0(
      "WT: ", n_nt_rnp, " vs ", n_nt_ctrl,
      " | NBN: ", n_nbn_rnp, " vs ", n_nbn_ctrl, " cells"
    ),
    x = "Log2FC (WT Stress Response: RNP vs CTRL in NT_NT)",
    y = "Log2FC (NBN Stress Response: RNP vs CTRL in NBN_NT)",
    color = ""
  )

ggsave(file.path(fig_dir, "Tier2_interaction_scatter.pdf"), p_scatter, width = 7, height = 7)

# --- Volcano-style plot of interaction scores ---
p_volcano <- ggplot(merged, aes(x = interaction_score, y = -log10(p_plot), color = class)) +
  geom_point(size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c(
    "Hyper-responsive in NBN" = "#E64B35",
    "Fails to respond in NBN" = "#4DBBD5",
    "Similar Response" = "grey70"
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = label_genes, aes(label = gene),
    size = 3.5, max.overlaps = 20, color = "black"
  ) +
  theme_classic(base_size = 14) +
  theme(aspect.ratio = 1) +
  labs(
    title = "Tier 2: Interaction Score Distribution",
    x = "Interaction Score (log2FC_NBN - log2FC_WT)",
    y = "-Log10(min P-value)", color = ""
  )

ggsave(file.path(fig_dir, "Tier2_interaction_volcano.pdf"), p_volcano, width = 7, height = 7)

# --- Summary ---
cat("\n=== Tier 2 Summary ===\n")
cat(
  "Vector A (WT stress): ", nrow(de_wt), " genes |",
  sum(de_wt$p_val_adj < 0.05), "significant (adj.P < 0.05)\n"
)
cat(
  "Vector B (NBN stress):", nrow(de_nbn), "genes |",
  sum(de_nbn$p_val_adj < 0.05), "significant (adj.P < 0.05)\n"
)
cat("Merged (full outer join):", nrow(merged), "genes\n")
cat("Hyper-responsive in NBN (score > 0.5):", sum(merged$class == "Hyper-responsive in NBN"), "\n")
cat("Fails to respond in NBN (score < -0.5):", sum(merged$class == "Fails to respond in NBN"), "\n")

cat("\nStep 03 complete.\n")
cat("Tables written to:", out_dir, "\n")
cat("Figures written to:", fig_dir, "\n")

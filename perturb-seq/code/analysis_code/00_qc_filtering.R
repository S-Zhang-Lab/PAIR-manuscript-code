##############################################################################
# Step 00: QC Filtering
# PAIR-Perturb-Seq Analysis Pipeline
#
# Input:  data/seu_prep.qs (unfiltered Seurat object, 38617 genes x 2552 cells)
# Output: data/seu_qc.qs (QC-filtered, re-processed Seurat object)
#         res/00_qc/ (QC plots and summary tables)
#
# Method: MAD-based adaptive thresholds (scverse best practice)
#   - log(nCount_RNA): median +/- 3*MAD
#   - log(nFeature_RNA): median +/- 3*MAD
#   - percent.mt: median + 3*MAD (upper bound only)
#
# No doublet detection: HTO demux from Cell Ranger multi already filters.
# No batch correction: all samples from same 10x lane.
##############################################################################

library(Seurat)
library(qs)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(42)

# --- Source config.R (sets DATA_ROOT, DATA_DIR, RES_DIR; see code/config.R) -
if (file.exists("code/config.R")) {
  source("code/config.R")                 # run from repo/
} else if (file.exists("../config.R")) {
  source("../config.R")                   # run from repo/code/analysis_code/
} else {
  stop("config.R not found. Run from repo/ or repo/code/analysis_code/.")
}
# ---------------------------------------------------------------------------

# --- Paths ---
base_dir  <- DATA_ROOT
data_dir  <- file.path(base_dir, "data")
out_dir <- file.path(OUTPUT_DIR, "00_qc")
fig_dir <- file.path(FIGURES_DIR, "00_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load ---
cat("Loading seu_prep.qs...\n")
seu <- qread(SEU_PREP)
cat("Input:", ncol(seu), "cells,", nrow(seu), "genes\n")

# --- Pre-QC summary ---
pre_qc_summary <- data.frame(
  metric = c("total_cells", "total_genes",
             "median_nCount", "median_nFeature", "median_pct_mt",
             "cells_nFeature_lt500", "cells_pct_mt_gt10", "cells_nCount_gt100k"),
  value = c(ncol(seu), nrow(seu),
            median(seu$nCount_RNA), median(seu$nFeature_RNA), round(median(seu$percent.mt), 3),
            sum(seu$nFeature_RNA < 500), sum(seu$percent.mt > 10), sum(seu$nCount_RNA > 100000))
)

# --- Compute MAD-based thresholds ---
log_counts   <- log1p(seu$nCount_RNA)
log_features <- log1p(seu$nFeature_RNA)
pct_mt       <- seu$percent.mt

# MAD thresholds (3 MAD from median)
count_median   <- median(log_counts)
count_mad      <- mad(log_counts)
feature_median <- median(log_features)
feature_mad    <- mad(log_features)
mt_median      <- median(pct_mt)
mt_mad         <- mad(pct_mt)

count_lo   <- exp(count_median - 3 * count_mad) - 1
count_hi   <- exp(count_median + 3 * count_mad) - 1
feature_lo <- exp(feature_median - 3 * feature_mad) - 1
feature_hi <- exp(feature_median + 3 * feature_mad) - 1
mt_hi      <- mt_median + 3 * mt_mad

cat("\n=== MAD-based QC Thresholds ===\n")
cat("nCount_RNA:   [", round(count_lo), ",", round(count_hi), "]\n")
cat("nFeature_RNA: [", round(feature_lo), ",", round(feature_hi), "]\n")
cat("percent.mt:   <=", round(mt_hi, 2), "\n")

# --- Apply filters ---
keep <- seu$nCount_RNA >= count_lo &
        seu$nCount_RNA <= count_hi &
        seu$nFeature_RNA >= feature_lo &
        seu$nFeature_RNA <= feature_hi &
        seu$percent.mt <= mt_hi

cat("\nCells passing QC:", sum(keep), "/", length(keep),
    "(removed:", sum(!keep), "=", round(100 * sum(!keep) / length(keep), 1), "%)\n")

# Tag the filter reason for removed cells
seu$qc_pass <- keep
filter_reasons <- rep("pass", ncol(seu))
filter_reasons[seu$nCount_RNA < count_lo]    <- "low_counts"
filter_reasons[seu$nCount_RNA > count_hi]    <- "high_counts"
filter_reasons[seu$nFeature_RNA < feature_lo] <- "low_features"
filter_reasons[seu$nFeature_RNA > feature_hi] <- "high_features"
filter_reasons[seu$percent.mt > mt_hi]        <- "high_mito"
seu$filter_reason <- filter_reasons

# --- QC Plots: Before filtering ---

# Violin plots
p_vln_before <- VlnPlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                         group.by = "filter_reason", pt.size = 0.1, ncol = 3) +
  plot_annotation(title = "QC Metrics by Filter Status (Before Filtering)")

ggsave(file.path(fig_dir, "qc_violin_by_filter_status.pdf"), p_vln_before,
       width = 14, height = 5)

# Scatter: nCount vs nFeature, colored by mito
df_scatter <- data.frame(
  nCount = seu$nCount_RNA,
  nFeature = seu$nFeature_RNA,
  pct_mt = seu$percent.mt,
  pass = seu$qc_pass
)

p_scatter1 <- ggplot(df_scatter, aes(x = nCount, y = nFeature, color = pct_mt)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(name = "% Mito") +
  geom_vline(xintercept = c(count_lo, count_hi), linetype = "dashed", color = "red") +
  geom_hline(yintercept = c(feature_lo, feature_hi), linetype = "dashed", color = "red") +
  scale_x_log10() +
  theme_classic() +
  labs(title = "QC Scatter: nCount vs nFeature", subtitle = "Red lines = MAD thresholds")

p_scatter2 <- ggplot(df_scatter, aes(x = nCount, y = pct_mt, color = pass)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "grey30", "FALSE" = "red"), name = "Pass QC") +
  geom_hline(yintercept = mt_hi, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(count_lo, count_hi), linetype = "dashed", color = "red") +
  scale_x_log10() +
  theme_classic() +
  labs(title = "QC Scatter: nCount vs % Mito")

p_scatter <- p_scatter1 / p_scatter2
ggsave(file.path(fig_dir, "qc_scatter_thresholds.pdf"), p_scatter, width = 8, height = 10)

# --- Cell counts per tag x treatment before/after ---
tag_treat_before <- as.data.frame(table(seu$assigned_tag, seu$treatment))
colnames(tag_treat_before) <- c("Tag", "Treatment", "Before")

seu_filt <- subset(seu, cells = colnames(seu)[keep])
tag_treat_after <- as.data.frame(table(seu_filt$assigned_tag, seu_filt$treatment))
colnames(tag_treat_after) <- c("Tag", "Treatment", "After")

cell_counts <- merge(tag_treat_before, tag_treat_after, by = c("Tag", "Treatment"), all = TRUE)
cell_counts$After[is.na(cell_counts$After)] <- 0
cell_counts$Removed <- cell_counts$Before - cell_counts$After
cell_counts$Pct_Removed <- round(100 * cell_counts$Removed / cell_counts$Before, 1)

# Add totals per tag
tag_totals <- cell_counts %>%
  group_by(Tag) %>%
  summarise(Total_Before = sum(Before), Total_After = sum(After), .groups = "drop") %>%
  mutate(Flag = ifelse(Total_After < 30, "LOW_POWER", ""))

write.csv(cell_counts, file.path(out_dir, "cell_counts_by_tag_treatment.csv"), row.names = FALSE)
write.csv(tag_totals, file.path(out_dir, "cell_counts_by_tag_totals.csv"), row.names = FALSE)

cat("\n=== Cell Counts by Tag (Before -> After) ===\n")
print(tag_totals)

# Flag low-power groups
low_power <- tag_totals$Tag[tag_totals$Total_After < 30]
if (length(low_power) > 0) {
  cat("\nWARNING: Tags with < 30 cells after QC:", paste(low_power, collapse = ", "), "\n")
}

# --- Re-process filtered object ---
cat("\nRe-processing filtered object...\n")

# Remove old reductions and re-process
seu_filt <- DietSeurat(seu_filt, layers = c("counts", "data"), dimreducs = NULL, graphs = NULL)

# Re-normalize
seu_filt <- NormalizeData(seu_filt, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features (2000)
seu_filt <- FindVariableFeatures(seu_filt, selection.method = "vst", nfeatures = 2000)

# Scale data (variable features only — more efficient than scaling all 38K genes)
seu_filt <- ScaleData(seu_filt)

# PCA
seu_filt <- RunPCA(seu_filt, npcs = 50, verbose = FALSE)

# Elbow plot
p_elbow <- ElbowPlot(seu_filt, ndims = 50) + ggtitle("Elbow Plot (Post-QC)")
ggsave(file.path(fig_dir, "elbow_plot.pdf"), p_elbow, width = 6, height = 4)

# Clustering and UMAP (use 20 PCs as in original)
seu_filt <- FindNeighbors(seu_filt, dims = 1:20)
seu_filt <- FindClusters(seu_filt, resolution = 0.4)
seu_filt <- RunUMAP(seu_filt, dims = 1:20)

# --- Post-QC plots ---

# UMAP by tag
p_umap_tag <- DimPlot(seu_filt, group.by = "assigned_tag", label = FALSE, pt.size = 0.5) +
  ggtitle("UMAP: Assigned Tag (Post-QC)") +
  theme(legend.text = element_text(size = 7))

# UMAP by treatment
p_umap_treat <- DimPlot(seu_filt, group.by = "treatment", pt.size = 0.5) +
  ggtitle("UMAP: Treatment (Post-QC)")

# UMAP by cluster
p_umap_clust <- DimPlot(seu_filt, group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP: Clusters (Post-QC)")

# UMAP by cell cycle
p_umap_phase <- DimPlot(seu_filt, group.by = "Phase", pt.size = 0.5) +
  ggtitle("UMAP: Cell Cycle Phase")

p_umap_all <- (p_umap_tag | p_umap_treat) / (p_umap_clust | p_umap_phase)
ggsave(file.path(fig_dir, "umap_post_qc.pdf"), p_umap_all, width = 16, height = 12)

# Post-QC violin plots
p_vln_after <- VlnPlot(seu_filt, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                        group.by = "seurat_clusters", pt.size = 0.1, ncol = 3) +
  plot_annotation(title = "QC Metrics by Cluster (Post-QC)")
ggsave(file.path(fig_dir, "qc_violin_post_qc.pdf"), p_vln_after, width = 14, height = 5)

# --- Save ---
cat("Saving filtered object...\n")
qsave(seu_filt, SEU_QC)

# --- Summary report ---
post_qc_summary <- data.frame(
  metric = c("total_cells", "total_genes",
             "median_nCount", "median_nFeature", "median_pct_mt",
             "n_clusters", "n_variable_features"),
  value = c(ncol(seu_filt), nrow(seu_filt),
            median(seu_filt$nCount_RNA), median(seu_filt$nFeature_RNA),
            round(median(seu_filt$percent.mt), 3),
            length(levels(seu_filt$seurat_clusters)),
            length(VariableFeatures(seu_filt)))
)

summary_text <- paste0(
  "=== QC Filtering Summary ===\n\n",
  "MAD Thresholds:\n",
  "  nCount_RNA:   [", round(count_lo), ", ", round(count_hi), "]\n",
  "  nFeature_RNA: [", round(feature_lo), ", ", round(feature_hi), "]\n",
  "  percent.mt:   <= ", round(mt_hi, 2), "\n\n",
  "Before QC: ", pre_qc_summary$value[1], " cells\n",
  "After QC:  ", ncol(seu_filt), " cells\n",
  "Removed:   ", ncol(seu) - ncol(seu_filt), " cells (",
  round(100 * (ncol(seu) - ncol(seu_filt)) / ncol(seu), 1), "%)\n\n",
  "Post-QC Medians:\n",
  "  nCount_RNA:   ", post_qc_summary$value[3], "\n",
  "  nFeature_RNA: ", post_qc_summary$value[4], "\n",
  "  percent.mt:   ", post_qc_summary$value[5], "\n\n",
  "Clusters: ", post_qc_summary$value[6], "\n",
  "Variable features: ", post_qc_summary$value[7], "\n\n",
  "Low-power tags (<30 cells): ",
  ifelse(length(low_power) > 0, paste(low_power, collapse = ", "), "None"), "\n"
)

writeLines(summary_text, file.path(out_dir, "qc_summary.txt"))
cat(summary_text)
cat("\nDone. Saved to:", SEU_QC, "\n")
cat("Tables written to:", out_dir, "\n")
cat("Figures written to:", fig_dir, "\n")

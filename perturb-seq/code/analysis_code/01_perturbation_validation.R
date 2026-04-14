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
seu <- qread(file.path(base_dir, "data", "seu_qc.qs"))
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

cat("\nStep 01 complete.\n")
cat("Tables written to:", out_dir, "\n")
cat("Figures written to:", fig_dir, "\n")

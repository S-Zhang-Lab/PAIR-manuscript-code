rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(tidyr)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")

NBN_subset <- c(
  "NBN_CRISPRa_NT_CasRx",
  "NBN_CRISPRa_PQ_CasRx",
  "NBN_CRISPRa_53BP1_CasRx",
  "NBN_CRISPRa_KU70_CasRx"
)

# Keep NBN subsets
seu <- subset(scRNA_seq, subset = assigned_tag %in% NBN_subset) 
meta_data <- seu@meta.data

### PAIR UMAP ### 
# Set identities
Idents(seu) <- "assigned_tag"

tag_levels <- NBN_subset

seu$assigned_tag <- factor(seu$assigned_tag, levels = tag_levels)
Idents(seu) <- "assigned_tag"

# ---- Soft pastel publication palette ----
tag_cols <- c(
  
  # NBN (slightly stronger pastel)
  "NBN_CRISPRa_NT_CasRx"     = "#5DA5DA",  # deeper blue
  "NBN_CRISPRa_PQ_CasRx"     = "#F28E2B",  # muted orange
  "NBN_CRISPRa_53BP1_CasRx"  = "#4EBA6F",  # deeper green
  "NBN_CRISPRa_KU70_CasRx"   = "#E15759"  # richer coral
)

# ---- UMAP ----
p_umap <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "assigned_tag",
  cols = tag_cols,
  pt.size = 0.7
) +
  theme_classic(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 10),
    axis.text    = element_text(size = 10),
    axis.title   = element_text(size = 12),
    title = element_blank()
  )

ggsave("./res/2026_0225/UMAP_assigned_tag_NBNsubset.pdf",
       p_umap,
       width = 7,
       height = 5,
       useDingbats = FALSE)



### treatment UMAP ### 
# Set identity
Idents(seu) <- "treatment"

# Make sure factor levels are ordered
seu$treatment <- factor(seu$treatment, levels = c("CTRL", "RNP"))

# Soft publication colors
treat_cols <- c(
  "CTRL" = "#7FB3D5",   # soft powder blue
  "RNP"  = "#E07A8D"    # muted rose
)

# UMAP
p_treat <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "treatment",
  cols = treat_cols,
  pt.size = 0.7
) +
  theme_classic(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 11),
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 10),
    title = element_blank()
  )

ggsave("./res/2026_0225/UMAP_treatment_NBNsubset.pdf",
       p_treat,
       width = 5.5,
       height = 5,
       useDingbats = FALSE)


### quantification of number of cells ###
meta <- seu@meta.data

df_tag_treat <- meta %>%
  count(assigned_tag, treatment, name = "n") %>%
  group_by(assigned_tag) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# soft-ish, darker treatment colors
treat_cols <- c("CTRL" = "#7FB3D5", "RNP" = "#E07A8D")

p1 <- ggplot(df_tag_treat, aes(x = assigned_tag, y = n, fill = treatment)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = treat_cols) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        legend.title = element_blank()) +
  labs(x = NULL, y = "Number of cells")

ggsave("./res/2026_0225/quantification_num_cells_NBNsubset.pdf",
       p1,
       width = 3,
       height = 5,
       useDingbats = FALSE)


### Gene expression sanity check ###
gene_list <- c("NBN", 
               "XRCC6", # "KU70", 
               "TP53BP1", # "53BP1", 
               "POLQ" # "PQ",
              )
seu$tag_treat <- paste(seu$assigned_tag, seu$treatment, sep = "_")

Idents(seu) <- "tag_treat"


p_dot <- DotPlot(
  seu,
  features = gene_list,
  group.by = "tag_treat",
  scale = FALSE,        # ← THIS removes z-scoring
  dot.scale = 6
) +
  scale_color_gradient(
    low = "#F7FBFF",
    high = "#E15759",
    limits = c(0, 3),      # cap at 2
    oob = scales::squish   # values >2 shown as max color
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )

ggsave("./res/2026_0225/sanity_gene_exp_dot_NBNsubset.pdf",
       p_dot,
       width = 6,
       height = 4,
       useDingbats = FALSE)

### cell expression heatmap ###
# -----------------------------
# 1. Genes of interest
# -----------------------------
gene_list <- c("NBN", "XRCC6", "TP53BP1", "POLQ")

# Check genes exist
gene_list <- gene_list[gene_list %in% rownames(seu)]
if (length(gene_list) == 0) stop("None of the genes found in object.")

# -----------------------------
# 2. Sample cells per assigned_tag
#    (recommended for clean figure)
# -----------------------------
set.seed(123)
n_per_group <- 200

cells_sampled <- meta %>%
  group_by(assigned_tag) %>%
  arrange(runif(n()), .by_group = TRUE) %>%  # shuffle within each group
  slice_head(n = n_per_group) %>%           # take up to n_per_group (if fewer, takes all)
  ungroup() %>%
  arrange(assigned_tag) %>%
  pull(cell)

# -----------------------------
# 3. Extract expression matrix
#    Using RNA log-normalized data
# -----------------------------
mat <- GetAssayData(
  seu,
  assay = "RNA",
  slot = "data"
)[gene_list, cells_sampled, drop = FALSE]

mat <- as.matrix(mat)

# -----------------------------
# 4. Cap expression at 2
# -----------------------------
mat[mat < 0] <- 0
mat[mat > 3] <- 3

# -----------------------------
# 5. Column annotation
# -----------------------------
anno_col <- meta[match(cells_sampled, meta$cell), "assigned_tag", drop = FALSE]
rownames(anno_col) <- cells_sampled

# --- your strategy colors (named) ---
tag_cols_base <- c(
  "NBN_CRISPRa_NT_CasRx"    = "#5DA5DA",  # deeper blue
  "NBN_CRISPRa_PQ_CasRx"    = "#F28E2B",  # muted orange
  "NBN_CRISPRa_53BP1_CasRx" = "#4EBA6F",  # deeper green
  "NBN_CRISPRa_KU70_CasRx"  = "#E15759"   # richer coral
)

# Levels actually present (and in your annotation)
tag_levels <- unique(anno_col$assigned_tag)

# Add colors for any tags not in tag_cols_base (to avoid missing-level errors)
missing_tags <- setdiff(tag_levels, names(tag_cols_base))
if (length(missing_tags) > 0) {
  extra_cols <- setNames(
    colorRampPalette(c("#9C7ED6", "#76B7B2", "#EDC948", "#B07AA1", "#59A14F"))(length(missing_tags)),
    missing_tags
  )
  tag_cols <- c(tag_cols_base, extra_cols)
} else {
  tag_cols <- tag_cols_base
}

# Keep only tags present (and in the same order, optional)
tag_cols <- tag_cols[tag_levels]

anno_colors <- list(
  assigned_tag = tag_cols
)



# -----------------------------
# 6. Save heatmap to PDF
# -----------------------------
pdf("./res/2026_0225/Heatmap_single_cell_assigned_tag_NBNsubset.pdf",
    width = 8,
    height = 3)

pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  color = colorRampPalette(c("#F7FBFF", "#FF9F9B", "#E15759"))(100),
  fontsize_row = 12
)

dev.off()


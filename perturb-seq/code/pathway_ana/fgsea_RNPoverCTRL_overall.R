rm(list = ls())
set.seed(7)
library(data.table) 
library(Seurat)
library(dplyr)
library(fgsea)
library(qs)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")
meta_data <- scRNA_seq@meta.data

# Convert loge to log2
scRNA_seq <- log2(exp(scRNA_seq[["RNA"]]$data)) %>% as.matrix()

all(colnames(scRNA_seq) == rownames(meta_data)) # TRUE

pwy <- c("c2")

RNP_exp <- scRNA_seq[, meta_data$treatment == "RNP"]
CTRL_exp <- scRNA_seq[, meta_data$treatment == "CTRL"]

# Calculate log2 fold changes
RNP_over_CTRL <- sort(rowMeans(RNP_exp) - rowMeans(CTRL_exp), decreasing = TRUE)

signature <- readRDS(paste0("./data/pathways/", pwy, ".rds"))

# Run fGSEA----
fgseaRes_RNP_over_CTRL <- fgsea(signature, RNP_over_CTRL, minSize=15, maxSize=500)

# only keep significant pathways
fwrite(fgseaRes_RNP_over_CTRL, paste0("./res/2025_1203/RNP_over_CTRL_overall.csv"))

pathway_embed <- qread("./data/pathway_embedding_metadata_umap.qs")
rownames(pathway_embed) <- pathway_embed$cell_id

rownames(fgseaRes_RNP_over_CTRL) <- fgseaRes_RNP_over_CTRL$pathway

overlap_pathways <- intersect(rownames(fgseaRes_RNP_over_CTRL), rownames(pathway_embed))

fgseaRes_RNP_over_CTRL <- as.data.frame(fgseaRes_RNP_over_CTRL)
fgseaRes_RNP_over_CTRL_sub <- fgseaRes_RNP_over_CTRL[overlap_pathways, ]

pathway_embed <- pathway_embed[overlap_pathways, ]

plot_tab <- cbind(fgseaRes_RNP_over_CTRL_sub, pathway_embed)

plot_tab_clean <- plot_tab %>%
  as.data.frame() %>%
  dplyr::select(padj, NES, umap_1, umap_2)  # drops any duplicated cols

q <- ggplot(plot_tab_clean, aes(x = umap_1, y = umap_2, color = NES)) +
  geom_point(size = 1, alpha = 0.9) +
  scale_color_gradient2(
    low = "#2166AC",     # blue for negative NES
    mid = "white",       # zero
    high = "#B2182B",    # red for positive NES
    midpoint = 0,
    name = "NES"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
ggsave("./res/2025_1203/fgsea_RNPoverCTRL_overall_umap_NES.png", plot = q, width = 6, height = 5)

# NES larger than 1
plot_tab_clean_filtered <- plot_tab_clean %>%
  filter(abs(NES) > 1)
q <- ggplot(plot_tab_clean_filtered, aes(x = umap_1, y = umap_2, color = NES)) +
  geom_point(size = 1, alpha = 0.9) +
  scale_color_gradient2(
    low = "#2166AC",     # blue for negative NES
    mid = "white",       # zero
    high = "#B2182B",    # red for positive NES
    midpoint = 0,
    name = "NES"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
ggsave("./res/2025_1203/fgsea_RNPoverCTRL_overall_umap_NES_larger1.png", plot = q, width = 6, height = 5)





# Create a categorical variable: significant vs non-significant
plot_tab_clean_sig <- plot_tab_clean %>%
  mutate(
    sig = padj < 0.05,
    NES_for_color = ifelse(sig, NES, NA)   # NA → grey
  )

# Plot
q <- ggplot(plot_tab_clean_sig, aes(x = umap_1, y = umap_2)) +
  # non-significant dots
  geom_point(
    data = subset(plot_tab_clean_sig, !sig),
    color = "grey80",
    size = 0.8,
    alpha = 0.9
  ) +
  # significant dots colored by NES
  geom_point(
    data = subset(plot_tab_clean_sig, sig),
    aes(color = NES_for_color),
    size = 1.0,
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "NES"
  ) +
  theme_classic(base_size = 14) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

# Save
ggsave(
  "./res/2025_1203/fgsea_RNPoverCTRL_overall_umap_NES_sig.png",
  plot = q,
  width  = 6,
  height = 5,
  dpi    = 300
)

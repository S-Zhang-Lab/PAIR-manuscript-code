rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(tidyr)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")

# remove MRE11_CRISPRa_NT_CasRx
seu <- subset(scRNA_seq, subset = assigned_tag != "MRE11_CRISPRa_NT_CasRx")
meta_data <- seu@meta.data

### PAIR UMAP ### 
# Set identities
Idents(seu) <- "assigned_tag"

tag_levels <- c(
  "NBN_CRISPRa_NT_CasRx",
  "NBN_CRISPRa_PQ_CasRx",
  "NBN_CRISPRa_53BP1_CasRx",
  "NBN_CRISPRa_KU70_CasRx",
  "NBN_CRISPRa_DNPKcs_CasRx",
  "NT_CRISPRa_NT_CasRx",
  "NT_CRISPRa_PQ_CasRx",
  "NT_CRISPRa_53BP1_CasRx",
  "NT_CRISPRa_KU70_CasRx",
  "NT_CRISPRa_DNPKcs_CasRx"
)

seu$assigned_tag <- factor(seu$assigned_tag, levels = tag_levels)
Idents(seu) <- "assigned_tag"

# ---- Soft pastel publication palette ----
tag_cols <- c(
  
  # NBN (slightly stronger pastel)
  "NBN_CRISPRa_NT_CasRx"     = "#5DA5DA",  # deeper blue
  "NBN_CRISPRa_PQ_CasRx"     = "#F28E2B",  # muted orange
  "NBN_CRISPRa_53BP1_CasRx"  = "#4EBA6F",  # deeper green
  "NBN_CRISPRa_KU70_CasRx"   = "#E15759",  # richer coral
  "NBN_CRISPRa_DNPKcs_CasRx" = "#9C7ED6",   # deeper lavender
  
  # NT (lighter version of same hues)
  "NT_CRISPRa_NT_CasRx"      = "#A1C9F4",
  "NT_CRISPRa_PQ_CasRx"      = "#FFB482",
  "NT_CRISPRa_53BP1_CasRx"   = "#8DE5A1",
  "NT_CRISPRa_KU70_CasRx"    = "#FF9F9B",
  "NT_CRISPRa_DNPKcs_CasRx"  = "#D0BBFF"
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

ggsave("./res/2026_0225/UMAP_assigned_tag_allPAIR.pdf",
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

ggsave("./res/2026_0225/UMAP_treatment_allPAIR.pdf",
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

ggsave("./res/2026_0225/quantification_num_cells.pdf",
       p1,
       width = 5.5,
       height = 5,
       useDingbats = FALSE)


### Gene expression sanity check ###
gene_list <- c("NBN", 
               "XRCC6", # "KU70", 
               "TP53BP1", # "53BP1", 
               "PRKDC", # "DNPKcs", 
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

ggsave("./res/2026_0225/sanity_gene_exp_dot_allPAIR.pdf",
       p_dot,
       width = 7,
       height = 5,
       useDingbats = FALSE)


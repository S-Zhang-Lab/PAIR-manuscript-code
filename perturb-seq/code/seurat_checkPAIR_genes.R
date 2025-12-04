rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(tidyr)
library(readxl)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")
meta_data <- scRNA_seq@meta.data

# ---- Input: your gene list ----
gene_list <- c("NBN", 
               "XRCC6", # "KU70", 
               "TP53BP1", # "53BP1", 
               "PRKDC", # "DNPKcs", 
               "PARP1", # "PQ", 
               "MRE11",
               "TP53") 

p <- VlnPlot(
  scRNA_seq,
  features = gene_list,
  group.by = "assigned_tag",
  pt.size = 0.5,
  combine = FALSE
)

# apply the theme to each violin plot
p <- lapply(
  p,
  function(x) x + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                        legend.position = "none")
)

# display
p <- patchwork::wrap_plots(p)

ggsave("./res/2025_1203/violin_plot_DNA_repair_genes_by_celltype_combined.png", plot = p, width = 20, height = 18)

### RNP ###
scRNA_seq_RNP <- subset(scRNA_seq, subset = treatment == "RNP")
p <- VlnPlot(
  scRNA_seq_RNP,
  features = gene_list,
  group.by = "assigned_tag",
  pt.size = 0.5,
  combine = FALSE
)

# apply the theme to each violin plot
p <- lapply(
  p,
  function(x) x + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                        legend.position = "none")
)

# display
p <- patchwork::wrap_plots(p)

ggsave("./res/2025_1203/violin_plot_DNA_repair_genes_by_celltype_RNP.png", plot = p, width = 20, height = 18)


### CTRL ###
scRNA_seq_CTRL <- subset(scRNA_seq, subset = treatment == "CTRL")
p <- VlnPlot(
  scRNA_seq_CTRL,
  features = gene_list,
  group.by = "assigned_tag",
  pt.size = 0.5,
  combine = FALSE
)

# apply the theme to each violin plot
p <- lapply(
  p,
  function(x) x + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                        legend.position = "none")
)

# display
p <- patchwork::wrap_plots(p)

ggsave("./res/2025_1203/violin_plot_DNA_repair_genes_by_celltype_CTRL.png", plot = p, width = 20, height = 18)

rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(tidyr)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")

# 1) Get canonical S and G2M gene sets 
cc.genes <- Seurat::cc.genes.updated.2019
s.genes  <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# 2) Score cell cycle (adds S.Score, G2M.Score, Phase to metadata)
# Make sure default assay is RNA (or the assay that contains your expression)
DefaultAssay(scRNA_seq) <- "RNA"

scRNA_seq <- CellCycleScoring(
  object = scRNA_seq,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)

# qsave(scRNA_seq, "./data/seu_prep.qs")  # save updated Seurat object

p <- DimPlot(scRNA_seq, reduction = "umap", group.by = "Phase", label = FALSE)
ggsave("./res/2025_1217/umap_cell_cycle_phase.png", plot = p, width = 6, height = 5)

p <- DimPlot(scRNA_seq, reduction = "umap", group.by = "treatment", label = FALSE)
ggsave("./res/2025_1217/umap_treatment.png", plot = p, width = 6, height = 5)


seu_sub <- subset(scRNA_seq, subset = assigned_tag == "NBN_CRISPRa_NT_CasRx")
p <- DimPlot(seu_sub, reduction = "umap", group.by = "Phase", label = FALSE)
ggsave("./res/2025_1217/umap_cell_cycle_phase_NBN_CRISPRa_NT_CasRx.png", plot = p, width = 6, height = 5)
p <- DimPlot(seu_sub, reduction = "umap", group.by = "treatment", label = FALSE)
ggsave("./res/2025_1217/umap_treatment_NBN_CRISPRa_NT_CasRx.png", plot = p, width = 6, height = 5)


meta_data_sub <- seu_sub@meta.data
meta_data_sub$treat_phase <- paste0(meta_data_sub$treatment, "_", meta_data_sub$Phase)
seu_sub@meta.data <- meta_data_sub

gene_list <- c("NBN", 
               "XRCC6", # "KU70", 
               "TP53BP1", # "53BP1", 
               "PRKDC", # "DNPKcs", 
               "POLQ", # "PQ", 
               "MRE11",
               "TP53",
               "RPA1", 
               "RPA2", 
               "RPA3",
               "RAD50",
               "BRCA1",
               "ATM",
               "PALB2") 


p <- VlnPlot(
  seu_sub,
  features = gene_list,
  group.by = "treat_phase",
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

ggsave("./res/2025_1217/violin_plot_DNA_repair_genes_NBN_NT.png", 
       plot = p, width = 22, height = 22)
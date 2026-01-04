rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(tidyr)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")
meta_data <- scRNA_seq@meta.data
# Create combined grouping
meta_data$group <- paste(meta_data$treatment, meta_data$assigned_tag, sep = "_")

scRNA_seq@meta.data <- meta_data

# ---- Input: your gene list ----
gene_list <- c(
  "NBN", 
  "XRCC6", # "KU70",
  "TP53BP1", # "53BP1",
  "PRKDC", # "DNPKcs",
  "POLQ" # "PQ",
)

cond <- c("CTRL", "RNP")
pos2_genes <- c("53BP1", "DNPKcs", "KU70", "PQ")
pos2_gene_symbol <- c("TP53BP1", "PRKDC", "XRCC6", "POLQ")

chosen_cond <- "CTRL"
chosen_pos2 <- "53BP1"
chosen_pos2_symbol <- "TP53BP1"

for (i in 1:length(cond)) {
  chosen_cond <- cond[i]
  
  for (j in 1:length(pos2_genes)){
    chosen_pos2 <- pos2_genes[j]
    chosen_pos2_symbol <- pos2_gene_symbol[j]
    
    # Construct the names of groups to plot
    chosen_groups <- c(
      paste0(chosen_cond, "_NT_CRISPRa_NT_CasRx"),
      paste0(chosen_cond, "_NBN_CRISPRa_NT_CasRx"),
      paste0(chosen_cond, "_NT_CRISPRa_", chosen_pos2, "_CasRx"),
      paste0(chosen_cond, "_NBN_CRISPRa_", chosen_pos2, "_CasRx")
    )
    
    Idents(scRNA_seq) <- "group"
    scRNA_seq_sub <- subset(scRNA_seq, idents = chosen_groups)
    
    # Make violin plots
    chosen_genes <- c("NBN", chosen_pos2_symbol)
    
    p <- VlnPlot(
      scRNA_seq_sub,
      features = chosen_genes,
      group.by = "group",
      pt.size = 0.1,
      ncol = length(chosen_genes),
      combine = FALSE
    )
    
    # apply the theme to each violin plot
    p <- lapply(
      p,
      function(x) x + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                            axis.title.x = element_blank(),
                            legend.position = "none")
    )
    
    # display
    p <- patchwork::wrap_plots(p)
    
    save_path <- paste0(
      "./res/2026_0104/",
      chosen_cond,
      "_NBN_",
      chosen_pos2_symbol,
      ".pdf"
    )
    
    ggsave(
      filename = save_path, # output name
      plot = p, # plot object
      width = 6,
      height = 6
    )
    
  }
}




rm(list = ls())
set.seed(7)
library(data.table) 
library(Seurat)
library(dplyr)
library(fgsea)
library(qs)
source("./code/pathway_ana/fcn_plot_pathway_embed.R")
pathway_embed <- qread("./data/pathway_embedding_metadata_umap.qs")

input_dir  <- "./res/2026_0225/pwy_results_RNP"
output_dir_overall <- "./res/2026_0225/pwy_results_RNP/overall"
output_dir_NES_larger1 <- "./res/2026_0225/pwy_results_RNP/NES_larger1"
output_dir_sig <- "./res/2026_0225/pwy_results_RNP/sig"

# Create output folder if it does not exist
if (!dir.exists(output_dir_overall)) {
  dir.create(output_dir_overall)
}

if (!dir.exists(output_dir_NES_larger1)) {
  dir.create(output_dir_NES_larger1)
}

if (!dir.exists(output_dir_sig)) {
  dir.create(output_dir_sig)
}

# List all CSV files in the directory
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# Loop over files
for (f in files) {
  
  # Extract file name (without path)
  fname <- basename(f)
  
  # Read CSV
  df <- read.csv(f, header = TRUE)
  
  # Apply your function
  result <- organize_plot_dat(df, pathway_embed)
  
  # Build output path
  # Example: RNP_over_CTRL_MRE11_CRISPRa_NT_CasRx_c2_processed.csv
  out_name_overall <- sub("\\.csv$", "_overall.png", fname)
  out_path_overall <- file.path(output_dir_overall, out_name_overall)
  
  out_name_NES_larger1 <- sub("\\.csv$", "_NES_larger1.png", fname)
  out_path_NES_larger1 <- file.path(output_dir_NES_larger1, out_name_NES_larger1)
  
  out_name_sig <- sub("\\.csv$", "_sig_padj.png", fname)
  out_path_sig <- file.path(output_dir_sig, out_name_sig)
  
  plot_overall(
    result,
    save_path = out_path_overall,
    width = 6,
    height = 5
  )
  
  plot_NES_larger1(
    result,
    save_path = out_path_NES_larger1,
    width = 6,
    height = 5
  )
  
  plot_sig(
    result,
    save_path = out_path_sig,
    width = 6,
    height = 5
  )
}


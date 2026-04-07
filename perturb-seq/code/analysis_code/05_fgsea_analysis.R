##############################################################################
# Step 05: fGSEA Pathway Analysis
# PAIR-Perturb-Seq Analysis Pipeline
#
# Input:  res/03_tier2/Tier2_interaction_table.csv
#         res/04_tier3/Tier3_*_DE.csv
# Output: res/05_gsea/ (enrichment tables, GSEA plots)
#
# Runs fgsea on:
#   1) Tier 2 interaction scores (ranked by interaction_score)
#   2) Tier 3 partner DE results (ranked by avg_log2FC per partner)
##############################################################################

library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

set.seed(42)

base_dir <- "."
res_dir <- file.path(base_dir, "res", "05_gsea")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load MSigDB Hallmark gene sets ---
cat("Loading MSigDB Hallmark gene sets...\n")
msig_df <- msigdbr(species = "Homo sapiens", collection = "H")

# Convert to list format for fgsea
pathways <- split(msig_df$gene_symbol, msig_df$gs_name)
cat("Loaded", length(pathways), "Hallmark pathways\n")

# Save pathway gene sets for reference
pathway_dir <- file.path(base_dir, "data", "pathways")
dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(pathways, file.path(pathway_dir, "hallmark_pathways.rds"))

# --- Helper: run fgsea and produce output ---
run_fgsea_analysis <- function(ranks, pathways, prefix, title_str, res_dir) {
  # Remove NAs and duplicates
  ranks <- ranks[!is.na(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]

  cat("Running fgsea for:", prefix, "(", length(ranks), "genes)...\n")

  fgsea_res <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = 15,
    maxSize = 500
  )

  fgsea_res <- fgsea_res %>% arrange(pval)

  n_sig <- sum(fgsea_res$padj < 0.05)
  cat("  Significant pathways (padj < 0.05):", n_sig, "\n")

  # Save full results
  # Convert leadingEdge list column to string for CSV export
  fgsea_export <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
  write.csv(fgsea_export, file.path(res_dir, paste0(prefix, "_fgsea_results.csv")),
    row.names = FALSE
  )

  # Bar plot of top pathways by NES
  top_pathways <- fgsea_res %>%
    filter(padj < 0.25) %>%
    arrange(NES)

  if (nrow(top_pathways) > 0) {
    # Clean pathway names
    top_pathways$pathway_clean <- gsub("HALLMARK_", "", top_pathways$pathway)
    top_pathways$pathway_clean <- gsub("_", " ", top_pathways$pathway_clean)
    top_pathways$pathway_clean <- factor(top_pathways$pathway_clean,
      levels = top_pathways$pathway_clean
    )
    top_pathways$direction <- ifelse(top_pathways$NES > 0, "Up", "Down")

    p_bar <- ggplot(top_pathways, aes(x = NES, y = pathway_clean, fill = direction)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#2166AC")) +
      geom_vline(xintercept = 0, color = "grey30") +
      theme_classic(base_size = 14) +
      theme(axis.text.y = element_text(size = 11)) +
      labs(
        title = title_str,
        subtitle = paste0(n_sig, " pathways with padj < 0.05"),
        x = "Normalized Enrichment Score (NES)",
        y = "", fill = "Direction"
      )

    ggsave(file.path(res_dir, paste0(prefix, "_NES_barplot.pdf")),
      p_bar,
      width = 10, height = max(5, nrow(top_pathways) * 0.35 + 2)
    )
  }

  # GSEA table plot for top pathways (using grid, not ggsave)
  top_n_paths <- head(fgsea_res[order(fgsea_res$pval), ]$pathway, 10)
  if (length(top_n_paths) > 0) {
    pdf(file.path(res_dir, paste0(prefix, "_gsea_table.pdf")), width = 12, height = 8)
    plotGseaTable(pathways[top_n_paths], ranks, fgsea_res, gseaParam = 0.5)
    dev.off()
  }

  return(fgsea_res)
}

# --- 1. Tier 2: Interaction Score GSEA ---
cat("\n=== Tier 2 Interaction Score GSEA ===\n")
tier2 <- read.csv(file.path(base_dir, "res", "03_tier2", "Tier2_interaction_table.csv"))

ranks_tier2 <- setNames(tier2$interaction_score, tier2$gene)
fgsea_tier2 <- run_fgsea_analysis(
  ranks_tier2, pathways,
  "Tier2_interaction",
  "Tier 2: GSEA on Interaction Scores (NBN vs WT Stress)",
  res_dir
)

# --- 2. Tier 3: Partner DE GSEA ---
cat("\n=== Tier 3 Partner DE GSEA ===\n")
partners <- c("TP53BP1", "XRCC6", "POLQ")
fgsea_tier3 <- list()

for (partner in partners) {
  cat("\n--- Partner:", partner, "---\n")
  de_file <- file.path(base_dir, "res", "04_tier3", paste0("Tier3_", partner, "_DE.csv"))

  if (!file.exists(de_file)) {
    cat("  WARNING: DE file not found:", de_file, "\n")
    next
  }

  de <- read.csv(de_file)
  ranks <- setNames(de$avg_log2FC, de$gene)

  fgsea_tier3[[partner]] <- run_fgsea_analysis(
    ranks, pathways,
    paste0("Tier3_", partner),
    paste0("Tier 3: GSEA - ", partner, " KD vs NT (RNP, NBN)"),
    res_dir
  )
}

# --- 3. Cross-partner comparison heatmap ---
cat("\n=== Cross-Partner GSEA Comparison ===\n")

# Collect NES values across partners
nes_list <- list()
for (partner in names(fgsea_tier3)) {
  res <- fgsea_tier3[[partner]]
  nes_list[[partner]] <- setNames(res$NES, res$pathway)
}

# Also add Tier 2 interaction
nes_list[["Tier2_Interaction"]] <- setNames(fgsea_tier2$NES, fgsea_tier2$pathway)

# Create NES matrix
all_pathways <- unique(unlist(lapply(nes_list, names)))
nes_mat <- matrix(0, nrow = length(all_pathways), ncol = length(nes_list))
rownames(nes_mat) <- all_pathways
colnames(nes_mat) <- names(nes_list)

for (nm in names(nes_list)) {
  nes_mat[names(nes_list[[nm]]), nm] <- nes_list[[nm]]
}

# Filter to pathways significant in at least one comparison
padj_mat <- matrix(1, nrow = length(all_pathways), ncol = length(nes_list))
rownames(padj_mat) <- all_pathways
colnames(padj_mat) <- names(nes_list)

for (nm in names(fgsea_tier3)) {
  res <- fgsea_tier3[[nm]]
  padj_mat[res$pathway, nm] <- res$padj
}
padj_mat[fgsea_tier2$pathway, "Tier2_Interaction"] <- fgsea_tier2$padj

sig_pathways <- rownames(padj_mat)[apply(padj_mat, 1, function(x) any(x < 0.1))]

if (length(sig_pathways) >= 2) {
  nes_sig <- nes_mat[sig_pathways, , drop = FALSE]
  rownames(nes_sig) <- gsub("HALLMARK_", "", rownames(nes_sig))
  rownames(nes_sig) <- gsub("_", " ", rownames(nes_sig))

  # Heatmap
  library(pheatmap)
  pdf(file.path(res_dir, "GSEA_cross_comparison_heatmap.pdf"),
    width = 8, height = max(5, length(sig_pathways) * 0.4 + 2)
  )
  pheatmap(nes_sig,
    color = colorRampPalette(c("#2166AC", "white", "#E64B35"))(100),
    breaks = seq(-3, 3, length.out = 101),
    cluster_cols = FALSE,
    main = "Hallmark GSEA: NES Across Comparisons",
    fontsize_row = 11,
    fontsize_col = 12
  )
  dev.off()

  cat("Cross-comparison heatmap generated with", length(sig_pathways), "pathways\n")
} else {
  cat("Fewer than 2 significant pathways across comparisons — skipping heatmap\n")
}

# Save NES matrix
write.csv(nes_mat, file.path(res_dir, "GSEA_NES_matrix.csv"))

cat("\nStep 05 complete. Results in:", res_dir, "\n")

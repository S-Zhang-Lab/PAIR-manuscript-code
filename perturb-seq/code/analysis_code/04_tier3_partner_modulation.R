##############################################################################
# Step 04: Tier 3 - Partner-Induced Modulation
# PAIR-Perturb-Seq Analysis Pipeline
#
# Input:  data/seu_qc.qs
# Output: res/04_tier3/ (DE tables, volcano plots per partner)
#
# Question: Within RNP-stressed, NBN-activated cells, does secondary
# suppression of a repair partner (TP53BP1, XRCC6, POLQ) rescue or
# aggravate the transcriptional phenotype?
#
# Comparisons (all within RNP-treated, NBN-CRISPRa cells):
#   RNP_NBN_TP53BP1 vs RNP_NBN_NT
#   RNP_NBN_XRCC6   vs RNP_NBN_NT
#   RNP_NBN_POLQ    vs RNP_NBN_NT
##############################################################################
rm(list = ls())
library(Seurat)
library(qs)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

set.seed(42)

# --- Source config.R (sets DATA_ROOT, RES_DIR, SEU_QC; see code/config.R) ---
if (file.exists("code/config.R")) {
  source("code/config.R")                 # run from repo/
} else if (file.exists("../config.R")) {
  source("../config.R")                   # run from repo/code/analysis_code/
} else {
  stop("config.R not found. Run from repo/ or repo/code/analysis_code/.")
}
# ---------------------------------------------------------------------------

res_dir <- file.path(RES_DIR, "04_tier3")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading seu_qc.qs...\n")
seu <- qread(SEU_QC)

# --- Subset to RNP-treated, NBN-CRISPRa cells ---
nbn_tags <- c(
  "NBN_CRISPRa_NT_CasRx", "NBN_CRISPRa_53BP1_CasRx",
  "NBN_CRISPRa_KU70_CasRx", "NBN_CRISPRa_PQ_CasRx"
)
seu_nbn_rnp <- subset(seu, treatment == "RNP" & assigned_tag %in% nbn_tags)
Idents(seu_nbn_rnp) <- "assigned_tag"

cat("RNP NBN cells:", ncol(seu_nbn_rnp), "\n")
cat("Breakdown:\n")
print(table(seu_nbn_rnp$assigned_tag))

# --- Define comparisons ---
comparisons <- list(
  "TP53BP1" = list(
    ident.1 = "NBN_CRISPRa_53BP1_CasRx",
    ident.2 = "NBN_CRISPRa_NT_CasRx",
    label   = "TP53BP1 KD"
  ),
  "XRCC6" = list(
    ident.1 = "NBN_CRISPRa_KU70_CasRx",
    ident.2 = "NBN_CRISPRa_NT_CasRx",
    label   = "XRCC6 KD"
  ),
  "POLQ" = list(
    ident.1 = "NBN_CRISPRa_PQ_CasRx",
    ident.2 = "NBN_CRISPRa_NT_CasRx",
    label   = "POLQ KD"
  )
)

# --- Run DE for each partner ---
all_de <- list()
volcano_plots <- list()

for (partner in names(comparisons)) {
  comp <- comparisons[[partner]]

  n1 <- sum(seu_nbn_rnp$assigned_tag == comp$ident.1)
  n2 <- sum(seu_nbn_rnp$assigned_tag == comp$ident.2)
  cat("\n--- ", comp$label, ": ", comp$ident.1, " (n=", n1,
    ") vs ", comp$ident.2, " (n=", n2, ") ---\n",
    sep = ""
  )

  de <- FindMarkers(
    seu_nbn_rnp,
    ident.1 = comp$ident.1,
    ident.2 = comp$ident.2,
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0,
    verbose = FALSE
  )
  de$gene <- rownames(de)
  de$partner <- partner
  de <- de %>% arrange(p_val_adj, desc(abs(avg_log2FC)))

  n_sig <- sum(de$p_val_adj < 0.05)
  n_sig_fc <- sum(de$p_val_adj < 0.05 & abs(de$avg_log2FC) > 0.5)
  cat("Genes tested:", nrow(de), "\n")
  cat("Significant (adj.P < 0.05):", n_sig, "\n")
  cat("Significant (adj.P < 0.05 & |FC| > 0.5):", n_sig_fc, "\n")

  # Save individual DE table
  write.csv(de, file.path(res_dir, paste0("Tier3_", partner, "_DE.csv")), row.names = FALSE)
  all_de[[partner]] <- de

  # --- Volcano plot ---
  de$class <- "Not Sig"
  de$class[de$p_val_adj < 0.05 & de$avg_log2FC > 0.5] <- paste0("Up in ", comp$label)
  de$class[de$p_val_adj < 0.05 & de$avg_log2FC < -0.5] <- paste0("Down in ", comp$label)

  min_nonzero_p <- min(de$p_val_adj[de$p_val_adj > 0], na.rm = TRUE) # use p_val_adj
  de$p_plot <- ifelse(de$p_val_adj == 0, min_nonzero_p * 0.1, de$p_val_adj)

  # Top genes to label
  top_up <- de %>%
    filter(grepl("^Up", class)) %>%
    top_n(10, wt = -p_val_adj)
  top_down <- de %>%
    filter(grepl("^Down", class)) %>%
    top_n(10, wt = -p_val_adj)
  if (nrow(top_up) + nrow(top_down) == 0) {
    top_up <- de %>%
      filter(avg_log2FC > 0) %>%
      top_n(5, wt = -p_val_adj)
    top_down <- de %>%
      filter(avg_log2FC < 0) %>%
      top_n(5, wt = -p_val_adj)
  }
  label_genes <- rbind(top_up, top_down)

  up_color <- "#E64B35"
  down_color <- "#4DBBD5"
  color_vals <- setNames(
    c(up_color, down_color, "grey70"),
    c(paste0("Up in ", comp$label), paste0("Down in ", comp$label), "Not Sig")
  )

  p <- ggplot(de, aes(x = avg_log2FC, y = -log10(p_plot), color = class)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = color_vals) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = label_genes, aes(label = gene),
      size = 3.5, max.overlaps = 20, color = "black"
    ) +
    theme_classic(base_size = 14) +
    theme(aspect.ratio = 1) +
    labs(
      title = paste0("Tier 3: ", comp$label, " vs NT (RNP, NBN)"),
      subtitle = paste0(
        n1, " vs ", n2, " cells | ",
        n_sig_fc, " genes (adj.P<0.05, |FC|>0.5)"
      ),
      x = "Log2 Fold Change", y = "-Log10(P-value)", color = ""
    )

  volcano_plots[[partner]] <- p
  ggsave(file.path(res_dir, paste0("Tier3_Volcano_", partner, ".pdf")),
    p,
    width = 7, height = 7
  )
}

# --- Combined volcano panel ---
p_combined <- wrap_plots(volcano_plots, ncol = 3) +
  plot_annotation(title = "Tier 3: Partner Modulation in RNP-Stressed NBN Cells")
ggsave(file.path(res_dir, "Tier3_Volcano_combined.pdf"), p_combined, width = 22, height = 8)

# --- Combined DE table ---
all_de_df <- do.call(rbind, all_de)
write.csv(all_de_df, file.path(res_dir, "Tier3_all_partners_DE.csv"), row.names = FALSE)

# --- Summary ---
cat("\n=== Tier 3 Summary ===\n")
for (partner in names(comparisons)) {
  de <- all_de[[partner]]
  sig <- sum(de$p_val_adj < 0.05 & abs(de$avg_log2FC) > 0.5)
  up <- sum(de$p_val_adj < 0.05 & de$avg_log2FC > 0.5)
  down <- sum(de$p_val_adj < 0.05 & de$avg_log2FC < -0.5)
  cat(partner, ": ", nrow(de), " genes tested | ", sig,
    " sig (", up, " up, ", down, " down)\n",
    sep = ""
  )
}

cat("\nStep 04 complete. Results in:", res_dir, "\n")

# =============================================================================
# A4 — Cross-Partner Gene Set Overlap (UpSet / Jaccard Analysis)
# Inspired by: Fielden et al. 2025 Nature — DDR synthetic lethality network
# =============================================================================
.libPaths(c("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library", .libPaths()))

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(stringr)
library(UpSetR)
library(fgsea)

# --- Source config.R (sets DATA_ROOT, RES_DIR, DATA_DIR; see code/config.R) --
if (file.exists("code/config.R")) {
  source("code/config.R")                 # run from repo/
} else if (file.exists("../config.R")) {
  source("../config.R")                   # run from repo/code/analysis_code/
} else {
  stop("config.R not found. Run from repo/ or repo/code/analysis_code/.")
}
# ---------------------------------------------------------------------------

# Output directories
out_dir <- file.path(OUTPUT_DIR, "09_cross_partner_overlap")
fig_dir <- file.path(FIGURES_DIR, "09_cross_partner_overlap")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load Tier3 DE results ──────────────────────────────────────────────────
tp53bp1 <- read.csv(file.path(OUTPUT_DIR, "04_tier3/Tier3_TP53BP1_DE.csv"))
xrcc6 <- read.csv(file.path(OUTPUT_DIR, "04_tier3/Tier3_XRCC6_DE.csv"))
polq <- read.csv(file.path(OUTPUT_DIR, "04_tier3/Tier3_POLQ_DE.csv"))

# Use lenient thresholds first; if very few genes, fall back
get_sig_genes <- function(df, fc_thresh = 0.5, pval_thresh = 0.05) {
  genes <- df$gene[df$p_val_adj < pval_thresh & abs(df$avg_log2FC) > fc_thresh]
  if (length(genes) < 20) {
    message("  Few genes at strict threshold, using p_val<0.05 & |FC|>0.3")
    genes <- df$gene[df$p_val < pval_thresh & abs(df$avg_log2FC) > 0.3]
  }
  unique(genes)
}

genes_tp53bp1 <- get_sig_genes(tp53bp1)
genes_xrcc6 <- get_sig_genes(xrcc6)
genes_polq <- get_sig_genes(polq)

cat("Significant DE genes:\n")
cat("  TP53BP1:", length(genes_tp53bp1), "\n")
cat("  XRCC6:  ", length(genes_xrcc6), "\n")
cat("  POLQ:   ", length(genes_polq), "\n")

# ── 2. UpSet Plot ─────────────────────────────────────────────────────────────
all_genes <- unique(c(genes_tp53bp1, genes_xrcc6, genes_polq))

upset_df <- data.frame(
  gene     = all_genes,
  TP53BP1  = as.integer(all_genes %in% genes_tp53bp1),
  XRCC6    = as.integer(all_genes %in% genes_xrcc6),
  POLQ     = as.integer(all_genes %in% genes_polq)
)

pdf(file.path(fig_dir, "A4_upset_plot.pdf"), width = 8, height = 5)
upset(
  upset_df,
  sets = c("TP53BP1", "XRCC6", "POLQ"),
  order.by = "freq",
  sets.bar.color = c("#E41A1C", "#377EB8", "#4DAF4A"),
  main.bar.color = "gray30",
  text.scale = 1.4,
  mb.ratio = c(0.55, 0.45),
  mainbar.y.label = "Number of DE Genes",
  sets.x.label = "Set Size"
)
dev.off()
cat("Saved: A4_upset_plot.pdf\n")

# ── 3. Jaccard Similarity Matrix ──────────────────────────────────────────────
jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) {
    return(0)
  }
  inter / union
}

partners <- list(TP53BP1 = genes_tp53bp1, XRCC6 = genes_xrcc6, POLQ = genes_polq)
pnames <- names(partners)
jmat <- matrix(NA, 3, 3, dimnames = list(pnames, pnames))
for (i in pnames) for (j in pnames) jmat[i, j] <- jaccard(partners[[i]], partners[[j]])

cat("\nJaccard similarity matrix:\n")
print(round(jmat, 3))
write.csv(jmat, file.path(out_dir, "A4_jaccard_matrix.csv"))

# Plot Jaccard heatmap
jdf <- as.data.frame(as.table(jmat)) %>% rename(Partner1 = Var1, Partner2 = Var2, Jaccard = Freq)

pdf(file.path(fig_dir, "A4_jaccard_heatmap.pdf"), width = 5, height = 4)
ggplot(jdf, aes(x = Partner1, y = Partner2, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Jaccard, 3)), size = 4) +
  scale_fill_gradient(low = "white", high = "#377EB8", limits = c(0, 1)) +
  labs(
    title = "Jaccard Similarity: DE Gene Overlap Across Partners",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 11, face = "bold")
  )
dev.off()
cat("Saved: A4_jaccard_heatmap.pdf\n")

# ── 4. Gene Class Assignment ──────────────────────────────────────────────────
# Core: in all 3; Partner-unique: in only 1; Pairwise: in exactly 2
upset_df <- upset_df %>%
  mutate(
    n_partners = TP53BP1 + XRCC6 + POLQ,
    gene_class = case_when(
      n_partners == 3 ~ "Core (all 3)",
      n_partners == 2 & TP53BP1 & XRCC6 ~ "TP53BP1+XRCC6",
      n_partners == 2 & TP53BP1 & POLQ ~ "TP53BP1+POLQ",
      n_partners == 2 & XRCC6 & POLQ ~ "XRCC6+POLQ",
      n_partners == 1 & TP53BP1 ~ "TP53BP1-unique",
      n_partners == 1 & XRCC6 ~ "XRCC6-unique",
      n_partners == 1 & POLQ ~ "POLQ-unique",
      TRUE ~ "Other"
    )
  )

class_summary <- upset_df %>%
  count(gene_class) %>%
  arrange(desc(n))
cat("\nGene class distribution:\n")
print(class_summary)
write.csv(upset_df, file.path(out_dir, "A4_gene_classes.csv"), row.names = FALSE)

# ── 5. fgsea on each gene class ───────────────────────────────────────────────
pathways <- readRDS(file.path(DATA_DIR, "pathways/hallmark_pathways.rds"))

# Build ranked gene list per class from combined FC
all_de <- bind_rows(
  tp53bp1 %>% select(gene, avg_log2FC) %>% mutate(partner = "TP53BP1"),
  xrcc6 %>% select(gene, avg_log2FC) %>% mutate(partner = "XRCC6"),
  polq %>% select(gene, avg_log2FC) %>% mutate(partner = "POLQ")
)

run_gsea_for_class <- function(class_name, gene_set, all_de_df) {
  gene_vec <- gene_set
  if (length(gene_vec) < 5) {
    message("  Skipping ", class_name, ": too few genes (", length(gene_vec), ")")
    return(NULL)
  }
  stats <- all_de_df %>%
    filter(gene %in% gene_vec) %>%
    group_by(gene) %>%
    summarise(stat = mean(avg_log2FC), .groups = "drop") %>%
    arrange(desc(stat))
  ranked <- setNames(stats$stat, stats$gene)
  set.seed(42)
  tryCatch(
    fgsea(pathways = pathways, stats = ranked, minSize = 5, maxSize = 500, nPerm = 1000) %>%
      filter(padj < 0.25) %>%
      arrange(NES) %>%
      mutate(gene_class = class_name),
    error = function(e) {
      message("  fgsea error for ", class_name, ": ", e$message)
      NULL
    }
  )
}

class_gene_sets <- list(
  "Core (all 3)"   = intersect(intersect(genes_tp53bp1, genes_xrcc6), genes_polq),
  "TP53BP1-unique" = setdiff(genes_tp53bp1, union(genes_xrcc6, genes_polq)),
  "XRCC6-unique"   = setdiff(genes_xrcc6, union(genes_tp53bp1, genes_polq)),
  "POLQ-unique"    = setdiff(genes_polq, union(genes_tp53bp1, genes_xrcc6))
)

gsea_results <- map2_dfr(
  names(class_gene_sets), class_gene_sets,
  ~ run_gsea_for_class(.x, .y, all_de)
)

if (nrow(gsea_results) > 0) {
  write.csv(gsea_results %>% select(-leadingEdge), file.path(out_dir, "A4_class_fgsea.csv"), row.names = FALSE)

  pdf(file.path(fig_dir, "A4_class_fgsea_NES.pdf"), width = 10, height = max(4, nrow(gsea_results) * 0.3 + 2))
  p <- ggplot(gsea_results, aes(x = NES, y = reorder(pathway, NES), fill = gene_class)) +
    geom_col(alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~gene_class, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title = "Hallmark GSEA by Gene Overlap Class",
      x = "Normalized Enrichment Score", y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"))
  print(p)
  dev.off()
  cat("Saved: A4_class_fgsea_NES.pdf\n")
} else {
  cat("No significant GSEA hits found for gene classes.\n")
}

# ── 6. Save summary markdown ──────────────────────────────────────────────────
n_core <- sum(upset_df$gene_class == "Core (all 3)")
n_tp_uniq <- sum(upset_df$gene_class == "TP53BP1-unique")
n_xr_uniq <- sum(upset_df$gene_class == "XRCC6-unique")
n_pq_uniq <- sum(upset_df$gene_class == "POLQ-unique")

md <- sprintf(
  "# A4 — Cross-Partner Gene Set Overlap
**Date:** %s

## Thresholds Used
- Primary: adj.P < 0.05, |log2FC| > 0.5
- Fallback (if <20 genes): p_val < 0.05, |log2FC| > 0.3

## Significant DE Genes
| Partner  | N genes |
|----------|---------|
| TP53BP1  | %d      |
| XRCC6    | %d      |
| POLQ     | %d      |

## Jaccard Similarity
```
%s
```

## Gene Classes
| Class            | N genes |
|------------------|---------|
| Core (all 3)     | %d      |
| TP53BP1-unique   | %d      |
| XRCC6-unique     | %d      |
| POLQ-unique      | %d      |

## Outputs
- `A4_upset_plot.pdf` — Intersection counts across partners
- `A4_jaccard_heatmap.pdf` — Pairwise similarity
- `A4_gene_classes.csv` — Per-gene class assignments
- `A4_class_fgsea.csv` / `A4_class_fgsea_NES.pdf` — Pathway enrichment per class
",
  Sys.Date(),
  length(genes_tp53bp1), length(genes_xrcc6), length(genes_polq),
  paste(capture.output(round(jmat, 3)), collapse = "\n"),
  n_core, n_tp_uniq, n_xr_uniq, n_pq_uniq
)
writeLines(md, file.path(out_dir, "A4_cross_partner_overlap_summary.md"))
cat("Done. Tables -> ", out_dir, " | Figures -> ", fig_dir, "\n")

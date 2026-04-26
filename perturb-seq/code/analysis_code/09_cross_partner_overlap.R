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
library(patchwork)
library(tibble)
library(ggrepel)

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

pdf(file.path(fig_dir, "A4_upset_plot.pdf"), width = 6, height = 4.2, onefile = FALSE)
upset(
  upset_df,
  sets = c("TP53BP1", "XRCC6", "POLQ"),
  order.by = "freq",
  sets.bar.color = c("#E41A1C", "#377EB8", "#4DAF4A"),
  main.bar.color = "gray30",
  text.scale = 1.3,
  mb.ratio = c(0.72, 0.28),
  point.size = 4.5,
  line.size = 1.2,
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

# ── 7. Pathway-Level UpSet: Down-regulated Hallmark Pathways ─────────────────
cat("\n--- Pathway UpSet: down-regulated Hallmark pathways ---\n")

fgsea_tp <- read.csv(file.path(OUTPUT_DIR, "05_fgsea/Tier3_TP53BP1_fgsea_results.csv"))
fgsea_xr <- read.csv(file.path(OUTPUT_DIR, "05_fgsea/Tier3_XRCC6_fgsea_results.csv"))
fgsea_pq <- read.csv(file.path(OUTPUT_DIR, "05_fgsea/Tier3_POLQ_fgsea_results.csv"))

padj_cut <- 0.25

down_tp <- fgsea_tp$pathway[fgsea_tp$NES < 0 & fgsea_tp$padj < padj_cut]
down_xr <- fgsea_xr$pathway[fgsea_xr$NES < 0 & fgsea_xr$padj < padj_cut]
down_pq <- fgsea_pq$pathway[fgsea_pq$NES < 0 & fgsea_pq$padj < padj_cut]

cat("Down-regulated pathways (padj <", padj_cut, "):\n")
cat("  TP53BP1:", length(down_tp), "\n")
cat("  XRCC6:  ", length(down_xr), "\n")
cat("  POLQ:   ", length(down_pq), "\n")

all_paths <- unique(c(down_tp, down_xr, down_pq))
if (length(all_paths) < 2) {
  cat("Too few pathways to build UpSet plot — skipping.\n")
} else {
  # Partner aesthetics
  partner_cols  <- c(TP53BP1 = "#A31621", XRCC6 = "#C07A2C", POLQ = "#6B3FA0")
  partner_names <- c("TP53BP1", "XRCC6", "POLQ")
  # y positions: TP53BP1=3 (top), XRCC6=2, POLQ=1 (bottom)
  partner_y <- c(TP53BP1 = 3L, XRCC6 = 2L, POLQ = 1L)

  # Membership table
  mem <- tibble(
    pathway = all_paths,
    TP53BP1 = all_paths %in% down_tp,
    XRCC6   = all_paths %in% down_xr,
    POLQ    = all_paths %in% down_pq
  ) %>%
    mutate(
      n_sets = TP53BP1 + XRCC6 + POLQ,
      set_id = paste0(
        ifelse(TP53BP1, "T", ""),
        ifelse(XRCC6,   "X", ""),
        ifelse(POLQ,    "P", "")
      )
    )

  # Intersection summary, ordered by count descending
  inter_sum <- mem %>%
    count(set_id, TP53BP1, XRCC6, POLQ, n_sets) %>%
    arrange(desc(n)) %>%
    mutate(
      x_pos = row_number(),
      bar_color = case_when(
        set_id == "T"   ~ partner_cols["TP53BP1"],
        set_id == "X"   ~ partner_cols["XRCC6"],
        set_id == "P"   ~ partner_cols["POLQ"],
        TRUE            ~ "gray40"
      ),
      x_label = case_when(
        set_id == "T"   ~ "Unique to\nTP53BP1",
        set_id == "X"   ~ "Unique to\nXRCC6",
        set_id == "P"   ~ "Unique to\nPOLQ",
        set_id == "TX"  ~ "TP53BP1\n+XRCC6",
        set_id == "TP"  ~ "TP53BP1\n+POLQ",
        set_id == "XP"  ~ "XRCC6\n+POLQ",
        set_id == "TXP" ~ "All\nThree",
        TRUE ~ set_id
      )
    )

  # ── Top panel: intersection bar chart ──────────────────────────────────────
  p_bars <- ggplot(inter_sum, aes(x = x_pos, y = n, fill = bar_color)) +
    geom_col(width = 0.55) +
    geom_text(aes(label = n), vjust = -0.45, size = 3.8, fontface = "bold",
              color = "gray20") +
    scale_fill_identity() +
    scale_x_continuous(
      breaks   = inter_sum$x_pos,
      labels   = inter_sum$x_label,
      position = "top",
      expand   = expansion(add = 0.7)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "Cross-Partner UpSet Plot",
      subtitle = "Down-regulated Hallmark pathways (NES < 0, padj < 0.25)",
      x = NULL, y = "Number of\nPathways"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x  = element_text(size = 8.5, color = "gray20", lineheight = 1.2),
      axis.line.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title   = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray45"),
      plot.margin  = margin(t = 8, r = 8, b = 0, l = 8)
    )

  # ── Dot matrix data ─────────────────────────────────────────────────────────
  dot_df <- expand.grid(
    x_pos   = inter_sum$x_pos,
    partner = partner_names,
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      inter_sum %>% select(x_pos, TP53BP1, XRCC6, POLQ),
      by = "x_pos"
    ) %>%
    mutate(
      active   = case_when(
        partner == "TP53BP1" ~ TP53BP1,
        partner == "XRCC6"   ~ XRCC6,
        partner == "POLQ"    ~ POLQ
      ),
      y_pos    = partner_y[partner],
      dot_fill = ifelse(active, partner_cols[partner], "gray88"),
      dot_size = ifelse(active, 5.5, 3.8)
    )

  # Connector lines for multi-partner intersections
  line_df <- inter_sum %>%
    filter(n_sets > 1) %>%
    mutate(
      y_min = pmin(
        ifelse(TP53BP1, partner_y["TP53BP1"], Inf),
        ifelse(XRCC6,   partner_y["XRCC6"],   Inf),
        ifelse(POLQ,    partner_y["POLQ"],     Inf)
      ),
      y_max = pmax(
        ifelse(TP53BP1, partner_y["TP53BP1"], -Inf),
        ifelse(XRCC6,   partner_y["XRCC6"],   -Inf),
        ifelse(POLQ,    partner_y["POLQ"],     -Inf)
      )
    )

  # ── Bottom panel: dot matrix ─────────────────────────────────────────────────
  stripe_df <- tibble(
    y     = 1:3,
    fill  = c("gray96", "white", "gray96")   # POLQ, XRCC6, TP53BP1
  )

  p_dots <- ggplot() +
    geom_rect(
      data = stripe_df,
      aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5, fill = fill),
      inherit.aes = FALSE
    ) +
    scale_fill_identity(guide = "none") +
    { if (nrow(line_df) > 0)
        geom_segment(
          data = line_df,
          aes(x = x_pos, xend = x_pos, y = y_min, yend = y_max),
          color = "gray35", linewidth = 1.4
        )
    } +
    geom_point(
      data  = dot_df,
      aes(x = x_pos, y = y_pos, fill = dot_fill, size = dot_size),
      shape = 21, color = "white", stroke = 0.4
    ) +
    scale_size_identity() +
    scale_x_continuous(
      breaks = inter_sum$x_pos,
      labels = inter_sum$n,
      expand = expansion(add = 0.7)
    ) +
    scale_y_continuous(
      breaks = 1:3,
      labels = rev(partner_names),
      limits = c(0.5, 3.5)
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y  = element_text(face = "bold", size = 11, color = "gray15"),
      axis.text.x  = element_text(size = 10, color = "gray30"),
      axis.line    = element_blank(),
      axis.ticks   = element_blank(),
      plot.margin  = margin(t = 0, r = 8, b = 8, l = 8)
    )

  # ── Combine and save ─────────────────────────────────────────────────────────
  p_pathway_upset <- p_bars / p_dots + plot_layout(heights = c(2.5, 1.5))

  n_cols <- nrow(inter_sum)
  pdf(file.path(fig_dir, "A4_pathway_upset_downreg.pdf"),
      width = max(6, n_cols * 1.4 + 2), height = 6)
  print(p_pathway_upset)
  dev.off()
  cat("Saved: A4_pathway_upset_downreg.pdf\n")

  # Save intersection table
  write.csv(
    mem %>% arrange(desc(n_sets), set_id),
    file.path(out_dir, "A4_pathway_intersections_downreg.csv"),
    row.names = FALSE
  )

  # ── 8. Grouped NES Heatmap (PRIMARY) + Dot Plot (SUPPLEMENT) ───────────────
  # Goal: name every pathway, group by intersection class, show NES + padj

  # Full per-partner NES/padj for every pathway in the union
  get_stats <- function(df, partner_label) {
    df %>%
      filter(pathway %in% all_paths) %>%
      transmute(pathway, partner = partner_label, NES, padj)
  }

  stats_df <- bind_rows(
    get_stats(fgsea_tp, "TP53BP1"),
    get_stats(fgsea_xr, "XRCC6"),
    get_stats(fgsea_pq, "POLQ")
  )

  # Ensure every (pathway, partner) pair exists (fill missing with NA)
  stats_df <- stats_df %>%
    complete(pathway = all_paths, partner = partner_names) %>%
    mutate(
      partner = factor(partner, levels = partner_names),
      sig_label = case_when(
        is.na(padj)      ~ "",
        padj < 0.001     ~ "***",
        padj < 0.01      ~ "**",
        padj < 0.05      ~ "*",
        padj < 0.25      ~ "\u00B7",   # middle dot for padj<0.25
        TRUE             ~ ""
      )
    )

  # Intersection class per pathway (inherited from `mem`)
  class_lookup <- mem %>%
    mutate(
      intersection_class = case_when(
        n_sets == 3                               ~ "Core (all 3)",
        set_id == "XP"                            ~ "XRCC6 + POLQ",
        set_id == "TX"                            ~ "TP53BP1 + XRCC6",
        set_id == "TP"                            ~ "TP53BP1 + POLQ",
        set_id == "T"                             ~ "TP53BP1-unique",
        set_id == "X"                             ~ "XRCC6-unique",
        set_id == "P"                             ~ "POLQ-unique",
        TRUE                                      ~ "Other"
      )
    ) %>%
    select(pathway, intersection_class, n_sets, set_id)

  # Class ordering: Core first, then pairwise, then partner-unique
  class_levels <- c(
    "Core (all 3)",
    "TP53BP1 + XRCC6", "TP53BP1 + POLQ", "XRCC6 + POLQ",
    "TP53BP1-unique", "XRCC6-unique", "POLQ-unique"
  )
  class_levels <- class_levels[class_levels %in% class_lookup$intersection_class]

  # Clean pathway labels (strip HALLMARK_ and replace _ with space)
  plot_df <- stats_df %>%
    left_join(class_lookup, by = "pathway") %>%
    mutate(
      intersection_class = factor(intersection_class, levels = class_levels),
      pathway_clean = gsub("^HALLMARK_", "", pathway),
      pathway_clean = gsub("_", " ", pathway_clean)
    )

  # Within each class, order pathways by mean NES (most negative at top)
  pathway_order <- plot_df %>%
    group_by(pathway_clean, intersection_class) %>%
    summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    arrange(intersection_class, mean_NES) %>%
    pull(pathway_clean)

  plot_df$pathway_clean <- factor(plot_df$pathway_clean, levels = rev(pathway_order))

  # Symmetric NES color scale
  nes_max <- max(abs(plot_df$NES), na.rm = TRUE)
  nes_lim <- c(-nes_max, nes_max)

  # ── PRIMARY: grouped NES heatmap ─────────────────────────────────────────────
  p_heat <- ggplot(plot_df, aes(x = partner, y = pathway_clean, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sig_label), size = 4.2, color = "gray10", vjust = 0.75) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = nes_lim,
      name     = "NES",
      na.value = "gray92"
    ) +
    facet_grid(
      rows   = vars(intersection_class),
      scales = "free_y",
      space  = "free_y",
      switch = "y"
    ) +
    labs(
      title    = "Down-regulated Hallmark pathways across repair-partner knockdowns",
      subtitle = "Grouped by intersection class; *** padj<0.001, ** <0.01, * <0.05, \u00B7 <0.25",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid         = element_blank(),
      axis.text.x        = element_text(face = "bold", size = 11, color = "gray15"),
      axis.text.y        = element_text(size = 9.5, color = "gray15"),
      strip.text.y.left  = element_text(angle = 0, face = "bold", size = 9.5,
                                        hjust = 1, color = "gray20"),
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "gray94", color = NA),
      panel.spacing.y    = unit(3, "pt"),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 9, color = "gray45"),
      legend.position    = "right"
    )

  n_paths <- length(unique(plot_df$pathway_clean))
  pdf(file.path(fig_dir, "A4_pathway_heatmap_downreg.pdf"),
      width = 8.5, height = max(4.5, n_paths * 0.32 + 2))
  print(p_heat)
  dev.off()
  cat("Saved: A4_pathway_heatmap_downreg.pdf\n")

  # ── SUPPLEMENT: dot plot (bubble) ────────────────────────────────────────────
  dot_plot_df <- plot_df %>%
    mutate(
      neglog10_padj = ifelse(is.na(padj), NA_real_, -log10(pmax(padj, 1e-12))),
      show_dot      = !is.na(NES)
    ) %>%
    filter(show_dot)

  p_dot <- ggplot(dot_plot_df,
                  aes(x = partner, y = pathway_clean,
                      size = neglog10_padj, fill = NES)) +
    geom_point(shape = 21, color = "gray30", stroke = 0.3) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = nes_lim, name = "NES"
    ) +
    scale_size_continuous(
      name   = expression(-log[10]~padj),
      range  = c(1.5, 9),
      breaks = c(1, 2, 4, 8)
    ) +
    facet_grid(
      rows   = vars(intersection_class),
      scales = "free_y",
      space  = "free_y",
      switch = "y"
    ) +
    labs(
      title    = "Down-regulated Hallmark pathways: effect size and significance",
      subtitle = "Dot size = -log10(padj); fill = NES",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major   = element_line(color = "gray92", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(face = "bold", size = 11, color = "gray15"),
      axis.text.y        = element_text(size = 9.5, color = "gray15"),
      strip.text.y.left  = element_text(angle = 0, face = "bold", size = 9.5,
                                        hjust = 1, color = "gray20"),
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "gray94", color = NA),
      panel.spacing.y    = unit(3, "pt"),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 9, color = "gray45"),
      legend.position    = "right",
      legend.box         = "vertical"
    )

  pdf(file.path(fig_dir, "A4_pathway_dotplot_downreg.pdf"),
      width = 9, height = max(4.5, n_paths * 0.32 + 2))
  print(p_dot)
  dev.off()
  cat("Saved: A4_pathway_dotplot_downreg.pdf\n")

  # ── 9. OPTION A: Slope / parallel-coordinates plot ──────────────────────────
  # Each pathway = one line connecting its NES across the 3 partners,
  # colored by intersection class. Reveals "consistent-across-all" vs
  # "partner-specific" as visual shapes rather than grid values.

  class_palette <- c(
    "Core (all 3)"     = "#111111",
    "TP53BP1 + XRCC6"  = "#6A3D9A",
    "TP53BP1 + POLQ"   = "#FF7F00",
    "XRCC6 + POLQ"     = "#1F78B4",
    "TP53BP1-unique"   = "#A31621",
    "XRCC6-unique"     = "#C07A2C",
    "POLQ-unique"      = "#6B3FA0"
  )
  class_palette <- class_palette[names(class_palette) %in% class_levels]

  slope_df <- plot_df %>%
    mutate(
      partner = factor(partner, levels = partner_names),
      is_core = intersection_class == "Core (all 3)"
    )

  # Endpoint labels: place on the right side for readability
  label_df <- slope_df %>%
    filter(partner == "POLQ") %>%
    mutate(lbl = pathway_clean)

  p_slope <- ggplot(slope_df,
                    aes(x = partner, y = NES,
                        group = pathway_clean,
                        color = intersection_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.4) +
    geom_line(aes(linewidth = is_core, alpha = is_core)) +
    geom_point(aes(size = is_core), shape = 16) +
    ggrepel::geom_text_repel(
      data          = label_df,
      aes(label = lbl),
      nudge_x       = 0.35,
      direction     = "y",
      hjust         = 0,
      segment.size  = 0.2,
      segment.color = "gray70",
      size          = 3,
      max.overlaps  = Inf,
      xlim          = c(3.2, NA)
    ) +
    scale_color_manual(values = class_palette, name = "Intersection class") +
    scale_linewidth_manual(values = c(`TRUE` = 1.3, `FALSE` = 0.55), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.75), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 2.6, `FALSE` = 1.8), guide = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.08, 0.55))) +
    labs(
      title    = "Trajectory of Hallmark NES across repair-partner knockdowns",
      subtitle = "Down-regulated pathways (padj < 0.25 in \u22651 partner). Core pathways bold.",
      x = NULL, y = "Normalized Enrichment Score (NES)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x   = element_text(face = "bold", size = 11, color = "gray15"),
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray45"),
      legend.position = "right"
    )

  pdf(file.path(fig_dir, "A4_pathway_slope_downreg.pdf"),
      width = 10, height = max(5, n_paths * 0.25 + 2))
  print(p_slope)
  dev.off()
  cat("Saved: A4_pathway_slope_downreg.pdf\n")

  # ── 10. OPTION B: Diverging lollipop, faceted by intersection class ─────────
  # One colored dot per (pathway, partner); rows grouped into intersection
  # facets. Direction and magnitude on one axis; facets carry the structure.

  partner_palette <- partner_cols  # already defined above

  lolli_df <- plot_df %>%
    filter(!is.na(NES)) %>%
    mutate(partner = factor(partner, levels = rev(partner_names)))  # for dodge

  # Within each class, order pathways by mean NES (most negative at bottom of facet)
  lolli_order <- lolli_df %>%
    group_by(pathway_clean, intersection_class) %>%
    summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    arrange(intersection_class, desc(mean_NES)) %>%
    pull(pathway_clean)
  lolli_df$pathway_clean <- factor(lolli_df$pathway_clean, levels = lolli_order)

  p_lolli <- ggplot(lolli_df,
                    aes(x = NES, y = pathway_clean, color = partner)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray70", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = NES,
                     y = pathway_clean, yend = pathway_clean),
                 position = position_dodge(width = 0.7),
                 linewidth = 0.6, alpha = 0.55) +
    geom_point(aes(size = -log10(pmax(padj, 1e-12))),
               position = position_dodge(width = 0.7)) +
    scale_color_manual(values = partner_palette, name = "Partner KD",
                       breaks = partner_names) +
    scale_size_continuous(
      name   = expression(-log[10]~padj),
      range  = c(1.5, 6),
      breaks = c(1, 2, 4, 8)
    ) +
    facet_grid(
      rows   = vars(intersection_class),
      scales = "free_y",
      space  = "free_y",
      switch = "y"
    ) +
    labs(
      title    = "NES by partner, grouped by intersection class",
      subtitle = "Down-regulated Hallmark pathways (padj < 0.25 in \u22651 partner)",
      x = "Normalized Enrichment Score (NES)", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 9.5, color = "gray15"),
      strip.text.y.left  = element_text(angle = 0, face = "bold", size = 9.5,
                                        hjust = 1, color = "gray20"),
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "gray94", color = NA),
      panel.spacing.y    = unit(4, "pt"),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 9, color = "gray45"),
      legend.position    = "right",
      legend.box         = "vertical"
    )

  pdf(file.path(fig_dir, "A4_pathway_lollipop_downreg.pdf"),
      width = 9.5, height = max(5, n_paths * 0.35 + 2.5))
  print(p_lolli)
  dev.off()
  cat("Saved: A4_pathway_lollipop_downreg.pdf\n")
}

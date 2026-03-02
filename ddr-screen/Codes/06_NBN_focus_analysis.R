## =============================================================================
## 06_NBN_focus_analysis.R
## NBN-focused visualization: enrichment bias + partner waterfall
## =============================================================================
##
## Purpose:
##   Generate two panels for main Figure 2 that bridge the screen-level
##   analysis to NBN-specific validation:
##
##   Panel A (Enrichment Bias):  For each high-frequency gene at position 1,
##     show the fraction of its significant PAIRs that are enriched in BFP+
##     cells, relative to the global baseline (101/217 = 47%).
##     Statistical assessment: per-gene binomial test, Fisher's exact test,
##     and a permutation test on the max z-score (10,000 iterations) that
##     corrects for post-hoc selection of the top-ranked gene. NBN has the
##     highest enrichment ratio (8/11 = 73%) and z-score, but none of the
##     tests reach conventional significance, consistent with limited power
##     from small per-gene sample sizes (5-19 PAIRs).
##
##   Panel B (NBN Waterfall):  All 11 DESeq2-significant PAIRs with NBN at
##     position 1, ranked by log2FC, with bars colored by the position 2
##     partner gene's DDR pathway membership (NHEJ vs HR vs DSB sensing vs
##     other). Visually demonstrates both the enrichment skew and the partner
##     coherence.
##
##   NOTE: All gene symbols use official HUGO nomenclature.
##     Common name → HUGO: DNA-PKcs → PRKDC, 53BP1 → TP53BP1,
##     KU70 → XRCC6, CtIP → RBBP8, XLF → NHEJ1.
##     MRE11 is classified as DSB sensing (MRN) alongside NBN and RAD50.
##
## Input:
##   - ./Results/DESeq2_primary_significant_PAIRs.csv
##
## Output:
##   - ./Figures/NBN_Enrichment_Bias.pdf
##   - ./Figures/NBN_Waterfall_Partners.pdf
##   - ./Figures/NBN_Combined_Panel.pdf
##   - ./Results/NBN_position1_PAIRs.csv
##   - ./Results/NBN_enrichment_bias_statistics.csv
## =============================================================================

# --- Load required libraries --------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- Official HUGO gene symbol mapping ----------------------------------------
# Maps common/colloquial names (as used in sgRNA identifiers) to official HUGO
gene_symbol_map <- c(
  "DNA-PKcs" = "PRKDC",
  "53BP1"    = "TP53BP1",
  "KU70"     = "XRCC6",
  "CtIP"     = "RBBP8",
  "XLF"      = "NHEJ1"
)
# All other names (NBN, MRE11, BRCA1, etc.) are already official HUGO symbols.

to_hugo <- function(x) {
  mapped <- gene_symbol_map[x]
  ifelse(is.na(mapped), x, mapped)
}

# --- DDR pathway annotations --------------------------------------------------
# Keyed by HUGO symbols. MRE11 classified as DSB sensing (MRN) with NBN/RAD50.
pathway_map <- data.frame(
  Gene = c("TP53BP1", "XRCC6", "NHEJ1", "XRCC4", "LIG4", "PRKDC",  # NHEJ
           "BRCA1", "BRCA2", "RAD51", "PALB2", "RBBP8",              # HR
           "NBN", "RAD50", "MRE11",                                    # DSB sensing (MRN)
           "POLQ", "LIG3",                                             # Alt-EJ / MMEJ
           "EXO1", "POLD3", "XRCC1",                                  # Other repair
           "ATM",                                                       # Signaling
           "NT"),                                                       # Non-targeting
  Pathway_group = c(rep("NHEJ", 6),
                    rep("HR", 5),
                    rep("DSB sensing (MRN)", 3),
                    rep("Alt-EJ / MMEJ", 2),
                    rep("Other repair", 3),
                    "DDR signaling",
                    "Non-targeting control"),
  stringsAsFactors = FALSE
)

# Color palette for pathways
pathway_colors <- c(
  "NHEJ"                   = "#2166AC",   # blue
  "HR"                     = "#B2182B",   # red
  "Alt-EJ / MMEJ"          = "#F4A582",   # salmon
  "DSB sensing (MRN)"      = "#7570B3",   # purple
  "DDR signaling"          = "#66A61E",   # green
  "Other repair"           = "#999999",   # gray
  "Non-targeting control"  = "#D9D9D9"    # light gray
)

# =============================================================================
# SECTION 1: Load data
# =============================================================================

deseq2_sig <- read.csv("./Results/DESeq2_primary_significant_PAIRs.csv")

cat("Loaded", nrow(deseq2_sig), "DESeq2-significant PAIRs\n")
cat("  Enriched:", sum(deseq2_sig$direction == "enriched"), "\n")
cat("  Depleted:", sum(deseq2_sig$direction == "depleted"), "\n")

# Parse gene names from sgRNA identifiers and convert to HUGO symbols
parse_pair <- function(sgRNA) {
  p1_raw <- sub("_CRISPRa.*", "", sgRNA)
  p2_raw <- sub(".*_CRISPRa_\\d+_", "", sgRNA)
  p2_raw <- sub("_CasRx_.*", "", p2_raw)
  data.frame(sgRNA = sgRNA,
             Gene_P1 = to_hugo(p1_raw),
             Gene_P2 = to_hugo(p2_raw),
             stringsAsFactors = FALSE)
}

sig_parsed <- parse_pair(deseq2_sig$sgRNA)
deseq2_sig <- merge(deseq2_sig, sig_parsed, by = "sgRNA")

# =============================================================================
# SECTION 2: Panel A — Enrichment ratio by gene at position 1
# =============================================================================
# For each gene that appears at position 1 in ≥ N significant PAIRs,
# compute: (# enriched) / (# total significant) and compare to baseline.

n_enriched_total <- sum(deseq2_sig$direction == "enriched")
n_sig_total      <- nrow(deseq2_sig)
baseline_ratio   <- n_enriched_total / n_sig_total

cat("\nBaseline enrichment ratio:", round(baseline_ratio, 3),
    "(", n_enriched_total, "/", n_sig_total, ")\n")

# Compute per-gene enrichment ratio at position 1
gene_p1_stats <- deseq2_sig %>%
  group_by(Gene_P1) %>%
  summarise(
    n_total    = n(),
    n_enriched = sum(direction == "enriched"),
    n_depleted = sum(direction == "depleted"),
    .groups    = "drop"
  ) %>%
  mutate(
    enrichment_ratio = n_enriched / n_total,
    gene_label = paste0(Gene_P1, "\n(", n_enriched, "/", n_total, ")")
  ) %>%
  # Filter: show genes with ≥ 5 significant PAIRs at position 1
  filter(n_total >= 5) %>%
  arrange(desc(enrichment_ratio))

# ---------------------------------------------------------------------------
# Formal statistical test: enrichment bias
# ---------------------------------------------------------------------------
#
# Statistical framework (three levels):
#
#   1. Per-gene exact binomial test (one-sided, greater):
#      H0: P(enriched) = baseline_ratio. Tests each gene individually.
#      Does not account for post-hoc selection of the top-ranked gene.
#
#   2. Per-gene Fisher's exact test (one-sided, greater):
#      Compares the gene's enriched/depleted split against all OTHER genes
#      pooled. Avoids including the gene in its own null.
#
#   3. Permutation test on max z-score (primary, corrects for selection):
#      Under H0 (direction labels are exchangeable), how often does the
#      maximum sample-size-normalized enrichment z-score across all genes
#      with >= min_n PAIRs exceed the observed maximum?
#      The z-score normalizes for unequal sample sizes:
#        z_g = (k_g - n_g * p0) / sqrt(n_g * p0 * (1 - p0))
#      This is the most appropriate test because we selected NBN as the
#      top-ranked gene, making it a max-statistic problem.
#
# Conclusion: No gene reaches significance by any framework, consistent
# with the limited power of small per-gene sample sizes (5-19 PAIRs).
# NBN is selected for validation based on the convergence of:
#   (a) highest enrichment ratio and z-score among all tested genes,
#   (b) biological coherence with MRN-mediated DSB sensing, and
#   (c) partner pathway diversity (Panel B).
# ---------------------------------------------------------------------------

# --- Per-gene tests (binomial + Fisher's exact) ---
gene_p1_stats$binom_p <- mapply(function(k, n) {
  binom.test(k, n, p = baseline_ratio, alternative = "greater")$p.value
}, gene_p1_stats$n_enriched, gene_p1_stats$n_total)

gene_p1_stats$fisher_p <- mapply(function(k, n) {
  a <- k                                    # gene enriched
  b <- n - k                                # gene depleted
  c <- n_enriched_total - k                 # rest enriched
  d <- (n_sig_total - n) - (n_enriched_total - k)  # rest depleted
  fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
}, gene_p1_stats$n_enriched, gene_p1_stats$n_total)

gene_p1_stats$binom_padj  <- p.adjust(gene_p1_stats$binom_p,  method = "BH")
gene_p1_stats$fisher_padj <- p.adjust(gene_p1_stats$fisher_p, method = "BH")

# --- Per-gene z-scores (for ranking and permutation test) ---
gene_p1_stats$z_score <- (gene_p1_stats$n_enriched -
                            gene_p1_stats$n_total * baseline_ratio) /
  sqrt(gene_p1_stats$n_total * baseline_ratio * (1 - baseline_ratio))

# --- Permutation test: max z-score across genes ---
set.seed(42)
n_perm <- 10000

# Direction vector and gene assignments (full dataset, not filtered)
dir_vec   <- deseq2_sig$direction == "enriched"
gene_vec  <- deseq2_sig$Gene_P1

# Identify genes with >= 5 PAIRs (the tested set)
gene_ns <- table(gene_vec)
tested_genes <- names(gene_ns[gene_ns >= 5])
n_tested <- length(tested_genes)

# Pre-compute per-gene indices into dir_vec
gene_idx_list <- lapply(tested_genes, function(g) which(gene_vec == g))
names(gene_idx_list) <- tested_genes

# Observed max z-score
obs_max_z      <- max(gene_p1_stats$z_score)
obs_max_z_gene <- gene_p1_stats$Gene_P1[which.max(gene_p1_stats$z_score)]

cat("\n--- Permutation test for max enrichment z-score ---\n")
cat("  Genes tested (>= 5 sig PAIRs at position 1):", n_tested, "\n")
cat("  Observed max z-score:", round(obs_max_z, 3),
    "(", as.character(obs_max_z_gene), ")\n")

perm_max_z <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  perm_dir <- sample(dir_vec)  # shuffle direction labels
  max_z <- -Inf
  for (g in tested_genes) {
    idx <- gene_idx_list[[g]]
    ng  <- length(idx)
    kg  <- sum(perm_dir[idx])
    zg  <- (kg - ng * baseline_ratio) / sqrt(ng * baseline_ratio * (1 - baseline_ratio))
    if (zg > max_z) max_z <- zg
  }
  perm_max_z[i] <- max_z
}

perm_p_max_z <- sum(perm_max_z >= obs_max_z) / n_perm

cat("  Permutation p-value (n=", n_perm, "): ", round(perm_p_max_z, 4), "\n", sep = "")
cat("  Null distribution: mean =", round(mean(perm_max_z), 3),
    ", 95th pct =", round(quantile(perm_max_z, 0.95), 3), "\n")

cat("\nGenes at position 1 with >= 5 significant PAIRs:\n")
cat("  (Binomial test H0: P(enriched) = ", round(baseline_ratio, 3), ")\n", sep = "")
print(gene_p1_stats[, c("Gene_P1", "n_total", "n_enriched", "enrichment_ratio",
                         "z_score", "binom_p", "fisher_p")])

# Build the enrichment ratio plot
gene_p1_stats$Gene_P1 <- factor(gene_p1_stats$Gene_P1,
  levels = gene_p1_stats$Gene_P1[order(gene_p1_stats$enrichment_ratio)])

# Highlight NBN
gene_p1_stats$is_NBN <- gene_p1_stats$Gene_P1 == "NBN"

# Baseline annotation label (positioned at top gene level on the discrete axis)
baseline_label_gene <- levels(gene_p1_stats$Gene_P1)[nrow(gene_p1_stats)]

p_ratio <- ggplot(gene_p1_stats,
    aes(x = enrichment_ratio, y = Gene_P1)) +
  # Baseline reference
  geom_vline(xintercept = baseline_ratio, linetype = "dashed",
             color = "gray40", linewidth = 0.6) +
  annotate("text",
           x = baseline_ratio + 0.02, y = baseline_label_gene,
           label = paste0("Baseline\n(", n_enriched_total, "/",
                          n_sig_total, " = ",
                          round(baseline_ratio * 100), "%)"),
           size = 3, color = "gray40", hjust = 0, vjust = 0) +
  # Points
  geom_segment(aes(x = baseline_ratio, xend = enrichment_ratio,
                   y = Gene_P1, yend = Gene_P1,
                   color = is_NBN),
               linewidth = 0.8, show.legend = FALSE) +
  geom_point(aes(fill = is_NBN), size = 4, shape = 21,
             color = "black", stroke = 0.5) +
  # Ratio labels (count only — no stars; see permutation test in console)
  geom_text(aes(label = paste0(n_enriched, "/", n_total)),
            hjust = -0.3, size = 3.2, fontface = "plain") +
  # Scales
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "#D62728"),
                    guide = "none") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#D62728")) +
  scale_x_continuous(
    labels = scales::percent_format(),
    limits = c(0.1, max(gene_p1_stats$enrichment_ratio) + 0.12),
    breaks = seq(0.1, 0.8, 0.1)
  ) +
  labs(
    title = "Enrichment bias among position 1 genes",
    subtitle = "Fraction of significant PAIRs enriched in BFP+ cells",
    x = "Enrichment ratio (enriched / total significant)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11)
  )

# Note: NBN is visually distinguished by the red point and segment.
# Vectorized face/color in element_text() is not officially supported in ggplot2
# and will produce a warning, so we keep the y-axis styling uniform here.

# Export enrichment bias statistics
stats_export <- gene_p1_stats[, c("Gene_P1", "n_total", "n_enriched", "n_depleted",
                                   "enrichment_ratio", "z_score",
                                   "binom_p", "binom_padj",
                                   "fisher_p", "fisher_padj")]
# Add permutation test result as metadata columns
stats_export$perm_max_z_obs  <- obs_max_z
stats_export$perm_max_z_gene <- as.character(obs_max_z_gene)
stats_export$perm_p_value    <- perm_p_max_z
stats_export$perm_n_iter     <- n_perm

write.csv(stats_export,
          "./Results/NBN_enrichment_bias_statistics.csv", row.names = FALSE)

print(p_ratio)
ggsave("./Figures/NBN_Enrichment_Bias.pdf", plot = p_ratio,
       width = 6, height = 5)

# =============================================================================
# SECTION 3: Panel B — Waterfall plot of NBN position 1 PAIRs
# =============================================================================

nbn_p1 <- deseq2_sig %>%
  filter(Gene_P1 == "NBN") %>%
  arrange(log2FoldChange)

cat("\nNBN at position 1:", nrow(nbn_p1), "significant PAIRs\n")
cat("  Enriched:", sum(nbn_p1$direction == "enriched"), "\n")
cat("  Depleted:", sum(nbn_p1$direction == "depleted"), "\n")

# Add pathway annotation for position 2 partner (using HUGO symbols)
nbn_p1 <- merge(nbn_p1, pathway_map, by.x = "Gene_P2", by.y = "Gene",
                all.x = TRUE)
nbn_p1$Pathway_group[is.na(nbn_p1$Pathway_group)] <- "Other"

# Create clean partner label: HUGO symbol (sgRNA info)
nbn_p1$partner_label <- paste0(nbn_p1$Gene_P2)

# Extract CasRx guide number for disambiguation
nbn_p1$casrx_guide <- sub(".*CasRx_", "", nbn_p1$sgRNA)
nbn_p1$full_label <- paste0(nbn_p1$Gene_P2, " (CasRx_", nbn_p1$casrx_guide, ")")

# Sort by LFC for waterfall
nbn_p1 <- nbn_p1 %>% arrange(log2FoldChange)
nbn_p1$full_label <- factor(nbn_p1$full_label,
                             levels = nbn_p1$full_label)

# Save to CSV
write.csv(nbn_p1[, c("sgRNA", "Gene_P1", "Gene_P2", "log2FoldChange",
                       "padj", "direction", "Pathway_group")],
          "./Results/NBN_position1_PAIRs.csv", row.names = FALSE)

# Build waterfall plot
p_waterfall <- ggplot(nbn_p1,
    aes(x = full_label, y = log2FoldChange, fill = Pathway_group)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  # Significance stars
  geom_text(
    aes(label = ifelse(padj < 0.001, "***",
                ifelse(padj < 0.01, "**",
                ifelse(padj < 0.05, "*", ""))),
        y = log2FoldChange + sign(log2FoldChange) * 0.8),
    size = 4, vjust = 0.5
  ) +
  scale_fill_manual(
    values = pathway_colors,
    name = "Position 2 partner\npathway"
  ) +
  scale_y_continuous(
    breaks = seq(-20, 20, 5),
    limits = c(min(nbn_p1$log2FoldChange) - 2,
               max(nbn_p1$log2FoldChange) + 2)
  ) +
  # Annotate enriched/depleted regions
  annotate("text", x = 0.6, y = max(nbn_p1$log2FoldChange) - 1,
           label = "Enriched in BFP+\n(DDR active)",
           size = 3.2, color = "gray30", hjust = 0, fontface = "italic") +
  annotate("text", x = 0.6, y = min(nbn_p1$log2FoldChange) + 1,
           label = "Depleted from BFP+",
           size = 3.2, color = "gray30", hjust = 0, fontface = "italic") +
  labs(
    title = "NBN activation (position 1): partner-specific DDR outcomes",
    subtitle = "11 DESeq2-significant PAIRs with NBN at CRISPRa position",
    x = "Position 2 partner (CasRx knockdown target)",
    y = expression(log[2]~"fold change (BFP+ / unsorted)")
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10)
  )

print(p_waterfall)
ggsave("./Figures/NBN_Waterfall_Partners.pdf", plot = p_waterfall,
       width = 7, height = 5)

# =============================================================================
# SECTION 4: Combined panel (Option C)
# =============================================================================

p_combined <- p_ratio + p_waterfall +
  plot_layout(widths = c(1, 1.4)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

ggsave("./Figures/NBN_Combined_Panel.pdf", plot = p_combined,
       width = 16, height = 6)

cat("\n--- Output files ---\n")
cat("  ./Figures/NBN_Enrichment_Bias.pdf\n")
cat("  ./Figures/NBN_Waterfall_Partners.pdf\n")
cat("  ./Figures/NBN_Combined_Panel.pdf\n")
cat("  ./Results/NBN_position1_PAIRs.csv\n")
cat("  ./Results/NBN_enrichment_bias_statistics.csv\n")

# --- Session info -------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           "./Results/06_NBN_sessionInfo.txt")

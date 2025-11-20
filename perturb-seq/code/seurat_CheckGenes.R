rm(list = ls())
set.seed(7)

library(dplyr)
library(Seurat)
library(qs)
library(ggplot2)
library(tidyr)
library(readxl)

# load the genes from excel
gene_tab <- read_excel("./info/updated 250918_information_final.xlsx", sheet = 5)
genes_to_plot <- gene_tab$gene

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")
meta_data <- scRNA_seq@meta.data

# ---- Input: your gene list ----
gene_list <- c("NBN", "XRCC6", # "KU70", 
               "TP53BP1", # "53BP1", 
               "PRKDC", # "DNPKcs", 
               "PARP1", # "PQ", 
               "MRE11") 

gene_list <- unique(c(gene_list, genes_to_plot))

# -----------------------------
# 2. Extract expression matrix
# -----------------------------
expr <- GetAssayData(scRNA_seq, slot = "data")   # log-normalized assay

# -----------------------------
# 3. Metadata
# -----------------------------
meta <- scRNA_seq@meta.data

# Create combined grouping
meta$group <- paste(meta$treatment, meta$assigned_tag, sep = "_")

# -----------------------------
# 4. Ensure genes exist
# -----------------------------
genes_present <- gene_list[gene_list %in% rownames(expr)]
if (length(genes_present) == 0) stop("None of the genes found in the object.")

# -----------------------------
# 5. Compute pct expression + mean expression
# -----------------------------
df_list <- lapply(genes_present, function(g) {
  
  tibble(
    gene = g,
    group = meta$group,
    expr = expr[g, ]
  ) %>%
    group_by(group, gene) %>%
    summarise(
      pct_exp = mean(expr > 0) * 100,   # % cells expressing gene
      avg_exp = mean(expr),            # average expression
      .groups = "drop"
    )
})

df <- bind_rows(df_list)

# -----------------------------
# 2. Extract expression matrix
# -----------------------------
expr <- GetAssayData(scRNA_seq, slot = "data")   # log-normalized assay

# -----------------------------
# 3. Metadata
# -----------------------------
meta <- scRNA_seq@meta.data

# Create combined grouping
meta$group <- paste(meta$treatment, meta$assigned_tag, sep = "_")

# -----------------------------
# 4. Ensure genes exist
# -----------------------------
genes_present <- gene_list[gene_list %in% rownames(expr)]
if (length(genes_present) == 0) stop("None of the genes found in the object.")

# -----------------------------
# 5. Compute pct expression + mean expression
# -----------------------------
df_list <- lapply(genes_present, function(g) {

  tibble(
    gene = g,
    group = meta$group,
    expr = expr[g, ]
  ) %>%
    group_by(group, gene) %>%
    summarise(
      pct_exp = mean(expr > 0) * 100,   # % cells expressing gene
      avg_exp = mean(expr),            # average expression
      .groups = "drop"
    )
})

df <- bind_rows(df_list)
df$gene <- factor(df$gene, levels = genes_present)

avg_exp_wide <- df %>%
  dplyr::select(gene, group, avg_exp) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = avg_exp
  )
write.csv(avg_exp_wide, "./res/2025_1119/genes_avg_exp_allgroups.csv", row.names = FALSE)

pct_exp_wide <- df %>%
  dplyr::select(gene, group, pct_exp) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = pct_exp
  )
write.csv(pct_exp_wide, "./res/2025_1119/genes_pct_exp_allgroups.csv", row.names = FALSE)



df$avg_exp[df$avg_exp >= 1.5] <- 1.5

chosen_groups <- c("CTRL_NT_CRISPRa_NT_CasRx",
                   "CTRL_NBN_CRISPRa_NT_CasRx",
                   "CTRL_NT_CRISPRa_53BP1_CasRx",
                   "CTRL_NBN_CRISPRa_53BP1_CasRx")

# chosen_groups <- c("CTRL_NT_CRISPRa_NT_CasRx", 
#                    "RNP_NT_CRISPRa_NT_CasRx")

df_sub <- df %>% filter(group %in% chosen_groups)
df_sub$group <- factor(df_sub$group, levels = chosen_groups)

# -----------------------------
# 6. Horizontal Dot Plot
# -----------------------------
p <- ggplot(df_sub, aes(x = gene, y = group)) +     # swapped axes
  geom_point(aes(size = pct_exp, color = avg_exp)) +
  
  scale_color_gradientn(
    colors = c("white", "#B2182B"),
    name = "Average\nExpression"
  ) +
  
  scale_size(
    range = c(1, 8),
    name = "% Expressing"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 15),
    legend.position = "top",
    plot.title = element_text(size = 16, face = "bold")
  ) +
  
  labs(
    x = "",
    y = "",
    title = ""
  )

ggsave(
  filename = "./res/2025_1120/genes_exp_NBN_53BP1.png",   # output name
  plot = p,                              # plot object
  width = 20, 
  height = 4, 
  dpi = 300
)

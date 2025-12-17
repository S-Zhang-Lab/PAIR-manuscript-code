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

# ---- Input: your gene list ----
gene_list <- c("NBN", "XRCC6", # "KU70", 
               "TP53BP1", # "53BP1", 
               "PRKDC", # "DNPKcs", 
               "POLQ", # "PQ", 
               "MRE11",
               "TP53",
               "RPA1", 
               "RPA2", 
               "RPA3") 

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


df$avg_exp[df$avg_exp >= 1.5] <- 1.5

chosen_groups <- c("CTRL_NT_CRISPRa_NT_CasRx",
                   "CTRL_NBN_CRISPRa_NT_CasRx",
                   "CTRL_NT_CRISPRa_53BP1_CasRx",
                   "CTRL_NBN_CRISPRa_53BP1_CasRx")

chosen_groups <- c("CTRL_NT_CRISPRa_NT_CasRx",
                   "CTRL_NBN_CRISPRa_NT_CasRx",
                   "CTRL_NT_CRISPRa_DNPKcs_CasRx",
                   "CTRL_NBN_CRISPRa_DNPKcs_CasRx")

chosen_groups <- c("CTRL_NT_CRISPRa_NT_CasRx",
                   "CTRL_NBN_CRISPRa_NT_CasRx",
                   "CTRL_NT_CRISPRa_KU70_CasRx",
                   "CTRL_NBN_CRISPRa_KU70_CasRx")

chosen_groups <- c("CTRL_NT_CRISPRa_NT_CasRx",
                   "CTRL_NBN_CRISPRa_NT_CasRx",
                   "CTRL_NT_CRISPRa_PQ_CasRx",
                   "CTRL_NBN_CRISPRa_PQ_CasRx")

chosen_groups <- c("RNP_NT_CRISPRa_NT_CasRx",
                   "RNP_NBN_CRISPRa_NT_CasRx",
                   "RNP_NT_CRISPRa_53BP1_CasRx",
                   "RNP_NBN_CRISPRa_53BP1_CasRx")

chosen_groups <- c("RNP_NT_CRISPRa_NT_CasRx",
                   "RNP_NBN_CRISPRa_NT_CasRx",
                   "RNP_NT_CRISPRa_DNPKcs_CasRx",
                   "RNP_NBN_CRISPRa_DNPKcs_CasRx")

chosen_groups <- c("RNP_NT_CRISPRa_NT_CasRx",
                   "RNP_NBN_CRISPRa_NT_CasRx",
                   "RNP_NT_CRISPRa_KU70_CasRx",
                   "RNP_NBN_CRISPRa_KU70_CasRx")

chosen_groups <- c("RNP_NT_CRISPRa_NT_CasRx",
                   "RNP_NBN_CRISPRa_NT_CasRx",
                   "RNP_NT_CRISPRa_PQ_CasRx",
                   "RNP_NBN_CRISPRa_PQ_CasRx")

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
  filename = "./res/2025_1217/RNP_genes_exp_NBN_PQ.png",   # output name
  plot = p,                              # plot object
  width = 10, 
  height = 4, 
  dpi = 300
)

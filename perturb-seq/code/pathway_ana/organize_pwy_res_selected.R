rm(list = ls())
set.seed(7)
library(data.table) 
library(Seurat)
library(dplyr)
library(fgsea)
library(qs)
library(ggplot2)

root_dir <- "./res/2025_1120/"

treatments <- c("CTRL", "RNP")

pwy_types <- c("base_excision_repair", "cell_cycle_checkpoint", "dsb_repair", "recombinational_repair")

NBN_53BP1_set <- c("NBN_NT", "NT_53BP1", "NBN_53BP1")
NBN_DNPKcs_set <- c("NBN_NT", "NT_DNPKcs", "NBN_DNPKcs")
NBN_KU70_set <- c("NBN_NT", "NT_KU70", "NBN_KU70")
NBN_PQ_set <- c("NBN_NT", "NT_PQ", "NBN_PQ")

chosen_set <- NBN_53BP1_set
chosen_pwy_type <- pwy_types[4]
treatment <- "RNP"

fig_save_path <- paste0(root_dir, "pwy_results_", treatment, "/", chosen_pwy_type, "/", 
                        chosen_set[3], "_dot.png")

tab_all <- c()
for (i in 1:length(chosen_set)) {
  pwy_res <- read.csv(paste0(root_dir, "pwy_results_", treatment, "/",
                             chosen_pwy_type, "/", chosen_set[i], "_over_NT_NT_c2.csv"))
  pwy_res$sources <- paste0(treatment, "_", chosen_set[i])
  tab_all <- rbind(tab_all, pwy_res)
}

plot_tab <- tab_all[, c("pathway", "NES", "padj", "sources")]
plot_tab$sources <- factor(plot_tab$sources, levels = paste0(treatment, "_", chosen_set))

# Make sure pathway is a factor (controls plotting order)
plot_tab$pathway <- factor(plot_tab$pathway, levels = rev(unique(plot_tab$pathway)))

df <- plot_tab

df2 <- df %>%
  mutate(
    sig = padj <= 0.05,
    size_plot = ifelse(sig, -log10(padj), 1),
    color_plot = NES
  )

p <- ggplot(df2, aes(x = sources, y = pathway)) +
  geom_point(aes(
    size = size_plot,
    color = color_plot
  )) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "NES",
    na.value = "gray80"      # nonsignificant
  ) +
  scale_size(range = c(2, 10), name = "-log10(padj)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90")
  ) +
  ggtitle(chosen_pwy_type)

ggsave(
  filename = fig_save_path,   # output file name
  plot = p,                           # plot object
  width = 10,                          # width in inches
  height = 6,                         # height in inches
  dpi = 300                           # publication-quality resolution
)

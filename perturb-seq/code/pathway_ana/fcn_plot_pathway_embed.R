library(data.table) 
library(dplyr)
library(fgsea)

organize_plot_dat <- function(pathway_enrich, pathway_embed) {
  rownames(pathway_enrich) <- pathway_enrich$pathway
  rownames(pathway_embed) <- pathway_embed$cell_id
  
  overlap_pathways <- intersect(rownames(pathway_enrich), rownames(pathway_embed))
  pathway_enrich <- pathway_enrich[overlap_pathways, ]
  pathway_embed <- pathway_embed[overlap_pathways, ]
  
  plot_tab <- cbind(pathway_enrich, pathway_embed)
  
  plot_tab_clean <- plot_tab %>%
    as.data.frame() %>%
    dplyr::select(padj, NES, umap_1, umap_2) %>%
    mutate(
      sig_padj = padj < 0.05,
      sig_NES  = abs(NES) > 1,
      NES_for_color_sig = ifelse(sig_padj, NES, NA),   # NA → grey
      NES_for_color_larger1 = ifelse(sig_NES, NES, NA)   # NA → grey
    )
  
  return(plot_tab_clean)
}

plot_overall <- function(plot_tab_clean,
                         save_path,
                         width=6,
                         height=5) {
  q <- ggplot(plot_tab_clean, aes(x = umap_1, y = umap_2, color = NES)) +
    geom_point(size = 1, alpha = 0.9) +
    scale_color_gradient2(
      low = "#2166AC",     # blue for negative NES
      mid = "white",       # zero
      high = "#B2182B",    # red for positive NES
      midpoint = 0,
      name = "NES"
    ) +
    theme_classic(base_size = 14) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      legend.position = "right",
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
  ggsave(save_path, plot = q, width = width, height = height, dpi = 300)
}

plot_NES_larger1 <- function(plot_tab_clean,
                             save_path,
                             width=6,
                             height=5){
  # Plot
  q <- ggplot(plot_tab_clean, aes(x = umap_1, y = umap_2)) +
    # Grey background points (|NES| <= 1)
    geom_point(
      data = subset(plot_tab_clean, !sig_NES),
      color = "grey80",
      size = 0.8,
      alpha = 0.9
    ) +
    # Colored points (NES > 1 or NES < -1)
    geom_point(
      data = subset(plot_tab_clean, sig_NES),
      aes(color = NES_for_color_larger1),
      size = 1.0,
      alpha = 0.9
    ) +
    scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    theme_classic(base_size = 14) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(
      legend.position = "right",
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
  
  # Save figure
  ggsave(
    save_path,
    plot = q,
    width  = width,
    height = height,
    dpi    = 300
  )
  
}

plot_sig <- function(plot_tab_clean,
                     save_path,
                     width=6,
                     height=5){
  # Plot
  q <- ggplot(plot_tab_clean, aes(x = umap_1, y = umap_2)) +
    # non-significant dots
    geom_point(
      data = subset(plot_tab_clean, !sig_padj),
      color = "grey80",
      size = 0.8,
      alpha = 0.9
    ) +
    # significant dots colored by NES
    geom_point(
      data = subset(plot_tab_clean, sig_padj),
      aes(color = NES_for_color_sig),
      size = 1.0,
      alpha = 0.9
    ) +
    scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    theme_classic(base_size = 14) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(
      legend.position = "right",
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
  
  # Save
  ggsave(
    save_path,
    plot = q,
    width  = width,
    height = height,
    dpi    = 300
  )
  
}



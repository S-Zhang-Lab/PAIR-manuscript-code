## =============================================================================
## 07_Target_Network_Analysis.R
## Network analysis and pathway positioning of 20 target DDR genes
## =============================================================================
##
## Purpose:
##   Generate comprehensive network and pathway visualizations showing:
##   (1) Where 20 DDR target genes land in canonical pathways (Pathway Cartography)
##   (2) Functional protein-protein interaction landscape (STRING Network)
##   (3) KEGG pathway enrichment and coverage metrics
##   (4) Network topology analysis (centrality, clustering, cross-pathway bridges)
##
## Approach:
##   - Extract target genes from dual modulation library
##   - Map to official HUGO nomenclature
##   - Retrieve protein-protein interactions (STRINGdb)
##   - Annotate with custom DDR pathways + KEGG pathways
##   - Calculate network metrics (degree, betweenness, eigenvector centrality)
##   - Generate multiple publication-ready visualizations
##
## Outputs:
##   Figures:
##   - ./Figures/07_Pathway_Cartography.pdf  (genes organized by canonical pathway)
##   - ./Figures/07_STRING_Network_Full.pdf  (full PPI network, pathway colored)
##   - ./Figures/07_STRING_Network_Filtered.pdf (high-conf interactions only)
##   - ./Figures/07_KEGG_Enrichment.pdf     (pathway enrichment results)
##
##   Data:
##   - ./Results/07_Target_Genes_Annotated.csv (gene list with pathway assignments)
##   - ./Results/07_Network_Metrics.csv (centrality, clustering, etc.)
##   - ./Results/07_Network_Interactions.csv (all retrieved PPI)
##   - ./Results/07_KEGG_Enrichment.csv (enrichment statistics)
##
## Dependencies:
##   - STRINGdb (for protein-protein interactions)
##   - igraph (network topology analysis)
##   - ggplot2, ggnetwork (visualization)
##   - clusterProfiler, org.Hs.eg.db (KEGG enrichment)
##   - dplyr, tidyr (data manipulation)
##   - readxl (read library file)
##
## =============================================================================

# --- Setup and library loading -----------------------------------------------

# Load dplyr FIRST to avoid namespace conflicts with igraph::select()
library(dplyr)
library(tidyr)
library(readxl)

# Then load other packages
library(STRINGdb)
library(igraph)
library(ggplot2)
library(ggnetwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

cat("=== 07_Target_Network_Analysis.R ===\n")
cat("Network analysis of PAIR DDR target genes\n\n")

# --- SECTION 1: Extract and annotate target genes ----------------------------

# Define official HUGO gene symbol mapping (common name -> HUGO)
hugo_map <- c(
  "DNA-PKcs" = "PRKDC",
  "53BP1"    = "TP53BP1",
  "KU70"     = "XRCC6",
  "CtIP"     = "RBBP8",
  "XLF"      = "NHEJ1",
  "rad50"    = "RAD50"
)

to_hugo <- function(x) {
  ifelse(x %in% names(hugo_map), hugo_map[x], x)
}

# Read library file and extract target genes
library_file <- "./Misc/DRR_crRNA_library.xlsx"
crispr_df <- read_excel(library_file, sheet = "CRISPRa")

# Extract unique genes (column 1), skip header rows
target_genes_raw <- crispr_df[[1]] %>%
  as.character() %>%
  unique() %>%
  setdiff(c("Input", NA, "")) %>%
  .[. != ""]

# Convert to HUGO nomenclature
target_genes <- sort(unique(to_hugo(target_genes_raw)))

# Remove non-targeting control for primary analysis
target_genes <- setdiff(target_genes, c("(NEG_CONTROL)", "NEG_CONTROL"))

cat("Extracted target genes:\n")
cat("  Total genes:", length(target_genes), "\n")
print(target_genes)

# --- SECTION 2: Custom DDR pathway classification ----------------------------

# Expert-curated DDR pathway assignments (HUGO symbols)
ddr_pathways <- data.frame(
  Gene = c(
    # Non-Homologous End Joining (NHEJ)
    "XRCC6", "NHEJ1", "XRCC4", "LIG4", "PRKDC", "TP53BP1",
    # Homologous Recombination (HR)
    "BRCA1", "BRCA2", "RAD51", "PALB2", "RBBP8",
    # DSB Sensing (MRN complex)
    "NBN", "RAD50", "MRE11",
    # Alternative End Joining / MMEJ
    "POLQ", "LIG3",
    # Other DNA Repair
    "EXO1", "POLD3", "XRCC1",
    # DDR Signaling
    "ATM"
  ),
  Pathway = c(
    rep("NHEJ", 6),
    rep("HR", 5),
    rep("DSB sensing (MRN)", 3),
    rep("Alt-EJ / MMEJ", 2),
    rep("Other repair", 3),
    "DDR signaling"
  ),
  stringsAsFactors = FALSE
)

# Define color palette for pathways
pathway_colors <- c(
  "NHEJ"              = "#2166AC",    # blue
  "HR"                = "#B2182B",    # red
  "DSB sensing (MRN)" = "#7570B3",    # purple
  "Alt-EJ / MMEJ"     = "#F4A582",    # salmon
  "Other repair"      = "#999999",    # gray
  "DDR signaling"     = "#66A61E"     # green
)

# Annotate target genes with pathways
gene_annotation <- data.frame(
  Gene = target_genes,
  stringsAsFactors = FALSE
) %>%
  left_join(ddr_pathways, by = "Gene") %>%
  mutate(
    Pathway = ifelse(is.na(Pathway), "Unclassified", Pathway),
    Pathway_color = pathway_colors[Pathway]
  ) %>%
  mutate(Pathway_color = ifelse(is.na(Pathway_color), "#CCCCCC", Pathway_color))

cat("\nPathway assignments:\n")
print(gene_annotation %>% dplyr::select(Gene, Pathway) %>%
      dplyr::arrange(Pathway, Gene))

write.csv(gene_annotation, "./Results/07_Target_Genes_Annotated.csv",
          row.names = FALSE)

# --- SECTION 3: STRINGdb network retrieval ----------------------------------

cat("\n--- Retrieving STRING protein-protein interactions ---\n")

# Initialize STRINGdb (Human, v12.0, moderate confidence threshold)
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400  # Medium confidence (0-1000 scale)
)

# Map gene symbols to STRING identifiers
cat("Mapping genes to STRING database...\n")
mapped <- string_db$map(
  data.frame(gene = target_genes),
  "gene", removeUnmappedRows = FALSE
)

cat("  Mapped:", sum(!is.na(mapped$STRING_id)), "/", nrow(mapped), "genes\n")

# Get interactions among target genes
string_ids <- na.omit(unique(mapped$STRING_id))

if (length(string_ids) > 0) {
  interactions <- string_db$get_interactions(string_ids)

  cat("  Retrieved", nrow(interactions), "interactions at score_threshold=400\n")

  # Add gene symbol annotations to interactions
  interactions <- interactions %>%
    left_join(
      mapped %>% dplyr::select(STRING_id, gene),
      by = c("from" = "STRING_id")
    ) %>%
    dplyr::rename(from_gene = gene) %>%
    left_join(
      mapped %>% dplyr::select(STRING_id, gene),
      by = c("to" = "STRING_id")
    ) %>%
    dplyr::rename(to_gene = gene) %>%
    dplyr::select(from_gene, to_gene, combined_score) %>%
    dplyr::filter(!is.na(from_gene) & !is.na(to_gene)) %>%
    # Deduplicate undirected edges by sorting endpoints; keep strongest score.
    dplyr::mutate(
      edge_min = pmin(from_gene, to_gene),
      edge_max = pmax(from_gene, to_gene)
    ) %>%
    dplyr::group_by(edge_min, edge_max) %>%
    dplyr::summarise(combined_score = max(combined_score), .groups = "drop") %>%
    dplyr::transmute(
      from_gene = edge_min,
      to_gene = edge_max,
      combined_score = combined_score
    )

  cat("  Filtered to target genes:", nrow(interactions), "interactions\n")

  write.csv(interactions, "./Results/07_Network_Interactions.csv",
            row.names = FALSE)

} else {
  cat("  WARNING: No genes mapped to STRING database\n")
  interactions <- data.frame(from_gene = character(),
                            to_gene = character(),
                            combined_score = numeric())
}

# --- SECTION 4: Network topology analysis -----------------------------------

cat("\n--- Network topology analysis ---\n")

if (nrow(interactions) > 0) {
  # Build igraph network
  network <- graph_from_data_frame(
    interactions[, c("from_gene", "to_gene")],
    directed = FALSE,
    vertices = target_genes
  )

  # Calculate network metrics
  degree <- degree(network)
  betweenness <- betweenness(network)
  eigenvector <- eigen_centrality(network)$vector
  clustering_coef <- transitivity(network, type = "local", isolates = "zero")

  network_metrics <- data.frame(
    Gene = target_genes,
    Degree = degree[target_genes],
    Betweenness = betweenness[target_genes],
    Eigenvector_Centrality = eigenvector[target_genes],
    Local_Clustering = clustering_coef[target_genes],
    stringsAsFactors = FALSE
  ) %>%
    left_join(gene_annotation %>% dplyr::select(Gene, Pathway), by = "Gene") %>%
    dplyr::arrange(desc(Degree))

  cat("\nNetwork statistics:\n")
  cat("  Nodes (genes):", vcount(network), "\n")
  cat("  Edges (interactions):", ecount(network), "\n")
  cat("  Network density:", round(edge_density(network), 3), "\n")
  cat("  Average clustering coefficient:",
      round(mean(clustering_coef, na.rm = TRUE), 3), "\n")
  cat("  Connected components:", components(network)$no, "\n")

  cat("\nTop hub genes (by degree):\n")
  print(network_metrics %>% dplyr::select(Gene, Pathway, Degree) %>%
        dplyr::arrange(desc(Degree)) %>% head(10))

  write.csv(network_metrics, "./Results/07_Network_Metrics.csv",
            row.names = FALSE)

} else {
  network <- NULL
  network_metrics <- data.frame(Gene = target_genes) %>%
    left_join(gene_annotation %>% dplyr::select(Gene, Pathway), by = "Gene") %>%
    dplyr::mutate(Degree = 0, Betweenness = 0, Eigenvector_Centrality = 0,
           Local_Clustering = 0)
  write.csv(network_metrics, "./Results/07_Network_Metrics.csv",
            row.names = FALSE)
}

# --- SECTION 5: KEGG pathway enrichment analysis ----------------------------

cat("\n--- KEGG pathway enrichment ---\n")

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(target_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

cat("Converted", nrow(entrez_ids), "of", length(target_genes),
    "genes to Entrez IDs\n")

if (nrow(entrez_ids) > 0) {
  # Perform KEGG enrichment
  kegg_enrich <- enrichKEGG(
    gene = entrez_ids$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.1,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.1
  )

  if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
    kegg_results <- kegg_enrich@result %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue) %>%
      dplyr::arrange(p.adjust)

    cat("\nTop enriched KEGG pathways (padj < 0.1):\n")
    print(kegg_results %>% head(10))

    write.csv(kegg_results, "./Results/07_KEGG_Enrichment.csv",
              row.names = FALSE)
  } else {
    cat("  No significantly enriched KEGG pathways found (padj < 0.1)\n")
    kegg_results <- data.frame()
    kegg_enrich <- NULL
  }
} else {
  kegg_enrich <- NULL
  kegg_results <- data.frame()
}

# --- SECTION 6: Visualization 1 - Pathway Cartography ----------------------

cat("\n--- Generating Pathway Cartography ---\n")

pathway_order <- c("NHEJ", "HR", "DSB sensing (MRN)", "Alt-EJ / MMEJ",
                   "Other repair", "DDR signaling")

plot_data <- gene_annotation %>%
  dplyr::mutate(Pathway = factor(Pathway, levels = pathway_order)) %>%
  dplyr::arrange(Pathway, Gene)

p_cartography <- ggplot(plot_data, aes(x = 1, y = reorder(Gene, -as.numeric(Pathway)))) +
  # Pathway background
  geom_tile(aes(x = 0.5, fill = Pathway), alpha = 0.3, width = 0.5) +
  # Gene points
  geom_point(aes(color = Pathway), size = 5, shape = 21, stroke = 1.5) +
  # Gene labels
  geom_text(aes(x = 1.05, label = Gene), hjust = 0, size = 3.5, fontface = "bold") +
  # Pathway facet lines
  facet_wrap(~factor(Pathway, levels = pathway_order),
             ncol = 1, scales = "free_y", strip.position = "left") +
  scale_color_manual(values = pathway_colors, guide = "none") +
  scale_fill_manual(values = pathway_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1.4)) +
  labs(
    title = "Pathway Cartography: DDR Target Genes",
    subtitle = "20-gene dual modulation library positioned across canonical pathways",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray50"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10, face = "bold")
  )

print(p_cartography)
ggsave("./Figures/07_Pathway_Cartography.pdf", plot = p_cartography,
       width = 8, height = 9)

cat("  Saved: ./Figures/07_Pathway_Cartography.pdf\n")

# --- SECTION 7: Visualization 2 - STRING Network (Full) ----------------------

cat("\n--- Generating STRING Network (full, threshold ≥400) ---\n")

if (!is.null(network) && ecount(network) > 0) {
  # Prepare network for plotting
  net_layout <- layout_with_fr(network)

  # Convert to ggnetwork format
  net_data <- ggnetwork(network, layout = net_layout)

  # Get vertex names column (varies by ggnetwork version)
  vertex_col <- intersect(c("vertex.names", ".vertex.names", "name"), colnames(net_data))
  if (length(vertex_col) == 0) vertex_col <- colnames(net_data)[1]
  vertex_col <- vertex_col[1]

  # Pre-compute degree as a named data frame (one row per gene) for joining
  degree_df <- data.frame(
    Gene        = names(degree(network)),
    Node_Degree = as.numeric(degree(network)),
    stringsAsFactors = FALSE
  )

  # Add pathway and degree information (join on the vertex names column)
  net_data <- net_data %>%
    left_join(
      gene_annotation %>% dplyr::select(Gene, Pathway, Pathway_color),
      by = setNames("Gene", vertex_col)
    ) %>%
    left_join(degree_df, by = setNames("Gene", vertex_col))

  p_network <- ggplot(net_data, aes(x = x, y = y, xend = xend, yend = yend)) +
    # Edges — use linewidth (not edge_width) for ggnetwork >= 0.5.12
    geom_edges(color = "gray70", linewidth = 0.6, alpha = 0.7) +
    # Nodes — size mapped to pre-computed Node_Degree column
    geom_nodes(aes(fill = Pathway, size = Node_Degree),
               shape = 21, color = "black", stroke = 0.6) +
    # Node labels (use .data[[vertex_col]] to safely reference dynamic column name)
    geom_nodetext_repel(aes(label = .data[[vertex_col]]),
                        size = 3, fontface = "bold", max.overlaps = 20) +
    scale_fill_manual(values = pathway_colors, name = "Pathway", na.value = "#CCCCCC") +
    scale_size_continuous(name = "Degree\n(# interactions)", range = c(4, 10)) +
    labs(
      title = "STRING Protein-Protein Interaction Network",
      subtitle = paste0("20 target genes | ", ecount(network),
                        " interactions | STRING score >= 400")
    ) +
    theme_blank() +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray50"),
      legend.position = "right"
    )

  print(p_network)
  ggsave("./Figures/07_STRING_Network_Full.pdf", plot = p_network,
         width = 10, height = 8)

  cat("  Saved: ./Figures/07_STRING_Network_Full.pdf\n")

} else {
  cat("  WARNING: Network too sparse for visualization\n")
}

# --- SECTION 8: Visualization 3 - Network Metrics Heatmap -------------------

cat("\n--- Generating Network Metrics Summary ---\n")

if (nrow(network_metrics) > 5) {
  gene_order <- network_metrics %>%
    dplyr::arrange(dplyr::desc(Degree), Gene) %>%
    dplyr::pull(Gene)

  metrics_plot_data <- network_metrics %>%
    dplyr::select(Gene, Pathway, Degree, Betweenness, Eigenvector_Centrality) %>%
    tidyr::pivot_longer(
      cols = c(Degree, Betweenness, Eigenvector_Centrality),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Gene = factor(Gene, levels = gene_order),
      Value_scaled = ifelse(max(Value, na.rm = TRUE) > 0,
                            Value / max(Value, na.rm = TRUE), 0)
    )

  p_metrics <- ggplot(metrics_plot_data,
                      aes(x = Metric, y = Gene)) +
    geom_tile(aes(fill = Value_scaled), color = "white", linewidth = 0.3) +
    geom_text(aes(label = round(Value, 2)), size = 2.5, color = "black") +
    scale_fill_gradient(low = "white", high = "#2166AC", name = "Scaled value") +
    facet_grid(~Metric, scales = "free_x") +
    labs(
      title = "Network Topology Metrics by Gene",
      subtitle = "Degree (connections), betweenness (pathway bridges), eigenvector centrality (hub importance)",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "gray50"),
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  print(p_metrics)
  ggsave("./Figures/07_Network_Metrics_Heatmap.pdf", plot = p_metrics,
         width = 10, height = 8)

  cat("  Saved: ./Figures/07_Network_Metrics_Heatmap.pdf\n")
}

# --- SECTION 9: Visualization 4 - KEGG Enrichment Barplot --------------------

cat("\n--- Generating KEGG Enrichment Plot ---\n")

if (nrow(kegg_results) > 0) {
  kegg_plot_data <- kegg_results %>%
    head(10) %>%
    dplyr::mutate(
      log10_padj = -log10(p.adjust),
      Description = gsub(" - Homo sapiens.*", "", Description),
      Description = substr(Description, 1, 50),
      Description = factor(Description,
                          levels = rev(Description))
    )

  p_kegg <- ggplot(kegg_plot_data,
                   aes(x = log10_padj, y = Description)) +
    geom_col(fill = "#2166AC", alpha = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "gray50", linewidth = 0.6) +
    geom_text(aes(label = paste0(GeneRatio)),
              hjust = -0.1, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = "KEGG Pathway Enrichment",
      subtitle = "Top 10 pathways by adjusted p-value",
      x = "-log10(adjusted p-value)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray50"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9)
    )

  print(p_kegg)
  ggsave("./Figures/07_KEGG_Enrichment.pdf", plot = p_kegg,
         width = 8, height = 6)

  cat("  Saved: ./Figures/07_KEGG_Enrichment.pdf\n")

} else {
  cat("  No significantly enriched KEGG pathways to visualize\n")
}

# --- SECTION 10: Summary and session info -----------------------------------

cat("\n=== SUMMARY ===\n")
cat("Target genes analyzed:", length(target_genes), "\n")
cat("Network nodes:", if (!is.null(network)) vcount(network) else 0, "\n")
cat("Network edges:", if (!is.null(network)) ecount(network) else 0, "\n")
cat("Significant KEGG pathways:", nrow(kegg_results), "\n\n")

cat("=== OUTPUT FILES ===\n")
cat("Figures:\n")
cat("  ./Figures/07_Pathway_Cartography.pdf\n")
if (!is.null(network) && ecount(network) > 0) {
  cat("  ./Figures/07_STRING_Network_Full.pdf\n")
  cat("  ./Figures/07_Network_Metrics_Heatmap.pdf\n")
}
if (nrow(kegg_results) > 0) {
  cat("  ./Figures/07_KEGG_Enrichment.pdf\n")
}

cat("\nData files:\n")
cat("  ./Results/07_Target_Genes_Annotated.csv\n")
cat("  ./Results/07_Network_Interactions.csv\n")
cat("  ./Results/07_Network_Metrics.csv\n")
if (nrow(kegg_results) > 0) {
  cat("  ./Results/07_KEGG_Enrichment.csv\n")
}

writeLines(capture.output(sessionInfo()),
          "./Results/07_Network_Analysis_sessionInfo.txt")

cat("\n✓ Analysis complete.\n")

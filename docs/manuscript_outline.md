# Figure → source-code map

Maps every panel of every code-derived figure in the manuscript to the script
that produces it and the output file you'll find it under. Schematic /
wet-lab panels (no code) are listed for completeness as "_schematic_" or
"_wet-lab_".

## Figure 2 — DDR combinatorial CRISPR screen

Source: [`ddr-screen/`](../ddr-screen/)

| Panel | Description | Source code | Output |
|---|---|---|---|
| 2a | PAIR-seq screening strategy | _schematic_ | — |
| 2b | KEGG Pathway Enrichment (top 10) | [`Codes/07_Target_Network_Analysis.R`](../ddr-screen/Codes/07_Target_Network_Analysis.R) | [`Figures/07_KEGG_Enrichment.pdf`](../ddr-screen/Figures/07_KEGG_Enrichment.pdf), [`Results/07_KEGG_Enrichment.csv`](../ddr-screen/Results/07_KEGG_Enrichment.csv) |
| 2c | STRING protein–protein interaction network | [`Codes/07_Target_Network_Analysis.R`](../ddr-screen/Codes/07_Target_Network_Analysis.R) | [`Figures/07_STRING_Network_Full.pdf`](../ddr-screen/Figures/07_STRING_Network_Full.pdf), [`Results/07_Network_Interactions.csv`](../ddr-screen/Results/07_Network_Interactions.csv), [`Results/07_Network_Metrics.csv`](../ddr-screen/Results/07_Network_Metrics.csv) |
| 2d | Heatmap of significant DE PAIRs (DESeq2) | [`Codes/03_DE-seq2.R`](../ddr-screen/Codes/03_DE-seq2.R) | [`Figures/DESeq2_Primary_Heatmap.pdf`](../ddr-screen/Figures/DESeq2_Primary_Heatmap.pdf), [`Results/DESeq2_primary_significant_PAIRs.csv`](../ddr-screen/Results/DESeq2_primary_significant_PAIRs.csv) |
| 2e | Volcano plot of significantly changed PAIRs | [`Codes/03_DE-seq2.R`](../ddr-screen/Codes/03_DE-seq2.R) | [`Figures/DEseq2_Volcano_Plot_with_Batch.pdf`](../ddr-screen/Figures/DEseq2_Volcano_Plot_with_Batch.pdf), [`Results/DEseq2_all_results.csv`](../ddr-screen/Results/DEseq2_all_results.csv) |
| 2f | Enrichment bias among position-1 genes | [`Codes/05_integration_analysis.R`](../ddr-screen/Codes/05_integration_analysis.R) | [`Results/gene_frequency_enriched.csv`](../ddr-screen/Results/gene_frequency_enriched.csv), [`Results/NBN_enrichment_bias_statistics.csv`](../ddr-screen/Results/NBN_enrichment_bias_statistics.csv) |
| 2g | NBN (position 1) partner enrichments waterfall | [`Codes/06_NBN_focus_analysis.R`](../ddr-screen/Codes/06_NBN_focus_analysis.R) | [`Figures/NBN_Waterfall_Partners.pdf`](../ddr-screen/Figures/NBN_Waterfall_Partners.pdf), [`Figures/NBN_Combined_Panel.pdf`](../ddr-screen/Figures/NBN_Combined_Panel.pdf), [`Results/NBN_position1_PAIRs.csv`](../ddr-screen/Results/NBN_position1_PAIRs.csv) |

## Figure 4 — Bidirectional single PAIR-perturb-seq

Source: [`perturb-seq/`](../perturb-seq/)

| Panel | Description | Source code | Output |
|---|---|---|---|
| 4a | scPAIR-seq workflow | _schematic_ | — |
| 4b | TapeStation library traces | _wet-lab_ | — |
| 4c | Target-gene expression dotplot by condition | [`code/analysis_code/01_perturbation_validation.R`](../perturb-seq/code/analysis_code/01_perturbation_validation.R) | `figures/01_validation/` |
| 4d | NBN's effect (ΔNBN) on DSB scatterplot | [`code/analysis_code/03_tier2_stress_interaction.R`](../perturb-seq/code/analysis_code/03_tier2_stress_interaction.R) | `figures/03_tier2/`, `output/03_tier2/` |
| 4e | Hallmark GSEA NES heatmap (4 comparisons) | [`code/analysis_code/05_fgsea_analysis.R`](../perturb-seq/code/analysis_code/05_fgsea_analysis.R) | `figures/05_fgsea/`, `output/05_fgsea/` |
| 4f | Down-regulated Hallmark pathways (effect size + significance) | [`code/analysis_code/05_fgsea_analysis.R`](../perturb-seq/code/analysis_code/05_fgsea_analysis.R) | `figures/05_fgsea/` |
| 4g | AUCell DDR module activity (HR / NHEJ / MMEJ) | [`code/analysis_code/07_tier5_aucell.R`](../perturb-seq/code/analysis_code/07_tier5_aucell.R) | `figures/07_tier5/`, `output/07_tier5/` |
| 4h | Candidate Rescue by Partner KD (% genes rescued) | [`code/analysis_code/08_phase2_epistasis.R`](../perturb-seq/code/analysis_code/08_phase2_epistasis.R) | `figures/08_phase2/`, `output/08_phase2/` |
| 4i | Number of DE genes per dual-perturbation arm | [`code/analysis_code/09_cross_partner_overlap.R`](../perturb-seq/code/analysis_code/09_cross_partner_overlap.R) | `figures/09_cross_partner_overlap/`, `output/09_cross_partner_overlap/` |

## Extended Data Fig 5 — perturb-seq supporting analyses

Source: [`perturb-seq/`](../perturb-seq/)

| Panel | Description | Source code | Output |
|---|---|---|---|
| 5a, 5b | scPAIR-seq library schematic + sequence | _schematic_ | — |
| 5c | UMAP coloured by DDR-PAIR assignment | [`code/analysis_code/00_qc_filtering.R`](../perturb-seq/code/analysis_code/00_qc_filtering.R) | `figures/00_qc/` |
| 5d | UMAP coloured by RNP treatment condition | [`code/analysis_code/00_qc_filtering.R`](../perturb-seq/code/analysis_code/00_qc_filtering.R) | `figures/00_qc/` |
| 5e | Pseudobulk PCA manifold | [`code/analysis_code/00_qc_filtering.R`](../perturb-seq/code/analysis_code/00_qc_filtering.R) | `figures/00_qc/` |
| 5f | Volcano: NBN_CRISPRa baseline (CTRL only) | [`code/analysis_code/02_tier1_nbn_baseline.R`](../perturb-seq/code/analysis_code/02_tier1_nbn_baseline.R) | `figures/02_tier1/`, `output/02_tier1/` |
| 5g | Volcanos: TP53BP1 / XRCC6 / POLQ KD vs NT (RNP, NBN) | [`code/analysis_code/04_tier3_partner_modulation.R`](../perturb-seq/code/analysis_code/04_tier3_partner_modulation.R) | `figures/04_tier3/`, `output/04_tier3/` |
| 5h | RNP-context residual scatter per partner | [`code/analysis_code/08_phase2_epistasis.R`](../perturb-seq/code/analysis_code/08_phase2_epistasis.R) | `figures/08_phase2/` |
| 5i | AUCell mean AUC per condition × treatment | [`code/analysis_code/07_tier5_aucell.R`](../perturb-seq/code/analysis_code/07_tier5_aucell.R) | `figures/07_tier5/` |
| 5j | Jaccard similarity of partner-KD DE gene sets | [`code/analysis_code/09_cross_partner_overlap.R`](../perturb-seq/code/analysis_code/09_cross_partner_overlap.R) | `figures/09_cross_partner_overlap/` |

See also [`perturb-seq/docs/PAIR_ASSIGNMENT_POLICY.md`](../perturb-seq/docs/PAIR_ASSIGNMENT_POLICY.md) for the rationale behind the per-cell PAIR identity assignment used in panels 5c-5d.

## Extended Data Fig 2 (related to Fig 1)

Wet-lab characterization (FACS gating, RT-qPCR, cloning schematics). No
analysis code is shipped for these panels.

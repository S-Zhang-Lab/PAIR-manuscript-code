# Manuscript outline

Map of each figure to the source code that generates it. Final figure numbering
will be confirmed at publication.

## Figure: DDR combinatorial screen

Source: [`ddr-screen/`](../ddr-screen/)

| Panel | Source code | Output |
|---|---|---|
| Volcano | [`ddr-screen/Codes/03_DE-seq2.R`](../ddr-screen/Codes/03_DE-seq2.R) | `ddr-screen/Figures/DEseq2_Volcano_Plot_with_Batch.pdf` |
| Heatmap (significant PAIRs) | [`ddr-screen/Codes/03_DE-seq2.R`](../ddr-screen/Codes/03_DE-seq2.R) | `ddr-screen/Figures/DESeq2_Primary_Heatmap.pdf` |
| edgeR cross-validation | [`ddr-screen/Codes/04_EdgeR.R`](../ddr-screen/Codes/04_EdgeR.R) | `ddr-screen/Figures/edgeR_Volcano_Plot_with_Batch.pdf` |
| NBN enrichment + waterfall | [`ddr-screen/Codes/06_NBN_focus_analysis.R`](../ddr-screen/Codes/06_NBN_focus_analysis.R) | `ddr-screen/Figures/NBN_Combined_Panel.pdf` |
| Network / pathway cartography | [`ddr-screen/Codes/07_Target_Network_Analysis.R`](../ddr-screen/Codes/07_Target_Network_Analysis.R) | `ddr-screen/Figures/07_*.pdf` |

See also [`ddr-screen/Docs/Figure_2_descriptions.md`](../ddr-screen/Docs/Figure_2_descriptions.md).

## Figures: Perturb-seq mechanistic follow-up

Source: [`perturb-seq/`](../perturb-seq/)

| Panel | Source code | Output |
|---|---|---|
| QC / sample composition | [`perturb-seq/code/analysis_code/00_qc_filtering.R`](../perturb-seq/code/analysis_code/00_qc_filtering.R) | `figures/00_qc/` |
| NBN baseline DE | [`perturb-seq/code/analysis_code/02_tier1_nbn_baseline.R`](../perturb-seq/code/analysis_code/02_tier1_nbn_baseline.R) | `output/02_tier1/` |
| Partner attenuation (IFNα/γ, MYC) | [`perturb-seq/code/analysis_code/05_fgsea_analysis.R`](../perturb-seq/code/analysis_code/05_fgsea_analysis.R) | `output/05_fgsea/` |
| Robustness (downsampling + low-count) | scripts 14 / 15 / 16 | `output/15_*`, `output/16_*` |

See also [`perturb-seq/docs/OBSERVATIONS_SUMMARY.md`](../perturb-seq/docs/OBSERVATIONS_SUMMARY.md)
and [`perturb-seq/docs/DATA_ANALYSIS_GUIDE.md`](../perturb-seq/docs/DATA_ANALYSIS_GUIDE.md).

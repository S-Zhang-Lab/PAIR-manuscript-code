# Manuscript outline

*Stub — fill in as the manuscript stabilizes.*

Map each figure to its source code, key result tables, and evidence tier. Use
this as the index reviewers and co-authors hit first.

## Figure 1 — TBD

| Panel | Source code | Output table | Evidence tier |
|---|---|---|---|
| 1A | _TBD_ | _TBD_ | _TBD_ |

## Figure 2 — DDR combinatorial screen

Source: [`ddr-screen/`](../ddr-screen/)

| Panel | Source code | Output | Evidence tier |
|---|---|---|---|
| Volcano | [`ddr-screen/Codes/03_DE-seq2.R`](../ddr-screen/Codes/03_DE-seq2.R) | `ddr-screen/Figures/DEseq2_Volcano_Plot_with_Batch.pdf` | L1 |
| Heatmap (217 PAIRs) | [`ddr-screen/Codes/03_DE-seq2.R`](../ddr-screen/Codes/03_DE-seq2.R) | `ddr-screen/Figures/DESeq2_Primary_Heatmap.pdf` | L1 |
| edgeR cross-validation | [`ddr-screen/Codes/04_EdgeR.R`](../ddr-screen/Codes/04_EdgeR.R) | `ddr-screen/Figures/edgeR_Volcano_Plot_with_Batch.pdf` | L2 |
| NBN enrichment + waterfall | [`ddr-screen/Codes/06_NBN_focus_analysis.R`](../ddr-screen/Codes/06_NBN_focus_analysis.R) | `ddr-screen/Figures/NBN_Combined_Panel.pdf` | L1 |
| Network / pathway cartography | [`ddr-screen/Codes/07_Target_Network_Analysis.R`](../ddr-screen/Codes/07_Target_Network_Analysis.R) | `ddr-screen/Figures/07_*.pdf` | L1 + L3 |

See also [`ddr-screen/Docs/Figure_2_descriptions.md`](../ddr-screen/Docs/Figure_2_descriptions.md).

## Figure 3+ — Perturb-seq mechanistic follow-up

Source: [`perturb-seq/`](../perturb-seq/)

| Panel | Source code | Output | Evidence tier |
|---|---|---|---|
| QC / sample composition | [`perturb-seq/code/analysis_code/00_qc_filtering.R`](../perturb-seq/code/analysis_code/00_qc_filtering.R) | `figures/00_qc/` | L1 |
| NBN baseline DE | [`perturb-seq/code/analysis_code/02_tier1_nbn_baseline.R`](../perturb-seq/code/analysis_code/02_tier1_nbn_baseline.R) | `output/02_tier1/` | L2 |
| Partner attenuation (IFNα/γ, MYC) | [`perturb-seq/code/analysis_code/05_fgsea_analysis.R`](../perturb-seq/code/analysis_code/05_fgsea_analysis.R) | `output/05_fgsea/` | L1 + L3 |
| Robustness (downsampling + low-count) | scripts 14/15/16 | `output/15_*`, `output/16_*` | L2 |

See also [`perturb-seq/docs/OBSERVATIONS_SUMMARY.md`](../perturb-seq/docs/OBSERVATIONS_SUMMARY.md)
and [`perturb-seq/docs/REVIEW_SUMMARY.md`](../perturb-seq/docs/REVIEW_SUMMARY.md).

# Code Overview — PAIR Perturb-seq Analysis

This folder contains all code used for sequencing alignment, PAIR barcode extraction, and downstream single-cell RNA-seq analysis.

---

## Folder Structure

```
code/
├── SeqAlignment/        # CellRanger alignment scripts and configuration files
├── PAIR_extraction/     # PAIR barcode extraction pipeline (Python + shell)
├── analysis_code/       # R scripts for downstream analysis (run in order)
├── organize_data.R      # Loads raw CellRanger output and assembles the initial Seurat object
├── seurat_prep.R        # Converts Ensembl IDs to HGNC symbols and prepares the Seurat object
└── prepare_pathways.R   # Downloads and saves MSigDB pathway gene sets
```

---

## SeqAlignment

Scripts for reference genome download and CellRanger multi pipeline to generate mRNA count matrices and demultiplex HTO-labeled samples.

---

## PAIR_extraction

Sequential pipeline to extract and quantify PAIR perturbation barcodes from raw sequencing reads.

| Script | Description |
|--------|-------------|
| `Step1_PAIR_extraction.sh` + `01_Extract_PAIR.py` | Extract PAIR reads from FASTQ files |
| `Step2_PAIR_Cellbc_mapping.sh` + `02_Cellbc_Whitelist_Mapping.py` | Map PAIR reads to cell barcodes |
| `Step3_PAIR_whitelist_mapping.sh` + `03_PAIR_Whitelist_Mapping.py` | Map to PAIR perturbation whitelist |
| `Step4_PAIR_UMI_mat.sh` + `04_Generate_PAIR_UMI_matrix.py` | Generate UMI count matrix |
| `Steps2to4.sh` | Convenience wrapper to run Steps 2–4 sequentially |

Reference files (`mRNA_barcodes.csv`, `PAIR_perturb_whitelist.csv`, etc.) are included in this folder.

---

## Analysis Code

R scripts for all downstream analyses. Run in numbered order. All scripts assume the working directory is the project root (`./`), with data in `./data/` and results written to `./res/`.

| Script | Description |
|--------|-------------|
| `00_check_and_install_packages.R` | Checks and installs all required R packages |
| `00_qc_filtering.R` | MAD-based QC filtering; produces `data/seu_qc.qs` |
| `01_perturbation_validation.R` | Validates CRISPRa (NBN activation) and CasRx (partner knockdown) efficiency |
| `02_tier1_nbn_baseline.R` | Differential expression: NBN-CRISPRa vs. NT control (baseline NBN effect) |
| `03_tier2_stress_interaction.R` | Interaction analysis: genes hyper- or hypo-responsive to stress in NBN-perturbed cells |
| `04_tier3_partner_modulation.R` | DE analysis for each dual perturbation (NBN + TP53BP1/XRCC6/POLQ) vs. NBN alone |
| `05_fgsea_analysis.R` | fGSEA pathway enrichment on MSigDB Hallmark gene sets for Tier 2 and Tier 3 results |
| `06_tier4_repair_modules.R` | Computes HR, NHEJ, and MMEJ module scores; characterizes repair pathway activity |
| `07_tier5_aucell.R` | AUCell-based per-cell pathway activity scoring |
| `08_phase2_epistasis.R` | Epistasis modeling: additive model residuals, pseudobulk PCA manifold, pathway tau scores |
| `09_cross_partner_overlap.R` | Cross-partner gene set overlap: UpSet plots, Jaccard similarity, per-class enrichment |
| `10_repair_profile_clustering.R` | Repair profile clustering using module scores (v1) |
| `10_repair_profile_clustering_v2.R` | Repair profile clustering using AUCell scores; per-cell ternary plots and heterogeneity analysis (v2, preferred) |
| `11_GxGxE_ctrl_epistasis.R` | GxGxE framework: epistasis in CTRL vs. RNP conditions; identifies rewired genes |
| `12_cell_heterogeneity.R` | Single-cell heterogeneity and bimodality analysis (Hartigan's dip test, synergy score UMAP) |
| `regenerate_figures.R` | Regenerates key figures from pre-computed CSVs without re-running Seurat |
| `regenerate_hugo_figures.R` | Updates figures to use standardized HUGO gene names |

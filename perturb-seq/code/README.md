# Code Overview — PAIR Perturb-seq Analysis

This folder contains all code used for sequencing alignment, PAIR barcode extraction, and downstream single-cell RNA-seq analysis.

---

## Folder Structure

```
code/
├── SeqAlignment/        # CellRanger alignment scripts and configuration files
├── PAIR_extraction/     # PAIR barcode extraction pipeline (Python + shell)
├── analysis_code/       # R scripts for downstream analysis (run in order)
│   └── pathway_embedding/   # Optional ProtTrans pathway-embedding pipeline
├── config.R             # Project configuration (sources config.local.R)
├── config.local.R       # Per-machine DATA_ROOT setting (gitignored)
├── organize_data.R      # Loads raw CellRanger output and assembles the initial Seurat object
├── seurat_prep.R        # Converts Ensembl IDs to HGNC symbols and prepares the Seurat object
└── prepare_pathways.R   # Downloads and saves MSigDB pathway gene sets
```

---

## Configuration

Every R script sources `code/config.R`, which in turn sources `code/config.local.R`
(per-machine, gitignored) to obtain `DATA_ROOT`. The following variables are
exported to the global environment:

| Variable | Role |
|---|---|
| `DATA_ROOT` | Project root — set in `config.local.R` per machine |
| `DATA_DIR` | `DATA_ROOT/data` |
| `OBJS_DIR` | `DATA_ROOT/data/objs` — Seurat objects |
| `OUTPUT_DIR` | `DATA_ROOT/output` — tabular results (CSV / RDS / MD) |
| `FIGURES_DIR` | `DATA_ROOT/figures` — plot files (PDF / PNG) |
| `RES_DIR` | `DATA_ROOT/res` — legacy, kept for backward compat |
| `SEU_RAW` | `OBJS_DIR/raw_seu_with_hto_pair_tags.qs` |
| `SEU_PREP` | `OBJS_DIR/seu_prep.qs` |
| `SEU_QC` | `OBJS_DIR/seu_qc.qs` |
| `EMBED_DIR` | `DATA_DIR/embedding` |
| `C2_EMBED` | `EMBED_DIR/c2_pathway_embeddings.qs` |
| `PROT_EMBED` | `EMBED_DIR/hs_ProtTrans_embed_All.rds` |

To set up: `cp code/config.R code/config.local.R` and edit `DATA_ROOT` in the
local copy. See `repo/README.md` → "Local setup" for full instructions.

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

R scripts for all downstream analyses. Run in numbered order from the repo root. Each
script sources `code/config.R`, reads inputs via `SEU_QC`/`OUTPUT_DIR`, and writes
tables to `OUTPUT_DIR/<step>/` and figures to `FIGURES_DIR/<step>/` (both under
`DATA_ROOT`).

| Script | Description |
|--------|-------------|
| `00_check_and_install_packages.R` | Checks and installs all required R packages |
| `00_qc_filtering.R` | MAD-based QC filtering; produces `data/objs/seu_qc.qs` |
| `01_perturbation_validation.R` | Validates CRISPRa (NBN activation) and CasRx (partner knockdown) efficiency; §6 adds PAIR assignment diagnostics (UMI margin histograms, multi-tag fraction) |
| `02_tier1_nbn_baseline.R` | Differential expression: NBN-CRISPRa vs. NT control (baseline NBN effect) |
| `03_tier2_stress_interaction.R` | Stress-response deviation scoring: genes hyper- or hypo-responsive to stress in NBN-perturbed cells (exploratory L2) |
| `04_tier3_partner_modulation.R` | Pooled cell-level DE for each dual perturbation (NBN + TP53BP1/XRCC6/POLQ) vs. NBN alone; results used as internal ranking |
| `05_fgsea_analysis.R` | fGSEA pathway enrichment on MSigDB Hallmark gene sets for Tier 2 and Tier 3 results |
| `06_tier4_repair_modules.R` | Computes HR, NHEJ, and MMEJ module scores; characterizes repair pathway activity |
| `07_tier5_aucell.R` | AUCell-based per-cell pathway activity scoring |
| `08_phase2_epistasis.R` | Descriptive combinatorial perturbation scoring: additive model residuals (Module A candidate classes), pseudobulk PCA manifold, pathway tau scores (exploratory L2) |
| `09_cross_partner_overlap.R` | Cross-partner gene set overlap: UpSet plots, Jaccard similarity, per-class enrichment |
| `10_repair_profile_clustering.R` | Repair profile clustering using module scores (v1) |
| `10_repair_profile_clustering_v2.R` | Repair profile clustering using AUCell scores; per-cell ternary plots and heterogeneity analysis (v2, preferred) |
| `11_GxGxE_ctrl_epistasis.R` | Condition-dependent residual pattern analysis (GxGxE): candidate condition-dependent residuals in CTRL vs. RNP; exploratory L2 |
| `12_cell_heterogeneity.R` | Cell heterogeneity and distribution diagnostics: Hartigan's dip test (unimodality check, low power), synergy-score UMAP |
| `13_linear_interaction_model.R` | Fold-change geometry / contribution analysis: `FC_AB ~ beta1*FC_A + beta2*FC_B`; reports NBN-weighted contribution fraction (A1 + A6) |
| `14_replicate_diagnostics.R` | **Robustness (critique §C.1/C.2):** replicate composition, per-replicate module scores, and RNP2-vs-RNP3 pathway direction sign-consistency check |
| `15_downsampling_stability.R` | **Robustness (critique §C.3):** 50-bootstrap downsampling of NBN+NT reference arm; 95% CI intervals for Hallmark pathway NES per partner |
| `16_low_count_sensitivity.R` | **Robustness (critique §C.4):** re-runs Tier 3 fGSEA at progressive cell-count cutoffs; confirms direction stability when low-count arms (POLQ, 132 cells) are excluded |
| `regenerate_figures.R` | Regenerates key figures from pre-computed CSVs without re-running Seurat |
| `regenerate_hugo_figures.R` | Updates figures to use standardized HUGO gene names |

### `pathway_embedding/` — optional ProtTrans pathway embedding pipeline

An independent analysis branch (not part of the numbered chain) that projects
MSigDB C2 pathways into a ProtTrans-gene-embedding space, clusters the
resulting pathway UMAP, then runs fgsea on each cluster against Tier 2 / Tier 3
DE ranks. Outputs land under `output/embedding/<run_id>/` and
`figures/embedding/<run_id>/`. Entry point:

```bash
Rscript code/analysis_code/pathway_embedding/run_embedding_standalone.R
```

See `docs/DATA_ANALYSIS_GUIDE.md` Part VI for the full pipeline description.

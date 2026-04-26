# PAIR-DDR

Investigation of gene-gene interactions in DNA damage response (DDR) pathways
by Programmable CRISPR Paired Sequencing (PAIR-seq). Code here generates the
panels in **Figure 2** of the manuscript. See
[`../docs/manuscript_outline.md`](../docs/manuscript_outline.md) for the
panel-by-panel script map.

## Background

### Library

Dual-regulatory crRNA cassettes targeting 20 DDR genes (10,525 PAIRs total
including non-targeting controls), packaged into lentivirus and used to
transduce 293T-PAIR cells (stably expressing Cas13d + CRISPRa).

The U6-driven PAIR RNA cassette layout:

1. `e1`: `AGGGCCTATTTCCCATGATTcgtctcacaccg`
2. `e2`: N20 — CRISPRa gRNA (position 1)
3. `e3`: scaffold linker
4. `e4`: N23 — Cas13d (CasRx) crRNA (position 2)
5. `e5`: `TTTTTTT`
6. `e6`: 14-bp barcode (Hamming distance ≥5, 2-bit error correction)
7. `e7`: `ctacagagacgcacttgtacttcagcggtc`

The DNA oligo library was designed using
[`Codes/01_Lib_assemble.R`](Codes/01_Lib_assemble.R).

### Experiment

293T-PAIR cells were transduced with the DDR PAIR library at MOI 0.2,
puromycin-selected for 7 days, then split into three dishes and transfected
with SaCas9/sgRNA + a BFP DNA donor to introduce a DSB. 4 days
post-transfection, BFP+ cells were FACS-sorted and harvested for amplicon
library construction.

**Groups:**
- Replicate 1: DDR1 (unsorted), DDR2 (BFP-sorted)
- Replicate 2: DDR3 (unsorted), DDR4 (BFP-sorted)
- DDR5: puromycin-selected reference (no DSB induction)

**Sequencing:** NovaSeq X Plus, PE150 (one lane). Date: 11/15/2024.

### Raw data

See [`../data/README.md`](../data/README.md) for raw FASTQ accessions and
processed-data download links.

## Pipeline

### FASTQ → count matrix (Python + shell)

`Codes/Fastq_processing_code/UTSW25_DDR.sh` orchestrates the full FASTQ
processing. Set `PROJECT_ROOT` to your local FASTQ download root before
running.

The pipeline: fastp QC → `UTSW25_DDR_BC_extraction.py` (extract DDR target
+ barcode) → `calculate_frequencies.py` → `01_Mapping.py` (k-mer error
correction, Hamming distance ≤1, count threshold ≥50) →
`02_assemble_count_matrix.py`.

The k-mer corrected count matrix
[`Misc/kmer20_HD1_count_matrix.csv`](Misc/kmer20_HD1_count_matrix.csv) is
already shipped in this repo, so the FASTQ steps are only needed if you
want to re-derive it from raw reads.

### Statistical analysis (R)

**DESeq2** is the primary statistical method (padj < 0.01, |log2FC| > 5,
pre-filter ≥3 zero counts). Among parametric tools available, DESeq2's
shrinkage-based dispersion estimation is most robust for the n=2 design.

**Gene-level frequency analysis** (script 05) identifies individual genes
that appear at unexpectedly high frequency among DESeq2-significant PAIRs
via a hypergeometric null model.

#### Run order

```r
# From RStudio with PAIR_DDR.Rproj open, or from the command line:
Rscript Codes/00_run_all.R
```

Or step-by-step:

| Script | Role |
|---|---|
| [`Codes/02_count_matrix_analysis.R`](Codes/02_count_matrix_analysis.R) | Exploratory CPM normalization and per-replicate log fold changes |
| [`Codes/03_DE-seq2.R`](Codes/03_DE-seq2.R) | **Primary** DESeq2 with batch correction; produces volcano + heatmap (Fig 2d, 2e) |
| [`Codes/05_integration_analysis.R`](Codes/05_integration_analysis.R) | Hit integration + gene-level hypergeometric frequency analysis (Fig 2f) |
| [`Codes/06_NBN_focus_analysis.R`](Codes/06_NBN_focus_analysis.R) | NBN-position-1 partner waterfall + enrichment-bias panels (Fig 2g) |
| [`Codes/06_final_plots.R`](Codes/06_final_plots.R) | Final visualization polish |
| [`Codes/07_Target_Network_Analysis.R`](Codes/07_Target_Network_Analysis.R) | STRINGdb network + KEGG enrichment + centrality metrics (Fig 2b, 2c) |

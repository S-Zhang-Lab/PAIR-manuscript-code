# Installation and environment setup

Tested on macOS (Apple Silicon) and Linux. The figure-generating scripts in
this repo run end-to-end once R, Python, and the listed packages are in place.

## R

- **R version:** ≥ 4.5

### Core packages (both subprojects)

```r
install.packages(c(
  "tidyverse", "ggplot2", "ggrepel", "patchwork", "data.table",
  "openxlsx", "readxl", "pheatmap"
))
```

### perturb-seq subproject

```r
install.packages(c("Seurat", "qs"))
BiocManager::install(c("fgsea", "AUCell"))
install.packages("msigdbr")
```

`Seurat` v5 is required. `qs` is used for fast object serialization (`.qs`
files in `data/objs/`).

### ddr-screen subproject

```r
BiocManager::install(c("DESeq2", "STRINGdb", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("igraph", "ggnetwork"))
```

## Python (only required for FASTQ processing)

- **Python version:** ≥ 3.9

```bash
pip install numpy pandas scipy biopython
```

The Python code is used for FASTQ-level processing (PAIR barcode extraction in
`perturb-seq/code/PAIR_extraction/` and DDR target/barcode extraction in
`ddr-screen/Codes/Fastq_processing_code/`). If you start from the pre-counted
matrices in `data/`, Python is not needed.

## Optional / system-level

- **Cell Ranger** (10x Genomics) — only needed if re-aligning raw FASTQs for
  the perturb-seq side. The repo ships the wrapper script
  `perturb-seq/code/SeqAlignment/mRNA_HTO_cellranger_count.sh`.
- **fastp** — used by the DDR FASTQ processing wrapper.

## Configuration

### perturb-seq

All scripts source `perturb-seq/code/config.R`, which sources
`perturb-seq/code/config.local.R` (per-machine, gitignored). To set up:

```bash
cp perturb-seq/code/config.example.R perturb-seq/code/config.local.R
```

Then edit `config.local.R` and set `DATA_ROOT` to the directory where you
downloaded the perturb-seq data archive (see `data/README.md`).

### ddr-screen

Scripts use paths relative to `ddr-screen/`. Open
`ddr-screen/PAIR_DDR.Rproj` in RStudio (or `cd ddr-screen` before sourcing
scripts) so the working directory is correct.

The FASTQ processing shell script `ddr-screen/Codes/Fastq_processing_code/UTSW25_DDR.sh`
uses `${PROJECT_ROOT:-.}` for paths. Set `PROJECT_ROOT` to your local FASTQ
download root before running.

The perturb-seq FASTQ-processing shell scripts under
`perturb-seq/code/PAIR_extraction/` and `perturb-seq/code/SeqAlignment/`
use the same `${PROJECT_ROOT:-.}` convention.

## Reproducing figures only

If you don't need to re-run upstream processing, the `*.csv` files in
`ddr-screen/Results/` plus the perturb-seq processed Seurat object
(`seu_qc.qs`, downloadable from the data archive) are sufficient. Each
figure-generating script lists its inputs at the top.

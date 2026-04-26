# Installation and environment setup

This repo contains analyses run on macOS (Apple Silicon) and Linux (UTSW HPC).
You should be able to reproduce the figure-generating steps on either platform
once R, Python, and the listed packages are in place.

## R

- **R version:** ≥ 4.5 (arm64 supported)

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
install.packages(c("msigdbr", "diptest"))
```

`Seurat` v5 is required. `qs` is used for fast object serialization (`.qs`
files in `data/objs/`).

### ddr-screen subproject

```r
BiocManager::install(c("DESeq2", "edgeR", "STRINGdb", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("igraph", "ggnetwork"))
```

For the MAGeCK secondary analysis (referenced in the screen but not the
primary call set), install MAGeCK separately:
<https://sourceforge.net/p/mageck/wiki/install/>.

## Python

- **Python version:** ≥ 3.9

```bash
pip install numpy pandas scipy biopython
```

The Python code is used for FASTQ-level processing (PAIR barcode extraction in
`perturb-seq/code/PAIR_extraction/` and DDR target/barcode extraction in
`ddr-screen/Codes/Fastq_processing_code/`). If you start from the pre-counted
matrices in `data/` you do not need Python.

## Optional / system-level

- **CellRanger** (10x Genomics) — only needed if re-aligning raw FASTQs for
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
scripts) so working directory is correct.

The FASTQ processing shell script `ddr-screen/Codes/Fastq_processing_code/UTSW25_DDR.sh`
contains a hardcoded UTSW HPC path used during the original run. To re-run it
yourself, edit the path at the top of the script to point to your local FASTQ
download.

## Reproducing figures only

If you don't need to re-run upstream processing, the `*.csv` and `*.qs` /
`*.rds` outputs already in the repo (or downloadable from the data archive)
are sufficient. Each figure-generating script will print which input files it
needs at the top.

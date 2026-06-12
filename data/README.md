# Data

Raw and bulky processed data are hosted externally because git is not the right
place for them. This file is the canonical pointer to where each dataset lives.

The preprint describing this dataset is on bioRxiv:
<https://www.biorxiv.org/content/10.64898/2026.05.08.722799v1>

> **Reviewers / readers:** if any link below is missing or returns a 404,
> please open an issue or contact the corresponding author. Final accession
> numbers (GEO, SRA, Zenodo) will be filled in at publication.

## perturb-seq (single-cell)

| Asset | Format | Where | Approx size |
|---|---|---|---|
| Raw FASTQ (mRNA + HTO + PAIR libraries) | `.fastq.gz` | _GEO accession TBD_ | _TBD_ |
| Cell Ranger output | standard 10x | _GEO accession TBD_ | _TBD_ |
| Final QC Seurat object (`seu_qc.qs`) | `qs` (R) | _Zenodo DOI TBD_ | 85 MB |
| Intermediate Seurat object (`seu_prep.qs`) | `qs` (R) | _Zenodo DOI TBD_ | 321 MB |
| Pathway embeddings | `rds` | _Zenodo DOI TBD_ | 134 MB |
| PAIR sparse UMI matrix | `mtx` | _Zenodo DOI TBD_ | 1.3 MB |

After download, set `DATA_ROOT` in `perturb-seq/code/config.local.R` to the
directory containing the unpacked archive. The directory layout the scripts
expect is:

```
$DATA_ROOT/
├── data/
│   ├── objs/                          (Seurat .qs files)
│   ├── PAIR_output/                   (PAIR sparse UMI matrix)
│   ├── embedding/                     (pathway embeddings)
│   ├── mRNA_raw/                      (Cell Ranger output)
│   ├── pathways/                      (MSigDB Hallmark gene sets)
│   └── reference/                     (Ensembl→HGNC mapping)
└── (output/ and figures/ are written by the analysis scripts)
```

## ddr-screen (combinatorial CRISPR)

| Asset | Format | Where | Approx size |
|---|---|---|---|
| Raw FASTQ (DDR1–DDR5) | `.fastq.gz` | _SRA / GEO accession TBD_ | _TBD_ |
| k-mer corrected count matrix | `.csv` | [`ddr-screen/Misc/kmer20_HD1_count_matrix.csv`](../ddr-screen/Misc/kmer20_HD1_count_matrix.csv) | 988 KB (in repo) |
| Oligo library design | `.csv` | [`ddr-screen/Misc/Final_oligos*.csv`](../ddr-screen/Misc/) | 12 MB (in repo) |
| DESeq2 result tables + NBN focus + network outputs | `.csv` | [`ddr-screen/Results/`](../ddr-screen/Results/) | < 5 MB (in repo) |
| Figures | `.pdf` | [`ddr-screen/Figures/`](../ddr-screen/Figures/) | 432 KB (in repo) |

The ddr-screen processed outputs are checked into the repo because they are
small and useful for reviewers. The raw FASTQs are needed only if you wish to
re-derive the count matrix from scratch.

## Replication scope

To regenerate every figure in the manuscript **without** re-running upstream
processing, you need:

- For ddr-screen: nothing extra — everything required is already in the repo.
- For perturb-seq: only `seu_qc.qs` (the final QC Seurat object) and the
  MSigDB Hallmark gene sets. These are sufficient for scripts
  `00_qc_filtering.R` through `09_cross_partner_overlap.R`.

To re-run from raw FASTQs end-to-end, you additionally need the Cell Ranger
output (perturb-seq) and the FASTQ files (both subprojects).

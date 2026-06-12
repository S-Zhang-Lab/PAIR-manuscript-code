# Data

Raw and bulky processed data are hosted externally because git is not the right
place for them. This file is the canonical pointer to where each dataset lives.

The preprint describing this dataset is on bioRxiv:
<https://www.biorxiv.org/content/10.64898/2026.05.08.722799v1>

## perturb-seq (single-cell)

### Processed data — Zenodo

All processed perturb-seq artifacts are deposited as a single Zenodo record:

> **Chen Chang & Siyuan Zhang (2026).** *Parallel Activation and Interference
> CRISPR (PAIR) with Sequencing Uncovers DNA Repair Networks Guiding Precision
> Cell Engineering — Processed Data.*
> [doi:10.5281/zenodo.20672625](https://doi.org/10.5281/zenodo.20672625)
> · canonical URL: <https://zenodo.org/records/20672625>
> · total 565.2 MB

Files in the deposit:

| File | Size | Used by |
|---|---|---|
| `seu_qc.qs` | 88.8 MB | All `0X_*.R` analysis scripts (final QC Seurat object, 2,271 cells) |
| `seu_prep.qs` | 336.3 MB | Upstream of QC; only needed if you want to re-run `00_qc_filtering.R` from scratch |
| `c2_pathway_embeddings.qs` | 43.1 MB | Pathway-embedding analyses |
| `hs_ProtTrans_embed_All.rds` | 94.9 MB | ProtTrans protein embeddings |
| `c2.rds` | 2.0 MB | MSigDB C2 pathway gene sets |
| `h.rds` | 28 KB | MSigDB Hallmark pathway gene sets (input to `05_fgsea_analysis.R`) |
| `PAIR_matched_sparse_UMI_matrix.mtx` | 60.2 KB | PAIR barcode × cell sparse UMI matrix |
| `PAIR_matched_sparse_UMI_matrix_rows.txt` | 105.9 KB | Row index for the `.mtx` |
| `PAIR_matched_sparse_UMI_matrix_columns.txt` | 246 B | Column index for the `.mtx` |

### Raw data

| Asset | Format | Where |
|---|---|---|
| Raw FASTQ (mRNA + HTO + PAIR libraries) | `.fastq.gz` | _GEO accession TBD_ |
| Cell Ranger output | standard 10x | _GEO accession TBD_ |

### Local setup

After downloading the Zenodo files, arrange them under your `DATA_ROOT`
(set in `perturb-seq/code/config.local.R`) like this:

```
$DATA_ROOT/
└── data/
    ├── objs/
    │   ├── seu_qc.qs                                 ← from Zenodo
    │   └── seu_prep.qs                               ← from Zenodo (optional)
    ├── PAIR_output/
    │   ├── PAIR_matched_sparse_UMI_matrix.mtx        ← from Zenodo
    │   ├── PAIR_matched_sparse_UMI_matrix_rows.txt   ← from Zenodo
    │   └── PAIR_matched_sparse_UMI_matrix_columns.txt ← from Zenodo
    ├── embedding/
    │   ├── c2_pathway_embeddings.qs                  ← from Zenodo
    │   ├── hs_ProtTrans_embed_All.rds                ← from Zenodo
    │   └── pathways/
    │       ├── c2.rds                                ← from Zenodo
    │       └── h.rds                                 ← from Zenodo
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
- For perturb-seq: only `seu_qc.qs` and the Hallmark gene sets (`h.rds`) from
  the Zenodo deposit. These are sufficient for scripts `00_qc_filtering.R`
  through `09_cross_partner_overlap.R`.

To re-run from raw FASTQs end-to-end, you additionally need the Cell Ranger
output (perturb-seq) and the FASTQ files (both subprojects).

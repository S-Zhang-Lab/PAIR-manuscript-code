# PAIR-manuscript-code

This repository contains exactly the analysis code that generates the figures
in the manuscript — no exploratory or robustness side-analyses.

| Subdirectory | Experiment | Figures |
|---|---|---|
| [`ddr-screen/`](ddr-screen/) | Combinatorial CRISPR screen across 20 DDR genes (10,525 dual-crRNA PAIRs) in 293T-PAIR cells, BFP-reporter readout | Fig 2 |
| [`perturb-seq/`](perturb-seq/) | Single-cell PAIR-perturb-seq of NBN-CRISPRa × partner KD (TP53BP1 / XRCC6 / POLQ) under DSB stress | Fig 4, Ext Fig 5 |

## Repository layout

```
PAIR-manuscript-code/
├── ddr-screen/         Combinatorial CRISPR screen (Fig 2)
│   ├── Codes/          R + Python pipeline (00–07) + FASTQ processing
│   ├── Results/        DESeq2 / NBN-focus / network result tables (CSV)
│   ├── Figures/        Manuscript Fig 2 panels (PDF)
│   └── Misc/           Oligo library design, k-mer corrected count matrix
├── perturb-seq/        Single-cell PAIR-perturb-seq (Fig 4, Ext Fig 5)
│   ├── code/
│   │   ├── analysis_code/  Numbered R pipeline (00–09)
│   │   ├── PAIR_extraction/  Python + shell for PAIR barcode extraction
│   │   └── SeqAlignment/   Cell Ranger alignment scripts
│   ├── docs/           PAIR assignment policy
│   └── info/           PAIR library structure schematic
├── docs/               Cross-experiment figure → source map
├── data/               Where to download raw + processed data (external)
├── INSTALL.md          R + Python environment setup
├── LICENSE             MIT
└── README.md           This file
```

## Quick start

See [`INSTALL.md`](INSTALL.md) for environment setup. In brief:

1. Install R ≥ 4.5 and the packages listed in `INSTALL.md`.
2. (Only for re-running FASTQ processing) install Python ≥ 3.9 and the
   packages listed in `INSTALL.md`.
3. Download the data — see [`data/README.md`](data/README.md) for accession
   numbers.
4. Configure paths:
   - **perturb-seq:** copy `perturb-seq/code/config.example.R` to
     `perturb-seq/code/config.local.R` and set `DATA_ROOT` to your local
     download path.
   - **ddr-screen:** scripts use relative paths; open
     `ddr-screen/PAIR_DDR.Rproj` in RStudio.
5. Run the numbered scripts in order (`00_*` first).

The figure-generating CSV/PDF outputs already in the repo (under
`ddr-screen/Results/`, `ddr-screen/Figures/`) are sufficient to inspect every
panel without re-running the upstream steps.

## Citation

A preprint describing this work is available on bioRxiv:

> Zhang Lab et al. (2026). *PAIR-seq: Programmable CRISPR Paired Sequencing for
> dissecting combinatorial DNA damage response.* bioRxiv.
> <https://www.biorxiv.org/content/10.64898/2026.05.08.722799v1>
> doi:10.64898/2026.05.08.722799

The peer-reviewed citation will be added here once the paper is published.
If you use this code or any of the analysis tables, please cite the preprint
above and (optionally) link this repository.

## Data availability

Raw sequencing data, processed Seurat objects, and other files too large for
git are hosted externally — see [`data/README.md`](data/README.md) for
download links and accessions.

## License

[MIT License](LICENSE) — code is free to reuse with attribution.

## Contact

Siyuan Zhang Lab, UT Southwestern Medical Center  
Lab GitHub: <https://github.com/S-Zhang-Lab>

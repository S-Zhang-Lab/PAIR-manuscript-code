# PAIR-manuscript-code

Code, analyses, and intermediate result tables supporting the PAIR (Programmable
CRISPR Paired Sequencing) manuscript from the **S. Zhang Lab, UT Southwestern
Medical Center**.

This repository contains two complementary experiments:

| Subdirectory | Experiment | Phenotype | Scale |
|---|---|---|---|
| [`ddr-screen/`](ddr-screen/) | Combinatorial CRISPR screen across 20 DDR genes (10,525 dual-crRNA PAIRs) in 293T-PAIR cells | BFP-reporter activation (DSB sensing + repair pathway engagement) | Bulk amplicon-seq, 5 samples |
| [`perturb-seq/`](perturb-seq/) | Single-cell perturb-seq of 11 PAIR lenti lines (NBN-CRISPRa × partner KD: TP53BP1 / XRCC6 / POLQ) under DSB stress | scRNA-seq transcriptome | 2,271 cells across 1 CTRL + 2 RNP replicates |

The DDR screen identifies which gene pairs gate DSB-sensing in a focused
combinatorial space; the perturb-seq follow-up characterizes the
transcriptional consequences of selected partner combinations under DSB stress.

## Repository layout

```
PAIR-manuscript-code/
├── ddr-screen/         Combinatorial CRISPR screen
│   ├── Codes/          Numbered R + Python pipeline (00–07) + FASTQ processing
│   ├── Results/        DESeq2/edgeR result tables, gene-frequency analysis, network outputs
│   ├── Figures/        Publication-ready PDFs (volcano, heatmap, NBN/MRE11/RBBP8 panels, network)
│   ├── Misc/           Oligo library design, k-mer count matrix, MAGeCK input
│   └── Docs/           Figure 2 descriptions, NBN focus interpretation, network README
├── perturb-seq/        Single-cell perturb-seq follow-up
│   ├── code/
│   │   ├── analysis_code/  Numbered R pipeline (00–16, including robustness 14/15/16)
│   │   ├── PAIR_extraction/  Python + shell for PAIR barcode extraction
│   │   └── SeqAlignment/   CellRanger alignment scripts
│   ├── docs/           Pipeline guide, observations summary, PAIR assignment policy
│   └── info/           PAIR library structure schematic
└── docs/               Cross-experiment documentation
```

## Quick start

See [`INSTALL.md`](INSTALL.md) for environment setup. In brief:

1. Install R ≥ 4.5 and the packages listed in `INSTALL.md`.
2. Install Python ≥ 3.9 with the packages listed in `INSTALL.md`.
3. Download the data — see [`data/README.md`](data/README.md) for accession
   numbers and host locations.
4. Configure paths:
   - **perturb-seq:** copy `perturb-seq/code/config.example.R` to
     `perturb-seq/code/config.local.R` and set `DATA_ROOT` to your local
     download path.
   - **ddr-screen:** scripts use relative paths from the `ddr-screen/`
     directory; no per-machine configuration needed.
5. Run the numbered scripts in order (`00_*` first). Each subproject's README
   has a script-by-script execution guide.

## Citation

> _Manuscript citation TBD — to be added at publication._

If you use this code or any of the analysis tables, please cite the manuscript
above and (optionally) link this repository.

## Data availability

Raw sequencing data, processed Seurat objects, pathway embeddings, and other
files too large for git are hosted externally — see
[`data/README.md`](data/README.md) for download links and accessions.

The smaller analysis outputs that this repo *does* track (e.g. DESeq2 result
tables under `ddr-screen/Results/`, k-mer corrected count matrix under
`ddr-screen/Misc/`) are sufficient to regenerate every figure in the manuscript
without re-running the upstream alignment / counting steps.

## License

[MIT License](LICENSE) — code is free to reuse with attribution.

## Contact

Siyuan Zhang Lab, UT Southwestern Medical Center  
Lab GitHub: <https://github.com/S-Zhang-Lab>

# PAIR Manuscript

Code, analyses, and documentation supporting the PAIR (Programmable CRISPR Paired
Sequencing) manuscript from the S. Zhang Lab. Two complementary experiments live
side-by-side here:

| Subdirectory | Experiment | Phenotype | Scale |
|---|---|---|---|
| [`ddr-screen/`](ddr-screen/) | Combinatorial CRISPR screen across 20 DDR genes (10,525 dual-crRNA PAIRs) in 293T-PAIR cells | BFP-reporter activation (DSB sensing + repair pathway engagement) | Bulk amplicon-seq, 5 samples |
| [`perturb-seq/`](perturb-seq/) | Single-cell perturb-seq of 11 PAIR lenti lines (NBN-CRISPRa × partner KD: TP53BP1 / XRCC6 / POLQ) under DSB stress | scRNA-seq transcriptome | 2,271 cells across 1 CTRL + 2 RNP replicates |

The DDR screen identifies which gene pairs gate DSB-sensing in a focused
combinatorial space; perturb-seq characterizes the transcriptional consequences
of selected partner combinations under DSB stress. Together they form the
manuscript's mechanistic story.

## Where to start

- **Just exploring?** Read [`ddr-screen/README.md`](ddr-screen/README.md) and
  [`perturb-seq/README.md`](perturb-seq/README.md) for each experiment's
  scope, design, and how to run scripts.
- **Manuscript-level docs:** [`docs/`](docs/) holds outline, data locations,
  and any cross-cutting writing that spans both experiments.
- **AI/Claude session?** Read [`CLAUDE.md`](CLAUDE.md) first — it carries the
  evidence-tier framework and prohibited-language rules that apply across the
  whole manuscript.

## Repo provenance

This repo was consolidated 2026-04-25 from two source repositories, each
preserved with full git history under its respective subdirectory:

- `ddr-screen/` ← `S-Zhang-Lab/PAIR-DDR` (now archived)
- `perturb-seq/` ← `S-Zhang-Lab/PAIR-perturb-seq_2025-10` (now archived)

Tags `ddr-screen/pre-consolidation` and `perturb-seq/pre-consolidation` mark
the tip of each source repo at the moment of merge. The active perturb-seq
refactor lives on branch `perturb-seq/refactor-codex-rerun`.

## Data

Raw and processed data are **not** stored in git. See
[`docs/data_locations.md`](docs/data_locations.md) for Box / HPC paths, and
each subproject's `code/config.local.R` (or equivalent) for the per-machine
setup.

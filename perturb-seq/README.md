# PAIR-perturb-seq_2025-10

This is the 2nd round of PAIR-perturb-seq experiment on DDR gene paired modulation conducted by Chen Chang.

## Summary

11 dishes of human 293T-PAIR cells were infected with 11 individual lenti-PAIR RNA (P21, P53, 54, 56, 57, 61-66). After 48 hours of infection, 10 µg/mL of puromycin was added to the culture medium for drug selection. After 24 hours of drug selection, when the untransfected control cells showed 0% survival, the 11 individual PAIR-line were mix together with similar cell density. One day before 10x experiment, split PAIR-line as two groups, each group has 3 biological replicates. Group#1 transfected RNP/sgRNA targeting AAVS1 locus as DDR stimulation group, groups#2 transfected RNP/NT as transfection only control group. Each individual samples were stained with one sepecific HTO antibody (see HTO part). Mix cell together for 10x sequencing.

## Groups and HTOs used

Two groups each of them contains 3 biological replicates:

Group #1 is RNP/sgRNA transfected group, label with TotalSeq™-C0254,C0255,C0256 separately.

Group #2 is RNP transfected control group, label with TotalSeq™-C0257,C0258,C0259 separately.

TotalSeq antibody HTO barcode index sequences can be found on BioLegend website.

## Genome

Human

## Overall experimental design

11 individual PAIR lenti infected PAIR cell lines were mix with similar density and separated as two groups with 3 biological replicates (total 6). The group #1 was transfected with RNP/sgRNA(group#1) and group#2 was transfected with RNP/NT (group#2) as control separately 24 hours prior 10x experiment. After confirming high cell viability (\>80%), cell suspensions were stained individually with a unique hashing antibody so that we could subsequently pool the samples for loading onto one lane on the 10x Chromium machine with a target cell recovery of 40,000. Later, three separate libraries were prepared consisting of those samples: transcriptome, hashtag, and PAIR RNA. The 3x libraries were then pooled at an appropriate ratio and sequenced on a NovaSeq X machine.

## Experiment date

9/24/25

## Platform

NovaSeq X Plus (PE150) on one lane

## RAW data location

Box: <https://utsw.box.com/s/msqf44ib1ozl3jw2w5skw3be7d6kinqn>

## PAIR RNA library structure

<https://benchling.com/s/seq-SQx0GzDjTDG89LR1RKFs?m=slm-wvbEzq2eRuiNKKWFKc3e>

![](info/PAIR_lib_structure.png)

## HTOs and PAIRs information

Please refer to xls info sheet [here](/info/250918_information.xlsx).

## Local setup

Raw and processed data are not tracked in this repo — they live in the Box
folder linked above. The scripts in `code/` and `code/analysis_code/` find
data via a small config file:

1. **Clone the repo.** If you keep the clone inside a cloud-synced folder
   (Box, Dropbox, iCloud), exclude `.git/` from sync to avoid packfile
   corruption.

2. **Point the code at your data.** Copy `code/config.R` to
   `code/config.local.R` and edit `DATA_ROOT` to the absolute path of your
   local data directory. `config.local.R` is gitignored so each user can
   have their own path.

   ```r
   # code/config.local.R
   DATA_ROOT <- "/path/to/PAIR-perturb-seq_2025-10_SZ"
   ```

3. **Data layout under `DATA_ROOT`.** The scripts expect:

   ```
   DATA_ROOT/
   ├── data/
   │   ├── objs/
   │   │   ├── raw_seu_with_hto_pair_tags.qs
   │   │   ├── seu_prep.qs           # produced by code/seurat_prep.R
   │   │   └── seu_qc.qs             # produced by code/analysis_code/00_qc_filtering.R
   │   ├── mRNA_raw/                 # raw 10x CellRanger mRNA output
   │   ├── HTO_raw/                  # raw 10x HTO counts
   │   ├── PAIR_output/              # PAIR sparse UMI matrices
   │   ├── infercnv/                 # CNV analysis inputs
   │   ├── pathways/
   │   │   └── hallmark_pathways.rds
   │   ├── embedding/                # ProtTrans + C2 pathway embeddings
   │   │   ├── c2_pathway_embeddings.qs
   │   │   ├── hs_ProtTrans_embed_All.rds
   │   │   └── pathways/c2.rds
   │   └── data_info/                # metadata xlsx files
   ├── output/                       # tabular results (CSV / RDS / MD), gitignored
   │   └── 00_qc/ 01_validation/ 02_tier1/ ... 13_linear_interaction_model/
   └── figures/                      # plot files (PDF / PNG), gitignored
       └── 00_qc/ 01_validation/ 02_tier1/ ... 13_linear_interaction_model/
   ```

   `output/` and `figures/` mirror each other by step name: every script writes
   tables under `output/<step>/` and plots under `figures/<step>/`.

4. **Install R packages.**
   `Rscript code/analysis_code/00_check_and_install_packages.R`
   (includes `Seurat`, `qs`, `fgsea`, `msigdbr`, `AUCell`, `UpSetR`, `ggtern`,
   `diptest` — required by scripts 07/09/10/12.)

5. **Run the pipeline.** Scripts assume you start in the repo root:

   ```bash
   cd repo
   Rscript code/analysis_code/00_qc_filtering.R          # seu_prep → seu_qc
   Rscript code/analysis_code/01_perturbation_validation.R
   Rscript code/analysis_code/02_tier1_nbn_baseline.R
   # … 03 through 13
   Rscript code/analysis_code/13_linear_interaction_model.R

   # Optional (independent branch): ProtTrans pathway embedding pipeline
   Rscript code/analysis_code/pathway_embedding/run_embedding_standalone.R
   ```

### Contributing: no hardcoded paths

New scripts **must not** hardcode user-specific paths. Always:

```r
source("code/config.R")  # or "../config.R" from code/analysis_code/
seu <- qread(SEU_QC)                         # not "/Users/me/.../seu_qc.qs"
out <- file.path(OUTPUT_DIR, "my_step")      # tables (CSV / RDS / MD)
fig <- file.path(FIGURES_DIR, "my_step")     # plots (PDF / PNG)
```

Pull requests that introduce `base_dir <- "/Users/…"`, `BASE <-
"/Users/…"`, or similar will be rejected. If you need a new derived path,
add it to `code/config.R` as a new variable. `RES_DIR` is still exported
for backward compatibility but new code should use `OUTPUT_DIR` /
`FIGURES_DIR`.

### Where to read next

- **[`ANALYSIS_WORKFLOW.md`](ANALYSIS_WORKFLOW.md)** — conceptual tier-1/2/3 analysis design.
- **[`code/README.md`](code/README.md)** — developer index of every script in `code/`.
- **[`docs/DATA_ANALYSIS_GUIDE.md`](docs/DATA_ANALYSIS_GUIDE.md)** — per-step pipeline reference (inputs, algorithms, outputs, citations).
- **[`docs/OBSERVATIONS_SUMMARY.md`](docs/OBSERVATIONS_SUMMARY.md)** — current findings grounded in the most recent pipeline run.

## =============================================================================
## 00_run_all.R
## Master script: run the complete PAIR-DDR analysis pipeline (scripts 02–05)
## =============================================================================
##
## Usage:
##   Open the PAIR_DDR.Rproj in RStudio, then:
##     source("Codes/00_run_all.R")
##
##   Or from the command line:
##     Rscript Codes/00_run_all.R
##
## Pipeline order:
##   02_count_matrix_analysis.R  →  Normalize counts, exploratory analysis,
##                                   generate MAGeCK input
##   03_DE-seq2.R                →  [Supplementary] DESeq2 differential
##                                   abundance (batch-corrected)
##   04_EdgeR.R                  →  [Supplementary] edgeR GLM differential
##                                   abundance (batch-corrected)
##   05_integration_analysis.R   →  [Primary] DESeq2-based hit integration,
##                                   cross-validation with edgeR/MAGeCK,
##                                   gene-frequency analysis with null model
##   06_NBN_focus_analysis.R     →  [Main Fig] NBN enrichment bias and
##                                   partner-specific waterfall panels
##
## Dependencies:
##   02 produces ./Results/kmer20_HD1_count_matrix_for_mageck.txt
##   03 and 04 both read that file and produce full result tables
##   05 reads DESeq2/edgeR results, plus MAGeCK output if available
##
## NOTE:
##   MAGeCK output is optional. If present in ./Results/MAGeCK/, it will be
##   included in the cross-validation. If absent, the integration script
##   runs with DESeq2 and edgeR only.
## =============================================================================

# --- Set working directory to project root ------------------------------------
# If sourced from RStudio with .Rproj open, getwd() is already the project root.
# If run via Rscript, detect the script location and set accordingly.

if (interactive()) {
  # Running in RStudio — .Rproj should set the correct wd
  cat("Working directory:", getwd(), "\n")
} else {
  # Running via Rscript — navigate to project root
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    script_dir <- dirname(normalizePath(script_path))
    project_root <- normalizePath(file.path(script_dir, ".."))
    setwd(project_root)
  }
  cat("Working directory set to:", getwd(), "\n")
}

# --- Pre-flight check: note MAGeCK status -------------------------------------
mageck_sgrna_file <- "./Results/MAGeCK/DDR_default.sgrna_summary.txt"
if (!file.exists(mageck_sgrna_file)) {
  cat("\n")
  cat("================================================================\n")
  cat("  NOTE: MAGeCK output not found\n")
  cat("================================================================\n")
  cat("  The integration script (step 4) will run with DESeq2 and\n")
  cat("  edgeR only. To include MAGeCK in cross-validation, place\n")
  cat("  output files in ./Results/MAGeCK/\n")
  cat("================================================================\n\n")
}

# --- Utility: timed sourcing --------------------------------------------------
run_step <- function(script_path, step_name) {
  cat("\n")
  cat("================================================================\n")
  cat(" ", step_name, "\n")
  cat("  Sourcing:", script_path, "\n")
  cat("================================================================\n")

  if (!file.exists(script_path)) {
    stop(paste("Script not found:", script_path))
  }

  t0 <- proc.time()
  source(script_path, local = new.env(parent = globalenv()))
  elapsed <- (proc.time() - t0)["elapsed"]

  cat("  Completed in", round(elapsed, 1), "seconds\n")
  cat("================================================================\n\n")
}

# --- Run pipeline -------------------------------------------------------------

pipeline_start <- proc.time()

run_step(
  "Codes/02_count_matrix_analysis.R",
  "Step 1/5: Count matrix normalization & MAGeCK formatting"
)

run_step(
  "Codes/03_DE-seq2.R",
  "Step 2/5: [Supplementary] DESeq2 differential abundance analysis"
)

run_step(
  "Codes/04_EdgeR.R",
  "Step 3/5: [Supplementary] edgeR differential abundance analysis"
)

run_step(
  "Codes/05_integration_analysis.R",
  "Step 4/5: [Primary] DESeq2 integration & cross-validation"
)

run_step(
  "Codes/06_NBN_focus_analysis.R",
  "Step 5/5: [Main Figure] NBN enrichment bias & partner waterfall"
)

# --- Summary ------------------------------------------------------------------

total_elapsed <- (proc.time() - pipeline_start)["elapsed"]

cat("\n")
cat("================================================================\n")
cat("  Pipeline complete\n")
cat("  Total runtime:", round(total_elapsed, 1), "seconds\n")
cat("================================================================\n")
cat("\nOutput locations:\n")
cat("  Results:  ./Results/\n")
cat("  Figures:  ./Figures/\n")
cat("  Session:  ./Results/*_sessionInfo.txt\n")

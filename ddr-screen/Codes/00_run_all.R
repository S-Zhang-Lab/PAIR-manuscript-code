## =============================================================================
## 00_run_all.R
## Master script: run the PAIR-DDR analysis pipeline end-to-end.
## =============================================================================
##
## Usage:
##   Open PAIR_DDR.Rproj in RStudio, then:
##     source("Codes/00_run_all.R")
##
##   Or from the command line:
##     Rscript Codes/00_run_all.R
##
## Pipeline order:
##   02_count_matrix_analysis.R  →  Normalize counts, exploratory analysis
##   03_DE-seq2.R                →  [Primary] DESeq2 differential abundance
##                                   (batch-corrected); produces volcano +
##                                   heatmap (Fig 2d, 2e)
##   05_integration_analysis.R   →  DESeq2-based hit integration and
##                                   gene-level frequency analysis with a
##                                   hypergeometric null model (Fig 2f)
##   06_NBN_focus_analysis.R     →  NBN enrichment bias and
##                                   partner-specific waterfall panels (Fig 2g)
##   07_Target_Network_Analysis.R → STRINGdb network + KEGG enrichment
##                                   (Fig 2b, 2c)
## =============================================================================

# --- Set working directory to project root ------------------------------------

if (interactive()) {
  cat("Working directory:", getwd(), "\n")
} else {
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

run_step("Codes/02_count_matrix_analysis.R",
         "Step 1/5: Count matrix normalization & exploratory analysis")
run_step("Codes/03_DE-seq2.R",
         "Step 2/5: DESeq2 differential abundance (Fig 2d, 2e)")
run_step("Codes/05_integration_analysis.R",
         "Step 3/5: Hit integration & gene-level frequency analysis (Fig 2f)")
run_step("Codes/06_NBN_focus_analysis.R",
         "Step 4/5: NBN enrichment bias & partner waterfall (Fig 2g)")
run_step("Codes/07_Target_Network_Analysis.R",
         "Step 5/5: STRINGdb network + KEGG enrichment (Fig 2b, 2c)")

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

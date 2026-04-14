# config.R — Project data-root configuration (committed template).
#
# USAGE:
#   1. Copy this file to `config.local.R` in the same directory.
#   2. Edit `config.local.R` to set DATA_ROOT to your local path.
#   3. All analysis scripts should source this file and reference paths
#      via the derived variables below (SEU_PREP, SEU_QC, RES_DIR, ...).
#
# config.local.R is gitignored so each user/machine can have its own path
# without polluting the repo.

# Locate config.local.R. Strategy: walk the source() call stack for a frame
# whose ofile ends in "config.R" and look for "config.local.R" next to it;
# if that fails, fall back to common relative paths from the current wd.
.find_local_cfg <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    f <- try(frames[[i]]$ofile, silent = TRUE)
    if (!inherits(f, "try-error") && !is.null(f) && nzchar(f)) {
      nf <- normalizePath(f, mustWork = FALSE)
      if (endsWith(nf, "config.R") && !endsWith(nf, "config.local.R")) {
        candidate <- file.path(dirname(nf), "config.local.R")
        if (file.exists(candidate)) return(candidate)
      }
    }
  }
  # Fallbacks for common cwd locations
  for (p in c(
    "code/config.local.R",           # cwd = repo/
    "config.local.R",                # cwd = repo/code/
    "../config.local.R",             # cwd = repo/code/analysis_code/
    "../../config.local.R"           # cwd = repo/code/analysis_code/<sub>/
  )) {
    if (file.exists(p)) return(normalizePath(p, mustWork = FALSE))
  }
  NULL
}
.local_cfg <- .find_local_cfg()

if (!is.null(.local_cfg) && file.exists(.local_cfg)) {
  source(.local_cfg, local = FALSE)
} else {
  stop(
    "config.local.R not found. Searched the source call stack and the ",
    "common locations (code/config.local.R, config.local.R, ",
    "../config.local.R, ../../config.local.R).\n",
    "Create one next to code/config.R with a single line, e.g.:\n",
    '  DATA_ROOT <- "/path/to/PAIR-perturb-seq_2025-10_SZ"\n'
  )
}

if (!exists("DATA_ROOT") || is.null(DATA_ROOT) || !nzchar(DATA_ROOT)) {
  stop("DATA_ROOT is not set. Edit config.local.R to define it.")
}
if (!dir.exists(DATA_ROOT)) {
  stop("DATA_ROOT does not exist on disk: ", DATA_ROOT)
}

# Derived paths — edit DATA_ROOT (in config.local.R) to move everything.
DATA_DIR    <- file.path(DATA_ROOT, "data")
OBJS_DIR    <- file.path(DATA_DIR, "objs")
RES_DIR     <- file.path(DATA_ROOT, "res")
OUTPUT_DIR  <- file.path(DATA_ROOT, "output")   # tabular results (CSV/TSV)
FIGURES_DIR <- file.path(DATA_ROOT, "figures")  # all plot files (PDF/PNG/SVG)

SEU_RAW    <- file.path(OBJS_DIR, "raw_seu_with_hto_pair_tags.qs")
SEU_PREP   <- file.path(OBJS_DIR, "seu_prep.qs")
SEU_QC     <- file.path(OBJS_DIR, "seu_qc.qs")

EMBED_DIR  <- file.path(DATA_DIR, "embedding")
C2_EMBED   <- file.path(EMBED_DIR, "c2_pathway_embeddings.qs")
PROT_EMBED <- file.path(EMBED_DIR, "hs_ProtTrans_embed_All.rds")

rm(.find_local_cfg, .local_cfg)

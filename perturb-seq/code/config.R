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

# Locate config.local.R next to this file.
.config_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
.local_cfg <- file.path(.config_dir, "config.local.R")

if (file.exists(.local_cfg)) {
  source(.local_cfg, local = FALSE)
} else {
  stop(
    "config.local.R not found at ", .local_cfg, "\n",
    "Create it with a single line, for example:\n",
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
DATA_DIR   <- file.path(DATA_ROOT, "data")
OBJS_DIR   <- file.path(DATA_DIR, "objs")
RES_DIR    <- file.path(DATA_ROOT, "res")

SEU_RAW    <- file.path(OBJS_DIR, "raw_seu_with_hto_pair_tags.qs")
SEU_PREP   <- file.path(OBJS_DIR, "seu_prep.qs")
SEU_QC     <- file.path(OBJS_DIR, "seu_qc.qs")

EMBED_DIR  <- file.path(DATA_DIR, "embedding")
C2_EMBED   <- file.path(EMBED_DIR, "c2_pathway_embeddings.qs")
PROT_EMBED <- file.path(EMBED_DIR, "hs_ProtTrans_embed_All.rds")

rm(.config_dir, .local_cfg)

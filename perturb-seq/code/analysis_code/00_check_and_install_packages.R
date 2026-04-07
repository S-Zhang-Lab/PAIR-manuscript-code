##############################################################################
# PAIR-Perturb-Seq Analysis Pipeline — Package Environment Setup
#
# Purpose : Audit all R packages required by scripts 00–13 and
#           regenerate_*.R. Report installation status, versions,
#           and install anything missing.
#
# Usage   : Rscript code/00_check_and_install_packages.R
#           Or source interactively in RStudio.
#
# R version tested : 4.5-arm64 (Apple Silicon)
# Library path     : /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
#
# Date    : 2026-03-08
##############################################################################

# ── 0. Set library path ───────────────────────────────────────────────────────
LIB <- "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"
.libPaths(c(LIB, .libPaths()))
cat("Library path set to:", LIB, "\n\n")

# ── 1. Full package manifest ───────────────────────────────────────────────────
# Each entry: list(pkg, source, scripts_used, purpose)
#   source: "CRAN" | "Bioconductor" | "base"
#
manifest <- list(

  # ── Core data wrangling (used in almost every script) ──────────────────────
  list(pkg = "dplyr",    source = "CRAN", min_ver = "1.1.0",
       scripts = "00–13, regenerate_*.R",
       purpose = "Data manipulation: filter, mutate, summarise, join"),

  list(pkg = "tidyr",    source = "CRAN", min_ver = "1.3.0",
       scripts = "01, 05–13, regenerate_*.R",
       purpose = "Data reshaping: pivot_wider, pivot_longer, nest"),

  list(pkg = "purrr",    source = "CRAN", min_ver = "1.0.0",
       scripts = "09–13",
       purpose = "Functional programming: map, map_dfr, map2_dfr"),

  list(pkg = "stringr",  source = "CRAN", min_ver = "1.5.0",
       scripts = "09–13",
       purpose = "String manipulation: str_detect, gsub wrappers"),

  list(pkg = "tibble",   source = "CRAN", min_ver = "3.2.0",
       scripts = "11–13",
       purpose = "Modern data frames; rownames_to_column()"),

  # ── Visualization ──────────────────────────────────────────────────────────
  list(pkg = "ggplot2",  source = "CRAN", min_ver = "3.4.0",
       scripts = "00–13, regenerate_*.R",
       purpose = "Core plotting framework for all figures"),

  list(pkg = "ggrepel",  source = "CRAN", min_ver = "0.9.0",
       scripts = "02, 03, 04, 11, 12, regenerate_*.R",
       purpose = "Non-overlapping text labels on scatter plots"),

  list(pkg = "ggtern",   source = "CRAN", min_ver = "3.4.0",
       scripts = "10",
       purpose = "Ternary plots for HR/NHEJ/MMEJ repair composition (A3)"),

  list(pkg = "patchwork", source = "CRAN", min_ver = "1.1.0",
       scripts = "00, 01, 04, 06, 12, regenerate_*.R",
       purpose = "Combining multiple ggplot panels into one figure"),

  list(pkg = "pheatmap", source = "CRAN", min_ver = "1.0.12",
       scripts = "05, 08, regenerate_*.R",
       purpose = "Hierarchical clustering heatmaps"),

  list(pkg = "grid",     source = "base", min_ver = NA,
       scripts = "05, regenerate_*.R",
       purpose = "Low-level graphics (base R package, always available)"),

  # ── Genomics / Single-cell ─────────────────────────────────────────────────
  list(pkg = "Seurat",   source = "CRAN", min_ver = "5.0.0",
       scripts = "00–04, 06–08, 11, 13, regenerate_*.R",
       purpose = "Single-cell RNA-seq analysis: QC, normalization, UMAP, DE"),

  list(pkg = "qs",       source = "CRAN", min_ver = "0.25.0",
       scripts = "00–04, 06–08, 11, 13, regenerate_hugo_figures.R",
       purpose = "Fast serialization: qread/qsave for Seurat objects"),

  list(pkg = "AUCell",   source = "Bioconductor", min_ver = "1.20.0",
       scripts = "07, 08",
       purpose = "Rank-based AUC per cell for pathway activity (Tier 5)"),

  list(pkg = "fgsea",    source = "Bioconductor", min_ver = "1.24.0",
       scripts = "05, 09, regenerate_*.R",
       purpose = "Fast pre-ranked gene set enrichment analysis (GSEA)"),

  list(pkg = "msigdbr",  source = "CRAN", min_ver = "7.5.0",
       scripts = "05, 07, regenerate_*.R",
       purpose = "MSigDB gene sets (Hallmark, C2, etc.) in R"),

  # ── Overlap / Intersection analysis ────────────────────────────────────────
  list(pkg = "UpSetR",   source = "CRAN", min_ver = "1.4.0",
       scripts = "09",
       purpose = "UpSet plots for multi-set intersection visualization (A4)"),

  list(pkg = "ComplexUpset", source = "CRAN", min_ver = "1.3.3",
       scripts = "(optional alternative to UpSetR)",
       purpose = "ggplot2-compatible UpSet plots with richer aesthetics"),

  # ── Statistical tests ──────────────────────────────────────────────────────
  list(pkg = "diptest",  source = "CRAN", min_ver = "0.76.0",
       scripts = "13",
       purpose = "Hartigan's dip test for bimodality / unimodality (A5)")
)

# ── 2. Audit installed packages ────────────────────────────────────────────────
cat("=============================================================\n")
cat(" PACKAGE AUDIT — PAIR-Perturb-Seq Pipeline\n")
cat("=============================================================\n\n")

installed_pkgs <- rownames(installed.packages(lib.loc = LIB))

status_tbl <- lapply(manifest, function(m) {
  is_inst <- m$pkg %in% installed_pkgs
  cur_ver  <- tryCatch(as.character(packageVersion(m$pkg)), error = function(e) NA_character_)
  meets_min <- if (!is.na(m$min_ver) && !is.na(cur_ver)) {
    tryCatch(utils::compareVersion(cur_ver, m$min_ver) >= 0, error = function(e) NA)
  } else TRUE

  list(
    Package   = m$pkg,
    Source    = m$source,
    Min_ver   = if (is.na(m$min_ver)) "any" else m$min_ver,
    Installed = is_inst,
    Version   = if (is.na(cur_ver)) "—" else cur_ver,
    OK        = is_inst && isTRUE(meets_min),
    Scripts   = m$scripts,
    Purpose   = m$purpose
  )
})

# Pretty-print
cat(sprintf("%-18s %-14s %-8s %-12s %-8s\n",
            "Package", "Source", "Min ver", "Installed", "Status"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (s in status_tbl) {
  status_icon <- if (isTRUE(s$OK)) "[OK]" else if (s$Installed) "[OLD]" else "[MISSING]"
  cat(sprintf("%-18s %-14s %-8s %-12s %-8s\n",
              s$Package, s$Source, s$Min_ver,
              if (s$Installed) s$Version else "not found",
              status_icon))
}
cat("\n")

# ── 3. Identify what needs installation ────────────────────────────────────────
missing_cran  <- Filter(function(s) !s$Installed && s$Source == "CRAN",  status_tbl)
missing_bioc  <- Filter(function(s) !s$Installed && s$Source == "Bioconductor", status_tbl)
outdated      <- Filter(function(s) s$Installed && !isTRUE(s$OK), status_tbl)
missing_all   <- c(missing_cran, missing_bioc)

if (length(missing_all) == 0 && length(outdated) == 0) {
  cat("All packages are installed and meet minimum version requirements.\n")
  cat("No installation needed.\n\n")
} else {
  if (length(outdated) > 0) {
    cat("Packages below minimum version:\n")
    for (s in outdated) cat(" -", s$Package, "(installed:", s$Version, "| required:", s$Min_ver, ")\n")
    cat("\n")
  }
  if (length(missing_all) > 0) {
    cat("Missing packages:\n")
    for (s in missing_all) cat(" -", s$Package, "(", s$Source, ")\n")
    cat("\n")
  }
}

# ── 4. Install missing / outdated packages ─────────────────────────────────────
# Set INSTALL_MISSING <- TRUE to actually install; FALSE for dry-run audit only
INSTALL_MISSING <- TRUE

if (INSTALL_MISSING && (length(missing_all) > 0 || length(outdated) > 0)) {
  cat("=============================================================\n")
  cat(" INSTALLING MISSING / OUTDATED PACKAGES\n")
  cat("=============================================================\n\n")

  # ── 4a. Ensure BiocManager is available ───────────────────────────────────
  if (!"BiocManager" %in% rownames(installed.packages())) {
    cat("Installing BiocManager (required for Bioconductor packages)...\n")
    install.packages("BiocManager", lib = LIB, repos = "https://cran.r-project.org")
  }
  library(BiocManager)

  # ── 4b. CRAN packages ─────────────────────────────────────────────────────
  cran_outdated   <- Filter(function(s) s$Source == "CRAN", outdated)
  to_install_cran <- unique(c(
    vapply(missing_cran,  function(s) s$pkg, character(1)),
    vapply(cran_outdated, function(s) s$pkg, character(1))
  ))
  if (length(to_install_cran) > 0) {
    cat("Installing from CRAN:", paste(to_install_cran, collapse = ", "), "\n")
    install.packages(to_install_cran,
                     lib   = LIB,
                     repos = "https://cran.r-project.org",
                     dependencies = TRUE)
  } else {
    cat("No CRAN packages to install.\n")
  }

  # ── 4c. Bioconductor packages ──────────────────────────────────────────────
  bioc_outdated   <- Filter(function(s) s$Source == "Bioconductor", outdated)
  to_install_bioc <- unique(c(
    vapply(missing_bioc,  function(s) s$pkg, character(1)),
    vapply(bioc_outdated, function(s) s$pkg, character(1))
  ))
  if (length(to_install_bioc) > 0) {
    cat("Installing from Bioconductor:", paste(to_install_bioc, collapse = ", "), "\n")
    BiocManager::install(to_install_bioc, lib = LIB, update = FALSE, ask = FALSE)
  } else {
    cat("No Bioconductor packages to install.\n")
  }

  cat("\nInstallation complete. Re-running audit...\n\n")

  # ── 4d. Post-install verification ──────────────────────────────────────────
  installed_pkgs2 <- rownames(installed.packages(lib.loc = LIB))
  still_missing <- setdiff(
    sapply(missing_all, `[[`, "pkg"),
    installed_pkgs2
  )
  if (length(still_missing) > 0) {
    cat("WARNING — Still not installed after attempt:\n")
    for (p in still_missing) cat(" -", p, "\n")
    cat("  Please install manually or check your internet connection.\n\n")
  } else {
    cat("All previously missing packages are now installed.\n\n")
  }
} else if (!INSTALL_MISSING && length(missing_all) > 0) {
  cat("Dry-run mode (INSTALL_MISSING = FALSE). Set to TRUE to install.\n\n")
}

# ── 5. Full package × script mapping table ────────────────────────────────────
cat("=============================================================\n")
cat(" PACKAGE × SCRIPT USAGE MAP\n")
cat("=============================================================\n\n")

cat(sprintf("%-18s %-14s %-10s  %s\n", "Package", "Source", "Version", "Scripts"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (s in status_tbl) {
  cat(sprintf("%-18s %-14s %-10s  %s\n",
              s$Package, s$Source, s$Version, s$Scripts))
}
cat("\n")

# ── 6. Purpose reference table ─────────────────────────────────────────────────
cat("=============================================================\n")
cat(" PACKAGE PURPOSE REFERENCE\n")
cat("=============================================================\n\n")

cat(sprintf("%-18s  %s\n", "Package", "Purpose"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (s in status_tbl) {
  cat(sprintf("%-18s  %s\n", s$Package, s$Purpose))
}
cat("\n")

# ── 7. Session info snapshot ───────────────────────────────────────────────────
cat("=============================================================\n")
cat(" SESSION INFO SNAPSHOT\n")
cat("=============================================================\n")
cat("R version  :", R.version$version.string, "\n")
cat("Platform   :", R.version$platform, "\n")
cat("Library    :", LIB, "\n")
cat("Date run   :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Quietly load all installed packages to confirm no load errors
pkg_load_status <- sapply(sapply(manifest, `[[`, "pkg"), function(p) {
  if (p == "grid") return("base — always available")
  if (!p %in% rownames(installed.packages(lib.loc = LIB))) return("NOT INSTALLED")
  tryCatch({
    suppressPackageStartupMessages(library(p, lib.loc = LIB,
                                           character.only = TRUE,
                                           quietly = TRUE))
    paste0("OK (", packageVersion(p), ")")
  }, error = function(e) paste0("LOAD ERROR: ", conditionMessage(e)))
})

cat("Package load test:\n")
for (nm in names(pkg_load_status)) {
  cat(sprintf("  %-18s %s\n", nm, pkg_load_status[nm]))
}
cat("\nDone.\n")

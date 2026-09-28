# =============================================================================
# 00_setup.R — shared paths, seeds and helpers
# =============================================================================
# Sourced at the top of every analysis script:
#
#   source("R/00_setup.R")
#
# It defines repo_root / data_dir / results_dir / tables_dir / figures_dir and
# a small set of helpers. No analysis happens here.
#
# Environment this work was run in: R 4.6.0 (2026-04-24 ucrt).
# See sessionInfo.txt for exact package versions.
# =============================================================================

# --- Repository root ---------------------------------------------------------
# Resolution order:
#   1. the REPO_ROOT environment variable, if set
#   2. if sourced from an interactive RStudio session, the project directory
#   3. otherwise the working directory
# Every script is expected to be run from the repository root, e.g.
#   Rscript R/03_RRA_discovery.R
if (!exists("repo_root")) {
  repo_root <- Sys.getenv("REPO_ROOT", unset = "")
  if (!nzchar(repo_root)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      repo_root <- rstudioapi::getActiveProject()
    } else {
      repo_root <- getwd()
    }
  }
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)
}

data_dir    <- file.path(repo_root, "data")
results_dir <- file.path(repo_root, "results")
tables_dir  <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")

invisible(lapply(c(tables_dir, figures_dir),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# --- Sources (accessions only; no data are stored in this repository) --------
# Liver discovery cohorts
liver_cohorts  <- c("GSE48452", "GSE89632", "GSE83452",
                    "GSE126848", "GSE130970", "GSE162694")
# Gastric discovery cohorts
gastric_cohorts <- c("GSE116312", "GSE27411")
# Held-out external validation cohorts (never used for candidate selection)
validation_cohorts <- c("GSE135251", "GSE153224")
# Single-cell / single-nucleus datasets
singlecell_datasets <- c("GSE202379", "GSE115469", "GSE134520")

discovery_cohorts <- c(liver_cohorts, gastric_cohorts)
all_cohorts       <- c(discovery_cohorts, validation_cohorts)

# --- Fixed analysis parameters ----------------------------------------------
# These are the values stated in the manuscript. They are repeated here so that
# no number in the analysis is silently retyped inside a script body.

# Robust rank aggregation
# NOTE: two background sizes are used in this study and they are NOT the same number.
#   - the discovery screen built its ranked lists and aggregated them with the runtime
#     size of the union of symbols across the six liver cohorts:
#         N <- length(unique(unlist(lapply(six, function(d) unique(d$symbol)))))
#     (see R/02_RRA_discovery.R)
#   - the robustness audit recomputed the aggregation with a fixed background of
#     20,000, which was selected after comparing three definitions: the six-cohort
#     common-gene count (11,982), per-cohort sizes, and 20,000.
# The manuscript quotes 20,000; that value belongs to the audit. Both are recorded here
# so the distinction survives into the repository.
RRA_BACKGROUND_N_RUNTIME <- NA_integer_  # resolved at run time from the cohort union
RRA_BACKGROUND_N_AUDIT   <- 20000L       # fixed background used for the LODO audit
RRA_SCORE_THRESHOLD      <- 0.01         # RRA score cutoff for the discovery screen
RRA_MIN_AGREEING_COHORTS <- 3L           # >= 3 of 6 liver cohorts, same direction
RRA_TOP_N_PER_DIRECTION  <- 200L         # up and down lists truncated to this size
CANDIDATE_MIN_AGREEING_COHORTS <- 4L     # >= 4 of 6, for the final 14 candidates

# Cross-organ selection
GASTRIC_DEG_P    <- 0.05
GASTRIC_DEG_LOGF <- 0.5                  # |log2FC| cutoff in the gastric cohorts
GASTRIC_SUPPORT_LOGF <- log2(1.5)        # audit-level support threshold
LODO_RANK_CUTOFF <- 200L                 # 'in200' retention criterion
LODO_P_CUTOFF    <- 0.01                 # 'lt001' retention criterion

# Single-cell pseudobulk
MIN_CELLS_PER_PATIENT_CELLTYPE <- 20L  # minimum cells for a patient-cell-type unit
PSEUDOBULK_N_GENES     <- 13L   # candidate genes tested
PSEUDOBULK_N_CELLTYPES <- 8L    # cell types tested
PSEUDOBULK_N_TESTS     <- PSEUDOBULK_N_GENES * PSEUDOBULK_N_CELLTYPES  # = 104
BONFERRONI_ALPHA       <- 0.05 / PSEUDOBULK_N_TESTS  # = 4.8077e-04

# Adjusted P-value threshold quoted in the manuscript
MANUSCRIPT_BONFERRONI_P <- 4.8e-4

# Multiple-testing families
N_EFFECT_HEATMAP_TESTS <- 140L  # 14 candidate genes x 10 cohorts (Fig 1)
N_VALIDATION_CANDIDATES <- 5L   # pre-specified candidates per validation cohort
N_MR_ESTIMATES         <- 69L

# --- Candidate genes ---------------------------------------------------------
candidates_14 <- c("IL32", "CDHR2", "ANXA4", "CADM2", "LGALS3",
                   "RPS6KA1", "GOLM1", "ANO10", "SLC6A16", "KIAA1958",
                   "FGA", "FGB", "LEPR", "TSPAN3")
candidates_main <- c("IL32", "CDHR2", "ANXA4", "CADM2", "LGALS3")
candidate_observation <- "RPS6KA1"
candidates_downgraded <- c("GOLM1", "ANO10", "SLC6A16", "KIAA1958",
                           "FGA", "FGB", "LEPR", "TSPAN3")

stopifnot(length(candidates_14) == 14L,
          length(candidates_main) == 5L,
          length(candidates_downgraded) == 8L,
          all(c(candidates_main, candidate_observation,
                candidates_downgraded) %in% candidates_14))

# --- Seeded randomness -------------------------------------------------------
# Every stochastic step in this repository must go through seed_with() so that
# the seed actually used is recorded in the log rather than living in a comment.
PERMUTATION_SEED  <- 20260921L   # <-- [CONFIRM] against the original script
PERMUTATION_N     <- 10000L      # <-- [CONFIRM] against the original script

seed_with <- function(seed, label = NULL) {
  if (!is.null(label)) message(sprintf("[seed] %s: set.seed(%s)", label, seed))
  set.seed(seed)
  invisible(seed)
}

# --- Helpers -----------------------------------------------------------------

# Assert that a dataset folder exists before a script tries to read from it.
require_dataset <- function(accession, hint = NULL) {
  path <- file.path(data_dir, accession)
  if (!dir.exists(path)) {
    stop(sprintf(paste0("Dataset folder not found: %s\n",
                        "  Download %s into data/ first - see data/datasets.md.%s"),
                 path, accession,
                 if (!is.null(hint)) paste0("\n  ", hint) else ""),
         call. = FALSE)
  }
  path
}

# Write a table into results/tables/ with a consistent name.
write_table <- function(x, filename, row.names = FALSE, ...) {
  stopifnot(grepl("\\.(csv|tsv|txt)$", filename))
  out <- file.path(tables_dir, filename)
  utils::write.table(x, out, sep = if (grepl("\\.tsv$", filename)) "\t" else ",",
                     row.names = row.names, quote = FALSE, ...)
  message("[out] ", out)
  invisible(out)
}

# Save a figure into results/figures/.
save_figure <- function(plot, filename, width = 180, height = 120, dpi = 300) {
  out <- file.path(figures_dir, filename)
  ggplot2::ggsave(out, plot = plot, width = width, height = height,
                  units = "mm", dpi = dpi)
  message("[out] ", out)
  invisible(out)
}

# Confirm the seed actually used in a stochastic analysis, and log it.
report_seed <- function(seed, n, what) {
  message(sprintf("[audit] %s: seed = %s, n = %s", what, seed, n))
  invisible(list(seed = seed, n = n, what = what))
}

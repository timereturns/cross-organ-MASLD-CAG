# =============================================================================
# 00_setup.R — shared paths, parameters, seeds and helpers
# =============================================================================
# Sourced at the top of every analysis script:
#
#   source("R/00_setup.R")
#
# It defines repo_root / data_dir / results_dir / tables_dir / figures_dir, the
# fixed analysis parameters, and the helper functions the pipeline shares. No
# analysis happens here.
#
# Environment this work was run in: R 4.6.0 (2026-04-24 ucrt).
# See sessionInfo.txt for exact package versions.
#
# Restored from the original working files; see docs/build_notes.md for what was
# changed and what was left alone.
# =============================================================================

# --- Repository root ---------------------------------------------------------
# Resolution order:
#   1. the REPO_ROOT environment variable, if set
#   2. if sourced from an interactive RStudio session, the project directory
#   3. otherwise the working directory
# Every script is expected to be run from the repository root, e.g.
#   Rscript R/02_RRA_discovery.R
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
cache_dir   <- file.path(repo_root, "cache")   # large intermediate objects

invisible(lapply(c(tables_dir, figures_dir, cache_dir),
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# --- Sources (accessions only; no data are stored in this repository) --------
# Liver discovery cohorts
liver_cohorts  <- c("Ahrens", "Arendt", "Lefebvre",
                    "Suppli", "Hoang", "Pantano")
# Accession <-> internal label
cohort_accession <- c(Ahrens = "GSE48452", Arendt = "GSE89632",
                      Lefebvre = "GSE83452", Suppli = "GSE126848",
                      Hoang = "GSE130970", Pantano = "GSE162694",
                      GSE116312 = "GSE116312", GSE27411 = "GSE27411")
# Gastric discovery cohorts (labelled by accession throughout)
gastric_cohorts <- c("GSE116312", "GSE27411")
# Held-out external validation cohorts (never used for candidate selection)
validation_cohorts <- c("GSE135251", "GSE153224")
# Single-cell / single-nucleus datasets
singlecell_datasets <- c("GSE202379", "GSE115469", "GSE134520")
# Fourth bulk cohort, deconvolution sensitivity analysis only
deconvolution_cohorts <- c("GSE126848", "GSE130970", "GSE162694", "GSE174478")

discovery_cohorts <- c(liver_cohorts, gastric_cohorts)

# --- Cohort sizes, asserted against the actual data -------------------------
# These are the counts the original audit hard-checked. Keep them as assertions,
# not as documentation: a grouping bug should fail loudly.
expected_case <- c(Ahrens = 26L, Arendt = 39L, Lefebvre = 104L,
                   Suppli = 31L, Hoang = 53L, Pantano = 112L)
expected_ctrl <- c(Ahrens = 28L, Arendt = 24L, Lefebvre = 44L,
                   Suppli = 26L, Hoang = 25L, Pantano = 31L)
expected_samples <- c(Ahrens = 54L, Arendt = 63L, Lefebvre = 148L,
                      Suppli = 57L, Hoang = 78L, Pantano = 143L,
                      GSE135251 = 216L, GSE174478 = 94L)

# --- Robust rank aggregation -------------------------------------------------
# NOTE: two background sizes appear in this study and they are NOT the same number.
#
#   - The discovery screen aggregated with N <- length(all_sym), the runtime size
#     of the union of symbols across the six liver cohorts, computed at run time.
#   - The robustness audit recomputed the aggregation with a fixed 20,000, after
#     comparing three definitions: the six-cohort common-gene count (11,982),
#     per-cohort gene counts, and 20,000. 20,000 was then chosen by hand.
#
# The manuscript quotes 20,000, which belongs to the audit. Both are recorded so
# the distinction survives into the repository. See manuscript_discrepancies.md.
RRA_BACKGROUND_N_RUNTIME <- NA_integer_  # resolved at run time from the cohort union
RRA_BACKGROUND_N_AUDIT   <- 20000L       # fixed background, LODO audit
RRA_SCORE_THRESHOLD      <- 0.01         # final score cutoff (first pass used 0.05)
RRA_SCORE_THRESHOLD_FIRST_PASS <- 0.05   # recorded for transparency
RRA_MIN_AGREEING_COHORTS <- 3L           # >= 3 of 6 liver cohorts, same direction
RRA_TOP_N_PER_DIRECTION  <- 200L

# --- Cross-organ candidate finalisation --------------------------------------
GASTRIC_DEG_P        <- 0.05
GASTRIC_DEG_LOGF     <- 0.5              # |log2FC| in the gastric cohorts
GASTRIC_SUPPORT_LOGF <- log2(1.5)        # audit-level support threshold
CANDIDATE_MIN_AGREEING_COHORTS <- 4L     # >= 4 of 6 liver cohorts, same direction

# --- Robustness audit --------------------------------------------------------
LODO_RANK_CUTOFF <- 200L                 # 'in200' retention criterion
LODO_P_CUTOFF    <- 0.01                 # 'lt001' retention criterion
META_METHOD      <- "REML"               # metafor::rma method

# --- Single-cell pseudobulk --------------------------------------------------
MIN_CELLS_PER_PATIENT_CELLTYPE <- 20L    # minimum cells for a patient-cell-type unit
FIBROSIS_STAGE_RULE <- "mode"            # mode of per-cell score; rounded mean was
                                         # computed as a comparator and rejected

# The Bonferroni family quoted in the manuscript: 13 genes x 8 cell types = 104.
# NOTE: the composition-adjustment marker panel in R/06 has 9 cell types and the
# candidate vectors elsewhere are 6, 13 and 14 long depending on the block. The
# 13 x 8 pair is the family behind the manuscript's threshold; see
# manuscript_discrepancies.md item 3.
PSEUDOBULK_N_GENES     <- 13L
PSEUDOBULK_N_CELLTYPES <- 8L
PSEUDOBULK_N_TESTS     <- PSEUDOBULK_N_GENES * PSEUDOBULK_N_CELLTYPES  # 104
BONFERRONI_ALPHA       <- 0.05 / PSEUDOBULK_N_TESTS                    # 4.8077e-04
MANUSCRIPT_BONFERRONI_P <- 4.8e-4

# --- Multiple-testing families ----------------------------------------------
N_EFFECT_HEATMAP_TESTS  <- 140L          # 14 candidate genes x 10 cohorts (Fig 1)
N_VALIDATION_CANDIDATES <- 5L            # pre-specified candidates per validation cohort
N_MR_ESTIMATES          <- 69L

# --- Deconvolution -----------------------------------------------------------
MUSIC_MARKERS_PER_CELLTYPE <- 50L        # top-N by specificity, per cell type
MUSIC_MARKER_PANEL_MAX     <- 400L       # 50 x 8 cell types before deduplication;
                                         # the true size is <= this and is recorded
                                         # at run time rather than asserted
MUSIC_COVARIATE_MIN_PROP   <- 0.02       # proportion above which a cell type enters
MUSIC_FORCED_COVARIATES    <- c("Hepatocytes", "Stellate")
BULK_GATE_MAX_READS        <- 50         # max value below which data are not raw counts

# --- Mendelian randomisation -------------------------------------------------
MR_F_STATISTIC_MIN <- 10                 # per-SNP instrument strength filter
MR_PRESTO_NB       <- 1000L              # run_mr_presso NbDistribution
MR_CLUMPING_R2     <- 0.001              # from the instrument directory name
                                         # eqtlclump_r2_0.001 -- window and panel are
                                         # not recorded in any surviving script

# --- Permutation: there is none ----------------------------------------------
# The manuscript Methods describes "20,000 rank-permutation replicates" and Table
# S02 is titled "Leave-one-dataset-out exact RRA (20,000 permutations)". Neither
# describes a permutation analysis. There is no Monte Carlo step anywhere in this
# pipeline:
#
#   * 20,000 is the RRA BACKGROUND GENE COUNT (RRA_BACKGROUND_N_AUDIT above), passed
#     as the denominator of the closed-form exact P-value;
#   * the leave-one-dataset-out audit recomputes that closed form once per dropped
#     cohort -- six iterations, and every output is a deterministic function of the
#     input ranks.
#
# Because nothing is resampled, there is no seed to set and no replicate count to
# record. The LODO table is fully reproducible from R/05_audit_LODO_meta.R.
#
# The only genuinely stochastic procedure in the project is the AUC bootstrap, below.
# See docs/manuscript_discrepancies.md item 1 for the three places the manuscript
# wording needs fixing.

# --- Bootstrap and jitter -----------------------------------------------------
# The AUC bootstrap for the GSE135251 validation: 2,000 resamples, percentile CI.
# This is the single seed that appears in the original working files
# (第一次预审稿文件处理.R, two blocks).
BOOTSTRAP_SEED <- 20260914L
BOOTSTRAP_N    <- 2000L

# Figure jitter positions. Added by this repository: the original used base-R
# stripchart(method = "jitter") with no seed, so point positions were not
# reproducible.
JITTER_SEED    <- 42L

# --- Candidate genes ---------------------------------------------------------
candidates_14 <- c("IL32", "CDHR2", "ANXA4", "CADM2", "LGALS3",
                   "RPS6KA1", "GOLM1", "ANO10", "SLC6A16", "KIAA1958",
                   "FGA", "FGB", "LEPR", "TSPAN3")
candidates_main        <- c("IL32", "CDHR2", "ANXA4", "CADM2", "LGALS3")
candidate_observation  <- "RPS6KA1"
candidates_downgraded  <- c("GOLM1", "ANO10", "SLC6A16", "KIAA1958",
                            "FGA", "FGB", "LEPR", "TSPAN3")

stopifnot(length(candidates_14) == 14L,
          length(candidates_main) == 5L,
          length(candidates_downgraded) == 8L,
          all(c(candidates_main, candidate_observation,
                candidates_downgraded) %in% candidates_14))

# Expected direction of change, from candidate_genes_final.csv. Hardcoded in the
# original figure script for the 5-gene panels; kept here so the two sources can
# be compared rather than silently diverging.
expected_direction <- c(IL32 = "up", CDHR2 = "up", ANXA4 = "up",
                        CADM2 = "dn", LGALS3 = "up", RPS6KA1 = "up",
                        GOLM1 = "up", ANO10 = "up", SLC6A16 = "dn",
                        KIAA1958 = "dn", FGA = "up", FGB = "up",
                        LEPR = "up", TSPAN3 = "up")

# --- Liver cell-type marker panel (composition adjustment, R/06) -------------
# 9 liver lineages. Used for the marker-score composition adjustment, which is a
# different analysis from the MuSiC reference (8 cell types) -- see
# manuscript_discrepancies.md item 4.
mk_liver <- list(
  Hepatocyte    = c("ALB", "CYP3A4", "CYP2E1", "HNF4A", "ASGR1", "TTR", "FABP1"),
  Cholangiocyte = c("KRT19", "KRT7", "SOX9", "SCTR"),
  Kupffer_macro = c("CD163", "CD68", "MARCO", "C1QA", "C1QB", "MS4A7"),
  T_cell        = c("CD3D", "CD3E", "CD2", "CD8A", "IL7R", "TRAC", "GZMA"),
  B_plasma      = c("MS4A1", "CD79A", "CD79B", "IGKC", "IGHG1", "MZB1"),
  NK            = c("KLRD1", "NKG7", "KLRF1"),
  Neutrophil    = c("CSF3R", "S100A8", "S100A9", "FCGR3B"),
  HSC_fibro     = c("ACTA2", "COL1A1", "COL3A1", "PDGFRB", "LUM", "DCN"),
  Endothelial   = c("PECAM1", "VWF", "CDH5", "ENG", "KDR")
)

# --- Mendelian randomisation gene list ---------------------------------------
# 10 of the 14 candidates have cis-eQTL instruments in the blood eQTL set.
mr_genes <- c("IL32", "GOLM1", "TSPAN3", "ANXA4", "RPS6KA1", "ANO10",
              "SLC6A16", "KIAA1958", "LGALS3", "CADM2")
mr_direction <- c(IL32 = "up", GOLM1 = "up", TSPAN3 = "up", ANXA4 = "up",
                  RPS6KA1 = "up", ANO10 = "up", SLC6A16 = "dn",
                  KIAA1958 = "dn", LGALS3 = "up", CADM2 = "dn")

# --- Seeded randomness -------------------------------------------------------
# Every stochastic step in this repository goes through seed_with() so that the
# seed actually used is recorded rather than living in a comment.
seed_with <- function(seed, label = NULL) {
  if (!is.null(label)) message(sprintf("[seed] %s: set.seed(%s)", label, seed))
  set.seed(seed)
  invisible(seed)
}

# --- Helpers -----------------------------------------------------------------

# Assert that a dataset directory exists before a script reads from it.
require_dataset <- function(accession, hint = NULL) {
  path <- file.path(data_dir, accession)
  if (!dir.exists(path)) {
    stop(sprintf(paste0("Dataset dir not found: %s\n",
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

# --- GEO utilities -----------------------------------------------------------
# These three readers appear, in near-identical form, in several of the original
# working files. They are defined once here.

# Read a series matrix expression table. Handles gzipped or plain files.
read_series_matrix_expr <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  L <- readLines(con, warn = FALSE)
  b <- grep("^!series_matrix_table_begin", L)
  e <- grep("^!series_matrix_table_end", L)
  if (length(b) != 1L || length(e) != 1L || e <= b) {
    stop("series matrix table markers not found in ", path, call. = FALSE)
  }
  dat <- read.table(text = L[(b + 1):(e - 1)], header = TRUE, sep = "\t",
                    check.names = FALSE, comment.char = "",
                    stringsAsFactors = FALSE)
  mat <- as.matrix(dat[, -1])
  rownames(mat) <- dat[, 1]
  storage.mode(mat) <- "double"
  mat
}

# Pull one !Sample_characteristics_ch1 field, named by GSM.
get_char <- function(path, tag) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  L <- readLines(con, warn = FALSE)
  ga  <- grep("^!Sample_geo_accession", L, value = TRUE)[1]
  gsm <- gsub('^"|"$', "", strsplit(ga, "\t")[[1]][-1])
  lines <- L[startsWith(L, paste0('!Sample_characteristics_ch1\t"', tag))]
  if (length(lines) == 0L) stop("characteristic not found: ", tag, call. = FALSE)
  x <- gsub('^"|"$', "", strsplit(lines[1], "\t")[[1]][-1])
  setNames(sub(paste0(tag, ": "), "", x, fixed = TRUE), gsm)
}

# GSM / title pairs.
read_titles <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  L <- readLines(con, warn = FALSE)
  ti <- gsub('^"|"$', "", strsplit(grep("^!Sample_title", L, value = TRUE)[1],
                                   "\t")[[1]][-1])
  ga <- gsub('^"|"$', "", strsplit(grep("^!Sample_geo_accession", L,
                                        value = TRUE)[1], "\t")[[1]][-1])
  data.frame(gsm = ga, title = ti, stringsAsFactors = FALSE)
}

# Full sample metadata table from a series matrix, one column per characteristics
# field, plus gsm.
meta_from_series <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  L <- readLines(con, warn = FALSE)
  gsm <- gsub('^"|"$', "", strsplit(grep("^!Sample_geo_accession", L,
                                         value = TRUE)[1], "\t")[[1]][-1])
  df <- data.frame(gsm = gsm, stringsAsFactors = FALSE)
  for (line in L[startsWith(L, "!Sample_characteristics_ch1")]) {
    x <- gsub('^"|"$', "", strsplit(line, "\t")[[1]][-1])
    df[[sub(":.*$", "", x[1])]] <- sub("^.*?: ", "", x)
  }
  df
}

# Print the value distribution of every characteristics field - used to work out
# each cohort's grouping rule.
summarize_chars <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  L <- readLines(con, warn = FALSE)
  for (line in grep("^!Sample_characteristics_ch1", L, value = TRUE)) {
    x <- gsub('^"|"$', "", strsplit(line, "\t")[[1]][-1])
    cat("\n==== field:", sub(":.*$", "", x[1]), "====\n")
    print(table(sub("^.*?: ", "", x)))
  }
  invisible(NULL)
}

# --- limma helpers -----------------------------------------------------------

# Array pipeline: design ~ 0 + grp, contrast case - control, probe -> symbol,
# keeping the most significant row per symbol. This is the original
# run_limma_fold(), unchanged apart from the annotation lookup being passed in.
run_limma_fold <- function(mat, grp, anno_db) {
  grp <- grp[colnames(mat)]
  if (anyNA(grp)) stop("group vector has NA after matching columns", call. = FALSE)
  grp <- factor(grp, levels = c("control", "case"))
  sym <- AnnotationDbi::mapIds(anno_db, keys = rownames(mat), column = "SYMBOL",
                               keytype = "PROBEID", multiVals = "first")
  keep <- !is.na(sym) & sym != "---"
  mat2 <- mat[keep, , drop = FALSE]
  sym2 <- sym[keep]
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  fit <- limma::eBayes(limma::contrasts.fit(
    limma::lmFit(mat2, design),
    limma::makeContrasts(case - control, levels = design)))
  tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
  tt$symbol <- sym2[rownames(tt)]
  tt <- tt[order(-abs(tt$t)), , drop = FALSE]
  tt[!duplicated(tt$symbol), , drop = FALSE]
}

# Same, for a platform whose annotation is supplied as a probe -> symbol vector
# (used for GPL6255, parsed from the platform annotation file).
run_limma_symbol <- function(mat, grp, sym) {
  grp <- grp[colnames(mat)]
  if (anyNA(grp)) stop("group vector has NA after matching columns", call. = FALSE)
  grp <- factor(grp, levels = c("control", "case"))
  keep <- !is.na(sym) & sym != "" & sym != "---"
  mat2 <- mat[keep, , drop = FALSE]
  sym2 <- sym[keep]
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  fit <- limma::eBayes(limma::contrasts.fit(
    limma::lmFit(mat2, design),
    limma::makeContrasts(case - control, levels = design)))
  tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
  tt$symbol <- sym2[rownames(tt)]
  tt <- tt[order(-abs(tt$t)), , drop = FALSE]
  tt[!duplicated(tt$symbol), , drop = FALSE]
}

# RNA-seq pipeline: voom on raw counts, then limma, then ID -> symbol.
#   id_keytype is "ENSEMBL" or "ENTREZID".
voom_limma <- function(cnt, grp, id_keytype = "ENSEMBL") {
  grp <- grp[colnames(cnt)]
  if (anyNA(grp)) stop("group vector has NA after matching columns", call. = FALSE)
  grp <- factor(grp, levels = c("control", "case"))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  v <- limma::voom(cnt, design)
  fit <- limma::eBayes(limma::contrasts.fit(
    limma::lmFit(v, design),
    limma::makeContrasts(case - control, levels = design)))
  tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
  ids <- sub("\\..*$", "", rownames(tt))
  sym <- AnnotationDbi::mapIds(hs_eg_db(), keys = ids,
                               keytype = id_keytype, column = "SYMBOL",
                               multiVals = "first")
  tt$symbol <- sym
  tt <- tt[!is.na(tt$symbol), , drop = FALSE]
  tt <- tt[order(-abs(tt$t)), , drop = FALSE]
  tt[!duplicated(tt$symbol), , drop = FALSE]
}

# Access the org.Hs.eg.db annotation object without attaching the package.
hs_eg_db <- function() {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db is required; see README section 3", call. = FALSE)
  }
  getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
}

# --- Closed-form robust rank aggregation -------------------------------------
# The exact implementation used for the audit and the LODO recomputation:
#
#     rho = min_k Beta(x_(k); k, n - k + 1)
#
# where x is the sorted vector of rank/N. This is the same quantity the
# RobustRankAggreg package computes; the audit used this function so that the
# per-gene ranking could be recomputed with an arbitrary denominator.
rra_rho <- function(ranks, denom) {
  r <- as.numeric(ranks)
  r <- r[!is.na(r)]
  if (length(r) < 2L) return(NA_real_)
  if (length(denom) == 1L) denom <- rep(denom, length(r))
  x <- sort(r / denom)
  n <- length(x)
  min(vapply(seq_len(n), function(k) pbeta(x[k], k, n - k + 1), numeric(1)))
}

# Build a gene x cohort matrix of integer ranks from a list of limma tables.
# direction = "up" ranks by descending t, "dn" by ascending t.
# Restricts to genes present in every cohort, matching the discovery screen.
build_rankmat <- function(dat_list, direction) {
  out <- lapply(dat_list, function(d) {
    r <- rank(if (direction == "up") -d$t else d$t, ties.method = "min")
    structure(as.integer(r), names = d$symbol)
  })
  cmn <- Reduce(intersect, lapply(out, names))
  do.call(cbind, lapply(out, function(x) x[cmn]))
}

# --- misc --------------------------------------------------------------------

# Confirm a stochastic analysis's settings and log them.
report_seed <- function(seed, n, what) {
  message(sprintf("[audit] %s: seed = %s, n = %s", what, seed, n))
  invisible(list(seed = seed, n = n, what = what))
}

message("[setup] repo_root = ", repo_root)
message("[setup] R ", getRversion(), " | ",
        length(rownames(installed.packages())), " packages installed")

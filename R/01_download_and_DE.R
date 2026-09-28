# =============================================================================
# 01_download_and_DE.R — cohort acquisition and per-cohort differential expression
# =============================================================================
# STATUS: RE-IMPLEMENTATION.
#
# The original working file for the liver side carried a third-party vendor notice
# asserting that use of the vendor's organised data and code to publish a paper
# without permission constitutes unauthorised use. It is therefore not
# redistributed with this repository, and this script is a fresh implementation of
# the same pipeline from the public GEO source files, using the parameters that are
# documented in docs/analysis_parameters.md and recorded in R/00_setup.R.
#
# Because this is a re-implementation rather than a restoration, it must be
# validated against the reported numbers before being relied on: the per-cohort
# limma outputs it produces should reproduce the tables in results/tables/, and
# the assertion at the end of the script checks the cohort sizes and group balances
# that the original audit hard-coded.
#
# What this script does
#   - downloads the series matrices and supplementary count files from GEO
#   - derives each cohort's case/control grouping from its own metadata field
#   - runs limma per cohort (arrays directly; RNA-seq via voom)
#   - collapses probes/IDs to gene symbols
#   - writes one limma table per cohort to cache/, plus a sample-level group table
#
# The gastric cohorts are included here as well, so that all eight discovery
# cohorts are produced in one place. Candidate selection happens in 02 and 03.
#
# Run:  Rscript R/01_download_and_DE.R
# =============================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi)
})

# -----------------------------------------------------------------------------
# 0. Download
# -----------------------------------------------------------------------------
# Files are only fetched when absent, so a re-run is cheap. Everything lands in
# data/<accession>/. See data/datasets.md for the accessions and their roles.

geo_base <- "https://ftp.ncbi.nlm.nih.gov/geo/series"

geo_files <- list(
  # --- arrays: liver ---
  Ahrens   = c("GSE48nnn/GSE48452/matrix/GSE48452_series_matrix.txt.gz"),
  Arendt   = c("GSE89nnn/GSE89632/matrix/GSE89632_series_matrix.txt.gz"),
  Lefebvre = c("GSE83nnn/GSE83452/matrix/GSE83452_series_matrix.txt.gz"),
  # --- RNA-seq: liver ---
  Suppli   = c("GSE126nnn/GSE126848/matrix/GSE126848_series_matrix.txt.gz",
               "GSE126nnn/GSE126848/suppl/GSE126848_Gene_counts_raw.txt.gz"),
  Hoang    = c("GSE130nnn/GSE130970/matrix/GSE130970_series_matrix.txt.gz",
               "GSE130nnn/GSE130970/suppl/GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"),
  Pantano  = c("GSE162nnn/GSE162694/matrix/GSE162694_series_matrix.txt.gz",
               "GSE162nnn/GSE162694/suppl/GSE162694_raw_counts.csv.gz"),
  # --- gastric ---
  GSE116312 = c("GSE116nnn/GSE116312/matrix/GSE116312_series_matrix.txt.gz"),
  GSE27411  = c("GSE27nnn/GSE27411/matrix/GSE27411_series_matrix.txt.gz",
                "GSE27nnn/GSE27411/suppl/GSE27411_non-normalized.txt.gz",
                "GSE27nnn/GSE27411/soft/GPL6255.annot.gz")
)

download_cohorts <- function(skip = FALSE) {
  if (skip) {
    message("[01] download skipped (skip_download = TRUE)")
    return(invisible(NULL))
  }
  for (cohort in names(geo_files)) {
    dir <- file.path(data_dir, cohort)
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    for (rel in geo_files[[cohort]]) {
      dest <- file.path(dir, basename(rel))
      if (file.exists(dest)) {
        message("[01] present: ", basename(dest))
        next
      }
      url <- file.path(geo_base, rel)
      message("[01] downloading ", basename(rel))
      utils::download.file(url, dest, mode = "wb", quiet = TRUE)
    }
  }
}

# -----------------------------------------------------------------------------
# 1. Grouping rules, one per cohort
# -----------------------------------------------------------------------------
# Each rule reads the cohort's own metadata field and returns a named vector of
# "case"/"control" over GSM ids. Rules that exclude samples are marked: the
# exclusions are part of the design, not incidental filtering.

group_ahrens <- function(dir) {
  f <- file.path(dir, "GSE48452_series_matrix.txt.gz")
  g    <- get_char(f, "group")
  surg <- get_char(f, "bariatric surgery")
  # exclude post-surgery samples: the comparison is baseline disease vs control
  keep <- is.na(surg) | surg != "after surgery"
  setNames(ifelse(g[keep] %in% c("Nash", "Steatosis"), "case", "control"),
           names(g)[keep])
}

group_arendt <- function(dir) {
  f <- file.path(dir, "GSE89632_series_matrix.txt.gz")
  dx <- get_char(f, "diagnosis")
  setNames(ifelse(dx %in% c("NASH", "SS"), "case", "control"), names(dx))
}

group_lefebvre <- function(dir) {
  f <- file.path(dir, "GSE83452_series_matrix.txt.gz")
  st <- get_char(f, "liver status")
  tm <- get_char(f, "time")
  # baseline only; the follow-up biopsy is a different comparison
  keep <- tm == "baseline" & st %in% c("NASH", "no NASH")
  setNames(ifelse(st[keep] == "NASH", "case", "control"), names(st)[keep])
}

group_suppli <- function(dir) {
  f <- file.path(dir, "GSE126848_series_matrix.txt.gz")
  dz <- get_char(f, "disease")
  # healthy and obese are both controls; NAFLD and NASH are cases
  setNames(ifelse(dz %in% c("healthy", "obese"), "control", "case"), names(dz))
}

group_hoang <- function(dir) {
  f <- file.path(dir, "GSE130970_series_matrix.txt.gz")
  fs <- get_char(f, "fibrosis stage")
  # stage 0 is the control group
  setNames(ifelse(fs == "0", "control", "case"), names(fs))
}

group_pantano <- function(dir) {
  f <- file.path(dir, "GSE162694_series_matrix.txt")
  ti <- read_titles(f)
  # title carries a stage suffix; *_N are the controls
  setNames(ifelse(grepl("_N$", sub(".* ", "", ti$title)), "control", "case"),
           ti$gsm)
}

group_gse116312 <- function(dir) {
  f <- file.path(dir, "GSE116312_series_matrix.txt.gz")
  ti <- read_titles(f)
  type <- sub(" biopsy.*$", "", ti$title)
  # chronic atrophic gastritis vs follicular gastritis; gastric cancer excluded
  g <- ifelse(type == "Chronic atrophic gastritis", "case",
              ifelse(type == "Follicular gastritis", "control", NA))
  setNames(g, ti$gsm)
}

group_gse27411 <- function(dir) {
  f <- file.path(dir, "GSE27411_series_matrix.txt.gz")
  ti <- read_titles(f)
  type <- sub(",.*$", "", ti$title)
  # atrophy vs H. pylori non-infected; the infected group is excluded
  g <- ifelse(type == "Atrophy", "case",
              ifelse(type == "Helicobacter pylori non-infected", "control", NA))
  setNames(g, ti$gsm)
}

group_fns <- list(Ahrens = group_ahrens, Arendt = group_arendt,
                  Lefebvre = group_lefebvre, Suppli = group_suppli,
                  Hoang = group_hoang, Pantano = group_pantano,
                  GSE116312 = group_gse116312, GSE27411 = group_gse27411)

# -----------------------------------------------------------------------------
# 2. Per-cohort differential expression
# -----------------------------------------------------------------------------

# --- arrays: probe-level matrices straight from the series matrix ------------
de_ahrens <- function(dir) {
  mat <- read_series_matrix_expr(file.path(dir, "GSE48452_series_matrix.txt.gz"))
  g <- group_ahrens(dir)
  mat <- mat[, names(g), drop = FALSE]
  # GPL11532, HuGene 1.1 ST
  run_limma_fold(mat, g, get_annotation_db("hugene11sttranscriptcluster.db"))
}

de_arendt <- function(dir) {
  mat <- read_series_matrix_expr(file.path(dir, "GSE89632_series_matrix.txt.gz"))
  g <- group_arendt(dir)
  mat <- mat[, names(g), drop = FALSE]
  # GPL14951, Illumina DASL
  run_limma_fold(mat, g, get_annotation_db("illuminaHumanWGDASLv4.db"))
}

de_lefebvre <- function(dir) {
  mat <- read_series_matrix_expr(file.path(dir, "GSE83452_series_matrix.txt.gz"))
  g <- group_lefebvre(dir)
  mat <- mat[, names(g), drop = FALSE]
  # GPL16686, HuGene 2.0 ST
  run_limma_fold(mat, g, get_annotation_db("hugene20sttranscriptcluster.db"))
}

# --- RNA-seq: raw counts, voom -----------------------------------------------
de_suppli <- function(dir) {
  cnt <- read.table(gzfile(file.path(dir, "GSE126848_Gene_counts_raw.txt.gz")),
                    header = TRUE, sep = "\t", row.names = 1,
                    check.names = FALSE)
  # count-file columns are zero-padded sample numbers; series matrix carries the
  # matching description, and the GSM ids come from the accession field
  sm <- file.path(dir, "GSE126848_series_matrix.txt.gz")
  dsc <- .read_field(sm, "!Sample_description")
  gsm <- .read_field(sm, "!Sample_geo_accession")
  colnames(cnt) <- gsm[match(sub("^0+", "", colnames(cnt)),
                             sub("^0+", "", dsc))]
  g <- group_suppli(dir)[colnames(cnt)]
  stopifnot(!anyNA(g))
  voom_limma(cnt, g, "ENSEMBL")
}

de_hoang <- function(dir) {
  cnt <- read.csv(gzfile(file.path(dir,
              "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz")),
                  row.names = 1, check.names = FALSE)
  ti <- read_titles(file.path(dir, "GSE130970_series_matrix.txt.gz"))
  colnames(cnt) <- ti$gsm[match(colnames(cnt), ti$title)]
  g <- group_hoang(dir)[colnames(cnt)]
  stopifnot(!anyNA(g))
  voom_limma(cnt, g, "ENTREZID")
}

de_pantano <- function(dir) {
  cnt <- read.csv(gzfile(file.path(dir, "GSE162694_raw_counts.csv.gz")),
                  row.names = 1, check.names = FALSE)
  ti <- read_titles(file.path(dir, "GSE162694_series_matrix.txt"))
  ti$suffix <- sub(".* ", "", ti$title)
  colnames(cnt) <- ti$gsm[match(colnames(cnt), ti$suffix)]
  g <- group_pantano(dir)[colnames(cnt)]
  stopifnot(!anyNA(g))
  voom_limma(cnt, g, "ENSEMBL")
}

# --- gastric ------------------------------------------------------------------
de_gse116312 <- function(dir) {
  mat <- read_series_matrix_expr(file.path(dir, "GSE116312_series_matrix.txt.gz"))
  g <- group_gse116312(dir)
  keep <- names(g)[!is.na(g) & names(g) %in% colnames(mat)]
  mat <- mat[, keep, drop = FALSE]
  # GPL6947 / GPL10558 family: use the platform annotation shipped with the series
  run_limma_fold(mat, g[keep], get_annotation_db("illuminaHumanv3.db"))
}

de_gse27411 <- function(dir) {
  mat <- read_series_matrix_expr(file.path(dir, "GSE27411_series_matrix.txt.gz"))
  g <- group_gse27411(dir)
  keep <- names(g)[!is.na(g) & names(g) %in% colnames(mat)]
  mat <- mat[, keep, drop = FALSE]
  # Expression values are Affymetrix raw intensities on GPL570; the series matrix
  # for this dataset is already on a log scale after RMA in the original record.
  # If the range indicates raw intensities, log2-transform as the original did:
  if (max(mat, na.rm = TRUE) > 100) {
    mat <- normalizeBetweenArrays(log2(mat - min(mat, na.rm = TRUE) + 1))
  }
  sym <- read_gpl_annotation(file.path(dir, "GPL6255.annot.gz"),
                             id_col = "ID", symbol_col = "Gene symbol")
  run_limma_symbol(mat, g[keep], sym[rownames(mat)])
}

de_fns <- list(Ahrens = de_ahrens, Arendt = de_arendt, Lefebvre = de_lefebvre,
               Suppli = de_suppli, Hoang = de_hoang, Pantano = de_pantano,
               GSE116312 = de_gse116312, GSE27411 = de_gse27411)

# -----------------------------------------------------------------------------
# 3. Run
# -----------------------------------------------------------------------------

run_all <- function(skip_download = FALSE) {
  download_cohorts(skip_download)

  results <- list()
  groups  <- list()

  for (cohort in names(de_fns)) {
    dir <- file.path(data_dir, cohort)
    if (!dir.exists(dir)) {
      message("[01] SKIP ", cohort, ": directory missing")
      next
    }
    message("[01] === ", cohort, " (", cohort_accession[[cohort]], ") ===")

    g <- group_fns[[cohort]](dir)
    groups[[cohort]] <- g
    tab <- table(g, useNA = "ifany")
    message("[01]   groups: ", paste(names(tab), tab, collapse = " / "))

    res <- de_fns[[cohort]](dir)
    results[[cohort]] <- res
    message("[01]   ", nrow(res), " genes after symbol collapse")

    saveRDS(res, file.path(cache_dir, paste0("limma_", cohort, ".rds")))
  }

  # --- assertions: the counts the original audit hard-coded -------------------
  # Split by platform, matching the two objects the original pipeline kept.
  arr <- intersect(names(results), c("Ahrens", "Arendt", "Lefebvre"))
  rna <- intersect(names(results), c("Suppli", "Hoang", "Pantano"))
  gas <- intersect(names(results), c("GSE116312", "GSE27411"))

  if (length(arr) == 3L) {
    saveRDS(results[arr], file.path(cache_dir, "chips_limma_results.rds"))
    for (nm in arr) {
      obs <- sum(groups[[nm]] == "case"); ctl <- sum(groups[[nm]] == "control")
      stopifnot(obs == expected_case[[nm]], ctl == expected_ctrl[[nm]])
      message("[01]   ", nm, " case/control ", obs, "/", ctl, " OK")
    }
  }
  if (length(rna) == 3L) {
    saveRDS(results[rna], file.path(cache_dir, "rnaseq_limma_results.rds"))
    for (nm in rna) {
      obs <- sum(groups[[nm]] == "case"); ctl <- sum(groups[[nm]] == "control")
      stopifnot(obs == expected_case[[nm]], ctl == expected_ctrl[[nm]])
      message("[01]   ", nm, " case/control ", obs, "/", ctl, " OK")
    }
  }
  if (length(gas) == 2L) {
    saveRDS(results[gas], file.path(cache_dir, "cag_limma_results.rds"))
    for (nm in gas) {
      message("[01]   ", nm, " case/control ",
              sum(groups[[nm]] == "case", na.rm = TRUE), "/",
              sum(groups[[nm]] == "control", na.rm = TRUE))
    }
  }

  # group vectors, named and aligned to the matrices
  saveRDS(groups, file.path(cache_dir, "all_groups.rds"))
  for (nm in names(groups)) {
    write_table(data.frame(gsm = names(groups[[nm]]), group = groups[[nm]]),
                paste0("meta_", nm, ".csv"))
  }

  message("[01] done: ", length(results), " cohorts in cache/")
  invisible(results)
}

# -----------------------------------------------------------------------------
# Internal helpers used above
# -----------------------------------------------------------------------------

# Read a raw series-matrix field by line prefix.
.read_field <- function(path, prefix) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  L <- readLines(con, warn = FALSE)
  line <- grep(paste0("^", prefix), L, value = TRUE)[1]
  if (is.na(line)) stop("field not found: ", prefix, call. = FALSE)
  gsub('^"|"$', "", strsplit(line, "\t")[[1]][-1])
}

# Load an annotation package's database object without attaching the package.
get_annotation_db <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("annotation package required but not installed: ", pkg,
         "\n  see README section 3", call. = FALSE)
  }
  getExportedValue(pkg, pkg)
}

# Parse a GEO platform annotation file into an id -> symbol vector.
read_gpl_annotation <- function(path, id_col = "ID", symbol_col = "Symbol") {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  raw <- readLines(con, warn = FALSE)
  i <- grep("^!platform_table_begin", raw)
  e <- grep("^!platform_table_end", raw)
  if (length(i) != 1L || length(e) != 1L || e <= i) {
    stop("platform table markers not found in ", path, call. = FALSE)
  }
  tab <- read.table(text = raw[(i + 1):(e - 1)], header = TRUE, sep = "\t",
                    comment.char = "#", quote = "", check.names = FALSE,
                    stringsAsFactors = FALSE)
  if (!all(c(id_col, symbol_col) %in% names(tab))) {
    stop("expected columns not found; got: ", paste(names(tab), collapse = ", "),
         call. = FALSE)
  }
  setNames(tab[[symbol_col]], tab[[id_col]])
}

# -----------------------------------------------------------------------------
if (!interactive()) run_all(skip_download = FALSE)

# =============================================================================
# 04_assets.R — rebuild the derived objects the later stages consume
# =============================================================================
# STATUS: RESTORATION (from the file the author called handover_assets.R).
#
# This is the bridge between the per-cohort analysis and everything downstream. It
# produces the objects that 05-09 read:
#
#   expr_arrays.rds        log2-scale array matrices, three liver cohorts
#   expr_rnaseq.rds        log2(CPM+1) RNA-seq matrices, three liver cohorts
#   GSE135251_data.rds     counts + log2(CPM+1) + sample metadata for validation
#   all_groups.rds         named case/control vectors, one per cohort
#   MASLD_RRA_{up,dn}_final_ranked.csv   RRA output + per-cohort rank columns
#
# The rank columns matter: the LODO audit in R/05 recomputes the aggregation from
# per-cohort ranks, and needs each candidate's rank in each cohort to do it.
#
# Two normalisation details worth knowing, both preserved from the original:
#   - arrays are taken from the series matrix as supplied (already log2 scale)
#   - RNA-seq counts go through log2(CPM + 1) here for the composition-adjustment
#     analysis; the voom differential expression in R/01 used raw counts
#
# Run:  Rscript R/04_assets.R                (after 01, 02)
# =============================================================================

source("R/00_setup.R")

# -----------------------------------------------------------------------------
# 1. Array matrices
# -----------------------------------------------------------------------------
# Straight from each series matrix, restricted to the samples that survived the
# cohort's grouping rule.

build_arrays <- function() {
  specs <- list(
    Ahrens   = list(file = "GSE48452_series_matrix.txt.gz", fn = group_ahrens),
    Arendt   = list(file = "GSE89632_series_matrix.txt.gz", fn = group_arendt),
    Lefebvre = list(file = "GSE83452_series_matrix.txt.gz", fn = group_lefebvre)
  )
  out <- list()
  for (nm in names(specs)) {
    dir <- file.path(data_dir, nm)
    if (!dir.exists(dir)) { message("[04] SKIP array ", nm); next }
    mat <- read_series_matrix_expr(file.path(dir, specs[[nm]]$file))
    g <- specs[[nm]]$fn(dir)
    mat <- mat[, intersect(names(g), colnames(mat)), drop = FALSE]
    out[[nm]] <- mat
    message("[04] ", nm, ": ", nrow(mat), " probes x ", ncol(mat), " samples")
  }
  out
}

# -----------------------------------------------------------------------------
# 2. RNA-seq matrices at log2(CPM + 1)
# -----------------------------------------------------------------------------
# NOTE: these are the count matrices re-normalised for the composition analysis.
# The voom-based differential expression in R/01 is a separate path and is not
# affected by this function.

log2_cpm <- function(cnt) log2(t(t(cnt) / colSums(cnt) * 1e6) + 1)

build_rnaseq <- function() {
  out <- list()

  dir <- file.path(data_dir, "Suppli")
  if (dir.exists(dir)) {
    cnt <- read.table(gzfile(file.path(dir, "GSE126848_Gene_counts_raw.txt.gz")),
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    sm <- file.path(dir, "GSE126848_series_matrix.txt.gz")
    con <- gzfile(sm, "rt"); L <- readLines(con, warn = FALSE); close(con)
    dsc <- gsub('^"|"$', "", strsplit(grep("^!Sample_description", L, value = TRUE)[1], "\t")[[1]][-1])
    gsm <- gsub('^"|"$', "", strsplit(grep("^!Sample_geo_accession", L, value = TRUE)[1], "\t")[[1]][-1])
    colnames(cnt) <- gsm[match(sub("^0+", "", colnames(cnt)), sub("^0+", "", dsc))]
    out$Suppli <- log2_cpm(cnt)
    message("[04] Suppli: ", nrow(out$Suppli), " genes x ", ncol(out$Suppli), " samples")
  }

  dir <- file.path(data_dir, "Hoang")
  if (dir.exists(dir)) {
    cnt <- read.csv(gzfile(file.path(dir,
                "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz")),
                    row.names = 1, check.names = FALSE)
    ti <- read_titles(file.path(dir, "GSE130970_series_matrix.txt.gz"))
    colnames(cnt) <- ti$gsm[match(colnames(cnt), ti$title)]
    out$Hoang <- log2_cpm(cnt)
    message("[04] Hoang: ", nrow(out$Hoang), " genes x ", ncol(out$Hoang), " samples")
  }

  dir <- file.path(data_dir, "Pantano")
  if (dir.exists(dir)) {
    cnt <- read.csv(gzfile(file.path(dir, "GSE162694_raw_counts.csv.gz")),
                    row.names = 1, check.names = FALSE)
    ti <- read_titles(file.path(dir, "GSE162694_series_matrix.txt"))
    colnames(cnt) <- ti$gsm[match(colnames(cnt), sub(".* ", "", ti$title))]
    out$Pantano <- log2_cpm(cnt)
    message("[04] Pantano: ", nrow(out$Pantano), " genes x ", ncol(out$Pantano), " samples")
  }

  out
}

# -----------------------------------------------------------------------------
# 3. GSE135251 — held-out validation cohort
# -----------------------------------------------------------------------------
# Per-sample count files, filtered to genes detected in >= 5 samples, then
# log2(CPM + 1). The 5-sample filter is part of the original pipeline and is kept.

build_gse135251 <- function() {
  dir <- require_dataset("GSE135251",
                         "needs GSE135251_counts/ (216 per-sample .gz files)")
  cnt_dir <- file.path(dir, "GSE135251_counts")
  if (!dir.exists(cnt_dir)) {
    stop("GSE135251_counts/ not found under ", dir, call. = FALSE)
  }
  fls <- list.files(cnt_dir, pattern = "\\.gz$", full.names = TRUE)
  message("[04] GSE135251 count files: ", length(fls), " (expect 216)")
  if (length(fls) != 216L) warning("[04] expected 216 count files, got ", length(fls))

  gsm_file <- sub("_.*$", "", basename(fls))
  first <- read.delim(gzfile(fls[1]), header = FALSE)
  cnt <- matrix(0, nrow = nrow(first), ncol = length(fls),
                dimnames = list(first[[1]], gsm_file))
  for (i in seq_along(fls)) {
    d <- read.delim(gzfile(fls[i]), header = FALSE)
    stopifnot(identical(d[[1]], rownames(cnt)))
    cnt[, i] <- d[[2]]
  }

  meta <- meta_from_series(file.path(dir, "GSE135251_series_matrix.txt.gz"))
  meta <- meta[match(colnames(cnt), meta$gsm), , drop = FALSE]
  meta$grp <- ifelse(meta$disease == "Control", "control", "case")

  keep <- rowSums(cnt > 0) >= 5
  message("[04] GSE135251 genes kept (>=5 non-zero): ", sum(keep), " of ", nrow(cnt))
  lcp <- log2(t(t(cnt[keep, , drop = FALSE]) /
                  colSums(cnt[keep, , drop = FALSE]) * 1e6) + 1)

  stopifnot(nrow(meta) == 216L,
            sum(meta$grp == "case") == 206L,
            sum(meta$grp == "control") == 10L)
  message("[04] GSE135251: 206 case / 10 control OK")

  list(cnt = cnt, lcp = lcp, meta = meta)
}

# -----------------------------------------------------------------------------
# 4. RRA rank columns
# -----------------------------------------------------------------------------
# Add each candidate's rank in each cohort to the RRA output, so that R/05 can
# recompute the aggregation with a cohort dropped without re-reading the
# expression matrices.

add_rank_columns <- function() {
  chips <- readRDS(file.path(cache_dir, "chips_limma_results.rds"))
  rnas  <- readRDS(file.path(cache_dir, "rnaseq_limma_results.rds"))
  six <- c(chips, rnas)
  if (length(six) != 6L) {
    stop("expected 6 liver cohorts, got ", length(six), call. = FALSE)
  }

  for (direction in c("up", "dn")) {
    path <- file.path(tables_dir, paste0("MASLD_RRA_", direction, "_final.csv"))
    if (!file.exists(path)) {
      stop("Run R/02_RRA_discovery.R first; missing ", basename(path), call. = FALSE)
    }
    csv <- read.csv(path, stringsAsFactors = FALSE)
    for (nm in names(six)) {
      ord <- six[[nm]]$symbol[order(if (direction == "up") -six[[nm]]$t else six[[nm]]$t)]
      csv[[paste0("rank_", nm)]] <- match(csv$Name, ord)
    }
    write_table(csv, paste0("MASLD_RRA_", direction, "_final_ranked.csv"))
    message("[04] added ", length(six), " rank columns to ", direction, " list")
  }
}

# -----------------------------------------------------------------------------
# 5. Run
# -----------------------------------------------------------------------------

arrays <- build_arrays()
if (length(arrays)) {
  saveRDS(arrays, file.path(cache_dir, "expr_arrays.rds"))
  # sample counts the original audit asserted
  for (nm in names(arrays)) {
    if (nm %in% names(expected_samples)) {
      obs <- ncol(arrays[[nm]])
      if (obs != expected_samples[[nm]]) {
        warning("[04] ", nm, ": ", obs, " samples, expected ",
                expected_samples[[nm]], call. = FALSE)
      }
    }
  }
}

rnaseq <- build_rnaseq()
if (length(rnaseq)) saveRDS(rnaseq, file.path(cache_dir, "expr_rnaseq.rds"))

# GSE135251
g135_path <- file.path(cache_dir, "GSE135251_data.rds")
if (dir.exists(file.path(data_dir, "GSE135251"))) {
  d135 <- build_gse135251()
  saveRDS(d135, g135_path)
} else {
  message("[04] GSE135251 not downloaded; skipping (needed by R/07)")
}

add_rank_columns()

# group vectors already written by R/01; assert they are present and named
grp_path <- file.path(cache_dir, "all_groups.rds")
if (file.exists(grp_path)) {
  grp <- readRDS(grp_path)
  for (nm in names(grp)) {
    stopifnot(!is.null(names(grp[[nm]])),
              length(grp[[nm]]) == length(names(grp[[nm]])))
  }
  message("[04] all_groups.rds: ", length(grp), " cohorts, all named OK")
} else {
  warning("[04] all_groups.rds missing - run R/01 first", call. = FALSE)
}

message("[04] done")

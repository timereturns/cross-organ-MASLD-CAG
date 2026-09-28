# =============================================================================
# 09_deconvolution_MuSiC.R — deconvolution sensitivity analysis
# =============================================================================
# STATUS: RESTORATION (from 第 2 步 v4 bulk 性质闸加正确模型 QC, plus the asset
# mapping stage from 正交去卷积敏感性分析).
#
# The question: is the candidate-gene differential expression still there once
# bulk cell-type composition is accounted for? R/06 answered this with marker
# scores. This answers it with a reference-based deconvolution, which is the
# stronger method and the reason both exist.
#
# THE HEADLINE RESULT IS A NEGATIVE ONE, and it is reported as such. Full-gene
# MuSiC fails its pre-registered QC: the snRNA reference and the polyA bulk have
# mismatched gene composition (reference hepatocyte theta_ALB ~0.0044 vs bulk ALB
# fraction 0.16-0.27, a 40-60x discrepancy), and the estimated proportions collapse
# onto hepatocytes -- 0.996 in GSE126848. A marker-restricted fallback is used
# instead, and the full-gene failure is recorded in the QC table rather than
# quietly dropped.
#
# Pipeline:
#   0. bulk-property gate       reject data that are not raw counts
#   1. reference signature      patient x cell-type mean pseudobulk (8 cell types)
#   2. full-gene MuSiC          markers = every shared gene        -> FAILS QC
#   3. marker panel             50 per cell type by specificity; candidates and
#                               MT-/MTRNR genes excluded, absence asserted
#   4. marker-restricted MuSiC  the accepted result
#   5. centered control         full-gene with centered = TRUE     -> degenerates
#   6. reconstruction QC        Y ~ sum_ct p_ct * S_ct * theta_ct
#   7. four-arm comparison      P0 effect, P0 adjusted, this model raw, and this
#                               model adjusted for the MuSiC proportions
#
# TWO DEFINITIONS OF theta IN THE ORIGINAL. The working file computed the reference
# fraction both as the mean of per-sample pseudobulk libraries and as the sum of
# raw counts over all cells, under the same variable name. This script uses the
# SUM OF COUNTS formulation, which is the reference signature as MuSiC defines it.
# See docs/manuscript_discrepancies.md item 15.
#
# ONE THING THIS SCRIPT REPORTS THAT THE ORIGINAL HARDCODED. The original's QC table
# contains several numbers transcribed from a console run log, and the marker panel
# size "400" is asserted only in prose. Here the panel size is printed and written
# out, so the table can be regenerated rather than trusted. See
# docs/manuscript_discrepancies.md items 13 and 14.
#
# Run:  Rscript R/09_deconvolution_MuSiC.R
# =============================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(Biobase)
  library(SingleCellExperiment)
  library(Matrix)
  library(limma)
})

if (!requireNamespace("MuSiC", quietly = TRUE)) {
  stop("MuSiC is required (v1.0.0 from github.com/xuranw/MuSiC); ",
       "see README section 3", call. = FALSE)
}

# The candidates examined here, and the six frozen ones excluded from the marker
# panel. Both are hand-typed literals in the original too.
GENES_TESTED <- c("CADM2", "ANXA4", "LGALS3")
CAND_EXCLUDE <- c("CDHR2", "ANXA4", "CADM2", "IL32", "LGALS3", "RPS6KA1")

# -----------------------------------------------------------------------------
# 0. Inputs and the bulk-property gate
# -----------------------------------------------------------------------------
bulk_path <- file.path(cache_dir, "B_bulks_mapped.rds")
ref_path  <- file.path(cache_dir, "B_reference_merged.rds")

if (!file.exists(bulk_path) || !file.exists(ref_path)) {
  stop("Deconvolution inputs missing from cache/:\n",
       "  B_bulks_mapped.rds      bulk cohorts, symbol-keyed\n",
       "  B_reference_merged.rds  the snRNA reference ($counts, $meta)\n",
       "  Both come from the single-cell preprocessing, which is not part of this\n",
       "  repository. See docs/repository_scope.md.", call. = FALSE)
}

mapped <- readRDS(bulk_path)
stopifnot(is.list(mapped), all(deconvolution_cohorts %in% names(mapped)))
message("[09] bulk cohorts: ", paste(names(mapped), collapse = ", "))

# --- the gate -----------------------------------------------------------------
# MuSiC needs raw counts. If the values look log-scaled the deconvolution is
# meaningless, and it fails quietly rather than loudly, so the gate is explicit
# and it stops rather than warns.
message("[09] ==== bulk property gate ====")
gate_ok <- TRUE
for (q in names(mapped)) {
  B <- mapped[[q]]
  message(sprintf("[09]   %-10s integer=%s  colSums median=%.3g  range [%.3g, %.3g]",
                  q, all(B == round(B)), median(colSums(B)), min(B), max(B)))
  if (max(B) < BULK_GATE_MAX_READS) gate_ok <- FALSE
}
if (!gate_ok) {
  stop("Bulk data look log-scaled, not raw counts (max < ", BULK_GATE_MAX_READS,
       ").\n  Deconvolution would be meaningless. Check how the bulk matrices ",
       "were built.", call. = FALSE)
}
message("[09]   gate passed: values are consistent with raw counts")

# --- reference signature ------------------------------------------------------
ref <- readRDS(ref_path)
ref_counts <- ref$counts
meta <- ref$meta
stopifnot(identical(names(meta), c("Patient.ID", "ct_merged")),
          nrow(meta) == ncol(ref_counts))

message("[09] mean library size per cell, by cell type (QC report only):")
print(round(tapply(Matrix::colSums(ref_counts), meta$ct_merged, mean)))

combo <- paste(meta$Patient.ID, meta$ct_merged, sep = "__")
combo_names <- names(table(combo))
message("[09] patient-celltype columns: ", length(combo_names))

patients_per_type <- tapply(meta$Patient.ID, meta$ct_merged,
                            function(x) length(unique(x)))
message("[09] patients per cell type (require >= 3):")
print(patients_per_type)
stopifnot(all(patients_per_type >= 3))

# Signature: mean over cells, per patient x cell-type column.
ref_sig <- matrix(0, nrow = nrow(ref_counts), ncol = length(combo_names),
                  dimnames = list(rownames(ref_counts), combo_names))
for (i in seq_along(combo_names)) {
  ref_sig[, i] <- Matrix::rowMeans(ref_counts[, combo == combo_names[i], drop = FALSE])
}

parts <- strsplit(combo_names, "__", fixed = TRUE)
ph <- data.frame(Patient.ID = vapply(parts, `[`, character(1), 1),
                 ct_merged  = vapply(parts, `[`, character(1), 2),
                 row.names  = combo_names)
sce <- SingleCellExperiment(assays = list(counts = ref_sig), colData = ph)
CT  <- sort(unique(ph$ct_merged))
message("[09] reference: ", nrow(sce), " genes x ", ncol(sce),
        " patient-celltype columns; ", length(CT), " cell types")
message("[09]   ", paste(CT, collapse = ", "))

# -----------------------------------------------------------------------------
# 1. MuSiC runner
# -----------------------------------------------------------------------------

run_music <- function(mapped, sce, ref_sig, CT, markers_arg = NULL,
                      centered = FALSE, label = "run") {
  prop_list <- list(); r2_list <- list()

  for (q in names(mapped)) {
    B <- mapped[[q]]
    common <- intersect(rownames(ref_sig), rownames(B))
    if (!is.null(markers_arg)) common <- intersect(common, markers_arg)
    if (length(common) < 50) {
      warning("[09] too few shared genes for ", q, " (", length(common), "); skipped",
              call. = FALSE)
      next
    }
    message("[09]   ", label, " / ", q, ": ", length(common), " shared genes")

    est <- MuSiC::music_prop(bulk.mtx = as.matrix(B[common, , drop = FALSE]),
                             sc.sce = sce[common, , drop = FALSE],
                             markers = common,
                             clusters = "ct_merged",
                             samples = "Patient.ID",
                             verbose = FALSE,
                             centered = centered)

    # MuSiC's return shape has changed between versions, so identify the
    # samples x celltypes matrix by shape rather than by name.
    mats <- est[vapply(est, is.matrix, logical(1))]
    hit <- names(mats)[vapply(mats, function(m)
      nrow(m) == ncol(B) && ncol(m) >= length(CT), logical(1))]
    pref <- if ("Est.prop" %in% hit) "Est.prop" else hit[grep("prop", hit)][1]
    if (length(pref) == 0 || is.na(pref)) {
      warning("[09] no proportion matrix returned for ", q, "; skipped", call. = FALSE)
      next
    }
    props <- est[[pref]][colnames(B), , drop = FALSE]
    props <- props[, intersect(CT, colnames(props)), drop = FALSE]
    prop_list[[q]] <- props

    r2name <- grep("r.squared", names(est), value = TRUE)[1]
    if (length(r2name) && !is.na(r2name)) {
      rv <- est[[r2name]]
      if (is.matrix(rv)) rv <- rv[, 1]
      r2_list[[q]] <- if (is.numeric(rv) && length(rv) == ncol(B)) rv else
        rep(NA_real_, ncol(B))
    } else {
      r2_list[[q]] <- rep(NA_real_, ncol(B))
    }
    message("[09]     mean proportions: ",
            paste(sprintf("%s=%.3f", colnames(props), colMeans(props)), collapse = " "))
  }
  list(prop = prop_list, r2 = r2_list)
}

# --- 1a. full-gene run (expected to fail QC) ---------------------------------
message("[09] ==== full-gene MuSiC (markers = all shared genes) ====")
full <- run_music(mapped, sce, ref_sig, CT, markers_arg = NULL, label = "full-gene")

# -----------------------------------------------------------------------------
# 2. Marker panel
# -----------------------------------------------------------------------------
# Specificity = a gene's fraction within a cell type, divided by its mean fraction
# across cell types. Top 50 per cell type.
#
# Candidates and mitochondrial genes are zeroed BEFORE ranking, and their absence
# from the panel is asserted. If a candidate were itself a marker of the cell type
# it is being tested in, the adjustment would be circular and the attenuation
# meaningless. This is the single most important safeguard in the analysis.

message("[09] ==== marker panel ====")

# theta: reference signature summed over cells of each type, then column-normalised
sig_ct <- do.call(cbind, lapply(CT, function(ct) {
  as.numeric(Matrix::rowSums(ref_counts[, meta$ct_merged == ct, drop = FALSE]))
}))
rownames(sig_ct) <- rownames(ref_counts)
colnames(sig_ct) <- CT
theta <- sweep(sig_ct, 2, colSums(sig_ct), "/")

spec0 <- theta / (rowMeans(theta) + 1e-12)

excl_nom <- unique(c(CAND_EXCLUDE,
                     grep("^MT-|^MTRNR", rownames(theta), value = TRUE)))
excl_in   <- excl_nom[excl_nom %in% rownames(theta)]
excl_skip <- setdiff(excl_nom, rownames(theta))
message("[09] exclusions: ", length(excl_nom), " named, ", length(excl_in),
        " present and zeroed, ", length(excl_skip), " absent from the reference")
if (length(excl_skip)) {
  # RPS6KA1 is expected here: it is not in the reference at all, so no
  # reference-side adjustment exists for it.
  message("[09]   absent: ", paste(excl_skip, collapse = ", "))
}

spec <- spec0
spec[excl_in, ] <- 0

markers <- unique(unlist(lapply(CT, function(ct)
  rownames(theta)[order(spec[, ct], decreasing = TRUE)[1:MUSIC_MARKERS_PER_CELLTYPE]])))

panel_size <- length(markers)
message("[09] marker panel: ", panel_size, " genes (",
        MUSIC_MARKERS_PER_CELLTYPE, " per cell type before deduplication; ",
        "manuscript states ", MUSIC_MARKER_PANEL_MAX, ")")
print(table(unlist(lapply(CT, function(ct)
  rownames(theta)[order(spec[, ct], decreasing = TRUE)[1:MUSIC_MARKERS_PER_CELLTYPE]]))))

# The anti-circularity assertion.
in_panel <- CAND_EXCLUDE[CAND_EXCLUDE %in% markers]
if (length(in_panel)) {
  stop("Candidate gene(s) present in the marker panel: ",
       paste(in_panel, collapse = ", "),
       "\n  The composition adjustment would be circular.", call. = FALSE)
}
message("[09] anti-circularity check passed: no candidate is in the marker panel")

write_table(data.frame(
  metric = c("panel_size", "per_celltype_target", "cell_types", "manuscript_states"),
  value = c(panel_size, MUSIC_MARKERS_PER_CELLTYPE, length(CT),
            MUSIC_MARKER_PANEL_MAX)),
  "B_music_marker_panel_size.csv")

# --- 2b. marker-restricted run (the accepted result) -------------------------
message("[09] ==== marker-restricted MuSiC (the accepted result) ====")
marker_run <- run_music(mapped, sce, ref_sig, CT, markers_arg = markers,
                        label = "marker-panel")

# --- 2c. centered control ----------------------------------------------------
message("[09] ==== centered control (full-gene, centered = TRUE) ====")
centered_run <- run_music(mapped, sce, ref_sig, CT, markers_arg = NULL,
                          centered = TRUE, label = "centered")

# -----------------------------------------------------------------------------
# 3. Reconstruction QC
# -----------------------------------------------------------------------------
# The model being checked is  Y ~ sum_ct p_ct * S_ct * theta_ct  on relative
# abundances. Reconstruction correlation and relative residual are reported per
# sample; the original used them to decide that the full-gene arm was unusable.

reconstruction_qc <- function(prop_list, label) {
  rows <- list()
  for (q in names(prop_list)) {
    P <- prop_list[[q]]
    B <- mapped[[q]]
    common <- intersect(rownames(ref_sig), rownames(B))
    Y <- sweep(B[common, , drop = FALSE], 2, colSums(B[common, , drop = FALSE]), "/")
    sig_frac <- sweep(sig_ct[common, CT, drop = FALSE], 2,
                      colSums(sig_ct[common, CT, drop = FALSE]), "/")
    Yhat <- sig_frac %*% t(P[, CT, drop = FALSE])
    cor_p <- sapply(seq_len(ncol(Y)), function(i) cor(Y[, i], Yhat[, i]))
    cor_s <- sapply(seq_len(ncol(Y)), function(i)
      cor(Y[, i], Yhat[, i], method = "spearman"))
    rel_resid <- sapply(seq_len(ncol(Y)), function(i)
      sqrt(mean((Y[, i] - Yhat[, i])^2) / mean(Y[, i]^2)))
    rows[[q]] <- data.frame(
      arm = label, cohort = q, sample = rownames(P),
      row_sum = rowSums(P), max_abs_dev_1 = abs(rowSums(P) - 1),
      n_neg = rowSums(P < -1e-8), min_prop = apply(P, 1, min),
      rec_cor_pearson = cor_p, rec_cor_spearman = cor_s,
      rel_resid = rel_resid, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

qc <- rbind(reconstruction_qc(full$prop, "full-gene"),
            reconstruction_qc(marker_run$prop, "marker-panel"))

message("[09] ==== QC ====")
message("[09] row-sum max deviation: ", signif(max(qc$max_abs_dev_1), 3))
message("[09] samples with negative proportions: ", sum(qc$n_neg > 0))
message("[09] reconstruction Pearson, median by arm and cohort:")
print(round(tapply(qc$rec_cor_pearson, list(qc$arm, qc$cohort), median), 3))
message("[09] relative residual, median by arm and cohort:")
print(round(tapply(qc$rel_resid, list(qc$arm, qc$cohort), median), 3))
write_table(qc, "B_music_qc.csv")

# --- the negative result, stated in numbers ----------------------------------
message("[09] ==== pre-registered QC verdict ====")
hep <- sapply(full$prop, function(P) mean(P[, "Hepatocytes"]))
message("[09] full-gene mean hepatocyte proportion (collapse indicator):")
print(round(hep, 3))
collapsed <- names(hep)[hep > 0.9]
message("[09] cohorts where full-gene proportions collapse onto hepatocytes: ",
        if (length(collapsed)) paste(collapsed, collapse = ", ") else "none")

# ALB discordance: the mechanism behind the collapse. Computed, not transcribed.
th_alb <- theta["ALB", "Hepatocytes"]
alb_rows <- do.call(rbind, lapply(names(mapped), function(q) {
  B <- mapped[[q]]
  f <- median(B["ALB", ] / colSums(B))
  data.frame(cohort = q, reference_theta_ALB = round(th_alb, 4),
             bulk_ALB_median = round(f, 4),
             ratio = round(f / th_alb, 1), stringsAsFactors = FALSE)
}))
message("[09] ALB discordance between reference and bulk:")
print(alb_rows)
write_table(alb_rows, "B_music_alb_discordance.csv")

centered_hep <- sapply(centered_run$prop, function(P)
  if ("Hepatocytes" %in% colnames(P)) mean(P[, "Hepatocytes"]) else NA_real_)
message("[09] centered control mean hepatocyte proportion: ",
        paste(sprintf("%s=%.3f", names(centered_hep), centered_hep), collapse = " "))

# -----------------------------------------------------------------------------
# 4. Four-arm comparison
# -----------------------------------------------------------------------------
# Arm 1  P0 per-cohort effect size (the discovery screen's estimate)
# Arm 2  P0 composition-adjusted summary (from R/06)
# Arm 3  this script's own model, unadjusted
# Arm 4  this script's own model, adjusted for the MuSiC proportions
#
# NOTE the identifiers: the original named these arm1_/arm2_/arm3_ and "arm 4" is
# arm3_logFC_adj, so the fourth column set carries an arm3 name. Kept as-is so the
# output lines up with the published Table S32, with a comment rather than a
# silent rename.

message("[09] ==== four-arm comparison ====")

es_path <- file.path(tables_dir, "p0_cohort_effect_sizes.csv")
mk_path <- file.path(tables_dir, "p0_composition_summary.csv")
if (!file.exists(es_path) || !file.exists(mk_path)) {
  stop("Arms 1 and 2 need p0_cohort_effect_sizes.csv (R/05) and ",
       "p0_composition_summary.csv (R/06).", call. = FALSE)
}
es <- read.csv(es_path, stringsAsFactors = FALSE)
mk <- read.csv(mk_path, stringsAsFactors = FALSE)

bulk_groups <- file.path(cache_dir, "all_groups.rds")
if (!file.exists(bulk_groups)) {
  stop("all_groups.rds needed for the group vectors. Run R/01.", call. = FALSE)
}
all_grp <- readRDS(bulk_groups)

# Endpoints differ by cohort: three are case/control, GSE174478 is a fibrosis
# stage. Arms 1 and 2 do not exist for GSE174478, because it was not part of P0.
fib_stage <- NULL
fib_path <- file.path(cache_dir, "GSE174478_fibrosis_stage.rds")

b4_rows <- list()
for (q in names(marker_run$prop)) {
  B <- mapped[[q]]
  P <- marker_run$prop[[q]]
  mu <- colMeans(P)

  # covariates: the two structural cell types always, plus any above the
  # proportion threshold; zero-variance columns dropped
  cov_ct <- union(MUSIC_FORCED_COVARIATES,
                  names(mu)[mu >= MUSIC_COVARIATE_MIN_PROP])
  Z <- P[, intersect(cov_ct, colnames(P)), drop = FALSE]
  Z <- Z[, apply(Z, 2, var) > 1e-6, drop = FALSE]

  endpoint <- NULL
  if (q %in% names(all_grp)) {
    g <- all_grp[[q]]
    ok <- intersect(names(g), colnames(B))
    y <- as.numeric(factor(g[ok])) - 1
    X <- B[, ok, drop = FALSE]
    endpoint <- "case_control"
    des_raw <- cbind(Int = 1, group = y)
  } else if (q == "GSE174478" && !is.null(fib_stage)) {
    y <- fib_stage[colnames(B)]
    X <- B
    endpoint <- "fibrosis_stage"
    des_raw <- cbind(Int = 1, endpoint = as.numeric(y))
  } else {
    message("[09]   ", q, ": no endpoint available; skipped from the four-arm table")
    next
  }

  have <- GENES_TESTED[GENES_TESTED %in% rownames(X)]
  if (!length(have)) { message("[09]   ", q, ": no candidates present; skipped"); next }

  Zc <- Z[colnames(X), , drop = FALSE]
  des_adj <- cbind(des_raw, Zc)
  fit_raw <- eBayes(lmFit(X[have, , drop = FALSE], des_raw))
  fit_adj <- eBayes(lmFit(X[have, , drop = FALSE], des_adj))

  for (gn in have) {
    e1 <- es[es$gene == gn & es$cohort == q, ]
    e2 <- mk[mk$gene == gn & mk$cohort == q, ]
    raw <- fit_raw$coef[gn, 2]; adj <- fit_adj$coef[gn, 2]
    att <- if (abs(raw) < 1e-8) NA_real_ else 1 - adj / raw
    b4_rows[[length(b4_rows) + 1L]] <- data.frame(
      cohort = q, gene = gn, endpoint = endpoint, n_samples = ncol(X),
      arm1_logFC_p0 = if (nrow(e1)) e1$logFC[1] else NA_real_,
      arm1_p_p0 = if (nrow(e1)) e1$P[1] else NA_real_,
      arm2_logFC_raw = if (nrow(e2)) e2$logFC_raw[1] else NA_real_,
      arm2_logFC_adj = if (nrow(e2)) e2$logFC_adj[1] else NA_real_,
      arm2_p_adj = if (nrow(e2)) e2$p_adj[1] else NA_real_,
      arm3_logFC_raw = raw, arm3_p_raw = fit_raw$p.value[gn, 2],
      arm3_logFC_adj = adj,   arm3_p_adj = fit_adj$p.value[gn, 2],
      att_music = att,
      dir_keep_music = sign(adj) == sign(raw),
      stringsAsFactors = FALSE)
  }
}

b4 <- do.call(rbind, b4_rows)
write_table(b4, "B_b4_compare.csv")

# NOTE: the original figure script appended `verdict` and `note` columns to this
# CSV by hand and wrote it back in place. Those two columns are curator
# annotations, not computed values. They are NOT recreated here; the published
# Table S32 contains them and the repository documents them as hand-entered.
# See docs/manuscript_discrepancies.md item 12.

message("[09] four-arm table: ", nrow(b4), " cohort-gene rows")
print(b4[, c("cohort", "gene", "arm1_logFC_p0", "arm3_logFC_raw",
             "arm3_logFC_adj", "att_music", "dir_keep_music")], digits = 3)

# -----------------------------------------------------------------------------
# 5. Persist
# -----------------------------------------------------------------------------
saveRDS(list(full = full, marker = marker_run, centered = centered_run,
             markers = markers, theta = theta, spec = spec,
             cell_types = CT, ref_sig = ref_sig),
        file.path(cache_dir, "music_results.rds"), compress = "xz")

# Long-format proportions of the accepted (marker-restricted) arm.
long <- do.call(rbind, lapply(names(marker_run$prop), function(q) {
  p <- marker_run$prop[[q]]
  data.frame(cohort = q,
             sample = rep(rownames(p), times = ncol(p)),
             cell_type = rep(colnames(p), each = nrow(p)),
             prop = as.vector(p), stringsAsFactors = FALSE)
}))
write_table(long, "B_music_marker_props.csv")

message("[09] done")

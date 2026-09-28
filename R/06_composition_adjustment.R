# =============================================================================
# 06_composition_adjustment.R — cell-composition correction and batch diagnostics
# =============================================================================
# STATUS: RESTORATION (from the file the author called P0_composition_batch_v2.R).
#
# The question this answers: is the candidate-gene differential expression just a
# reflection of changed cell-type composition? Six liver cohorts change their
# hepatocyte/immune/stellate balance in disease, and a bulk difference can follow
# from that alone.
#
# The approach is deliberately cheap and transparent: score each cell type from a
# small curated marker panel, add the three most variable scores as covariates to
# the limma model, and ask how much of the effect survives.
#
#   1. Marker-score composition  z-score each gene across samples, average the
#      markers per cell type. Requires >= 2 markers present.
#   2. Covariates  the top 3 cell types by variance across samples. Variance, not
#      group difference, so the choice is not informed by the outcome.
#   3. Two models  raw (~ group) and adjusted (~ group + score1 + score2 + score3).
#   4. Attenuation  1 - logFC_adj / logFC_raw, plus a direction-flip flag.
#   5. Batch diagnostic  PCA on variance-filtered genes, first 3 PCs regressed on
#      group, so that a cohort whose disease effect lives in PC1 is visible.
#
# TWO THINGS TO KNOW ABOUT THIS SCRIPT
#
#   * It uses NINE liver cell types (see mk_liver in R/00_setup.R). The MuSiC
#     deconvolution in R/09 uses EIGHT. Both are correct for their own analysis,
#     but the manuscript's single "8 cell types" phrase does not cover this one -
#     see docs/manuscript_discrepancies.md item 4.
#
#   * Composition scores are correlations of z-scored marker averages, i.e. a
#     surrogate. This is why the MuSiC arm exists as an orthogonal check rather
#     than a repeat.
#
# Run:  Rscript R/06_composition_adjustment.R   (after 04)
# =============================================================================

source("R/00_setup.R")

suppressPackageStartupMessages(library(limma))

# Candidate genes examined here. The original used six, not fourteen: the five
# retained genes plus the observation-list gene.
cand <- c("CDHR2", "ANXA4", "CADM2", "IL32", "LGALS3", "RPS6KA1")

# -----------------------------------------------------------------------------
# 1. Load matrices and groups
# -----------------------------------------------------------------------------
arr_path <- file.path(cache_dir, "expr_arrays.rds")
rna_path <- file.path(cache_dir, "expr_rnaseq.rds")
grp_path <- file.path(cache_dir, "all_groups.rds")

if (!file.exists(grp_path)) {
  stop("all_groups.rds missing. Run R/01 first.", call. = FALSE)
}
grp <- readRDS(grp_path)

arrays <- if (file.exists(arr_path)) readRDS(arr_path) else list()
rnaseq <- if (file.exists(rna_path)) readRDS(rna_path) else list()
mats <- c(arrays, rnaseq)

if (!length(mats)) {
  stop("No expression matrices in cache/. Run R/04_assets.R first.", call. = FALSE)
}
message("[06] cohorts with matrices: ", paste(names(mats), collapse = ", "))

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------

# Group vector, aligned to a matrix's columns, with the balance asserted.
get_group <- function(g_obj, samples, nm) {
  v <- g_obj[samples]
  if (anyNA(v)) stop(nm, ": GSM names do not match matrix columns", call. = FALSE)
  message("[06]   ", nm, ": ", paste(names(table(v)), table(v), collapse = " / "))
  if (nm %in% names(expected_case)) {
    stopifnot(sum(v == "case") == expected_case[[nm]],
              sum(v == "control") == expected_ctrl[[nm]])
  }
  factor(ifelse(v == "case", "case", "ctrl"), levels = c("ctrl", "case"))
}

# Marker-score composition. Genes are z-scored across samples, then averaged per
# cell type. Cell types with fewer than two markers present return NA rather than
# a single-marker average, which would be noise.
score_mk <- function(mat, mk) {
  rn <- toupper(rownames(mat))
  z  <- t(scale(t(mat)))
  sapply(names(mk), function(ct) {
    g <- mk[[ct]][toupper(mk[[ct]]) %in% rn]
    if (length(g) < 2) rep(NA_real_, ncol(mat))
    else colMeans(z[g, , drop = FALSE], na.rm = TRUE)
  })
}

# Marker hit diagnostics: if a platform is probe-level and unmapped, every cell
# type returns NA and the whole analysis is vacuous. Report before proceeding.
report_marker_hits <- function(mat, mk, nm) {
  mm <- sapply(mk, function(gs) sum(toupper(gs) %in% toupper(rownames(mat))))
  message("[06]   ", nm, " marker hits: ",
          paste0(names(mm), "=", mm, collapse = ", "))
  if (sum(mm >= 2) < 3) {
    warning("[06] ", nm, ": fewer than 3 cell types have >= 2 markers; ",
            "the adjustment for this cohort is not meaningful", call. = FALSE)
  }
  invisible(mm)
}

# -----------------------------------------------------------------------------
# 3. Per-cohort composition adjustment
# -----------------------------------------------------------------------------
res  <- list()
summ_rows <- list()

for (nm in names(mats)) {
  if (!nm %in% names(grp)) { warning("[06] ", nm, " has no group vector; skipped"); next }
  M <- mats[[nm]]
  g <- get_group(grp[[nm]], colnames(M), nm)
  report_marker_hits(M, mk_liver, nm)

  sc <- score_mk(M, mk_liver)
  ok <- apply(sc, 2, function(x) !all(is.na(x)))
  sc <- sc[, ok, drop = FALSE]
  if (ncol(sc) < 1) { warning("[06] ", nm, ": no usable composition scores; skipped"); next }

  vv <- apply(sc, 2, var, na.rm = TRUE)
  top3 <- names(sort(vv, decreasing = TRUE, na.last = NA))[seq_len(min(3, ncol(sc)))]
  message("[06]   ", nm, " covariates: ", paste(top3, collapse = ", "))

  # (a) composition differs between groups?
  comp_p <- apply(sc, 2, function(x) wilcox.test(x ~ g)$p.value)

  # (b) candidate expression vs each composition score
  have <- cand[cand %in% toupper(rownames(M))]
  corr <- sapply(have, function(gn) {
    x <- M[toupper(rownames(M)) == gn, , drop = FALSE][1, ]
    sapply(colnames(sc), function(ct)
      cor(x, sc[, ct], method = "spearman", use = "complete.obs"))
  })

  # (c) raw vs adjusted model
  des_raw <- model.matrix(~ g)
  des_adj <- if (length(top3) >= 1) {
    model.matrix(as.formula(paste("~ g +",
      paste0("sc[, top3[", seq_along(top3), "]]", collapse = " + "))))
  } else des_raw
  fit_raw <- eBayes(lmFit(M[have, , drop = FALSE], des_raw))
  fit_adj <- eBayes(lmFit(M[have, , drop = FALSE], des_adj))

  # (d) batch diagnostic. Variance-filter first: a constant gene breaks prcomp
  # under scaling, and the original had to patch this.
  Mfull <- M
  keep <- apply(Mfull, 1, var) > 0
  message("[06]   ", nm, " PCA genes kept: ", sum(keep), "/", nrow(Mfull))
  pc <- prcomp(t(Mfull[keep, , drop = FALSE]), center = TRUE, scale. = TRUE)$x[, 1:3, drop = FALSE]
  pc_group_p <- apply(pc, 2, function(x) summary(lm(x ~ g))$coef[2, 4])

  res[[nm]] <- list(comp_diff_p = comp_p, corr = corr, top3 = top3,
                    raw_logFC = fit_raw$coef[, 2], raw_p = fit_raw$p.value[, 2],
                    adj_logFC = fit_adj$coef[, 2], adj_p = fit_adj$p.value[, 2],
                    pc_group_p = pc_group_p)
}

saveRDS(res, file.path(cache_dir, "composition_res.rds"))

# -----------------------------------------------------------------------------
# 4. Summary table
# -----------------------------------------------------------------------------
for (nm in names(res)) {
  r <- res[[nm]]
  C <- r$corr
  gn_list <- names(r$raw_logFC)
  for (gn in gn_list) {
    # corr may be gene x celltype or celltype x gene depending on how many
    # candidates were present; identify the orientation rather than assume it.
    if (gn %in% rownames(C))      { v <- C[gn, ]; ct_names <- colnames(C) }
    else if (gn %in% colnames(C)) { v <- C[, gn]; ct_names <- rownames(C) }
    else next
    topct <- ct_names[which.max(abs(v))]

    raw <- r$raw_logFC[[gn]]
    adj <- r$adj_logFC[[gn]]
    # Attenuation is undefined if the raw effect is ~0: report NA rather than Inf.
    att <- if (abs(raw) < 1e-8) NA_real_ else 1 - adj / raw

    summ_rows[[length(summ_rows) + 1L]] <- data.frame(
      cohort = nm, gene = gn,
      top3_cov = paste(r$top3, collapse = ";"),
      logFC_raw = round(raw, 3),
      logFC_adj = round(adj, 3),
      p_raw = signif(r$raw_p[[gn]], 2),
      p_adj = signif(r$adj_p[[gn]], 2),
      attenuation = round(att, 2),
      direction_flip = sign(raw) != sign(adj),
      top_cor_celltype = topct,
      top_cor = round(v[which.max(abs(v))], 3),
      pc1_p = signif(r$pc_group_p[1], 2),
      pc2_p = signif(r$pc_group_p[2], 2),
      pc3_p = signif(r$pc_group_p[3], 2),
      stringsAsFactors = FALSE)
  }
}
summ <- do.call(rbind, summ_rows)
write_table(summ, "p0_composition_summary.csv")

# Composition difference table
comp_diff <- do.call(rbind, lapply(names(res), function(nm) {
  data.frame(cohort = nm,
             celltype = names(res[[nm]]$comp_diff_p),
             wilcox_p = signif(res[[nm]]$comp_diff_p, 2),
             stringsAsFactors = FALSE)
}))
write_table(comp_diff, "p0_composition_compdiff.csv")

# --- interpretation ----------------------------------------------------------
# Report the headline: how many candidate-cohort pairs flip direction after
# adjustment, and how many attenuate by more than half. Both are the numbers a
# reviewer will want.
n_flip <- sum(summ$direction_flip, na.rm = TRUE)
n_att50 <- sum(summ$attenuation > 0.5, na.rm = TRUE)
message("[06] candidate-cohort pairs: ", nrow(summ))
message("[06]   direction flips after adjustment: ", n_flip)
message("[06]   attenuated by more than 50%:      ", n_att50)
message("[06]   cohort PC1 associated with group:  ",
        sum(summ$pc1_p < 0.05, na.rm = TRUE), " of ", nrow(summ), " rows")

print(summ[order(-abs(summ$attenuation)), c("cohort", "gene", "logFC_raw",
      "logFC_adj", "attenuation", "direction_flip", "top_cor_celltype")][1:min(15, nrow(summ)), ])
message("[06] done")

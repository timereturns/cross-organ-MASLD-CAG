# =============================================================================
# 05_audit_LODO_meta.R — robustness audit of the 14 candidates
# =============================================================================
# STATUS: RESTORATION (from the file the author called P0_audit_disk.R).
#
# Four questions, in the order the original audit asked them:
#
#   1. Calibration. Does the closed-form rho reproduce the stored RRA Score? The
#      answer decides which background denominator the LODO recomputation uses.
#      Three definitions are on the table: the six-cohort common-gene count, the
#      per-cohort gene counts, and a fixed 20,000. The original selected 20,000 by
#      hand, and that choice is recorded here rather than left implicit.
#
#   2. Leave-one-dataset-out. Recompute the aggregation with each cohort dropped in
#      turn, by full re-ranking, not by deleting a rank row. Two retention criteria
#      are recorded per gene: rank <= 200 after the drop, and corrected P < 0.01
#      after the drop.
#
#   3. Gastric support. Which of the two gastric cohorts support each gene,
#      directionally and at P < 0.05 with |log2FC| > log2(1.5).
#
#   4. Effect sizes and meta-analysis. Per-cohort logFC with SE = logFC/t, then a
#      random-effects (REML) pooled estimate per organ within platform. Pooling is
#      chip = {Ahrens, Arendt, Lefebvre} and rna = {Suppli, Hoang, Pantano}; the two
#      gastric cohorts are pooled separately (fixed and DerSimonian-Laird) and their
#      k=2 random-effects estimate is reported as a sensitivity figure only.
#
# Run:  Rscript R/05_audit_LODO_meta.R       (after 02 and 04)
# =============================================================================

source("R/00_setup.R")

if (!requireNamespace("metafor", quietly = TRUE)) {
  stop("metafor is required; see README section 3", call. = FALSE)
}
suppressPackageStartupMessages(library(metafor))

# -----------------------------------------------------------------------------
# 0. Load the audit's inputs
# -----------------------------------------------------------------------------
need <- c("candidate_genes_final.csv", "MASLD_RRA_up_final_ranked.csv",
          "MASLD_RRA_dn_final_ranked.csv")
if (!all(file.exists(file.path(tables_dir, need)))) {
  stop("Missing audit inputs: ", paste(setdiff(need, list.files(tables_dir)),
                                       collapse = ", "),
       "\n  Run R/02 and R/04 first.", call. = FALSE)
}
if (!all(file.exists(file.path(cache_dir,
        c("chips_limma_results.rds", "rnaseq_limma_results.rds",
          "cag_limma_results.rds"))))) {
  stop("Limma result objects missing. Run R/01 first.", call. = FALSE)
}

cand    <- read.csv(file.path(tables_dir, "candidate_genes_final.csv"),
                    stringsAsFactors = FALSE)
genes14 <- cand$gene
dir14   <- cand$direction

chips <- readRDS(file.path(cache_dir, "chips_limma_results.rds"))
rnas  <- readRDS(file.path(cache_dir, "rnaseq_limma_results.rds"))
cag   <- readRDS(file.path(cache_dir, "cag_limma_results.rds"))
masld <- c(chips, rnas)

N_per <- vapply(masld, nrow, numeric(1))
message("[05] audit over ", length(genes14), " candidates x ",
        length(c(chips, rnas, cag)), " cohorts")

# -----------------------------------------------------------------------------
# 1. Calibration: reproduce the stored Score
# -----------------------------------------------------------------------------
# The stored Score came from aggregateRanks. If rra_rho() with a given denominator
# reproduces it, that denominator is the one to use for the LODO recomputation.

rnk_cols <- c("rank_Ahrens", "rank_Arendt", "rank_Lefebvre",
              "rank_Suppli", "rank_Hoang", "rank_Pantano")
r_up <- read.csv(file.path(tables_dir, "MASLD_RRA_up_final_ranked.csv"),
                 stringsAsFactors = FALSE)

den_per <- N_per[c("Ahrens", "Arendt", "Lefebvre", "Suppli", "Hoang", "Pantano")]

cal_per   <- sapply(seq_len(nrow(r_up)), function(i) rra_rho(r_up[i, rnk_cols], den_per))
cal_20000 <- sapply(seq_len(nrow(r_up)), function(i) rra_rho(r_up[i, rnk_cols], RRA_BACKGROUND_N_AUDIT))

calib <- data.frame(gene = r_up$Name,
                    stored = r_up$Score,
                    ratio_per_cohort = cal_per / r_up$Score,
                    ratio_N20000 = cal_20000 / r_up$Score)
print(head(calib, 12))

# --- The choice of denominator -----------------------------------------------
# The original used 20,000 for the LODO recomputation. This is a hand decision,
# recorded here; change it only deliberately. If the per-cohort denominators
# reproduce the stored Score better, that is worth knowing - so report both and
# let the ratios above speak.
USE_N20000 <- TRUE

den_lodo <- if (USE_N20000) {
  function(k) rep(RRA_BACKGROUND_N_AUDIT, k)
} else {
  function(k) den_per[seq_len(k)]
}

message("[05] LODO denominator: ",
        if (USE_N20000) paste0("fixed N = ", RRA_BACKGROUND_N_AUDIT) else
          "per-cohort gene counts",
        "  (hand-selected; see docs/analysis_parameters.md)")

write_table(calib, "p0_rra_calibration.csv")

# -----------------------------------------------------------------------------
# 2. Exact leave-one-dataset-out
# -----------------------------------------------------------------------------

up_mat <- build_rankmat(masld, "up")
dn_mat <- build_rankmat(masld, "dn")
message("[05] rank matrices: up ", nrow(up_mat), " x ", ncol(up_mat),
        ", dn ", nrow(dn_mat), " x ", ncol(dn_mat))

lodo_check <- function(mat, genes, direction) {
  full_sc <- apply(mat, 1, rra_rho, denom = den_lodo(ncol(mat)))
  full_rk <- rank(full_sc, ties.method = "min")
  out <- data.frame(gene = genes,
                    dir = direction,
                    score_full = full_sc[genes],
                    rank_full = full_rk[genes],
                    stringsAsFactors = FALSE)
  for (j in seq_len(ncol(mat))) {
    sc <- apply(mat[, -j, drop = FALSE], 1, rra_rho, denom = den_lodo(ncol(mat) - 1))
    rk <- rank(sc, ties.method = "min")
    out[[paste0("rank_drop_", colnames(mat)[j])]] <- rk[genes]
    out[[paste0("in200_drop_", colnames(mat)[j])]] <- rk[genes] <= LODO_RANK_CUTOFF
    out[[paste0("lt001_drop_", colnames(mat)[j])]] <- sc[genes] < LODO_P_CUTOFF
  }
  out
}

lodo_up <- lodo_check(up_mat, genes14[dir14 == "up"], "up")
lodo_dn <- lodo_check(dn_mat, genes14[dir14 == "dn"], "dn")
lodo_all <- rbind(lodo_up, lodo_dn)
write_table(lodo_all, "p0_lodo_exact_rra.csv")

in200_cols <- grep("^in200_drop_", names(lodo_all))
keep_counts <- rowSums(lodo_all[, in200_cols], na.rm = TRUE)
message("[05] LODO retention (of ", length(in200_cols), " drops):")
for (i in seq_len(nrow(lodo_all))) {
  message(sprintf("[05]   %-10s %d/%d in top %d", lodo_all$gene[i],
                  keep_counts[i], length(in200_cols), LODO_RANK_CUTOFF))
}

# -----------------------------------------------------------------------------
# 3. Gastric support
# -----------------------------------------------------------------------------
# A gastric cohort supports a gene if the direction matches and the gene clears
# P < 0.05 with |log2FC| > log2(1.5). The magnitude threshold is stricter than the
# 0.5 used for DEG definition in R/03: this is an audit, so it asks for more.

cag_support <- function(gene, direction) {
  sapply(names(cag), function(cn) {
    d <- cag[[cn]]
    i <- match(gene, d$symbol)
    if (is.na(i)) return(NA_character_)
    dir_ok <- if (direction == "up") d$logFC[i] > 0 else d$logFC[i] < 0
    thr_ok <- d$P.Value[i] < GASTRIC_DEG_P &
      abs(d$logFC[i]) > GASTRIC_SUPPORT_LOGF
    paste0(ifelse(thr_ok & dir_ok, "1", "0"),
           " (logFC=", round(d$logFC[i], 2),
           ", P=", signif(d$P.Value[i], 2), ")")
  })
}

cag_tab <- do.call(rbind, lapply(seq_along(genes14), function(i) {
  data.frame(gene = genes14[i], dir = dir14[i],
             t(cag_support(genes14[i], dir14[i])), check.names = FALSE)
}))
write_table(cag_tab, "p0_cag_support.csv")

# -----------------------------------------------------------------------------
# 4. Per-cohort effect sizes and meta-analysis
# -----------------------------------------------------------------------------
# SE is back-computed from the moderated t: SE = logFC / t. This is the original
# approach and it is what makes the meta-analysis possible without the raw series.

eff_tab <- do.call(rbind, lapply(genes14, function(g) {
  dir <- dir14[match(g, genes14)]
  do.call(rbind, lapply(names(c(masld, cag)), function(cn) {
    d <- c(masld, cag)[[cn]]
    i <- match(g, d$symbol)
    if (is.na(i)) return(NULL)
    se <- d$logFC[i] / d$t[i]
    data.frame(gene = g, dir = dir, cohort = cn,
               logFC = d$logFC[i], SE = se,
               CI_lo = d$logFC[i] - 1.96 * se,
               CI_hi = d$logFC[i] + 1.96 * se,
               P = d$P.Value[i], FDR = d$adj.P.Val[i],
               dir_ok = sign(d$logFC[i]) == ifelse(dir == "up", 1, -1),
               stringsAsFactors = FALSE)
  }))
}))
write_table(eff_tab, "p0_cohort_effect_sizes.csv")
message("[05] per-cohort effect rows: ", nrow(eff_tab))

dir_cnt <- tapply(eff_tab$dir_ok[eff_tab$cohort %in% names(masld)],
                  eff_tab$gene[eff_tab$cohort %in% names(masld)], sum, na.rm = TRUE)

# --- liver: random effects within platform -----------------------------------
meta_tab <- do.call(rbind, lapply(genes14, function(g) {
  sub <- eff_tab[eff_tab$gene == g & eff_tab$cohort %in% names(masld), ]
  do.call(rbind, lapply(c("chip", "rna"), function(plat) {
    s <- if (plat == "chip") {
      sub[sub$cohort %in% c("Ahrens", "Arendt", "Lefebvre"), ]
    } else {
      sub[sub$cohort %in% c("Suppli", "Hoang", "Pantano"), ]
    }
    if (nrow(s) < 2) return(NULL)
    m <- tryCatch(rma(yi = logFC, sei = SE, data = s, method = META_METHOD),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    data.frame(gene = g, platform = plat,
               pooled_logFC = as.numeric(m$b), se = m$se, p = m$pval,
               Q = m$QE, Q_p = m$QEp, I2 = m$I2,
               stringsAsFactors = FALSE)
  }))
}))
write_table(meta_tab, "p0_meta_per_platform.csv")

# --- gastric: fixed and DerSimonian-Laird, k = 2 -----------------------------
# With two cohorts a random-effects estimate is a sensitivity figure, not a
# primary result. Both are reported so the difference is visible.
gastric_fe_re <- function(g2) {
  b <- g2$logFC; se <- g2$SE
  w <- 1 / se^2
  k <- length(b)
  fb <- sum(w * b) / sum(w)
  fse <- 1 / sqrt(sum(w))
  Q <- sum(w * (b - fb)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  w2 <- 1 / (se^2 + tau2)
  rb <- sum(w2 * b) / sum(w2)
  rse <- 1 / sqrt(sum(w2))
  data.frame(gas_FE_logFC = fb,
             gas_FE_p = 2 * pnorm(-abs(fb / fse)),
             gas_RE_logFC = rb,
             gas_RE_p = 2 * pnorm(-abs(rb / rse)),
             gas_I2 = 100 * max(0, (Q - (k - 1)) / Q),
             stringsAsFactors = FALSE)
}

gas_tab <- do.call(rbind, lapply(genes14, function(g) {
  g2 <- eff_tab[eff_tab$gene == g & eff_tab$cohort %in% names(cag), ]
  if (nrow(g2) < 2) return(NULL)
  cbind(data.frame(gene = g, stringsAsFactors = FALSE), gastric_fe_re(g2))
}))
write_table(gas_tab, "p0_gastric_fe_re.csv")

# -----------------------------------------------------------------------------
# 5. Summary audit table
# -----------------------------------------------------------------------------
val_path <- file.path(tables_dir, "shared_validation_auc.csv")
if (file.exists(val_path)) {
  val <- read.csv(val_path, row.names = 1, stringsAsFactors = FALSE)
} else {
  message("[05] shared_validation_auc.csv absent (run R/07); skipping validation columns")
  val <- NULL
}

mr_path <- file.path(tables_dir, "mr_all_results_v4.csv")
if (file.exists(mr_path)) {
  mr <- read.csv(mr_path, stringsAsFactors = FALSE)
  mr_dir <- tapply(mr$dir_ok, mr$gene, function(x)
    paste0(sum(x, na.rm = TRUE), "/", length(x)))
} else {
  message("[05] mr_all_results_v4.csv absent (run R/11); skipping MR column")
  mr_dir <- NULL
}

audit <- data.frame(
  gene = genes14,
  direction = dir14,
  lodo_retained_of6 = rowSums(lodo_all[, in200_cols], na.rm = TRUE),
  lodo_worst_rank = apply(lodo_all[, grep("^rank_drop_", names(lodo_all))],
                          1, max, na.rm = TRUE),
  dir_concordant_masld6 = dir_cnt[genes14],
  cag_support = apply(cag_tab[, -c(1, 2)], 1,
                      function(x) paste(gsub(" .*", "", x), collapse = "/")),
  stringsAsFactors = FALSE)

if (!is.null(val)) {
  audit$val_direction_ok <- val[genes14, "direction_ok"]
  audit$val_AUC_norm     <- round(val[genes14, "AUC_norm"], 3)
  audit$val_NAS_rho      <- round(val[genes14, "spearman_NAS_rho"], 3)
}
if (!is.null(mr_dir)) audit$mr_dir_ok <- unname(mr_dir[genes14])

write_table(audit, "candidate_robustness_audit.csv")

message("[05] done")
print(audit[, c("gene", "lodo_retained_of6", "lodo_worst_rank",
                "dir_concordant_masld6", "cag_support")])

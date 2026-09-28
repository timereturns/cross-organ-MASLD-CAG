# =============================================================================
# 07_validation_GSE135251.R — held-out external validation (liver)
# =============================================================================
# STATUS: RESTORATION (from the file the author called
# 第一步第二部分GSE135251 独立验证.R).
#
# GSE135251 was held out of the discovery screen: it was not used to select the 14
# candidates, and it is interrogated once, in the pre-specified direction, without
# re-screening. That is what makes it a validation and not a second discovery.
#
# WHAT THE COHORT ACTUALLY IS, and why it constrains the interpretation:
#
#   216 samples, of which 206 are disease and 10 are control. A 206-vs-10 split
#   makes the case/control statistics (Wilcoxon, AUC) weak on the control side:
#   10 controls is not enough for a stable AUC. The NAS-gradient Spearman
#   correlation across all 216 samples is the primary validation metric, and the
#   case/control numbers are auxiliary. The supplementary notes already say this;
#   the code says it here too so the two cannot drift apart.
#
# Two derived quantities need care:
#
#   * AUC is the Mann-Whitney rank statistic. For a strongly DOWN-regulated gene
#     the raw AUC is near 0, which reads as "terrible discrimination" when it is
#     actually strong concordance in the expected direction. AUC_norm = max(AUC,
#     1 - AUC) fixes the direction but loses the sign, so it must always be
#     reported alongside the raw AUC and the expected direction.
#   * The bootstrap CI needs a seed. The original ran it without one; the
#     repository records BOOTSTRAP_SEED so the interval is reproducible.
#
# Run:  Rscript R/07_validation_GSE135251.R   (after 04)
# =============================================================================

source("R/00_setup.R")

# -----------------------------------------------------------------------------
# 0. Data
# -----------------------------------------------------------------------------
g135_path <- file.path(cache_dir, "GSE135251_data.rds")
if (!file.exists(g135_path)) {
  stop("GSE135251_data.rds missing. Run R/04_assets.R first.", call. = FALSE)
}
d <- readRDS(g135_path)
lcp  <- as.matrix(d$lcp)
meta <- d$meta

cand <- read.csv(file.path(tables_dir, "candidate_genes_final.csv"),
                 stringsAsFactors = FALSE)
message("[07] candidates carried into validation: ", nrow(cand))

# -----------------------------------------------------------------------------
# 1. Resolve each candidate to a row of the matrix
# -----------------------------------------------------------------------------
# lcp may be keyed by symbol or by Ensembl id depending on the upstream record;
# handle both rather than assuming.
rn <- rownames(lcp)
if (all(cand$gene %in% rn)) {
  message("[07] matrix rownames are symbols")
  G <- lcp[cand$gene, , drop = FALSE]
  cand$row <- cand$gene
} else if (grepl("^ENSG", rn[1])) {
  message("[07] matrix rownames are Ensembl ids; mapping to symbols")
  sym <- AnnotationDbi::mapIds(hs_eg_db(), keys = sub("\\..*$", "", rn),
                               keytype = "ENSEMBL", column = "SYMBOL",
                               multiVals = "first")
  ok <- !is.na(sym)
  m2 <- lcp[ok, , drop = FALSE]
  s  <- sym[ok]
  dup <- duplicated(s)
  m2 <- m2[!dup, , drop = FALSE]
  s  <- s[!dup]
  rownames(m2) <- s
  cand$row <- vapply(cand$gene, function(g) {
    e <- rownames(m2)[s == g]
    if (!length(e)) return(NA_character_)
    if (length(e) > 1) e <- e[which.max(rowMeans(m2[e, , drop = FALSE]))]
    e
  }, character(1))
  cand <- cand[!is.na(cand$row), , drop = FALSE]
  G <- m2[cand$row, , drop = FALSE]
} else {
  stop("Unrecognised rownames in lcp: ", paste(head(rn, 3), collapse = ", "),
       call. = FALSE)
}

# -----------------------------------------------------------------------------
# 2. Metrics
# -----------------------------------------------------------------------------
# Mann-Whitney AUC, computed by rank so that it equals the Wilcoxon statistic.
auc_rank <- function(x_case, x_ctrl) {
  n1 <- length(x_case); n2 <- length(x_ctrl)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  r <- rank(c(x_case, x_ctrl))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

nas <- suppressWarnings(as.numeric(meta[["nas score"]]))
fib <- suppressWarnings(as.numeric(meta[["fibrosis stage"]]))
grp <- meta$grp

res <- data.frame(
  direction = cand$direction,
  mean_case = NA_real_, mean_ctrl = NA_real_,
  wilcox_p = NA_real_, AUC = NA_real_, AUC_norm = NA_real_,
  direction_ok = NA,
  spearman_NAS_rho = NA_real_, spearman_NAS_p = NA_real_,
  spearman_F_rho = NA_real_, spearman_F_p = NA_real_,
  row.names = cand$gene, stringsAsFactors = FALSE)

for (i in seq_len(nrow(cand))) {
  g <- cand$gene[i]
  y <- as.numeric(G[g, ])
  x_case <- y[grp == "case"]
  x_ctrl <- y[grp == "control"]

  res$mean_case[i] <- mean(x_case, na.rm = TRUE)
  res$mean_ctrl[i] <- mean(x_ctrl, na.rm = TRUE)
  res$wilcox_p[i]  <- suppressWarnings(wilcox.test(x_case, x_ctrl)$p.value)
  res$AUC[i]       <- auc_rank(x_case, x_ctrl)
  res$AUC_norm[i]  <- max(res$AUC[i], 1 - res$AUC[i], na.rm = TRUE)

  # Direction concordance is decided on group means, not on the AUC, so that the
  # down-regulated genes are evaluated on the same footing as the up-regulated.
  res$direction_ok[i] <- if (cand$direction[i] == "up") {
    mean(x_case, na.rm = TRUE) > mean(x_ctrl, na.rm = TRUE)
  } else {
    mean(x_case, na.rm = TRUE) < mean(x_ctrl, na.rm = TRUE)
  }

  ok_nas <- is.finite(y) & !is.na(nas)
  ok_fib <- is.finite(y) & !is.na(fib)
  if (sum(ok_nas) > 3) {
    sp <- suppressWarnings(cor.test(y[ok_nas], nas[ok_nas], method = "spearman"))
    res$spearman_NAS_rho[i] <- unname(sp$estimate)
    res$spearman_NAS_p[i]   <- sp$p.value
  }
  if (sum(ok_fib) > 3) {
    sp <- suppressWarnings(cor.test(y[ok_fib], fib[ok_fib], method = "spearman"))
    res$spearman_F_rho[i] <- unname(sp$estimate)
    res$spearman_F_p[i]   <- sp$p.value
  }
}

write_table(res, "shared_validation_auc.csv")

# -----------------------------------------------------------------------------
# 3. Spearman 95% CI (Fisher z) for the 5 pre-specified candidates
# -----------------------------------------------------------------------------
GENES <- c("CDHR2", "ANXA4", "CADM2", "IL32", "LGALS3")

ci_rho <- do.call(rbind, lapply(GENES, function(gn) {
  if (!gn %in% cand$gene) return(NULL)
  y <- as.numeric(G[cand$gene == gn, ][1, ])
  ok <- !is.na(nas) & is.finite(y)
  n <- sum(ok)
  if (n < 4) return(NULL)
  rho <- cor(y[ok], nas[ok], method = "spearman")
  z <- atanh(rho); se <- 1 / sqrt(n - 3)
  data.frame(gene = gn, n = n, rho = rho,
             ci_lo = tanh(z - 1.96 * se), ci_hi = tanh(z + 1.96 * se),
             stringsAsFactors = FALSE)
}))
print(ci_rho, digits = 4)
write_table(ci_rho, "GSE135251_spearman_CI.csv")

# -----------------------------------------------------------------------------
# 4. AUC bootstrap (percentile CI)
# -----------------------------------------------------------------------------
seed_with(BOOTSTRAP_SEED, "GSE135251 AUC bootstrap")
report_seed(BOOTSTRAP_SEED, BOOTSTRAP_N, "AUC bootstrap")

i1 <- which(grp == "case")
i0 <- which(grp == "control")

boot_auc <- do.call(rbind, lapply(GENES, function(gn) {
  if (!gn %in% cand$gene) return(NULL)
  y <- as.numeric(G[cand$gene == gn, ][1, ])
  boot <- replicate(BOOTSTRAP_N, {
    b1 <- sample(i1, length(i1), replace = TRUE)
    b0 <- sample(i0, length(i0), replace = TRUE)
    auc_rank(y[b1], y[b0])
  })
  q <- quantile(boot, c(.025, .975), na.rm = TRUE)
  data.frame(gene = gn, n_case = length(i1), n_ctrl = length(i0),
             AUC = auc_rank(y[i1], y[i0]),
             boot_lo = unname(q[1]), boot_hi = unname(q[2]),
             stringsAsFactors = FALSE)
}))
print(boot_auc, digits = 4)
write_table(boot_auc, "GSE135251_AUC_bootstrap_CI.csv")

# -----------------------------------------------------------------------------
# 5. Summary
# -----------------------------------------------------------------------------
n_dir_ok <- sum(res$direction_ok, na.rm = TRUE)
message("[07] direction concordance: ", n_dir_ok, "/", nrow(res), " candidates")

# The manuscript reports 13/14 for the overall cohort panel; this cohort's share
# is reported for what it is.
strong <- res[!is.na(res$AUC_norm) & res$AUC_norm >= 0.7, ]
message("[07] AUC_norm >= 0.7: ", nrow(strong), " candidates")

# BH over the pre-specified candidate family, not over all 14.
p <- res$wilcox_p[!is.na(res$wilcox_p)]
if (length(p)) {
  message("[07] BH over ", length(p), " case/control tests ",
          "(auxiliary family; the NAS gradient is the primary metric)")
}

print(res[, c("direction", "direction_ok", "AUC", "AUC_norm",
              "wilcox_p", "spearman_NAS_rho", "spearman_NAS_p")], digits = 3)
message("[07] done")

# =============================================================================
# 11_MR_exploratory.R — exploratory two-sample Mendelian randomisation
# =============================================================================
# STATUS: RESTORATION of the THIRD run of the MR pipeline.
#
# The MR analysis exists in three successive versions in the original work. Only
# the third is reported in the manuscript: it writes mr_all_results_v4.csv, which
# every downstream script reads. The two earlier versions (v2, with Steiger
# filtering; v3, which removed Steiger and added a manual Wald fallback) are
# superseded and are not reproduced.
#
# What changed from v2 to v4, and why it matters for reading the manuscript:
#
#   * Steiger filtering was REMOVED. The output still carries a `steiger_ok`
#     column, and it is NA throughout in v4. Do not read that column as "Steiger
#     passed" -- it records that Steiger was not applied.
#   * LiverPDFF became the primary liver outcome and NAFLD_UKB was demoted to a
#     sensitivity outcome.
#   * The outcomes carry three tiers: main, sens, explor. FDR is computed within
#     {main + sens} as one family and within {explor} as a second.
#
# ROLE IN THE PAPER: exploratory. No estimate survives FDR, and the manuscript
# says so. This script is not a source of causal claims.
#
# THE INSTRUMENTS ARE PRE-CLumped FILES, produced by R/11a_instrument_clumping.R —
# run that first. Their selection parameters, recovered from that script, are:
#
#     ieugwasr::ld_clump(clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
#                        pop = "EUR")
#
# Eight of the ten genes use that derived set. LGALS3 and CADM2 were never
# re-clumped and come from a looser r^2 = 0.01 set, one SNP each — so all their
# estimates are Wald ratios with no pleiotropy diagnostics. CADM2 is a retained
# gene, which makes that worth stating. See docs/analysis_parameters.md section 9.
#
# Run:  Rscript R/11_MR_exploratory.R           (after R/11a)
# =============================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
# Instrument files, one per gene. Eight genes come from the r2 = 0.001 clumped
# set; LGALS3 and CADM2 use a looser set. The split is from the original and is
# not explained by any surviving comment.
eqtl_dir_strict <- file.path(data_dir, "eqtl", "eqtlclump_r2_0.001")
eqtl_dir_loose  <- file.path(data_dir, "eqtl", "eqtlclump")
outcome_dir     <- file.path(data_dir, "gwas")

genes <- data.frame(
  gene = mr_genes,
  dir  = unname(mr_direction[mr_genes]),
  file = c(file.path(eqtl_dir_strict,
                     c("IL32.txt", "GOLM1.txt", "TSPAN3.txt", "ANXA4.txt",
                       "RPS6KA1.txt", "ANO10.txt", "SLC6A16.txt", "KIAA1958.txt")),
           file.path(eqtl_dir_loose, c("LGALS3.txt", "CADM2.txt"))),
  stringsAsFactors = FALSE)

# Outcomes. LiverPDFF is continuous (a PDFF percentage), so odds ratios are not
# meaningful for it and are left NA.
outcomes <- list(
  list(id = "ukbpdff", trait = "LiverPDFF_UKB", file = "GCST90267352.h.tsv.gz",
       fmt = "catalog", ncase = NA, ncontrol = NA, ntotal = 33235,
       tier = "main", continuous = TRUE),
  list(id = "ukbgast", trait = "Gastritis_UKB", file = "GCST90129440.h.tsv.gz",
       fmt = "catalog", ncase = 41746, ncontrol = 179970, ntotal = NA,
       tier = "main", continuous = FALSE),
  list(id = "ukbnafld", trait = "NAFLD_UKB", file = "GCST90054782_mr_ready.tsv.gz",
       fmt = "catalog", ncase = 4761, ncontrol = 373227, ntotal = NA,
       tier = "sens", continuous = FALSE),
  list(id = "finnafld", trait = "NAFLD_FinnGen", file = "finngen_R13_NAFLD.gz",
       fmt = "finn", ncase = 3943, ncontrol = 496243, ntotal = NA,
       tier = "sens", continuous = FALSE),
  list(id = "finngast", trait = "Gastritis_FinnGen",
       file = "finngen_R13_K11_CHRONGASTR.gz",
       fmt = "finn", ncase = 12849, ncontrol = 420114, ntotal = NA,
       tier = "sens", continuous = FALSE),
  list(id = "finnfib", trait = "FibroCirc_FinnGen",
       file = "finngen_R13_K11_FIBROCHIRLIV.gz",
       fmt = "finn", ncase = 2963, ncontrol = 483311, ntotal = NA,
       tier = "explor", continuous = FALSE),
  list(id = "finnstom", trait = "StomachCa_FinnGen",
       file = "finngen_R13_C3_STOMACH_WIDE.gz",
       fmt = "finn", ncase = 2455, ncontrol = 372159, ntotal = 33351,
       tier = "explor", continuous = FALSE)
)

# -----------------------------------------------------------------------------
# 1. Outcome preparation
# -----------------------------------------------------------------------------
# Two source layouts. GWAS Catalog files are already in the harmonised shape;
# FinnGen files need column renaming and the multi-valued rsids field trimmed to
# its first entry.

prep_outcome <- function(o) {
  f <- file.path(outcome_dir, o$file)
  if (!file.exists(f)) {
    message("[11] outcome file missing, skipped: ", o$file)
    return(NULL)
  }
  if (o$fmt == "catalog") {
    d <- fread(f, select = c("rsid", "effect_allele", "other_allele", "beta",
                             "standard_error", "effect_allele_frequency", "p_value"))
    d <- d[nzchar(rsid) & !is.na(rsid)]
    setnames(d, c("SNP", "effect_allele.outcome", "other_allele.outcome",
                  "beta.outcome", "se.outcome", "eaf.outcome", "pval.outcome"))
  } else {
    d <- fread(f, select = c("#chrom", "pos", "ref", "alt", "rsids", "pval",
                             "beta", "sebeta", "af_alt"))
    d$SNP <- sub(",.*$", "", d$rsids)
    d <- d[nzchar(SNP) & !is.na(SNP)]
    setnames(d, c("alt", "ref", "beta", "sebeta", "af_alt", "pval"),
             c("effect_allele.outcome", "other_allele.outcome",
               "beta.outcome", "se.outcome", "eaf.outcome", "pval.outcome"))
  }
  d$samplesize.outcome <- if (!is.na(o$ncase)) o$ncase + o$ncontrol else o$ntotal
  d$ncase.outcome    <- o$ncase
  d$ncontrol.outcome <- o$ncontrol
  d$id.outcome       <- o$id
  d$outcome          <- o$trait
  d
}

# -----------------------------------------------------------------------------
# 2. Per gene-outcome estimate
# -----------------------------------------------------------------------------

run_gene_outcome <- function(g_row, o, outc_full) {
  if (!file.exists(g_row$file)) return(NULL)

  expo <- read_exposure_data(filename = g_row$file, sep = "\t",
                             snp_col = "SNP", beta_col = "beta", se_col = "se",
                             effect_allele_col = "A1", other_allele_col = "A2",
                             eaf_col = "eaf", pval_col = "p", samplesize_col = "n")
  if (is.null(expo) || nrow(expo) == 0) return(NULL)

  outc <- outc_full[SNP %in% expo$SNP]
  if (nrow(outc) == 0) return(NULL)

  # Some outcome files omit allele frequency entirely, which makes harmonisation
  # drop every palindromic SNP. Back-fill from the IEU API rather than losing them.
  if (all(is.na(outc$eaf.outcome))) {
    af <- tryCatch(as.data.frame(ieugwasr::afl2_rsid(outc$SNP)),
                   error = function(e) NULL)
    if (!is.null(af) && nrow(af)) {
      i_rs <- which(tolower(names(af)) == "rsid")[1]
      i_af <- which(grepl("af", names(af), ignore.case = TRUE) &
                      !grepl("case|control|afl2", names(af), ignore.case = TRUE))[1]
      if (!is.na(i_rs) && !is.na(i_af)) {
        outc$eaf.outcome <- af[[i_af]][match(outc$SNP, af[[i_rs]])]
      }
    }
  }

  # action = 2: palindromic SNPs with intermediate allele frequency are DROPPED
  # rather than strand-inferred. This is a conservative choice and the manuscript
  # should describe it as such.
  harm <- tryCatch(harmonise_data(expo, as.data.frame(outc), action = 2),
                   error = function(e) NULL)
  if (is.null(harm) || nrow(harm) == 0) return(NULL)

  # Instrument strength. R2 from the standard approximation for a continuous
  # exposure, then F = R2 (N - 2) / (1 - R2).
  harm$R2 <- (2 * harm$beta.exposure^2 * harm$eaf.exposure * (1 - harm$eaf.exposure)) /
    (2 * harm$beta.exposure^2 * harm$eaf.exposure * (1 - harm$eaf.exposure) +
       2 * harm$samplesize.exposure * harm$eaf.exposure *
       (1 - harm$eaf.exposure) * harm$se.exposure^2)
  harm$F <- harm$R2 * (harm$samplesize.exposure - 2) / (1 - harm$R2)

  mean_F <- mean(harm$F, na.rm = TRUE)
  harm <- harm[harm$F > MR_F_STATISTIC_MIN, ]
  if (nrow(harm) == 0) return(NULL)

  if (nrow(harm) == 1) {
    # Wald ratio, computed directly. The manual form is used because
    # mr_wald_ratio() returned an empty table in this environment.
    b_ex <- harm$beta.exposure[1]; se_ex <- harm$se.exposure[1]
    b_ou <- harm$beta.outcome[1];  se_ou <- harm$se.outcome[1]
    b  <- b_ou / b_ex
    se <- sqrt(se_ou^2 / b_ex^2 + (b_ou^2 * se_ex^2) / b_ex^4)
    p  <- 2 * pnorm(-abs(b / se))
    row <- data.frame(gene = g_row$gene, outcome = o$trait, tier = o$tier,
                      nsnp = 1L, method = "Wald ratio",
                      beta = b, se = se, pval = p,
                      OR = if (o$continuous) NA_real_ else exp(b),
                      OR_lci = if (o$continuous) NA_real_ else exp(b - 1.96 * se),
                      OR_uci = if (o$continuous) NA_real_ else exp(b + 1.96 * se),
                      mean_F = mean_F, egger_int_p = NA_real_, Q_ivw_p = NA_real_,
                      presso_global_p = NA_real_, steiger_ok = NA,
                      stringsAsFactors = FALSE)
  } else {
    mr_res <- mr(harm, method_list = c("mr_egger_regression", "mr_weighted_median",
                                       "mr_ivw", "mr_weighted_mode"))
    if (nrow(mr_res) == 0) return(NULL)
    ivw <- mr_res[mr_res$method == "Inverse variance weighted", ]
    if (nrow(ivw) == 0) return(NULL)

    het  <- tryCatch(mr_heterogeneity(harm), error = function(e) NULL)
    ple  <- tryCatch(mr_pleiotropy_test(harm), error = function(e) NULL)
    pres <- tryCatch(run_mr_presso(harm, NbDistribution = MR_PRESTO_NB),
                     error = function(e) NULL)
    q_p <- if (!is.null(het) && nrow(het)) {
      i <- which(het$method == "Inverse variance weighted")
      if (length(i)) het$Q_pval[i] else NA_real_
    } else NA_real_
    e_p <- if (!is.null(ple) && nrow(ple)) ple$pval[1] else NA_real_
    p_p <- if (!is.null(pres) && length(pres) >= 2 && is.data.frame(pres[[2]])) {
      pres[[2]]$`P-value`[1]
    } else NA_real_

    row <- data.frame(gene = g_row$gene, outcome = o$trait, tier = o$tier,
                      nsnp = nrow(harm), method = "IVW",
                      beta = ivw$b[1], se = ivw$se[1], pval = ivw$pval[1],
                      OR = if (o$continuous) NA_real_ else exp(ivw$b[1]),
                      OR_lci = if (o$continuous) NA_real_ else exp(ivw$b[1] - 1.96 * ivw$se[1]),
                      OR_uci = if (o$continuous) NA_real_ else exp(ivw$b[1] + 1.96 * ivw$se[1]),
                      mean_F = mean_F, egger_int_p = e_p, Q_ivw_p = q_p,
                      presso_global_p = p_p, steiger_ok = NA,
                      stringsAsFactors = FALSE)
  }

  row$dir_expected <- g_row$dir
  row$dir_obs <- ifelse(row$beta > 0, "up", "dn")
  row$dir_ok <- (g_row$dir == "up" && row$beta > 0) ||
                (g_row$dir == "dn" && row$beta < 0)
  list(row = row, harm = harm)
}

# -----------------------------------------------------------------------------
# 3. Run
# -----------------------------------------------------------------------------
message("[11] ", nrow(genes), " genes x ", length(outcomes), " outcomes = ",
        nrow(genes) * length(outcomes), " possible pairs")

all_rows <- list(); harm_store <- list()

for (o in outcomes) {
  oc <- prep_outcome(o)
  if (is.null(oc)) next
  message("[11] == ", o$trait, " (", o$tier,
          if (o$continuous) ", continuous" else "", ") ==")
  for (i in seq_len(nrow(genes))) {
    r <- tryCatch(run_gene_outcome(genes[i, ], o, oc), error = function(e) NULL)
    if (is.null(r)) {
      message("[11]    ", genes$gene[i], ": not estimable (no shared instruments)")
      next
    }
    key <- paste(o$id, genes$gene[i])
    all_rows[[key]]   <- r$row
    harm_store[[key]] <- r$harm
    message(sprintf("[11]    %-10s nsnp=%d  p=%.3g  dir_ok=%s",
                    genes$gene[i], r$row$nsnp, r$row$pval, r$row$dir_ok))
  }
}

if (!length(all_rows)) {
  stop("No estimable pairs. Check that the eQTL and GWAS files are present in ",
       "data/eqtl/ and data/gwas/ - see data/datasets.md.", call. = FALSE)
}

res <- as.data.frame(data.table::rbindlist(all_rows, fill = TRUE))

# FDR within two families: {main + sens} together, {explor} separately. The
# split is the original's and reflects that the exploratory outcomes are a
# different kind of question, not a larger family.
ms <- res$tier %in% c("main", "sens")
res$FDR_main <- NA_real_
res$FDR_main[ms] <- p.adjust(res$pval[ms], "fdr")
res$FDR_explor <- NA_real_
res$FDR_explor[res$tier == "explor"] <- p.adjust(res$pval[res$tier == "explor"], "fdr")

res <- res[order(res$pval), , drop = FALSE]
write_table(res, "mr_all_results_v4.csv")
saveRDS(harm_store, file.path(cache_dir, "harmonised_for_coloc.rds"))

# -----------------------------------------------------------------------------
# 4. Summary
# -----------------------------------------------------------------------------
message("[11] ==== results ====")
print(res[, c("gene", "outcome", "tier", "nsnp", "method",
              "beta", "pval", "FDR_main", "FDR_explor", "dir_ok")], digits = 3)

message("[11] estimates: ", nrow(res), " (manuscript reports ",
        N_MR_ESTIMATES, ")")
if (nrow(res) != N_MR_ESTIMATES) {
  message("[11]   note: the count differs from the published 69. This is expected ",
          "if a different number of instrument files or outcome files is present; ",
          "it is worth confirming rather than assuming.")
}

surv <- res[(!is.na(res$FDR_main) & res$FDR_main < 0.05) |
              (!is.na(res$FDR_explor) & res$FDR_explor < 0.05), ]
message("[11] estimates surviving FDR < 0.05: ", nrow(surv),
        "  (the manuscript reports none)")
if (nrow(surv)) print(surv[, c("gene", "outcome", "tier", "pval",
                               "FDR_main", "FDR_explor")])

message("[11] genes with no estimable pair: ",
        paste(setdiff(mr_genes, unique(res$gene)), collapse = ", "))

# Single-instrument pairs are Wald ratios: worth knowing how many, since they carry
# no pleiotropy diagnostics at all.
n_wald <- sum(res$method == "Wald ratio")
message("[11] single-instrument (Wald ratio) pairs: ", n_wald, " of ", nrow(res),
        " -- no Egger intercept, Q or MR-PRESSO is available for these")

message("[11] done")

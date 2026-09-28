# =============================================================================
# 11a_instrument_clumping.R — cis-eQTL instrument selection and QC
# =============================================================================
# STATUS: RESTORATION (from clump自查 - 8循环重跑).
#
# RUN ORDER: this precedes R/11_MR_exploratory.R, which reads the instrument files
# it produces. It is numbered 11a rather than 10a so that it sits beside the
# analysis it feeds.
#
# This is the step that produces the exposure instrument files the MR analysis
# reads. It exists because the parameters it records were previously only implicit
# in a directory name.
#
# WHAT IT FOUND, AND WHY IT IS WORTH HAVING IN THE REPOSITORY
#
#   1. The r2 = 0.01 instrument set was PURCHASED from a third party; the
#      r2 = 0.001 set was produced here by re-clumping it. So the stricter set is
#      this study's own work, and that is the set eight of the ten genes use.
#
#   2. TSPAN3's re-clump returned ZERO SNPs. The script has an explicit guard for
#      that and stops rather than writing an empty file. It was then re-run and, in
#      a second pass, one surviving LD pair still exceeded r2 = 0.001, so the
#      lower-F SNP of that pair was dropped by hand, guarded by
#      stopifnot(nrow(kept2) == nrow(d) - 1). Final state: 4 SNPs, max r2 <= 0.001.
#      The consequence is already recorded in the supplementary notes
#      ("ukbgast TSPAN3: 2 SNPs"); the repair itself is only visible here.
#
#   3. LGALS3 and CADM2 were never re-clumped. They are read from the purchased
#      r2 = 0.01 directory and each contributes ONE SNP, so every LGALS3 and CADM2
#      estimate in the MR table is a Wald ratio with no Egger intercept, no Q
#      statistic and no MR-PRESSO. CADM2 is one of the five retained genes, so this
#      belongs next to any MR-based support for it.
#
#   4. Two different F-statistic definitions are in play. This script records
#      instrument strength as F = (beta/se)^2, the simple univariate form. The MR
#      run computes F from R2 and is sample-size aware. They disagree. See
#      docs/analysis_parameters.md section 9.
#
# Requires internet access: ld_clump and ld_matrix call the IEU OpenGWAS API.
#
# Run:  Rscript R/11a_instrument_clumping.R     (before R/11)
# =============================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(ieugwasr)
  library(data.table)
})

# --- paths --------------------------------------------------------------------
# RAW holds the purchased r2 = 0.01 files; OUT receives the study's own r2 = 0.001
# set. Both are inputs, not repository contents: the purchased set is third-party
# and the derived set is produced by this script.
RAW <- file.path(data_dir, "eqtl", "eqtlclump_purchased_r2_0.01")
OUT <- file.path(data_dir, "eqtl", "eqtlclump_r2_0.001")

# --- the clumping parameters, now recorded as values rather than a directory name
CLUMP_KB <- 10000        # window, kb
CLUMP_R2 <- MR_CLUMPING_R2  # 0.001, from R/00_setup.R
CLUMP_P  <- 5e-8         # instrument selection p-value
CLUMP_POP <- "EUR"       # reference panel population

message("[15] clumping: kb = ", CLUMP_KB, ", r2 = ", CLUMP_R2,
        ", p = ", CLUMP_P, ", pop = ", CLUMP_POP)

# Eight genes are re-clumped. LGALS3 and CADM2 are deliberately excluded: they stay
# on the looser purchased set, one SNP each.
genes_reclumped <- c("IL32", "ANXA4", "RPS6KA1", "TSPAN3", "GOLM1", "ANO10",
                     "SLC6A16", "KIAA1958")
genes_loose <- c("LGALS3", "CADM2")

if (!dir.exists(RAW)) {
  stop("Purchased instrument directory not found:\n  ", RAW, "\n",
       "  It is third-party and is not redistributed with this repository.\n",
       "  See docs/repository_scope.md.", call. = FALSE)
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Re-clump to r2 < 0.001
# -----------------------------------------------------------------------------
clump_one <- function(g) {
  f <- file.path(RAW, paste0(g, ".txt"))
  if (!file.exists(f)) { message("[15] ", g, ": source file missing, skipped"); return(NULL) }

  dat <- fread(f)
  cl <- ieugwasr::ld_clump(data.frame(rsid = dat$SNP, pval = dat$p),
                           clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2,
                           clump_p = CLUMP_P, pop = CLUMP_POP)

  kept <- dat[SNP %in% cl$rsid]
  if (nrow(kept) == 0) {
    # The original stopped here rather than writing an empty file. That guard is
    # why the TSPAN3 problem surfaced instead of turning into a silent zero-SNP
    # instrument set.
    warning("[15] ", g, ": clumping returned 0 SNPs - NOT written",
            "\n  Review the source file and the API response before continuing.",
            call. = FALSE)
    return(NULL)
  }

  write.table(kept, file.path(OUT, paste0(g, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("[15] %-10s %d -> %d SNPs", g, nrow(dat), nrow(kept)))
  kept
}

message("[15] ==== re-clumping ", length(genes_reclumped), " genes ====")
for (g in genes_reclumped) clump_one(g)

# -----------------------------------------------------------------------------
# 2. Independence and strength check
# -----------------------------------------------------------------------------
# max r2 over the upper triangle of the LD matrix, and the F distribution. This is
# the check that found TSPAN3's remaining over-threshold pair.

check_one <- function(g) {
  f <- file.path(OUT, paste0(g, ".txt"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  fstat <- (d$beta / d$se)^2

  if (nrow(d) >= 2) {
    ld <- tryCatch(ieugwasr::ld_matrix(d$SNP, pop = CLUMP_POP,
                                        with_alleles = TRUE),
                   error = function(e) NULL)
    if (is.null(ld)) {
      ld_msg <- "LD matrix unavailable"
      max_r2 <- NA_real_
    } else {
      r2 <- ld[upper.tri(ld)]^2
      max_r2 <- max(r2, na.rm = TRUE)
      ld_msg <- sprintf("%d pairs, max r2 = %.5f %s", length(r2), max_r2,
                        if (isTRUE(max_r2 <= CLUMP_R2)) "OK" else "OVER THRESHOLD")
    }
  } else {
    ld_msg <- "single instrument, no LD check applicable"
    max_r2 <- NA_real_
  }

  message(sprintf("[15] %-10s %d SNPs | median F = %.1f  min F = %.1f | %s",
                  g, nrow(d), median(fstat), min(fstat), ld_msg))
  data.frame(gene = g, n_snp = nrow(d),
             median_F = round(median(fstat), 1), min_F = round(min(fstat), 1),
             max_r2 = if (is.na(max_r2)) NA_real_ else round(max_r2, 5),
             check = ld_msg, stringsAsFactors = FALSE)
}

message("[15] ==== independence and strength ====")
chk <- do.call(rbind, lapply(genes_reclumped, check_one))

# -----------------------------------------------------------------------------
# 3. Repair: drop the lower-F SNP of an over-threshold pair
# -----------------------------------------------------------------------------
# The original did this by hand for TSPAN3, guarded by an assertion that exactly
# one row was removed. Kept as a function so the operation is auditable rather
# than a one-off edit in a session.

repair_over_threshold <- function(g) {
  f <- file.path(OUT, paste0(g, ".txt"))
  d <- fread(f)
  if (nrow(d) < 2) { message("[15] ", g, ": fewer than 2 SNPs, nothing to repair"); return(invisible(NULL)) }

  ld <- ieugwasr::ld_matrix(d$SNP, pop = CLUMP_POP, with_alleles = TRUE)
  snps_bare <- sub("_[A-Z]+_[A-Z]+$", "", rownames(ld))   # strip the allele suffix
  r2mat <- ld^2
  idx <- which(r2mat > CLUMP_R2 & upper.tri(r2mat), arr.ind = TRUE)
  message("[15] ", g, ": over-threshold pairs = ", nrow(idx))
  if (nrow(idx) < 1) return(invisible(NULL))

  pair_bare <- snps_bare[idx[1, ]]
  fstat <- (d$beta / d$se)^2
  names(fstat) <- d$SNP
  drop <- pair_bare[which.min(fstat[pair_bare])]
  stopifnot(nzchar(drop))

  message("[15] ", g, ": dropping ", drop, " (F = ",
          round(min(fstat[pair_bare]), 1), ")")
  kept2 <- d[SNP != drop, ]
  stopifnot(nrow(kept2) == nrow(d) - 1)   # must lose exactly one row

  write.table(kept2, file.path(OUT, paste0(g, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  d2 <- fread(file.path(OUT, paste0(g, ".txt")))
  ld2 <- ieugwasr::ld_matrix(d2$SNP, pop = CLUMP_POP, with_alleles = TRUE)
  r2v <- ld2[upper.tri(ld2)]^2
  f2 <- (d2$beta / d2$se)^2
  message(sprintf("[15]   repaired: %d SNPs | max r2 = %.5f %s | median F = %.1f  min F = %.1f",
                  nrow(d2), max(r2v, na.rm = TRUE),
                  if (isTRUE(max(r2v, na.rm = TRUE) <= CLUMP_R2)) "OK" else "STILL OVER",
                  median(f2), min(f2)))
  invisible(NULL)
}

over <- chk$gene[!is.na(chk$max_r2) & chk$max_r2 > CLUMP_R2]
if (length(over)) {
  message("[15] ==== repairing: ", paste(over, collapse = ", "), " ====")
  for (g in over) repair_over_threshold(g)
} else {
  message("[15] no over-threshold pairs; nothing to repair")
}

# -----------------------------------------------------------------------------
# 4. The two genes that stay on the looser set
# -----------------------------------------------------------------------------
message("[15] ==== genes left on the purchased r2 = 0.01 set ====")
loose <- do.call(rbind, lapply(genes_loose, function(g) {
  f <- file.path(RAW, paste0(g, ".txt"))
  if (!file.exists(f)) { message("[15] ", g, ": missing"); return(NULL) }
  d <- fread(f)
  n <- nrow(d)
  message(sprintf("[15] %-10s %d SNP | F = %.1f | p = %.2e", g, n,
                  (d$beta[1] / d$se[1])^2, d$p[1]))
  if (n > 1) {
    message("[15]   note: ", g, " has ", n,
            " instruments on the loose set, so it is not a Wald ratio")
  } else {
    message("[15]   single instrument -> every ", g,
            " estimate is a Wald ratio, with no Egger intercept, Q or MR-PRESSO")
  }
  data.frame(gene = g, set = "purchased r2=0.01", n_snp = n,
             stringsAsFactors = FALSE)
}))

# -----------------------------------------------------------------------------
# 5. Summary
# -----------------------------------------------------------------------------
instr <- rbind(
  data.frame(gene = chk$gene, set = "derived r2=0.001", n_snp = chk$n_snp,
             median_F = chk$median_F, min_F = chk$min_F, max_r2 = chk$max_r2,
             stringsAsFactors = FALSE),
  loose)
instr <- instr[order(instr$gene), ]
write_table(instr, "p3_instrument_summary.csv")
print(instr)

message("[15] F here is the univariate (beta/se)^2 form. The MR run in R/11 uses ",
        "the sample-size-aware R2 form; the two disagree. See ",
        "docs/analysis_parameters.md section 9.")

# Cross-check against what the MR script will see.
mr_genes_present <- intersect(mr_genes, instr$gene)
message("[15] genes with instruments available to R/11: ",
        length(mr_genes_present), " of ", length(mr_genes),
        "  (absent: ", paste(setdiff(mr_genes, instr$gene), collapse = ", "), ")")

message("[15] done")

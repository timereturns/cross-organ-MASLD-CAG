# =============================================================================
# 03_cross_organ_candidates.R — gastric concordance and candidate finalisation
# =============================================================================
# STATUS: RESTORATION.
#
# This is the step that turns the liver discovery list into the 14 shared
# candidates. The logic is the original, extracted from the gastric download file,
# where it sat after the gastric limma runs.
#
# The filter chain, in the original order:
#   1. gastric DEGs per cohort: P.Value < 0.05 and |logFC| > 0.5
#   2. concordant genes: in the liver RRA list AND a gastric DEG, direction-matched
#   3. require significance in BOTH gastric cohorts, same direction
#   4. require >= 4 of the 6 liver cohorts to agree in direction
#   5. write candidate_genes_final.csv
#
# Steps 3 and 4 are the strongest internal guards in the whole screen: a gene has
# to survive an independent organ twice over, in the same direction both times.
#
# Run:  Rscript R/03_cross_organ_candidates.R   (after 01 and 02)
# =============================================================================

source("R/00_setup.R")

cache_cag <- file.path(cache_dir, "cag_limma_results.rds")
if (!file.exists(cache_cag)) {
  stop("Gastric limma results not found. Run R/01_download_and_DE.R first.",
       call. = FALSE)
}
cag <- readRDS(cache_cag)   # GSE116312, GSE27411
message("[03] gastric cohorts: ", paste(names(cag), collapse = ", "))

up_path <- file.path(tables_dir, "MASLD_RRA_up_final.csv")
dn_path <- file.path(tables_dir, "MASLD_RRA_dn_final.csv")
if (!all(file.exists(up_path, dn_path))) {
  stop("RRA output not found. Run R/02_RRA_discovery.R first.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 1. Gastric DEGs
# -----------------------------------------------------------------------------
deg <- function(d, p_cut = GASTRIC_DEG_P, fc_cut = GASTRIC_DEG_LOGF) {
  list(up = d$symbol[d$P.Value < p_cut & d$logFC >  fc_cut],
       dn = d$symbol[d$P.Value < p_cut & d$logFC < -fc_cut])
}
d_gas <- lapply(cag, deg)
for (nm in names(d_gas)) {
  message("[03] ", nm, ": ", nrow(cag[[nm]]), " genes tested, ",
          length(d_gas[[nm]]$up), " up / ", length(d_gas[[nm]]$dn), " dn")
}

# -----------------------------------------------------------------------------
# 2. Concordance with the liver discovery list
# -----------------------------------------------------------------------------
masld_up <- read.csv(up_path, stringsAsFactors = FALSE)$Name
masld_dn <- read.csv(dn_path, stringsAsFactors = FALSE)$Name

# a gene is concordant if it is a gastric DEG in the same direction in the
# respective cohorts
gas_up <- unique(unlist(lapply(d_gas, `[[`, "up")))
gas_dn <- unique(unlist(lapply(d_gas, `[[`, "dn")))

conc_up <- masld_up[masld_up %in% gas_up]
conc_dn <- masld_dn[masld_dn %in% gas_dn]

# discordant sets: recorded because a large discordant count would suggest the two
# organs are not sharing a program at all
disc_1 <- masld_up[masld_up %in% gas_dn]     # MASLD up, CAG down
disc_2 <- masld_dn[masld_dn %in% gas_up]     # MASLD down, CAG up

message("[03] concordant: ", length(conc_up), " up / ", length(conc_dn), " dn")
message("[03] discordant: ", length(disc_1), " up/dn, ", length(disc_2), " dn/up")

write_table(data.frame(gene = c(conc_up, conc_dn),
                       direction = c(rep("up_up", length(conc_up)),
                                     rep("dn_dn", length(conc_dn)))),
            "shared_genes_concordant.csv")

# -----------------------------------------------------------------------------
# 3. Require significance in BOTH gastric cohorts, same direction
# -----------------------------------------------------------------------------
# With only two gastric cohorts this is the strongest available gastric guard.
both_up <- intersect(d_gas[["GSE116312"]]$up, d_gas[["GSE27411"]]$up)
both_dn <- intersect(d_gas[["GSE116312"]]$dn, d_gas[["GSE27411"]]$dn)

final_up <- conc_up[conc_up %in% both_up]
final_dn <- conc_dn[conc_dn %in% both_dn]

message("[03] after both-gastric-cohort filter: ", length(final_up), " up / ",
        length(final_dn), " dn")

# -----------------------------------------------------------------------------
# 4. Require >= 4 of 6 liver cohorts to agree in direction
# -----------------------------------------------------------------------------
chips <- readRDS(file.path(cache_dir, "chips_limma_results.rds"))
rnas  <- readRDS(file.path(cache_dir, "rnaseq_limma_results.rds"))
six <- c(chips, rnas)

dir_tab <- function(genes) {
  sapply(genes, function(g) {
    tv <- vapply(six, function(d) {
      i <- which(d$symbol == g)
      if (length(i)) d$t[i[1]] else NA_real_
    }, numeric(1))
    tv <- tv[!is.na(tv)]
    c(up = sum(tv > 0), dn = sum(tv < 0), n = length(tv))
  })
}

du <- dir_tab(final_up)
dd <- dir_tab(final_dn)

keep_up <- final_up[if (length(final_up)) du["up", ] >= CANDIDATE_MIN_AGREEING_COHORTS else logical(0)]
keep_dn <- final_dn[if (length(final_dn)) dd["dn", ] >= CANDIDATE_MIN_AGREEING_COHORTS else logical(0)]

message("[03] FINAL candidates: ", length(keep_up), " up / ", length(keep_dn), " dn")

# -----------------------------------------------------------------------------
# 5. Output
# -----------------------------------------------------------------------------
cand <- data.frame(
  gene = c(keep_up, keep_dn),
  direction = c(rep("up", length(keep_up)), rep("dn", length(keep_dn))),
  stringsAsFactors = FALSE)

write_table(cand, "candidate_genes_final.csv")

# --- assertions --------------------------------------------------------------
# The published screen returned 14 candidates. This is not a rule of nature, so
# treat a different number as a signal to investigate rather than as a pass: it
# means either the data changed or a threshold moved.
message("[03] candidate count: ", nrow(cand),
        "  (the published screen returned 14)")

# The five retained genes and the eight downgraded ones must all be in the list;
# if they are not, the discovery chain no longer reproduces and nothing downstream
# is comparable.
known <- c(candidates_main, candidate_observation, candidates_downgraded)
missing <- setdiff(known, cand$gene)
if (length(missing)) {
  warning("[03] candidates expected from the published screen are absent: ",
          paste(missing, collapse = ", "),
          "\n  Investigate before using the downstream results.", call. = FALSE)
} else {
  message("[03] all 14 published candidates recovered")
}

print(cand)
message("[03] done")

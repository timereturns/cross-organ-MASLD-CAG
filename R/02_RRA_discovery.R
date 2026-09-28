# =============================================================================
# 02_RRA_discovery.R — robust rank aggregation discovery screen (liver)
# =============================================================================
# STATUS: RESTORATION.
#
# The analysis logic below is the original discovery screen, unchanged. What was
# removed is the RStudio-session scaffolding around it. See docs/build_notes.md.
#
# The screen, in order:
#   1. rank every liver cohort's genes by moderated t, separately for up and down
#   2. aggregate the six ranked lists with RobustRankAggreg::aggregateRanks
#      using the runtime size of the union of symbols as the background
#   3. keep genes with Score < 0.01 that are also directionally consistent in at
#      least 3 of the 6 cohorts
#   4. truncate to the top 200 in each direction
#
# NOTE ON THE BACKGROUND SIZE. Step 2 uses N <- length(all_sym), computed at run
# time, NOT a hardcoded 20,000. The manuscript quotes 20,000, which is the fixed
# background used by the LATER audit in 05_audit_LODO_meta.R, where three
# definitions were compared and 20,000 was chosen by hand. Both values are
# recorded; this script reports the one it actually used.
#
# Run:  Rscript R/02_RRA_discovery.R       (after 01)
# =============================================================================

source("R/00_setup.R")

if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) {
  stop("RobustRankAggreg is required; see README section 3", call. = FALSE)
}
library(RobustRankAggreg)

# -----------------------------------------------------------------------------
# 1. Load the liver cohort limma tables
# -----------------------------------------------------------------------------

chips_path <- file.path(cache_dir, "chips_limma_results.rds")
rnas_path  <- file.path(cache_dir, "rnaseq_limma_results.rds")

if (!all(file.exists(chips_path, rnas_path))) {
  stop("Limma results not found in cache/. Run R/01_download_and_DE.R first.",
       call. = FALSE)
}

chips <- readRDS(chips_path)   # Ahrens, Arendt, Lefebvre
rnas  <- readRDS(rnas_path)    # Suppli, Hoang, Pantano
six <- c(chips, rnas)
message("[02] cohorts: ", paste(names(six), collapse = ", "))

# -----------------------------------------------------------------------------
# 2. Ranked lists and background size
# -----------------------------------------------------------------------------
# One ranked symbol list per cohort per direction. unique() collapses symbols that
# survived symbol mapping more than once; the original screen did the same.

all_sym <- unique(unlist(lapply(six, function(d) unique(d$symbol))))
N <- length(all_sym)

message("[02] background size N = ", N,
        "  (runtime union of symbols across the six cohorts)")
message("[02] this is NOT the 20,000 used for the audit in R/05")

rank_up <- lapply(six, function(d) unique(d$symbol[order(-d$t)]))
rank_dn <- lapply(six, function(d) unique(d$symbol[order(d$t)]))

rra_up <- aggregateRanks(glist = rank_up, N = N)
rra_dn <- aggregateRanks(glist = rank_dn, N = N)

# -----------------------------------------------------------------------------
# 3. Threshold, direction filter, truncation
# -----------------------------------------------------------------------------

hub_up <- rra_up[rra_up$Score < RRA_SCORE_THRESHOLD, ]
hub_dn <- rra_dn[rra_dn$Score < RRA_SCORE_THRESHOLD, ]
message("[02] Score < ", RRA_SCORE_THRESHOLD, ": up ", nrow(hub_up),
        " / dn ", nrow(hub_dn))

# Directional consistency: a gene must be on the same side of zero in at least
# RRA_MIN_AGREEING_COHORTS of the six, by the sign of its moderated t.
dir_up <- lapply(six, function(d) d$symbol[d$t > 0])
dir_dn <- lapply(six, function(d) d$symbol[d$t < 0])

in_at_least_3_up <- names(which(table(unlist(dir_up)) >= RRA_MIN_AGREEING_COHORTS))
in_at_least_3_dn <- names(which(table(unlist(dir_dn)) >= RRA_MIN_AGREEING_COHORTS))

hub_up_final <- hub_up[hub_up$Name %in% in_at_least_3_up, , drop = FALSE]
hub_dn_final <- hub_dn[hub_dn$Name %in% in_at_least_3_dn, , drop = FALSE]

hub_up_top <- head(hub_up_final[order(hub_up_final$Score), , drop = FALSE],
                   RRA_TOP_N_PER_DIRECTION)
hub_dn_top <- head(hub_dn_final[order(hub_dn_final$Score), , drop = FALSE],
                   RRA_TOP_N_PER_DIRECTION)

message("[02] final lists: up ", nrow(hub_up_top), " / dn ", nrow(hub_dn_top))

# -----------------------------------------------------------------------------
# 4. Sensitivity check: drop Arendt and re-aggregate
# -----------------------------------------------------------------------------
# Recorded because the original analysis ran it. If the top-50 overlap is low, the
# discovery screen is being driven by a single cohort.

six_noAr <- six[setdiff(names(six), "Arendt")]
rra_noAr <- aggregateRanks(
  glist = lapply(six_noAr, function(d) unique(d$symbol[order(-d$t)])), N = N)

top50_full <- head(hub_up_top$Name, 50)
top50_noAr <- head(rra_noAr$Name, 50)
overlap_50  <- length(intersect(top50_full, top50_noAr))

message("[02] sensitivity (Arendt removed): top-50 overlap ", overlap_50, "/50")

write_table(data.frame(
  check = "top50 overlap, full vs Arendt-removed",
  overlap = overlap_50, of = 50L,
  background_N = N,
  score_threshold = RRA_SCORE_THRESHOLD,
  min_agreeing_cohorts = RRA_MIN_AGREEING_COHORTS),
  "rra_sensitivity.csv")

# -----------------------------------------------------------------------------
# 5. Output
# -----------------------------------------------------------------------------
# The ranked lists are written with a Score column; downstream scripts add the
# per-cohort rank columns (see R/04_assets.R).

write_table(hub_up_top, "MASLD_RRA_up_final.csv")
write_table(hub_dn_top, "MASLD_RRA_dn_final.csv")

message("[02] top 10 up:  ", paste(head(hub_up_top$Name, 10), collapse = ", "))
message("[02] top 10 dn:  ", paste(head(hub_dn_top$Name, 10), collapse = ", "))
message("[02] done")

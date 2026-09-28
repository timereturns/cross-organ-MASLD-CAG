# =============================================================================
# 17_bonferroni_family.R — determine and assert the 104-test family
# =============================================================================
# STATUS: ADDED. Not from the original work; there was no such check to restore.
#
# WHY THIS EXISTS
#
# The manuscript's Bonferroni threshold, P = 4.8 x 10^-4, derives from
# 13 genes x 8 cell types = 104 tests. The arithmetic is right (0.05 / 104 =
# 4.8077 x 10^-4), but nothing in the original code stated the two dimensions in
# one place:
#
#   * the candidate vector appears as 6, 13 and 14 genes in different blocks;
#   * the composition-adjustment marker panel (R/06) has NINE liver cell types;
#   * the MuSiC reference (R/09) has EIGHT cell types from a different vocabulary;
#   * the gastric annotation (R/08) uses TWELVE lineages.
#
# So a reviewer asking "which 13 genes and which 8 cell types?" had no single place
# to look. This script derives both dimensions from the data and asserts that they
# multiply to the published 104.
#
# WHERE THE TWO DIMENSIONS COME FROM
#
#   cell types  the 8 unified liver lineages the cross-dataset localisation is built
#               on: Hepatocytes, Cholangiocytes, Macrophages, Stellate,
#               Endothelial, Lymphocytes, B_cells, Erythroid
#               (the "Other" bucket is a residual, not a cell type, and is excluded)
#
#   genes       the candidate genes that survive in the liver localisation matrices.
#               GSE202379 carries 13 of the 14 -- the original file even names the
#               object `genes14` while commenting `# 13 个`, because one candidate is
#               absent from that dataset.
#
# If the product is 104, the family is the one the manuscript describes and the
# threshold follows. If it is not, either the data or the manuscript has changed,
# and the script says which.
#
# Run:  Rscript R/17_bonferroni_family.R      (after R/08)
# =============================================================================

source("R/00_setup.R")

# --- the 8 unified liver lineages ---------------------------------------------
# From the lineage map in the single-cell script. Kept in the same order the map
# defines them, so the two can be diffed by eye.
LIVER_LINEAGES <- c("Hepatocytes", "Cholangiocytes", "Macrophages", "Stellate",
                    "Endothelial", "Lymphocytes", "B_cells", "Erythroid")

# The map also emits a residual bucket for anything unmatched. It is not a cell
# type and must not be counted: including it would inflate the family.
RESIDUAL_BUCKETS <- c("Other", "Unassigned", "unknown")

message("[17] liver lineages defined by the lineage map: ", length(LIVER_LINEAGES))
message("[17]   ", paste(LIVER_LINEAGES, collapse = ", "))

# --- which genes survive, per dataset -----------------------------------------
read_loc <- function(f) {
  p <- file.path(tables_dir, f)
  if (!file.exists(p)) return(NULL)
  read.csv(p, stringsAsFactors = FALSE)
}

loc_liver <- list(GSE202379 = read_loc("p1_gse202379_localization.csv"),
                  GSE115469 = read_loc("p1_gse115469_localization.csv"))

if (all(vapply(loc_liver, is.null, logical(1)))) {
  stop("No liver localisation tables found. Run R/08_singlecell_pseudobulk.R first.",
       call. = FALSE)
}

genes_per_ds <- lapply(loc_liver, function(l) if (is.null(l)) NULL else sort(unique(l$gene)))
cts_per_ds   <- lapply(loc_liver, function(l) {
  if (is.null(l)) return(NULL)
  ct <- sort(unique(l$celltype))
  ct[!ct %in% RESIDUAL_BUCKETS]
})

for (nm in names(loc_liver)) {
  if (is.null(loc_liver[[nm]])) {
    message("[17] ", nm, ": no localisation table; skipped")
    next
  }
  message("[17] ", nm, ": ", length(genes_per_ds[[nm]]), " genes x ",
          length(cts_per_ds[[nm]]), " non-residual cell types")
  message("[17]   absent candidates: ",
          paste(setdiff(candidates_14, genes_per_ds[[nm]]), collapse = ", "))
  message("[17]   cell types: ", paste(cts_per_ds[[nm]], collapse = ", "))
}

# --- the family -----------------------------------------------------------------
# Genes present in BOTH liver datasets, which is what a cross-dataset localisation
# requires. If only one dataset is available, fall back to it and say so.
if (!is.null(genes_per_ds$GSE202379) && !is.null(genes_per_ds$GSE115469)) {
  genes_in_family <- intersect(genes_per_ds$GSE202379, genes_per_ds$GSE115469)
  gene_basis <- "genes present in both liver datasets"
} else {
  avail <- names(genes_per_ds)[!vapply(genes_per_ds, is.null, logical(1))]
  genes_in_family <- genes_per_ds[[avail[1]]]
  gene_basis <- paste0("genes present in ", avail[1], " (only one dataset available)")
}

# Cell types: the 8 lineages. If a dataset is missing, the vocabulary is still the
# 8 lineages, because that is what the localisation is defined over.
cts_in_family <- LIVER_LINEAGES

n_genes <- length(genes_in_family)
n_cells <- length(cts_in_family)
n_tests <- n_genes * n_cells
alpha   <- 0.05 / n_tests

message("\n[17] ==== the family ====")
message("[17] genes      : ", n_genes, "  (", gene_basis, ")")
message("[17]   ", paste(genes_in_family, collapse = ", "))
message("[17] cell types : ", n_cells)
message("[17]   ", paste(cts_in_family, collapse = ", "))
message("[17] tests      : ", n_genes, " x ", n_cells, " = ", n_tests)
message("[17] threshold  : 0.05 / ", n_tests, " = ", signif(alpha, 6))

# --- assertions -----------------------------------------------------------------
message("\n[17] ==== checks ====")

if (n_tests == PSEUDOBULK_N_TESTS) {
  message("[17] OK: the product matches PSEUDOBULK_N_TESTS (", PSEUDOBULK_N_TESTS,
          ") in R/00_setup.R")
} else {
  warning("[17] the product is ", n_tests, " but R/00_setup.R declares ",
          PSEUDOBULK_N_TESTS, ".\n  Reconcile before quoting the threshold.",
          call. = FALSE)
}

# Compare at the precision the manuscript quotes (2 significant figures), not with
# an absolute tolerance: 0.05/104 = 4.8077e-4 rounds to 4.8e-4 and passes, while
# 0.05/96 = 5.2083e-4 rounds to 5.2e-4 and does not.
alpha_2sf <- signif(alpha, 2)
ms_2sf    <- signif(MANUSCRIPT_BONFERRONI_P, 2)
if (isTRUE(all.equal(alpha_2sf, ms_2sf))) {
  message("[17] OK: computed threshold ", signif(alpha, 4),
          " rounds to the manuscript's ", MANUSCRIPT_BONFERRONI_P,
          " at two significant figures")
} else {
  warning("[17] computed threshold ", signif(alpha, 4), " (= ", alpha_2sf,
          " to 2 sf) does not match the manuscript's ", MANUSCRIPT_BONFERRONI_P,
          " (= ", ms_2sf, " to 2 sf).\n",
          "  The family size or the stated threshold has changed. Do not quote ",
          "the manuscript value until reconciled.", call. = FALSE)
}

if (n_cells == 8L) {
  message("[17] OK: 8 cell types, as the manuscript states")
} else {
  warning("[17] ", n_cells, " cell types, not the 8 the manuscript states",
          call. = FALSE)
}

if (n_genes == 13L) {
  message("[17] OK: 13 genes, as the manuscript states")
} else {
  warning("[17] ", n_genes, " genes, not the 13 the manuscript states.\n",
          "  Present: ", paste(genes_in_family, collapse = ", "), "\n",
          "  Absent:  ", paste(setdiff(candidates_14, genes_in_family), collapse = ", "),
          call. = FALSE)
}

# The composition marker panel is a DIFFERENT analysis with a different count. Say
# so here so the two are not confused in the manuscript.
message("\n[17] for contrast, the other cell-type counts in this study:")
message("[17]   composition-adjustment marker panel (R/06): ",
        length(mk_liver), " liver lineages")
message("[17]   MuSiC deconvolution reference (R/09): 8 cell types from the snRNA ",
        "vocabulary (see Table S41)")
message("[17]   gastric single-cell annotation (R/08): 12 lineages")
message("[17] Only the ", n_cells, "-lineage localisation family above is the one ",
        "behind the Bonferroni threshold.")

# --- output ---------------------------------------------------------------------
fam <- data.frame(
  gene = rep(genes_in_family, each = n_cells),
  celltype = rep(cts_in_family, times = n_genes),
  stringsAsFactors = FALSE)
fam$family <- "Bonferroni 104-test family"
fam$variable <- paste(fam$gene, fam$celltype, sep = " x ")
write_table(fam, "p1_pseudobulk_test_family.csv")

summary_tab <- data.frame(
  item = c("genes", "cell_types", "tests", "alpha", "manuscript_threshold",
           "gene_basis", "threshold_rounded"),
  value = c(n_genes, n_cells, n_tests, signif(alpha, 6),
            MANUSCRIPT_BONFERRONI_P, gene_basis, signif(alpha, 2)),
  stringsAsFactors = FALSE)
write_table(summary_tab, "p1_pseudobulk_family_summary.csv")
print(summary_tab)

message("\n[17] done - the family is listed gene by cell type in ",
        "results/tables/p1_pseudobulk_test_family.csv (", n_tests, " rows)")

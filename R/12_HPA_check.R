# =============================================================================
# 12_HPA_check.R — Human Protein Atlas protein-level annotation check
# =============================================================================
# STATUS: RESTORATION, with one addition.
#
# This produces the protein-level verification table for the five retained
# candidates (Table S27). The values are CURATOR-ASSIGNED: protein levels,
# cell-type strings, antibody identifiers and reliability grades were read off the
# Human Protein Atlas by hand and typed into a table. Nothing here queries HPA
# programmatically, and there is no per-gene retrieval record.
#
# That is a legitimate way to build an annotation table -- HPA's own reliability
# grades are the product of human curation too. What matters is that a reader can
# tell it apart from a computed result.
#
# THREE THINGS THE ORIGINAL DID NOT SETTLE, all marked below:
#
#   1. VERSION. The manuscript states HPA version 25.0. That number appears
#      nowhere in the code. It is carried here as a variable so the two can be
#      reconciled, but it has not been verified against a retrieval record.
#
#   2. ACCESS DATE. The manuscript says 5 September 2026. The supplementary NOTES
#      say 2026-09-13. Both are recorded below and the script reports the
#      disagreement rather than picking one.
#
#   3. ONE CELL CONTRADICTS ITSELF. `stomach_cells` is the constant string
#      "Glandular cells" for all five genes, including CADM2 whose own
#      `stomach_level` is "Not detected". The script now flags any gene whose
#      cell-type string asserts a localisation while its level says not detected,
#      so the contradiction surfaces at run time instead of sitting unnoticed in a
#      submitted table. The value itself is NOT changed: silently "fixing" curated
#      data would alter a published table.
#
# Run:  Rscript R/12_HPA_check.R
# =============================================================================

source("R/00_setup.R")

# --- the two contested metadata fields ---------------------------------------
HPA_VERSION      <- "25.0"          # manuscript; unverified against a record
HPA_ACCESS_MANUSCRIPT <- "2026-09-05"
HPA_ACCESS_NOTES      <- "2026-09-13"

if (HPA_ACCESS_MANUSCRIPT != HPA_ACCESS_NOTES) {
  warning("[12] HPA access date disagrees between the manuscript and the ",
          "supplementary notes:\n",
          "  manuscript: ", HPA_ACCESS_MANUSCRIPT, "\n",
          "  notes:      ", HPA_ACCESS_NOTES, "\n",
          "  Reconcile before submission.", call. = FALSE)
}

# --- the curated table --------------------------------------------------------
# Verbatim from the original script. Field meanings:
#   liver_level / stomach_level   HPA's ordinal level for that tissue
#   liver_cells / stomach_cells   the cell types HPA attributes the signal to
#   specificity                   HPA tissue-specificity category
#   antibody                      the antibody identifiers behind the call
#   reliability                   HPA's own reliability grade for those antibodies
#   manuscript_verdict            the author's interpretation, NOT an HPA field
hpa <- data.frame(
  gene = c("IL32", "ANXA4", "CDHR2", "CADM2", "LGALS3"),
  liver_level = c("Low", "Medium", "Medium (hepatocytes)", "Not detected",
                  "Low (cholangiocytes)"),
  liver_cells = c("Cholangiocytes Low; Hepatocytes Low",
                  "Cholangiocytes Medium-High; Hepatocytes Medium",
                  "Hepatocytes Medium; Cholangiocytes Not detected",
                  "Not detected (bulk IHC)",
                  "Cholangiocytes Low; myeloid cells by scRNA"),
  stomach_level = c("High/Medium", "Medium (High with 2 antibodies)", "Medium",
                    "Not detected (raw Low)", "High"),
  stomach_cells = rep("Glandular cells", 5),
  specificity = c("Immune-related", "GI-tract enriched", "Low tissue specificity",
                  "Brain enriched", "Low tissue specificity"),
  antibody = c("HPA029397; CAB030029",
               "HPA007393; CAB005076; CAB017560",
               "HPA012569; HPA017053",
               "HPA010024",
               "HPA003162; CAB005191"),
  reliability = c("Uncertain (IHC); WB Enhanced",
                  "Enhanced (orthogonal, 3 antibodies)",
                  "Supported / Enhanced",
                  "Enhanced (IHC, but ND in liver)",
                  "Enhanced (HPA003162); Supported (CAB005191)"),
  manuscript_verdict = c("limited: unvalidated IHC, pending annotation",
                         "strongest support: dual epithelial protein",
                         "good support: hepatocyte + gastric epithelial",
                         "negative-as-predicted: low-abundance stellate gene",
                         "gastric epithelial high; liver myeloid by scRNA"),
  stringsAsFactors = FALSE)

stopifnot(nrow(hpa) == 5L,
          all(hpa$gene %in% candidates_main))

# --- provenance columns -------------------------------------------------------
hpa$hpa_version <- HPA_VERSION
hpa$access_date <- HPA_ACCESS_MANUSCRIPT
hpa$provenance  <- "curator-assigned; no programmatic HPA query"

# --- consistency flag ---------------------------------------------------------
# A gene whose level says "Not detected" cannot also have a cell type attributed
# to it. This is not a formatting nitpick: it is the kind of pair a reviewer reads
# as carelessness.
not_detected <- grepl("^Not detected", hpa$stomach_level, ignore.case = TRUE)
has_cells    <- nzchar(hpa$stomach_cells) &
  !grepl("^(Not detected|not detected|ND)", hpa$stomach_cells)
hpa$flag_stomach_contradiction <- not_detected & has_cells

# --- antibody inventory -------------------------------------------------------
ab <- do.call(rbind, lapply(seq_len(nrow(hpa)), function(i) {
  ids <- trimws(strsplit(hpa$antibody[i], ";")[[1]])
  data.frame(gene = hpa$gene[i], antibody = ids,
             stringsAsFactors = FALSE)
}))
message("[12] antibodies recorded: ", nrow(ab), " across ",
        length(unique(ab$gene)), " genes")

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
message("[12] ==== HPA check ====")
print(hpa[, c("gene", "liver_level", "stomach_level", "specificity",
              "reliability")])

flagged <- hpa$gene[hpa$flag_stomach_contradiction]
if (length(flagged)) {
  warning("[12] stomach_level and stomach_cells contradict each other for: ",
          paste(flagged, collapse = ", "), "\n",
          "  stomach_cells is 'Glandular cells' for every gene, including genes ",
          "whose stomach level is 'Not detected'.\n",
          "  The value is left as curated. Either remove the cell-type string for ",
          "those genes or record why the signal is attribution-free.",
          call. = FALSE)
  message("[12] flagged genes: ", paste(flagged, collapse = ", "))
}

# Note the reliability grade that cuts against its own verdict.
il32 <- hpa[hpa$gene == "IL32", ]
message("[12] IL32 reliability is '", il32$reliability,
        "' while the verdict is '", il32$manuscript_verdict, "'",
        " -- these are consistent only if read as 'protein evidence exists but ",
        "the localisation does not settle the question'")

message("[12] version ", HPA_VERSION, " is stated in the manuscript but was not ",
        "recorded during curation; verify it against the site before submission")

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
write_table(hpa, "p4_HPA_protein_check.csv")
write_table(ab,  "p4_HPA_antibodies.csv")

message("[12] done")

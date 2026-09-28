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
# THREE THINGS THE ORIGINAL DID NOT SETTLE, all handled below:
#
#   1. VERSION. The manuscript cited HPA version 25.0; the author has confirmed
#      nowhere in the code or in the generated table. Worse, it cannot have been
#      the version in force: HPA released 25.1 on 2026-05-25, and the annotations
#      were made on 2026-09-14. So the release actually open during curation was
#      25.1, and the citation names the wrong one. See R/16_verify_hpa_version.R
#      for the check that settles whether it changes any of the five genes.
#
#   2. ACCESS DATE. Three values have been in play - the manuscript says
#      2026-09-05, the supplementary NOTES say 2026-09-13, and the author has since
#      confirmed 2026-09-14. The confirmed date is recorded as such below, and the
#      script reports the disagreement rather than quietly picking one.
#
#   3. ONE CELL CONTRADICTS ITSELF. `stomach_cells` is the constant string
#      "Glandular cells" for all five genes, including CADM2 whose own
#      `stomach_level` is "Not detected". The script flags any gene whose
#      cell-type string asserts a localisation while its level says not detected,
#      so the contradiction surfaces at run time instead of sitting unnoticed in a
#      submitted table. The value itself is NOT changed: silently "fixing" curated
#      data would alter a published table.
#
#      Note in the table's favour: "Glandular cells" IS a standard HPA cell-type
#      group value, appearing verbatim in HPA's own search vocabulary. The
#      attribution is HPA-consistent usage; the problem is that it is applied to a
#      gene whose level reads not detected.
#
# Run:  Rscript R/12_HPA_check.R
# =============================================================================

source("R/00_setup.R")

# --- contested metadata -------------------------------------------------------
# The version cited in the manuscript versus the version actually in force on the
# confirmed access date. Both are carried so the discrepancy is visible rather than
# resolved by assumption.
HPA_VERSION            <- "25.1"          # confirmed by the author; was cited as 25.0
HPA_VERSION_IN_FORCE   <- "25.1"          # released 2026-05-25; live on the access date
HPA_ACCESS_MANUSCRIPT  <- "2026-09-05"    # manuscript Data availability text
HPA_ACCESS_NOTES       <- "2026-09-13"    # supplementary NOTES sheet
HPA_ACCESS_CONFIRMED   <- "2026-09-14"    # author-confirmed actual access date

if (HPA_ACCESS_MANUSCRIPT != HPA_ACCESS_NOTES) {
  warning("[12] HPA access date has three values in play:\n",
          "  manuscript:          ", HPA_ACCESS_MANUSCRIPT, "\n",
          "  supplementary NOTES: ", HPA_ACCESS_NOTES, "\n",
          "  author-confirmed:    ", HPA_ACCESS_CONFIRMED, "\n",
          "  Use the confirmed date in both places.", call. = FALSE)
}

if (HPA_VERSION != HPA_VERSION_IN_FORCE) {
  warning("[12] HPA version cited (", HPA_VERSION, ") is not the version that was ",
          "live on the access date (", HPA_VERSION_IN_FORCE, ", released ",
          "2026-05-25).\n",
          "  Write ", HPA_VERSION_IN_FORCE, ", then compare the five genes across ",
          "both releases at vX.proteinatlas.org\n",
          "  to confirm the annotation is unchanged. ", 
          "See R/16_verify_hpa_version.R.", call. = FALSE)
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
# The version is recorded as the one IN FORCE on the access date, not the one the
# manuscript cites, so the output states what the annotation actually reflects.
# The cited version is kept alongside it so the discrepancy travels with the table.
hpa$hpa_version_used  <- HPA_VERSION_IN_FORCE
hpa$hpa_version_cited <- HPA_VERSION
hpa$access_date       <- HPA_ACCESS_CONFIRMED
hpa$provenance        <- "curator-assigned; no programmatic HPA query"

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

message("[12] version: manuscript cites ", HPA_VERSION,
        "; the release in force on ", HPA_ACCESS_CONFIRMED, " was ",
        HPA_VERSION_IN_FORCE, " (released 2026-05-25).")
message("[12]   the output records BOTH, so the table states what it reflects ",
        "and what it is cited as.")
message("[12]   next: compare the five genes across v25 and 25.1 at ",
        "vX.proteinatlas.org. If unchanged, the citation is the only fix.")
message("[12]   run R/16_verify_hpa_version.R to record what the site reports.")

# -----------------------------------------------------------------------------
# Write
# -----------------------------------------------------------------------------
write_table(hpa, "p4_HPA_protein_check.csv")
write_table(ab,  "p4_HPA_antibodies.csv")

message("[12] done")

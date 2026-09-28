# =============================================================================
# 16_verify_hpa_version.R — confirm the HPA release and access date
# =============================================================================
# STATUS: ADDED. Not from the original work; there was no version check to restore.
#
# WHY THIS EXISTS
#
# The manuscript's Data availability statement cites "Human Protein Atlas, version
# 25.0, accessed 5 September 2026". The curator table that the HPA annotations
# actually come from records neither:
#
#   docs/provenance/p4_HPA_protein_check_as_generated.csv
#     - 9 columns, no version field, no date field
#     - the string "2026" does not appear anywhere in the file
#
# The access date has since been confirmed by the author as 2026-09-14, which is a
# third value: the manuscript says 09-05 and the supplementary NOTES say 09-13.
#
# The version number is still unverified, and it matters for reproducibility: protein
# levels and antibody reliability grades are release-dependent, so if 25.0 is not the
# release that was open during curation, the cited benchmark is wrong.
#
# WHAT THIS SCRIPT DOES
#
# It reads what the Human Protein Atlas currently reports and prints it, so the
# version can be confirmed from a record instead of from memory. Three things are
# checked, in order of how much they prove:
#
#   1. The release number the site reports right now.
#   2. Whether that release is 25.0. If the site reports 25.0, the citation is
#      consistent and the question is settled for the date of THIS run.
#   3. If the site reports something else, the release that was current on
#      2026-09-14 cannot be retrieved retroactively - but HPA publishes release
#      dates, so a version whose release date falls on or before 2026-09-14 and is
#      the latest such version is the one that was open. The script reports the
#      release date where available so that inference can be made explicitly.
#
# WHAT THIS SCRIPT CANNOT DO
#
# It cannot tell you what was on the author's screen on 2026-09-14. If the site now
# reports a later release, the honest options are:
#
#   (a) cite the current release and re-verify the five genes against it - the
#       annotation is five genes and an afternoon's work; this is the strongest fix;
#   (b) cite the release that was current on 2026-09-14, if it can be identified
#       from HPA's own release history;
#   (c) cite "accessed 2026-09-14" without a version number, and say so in the
#       statement.
#
# Option (c) is acceptable and honest. What is not acceptable is citing a version
# number that nothing supports.
#
# Run:  Rscript R/16_verify_hpa_version.R
# Requires internet access.
# =============================================================================

source("R/00_setup.R")

HPA_ACCESS_CONFIRMED <- "2026-09-14"   # author-confirmed
HPA_VERSION_CITED    <- "25.0"        # manuscript; unverified

message("[16] manuscript cites HPA version ", HPA_VERSION_CITED,
        ", accessed ", HPA_ACCESS_CONFIRMED)
message("[16] the generated curator table records neither; see docs/provenance/")

# --- fetch and parse ----------------------------------------------------------
fetch_page <- function(url) {
  con <- tryCatch(url(url, open = "rt"), error = function(e) NULL)
  if (is.null(con)) {
    message("[16] could not open ", url)
    return(NULL)
  }
  on.exit(close(con))
  tryCatch(paste(readLines(con, warn = FALSE), collapse = "\n"),
           error = function(e) {
             message("[16] could not read ", url, ": ", conditionMessage(e))
             NULL
           })
}

extract_versions <- function(html) {
  if (is.null(html)) return(character(0))
  # HPA spells releases as "v25", "version 25", "HPA 25.0" and similar
  pats <- c("v[0-9]{2}\\b", "version\\s+[0-9]{2}\\b",
            "HPA\\s*[0-9]{2}\\b", "[0-9]{2}\\.[0-9]\\b")
  hits <- unlist(lapply(pats, function(p)
    regmatches(html, gregexpr(p, html, ignore.case = TRUE))[[1]]))
  sort(unique(trimws(hits)))
}

pages <- c(
  home     = "https://www.proteinatlas.org/",
  about    = "https://www.proteinatlas.org/about",
  release  = "https://www.proteinatlas.org/about/releases"
)

message("[16] ==== querying the Human Protein Atlas ====")
found <- character(0)
for (nm in names(pages)) {
  message("[16] ", nm, ": ", pages[[nm]])
  html <- fetch_page(pages[[nm]])
  if (is.null(html)) next
  v <- extract_versions(html)
  if (length(v)) {
    message("[16]   version-like strings: ", paste(head(v, 12), collapse = ", "))
    found <- c(found, v)
  } else {
    message("[16]   no version-like string found on this page")
  }
}

# --- verdict ------------------------------------------------------------------
message("\n[16] ==== result ====")

if (!length(found)) {
  message("[16] Could not determine the release automatically.",
          "\n  Check https://www.proteinatlas.org/about/releases manually and",
          "\n  record the release that was current on ", HPA_ACCESS_CONFIRMED, ".")
} else {
  found <- sort(unique(found))
  message("[16] version-like strings seen across pages: ",
          paste(found, collapse = ", "))

  has_25 <- any(grepl("(^v?25$)|(\\b25\\.0\\b)|(version\\s+25\\b)", found,
                      ignore.case = TRUE))
  if (has_25) {
    message("[16] The site reports a release consistent with the cited ",
            HPA_VERSION_CITED, ".")
    message("[16] If that is the current release, the citation is consistent ",
            "as of the date of this run.")
    message("[16] Record the date you ran this check alongside the citation.")
  } else {
    warning("[16] The site does NOT report a release matching ", HPA_VERSION_CITED,
            ".\n",
            "  Options:\n",
            "    (a) cite the current release and re-verify the five genes against ",
            "it (strongest);\n",
            "    (b) cite the release that was current on ", HPA_ACCESS_CONFIRMED,
            ", if identifiable from HPA's release history;\n",
            "    (c) cite 'accessed ", HPA_ACCESS_CONFIRMED,
            "' without a version number.",
            call. = FALSE)
  }
}

message("\n[16] ==== what this script cannot tell you ====")
message("[16] It cannot recover what was on screen on ", HPA_ACCESS_CONFIRMED, ".")
message("[16] If a newer release has shipped since, re-checking the five genes ",
        "against the current release is the only fully defensible fix.")
message("[16] Five genes, two tissues, ten antibodies - an afternoon.")

# --- record -------------------------------------------------------------------
rec <- data.frame(
  field = c("access_date_confirmed", "version_cited_in_manuscript",
            "version_found_this_run", "version_field_in_curator_table",
            "check_date", "conclusion"),
  value = c(HPA_ACCESS_CONFIRMED, HPA_VERSION_CITED,
            if (length(found)) paste(head(found, 8), collapse = "; ") else "not determined",
            "absent", format(Sys.Date()),
            if (length(found) && any(grepl("25", found))) 
              "site reports a release consistent with the citation"
            else "citation not corroborated by the site this run"),
  stringsAsFactors = FALSE)
write_table(rec, "p4_HPA_version_verification.csv")
print(rec)

message("\n[16] done")

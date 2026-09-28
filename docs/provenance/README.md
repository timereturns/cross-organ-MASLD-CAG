# Provenance artifacts

Output files that were produced during the study and were obtained for verification.
They are kept here because they let a reader check what the repository's scripts
reproduce, independently of the scripts themselves.

These are **not** analysis inputs. The pipeline regenerates its own outputs into
`results/tables/`. What is here is the evidence that the reconstruction matches what
was actually generated.

| File | Original name | sha256 (prefix) | What it is |
|---|---|---|---|
| `p4_HPA_protein_check_as_generated.csv` | `p4_HPA_protein_check.csv` | `a6c27264` | The Human Protein Atlas curator table as written by the original script |

---

## `p4_HPA_protein_check_as_generated.csv`

The curator-assigned HPA protein table behind Table S27. Produced by the original
`落盘 最终核对表 2026.9.13 15 12.R`, which built it as a `data.frame` literal and
wrote it out. Restored as `R/12_HPA_check.R`.

### What it confirms

The header is exactly nine columns, and the five data rows are reproduced in
`R/12_HPA_check.R` character for character:

```
gene,liver_level,liver_cells,stomach_level,stomach_cells,
specificity,antibody,reliability,manuscript_verdict
```

So the table in the repository is the table that was generated, not a paraphrase of it.

### What it does not contain — and this matters

Searched for and **not found** anywhere in the file:

| Search | Result |
|---|---|
| `version` | absent |
| `25.0` | absent |
| `2026-09-05` | absent |
| `2026-09-13` | absent |
| `2026` | **absent** |

The file has no version field and no access date, and does not record a year at all.

The manuscript's Data availability statement originally cited **Human Protein Atlas, version 25.0** (since corrected to **25.1**, the release in force on the access date)
and an access date of **5 September 2026**. The supplementary notes say
**2026-09-13**. Neither date, and not the version number, appears in the artifact.

This does not mean the manuscript is wrong — HPA version numbers and access dates are
normally tracked outside the analysis output. It means the claim is **not corroborated by
the generated table**, so it has to be confirmed another way. Two concrete actions:

1. Check the current HPA release number against the site, and against whatever notes or
   browser history exist from the curation session.
2. Choose one access date. The manuscript and the supplementary notes currently disagree,
   and an editor comparing the two will notice.

### Two related points

**The table is hand-curated, and that is stated rather than implied.** Every protein
level, cell-type string, antibody identifier and reliability grade was read off HPA and
typed in. HPA's own reliability grades are themselves the product of human curation, so
this is a legitimate way to build an annotation table — what matters is that a reader can
tell it apart from a computed result. The repository now adds a `provenance` column
saying so, and `R/12_HPA_check.R` prints it.

**One cell contradicts itself.** `stomach_cells` is the constant string `"Glandular cells"`
for all five genes, including CADM2, whose own `stomach_level` is
`"Not detected (raw Low)"`. A gene cannot have a cell type attributed to it while its
level reads not detected. The reconstruction **does not change the value** — silently
"fixing" curated data would alter a published table — but it flags the row at run time so
the inconsistency is visible rather than discovered. See
`../manuscript_discrepancies.md` item 8.

# Repository build notes

How the organised `R/` scripts relate to the original working files, and what
reconstruction involved. For the file-by-file provenance table see
`original_scripts_manifest.md`.

Source material: **23 files, ~8,400 lines**, arriving in two passes. The second pass added
`clump自查 - 8循环重跑9.12.R`, which resolved the MR clumping parameters, and `test.R`,
an unrelated `mtcars` practice script. Three files in the second pass were byte-identical
to files already analysed.

Output: **16 scripts, ~4,300 lines.**

The reduction is almost entirely removed scaffolding, not removed analysis. The largest
single contributor is the figure script: ~1,950 lines containing four superseded revisions
and one abandoned block became ~600 lines containing one version of each figure.

---

## Two kinds of script in this repository

**(a) Faithful restorations.** Most of the pipeline is restored, not rewritten. The
original working files contain the real analysis code; it was extracted, the
RStudio-session scaffolding was removed, and the analysis logic is unchanged. Where a
block depended on objects left in a previous session, the input reads were made explicit
— but the statistics are the ones that ran.

**(b) Re-implementations.** Where the original is unavailable or cannot be redistributed,
the stage is implemented from the public inputs and the documented parameters. **These are
marked in the script header**, so a reader is never left guessing which they are looking
at.

The distinction matters for a reviewer: a restored script can be checked line by line
against the analysis history, whereas a re-implementation has to be validated against the
reported numbers. The header of each affected script says which it is.

---

## Common transformations applied to every restored script

| Transformation | Reason |
|---|---|
| `setwd("E:/生信狂人/…")` removed | 100+ machine-specific absolute paths; one file repeats the same `setwd()` 18 times. Replaced by `R/00_setup.R` and repository-relative paths. |
| Runtime `install.packages()` / `BiocManager::install()` removed | Several blocks install packages mid-analysis, occasionally from a hardcoded CRAN mirror. Dependencies are declared in the README instead. |
| Blocks separated by blank lines merged by analysis stage | The originals are concatenations of RStudio sessions, not pipelines. |
| Superseded revisions dropped | The figure script alone holds four superseded revisions and one block with no output. The MR pipeline exists in three versions. |
| In-session object dependencies made explicit | Some blocks — including the one producing the submitted Fig 3 — reference objects from an earlier session and fail in a fresh one. |
| `set.seed()` added only where a stochastic step previously had **no** seed | One such step was found: a base-R `stripchart(method = "jitter")` whose point positions were irreproducible. The value added is recorded in `00_setup.R`. No existing seed was changed. |

## What was deliberately **not** changed

- **Hand-entered values stay where they are.** The G and H dimensions of the scoring
  matrix, the Human Protein Atlas table, two columns of the four-arm comparison table and
  several rows of the deconvolution QC table are curator-assigned, not computed. They
  remain, and are marked as hand-entered in `manuscript_discrepancies.md`. Undocumented
  manual values are a reproducibility problem; documented ones are a curation record.
  Replacing them with computed values would change the published results, which is not a
  decision for this repository to make.
- **Thresholds that were tightened during the analysis stay tightened.** The RRA score
  cutoff moved from 0.05 to 0.01 in the original code. The repository runs the final
  value and records the sequence in the parameter file.
- **The `theta` definition used for the published tables is the one kept.** Two
  definitions coexist in the working files; see `manuscript_discrepancies.md` §16 for
  which was identified and how.

---

## Verification approach

Each script is syntactically checked with `Rscript -e 'parse(file = …)'`. Where a script
can be executed without its large input data, it is run and its assertions exercised —
`00_setup.R` and the parameter constants are covered this way.

Full end-to-end re-running requires the datasets listed in `data/datasets.md` and is not
attempted in this environment. The scripts therefore carry their own internal consistency
assertions (the originals were unusually good about this — expected sample counts, group
balances, locked ρ values — and those assertions are preserved) so that a user with the
data gets a loud failure rather than a silent wrong answer.

## Known gaps

Listed in `repository_scope.md` §"Not included, and why", with the reason for each. The
one that affects the manuscript's wording is the permutation test.

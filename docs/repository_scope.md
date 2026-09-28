# Repository scope

This repository accompanies the manuscript:

> *Cross-organ transcriptomic integration reveals shared cell-state patterns between
> metabolic dysfunction-associated steatotic liver disease and chronic atrophic
> gastritis* — submitted to Scientific Reports.

It contains the analysis code for the study. This file states plainly what is included,
what is not, and why — so that a reviewer can tell the difference between a gap and a
deliberate boundary.

---

## Included

| Stage | Script |
|---|---|
| Cohort-level differential expression, liver and gastric | `R/01_download_and_DE.R` |
| Robust rank aggregation discovery screen | `R/02_RRA_discovery.R` |
| Cross-organ concordance and candidate finalisation | `R/03_cross_organ_candidates.R` |
| Asset reconstruction (expression matrices, groups, rank columns) | `R/04_assets.R` |
| Robustness audit: leave-one-dataset-out, per-cohort effects, meta-analysis | `R/05_audit_LODO_meta.R` |
| Cell-composition adjustment and batch diagnostics | `R/06_composition_adjustment.R` |
| External validation, held-out cohort GSE135251 | `R/07_validation_GSE135251.R` |
| Single-cell / single-nucleus localisation and pseudobulk | `R/08_singlecell_pseudobulk.R` |
| Deconvolution sensitivity analysis (MuSiC) | `R/09_deconvolution_MuSiC.R` |
| Candidate adjudication: 8-dimension scoring | `R/10_candidate_scoring.R` |
| Exploratory two-sample Mendelian randomisation | `R/11_MR_exploratory.R` |
| Human Protein Atlas annotation check | `R/12_HPA_check.R` |
| Figure generation: Figs 1–4, S1, S2 | `R/13_figures.R` |
| Figure 5 redraw | `R/14_fig5_redraw.R` |
| Analysed result tables (the numbers cited in the paper) | `results/tables/` |

---

## Not included, and why

### 1. Third-party derived code — deliberate exclusion

One working script used during the study carried a vendor notice asserting that use of
the vendor's organised data and code to publish a paper without permission constitutes
unauthorised use. The author's position is that the analysis code is their own work and
the notice is boilerplate attached to a data package; that position has not been
independently confirmed.

**Any part of that file that may derive from a third-party package is therefore not
redistributed here**, and no licence is asserted over it. The corresponding stage is
covered by `R/01_download_and_DE.R`, which is a **fresh implementation** of the same
pipeline from the public GEO source files rather than a copy of the original file. It is
marked as a re-implementation in its own header, and its per-cohort outputs should be
checked against the reported numbers before it is relied on.

This is a scoping decision, not a claim that the excluded material is unavailable. The
author can provide it to the editors and reviewers on request, as the Nature Portfolio
code policy provides.

### 2. Resampling — there is none

This section previously listed a permutation analysis as not included. **That was based on
a misreading and has been corrected: the analysis contains no Monte Carlo step.**

The manuscript Methods refers to *"20,000 rank-permutation replicates"* and Table S02 is
titled *"Leave-one-dataset-out exact RRA (20,000 permutations)"*. Both describe **the same
deterministic procedure**:

- `R/05_audit_LODO_meta.R` recomputes the robust rank aggregation with each of the six
  discovery cohorts removed in turn, using the closed-form P-value implementation
  `min_k pbeta(x_(k); k, n − k + 1)`;
- it produces `results/tables/p0_lodo_exact_rra.csv`, with per-cohort retention columns
  `in200_drop_*` (rank ≤ 200 after dropping that cohort) and `lt001_drop_*` (P < 0.01
  after dropping that cohort);
- **20,000 is the background gene count `N`, not an iteration count.** In the code it is
  the denominator passed to the closed-form function.

Nothing is resampled, so nothing was lost and there is no seed to quote. **The LODO table
is fully reproducible from this repository**, which is a stronger position than having a
seed would be.

The only defect is the label, in three places: the Methods text, the Table S02 title, and
the Code availability statement's promise of "the random seeds used for permutation
testing". None should refer to permutation. See `manuscript_discrepancies.md` item 1 for
the specific rewording.

The one genuinely stochastic procedure in the project is the **AUC bootstrap** for the
GSE135251 validation (2,000 resamples, seed `20260914`), which is fully implemented in
`R/07_validation_GSE135251.R`.

### 3. Preprocessing of single-cell datasets

The scripts that convert the raw single-cell inputs into the compact analysis objects
(`p1_gse202379_compact.rds`) and the deconvolution reference
(`B_reference_merged.rds`) were run in earlier sessions and are not reproduced here. The
derived objects' construction is documented in `docs/analysis_parameters.md`, and the
scripts that consume them are included.

### 4. GSE153224 external validation

The validation statistics for GSE153224 are reported in the manuscript and its
supplementary tables. The script that computed them is not included. The GSE135251
held-out validation — the other half of the external validation — **is** included in
full (`R/07_validation_GSE135251.R`).

### 5. Raw expression data

No expression data are redistributed. All datasets are public; `data/datasets.md` lists
every accession, platform and role. Several of these matrices are tens of gigabytes and
the GEO terms of use do not permit wholesale re-hosting.

### 6. Interactive session leftovers

The original working files contained several blocks that were exploratory rather than
analytical: package installation, filesystem probes, console output read back into the
next session, and superseded drafts of analyses. These are not included. Where a
superseded analysis matters for understanding a decision — for example the three
successive versions of the Mendelian randomisation pipeline, of which only the last
produced the reported results — the decision is documented in
`docs/analysis_parameters.md`.

---

## What a reviewer can verify from this repository

1. That the 14 candidate genes are the output of a stated, reproducible screen
   (`R/02`, `R/03`), including the thresholds applied.
2. That every reported parameter is declared in `R/00_setup.R` rather than buried in a
   script body, and tabulated in `docs/analysis_parameters.md`.
3. That the robustness audit's retention criteria and meta-analytic model are explicit
   (`R/05`).
4. That the Bonferroni threshold quoted in the manuscript, P = 4.8 × 10⁻⁴, follows
   arithmetically from 13 genes × 8 cell types = 104 tests.
5. Which figure came from which script (`figure_map.md`).
6. Which numbers in the supplementary tables are computed by code and which are
   curator annotations transcribed from a run log — see
   `docs/manuscript_discrepancies.md` §14.

## Known limitations, stated rather than hidden

`docs/manuscript_discrepancies.md` records every point at which the code and the
manuscript do not line up, or where two places in the code disagree. That document is
part of this repository on purpose. A reviewer who finds one of those points and sees
that it was already documented, with its resolution, is in a very different position
from one who finds it unaided.

# Manuscript availability statements — ready to paste

Two statements for the manuscript, plus the edits they require elsewhere. Written to
match the Nature Portfolio policy (*Reporting standards and availability of data,
materials, code and protocols*), which *Scientific Reports* follows.

**Placement:** in this order, at the end of the main text —

```
… Discussion / Conclusions
Data availability          ← statement 1
Code availability          ← statement 2
References
```

Both headings must appear verbatim as `Data availability` and `Code availability`.

---

## ⚠️ Three placeholders to fill before pasting

| Placeholder | Replace with | Where it comes from |
|---|---|---|
| `ACCOUNT/REPO` | your GitHub account and repository name | step 1 of the GitHub → Zenodo walkthrough |
| `10.5281/zenodo.XXXXXXX` | the **version DOI** | Zenodo, after the v1.0.0 release is published |
| *(nothing else)* | | |

**Use the version DOI, not the concept DOI.** Zenodo issues two: the version DOI pins the
v1.0.0 snapshot, which is what the analyses were run with, and the concept DOI always
points at the newest version. The manuscript should cite the snapshot.

---

## 1. Data availability

> **Data availability**
>
> All transcriptomic datasets analysed in this study are publicly available in the NCBI
> Gene Expression Omnibus under accessions GSE48452, GSE89632, GSE83452, GSE126848,
> GSE130970 and GSE162694 (liver discovery cohorts), GSE116312 and GSE27411 (gastric
> discovery cohorts), GSE135251 and GSE153224 (external validation cohorts), GSE174478
> (additional liver bulk cohort used for the deconvolution sensitivity analysis), and
> GSE202379, GSE115469 and GSE134520 (single-cell and single-nucleus datasets). Genetic
> instruments were obtained from eQTLGen (https://www.eqtlgen.org) and the GTEx Portal
> (v10; https://gtexportal.org); outcome genome-wide association summary statistics were
> obtained from OpenGWAS/MR-Base (https://gwas.mrcieu.ac.uk). Normal-tissue protein
> annotations were obtained from the Human Protein Atlas, version 25.1
> (https://www.proteinatlas.org; accessed 14 September 2026). No new data were generated
> in this study. Original publications describing the datasets analysed here are cited in
> the reference list.

### What changed from the current draft, and why

| Change | Reason |
|---|---|
| **GSE174478 added** | It is the fourth bulk cohort in the deconvolution sensitivity analysis (Tables S32–S36, Fig S2) and appears in the supplementary NOTES, but the current statement does not list it. Nature Portfolio expects every analysed dataset to be identifiable. |
| **HPA version 25.0 → 25.1** | HPA released 25.1 on 2026-05-25. The annotations were made on 2026-09-14, so 25.1 was the release in force. See `../manuscript_discrepancies.md` item 7. |
| **HPA access date → 14 September 2026** | Author-confirmed. The draft says 5 September 2026; the supplementary NOTES say 13 September 2026. Both need to become 14 September. |
| *"Original publications … are cited in the reference list"* added | Nature Portfolio asks that previously deposited datasets be cited as formal references. The cohorts' source papers are already in the reference list, so this sentence closes the gap without adding 14 entries. |

### Also delete this sentence from the existing Data availability paragraph

> ~~The analysis code is available from the corresponding author on reasonable request.~~

"On reasonable request" is precisely the formulation the Nature Portfolio code policy does
not accept for code central to the conclusions, and the Code availability statement below
supersedes it.

---

## 2. Code availability

> **Code availability**
>
> The R scripts used for cohort-level differential expression, robust rank aggregation,
> the robustness audit (leave-one-dataset-out and effect-size meta-analysis),
> cell-composition adjustment, external validation, patient-level single-cell pseudobulk
> analysis, MuSiC deconvolution, candidate adjudication, the exploratory Mendelian
> randomisation, the Human Protein Atlas annotation check and figure generation are
> available at https://github.com/ACCOUNT/REPO and archived at
> https://doi.org/10.5281/zenodo.XXXXXXX (version v1.0.0). All analytical parameters are
> declared in a single configuration script and documented in the repository, including
> the robust rank aggregation background size, the closed-form P-value implementation,
> the minimum-cell threshold applied in the pseudobulk analysis, the Bonferroni-corrected
> significance threshold (P = 4.8 × 10⁻⁴) and the random seed used for the validation
> bootstrap. The scripts require no proprietary code.

### What changed from the current draft, and why

**The permutation clause is removed.** The draft ends:

> ~~…and the random seeds used for permutation testing.~~

There is no permutation testing in this study. The 20,000 is the robust rank aggregation
background gene count, and the leave-one-dataset-out audit is computed in closed form and
is fully deterministic — six iterations, one per dropped cohort, with no resampling. The
Methods and the Supplementary Table S02 title also describe this as a permutation analysis
and need the same correction. See `../manuscript_discrepancies.md` item 1.

The replacement clause cites the seed that **does** exist: the 2,000-resample percentile
bootstrap for the GSE135251 AUC confidence intervals (`set.seed(20260914)`).

**Two further points of wording:**

- *"All analytical parameters were fixed in advance"* has been replaced with *"All
  analytical parameters are declared in a single configuration script and documented in
  the repository."* The repository records that several thresholds were adjusted during
  the analysis — the RRA score cutoff moved from 0.05 to 0.01, and the LODO background of
  20,000 was selected after comparing three definitions. Claiming they were fixed in
  advance would be inaccurate and checkable. Claiming they are declared and documented is
  accurate and is what a reviewer actually needs.

- *"The scripts require no proprietary code"* was added because the code policy allows an
  editor to decline a paper if important code is unavailable. A one-line statement that
  nothing is blocked behind a licence is cheap insurance.

---

## 3. One item to confirm before submission

The Code availability statement points at the repository. The repository documents the
following, and each should be consistent with what the manuscript says:

| Repository documents | Manuscript should say |
|---|---|
| RRA discovery background = runtime union of symbols; audit background = 20,000 | whichever the Methods describes — see item 3 in `../manuscript_discrepancies.md` |
| Bonferroni family = 13 genes × 8 cell types, the 8 being the unified liver lineage vocabulary | name the vocabulary, so "8 cell types" is checkable |
| Deconvolution: full-gene MuSiC failed pre-registered QC; marker-restricted fallback used, ≤ 400 markers at 50 per cell type | consistent with Tables S32–S36 |
| MR: 69 estimates, none surviving FDR, exploratory only | consistent with the Results |
| GSE153224 validation script not included in the repository | the statement points at the repository as a whole; if an editor asks for that step specifically, supply it directly |

---

## 4. Reviewer-facing note

Nothing in either statement claims the repository is complete. The Code availability
statement says where the scripts are, what they cover, and that the parameters are
documented. `docs/repository_scope.md` states plainly which stages are not included and
why — the preprocessing of the single-cell objects, the GSE153224 validation script, and
any third-party-derived material. That is a scope statement, not a hidden gap, and it is
better that a reviewer reads it from the repository than infers it.

## 5. Submission checklist

Carried over from the handover note, unchanged:

- **STROBE-MR checklist** must be uploaded as a separate file. *Scientific Reports* states
  that submissions without it will not be considered.
- Repository must be **public** before the DOI exists, so that the link in the Code
  availability statement resolves for reviewers.

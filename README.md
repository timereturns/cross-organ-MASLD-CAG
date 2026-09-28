# cross-organ-MASLD-CAG

Analysis code for the manuscript:

> **Cross-organ transcriptomic integration reveals shared cell-state patterns between metabolic dysfunction-associated steatotic liver disease and chronic atrophic gastritis**
> Cai J, Sentalati B, Fan Z, Chen Y\*, Zeng B\*
> Xinjiang Medical University, School of Traditional Chinese Medicine
> Submitted to *Scientific Reports*

The study integrates six liver and two gastric transcriptomic cohorts, identifies
candidate genes shared between MASLD and chronic atrophic gastritis by robust rank
aggregation, adjudicates them through a robustness audit, replicates them in two held-out
cohorts, localises them by single-cell and single-nucleus analysis, and reports an
exploratory two-sample Mendelian randomisation. This repository contains the R code for
that pipeline, the analysed result tables, and a full statement of the analysis
parameters.

**Current release:** v1.0.0 · **DOI:** see the Zenodo badge above (added at release)

> **A note on the manuscript's wording.** The Methods section describes "20,000
> rank-permutation replicates". There is no permutation step in this analysis: 20,000 is
> the robust rank aggregation background gene count, and the leave-one-dataset-out audit
> is computed in closed form. The statistics are correct and fully reproducible here; the
> label is not. See §6 and `docs/manuscript_discrepancies.md` item 1.

---

## 1. What is here, and what is not

Read `docs/repository_scope.md` before assuming a gap is an omission. In summary:

**Included** — the differential-expression pipeline, the rank-aggregation discovery
screen, the cross-organ candidate finalisation, the leave-one-dataset-out robustness
audit, the composition adjustment, the GSE135251 held-out validation, the single-cell
pseudobulk analysis, the MuSiC deconvolution sensitivity analysis, the candidate
adjudication, the exploratory MR, the Human Protein Atlas check, and the figure code.

**Not included** — raw expression data (public, and too large to re-host); the
single-cell preprocessing that produced the compact analysis objects; the GSE153224
validation script; and the standalone permutation script (see §6).

---

## 2. Running order

Every script sources `R/00_setup.R`, which defines the repository root, the output
directories and all analysis parameters. Run from the repository root:

```bash
Rscript R/01_download_and_DE.R
Rscript R/02_cross_organ_candidates.R
Rscript R/03_assets.R
Rscript R/04_audit_LODO_meta.R
Rscript R/05_composition_adjustment.R
Rscript R/06_validation_GSE135251.R
Rscript R/07_singlecell_pseudobulk.R
Rscript R/10_candidate_scoring.R
Rscript R/11a_instrument_clumping.R # before 11: produces the exposure instrument files
Rscript R/11_MR_exploratory.R       # before 10 if regenerating dimension F
Rscript R/12_HPA_check.R
Rscript R/13_figures.R
Rscript R/14_fig5_redraw.R          # Fig 5, from the authoritative stage table
Rscript R/16_verify_hpa_version.R   # HPA release check (needs internet; optional)
Rscript R/17_bonferroni_family.R    # derives and asserts the 104-test family
```

| Script | Depends on |
|---|---|
| `01` | internet access; writes `data/<accession>/` |
| `02` | `01` |
| `03` | `01`, `02` |
| `04` | `01`, `02` |
| `05` | `03`, `04` |
| `06` | `04` |
| `07` | `04` |
| `08` | `data/GSE202379`, `data/GSE115469`, `data/GSE134520` |
| `09` | `04` and the deconvolution assets (see §6) |
| `10` | `03`, `05`, `07`, `08`, `11` |
| `11a` | internet access (IEU OpenGWAS API); the purchased r² = 0.01 instrument set |
| `11` | `11a`; the GWAS outcome files |
| `12` | nothing (curated table) |
| `13` | `05`, `07`, `08`, `09` |
| `14` | `08` and the patient covariate table |

Note the ordering constraint around the MR: run `11a` before `11`, and `11` before `10`
if you are regenerating the scoring matrix, because dimension F of that matrix reads the
MR results.

A full re-run needs roughly 200 GB of free disk for the raw inputs, most of it the
GSE135251 per-sample count files and the three single-cell datasets.

Set `REPO_ROOT` if you run the scripts from elsewhere:

```r
Sys.setenv(REPO_ROOT = "/path/to/cross-organ-MASLD-CAG")
source("R/00_setup.R")
```

---

## 3. Environment

| | |
|---|---|
| R | **4.6.0** (2026-04-24 ucrt) — the version stated in the manuscript |
| Platform | x86_64-w64-mingw32/x64 (Windows 11) |
| Locale | Chinese (Simplified)\_China.utf8 · tz Asia/Shanghai |

| Package | Version | Used for |
|---|---|---|
| `limma` | 3.68.5 | differential expression |
| `RobustRankAggreg` | 1.2.1 | rank aggregation |
| `metafor` | 5.0-1 | random-effects meta-analysis |
| `MuSiC` | 1.0.0 | deconvolution |
| `Seurat` | 5.5.1 | single-cell handling |
| `TwoSampleMR` | 0.7.9 | Mendelian randomisation |
| `pheatmap` | 1.0.13 | heatmaps |
| `ComplexHeatmap` | 2.28.0 | heatmaps |
| `SingleCellExperiment` | 1.34.0 | MuSiC input |
| `Biobase` | 2.72.0 | MuSiC input |
| `edgeR` | 4.10.5 | RNA-seq normalisation |

`MuSiC` is **not on Bioconductor**; v1.0.0 was installed from
`github.com/xuranw/MuSiC` using the `bulk.mtx` + `SingleCellExperiment` interface.

Annotation packages: `org.Hs.eg.db`, `hugene11sttranscriptcluster.db`,
`hugene20sttranscriptcluster.db`, `illuminaHumanWGDASLv4.db`,
`hugene10sttranscriptcluster.db`, `illuminaHumanv3.db`.

Also required for parts of the pipeline: `ieugwasr`, `coloc`, `openxlsx`, `nnls`,
`quadprog`, `patchwork`, `ggtext`, `readxl`, `png`, `dplyr`.

Full `sessionInfo()` output: `sessionInfo.txt`. Package inventory at capture time:
`docs/installed_packages.csv`.

---

## 4. Data

No expression data are redistributed. Everything is public; `data/datasets.md` lists
every accession, platform, role and download location.

| Role | Accessions |
|---|---|
| Liver discovery (6) | GSE48452, GSE89632, GSE83452, GSE126848, GSE130970, GSE162694 |
| Gastric discovery (2) | GSE116312, GSE27411 |
| External validation (2) | GSE135251, GSE153224 |
| Single-cell / single-nucleus (3) | GSE202379, GSE115469, GSE134520 |
| **Deconvolution bulk cohort (4th)** | **GSE174478** |

Non-GEO sources: eQTLGen (`https://www.eqtlgen.org`), GTEx v10
(`https://gtexportal.org`), OpenGWAS / MR-Base (`https://gwas.mrcieu.ac.uk`),
eQTL Catalogue tabix paths (for the colocalisation follow-up), and the Human Protein
Atlas v25.1 (`https://www.proteinatlas.org`).

Layout expected by the scripts:

```
data/
├── GSE48452/GSE48452_series_matrix.txt.gz
├── GSE89632/GSE89632_series_matrix.txt.gz
├── GSE83452/GSE83452_series_matrix.txt.gz
├── GSE126848/GSE126848_series_matrix.txt.gz
│             GSE126848_Gene_counts_raw.txt.gz
├── GSE130970/GSE130970_series_matrix.txt.gz
│             GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz
├── GSE162694/GSE162694_series_matrix.txt  GSE162694_raw_counts.csv.gz
├── GSE116312/GSE116312_series_matrix.txt.gz  GPL6255.annot.gz
├── GSE27411/GSE27411_series_matrix.txt.gz  GSE27411_non-normalized.txt.gz
├── GSE135251/GSE135251_series_matrix.txt.gz  GSE135251_counts/  (216 .gz files)
├── GSE153224/GSE153224_mRNA_Expression_Profiling.xlsx
├── GSE202379/  GSE115469/  GSE134520/  GSE174478/
└── gwas/
    ├── GCST90267352.h.tsv.gz    LiverPDFF_UKB        (main outcome)
    ├── GCST90129440.h.tsv.gz    Gastritis_UKB        (main outcome)
    ├── GCST90054782_mr_ready.tsv.gz  NAFLD_UKB       (sensitivity)
    ├── finngen_R13_NAFLD.gz
    ├── finngen_R13_K11_CHRONGASTR.gz
    ├── finngen_R13_K11_FIBROCHIRLIV.gz
    └── finngen_R13_C3_STOMACH_WIDE.gz
```

| Source | Where | Accessed |
|---|---|---|
| GEO accessions | https://ftp.ncbi.nlm.nih.gov/geo/series/ | _(fill in per cohort)_ |
| eQTLGen cis-eQTLs | https://www.eqtlgen.org | _(fill in)_ |
| GTEx v10 | https://gtexportal.org | _(fill in)_ |
| OpenGWAS / MR-Base | https://gwas.mrcieu.ac.uk | _(fill in)_ |
| Human Protein Atlas | https://www.proteinatlas.org | **2026-09-14** (author-confirmed); version **25.1** — see §6 |

---

## 5. Key parameters

The complete table is `docs/analysis_parameters.md`. The values a reviewer is most
likely to ask about:

| Analysis | Parameter | Value |
|---|---|---|
| **RRA discovery** | package | `RobustRankAggreg` 1.2.1 |
| | background size | **runtime union of symbols across the six liver cohorts** (`N <- length(all_sym)`); **20,000** was used for the audit recomputation (§6) |
| | P-value implementation | closed form, `min_k pbeta(x_(k); k, n − k + 1)` |
| | score threshold | < 0.01, after ≥ 3 of 6 cohorts agreeing in direction |
| | carried forward | top 200 up + top 200 down |
| **Cross-organ finalisation** | gastric DEG rule | `P.Value < 0.05 & abs(logFC) > 0.5` |
| | gastric requirement | significant in **both** gastric cohorts, same direction |
| | liver requirement | **≥ 4 of 6** cohorts in the stated direction |
| **Robustness audit** | LODO | closed-form exact RRA `min_k pbeta(x_(k); k, n−k+1)`, recomputed with each of 6 cohorts dropped |
| | retention | rank ≤ 200, or P < 0.01, after the drop |
| | DT side correction | none — descriptive |
| | **there is no permutation step** | 20,000 is the RRA **background N**, not a replicate count — see §6 |
| **Meta-analysis** | model | `metafor::rma(yi = logFC, sei = SE, method = "REML")` |
| | SE derivation | `SE = logFC / t` |
| | pooling unit | within organ, within platform (chip vs RNA-seq) |
| **Composition adjustment** | cell types | 9 liver lineages, marker-based z-scores |
| | covariates | top 3 cell types by variance |
| | model | `limma ~ group + score_1 + score_2 + score_3` |
| **Pseudobulk** | minimum cells | **≥ 20 cells per patient × cell type** |
| | normalisation | CPM per patient–cell-type unit |
| | fibrosis stage | **mode** across that patient's cells |
| **Bonferroni** | family | 13 genes × 8 cell types = **104 tests** |
| | threshold | **P = 4.8 × 10⁻⁴** (0.05 / 104) |
| **MuSiC** | version | 1.0.0 from GitHub |
| | reference | GSE202379, 8 cell types, patient × cell-type mean pseudobulk |
| | marker panel | specificity `theta / rowMeans(theta)`, top 50 per cell type, ≤ 400 genes total |
| | anti-circularity | 6 candidates and MT-/MTRNR genes zeroed before ranking; absence hard-asserted |
| | bulk cohorts | GSE126848, GSE130970, GSE162694, GSE174478 |
| | candidates tested | CADM2, ANXA4, LGALS3 |
| **MR** | instrument filter | per-SNP F > 10 |
| | **clumping** | **window 10,000 kb · r² < 0.001 · p < 5×10⁻⁸ · EUR panel**, via `ieugwasr::ld_clump` |
| | independence check | `ld_matrix`, max r² asserted ≤ 0.001 |
| | harmonisation | `harmonise_data(action = 2)` — palindromic SNPs dropped |
| | Steiger | **not applied** in the reported version |
| | methods | Wald ratio (1 SNP) or IVW |
| | MR-PRESSO | `NbDistribution = 1000` |
| | estimates | 69; **none survived FDR** |
| **Multiple testing** | Fig 1 | BH over 140 gene–cohort tests |
| | validation | BH over the 5 pre-specified candidates |
| | pseudobulk | Bonferroni over 104 tests |

---

## 6. Known limitations

Stated here rather than left to be discovered. `docs/manuscript_discrepancies.md` gives
the full record with code locations.

**"20,000 permutations" is a misnomer — there is no permutation analysis.** The manuscript
Methods and the Table S02 title both describe the leave-one-dataset-out audit as a
permutation test with 20,000 replicates. In fact **20,000 is the robust rank aggregation
background size**, and the audit uses the closed-form exact P-value
`min_k pbeta(x_(k); k, n − k + 1)`, run once per dropped cohort — six runs, not 20,000.
Nothing is resampled, so nothing was lost and there is no seed to quote. **The LODO table
is fully reproducible from this repository.**

The label needs fixing in three places: the Methods wording, the Table S02 title, and the
Code availability statement's promise of "the random seeds used for permutation testing".
See `docs/manuscript_discrepancies.md` item 1.

The one genuinely stochastic procedure is the **AUC bootstrap** for GSE135251 (2,000
resamples, seed `20260914`), implemented in full in `R/07_validation_GSE135251.R`.

**Two background sizes appear in the RRA.** The discovery screen used
`N <- length(all_sym)`, the runtime size of the union of symbols across the six liver
cohorts. The value **20,000** belongs to the later audit recomputation, where three
definitions were compared and 20,000 was selected by hand. Both are documented; the
manuscript should be read with that distinction in mind.

**Several values in the supplementary tables are transcribed rather than computed.**
The deconvolution QC table (Table S35) contains numbers recorded from an interactive run
log, and two columns of the four-arm comparison table (`verdict`, `note`) were appended
by hand. `docs/manuscript_discrepancies.md` §13–14 lists which rows are computed and
which are transcripts.

**The Human Protein Atlas citation names the wrong release.** The manuscript cites
version 25.0. HPA released **25.1 on 2026-05-25**, and the annotations were made on the
author-confirmed access date of **2026-09-14** — so the release in force during curation
was 25.1. Fix the citation, then compare the five genes across `v25` and 25.1 at
`vX.proteinatlas.org`: if the tissue levels, cell-type attributions and antibody
reliability grades are unchanged, the citation is the only fix needed. The access date
also needs correcting — the manuscript says 2026-09-05 and the supplementary notes say
2026-09-13, against a confirmed 2026-09-14. `R/16_verify_hpa_version.R` records the check.
See `docs/manuscript_discrepancies.md` item 7.

**The Human Protein Atlas table is manually curated.** `R/12_HPA_check.R` holds the
protein levels, cell-type strings, antibody IDs and reliability grades as typed values;
there is no programmatic HPA query. That is a legitimate way to build an annotation table —
HPA's own reliability grades are human curation too — and the output now carries provenance
columns saying so.

**Some single-cell preprocessing is not reproduced here.** The compact analysis object
for GSE202379 and the deconvolution reference were produced in earlier sessions; the
scripts that consume them are included, the scripts that build them are not.

**The GSE153224 validation script is not included.** Its results are reported in the
manuscript and supplements.

**Figure numbering in the scripts differs from the manuscript.** See `figure_map.md`;
three figures were renamed by hand and Fig 5 came from a second script.

**One further accession is analysed that the Data availability statement does not
list:** GSE174478, the fourth bulk cohort in the deconvolution sensitivity analysis.
See `docs/manuscript_discrepancies.md`.

---

## 7. Repository contents

```
├── README.md                      this file
├── LICENSE                        MIT, code only
├── sessionInfo.txt                R 4.6.0 and every package version
├── figure_map.md                  manuscript figure -> script -> output file
├── R/
│   ├── 00_setup.R                 paths, parameters, seeds, shared helpers
│   ├── 01_download_and_DE.R       GEO download; limma per cohort  [re-implemented]
│   ├── 02_RRA_discovery.R         rank aggregation discovery screen
│   ├── 03_cross_organ_candidates.R  gastric concordance; the 14 candidates
│   ├── 04_assets.R                expression matrices; group vectors; rank columns
│   ├── 05_audit_LODO_meta.R       LODO; effect sizes; meta-analysis
│   ├── 06_composition_adjustment.R  marker-score adjustment; batch QC
│   ├── 07_validation_GSE135251.R  held-out validation; AUC; bootstrap CI
│   ├── 08_singlecell_pseudobulk.R GSE202379 / GSE115469 / GSE134520
│   ├── 09_deconvolution_MuSiC.R   bulk gate; full-gene MuSiC; marker panel;
│   │                              four-arm comparison
│   ├── 10_candidate_scoring.R     8-dimension adjudication
│   ├── 11a_instrument_clumping.R  cis-eQTL instrument selection and independence QC
│   ├── 11_MR_exploratory.R        two-sample MR (the reported v4 configuration)
│   ├── 12_HPA_check.R             protein annotation table
│   ├── 13_figures.R               Figures 1-4, S1, S2
│   ├── 14_fig5_redraw.R           Figure 5, from the authoritative stage table
│   ├── 16_verify_hpa_version.R    verifies the cited HPA release against the site
│   └── 17_bonferroni_family.R     derives and asserts the 104-test family
├── data/
│   └── datasets.md                accessions, platforms, dates
├── results/
│   ├── tables/                    analysed result tables cited in the paper
│   └── figures/                   submitted figure files
├── cache/                         large intermediate objects (not committed)
└── docs/
    ├── analysis_parameters.md     the complete parameter set
    ├── build_notes.md             how these scripts relate to the originals
    ├── code_inventory_raw.md      audit of the original working files
    ├── manuscript_discrepancies.md  code vs manuscript, item by item
    ├── original_scripts_manifest.md  provenance of each original file
    ├── repository_scope.md        what is in, what is out, and why
    └── provenance/                generated artifacts kept for verification
        ├── README.md
        └── p4_HPA_protein_check_as_generated.csv
```

Only `R/01_download_and_DE.R` is a re-implementation; every other script is a
restoration of the original analysis logic. The distinction is stated in each
script's header, and `docs/build_notes.md` explains what reconstruction involved.

---

## 8. Figure and table provenance

`figure_map.md` maps every manuscript figure to the script and output file that
produced it, including the four superseded revisions left in the figure script and the
one abandoned block.

`docs/analysis_parameters.md` §13 maps the manuscript's claims to the code that
generates them, and lists every point where the two need reconciling.

## 9. Citation

If you use this code, please cite the manuscript and the archived release:

```bibtex
@software{cai2026crossorgan,
  author  = {Cai, Junpu and Sentalati, Bizha and Fan, Zhengyang and
             Chen, Yanzhu and Zeng, Binfang},
  title   = {Analysis code for "Cross-organ transcriptomic integration reveals
             shared cell-state patterns between metabolic dysfunction-associated
             steatotic liver disease and chronic atrophic gastritis"},
  year    = {2026},
  version = {v1.0.0},
  doi     = {10.5281/zenodo.XXXXXXX},
  url     = {https://github.com/ACCOUNT/REPO}
}
```

## 10. Licence

MIT for the code — see `LICENSE`. The licence does not extend to any third-party
expression data analysed with it; all datasets remain subject to the terms of use of
their repositories, and none are redistributed here.

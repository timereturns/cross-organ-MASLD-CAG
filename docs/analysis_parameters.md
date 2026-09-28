# Analysis parameters

Every parameter a reader, reviewer or editor would need in order to reproduce the
analysis. Values here were read out of the working scripts, not restated from the
manuscript; where the two differ the difference is flagged and the reconciliation is
listed at the end.

Environment: R **4.6.0** (2026-04-24 ucrt). Full package versions: `sessionInfo.txt`.

Items still needing author confirmation are marked **`[CONFIRM]`**. Items the code
settles but where the manuscript disagrees are marked **`[MANUSCRIPT]`**.

---

## 1. Cohort-level differential expression

| Item | Value | Source |
|---|---|---|
| Method | `limma` moderated *t*-statistic, per cohort | `第一步数据下载（MASLD侧）`, `第一步数据下载（CAG）` |
| Array design | `model.matrix(~ 0 + grp)`, contrast `case − control` | `run_limma_fold` / `run_limma_symbol` |
| RNA-seq design | `voom` then limma — `res126 <- voom_limma(h126, grp126, "ENSEMBL")` | `第一步数据下载（MASLD侧）` line 599 |
| Multiple testing | Benjamini–Hochberg, per cohort, on the full gene table | `topTable(..., number = Inf)` |
| Probe → symbol collapse | `mapIds(annotation_db, keytype = "PROTEID")`, then **keep the most significant row per symbol** (`tt[!duplicated(tt$symbol), ]` after sorting by `-abs(t)`) | `run_limma_fold` line 119–120 |
| Covariates | **None.** The design matrices contain only the group term | both DE scripts |
| Genes entering integration | Genes present in all six liver cohorts — `Reduce(intersect, names)` | `磁盘层执行脚本` line 61 |

Annotation databases used: `hugene11sttranscriptcluster.db` (Ahrens/GSE48452),
`illuminaHumanWGDASLv4.db` (Arendt/GSE89632), `hugene20sttranscriptcluster.db`
(Lefebvre/GSE83452), `org.Hs.eg.db` (RNA-seq cohorts and gastric),
`hugene10sttranscriptcluster.db` + `illuminaHumanv3.db` (gastric).

### Cohort sizes and grouping rules

| Cohort | Accession | Samples | Case / control | Grouping rule |
|---|---|---|---|---|
| Ahrens | GSE48452 | 54 | 26 / 28 | `group ∈ {Nash, Steatosis}` → case; **`bariatric surgery == "after surgery"` samples excluded** |
| Arendt | GSE89632 | 63 | 39 / 24 | `diagnosis ∈ {NASH, SS}` → case |
| Lefebvre | GSE83452 | 148 | 104 / 44 | **`time == "baseline"` only**; `liver status == NASH` → case |
| Suppli | GSE126848 | 57 | 31 / 26 | `disease ∈ {healthy, obese}` → control |
| Hoang | GSE130970 | 78 | 53 / 25 | `fibrosis stage == "0"` → control |
| Pantano | GSE162694 | 143 | 112 / 31 | title suffix `_N` → control |
| GSE135251 | validation | 216 | 206 / 10 | `disease == "Control"` → control |
| GSE116312 | gastric | 10 retained | 3 / 7 | CAG → case; **follicular gastritis → control; gastric cancer (3) excluded** |
| GSE27411 | gastric | 12 retained | 6 / 6 | Atrophy → case; **H. pylori non-infected → control; H. pylori-infected (6) excluded** |

Group vectors carry names and are asserted against expected counts
(`stopifnot(sum(v == "case") == exp_case[[nm]], ...)`), so a mis-assignment fails loudly.

---

## 2. Robust rank aggregation — the discovery screen

| Item | Value |
|---|---|
| Package | `RobustRankAggreg` **1.2.1** |
| Call | `aggregateRanks(glist = rank_up, N = N)` |
| Ranking input | full per-cohort symbol ranking by moderated *t*; `rank_up` uses `order(-t)`, `rank_dn` uses `order(t)`; duplicates removed with `unique()` before ranking |
| **Background size, discovery run** | **`N <- length(all_sym)`** — the runtime size of the union of symbols across the six liver cohorts. **Not** a hardcoded 20,000. |
| **Background size, audit recomputation** | **20,000**, selected by hand: `use_N20000 <- TRUE  # 校准后手工指定` |
| Third background size examined | 11,982 = the six-cohort common-gene count, plus per-cohort sizes — compared in `P0组成校正` lines 400–426 |
| P-value implementation | Closed form: `min_k pbeta(x_(k); k, n − k + 1)` where `x = sort(rank / N)` — `磁盘层执行脚本` line 36 |
| Score threshold, first pass | `Score < 0.05` |
| Score threshold, final | `Score < 0.01` |
| Direction-consistency filter | retained only if **≥ 3 of 6** liver cohorts agree in direction |
| Truncation | top **200** up + top **200** down → `MASLD_RRA_up_final.csv`, `MASLD_RRA_dn_final.csv` |
| Sensitivity check | full ranking re-run with Arendt removed; top-50 overlap reported |
| Random seed | **None needed** — RRA is deterministic and no subsampling is applied |

**`[MANUSCRIPT]`** The manuscript states the background size as 20,000. The discovery
screen used a runtime value. See item 3.

### Cross-organ candidate finalisation

| Step | Rule | Source |
|---|---|---|
| Gastric DEG definition | `P.Value < 0.05 & abs(logFC) > 0.5`, per gastric cohort | `第一步数据下载（CAG）` line 280 |
| Concordance with liver | MASLD RRA list ∩ gastric DEGs, direction-matched | line 290–293 |
| Gastric internal requirement | gene must be significant in **both** gastric cohorts, same direction (`intersect(d312$up, d27$up)`) | line 356–357 |
| Liver internal requirement | **≥ 4 of 6** liver cohorts agree in direction | line 375–376 |
| Output | `candidate_genes_final.csv`, 14 genes with direction | line 380–382 |
| Gastric supporting threshold (audit) | `P.Value < 0.05 & abs(logFC) > log2(1.5)` | `磁盘层执行脚本` line 91 |

---

## 3. Robustness audit

### 3.1 Leave-one-dataset-out

| Item | Value |
|---|---|
| Procedure | Full re-ranking with each of the 6 liver cohorts dropped in turn; the aggregation is recomputed on the reduced set |
| Implementation | `apply(mat[, -j], 1, rra_rho, denom = den_lodo(ncol(mat) - 1))` — a genuine recomputation, not rank dropping |
| Denominator | `use_N20000` → `rep(20000, k)` |
| Retention criterion A (`in200`) | rank ≤ **200** after dropping that cohort |
| Retention criterion B (`lt001`) | corrected P < **0.01** after dropping that cohort |
| Output | `p0_lodo_exact_rra.csv`, with `rank_drop_*`, `in200_drop_*`, `lt001_drop_*` per cohort |
| Alternative-background LODO | `p0_lodo_rra_percohortN.csv`, using per-cohort gene counts as denominators |
| Package cross-check | `p0_lodo_aggregateRanks.csv`, recomputed with the `RobustRankAggreg` package rather than the closed-form function |

### 3.2 ✅ RESOLVED — there is no permutation test, and none is needed

**This section previously described a missing permutation script. That was wrong.** The
analysis contains no Monte Carlo step at all, and nothing has been lost.

The LODO audit is **deterministic and closed-form**:

```r
# 磁盘层执行脚本2026.9.12 16 56.R, lines 31-37
rra_rho <- function(r, denom) {
  r <- as.numeric(r); r <- r[!is.na(r)]
  if (length(r) < 2) return(NA_real_)
  if (length(denom) == 1) denom <- rep(denom, length(r))
  x <- sort(r / denom); n <- length(x)
  min(vapply(seq_len(n), function(k) pbeta(x[k], k, n - k + 1), numeric(1)))
}
```

Every quantity it produces is a deterministic function of the input ranks. There is no
resampling, so there is nothing to seed, and **the LODO table is fully reproducible from
the code in this repository** — which is a stronger position than having a seed to quote.

| Item | Value |
|---|---|
| Monte Carlo steps in the pipeline | **none** |
| Resampling anywhere in the original files | only `第一次预审稿文件处理.R` L255–257 and L408–409, inside the **AUC bootstrap** (`B <- 2000`) |
| Seed attached to that bootstrap | `set.seed(20260914)` (L206, L352) |
| `20000` in the LODO code | the **background gene count**, passed as the denominator at L45 and L52 — not an iteration count |
| Retention tallies | `in200_drop_*` (L76) and `lt001_drop_*` (L77), counting LODO runs, not resamples |

### ⚠️ The "20,000 permutations" label is a misnomer — and it is in the manuscript

The manuscript's Methods describes **"20,000 rank-permutation replicates"**, and
Supplementary Table S02 is titled *"Leave-one-dataset-out exact RRA (20,000
permutations)"*.

Both are describing the LODO audit. **Neither is performing permutations.** The 20,000 is
the RRA background size `N`, and the procedure is the closed-form exact aggregation, run
once per dropped cohort — six runs, not 20,000.

This is a labelling problem, not an analysis problem: the statistics are correct and
reproducible. But it must be fixed in three places, because a reviewer who looks for a
permutation step will not find one and will reasonably conclude the code is missing:

1. **Methods** — replace "20,000 rank-permutation replicates" with a description of the
   LODO audit: the exact robust rank aggregation, recomputed with each discovery cohort
   removed in turn, using a background of 20,000 genes.
2. **Table S02 title** — *"Leave-one-dataset-out exact RRA (20,000 permutations)"* →
   something like *"Leave-one-dataset-out exact RRA (background N = 20,000)"*.
3. **Code availability statement** — remove "the random seeds used for permutation
   testing". No permutation testing is performed, so this promise cannot be honoured and
   does not need to be. Replace it with the parameters that do matter: the RRA background
   size, the closed-form P-value implementation, the minimum-cell pseudobulk threshold and
   the Bonferroni threshold — all of which are documented and reproducible.

The full statement of what happened is in `manuscript_discrepancies.md` item 2.

### 3.3 Effect-size meta-analysis

| Item | Value |
|---|---|
| Model | Random effects, **REML** |
| Call | `metafor::rma(yi = logFC, sei = SE, data = s, method = "REML")` |
| Effect size | Per-cohort `logFC` from limma |
| Standard error | **`SE = logFC / t`** (back-computed from the moderated *t*) |
| Heterogeneity | `I2` from the `rma` object |
| Pooling unit | **Within platform, within organ**: `chip` = {Ahrens, Arendt, Lefebvre}; `rna` = {Suppli, Hoang, Pantano} |
| Minimum cohorts | pooled only if `nrow(s) >= 2` |
| Gastric pooling (audit) | separate fixed-effect and DerSimonian–Laird random-effect estimates over the 2 gastric cohorts, with I²; RE with k = 2 reported as a sensitivity estimate |
| Output | `p0_meta_per_platform.csv` |

### 3.4 Alternative cell-composition adjustment

| Item | Value |
|---|---|
| Purpose | Check whether candidate-gene differential expression survives adjustment for cell-type composition |
| Marker panel | **9 liver cell types**, hardcoded in `mk_liver`: Hepatocyte (ALB, CYP3A4, CYP2E1, HNF4A, ASGR1, TTR, FABP1), Cholangiocyte (KRT19, KRT7, SOX9, SCTR), Kupffer_macro (CD163, CD68, MARCO, C1QA, C1QB, MS4A7), T_cell (CD3D, CD3E, CD2, CD8A, IL7R, TRAC, GZMA), B_plasma (MS4A1, CD79A, CD79B, IGKC, IGHG1, MZB1), NK (KLRD1, NKG7, KLRF1), Neutrophil (CSF3R, S100A8, S100A9, FCGR3B), HSC_fibro (ACTA2, COL1A1, COL3A1, PDGFRB, LUM, DCN), Endothelial (PECAM1, VWF, CDH5, ENG, KDR) |
| Composition score | z-score each gene across samples, then mean z per cell type: `score_mk()` = `colMeans(z[g, ])`, requiring ≥ 2 markers present |
| Covariates chosen | the **top 3 cell types by variance** across samples: `top3 <- names(sort(vv, decreasing = TRUE))[1:3]` |
| Adjusted model | `model.matrix(~ g + sc[, top3[1]] + sc[, top3[2]] + sc[, top3[3]])`, fitted with `limma::lmFit` + `eBayes` |
| Raw comparison | same genes with `model.matrix(~ g)` only |
| Batch diagnostic | PCA on variance-filtered genes (`var > 0`), first 3 PCs regressed on group; per-sample PCA group P reported |
| Reported attenuation | `attenuation = 1 − logFC_adj / logFC_raw`; `direction_flip` = sign change |
| Genes/cell types analysed | **6** candidates × up to 9 cell types × the cohorts where the matrix was available |
| Outputs | `p0_composition_summary.csv` (arrays), `p0_composition_summary_rna.csv`, `p0_composition_compdiff_rna.csv`, `p0_composition_res.rds`, `p0_composition_res_rna.rds` |

`[MANUSCRIPT]` The Bonferroni family is stated as 13 genes × 8 cell types = 104. The
code's marker panel has **9** cell types and the candidate vectors vary between 6, 13
and 14 genes depending on the block. See item 4.

---

## 4. Candidate adjudication — the 8-dimension scoring matrix

`P2 开始` computes a score per gene across 8 dimensions (A–H), then ranks.

| Dimension | Meaning | How obtained |
|---|---|---|
| A `LODO_MASLD` | LODO stability | `in200` count: 6 → 2, 5 → 1.5, 4 → 1, else 0 |
| B `CAG` | gastric support | cohorts with `thr_ok & dir_ok`: 2 → 2, 1 → 1, 0 → 0 |
| C `Validation` | GSE135251 | direction must agree, then AUC_norm ≥ 0.9 → 2; ≥ 0.8 → 1.5; ≥ 0.7 → 1; ≥ 0.6 → 0.5; else 0.25 |
| D `CellSpecificity` | single-cell specificity | max over the 3 datasets; `pct_expr ≥ 0.25 & ratio ≥ 2` → 2; `pct ≥ 0.20` → 1.5; `pct ≥ 0.10` → 1; else 0.5 |
| E `CrossOrgan` | lineage conservation | liver two-dataset agreement + cross-organ agreement, four levels 2 / 1.5 / 1 / 0.5 |
| F `GeneticSupport` | MR direction support | count of `dir_ok` estimates: ≥5 → 2, 4 → 1.5, 3 → 1, 2 → 0.5, else 0; **default 0.5 when the gene has no MR row** |
| G `Novelty` | ⚠️ **hardcoded literal** | author-assigned values typed in the script, lines 105–107 |
| H `Verifiability` | ⚠️ **hardcoded literal** | author-assigned values typed in the script, lines 110–113 |
| Total | `rowSums(A:H)` | asserted equal to the stored `total` by `stopifnot` |

| Final classification | Genes |
|---|---|
| Retained (main) | **IL32, CDHR2, ANXA4, CADM2, LGALS3** (5) |
| Observation list | **RPS6KA1** (1) |
| Downgraded | **GOLM1, ANO10, SLC6A16, KIAA1958, FGA, FGB, LEPR, TSPAN3** (8) |
| Output | `p2_scoring_matrix.csv`, published as Table S11 |

`[CONFIRM]` The cut-points that separate "retained" from "observation" from "downgraded"
are **not** stated as a rule in the code — the script produces a ranked total and the
grouping is authorial. State the rule explicitly, because Table S11 is titled a
"transparent 8-dimension scoring matrix" and a reviewer will look for the decision
boundary. **Two of the eight dimensions (G, H) are hand-assigned**, and three total
values circulate in the files for IL32 and LGALS3 (see item 9).

---

## 5. External validation (held-out cohorts)

| Item | Value |
|---|---|
| Cohort | GSE135251 (liver, 216 samples: 206 disease / 10 control) |
| Procedure | **No re-screening.** The candidate genes were tested directly against the stored direction |
| Normalisation | genes with ≥ 5 non-zero samples → `log2(CPM + 1)` |
| Gene→Ensembl resolution | `org.Hs.eg.db` symbol lookup; if several Ensembl IDs map to one symbol, the **highest mean expression** row is kept |
| Direction concordance | `mean(case)` vs `mean(control)` compared with the stored expected direction |
| AUC | Mann–Whitney rank statistic — `AUC = (sum(ranks of cases) − n1(n1+1)/2) / (n1·n2)` |
| AUC normalisation | `AUC_norm = pmax(AUC, 1 − AUC)` — so a strongly **down**-regulated gene also scores near 1 |
| Association with severity | Spearman ρ against NAS score and against fibrosis stage, full 216 samples |
| AUC confidence interval | **2,000-resample percentile bootstrap** — `B <- 2000` |
| Spearman CI | Fisher z, `se = 1/sqrt(n − 3)` |
| Multiple testing | BH over the **5 pre-specified candidates** |
| Output | `shared_validation_auc.csv`, `GSE135251_spearman_CI.csv`, `GSE135251_AUC_bootstrap_CI.csv` |
| Secondary cohort | GSE153224 (gastric, 10 pooled libraries from 40 patients, 4 patients per pool, 5 vs 5 at library level) — **script not supplied**; results exist as `candidate_effects.csv` |

`[MANUSCRIPT]` Note that GSE135251's case/control split is 206 vs 10, so the AUC and
Wilcoxon statistics are auxiliary; the NAS-gradient Spearman correlation (n = 216) is
the primary validation metric. The supplementary NOTES state this explicitly. The
manuscript should be checked for the same framing.

Note on the AUC_norm convention: for a down-regulated gene such as CADM2 the raw AUC is
0.082. The `AUC_norm` transformation turns this into 0.918, which is why CADM2 scores
highly in dimension C of the scoring matrix. This is defensible for a *directional*
concordance metric but must be described as such, or a reviewer will read 0.918 as
discrimination performance.

---

## 6. Single-cell / single-nucleus analysis

| Item | Value |
|---|---|
| Datasets | GSE202379 (liver snRNA, 47 patients), GSE115469 (liver scRNA, cell-type confirmation), GSE134520 (stomach snRNA, 9 patients) |
| Aggregation | patient-level pseudobulk: summed counts per `Patient.ID × cell type` |
| Normalisation | CPM against the summed library size of that patient–cell-type unit: `pb / tot * 1e6` |
| **Minimum cells per patient–cell-type** | **≥ 20** |
| Localisation metric | mean expression per cell type + `pct_expr` (fraction of cells with expression > 0); `log1p` for count-based GSE134520 |
| Fibrosis stage per patient | **mode** of the per-cell `Fibrosis.score..F0.4.` annotation; the **rounded mean** was computed as a comparator and both were tested against three stored ρ values; the mode reproduced them and was adopted |
| Stage counts | F0–F4 = **5 / 9 / 12 / 12 / 9** patients |
| Locked fibrosis correlations | IL32 / Hepatocytes **ρ = +0.639**; CADM2 / Stellate **ρ = −0.761**; ANXA4 / Cholangiocytes **ρ = −0.585** |
| Patients per analysis pair | IL32/Hepatocytes **n = 47**; CADM2/Stellate **n = 37**; ANXA4/Cholangiocytes **n = 30** |
| Correlation test | Spearman, `exact = FALSE` (t asymptotic approximation) |
| GSE134520 grouping | case/control by lesion: NAG vs CAG, requiring ≥ 2 patients per side and ≥ 20 cells; Wilcoxon |
| GSE115469 lineage mapping | `Hepatocyte*`→Hepatocytes; Cholangiocytes; `Macrophage*`→Macrophages; `Hepatic_Stellate_Cells`→Stellate; `LSEC|Endothelial`; `T_Cells|NK-like`→Lymphocytes; `B_Cells|Plasma`→B_cells; `Erythroid`; else Other |
| GSE134520 cell-type annotation | marker-score `which.max` over 12 gastric lineages (G_epithelium, Pit, Neck, Parietal, Chief, Enteroendocrine, Fibroblast, Endothelial, Myeloid, T_NK, B_plasma, Mast) |
| Outputs | `p1_gse202379_localization.csv`, `p1_gse202379_pseudobulk_cpm.csv`, `p1_gse202379_fibrosis_cor.csv`, `p1_gse202379_fibrosis_stage_counts.csv`, `p1_gse115469_*`, `p1_gse134520_*` |

### The Bonferroni family — 13 genes × 8 cell types = 104

This is the family behind the manuscript's threshold, and **both dimensions are recoverable
from the code.** `R/17_bonferroni_family.R` derives them from the localisation tables and
asserts the product.

**The 8 cell types** — the unified liver lineage vocabulary the cross-dataset localisation
is built on (lineage map in the single-cell script, L129–151):

```
Hepatocytes · Cholangiocytes · Macrophages · Stellate ·
Endothelial · Lymphocytes · B_cells · Erythroid
```

The lineage map also emits an `"Other"` residual bucket for unmatched labels. **It is a
residual, not a cell type, and must not be counted** — including it would make the family
9 × 13 = 117 and change the threshold.

**The 13 genes** — the candidates surviving in the liver localisation matrices. GSE202379
carries 13 of the 14. **RPS6KA1 is the likely exclusion**, corroborated independently by
Supplementary NOTE 17 (*"RPS6KA1 is absent from the deconvolution reference gene rows"*)
and by the `excl_skip` mechanism in the deconvolution code, which exists precisely to
handle a candidate missing from the reference.

### ⚠️ Four cell-type counts appear in this study — none of them is wrong

| Analysis | Count | Vocabulary |
|---|---|---|
| Cross-dataset liver localisation (the 104-test family) | **8** | the unified lineage vocabulary above |
| Composition-adjustment marker panel (`R/06`) | **9** | Hepatocyte, Cholangiocyte, Kupffer_macro, T_cell, B_plasma, NK, Neutrophil, HSC_fibro, Endothelial |
| MuSiC reference (`R/09`) | **8** | the snRNA vocabulary; per Table S41 B-cell 1/2 merged, 99,687 cells |
| Gastric single-cell annotation (`R/08`) | **12** | gastric lineages including Pit, Neck, Chief, Parietal |

Each is correct for its own analysis. The manuscript should name the vocabulary wherever
the number appears, rather than letting a bare "8 cell types" stand unqualified for the
whole paper.

---

## 7. Human Protein Atlas annotation check

| Item | Value |
|---|---|
| **Version cited in the manuscript** | **25.0** — see the timeline below; the release in force on the access date was 25.1 |
| **Version actually current on the access date** | **25.1**, released **2026-05-25** |
| **Accessed** | **2026-09-14** — author-confirmed. The manuscript says 2026-09-05; the supplementary NOTES say 2026-09-13 |
| Genes | IL32, ANXA4, CDHR2, CADM2, LGALS3 (the 5 retained candidates) |
| Method | **Manual curation, hand-typed into a `data.frame`.** No programmatic HPA query. |
| Antibodies recorded | IL32 HPA029397, CAB030029 · ANXA4 HPA007393, CAB005076, CAB017560 · CDHR2 HPA012569, HPA017053 · CADM2 HPA010024 · LGALS3 HPA003162, CAB005191 — **10 identifiers** |
| Fields per gene | liver level, liver cell types, stomach level, stomach cell types, tissue specificity, antibody IDs, reliability grade, author verdict — **9 columns** |
| Output | `p4_HPA_protein_check.csv`, published as Table S27 |
| Reliability vocabulary | HPA grades: Enhanced / Supported / Approved / Uncertain / Enhanced (orthogonal) |

### ⚠️ Version versus access date — the citation does not match the date

Checked against HPA's own release history:

| Release | Date | Source |
|---|---|---|
| version 25 | 2025-11-11 | announced at HUPO, Toronto |
| **version 25.1** | **2026-05-25** | HPA release history — *"Twenty-fifth major release"* |
| current release at the time of checking | 25.1 | HPA home page |

The annotations were made on **2026-09-14**. By then 25.1 had been live for nearly four
months, so whatever was on screen during curation was **25.1** — the manuscript's
"version 25.0" does not name it.

**Two readings, needing different fixes:**

- If "25.0" is shorthand for the version 25 *family* — 25.0 and 25.1 share the major
  number — the intent is defensible but the citation is imprecise. Write **25.1**.
- If "25.0" was meant literally as `v25.proteinatlas.org`, the citation points at a
  release superseded before the access date, and the five genes need re-checking against
  25.1.

HPA hosts older releases at `vX.proteinatlas.org`, so the two can be compared directly.
**That comparison is the decisive check:** if the five genes' tissue levels, cell-type
attributions and antibody reliability grades are unchanged between 25.0 and 25.1, the
citation is corrected to 25.1 with no analysis change needed. `R/16_verify_hpa_version.R`
documents the procedure and prints what the site currently reports.

### What the curator table does and does not record — verified

The **generated artifact** was obtained and inspected
(`docs/provenance/p4_HPA_protein_check_as_generated.csv`, sha256 `a6c27264…`). It has
**exactly 9 columns**, matching the eight HPA fields plus the author verdict. It contains:

- **no version number** — the string `25.0` does not appear
- **no access date** — none of the three candidate dates appears
- **no year at all** — the string `2026` does not appear anywhere in the file

So the manuscript's version and any access date come from outside this artifact. HPA
versions are normally tracked outside an analysis output, so this is not an error in
itself — but the citation is uncorroborated by the artifact, and one element of it is now
known to be inconsistent with the confirmed access date.

The reconstruction in `R/12_HPA_check.R` reproduces the artifact's 9 columns exactly, then
adds four columns — `hpa_version`, `access_date`, `provenance` and
`flag_stomach_contradiction` — so the provenance of every value is explicit.

### One point in the table's favour

`"Glandular cells"` is **not** an ad-hoc string. It is a **standard HPA cell-type group
value**, appearing verbatim in HPA's own search vocabulary alongside groups such as
hepatocytes and smooth muscle cells. Using it as a cell-type attribution is HPA-consistent
usage, not a typo.

### The remaining problem is narrower

`"Glandular cells"` is applied uniformly to all five genes, including CADM2, whose own
`stomach_level` reads `"Not detected (raw Low)"`. A gene cannot be attributed to a cell
type while its level reads not detected. The value is **not changed** here — silently
editing curated data would alter a published table — but the script flags the row at run
time. See `manuscript_discrepancies.md` item 8.

---

## 8. MuSiC deconvolution (sensitivity analysis)

| Item | Value |
|---|---|
| Package | `MuSiC` **1.0.0**, installed from `github.com/xuranw/MuSiC` (no longer on Bioconductor) |
| Interface | `music_prop(bulk.mtx = …, sc.sce = …, clusters = "ct_merged", samples = "Patient.ID")` with a `SingleCellExperiment` |
| Reference | `B_reference_merged.rds`; signature = **mean over cells per `Patient.ID × ct_merged`** combination |
| Reference dataset | **GSE202379**, liver snRNA |
| Cell types | **8** (per the script's own prose); the true names come from `sort(unique(meta$ct_merged))` |
| **Marker panel construction** | specificity score `spec = theta / (rowMeans(theta) + 1e-12)`, where `theta` is the column-normalised reference signature; take the **top 50 per cell type**; `unique()` across cell types |
| **Anti-circularity** | the 6 frozen candidates **and** `^MT-`/`^MTRNR` genes are zeroed in `spec` before ranking, then `stopifnot(!any(CAND %in% markers))` |
| Panel size | **≤ 400** (top 50 × 8 cell types, deduplicated). "400" is prose in the table title; the true count is never printed or stored |
| Candidate exclusion list | `CDHR2, ANXA4, CADM2, IL32, LGALS3, RPS6KA1` (hand-typed literal) |
| Bulk cohorts | **4**: GSE126848 (57), GSE130970 (78), GSE162694 (143), **GSE174478 (94)** |
| Full-gene run (failed) | `markers = all shared genes` |
| Marker-panel run (accepted) | `markers = intersect(bulk, reference, markers)` |
| Centered control | `centered = TRUE`, full-gene, GSE126848 only |
| Candidates examined | **3**: CADM2, ANXA4, LGALS3 |
| Adjusted model | per cohort, `limma` on the 3 candidates with endpoint + selected MuSiC proportions as covariates |
| Covariate selection | `union(c("Hepatocytes","Stellate"), names(mu)[mu >= 0.02])`; zero-variance covariates dropped at `var > 1e-6` |
| Endpoint | GSE126848/GSE130970/GSE162694: case/control 0-1 · **GSE174478: numeric `fibrosis_stage` 0–4** |
| Attenuation | `1 − arm3_logFC_adj / arm3_logFC_raw` |
| Outputs | `B_music_props.csv`, `B_music_props_list.rds`, `B_music_qc.csv`, `B_music_marker_props_list.rds`, `B_b4_compare.csv`, `B_b4_pheno.csv` → Tables S32–S36 |

### Pre-registered QC thresholds

| Check | Threshold | Result |
|---|---|---|
| Bulk data are raw counts | `max(B) < 50` on any cohort fails the gate | passed |
| Proportion row sums | `max(abs(rowSums(P) − 1)) < 1e-6` | passed |
| Non-negative proportions | `P > -1e-8` | passed |
| Reconstruction | Pearson / Spearman on `Y ≈ Σ p_ct · S_ct · θ_ct` | **full-gene failed** |
| Proportion collapse | mean hepatocyte proportion | **full-gene 0.996 on GSE126848** |

### Recorded negative results (all reported in the manuscript's supplementary tables)

1. **Full-gene MuSiC failed pre-registered QC** — reference/vs/bulk gene-composition
   mismatch. Reference hepatocyte `theta_ALB ≈ 0.0044` vs bulk ALB fraction 0.16–0.27,
   a ~40–60× discrepancy, driving proportion collapse (GSE126848 hepatocytes 0.996).
2. **A second, independent mismatch axis**: `MALAT1`/`NEAT1` are top genes in every
   reference cell type (5–15%) but absent from the bulk top-10 — nuclear capture vs polyA.
3. **The centered control also degenerated**: `centered = TRUE` gave
   `r.squared.full = 0.014` with B-cells 0.667 / Lymphocytes 0.300.
4. **Plain NNLS was degenerate** — solutions did not sum to 1 without renormalisation;
   fraction-space R² medians 0.092 / 0.062 / 0.037 / 0.030.
5. **CIBERSORTx was never run** — installation infeasible in the analysis environment.
6. **RPS6KA1 is absent from the deconvolution reference gene rows**, so no reference-side
   adjustment exists for it.
7. **GSE174478 is not part of the P0 training set**, and its endpoint is a numeric stage
   rather than case/control, so comparison arms 1 and 2 are `NA` for it.

---

## 9. Exploratory Mendelian randomisation

**The reported version is the third run** (`eQTLGen_eqtl 主分析第三次跑`), which writes
`mr_all_results_v4.csv` — the file every downstream script reads. The first two runs are
superseded and should not be presented as the analysis.

| Item | Value |
|---|---|
| Package | `TwoSampleMR` **0.7.9**, with `data.table`, `ieugwasr`, `coloc` |
| Exposure | 10 genes with cis-eQTL instruments: IL32, GOLM1, TSPAN3, ANXA4, RPS6KA1, ANO10, SLC6A16, KIAA1958 (r² = 0.001 clumped set), LGALS3, CADM2 (**looser purchased set, 1 SNP each**) |
| **Instrument p threshold** | **`clump_p = 5e-8`** — applied during clumping, not re-applied in the MR scripts |
| **LD clumping** | **`ieugwasr::ld_clump(..., clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8, pop = "EUR")`** — 10,000 kb window, r² < 0.001, European reference panel |
| **Clumping independence check** | `ieugwasr::ld_matrix(snps, pop = "EUR", with_alleles = TRUE)`, then `max(r²)` over the upper triangle; asserted ≤ 0.001 |
| **The 0.01 source directory was purchased** | `eqtlclump/` is described in the original as "买来的 0.01 文件所在目录" (the directory holding the purchased r² = 0.01 files); the r² = 0.001 set was produced from it by re-clumping (`clump自查`) |
| Harmonisation | `harmonise_data(..., action = 2)` — palindromic SNPs **dropped**, not strand-inferred |
| **F statistic, as recorded during clumping** | **`F = (beta / se)²`** per SNP, univariate; median and minimum reported per gene |
| **F statistic, as implemented in the MR run** | `R2 = 2·β²·EAF·(1−EAF) / (2·β²·EAF·(1−EAF) + 2·N·EAF·(1−EAF)·SE²)`, `F = R2·(N−2)/(1−R2)`; per-SNP `F > 10` filter |
| Steiger filtering | **removed** in the v3 revision and absent from v4 — `steiger_ok` is retained as a column but is `NA` throughout |
| Methods | single instrument: **Wald ratio** (manual formula). Multiple instruments: **IVW** (`mr_ivw`) alongside Egger, weighted median and weighted mode |
| Heterogeneity | `mr_heterogeneity` → Cochran's Q, IVW row reported |
| Pleiotropy | `mr_pleiotropy_test` → Egger intercept p |
| MR-PRESSO | `run_mr_presso(harm, NbDistribution = 1000)` → global p |
| Multiple testing | BH **within tier**: {main + sens} pooled as one family, {explor} as a separate family |
| Estimates | **69** |
| Result | **No estimate survived FDR** |
| Follow-up (not in the MR scripts) | `P3 脚本（CADM2 × LiverPDFF × GTEx Liver）` runs `coloc` on CADM2 × LiverPDFF against GTEx liver eQTLs from the eQTL Catalogue |
| Output | `mr_all_results_v4.csv`, `harmonised_for_coloc.rds` |

### ⚠️ Two F-statistic definitions are in play

The clumping self-check records instrument strength as **`F = (beta/se)²`**, the simple
univariate form. The MR run itself computes F from **R²**, using

```
R2 = 2·β²·EAF·(1−EAF) / (2·β²·EAF·(1−EAF) + 2·N·EAF·(1−EAF)·SE²)
F  = R2·(N−2) / (1−R2)
```

These give different numbers for the same SNP — the second is sample-size aware, the first
is not. `mean_F` in `mr_all_results_v4.csv` comes from the second. If the manuscript quotes
a median or minimum F, it is worth checking which definition produced it.

### ⚠️ TSPAN3: an instrument set that was repaired by hand

The clumping self-check records a real problem and its fix:

1. Re-clumping TSPAN3 at r² < 0.001 returned **0 SNPs**. The script contains an explicit
   guard for this — `if (nrow(kept) == 0) stop("clump 结果为空，停下检查")` — and the
   file was left empty.
2. Re-run with the empty-file gate `stopifnot(nrow(d) >= 2)`.
3. The repair: the remaining SNPs still exceeded r² = 0.001 in one pair, so the
   **lower-F SNP of that pair was dropped by hand** (`drop <- pair_bare[which.min(fstat[pair_bare])]`),
   with `stopifnot(nrow(kept2) == nrow(d) - 1)` guarding the edit.
4. Final state: 4 SNPs, max r² ≤ 0.001.

The supplementary NOTES already record the consequence — *"ukbgast TSPAN3: 2 SNPs"* — so
the reduced instrument count is documented. But note that TSPAN3 is one of the eight
**downgraded** genes, and dimension F of the scoring matrix scores it 0.5 on a single
agreeing estimate. The instrument repair is worth mentioning in the README rather than
leaving a reviewer to infer it from a 4-SNP file.

### LGALS3 and CADM2 use a different instrument source

Both are read from the **looser purchased set** (`eqtlclump/`, r² = 0.01), not the
r² = 0.001 set, and each contributes **one SNP**. That means both pairs are Wald ratios
with no Egger intercept, no Q statistic and no MR-PRESSO. Since CADM2 is one of the five
retained genes, this should be stated wherever CADM2's MR-support is discussed.

---

## 10. Multiple-testing families

| Analysis | Family | Correction | Threshold |
|---|---|---|---|
| Fig 1 effect-size heatmap | **140** gene–cohort tests (14 candidates × 10 cohorts) | BH, within the heatmap | — |
| Discovery RRA | 200 up + 200 dn carried forward | RRA score threshold only | Score < 0.01 |
| LODO audit | 14 genes × 6 drop-iterations = **84 rank comparisons**, plus 14 × 6 direction tallies | none — descriptive | rank ≤ 200 or P < 0.01 |
| External validation | **5** pre-specified candidates per cohort | BH | — |
| Single-cell pseudobulk | **104** tests (13 genes × 8 cell types) | **Bonferroni** | **P = 4.8 × 10⁻⁴** |
| Composition adjustment | 6 candidates × up to 9 cell types | none — sensitivity | — |
| MR | main+sens one family; explor a second | BH, per family | — |
| AUC bootstrap CI | 5 genes × 2000 resamples | percentile CI, not a test | — |

Arithmetic check: 0.05 / 104 = 4.8077 × 10⁻⁴ → **4.8 × 10⁻⁴** ✓ consistent with the
manuscript.

Note there is **no permutation family**. The pipeline performs no resampling-based
testing, so there is no multiplicity correction across resamples. See §3.2.

---

## 11. Software environment

| Tool | Version |
|---|---|
| R | **4.6.0** (2026-04-24 ucrt) |
| Platform | x86_64-w64-mingw32/x64, Windows 11 x64 |
| `limma` | 3.68.5 |
| `Seurat` | 5.5.1 |
| `pheatmap` | 1.0.13 |
| `TwoSampleMR` | 0.7.9 |
| `MuSiC` | 1.0.0 (GitHub) |
| `RobustRankAggreg` | 1.2.1 |
| `metafor` | 5.0-1 |
| `edgeR` | 4.10.5 |
| `ComplexHeatmap` | 2.28.0 |
| `SingleCellExperiment` | 1.34.0 |
| `Biobase` | 2.72.0 |
| `Matrix` | 1.7-6 |
| `nnls` | 1.6 |
| `org.Hs.eg.db` | present |
| Annotation packages | `hugene11sttranscriptcluster.db`, `hugene20sttranscriptcluster.db`, `illuminaHumanWGDASLv4.db`, `hugene10sttranscriptcluster.db`, `illuminaHumanv3.db` |
| Locale | Chinese (Simplified)_China.utf8 |
| Time zone | Asia/Shanghai |

Full list: `sessionInfo.txt`. Package inventory at build time: `docs/installed_packages.csv`.

Other packages used for parts of the pipeline but not installed in this environment at
capture time: `GEOquery`, `sva`, `DESeq2`, `clusterProfiler`, `WGCNA`, `preprocessCore`,
`biomaRt`, `randomForest`, `coloc`, `ieugwasr`, `openxlsx`, `quadprog`, `TOAST`,
`MCMCpack`, `SingleR`, `celldex`.

---

## 12. Script execution order

```
00_setup.R                      paths, parameters, seeds
01_download_and_DE.R            GEO download; limma per cohort; RRA discovery screen
02_cross_organ_candidates.R     gastric DEGs; concordance; candidate_genes_final.csv
03_assets.R                     expression matrices, group vectors, RRA rank columns
04_audit_LODO_meta.R            LODO exact RRA; CAG support; effect sizes; meta
05_composition_adjustment.R     marker-score composition adjustment; batch diagnostics
06_validation_GSE135251.R       held-out liver validation; AUC; bootstrap CI
07_singlecell_pseudobulk.R      GSE202379 / GSE115469 / GSE134520 handling
08_deconvolution_MuSiC.R        bulk gate; full-gene MuSiC; marker panel; four-arm table
09_candidate_scoring.R          8-dimension adjudication
10_MR_exploratory.R             two-sample MR (v4 configuration)
11_HPA_check.R                  protein-level annotation table
12_figures.R                    Figures 1–5, S1
```

Dependency notes: `03` must run before `04`–`09`; `04` before `09` (arm 1 and arm 2);
`07` before `12` (Fig 3, Fig 5, Fig S1); `08` before `09` (dimension F uses the MR table,
which needs `10`; run `10` before `09` if regenerating the scoring matrix).

---

## 13. Reconciliation items — code vs manuscript

Full detail in `docs/manuscript_discrepancies.md`. Summary:

1. **RRA background size.** Manuscript: 20,000. Discovery screen: `length(all_sym)`, a
   runtime union. 20,000 belongs to the *audit* recomputation and was hand-selected.
   → Fix the manuscript wording or document both numbers.
2. **"20,000 permutations" is a misnomer.** The manuscript Methods and the Table S02
   title both describe the LODO audit as a permutation test with 20,000 replicates. There
   is no permutation step: the 20,000 is the RRA background size `N` and the procedure is
   closed-form, run once per dropped cohort. **The analysis is correct and fully
   reproducible; only the label is wrong.** Fix the Methods wording, the Table S02 title,
   and the Code availability statement. See §3.2.
3. **Bonferroni family.** Manuscript: 13 genes × 8 cell types = 104. Code: the marker
   panel has 9 cell types; candidate vectors are 6, 13 or 14 depending on the block.
   → State the exact 13 genes and 8 cell types in the README.
4. **Cell-type count — RESOLVED.** Four counts appear (8 / 9 / 8 / 12) and each is
   legitimate for its own analysis. The one behind the Bonferroni family is the 8-lineage
   cross-dataset liver localisation, named in section 6. The manuscript should name the
   vocabulary wherever the number appears.
5. **HPA access date.** Manuscript 2026-09-05; code NOTES 2026-09-13.
6. **HPA version.** 25.0 appears in the manuscript only, never in code.
7. **HPA table provenance.** Hand-typed; one cell self-contradictory; no per-gene
   retrieval record.
8. **Scoring matrix.** Dimensions G and H are hand-assigned; the retained/observation/
   downgraded boundary is not a coded rule; three different total values for IL32 and
   LGALS3 appear across the files.
9. **GSE174478.** Used as the fourth deconvolution cohort but **absent from the
   manuscript's Data availability list** of 13 accessions. See `docs/manuscript_discrepancies.md`.
10. **13 vs 14 genes.** `genes14 <- cc$genes_present # 13 个` in the single-cell script.
11. **Fig S2.** Not referenced by any script; the deconvolution output is tables only.
12. **Two F-statistic definitions.** The clumping record uses `(beta/se)²`; the MR run
    uses the sample-size-aware R² form. See §9. Whichever the manuscript quotes, name it.
13. **TSPAN3's instrument set was repaired by hand.** Re-clumping returned zero SNPs and
    one over-threshold LD pair was resolved by dropping the lower-F SNP. The reduced
    instrument count is already recorded in the supplementary notes, but the repair
    itself is only visible in the clumping self-check.
14. **LGALS3 and CADM2 use a single instrument each**, from the looser r² = 0.01 set.
    All their estimates are Wald ratios with no pleiotropy diagnostics. CADM2 is a
    retained gene, so this belongs next to any MR-based support for it.

---

## Resolved during repository build

- **"Permutation testing" — it does not exist.** The pipeline performs no Monte Carlo
  step. The LODO audit is deterministic and closed-form, and is fully reproducible from
  the code here. The label originated in the Methods text and the Table S02 title. See
  §3.2 and `manuscript_discrepancies.md` item 2.
- **MR clumping parameters** — supplied by the clumping self-check after the first
  inventory: `clump_kb = 10000`, `clump_r2 = 0.001`, `clump_p = 5e-8`, `pop = "EUR"`.
  See §9.
- **Marker-restricted MuSiC** — present in the bulk-gate file; identified on the second
  read. See `code_inventory_raw.md` §3.2.
- **HPA version and access date** — the generated table was obtained and searched: it
  records neither, and contains no year at all. The manuscript's version and both
  candidate dates come from outside the artifact. Not an error, but it needs confirming
  because the manuscript and the supplementary notes disagree with each other. See §7.

# Manuscript-relevant discrepancies found in the code

These are things the code says that the manuscript does not, or where two places in the
code disagree. Each one is a question a reviewer could ask. None of them are hidden —
they are all answerable — but they must be answered deliberately rather than discovered.

Status key: **OPEN** = needs author confirmation before submission ·
**RESOLVED** = code settles it · **DOC** = must be documented in README, not an error.

---

## 1. "20,000 permutations" — a misnomer. There is no permutation analysis — OPEN

**Corrected after a targeted re-audit.** An earlier version of this item described a missing
permutation script. That was wrong: **the analysis contains no Monte Carlo step at all**, so
nothing is missing.

The manuscript Methods describes **"20,000 rank-permutation replicates"**, and Supplementary
Table S02 is titled *"Leave-one-dataset-out exact RRA (20,000 permutations)"*. Both describe
the same procedure — and **that procedure performs no permutations.**

### What the code actually does

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

```
line 45   cal_200 <- rra_rho(r_up[i, rnk_cols], 20000)     # background = 20,000
line 51   use_N20000 <- TRUE
line 52   den_lodo <- function(k) rep(20000, k)
line 67   full_sc <- apply(mat, 1, rra_rho, denom = den_lodo(ncol(mat)))
line 72   sc <- apply(mat[, -j], 1, rra_rho, denom = den_lodo(ncol(mat) - 1))
line 76   out[[paste0("in200_drop_", ...)]] <- rk[genes] <= 200
line 77   out[[paste0("lt001_drop_", ...)]] <- sc[genes] < 0.01
```

`20000` is the **background gene count `N`**, passed as the denominator of a closed-form
Beta-distribution p-value. The audit drops one cohort at a time and recomputes —
**six iterations, not 20,000.** Every output is a deterministic function of the input ranks.

### Full-tree search for a resampling step

Across all original files, searching `set.seed`, `permut`, `置换`, `sample(`, `replicate(`,
`runif`, `rnorm`, `rmultinom`:

| Location | What it is |
|---|---|
| `第一次预审稿文件处理.R` L206, L352 | `set.seed(20260914)` |
| `第一次预审稿文件处理.R` L255–257 | `replicate(B, { b1 <- sample(i1, …); b0 <- sample(i0, …) })` with `B <- 2000` — the **AUC bootstrap** (Table S39) |
| `第一次预审稿文件处理.R` L408–409 | the same bootstrap, repeated for a second cohort |
| three MR scripts, e.g. `eQTLGen… 第三次跑` L126 | `run_mr_presso(..., NbDistribution = 1000)` — **MR-PRESSO** |
| `第一次预审稿文件处理.R` L442 | the title string itself |

**No permutation loop exists, because none was needed.** `set.seed(20260914)` belongs to the
AUC bootstrap; it was never a permutation seed.

### Why this is good news, not a gap

1. **The LODO table is fully reproducible from this repository.** A closed-form procedure
   needs no seed. Quoting a seed would have been the weaker position.
2. **No analysis is missing.** The earlier "script could not be recovered" framing was
   wrong; `repository_scope.md` §2, `code_inventory_raw.md` §3.1 and the README have all
   been corrected.
3. The defect is **wording**, not results. No statistic changes.

### But the wording must be fixed, in three places

| Where | Current | Should be |
|---|---|---|
| **Methods** | "20,000 rank-permutation replicates" | the exact robust rank aggregation, recomputed with each discovery cohort removed in turn, at a background of 20,000 genes |
| **Table S02 title** | "Leave-one-dataset-out exact RRA (20,000 permutations)" | "… (background N = 20,000)" |
| **Code availability** | "…and the random seeds used for permutation testing." | remove the clause — no permutation testing is performed. Replace with parameters that do matter: RRA background size, closed-form P-value implementation, minimum-cell pseudobulk threshold, Bonferroni threshold — all documented in this repository |

**Why this matters more than it looks.** A reviewer who reads "permutation testing", then
opens the repository and finds no permutation step, must choose between two readings: the
code is incomplete, or the manuscript describes the method wrongly. The first is a
rejection risk. Leaving the phrase in converts a labelling slip into an apparent
reproducibility failure — and the fix costs one sentence.

## 3. RRA background size — OPEN

The manuscript states the robust rank aggregation background size as **20,000 genes**.

The discovery screen does not use 20,000:

```r
# 第一步数据下载（MASLD侧 差异分析已完成）.R, line 613
N <- length(all_sym)   # 用六队列 symbol 并集当背景，别用拍脑袋的 18000
```

`N` is the runtime size of the union of symbols across the six liver cohorts. A later
audit block (`P0组成校正`, lines 400–426) compared three definitions —
`common = 11982`, per-cohort sizes, and `20000` — and then 20,000 was **selected by
hand** for the LODO recomputation:

```r
# 磁盘层执行脚本2026.9.12 16 56.R, line 51
use_N20000 <- TRUE   # 校准后手工指定
```

→ Confirm which value the manuscript is describing. If the recovered 14 candidates came
from the discovery run, the background was the runtime union, and the manuscript's
"20,000" describes the audit re-computation only. The README should state both numbers
and say which analysis each belongs to.

## 4. Cell-type count: 8 or 9 — OPEN

The Bonferroni threshold is stated as **P = 4.8 × 10⁻⁴**, derived from 13 genes ×
**8 cell types** = 104 tests. Arithmetic checks out: 0.05 / 104 = 4.8077 × 10⁻⁴. ✓

But the code uses different numbers of cell types in different places:

| Where | Count | List |
|---|---|---|
| `P0组成校正` `mk_liver` marker list | **9** | Hepatocyte, Cholangiocyte, Kupffer_macro, T_cell, B_plasma, NK, Neutrophil, HSC_fibro, Endothelial |
| `P2 开始` / `P0组成校正` candidate vector | **6** | IL32, CDHR2, LGALS3, ANXA4, CADM2, RPS6KA1 |
| `磁盘层执行脚本` `cag_support`, `eff_tab`, `meta_tab`, `audit` | **14** | all candidates |
| `正交去卷积…` MuSiC reference (`ct_merged`) | reported **8** | per Table S41 note, B-cell 1/2 merged, 99,687 cells |
| `P2 开始` `genes14` | **14** | |

→ Confirm the exact 8 cell types and the exact 13 genes behind the 104 tests. The
README must state them explicitly, because "13 × 8" is a checkable claim and the code
does not contain a single place where both numbers appear together.

## 5. Fibrosis pseudobulk threshold — OPEN

The manuscript describes a "pre-specified pseudobulk threshold" of ≥20 cells per
patient and cell type.

In the code the 20-cell rule appears as an **inline literal**, in several variants:

```r
sel <- pb_df$n_cells >= 20                     # 跑 compact 版三件套, line 48
sub <- pb_df[pb_df$lesion %in% c("NAG","CAG") & pb_df$n_cells >= 20, ]   # line 377
okF <- st & !is.na(meta$Fibrosis.score..F0.4.)  # 第一次预审稿文件处理, line 181
ok  <- ncell$Stellate >= 20                     # 第一次预审稿文件处理, line 744
stage_counts <- function(gene, ct, pb, ncell, min_cells = 20)   # line 790
```

→ It is a consistent value, but it is retyped rather than declared once. The repository
centralises it in `R/00_setup.R` as `MIN_CELLS_PER_PATIENT_CELLTYPE <- 20L`. Worth
noting in the README that the value is uniform across analyses.

## 6. 13 vs 14 genes in the single-cell localisation — DOC

`跑 compact 版三件套` line 8 reads `genes14 <- cc$genes_present # 13 个` — the variable
is named `genes14` but holds **13** genes. The comment flags it, so this appears to be
known. Later `genes14` vectors in the same file and in `P2 开始` do contain 14.

→ Probably one candidate (likely RPS6KA1) was absent from the GSE202379 compact object.
The README should state which gene was missing and from which dataset, rather than
leaving a variable whose name contradicts its contents.

## 7. Human Protein Atlas: version and access date — OPEN, now with evidence

The version and access date claimed in the manuscript appear **nowhere in the generated
artifact**. The original `p4_HPA_protein_check.csv` was obtained and searched
(kept at `docs/provenance/p4_HPA_protein_check_as_generated.csv`, sha256 `a6c27264…`):

| Search | Result |
|---|---|
| `version` | absent |
| `25.0` | absent |
| `2026-09-05` | absent |
| `2026-09-13` | absent |
| `2026` | **absent** |

The file has nine columns — the eight HPA fields plus the author verdict — and records
neither a version nor a retrieval date nor even a year.

Two conflicting dates are nonetheless recorded elsewhere:

| Source | Date |
|---|---|
| Manuscript Data availability text | **5 September 2026** |
| `落盘 最终核对表` NOTES sheet (both variants) | **2026-09-13** |

→ Confirm the HPA release number and pick one access date. This is not an error in the
manuscript — version numbers are normally tracked outside the output file — but the claim
is currently unverifiable from the repository, and the two dates contradict each other.

## 8. HPA protein table is hand-entered — OPEN

`落盘 最终核对表` builds `p4_HPA_protein_check.csv` as a `data.frame` literal: every
protein level, cell-type string, antibody ID and reliability grade is typed by hand.
Nothing is queried from HPA programmatically.

Two specific weaknesses in that table:

- `stomach_cells = rep("Glandular cells", 5)` applies the same cell-type string to all
  five genes, **including CADM2**, whose own `stomach_level` is `"Not detected (raw
  Low)"`. Those two cells contradict each other.
- The `liver_level` and `liver_cells` columns mix granularities and evidence modalities:
  `"Low"` next to `"Cholangiocytes Low; Hepatocytes Low"`; and strings like
  `"Not detected (bulk IHC)"` and `"myeloid cells by scRNA"` embed the method inside a
  cell that is otherwise method-free.

→ Either re-derive the table from HPA with a script, or annotate each cell with its
provenance. Reviewers accept curated annotation tables; they do not accept tables that
contradict themselves.

## 9. Scoring matrix: two hard-coded dimensions — OPEN

The 8-dimension scoring matrix (`P2 开始`, `p2_scoring_matrix.csv`) computes dimensions
A–F from data but **G (`G_Novelty`) and H (`H_Verifiability`) are literals typed by the
author** (lines 105–113).

Separately, three different totals for IL32 and LGALS3 circulate in the files:

| Value | Source |
|---|---|
| IL32 total 13.5, LGALS3 total 12.25 | quoted as *superseded draft* text inside generated notes |
| IL32 total 13.0, LGALS3 total 11.75 | hand-typed note 12, and stated as the final values |
| computed `p2$total` | asserted equal to the A–H sum by `stopifnot` |

→ The final values are the ones to publish, but the repository should carry a note
explaining that G and H are author-assigned expert scores, not derived quantities. A
reviewer seeing a "transparent 8-dimension scoring matrix" (Table S11) will otherwise
assume all eight dimensions are computed.

## 10. Contradiction between two Supplementary NOTES variants — OPEN

`落盘 最终核对表` contains two successive versions of the workbook's NOTES sheet that
say **opposite** things about the single not-estimable MR pair:

- Block 2 variant, note 9: *"empty/not-estimable pairs recorded as NA rows"*
- Block 3 variant, note 9: *"The single not-estimable pair (ANXA4 x Gastritis_UKB) is
  absent from the harmonised list"*

Both describe `TableS24`. Whichever block ran last decides the published text.

→ Confirm the correct statement and delete the stale variant. This is a live
inconsistency inside a submitted supplementary file, not merely in the code.

## 11. Audit checks that do not assert — OPEN

`结果部分确认2026.9.14 19 23.R` performs **no** numeric verification: it prints values
for the author to eyeball. It contains no `stopifnot`. One of its checks has a stated
expected value (chip matrix dimensions "should be 54/63/148") that is printed but never
asserted, and its GEO sample count is computed as
`length(grep("^!Sample_title", txt))` — a **line** count, not a sample count, so the
numbers it prints for GSE48452/GSE83452/GSE89632 are wrong by construction.

Also, the same file selects rows by `d$gene` in two checks and by `d$X` in a third, on
the same CSV — one of the two must silently return nothing.

→ This file is audit scaffolding, not a deliverable. It should be excluded from the
published repository, or replaced by a script with real assertions. The pass/fail checks
that matter are already asserted in `磁盘层执行脚本` and `第一次预审稿文件处理`.

## 12. Third-party code provenance — OPEN, BLOCKING

See `docs/code_inventory_raw.md` §4. The file that performs the liver-side download,
limma differential expression and RRA discovery carries a vendor notice threatening
retraction for unauthorised use. Must be resolved by the author before publication.

**Author decision (2026-09-28):** treat the repository as the author's own code only.
The affected script is held out of the public repository, and the scope limitation is
stated in the README. See `docs/repository_scope.md`.

## 13. `B_b4_compare.csv` was edited outside R — OPEN, CRITICAL

This is the single biggest reproducibility hazard found in the deconvolution code.

The block that builds Table S32 selects columns that the generating code **never
creates**:

```r
# 第 2 步 v4 …, line 1249–1250
nm <- c(..., verdict = "verdict", note = "note", ...)
s32 <- b4[, names(nm)]
```

But the data frame written at lines 1133–1145 — and saved to `B_b4_compare.csv` at
line 1155 — has **no `verdict` column and no `note` column**. Unless the CSV was edited
by hand in a spreadsheet between those two points, `b4[, names(nm)]` errors out with
"undefined columns selected".

Corroborating evidence that it *was* hand-edited: line 1194–1195 prints the column names
of the re-read CSV under the comment `# ---- B_b4_compare 当前列名（S32 取材）----`, and
NOTE 18 asserts that GSE126848 "is flagged deconvolution-sensitive in Table S32" — a flag
that only exists if `verdict` was added manually.

→ The published Table S32 therefore contains at least one column that is **not**
machine-generated from the pipeline. Either add `verdict`/`note` to the generating code
so the CSV is reproducible, or document explicitly in the README that these two columns
are curator annotations appended after the analysis.

## 14. Table S35 rows are transcribed from a console log — OPEN

Several rows of the deconvolution QC table are **typed into the script as literals**,
with the source labelled "recorded from run log":

```r
# lines 1297–1299
add("Full-gene NNLS", q, "R2 in fraction space (median)",
    c(GSE126848 = 0.092, GSE130970 = 0.062, GSE162694 = 0.037, GSE174478 = 0.030)[q],
    "recorded from run log (diagnostic D3)")

# lines 1312–1314
add("Marker-panel MuSiC (fallback)", q, "r.squared.full (median)",
    c(GSE126848 = 0.027, GSE130970 = 0.348, GSE162694 = 0.241, GSE174478 = 0.172)[q],
    "recorded from run log (E1)")

# lines 1316–1317
add("Centered control", "GSE126848", "r.squared.full", 0.014,
    "…proportions degenerated (B-cells 0.667, Lymphocytes 0.300); recorded from run log")
```

These values never touch a computed object, so **re-running the analysis would produce
identical output even if the numbers had changed.** The author's own comment at line 1270
acknowledges the split: `# ---- S35：QC 记录（可算的现算，须重跑 music_prop 的用 run-log 值并标注）`.

The same applies to the whole NOTES 13–18 block (lines 1347–1353), which carries
quantitative claims — `theta_ALB ~0.0044`, bulk ALB `0.16-0.27`, `~40-60-fold`, hepatocyte
proportion `0.996`, `400` markers, `Pearson 0.18, Spearman 0.48` — as hand-written prose.

→ Table S35's *mechanism* (ALB discordance, proportion collapse) **is** computed in-script
at lines 1277–1285; only the summary numbers are transcribed. The README should mark which
Table S35 rows are computed and which are transcribed, so a reviewer reading the QC record
knows what re-running would verify.

## 15. Marker panel size "400" is unverifiable — OPEN

Table S34's title and NOTE 15 both state **400 specificity-selected genes**. The code
builds the panel as `unique(unlist(top-50-per-cell-type))` (line 344), so duplicates
across cell types are collapsed and the true size is **≤ 400**. The value is never
printed, stored or asserted — only written as prose.

The code also contains an unexplained literal `"108 集"` in a `cat()` string at line 331,
whose meaning is not derivable from the file and appears to be a stale number from an
earlier run.

→ Record the actual `length(markers)` in the repository, or restate Table S34's title to
match the true count. This is a one-line fix with a checkable number attached.

## 16. Two conflicting definitions of `theta` in one file — OPEN

`第 2 步 v4 …` computes the reference signature fraction two different ways, under the
same variable name:

| Definition | Where | Formula |
|---|---|---|
| Mean of per-sample pseudobulk libraries | lines 44–48 | `rowMeans(ref_sig[, ph$ct_merged == ct])` |
| Sum of raw counts over all cells of that type | lines 209–212, 319–322, 1224–1227 | `Matrix::rowSums(cts[, idx])` |

Both are named `sig_frac`/`theta`. In an interactive session, whichever block ran last
determines the value. The two are not equivalent when cell numbers differ between
patients.

→ Determine which definition produced the published Tables S33–S35, and state it once in
the README. This is exactly the class of ambiguity the "documented in the repository"
promise in the Code availability statement is meant to cover.

## 17. Four-arm comparison: no arm 4 identifier, and Fig S2 is unreferenced — DOC

- The code defines `arm1_*`, `arm2_*`, `arm3_*` only. **No `arm4*` identifier exists.**
  The "four-arm" phrasing appears only as prose, and the fourth arm is
  `arm3_logFC_adj` — this script's own limma model adjusted for MuSiC proportions, as
  distinct from `arm3_logFC_raw`. The four arms are thus: (1) P0 cohort effect,
  (2) P0 composition-adjusted summary, (3) this script unadjusted, (4) this script adjusted.
- **Figure S2 is not referenced anywhere in the deconvolution code.** The strings
  `Fig`, `Figure`, `图` do not appear in the file; only Tables S32–S36.
- `GSE174478` is **not** part of the P0 training set, and its endpoint is numeric
  `fibrosis_stage` rather than case/control, so arms 1 and 2 are `NA` for it. NOTE 16
  states this, and the README should too, so the NA cells are not read as failures.

→ All three are documentation items. The "four-arm" label is defensible once the four
estimates are enumerated.

## 18. Supplementary workbook is built by positional indices — DOC

`Supplementary_Tables_v2.xlsx` is assembled by assuming the exact shape of a workbook
built in a different file:

```r
writeData(wb, "INDEX", idx_new, startRow = 33, ...)   # assumes exactly 32 INDEX rows
writeData(wb, "NOTES", nt_new, startRow = 14, ...)    # assumes exactly 13 NOTES rows
names(wb2)[33:38]                                     # assumes new sheets land at 33–38
```

Any change to the upstream workbook breaks these silently or with an error. Also,
`第 2 步 v4 …` reads other working scripts by timestamped filename (lines 682–688, 733,
783, 958) — including a probe to confirm that `正交去卷积敏感性分析MuSiC 主跑…R` really was
empty (`# ---- L0. 确认"备好脚本"第二个是否真空 ----`, line 732).

→ The published repository should not depend on positional workbook state. The README
should state the expected sheet and row counts as preconditions.

## 19. GSE174478 is analysed but missing from Data availability — OPEN

The deconvolution sensitivity analysis uses **four** bulk cohorts:

```r
mapped  # GSE126848, GSE130970, GSE162694, GSE174478
```

Confirmed by `stopifnot(all(c("GSE126848","GSE130970","GSE162694","GSE174478") %in%
names(bulks)))` in the deconvolution asset scripts.

**GSE174478 does not appear in the manuscript's Data availability statement**, which
lists 13 accessions (GSE48452, GSE89632, GSE83452, GSE126848, GSE130970, GSE162694,
GSE116312, GSE27411, GSE135251, GSE153224, GSE202379, GSE115469, GSE134520). It also
does not appear in `data/datasets.md` as originally drafted, because it was not in the
project handover's dataset table either.

It *is* mentioned in the supplementary NOTES as the fibrosis-gradient cohort
(NOTE 16: "GSE174478 was not part of the P0 training set … its endpoint is numeric
fibrosis stage 0-4"), and it is the subject of note 5 about GTEx V10 in Table S32. So the
analysis is reported; only the data-availability list is incomplete.

→ Add GSE174478 to the Data availability statement. Nature Portfolio's policy expects
every dataset analysed to be identifiable, and this one is analysed, reported in the
supplementary tables, and drawn in Figure S2 (the dagger-marked fibrosis-gradient arm).

Two related framing points from the same NOTES entry, worth checking in the manuscript:

- GSE174478 was **not** part of the P0 training set, so comparison arms 1 and 2 are `NA`
  for it. The main text should not present it as a discovery cohort.
- Its endpoint is a **numeric fibrosis stage**, not case/control, so its effect estimate
  is on a different scale from the other three cohorts. Figure S2 shows this with a
  dagger rather than labelling the axis as stage units.

## 20. Figure renames were done outside R — DOC

Three figures were renamed by hand between the script's output filename and the
submitted filename: `fig1_localization_heatmap.png` → Fig 3,
`fig3_cross_organ_schematic.pdf` → Fig 4, and `fig2_fibrosis_boxplots.pdf` → Fig 5
(superseded). No `file.rename()` or `file.copy()` call exists in any supplied script.

The consequence is that a reader who opens the figure script and sees `fig1_...` will
reasonably assume it is Figure 1. It is Figure 3.

→ This is fully resolved by `figure_map.md`. Flagging it here because it is the kind of
mismatch that looks like carelessness until it is documented, at which point it looks
like transparency.

---

## Items the code settles (no action needed)

| Item | Established value |
|---|---|
| RRA P-value implementation | `min_k pbeta(x_(k); k, n − k + 1)` — `磁盘层执行脚本` line 36, matches the manuscript exactly ✓ |
| Meta-analysis | `metafor::rma(yi = logFC, sei = SE, method = "REML")`, grouped chip vs RNA — line 126 ✓ |
| LODO procedure | full re-ranking with each of 6 cohorts dropped; retention criteria `in200` (rank ≤ 200) and `lt001` (P < 0.01) — lines 71–78 ✓ |
| Per-cohort SE | `SE = logFC / t` — line 108 ✓ |
| Gastric DEG rule | `P.Value < 0.05 & abs(logFC) > 0.5` — `第一步数据下载（CAG）` line 280 ✓ |
| Gastric supporting threshold | `P.Value < 0.05 & abs(logFC) > log2(1.5)` — `磁盘层执行脚本` line 91 ✓ |
| AUC definition | Mann–Whitney rank statistic, normalised as `pmax(AUC, 1 − AUC)` — `第一步第二部分GSE135251` line 222 ✓ |
| AUC bootstrap | 2,000 resamples, percentile CI — `第一次预审稿文件处理` line 251 ✓ |
| Spearman CI | Fisher z, `se = 1/sqrt(n − 3)` — line 230 ✓ |
| MR instrument filter | per-SNP `F > 10`, R² from beta/EAF/samplesize — `eQTLGen… 第三次跑` lines 97–102 ✓ |
| MR method | Wald ratio for 1 SNP (manual formula), otherwise IVW of 4 methods — lines 105–137 ✓ |
| MR clumping | instruments from a directory named `eqtlclump_r2_0.001` — i.e. r² = 0.001 — all MR scripts ✓ |
| MR MR-PRESSO | `NbDistribution = 1000` — line 126 ✓ |
| CVD risk / Steiger | **removed** in v3 and absent from v4 — v4 is the reported version ✓ |
| GSE135251 pipeline | counts → genes with ≥5 non-zero samples → `log2(CPM + 1)` — `补缺口脚本` lines 106–107 ✓ |
| GSE135251 groups | 206 disease vs 10 control ✓ |
| Fibrosis stage per patient | **mode** of per-cell `Fibrosis.score..F0.4.` — `第一次预审稿文件处理` line 731 ✓ |
| Fig 5 locked values | IL32/Hepatocytes ρ = +0.639; CADM2/Stellate ρ = −0.761; ANXA4/Cholangiocytes ρ = −0.585; stage counts F0–F4 = 5/9/12/12/9 ✓ |
| GSE153224 design | 10 pooled libraries, 4 patients per pool, statistics at library level (5 vs 5) ✓ |
| MuSiC installation | v1.0.0 from `github.com/xuranw/MuSiC`, not Bioconductor ✓ |
| Sample sizes (6 liver cohorts) | Ahrens 26/28, Arendt 39/24, Lefebvre 104/44, Suppli 31/26, Hoang 53/25, Pantano 112/31 ✓ |
| Chip matrix dimensions | Ahrens 54, Arendt 63, Lefebvre 148 samples ✓ |

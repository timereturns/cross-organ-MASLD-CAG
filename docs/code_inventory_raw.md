# Code inventory: what was received, and what is missing

Audit date: 2026-09-28. Auditor: automated code review pass over 21 author-supplied
`.R` working files (~8,300 lines).

The author ran these analyses interactively in RStudio. The files are **working logs**:
each one contains several RStudio sessions' worth of code separated by blank lines,
with blocks that were superseded mid-file, hand-substituted values, and output that
was read back into the next session. They are not a clean pipeline. This document maps
what exists, then lists what is missing and what is a reproducibility risk.

---

## 1. What was received

| Original filename | Lines | Role | Status |
|---|---|---|---|
| `第一步数据下载（MASLD侧 差异分析已完成）.R` | 737 | Liver cohort download + limma DE + **RRA discovery** | ✅ load-bearing |
| `第一步数据下载（CAG）.R` | 393 | Gastric cohort download + limma DE + **candidate finalisation** | ✅ load-bearing |
| `补缺口脚本9.12 16 10….R` (`handover_assets.R`) | 249 | Rebuilds expression matrices, groups, GSE135251 data, RRA rank columns | ✅ load-bearing |
| `磁盘层执行脚本2026.9.12 16 56.R` (`P0_audit_disk.R`) | 180 | **LODO exact RRA + CAG support + per-cohort effects + meta** | ✅ load-bearing |
| `P0组成校正2026.9.12 17 14.R` | 714 | Composition adjustment, batch diagnosis, RRA re-derivation | ✅ load-bearing |
| `第一次预审稿文件处理.R` | 1087 | Post-review evidence: bootstrap CIs, core matrix, Supplementary workbook build | ✅ load-bearing |
| `正交去卷积敏感性分析2026.9.14 05 16.R` | 679 | Gene-symbol mapping for deconvolution; MuSiC **full-gene** run | ✅ load-bearing |
| `第 2 步 v4 bulk 性质闸加正确模型 QC….R` | 1365 | Bulk gate + full-gene MuSiC + **marker-restricted fallback** + Tables S32–S36 | ✅ load-bearing |
| `P1 主图绘制 2026.9.12 19 55.R` | 1580 | Figure generation (many superseded drafts) | ✅ figures |
| `P1 执行序 2026.9.12 18 34.R` | 29 | Package installation only | ❌ scaffolding |
| `P2 开始 2026.9.12 22 04.R` | 132 | 8-dimension candidate scoring matrix (A–H) | ✅ load-bearing |
| `P3 脚本（CADM2 × LiverPDFF × GTEx Liver）…R` | 228 | Colocalisation follow-up | ⚠️ exploratory |
| `eQTLGen_eqtl 结局筛选与确定….R` | 70 | Outcome GWAS file reconnaissance + format conversion | ⚠️ exploratory |
| `eQTLGen_eqtl 主分析2026.9.12 14 38.R` | 250 | MR v2 (has `steiger_filtering`) | ❌ superseded |
| `eQTLGen_eqtl 主分析第二次跑….R` | 147 | MR v3 (Steiger removed, manual Wald) | ❌ superseded |
| `eQTLGen_eqtl 主分析第三次跑….R` | 168 | **MR v4 — the version reported in the paper** (`mr_all_results_v4.csv`) | ✅ load-bearing |
| `跑 compact 版三件套2026.9.12 19 21.R` | 448 | GSE202379 / GSE115469 / GSE134520 single-cell handling | ✅ load-bearing |
| `第一步第二部分GSE135251 独立验证.R` | 406 | GSE135251 held-out validation | ✅ load-bearing |
| `预审稿后图片重绘2026.9.20 19 36.R` | 206 | Fig 5 v3 redraw (fibrosis stages) | ✅ figures |
| `落盘 最终核对表 2026.9.13 15 12.R` | 184 | HPA check table + final checklist | ✅ load-bearing |
| `结果部分确认2026.9.14 19 23.R` | 82 | Results-number spot checks | ⚠️ audit only |
| `正交去卷积敏感性分析MuSiC 主跑2026.9.14 05 27.R` | **0** | ⚠️ **EMPTY FILE** | ❌ missing |

Note: three files are successive runs of the same MR analysis (v2 → v3 → v4). **Only
v4 corresponds to the manuscript.** v2 and v3 are historical and are archived, not
presented as the analysis.

---

## 2. The discovery chain, reconstructed

The 14 candidate genes were derived by this route:

```
6 liver cohorts: limma per cohort  →  rank by t-statistic (up / dn separately)
        ↓
RobustRankAggreg::aggregateRanks(glist, N = length(union of symbols))
        ↓
filter  Score < 0.01  and  ≥3 of 6 cohorts in the same direction
        ↓
top 200 up  +  top 200 dn            →  MASLD_RRA_up_final.csv / _dn_final.csv
        ↓
2 gastric cohorts: limma per cohort; DEG = P.Value < 0.05 and |logFC| > 0.5
        ↓
MASLD RRA list ∩ gastric DEGs        →  shared_genes_concordant.csv
        ↓
require direction agreement in BOTH gastric cohorts
require ≥4 of 6 liver cohorts in the stated direction
        ↓
                                     →  candidate_genes_final.csv  (14 genes)
```

Source: `第一步数据下载（MASLD侧 差异分析已完成）.R` lines 605–718, and
`第一步数据下载（CAG）.R` lines 279–382.

### ⚠️ Background size discrepancy — must be resolved

The two sizes appear in the code and they are not the same number:

| Value | Where | What it is |
|---|---|---|
| `N = length(all_sym)` | discovery script, line 613 | union of symbols across the six liver cohorts, computed at run time |
| `N = 11982` | `P0组成校正`, line 402 (`N_opt$common`) | the six-cohort common-gene count |
| `N = 20000` | `磁盘层执行脚本`, line 45 | **the audited choice**, used for the LODO table |

The discovery run used the **runtime union**, not a hard-coded 20,000. The value
**20,000** enters later, in the robustness audit, where three definitions were
compared and 20,000 was manually selected (`磁盘层执行脚本` line 51:
`use_N20000 <- TRUE  # 校准后手工指定`). The comment records a **hand calibration**,
not an automatic criterion — this needs to be documented explicitly, because a
reviewer comparing the manuscript's "20,000" against a script that computed
`N <- length(all_sym)` will find a contradiction.

---

## 3. Missing code

### 3.1 ⛔ CRITICAL — permutation testing

The manuscript Introduction refers to permutation testing, the Code availability
statement promises "the random seeds used for permutation testing", and Supplementary
Table S02 is titled **"Leave-one-dataset-out exact RRA (20,000 permutations)"**.

**No permutation loop exists in any of the 21 files.** A full-tree search for
`set.seed`, `sample(`, `replicate(`, `置换`, `permut`, `nperm` returns only:

- `第一次预审稿文件处理.R` lines 206 and 352: `set.seed(20260914)`
- `第一次预审稿文件处理.R` lines 251 and 404: `B <- 2000` — these drive the
  **AUC bootstrap** (Table S39), a different analysis
- `run_mr_presso(..., NbDistribution = 1000)` in the three MR scripts — MR-PRESSO,
  also unrelated
- the string `"20,000 permutations"` in a table title at line 442

So the seed used by the permutation test is recoverable
(**`set.seed(20260914)` appears at the top of the post-review blocks**), but the
permutation code itself and its iteration count are not present. The permission
setting must be confirmed against the script that actually produced Table S02.

**Action required:** supply the script that generated `p0_lodo_exact_rra.csv` /
Table S02, or confirm that `set.seed(20260914)` with 20,000 permutations is correct.

### 3.2 ✅ RESOLVED — marker-restricted MuSiC is present after all

**Correction to an earlier draft of this inventory.** The marker-restricted fallback is
not missing: it lives in `第 2 步 v4 bulk 性质闸加正确模型 QC….R` (1365 lines, not
1138 as first counted — the file tail was missed on the first pass).

The code behind Tables S32–S36 is:

| Element | Location | Value |
|---|---|---|
| Candidate exclusion list | line 324 | `CAND <- c("CDHR2","ANXA4","CADM2","IL32","LGALS3","RPS6KA1")` |
| Specificity score | line 328 | `spec = theta / (rowMeans(theta) + 1e-12)` |
| Marker selection | lines 344–349 | top **50 per cell type** by specificity |
| Anti-circularity assertion | line 347 | `stopifnot(!any(CAND %in% markers))` |
| Mitochondrial exclusion | line 333 | `grep("^MT-|^MTRNR", …)` |
| Reference | lines 353–363 | patient × cell-type **mean** pseudobulk, GSE202379, 8 cell types |
| Full-gene run (the one that failed) | lines 63–66 | `markers = common` (all shared genes) |
| Marker-panel run (the accepted one) | lines 372–375 | `markers = intersect(intersect(bulk, ref), markers)` |
| Centered control | lines 390–393 | `centered = TRUE`, GSE126848 only |

**Two caveats that must be documented rather than fixed silently:**

1. **The "400" is prose, not a computed value.** The panel is built with
   `unique(unlist(...))` over top-50 lists, so any gene that is top-50 for more than one
   cell type is counted once. The true panel size is therefore **≤ 400** and is never
   printed, stored, or asserted. Table S34's title *"400 specificity-selected genes"* and
   NOTE 15's *"400 specificity-selected markers"* are hardcoded strings.
2. **`excl_skip` explains RPS6KA1.** NOTE 17 says RPS6KA1 is absent from the reference
   gene rows, so it is skipped from the exclusion list — meaning it was never in the
   reference to begin with, and no reference-side adjustment exists for it. Consistent.

The 0-byte file `正交去卷积敏感性分析MuSiC 主跑2026.9.14 05 27.R` remains empty, but the
script it was apparently meant to become was written into the bulk-gate file instead.
Nothing is lost.


### 3.3 Missing upstream scripts and derived assets

| Missing item | Needed for | Referenced by |
|---|---|---|
| `GSE153224` external validation script | Fig 2B, Table S25 | `P1 主图绘制` line 936 reads `GSE153224_mRNA_Expression_Profiling.xlsx`; source dir `MASLD_CAG_followup_v5/` never supplied |
| GSE202379 preprocessing → `p1_gse202379_compact.rds` | every GSE202379 analysis | `跑 compact 版三件套` line 6; the compact object is read, never built |
| `B_reference_merged.rds` construction | MuSiC reference | read at `正交去卷积…` line 32; build script absent |
| `B_bulks_loaded.rds` construction | MuSiC bulk arm | read at line 41; build script absent |
| `GSE202379_患者级协变量与纤维化分期表.csv` | Fig 5 | `预审稿后图片重绘` line 20; source absent |
| Per-cohort limma result objects | everything downstream | `chips_limma_results.rds`, `rnaseq_limma_results.rds`, `cag_limma_results.rds` are read; the producing blocks are partly present in the download scripts but not as a single reproducible step |

### 3.4 Derived result files that are inputs to other scripts

These are read by supplied scripts but were not supplied. They are small CSV files and
should be committed to the repository as analysis outputs:

`candidate_genes_final.csv` · `MASLD_RRA_up_final_ranked.csv` ·
`MASLD_RRA_dn_final_ranked.csv` · `shared_genes_concordant.csv` ·
`shared_validation_auc.csv` · `p0_lodo_exact_rra.csv` · `p0_cag_support.csv` ·
`p0_cohort_effect_sizes.csv` · `p0_meta_per_platform.csv` ·
`p2_scoring_matrix_final.csv` · `p1_gse202379_localization.csv` ·
`p1_gse202379_pseudobulk_cpm.csv` · `p1_gse202379_fibrosis_cor.csv` ·
`p1_gse115469_localization.csv` · `p1_gse134520_localization.csv` ·
`p1_gse134520_pseudobulk_cpm.csv` · `mr_all_results_v4.csv`

---

## 4. ⛔ Blocking issue: third-party code provenance

`第一步数据下载（MASLD侧 差异分析已完成）.R` opens with this notice, verbatim:

```
#数据与代码声明
#如果没有购买SCI狂人团队或者生信狂人团队的正版会员
#没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
#如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿
####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
```

Translation: *"Data and code statement — if you have not purchased a genuine membership
from the SCI Kuangren / Bioinformatics Kuangren team, and have not obtained our
permission, using the data and code we have organised to publish a paper constitutes
unauthorised use; if we discover your paper used our data and code, we will use every
means to have your paper retracted."*

This header is attached to **the script that performs the liver-side download, limma
differential expression and the RRA discovery screen** — the foundation of the paper's
main results.

**Why this blocks publication of the repository:**

1. A public GitHub repository plus a Zenodo DOI under an MIT licence asserts the right
   to redistribute and relicense. If any part of this file is vendor-supplied, that
   assertion is not the author's to make.
2. The notice threatens retraction. Publishing the file publicly, with the notice either
   attached or removed, converts a dormant risk into an active one — and the manuscript
   will point reviewers directly at this repository.
3. Scientific Reports is a Nature Portfolio journal and takes research-integrity
   complaints seriously. A vendor claiming licence breach against a repository cited in
   the Code availability statement is a materially worse outcome than a slow submission.

**This must be resolved by the author before any push or upload.** Not an AI decision.
The available routes are:

- **(a)** Establish that the notice is boilerplate appended to a member-only data
  package, and that the analysis code below it is the author's own writing — get this
  in writing from the vendor if possible.
- **(b)** Rewrite the affected analysis from scratch so that the repository contains
  only the author's own code, and document the re-implementation.
- **(c)** Scope the repository to exclude the affected stages and state the scope
  limitation in the Code availability statement.

**Interim measure applied:** the offending header is quoted here for the record, and
the affected file is flagged. Nothing has been pushed anywhere.

---

## 5. Other reproducibility risks to disclose in the README

| Risk | Detail | Where |
|---|---|---|
| Hand-calibrated parameter | `use_N20000 <- TRUE  # 校准后手工指定` | `磁盘层执行脚本` line 51 |
| Hand-entered scoring weights | The G/H dimensions of the candidate scoring matrix are hard-coded numeric vectors typed by the author | `P2 开始` lines 105–113 |
| Hand-entered HPA table | Protein levels, antibody IDs, reliability grades typed as a `data.frame` literal | `落盘 最终核对表` lines 2–20 |
| Hand-copied effect size | `b_exp <- -0.3213241; se_exp <- 0.008639812` pasted from a previous session | `P0组成校正` lines 267–268 |
| Two stage definitions compared, one chosen post hoc | `stage_mode` vs `stage_mean`; the code picks whichever reproduces three stored ρ values | `第一次预审稿文件处理` lines 729–772 |
| Comparator definition is "mode OR rounded mean" | Selection rule is reproduction of stored numbers, not a pre-registered criterion | same |
| Threshold chosen after inspection | `Score < 0.05` → `Score < 0.01`, and top-200 truncation, appear as successive tightenings | `第一步数据下载（MASLD侧）` lines 621–706 |
| Superseded analyses in-file | Multiple abandoned versions of the same block remain in the working files | throughout |
| Platform-specific paths | 100+ occurrences of `E:/生信狂人/...` and `C:/Users/LEGION/...` | every file |

None of these are fatal. All of them are the kind of thing a careful reviewer asks
about, and all are answerable as long as the repository documents them rather than
leaving them to be discovered. The README drafts in this directory do so.

---

## 6. Summary of required author input

Blocking, in priority order:

1. **Third-party code provenance** — resolve before anything is published (§4).
   Author decision: the affected script is held out of the public repository.
2. **Permutation test script**, or confirmation that the seed is `20260914` and the
   count is 20,000 (§3.1). Author reports this script is not available.
3. **`B_b4_compare.csv`** — must contain the `verdict` and `note` columns, or the
   TableS32 block fails. See `manuscript_discrepancies.md` §12.
4. **The derived result CSVs** listed in §3.4.
5. GSE153224, GSE202379-preprocessing and deconvolution-asset scripts (§3.3), if
   they exist.

Resolved since the first pass: marker-restricted MuSiC (§3.2) is present.


# Original working files — manifest and provenance

The analysis was carried out interactively in RStudio. The code arrived as 21 working
`.R` files, ~8,300 lines, named with the author's own conventions (Chinese working names
with RStudio autosave timestamps appended). They are logs of many sessions rather than a
pipeline: each file contains several blocks separated by blank lines, each block
typically opening with `setwd()` and device cleanup, with superseded revisions left in
place.

This directory records what those files are, so that the organised `R/` scripts can be
traced back to their sources. The working files themselves are **not** redistributed
here — see the provenance note below.

Column meanings: **role** = where the file sits in the pipeline; **status** = whether the
organised repository depends on it.

---

## Inventory

| # | Original filename | Lines | Role | Organised script | Status |
|---|---|---|---|---|---|
| 1 | `第一步数据下载（MASLD侧 差异分析已完成）.R` | 737 | Liver cohort download, limma DE, **RRA discovery** | `01`, `02`, `04` | ⛔ **excluded — see provenance** |
| 2 | `第一步数据下载（CAG）.R` | 393 | Gastric download, limma DE, **candidate finalisation** | `01`, `03` | ✅ |
| 3 | `补缺口脚本9.12 16 10（现在跑 一次补齐所有缺失资产）.R` | 249 | Rebuilds expression matrices, groups, GSE135251 data, RRA rank columns | `04` | ✅ |
| 4 | `磁盘层执行脚本2026.9.12 16 56.R` | 180 | **LODO exact RRA**, CAG support, per-cohort effects, meta-analysis | `05` | ✅ |
| 5 | `P0组成校正2026.9.12 17 14.R` | 714 | Marker-score composition adjustment, batch diagnostics, RRA re-derivation | `06` | ✅ |
| 6 | `第一步第二部分GSE135251 独立验证.R` | 406 | GSE135251 held-out validation | `07` | ✅ |
| 7 | `跑 compact 版三件套2026.9.12 19 21.R` | 448 | GSE202379 / GSE115469 / GSE134520 single-cell handling | `08` | ✅ |
| 8 | `正交去卷积敏感性分析2026.9.14 05 16.R` | 679 | Gene-symbol mapping for deconvolution; full-gene MuSiC | `09` | ✅ |
| 9 | `第 2 步 v4 bulk 性质闸加正确模型 QC 2026.9.14 14 34.R` | 1365 | Bulk gate, full-gene MuSiC, **marker-restricted fallback**, Tables S32–S36 | `09` | ✅ |
| 10 | `P2 开始 2026.9.12 22 04.R` | 132 | 8-dimension candidate scoring matrix | `10` | ✅ |
| 11 | `eQTLGen_eqtl 主分析第三次跑2026.9.12 15 10.R` | 168 | **MR v4 — the reported version** | `11` | ✅ |
| 12 | `eQTLGen_eqtl 主分析2026.9.12 14 38.R` | 250 | MR v2 (with Steiger filtering) | — | ❌ superseded |
| 13 | `eQTLGen_eqtl 主分析第二次跑2026.9.12 14 56.R` | 147 | MR v3 (Steiger removed, manual Wald) | — | ❌ superseded |
| 14 | `eQTLGen_eqtl 结局筛选与确定2026.9.12 13 56.R` | 70 | Outcome GWAS reconnaissance and format conversion | `11` (partial) | ⚠️ exploratory |
| 15 | `P3 脚本（CADM2 × LiverPDFF × GTEx Liver）2026.9.12 22 13.R` | 228 | Colocalisation follow-up, eQTL Catalogue | — | ⚠️ exploratory |
| 16 | `落盘 最终核对表 2026.9.13 15 12.R` | 281 | HPA check table; supplementary workbook assembly | `12` | ✅ |
| 17 | `P1 主图绘制 2026.9.12 19 55.R` | 1958 | Figure generation, 18 blocks, many superseded | `13` | ✅ |
| 18 | `预审稿后图片重绘2026.9.20 19 36.R` | 206 | Fig 5 v3 redraw with locked-value assertions | `14` | ✅ |
| 19 | `P1 执行序 2026.9.12 18 34.R` | 29 | Package installation only | — | ❌ scaffolding |
| 20 | `结果部分确认2026.9.14 19 23.R` | 96 | Results-number spot checks; prints only, no assertions | — | ❌ scaffolding |
| 21 | `正交去卷积敏感性分析MuSiC 主跑2026.9.14 05 27.R` | **0** | Empty file | — | ❌ empty |

## Second submission (2026-09-28, later)

Five files arrived after the first inventory. Three were **byte-identical** to files
already analysed — `第一步第二部分GSE135251 独立验证.R` (sha256 `5c809d4c…`),
`第一步数据下载（CAG）.R` (`60327c63…`) and
`第一步数据下载（MASLD侧 差异分析已完成）.R` (`83249b81…`) — and required no new work.

The two new files:

### `clump自查 - 8循环重跑9.12.R` — 146 lines — ✅ genuinely useful

It records the MR instrument selection parameters that were previously only implicit
in a directory name:

```r
ieugwasr::ld_clump(data.frame(rsid = dat$SNP, pval = dat$p),
                   clump_kb = 10000, clump_r2 = 0.001,
                   clump_p = 5e-8, pop = "EUR")
```

window **10,000 kb** · **r² < 0.001** · **p < 5×10⁻⁸** · **European panel**.

It also documents three things that matter for reading the MR results:

1. **The r² = 0.01 instrument set was purchased.** The comment identifies
   `eqtlclump/` as *"买来的 0.01 文件所在目录"* — the directory holding the purchased
   r² = 0.01 files — and the r² = 0.001 set was derived from it here. A second
   instance of third-party content in the pipeline, though unlike file 1 the derived
   set is the study's own work.
2. **TSPAN3's instruments were repaired by hand.** Re-clumping returned **zero** SNPs
   (guarded by `if (nrow(kept) == 0) stop("clump 结果为空，停下检查")`), then one
   surviving LD pair still exceeded r² = 0.001 and the lower-F SNP was dropped,
   guarded by `stopifnot(nrow(kept2) == nrow(d) - 1)`.
3. **LGALS3 and CADM2 use a single instrument each**, from the looser purchased set,
   so every estimate for them is a Wald ratio with no Egger intercept, no Q and no
   MR-PRESSO. **CADM2 is one of the five retained genes**, which makes this worth
   stating wherever its MR support is discussed.

A fourth item: this script computes instrument strength as `F = (beta/se)²`, the
univariate form, while the MR run uses the sample-size-aware R² form. The two
disagree, and `analysis_parameters.md` §9 now records both.

### `test.R` — 25 lines — ❌ not useful

A practice script on the built-in `mtcars` dataset, containing a malformed code fence
(`` ` ``) inside the R source at line 11. No relation to this project.

---

Line counts are as read during the audit; some differ from the file's byte-derived count
because of trailing blank lines.

---

## ⛔ Provenance note — file 1 is excluded

`第一步数据下载（MASLD侧 差异分析已完成）.R` opens with a vendor notice stating that use
of the vendor's organised data and code to publish a paper without permission constitutes
unauthorised use, and threatening retraction. The notice is attached to the file that
performs the liver-side download, limma differential expression and the rank-aggregation
discovery screen.

The author's position is that the analysis code is their own work and the notice is
boilerplate attached to a data package; that position has not been independently
confirmed. **The file is therefore not redistributed here, and no licence is asserted
over any part of it.**

The corresponding pipeline stage is implemented in `R/01_download_and_DE.R` and
`R/02_RRA_discovery.R` from the public GEO source files, and the resulting candidates and
rankings are reported in `results/tables/`. The author can provide the original file to
the editors and reviewers on request.

---

## Notes on the reorganised scripts

1. **Blocks were merged by analysis stage, not by file.** Three of the original files
   each span several pipeline stages — file 1 alone covers download, differential
   expression, rank aggregation and the LODO audit — so the organised scripts do not map
   one-to-one onto the originals. The table above records the mapping.

2. **`setwd()` calls were removed and replaced by `R/00_setup.R`.** The originals contain
   100+ occurrences of machine-specific absolute paths including
   `E:/生信狂人/190_快速筛选主角基因/key_gene`, and one file repeats the same `setwd()` 18
   times. The organised scripts resolve everything relative to the repository root.

3. **Runtime `install.packages()` calls were removed.** Several original blocks install
   packages mid-analysis, including Bioconductor and GitHub sources and a hardcoded CRAN
   mirror. Dependencies are now declared in the README.

4. **Superseded revisions were dropped.** The figure script alone contains four
   superseded revisions and one abandoned block with no output. The MR pipeline exists in
   three versions; only the third produced the reported results. Supersession is recorded
   in `docs/analysis_parameters.md`.

5. **Blocks that depend on in-session objects were made self-contained.** Several blocks
   — including the one that produces the submitted Fig 3 — reference objects created in
   an earlier RStudio session and fail in a fresh session. These were rewritten to load
   their inputs explicitly.

6. **Hand-entered values were preserved, not silently replaced.** Where a number was
   typed by hand rather than computed — the G and H dimensions of the scoring matrix, the
   Human Protein Atlas table, two columns of the four-arm comparison table, several rows
   of the deconvolution QC table — it remains, and
   `docs/manuscript_discrepancies.md` marks it as hand-entered. Undocumented manual
   values are a reproducibility problem; documented ones are a curation record. The
   distinction is what matters.

7. **Scripts that could not be recovered are listed in `docs/repository_scope.md`** with
   the reason and the scope of the gap.

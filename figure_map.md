# Figure map

Which script produced which figure, and where the file went.

> ⚠️ **Read this first.** The working scripts number their figures **differently** from
> the manuscript, and the renames were done by hand outside R (no `file.rename()` call
> exists anywhere). Never take a filename from a script at face value.
>
> | Script's own filename | Manuscript figure |
> |---|---|
> | `fig1_localization_heatmap.png` | **Fig 3** |
> | `fig3_cross_organ_schematic.pdf` | **Fig 4** |
> | `fig2_fibrosis_boxplots.pdf` | **Fig 5** (superseded — see below) |
> | `Fig1_cross_cohort_audit.pdf` | **Fig 1** |
> | `Fig2_external_validation.pdf` | **Fig 2** |
> | `figS1_cag_nag_bars.pdf` | **Fig S1** |
> | `figS_B4_attenuation.pdf` | **Fig S2** |

---

## Main figures

| Manuscript | Script (repository) | Output written by the script | Evidence for the mapping |
|---|---|---|---|
| **Fig 1** — cross-cohort consistency of the 14 candidates across 10 cohorts | `R/12_figures.R` — final block, `★ 修` | `Fig1_cross_cohort_audit.pdf` + `.png` (300 dpi) | 14 × 10 heatmap; BH over the flattened matrix gives 140 tests, matching the legend. `TSPAN3` moved to *Demoted* and `gaps_row = c(1,3,4,5,6)`, `gaps_col = 7` in the final block only. |
| **Fig 2** — replication in the two held-out cohorts | `R/12_figures.R` — `★ v4 修复版` block | `Fig2_external_validation.pdf` + `.png` | Panel A GSE135251 (n = 216), Panel B GSE153224 (10 pooled libraries, 5 vs 5). The legend's "arrows next to the gene names" exist only in the `ggtext` version, which is the last block in the file. |
| **Fig 3** — cell-type localisation in three snRNA-seq datasets | `R/12_figures.R` — `统一配色` block | `fig1_localization_heatmap.png` → renamed by hand to `Fig3_localization_heatmap.png` | 14 genes × cell types of L1/L2/St, within-dataset z-scored; grey = gene not captured. On-disk `Fig3_localization_heatmap.png` is **86,727 B, 2026-09-12 20:58**, byte-identical in size and timestamp to the submitted `Fig 3.png`. |
| **Fig 4** — cross-organ schematic | `R/12_figures.R` — `v6` block | `fig3_cross_organ_schematic.pdf` → renamed by hand to `Fig 4.pdf` | The v6 node label `"bulk-detectable,\nlow abundance"` and the annotation `"(cell of origin unresolved)"` match the legend's "dashed curve and dashed box". |
| **Fig 5** — candidates across fibrosis stages | ⚠️ **`R/14_fig5_redraw.R`**, *not* the main figure script | `fig5_fibrosis_boxplots.png` → submitted as `Fig 5.png` | Submitted `Fig 5.png` is **315,435 B, 2026-09-20 20:01:22** — byte-identical in size and timestamp to the output of the later redraw script. |

### Fig 5 needs its own note

The main figure script **does** contain a fibre-stage block (a base-R version and a
ggplot version), but neither is the submitted figure. The submitted Fig 5 was produced
later by `预审稿后图片重绘2026.9.20 19 36.R`, which:

- re-reads the same `p1_gse202379_compact.rds` and `p1_gse202379_pseudobulk_cpm.csv`,
- **additionally** reads `GSE202379_患者级协变量与纤维化分期表.csv` — the authoritative
  stage table, which is a different source from the modal stage computed in the main
  script,
- uses `ppcor` for the partial-correlation check,
- asserts the locked values: stage counts F0–F4 = 5 / 9 / 12 / 12 / 9, and the three
  lock-in ρ values (IL32/Hepatocytes +0.639, CADM2/Stellate −0.761,
  ANXA4/Cholangiocytes −0.585), failing if they do not reproduce,
- writes `cairo_pdf` for embedded fonts, for submission safety.

That is why Fig 5 is listed against a separate script: the two stage definitions were
compared post hoc and the one that reproduced the stored numbers was adopted. This is
documented rather than hidden — see `docs/manuscript_discrepancies.md`.

## Supplementary figures

| Manuscript | Script | Output |
|---|---|---|
| **Fig S1** — GSE134520 gastric pseudobulk, descriptive | `R/12_figures.R` — `FigS1 优化版` block | `figS1_cag_nag_bars.pdf` (4 facets, one per gene, cell types G_epithelium / Pit / Neck / Chief / Myeloid, log2FC labels above bars) |
| **Fig S2** — MuSiC-adjusted effect attenuation | `R/12_figures.R` — `补充图 优化版 v2` block | `figS_B4_attenuation.pdf` + `.png`. On-disk `.png` is **197,981 B, 2026-09-14 16:03:28** — byte-identical to the submitted `Fig S2.png`. 50 % attenuation reference line, grey dashed boxes for no-signal and flipped cells, dagger marking the fibrosis-gradient arm. |

⚠️ **The Fig S2 block rewrites its own input.** It reads
`B_b4_compare.csv`, appends two hand-typed columns (`verdict`, `note`), and writes the
file back in place. See `docs/manuscript_discrepancies.md` §12 — this is
**the** reproducibility hazard in the repository, and the reason those two columns exist
in the published Table S32.

---

## Figures the script draws but which were superseded

Kept here so nobody re-derives them by accident:

| Superseded output | Superseded by |
|---|---|
| `fig1_localization_heatmap.png` — first version, 5-stop RdBu | the `统一配色` block |
| `fig3_cross_organ_schematic.pdf` — base-R version | the ggplot2 `桑基风格` version, then `v6` |
| `figS1_cag_nag_bars.pdf` — base-R version | the faceted ggplot version |
| `fig2_fibrosis_boxplots.pdf` — base-R, 6 panels | the ggplot `优化版`, which is itself superseded by `R/14_fig5_redraw.R` |
| `Fig1_cross_cohort_audit.pdf` — two earlier revisions | the final `★ 修` block |
| `Fig2_external_validation.pdf` — base-R version, then ggplot v2 | the `★ v4 修复版` block |

The file also contains one **abandoned block with no output at all** — a third copy of
the Fig 2 header that defines its inputs and helper functions and then stops without
drawing anything.

---

## Figures that no script references

**Figure S2 is not referenced in the deconvolution scripts.** The strings `Fig`,
`Figure` and `图` do not appear in the MuSiC code at all; that pipeline produces tables
S32–S36 only. The figure was assembled from `B_b4_compare.csv` inside the main figure
script. Anyone auditing the deconvolution path should therefore read the Fig S2 block of
`R/12_figures.R` alongside the MuSiC script.

---

## Figure parameters worth knowing

| Item | Value |
|---|---|
| Heatmap clustering | **None.** Every `pheatmap` call sets `cluster_rows = FALSE, cluster_cols = FALSE`; row and column orders are hand-specified. No distance metric or linkage is used anywhere. |
| Fig 3 / Fig 1 z-scoring | Within dataset, across cell types (Fig 3); per gene across cohorts (Fig 1) |
| Fig 1 colour limits | Clipped at the **98th percentile** of `abs(effect)` — values above are saturated, so the legend shows ±mx, not the full range |
| Fig 3 constant rows | Set to 0 rather than dropped |
| Missing genes | Rendered as `NA` in grey — **never imputed, never dropped** |
| Fig 5 y value | `log2(CPM + 1)` of patient-level pseudobulk |
| Fig 2 panel B value | `log2(FPKM + 1)` mean difference, CAG vs CNAG |
| GSE153224 grouping | From a **filename convention** — FPKM columns beginning `W-` are CAG |
| Star thresholds | `Q < 0.001 ***`, `< 0.01 **`, `< 0.05 *`; `NA → "n.d."` |
| BH family | 140 = 14 genes × 10 cohorts, over the flattened matrix |
| Fig S2 attenuation | `1 − adj/raw`; cells with `att > 1` or `verdict == "direction_flipped"` are **dropped from the plot** because they are sign flips, not attenuations |
| Jitter seed | `position_jitter(..., seed = 42)` in the ggplot fibrosis panel. The base-R `stripchart(method = "jitter")` has **no seed** — its point positions are not reproducible. |
| Output format | 300 dpi PNG + `cairo_pdf` for font embedding |

---

## Figure files in this repository

`results/figures/` holds the figure files that were **submitted with the manuscript**.
The scripts regenerate them into that directory when the upstream result tables are
present.

Note that three of the seven figures were renamed by hand between the script's output
name and the submitted name, and one figure (Fig 5) came from a second script. The
mapping above is derived from the manuscript figure legends, the project handover note,
and byte-size/timestamp identity between the on-disk script outputs and the submitted
files — not from any rename call in the code, because none exists.

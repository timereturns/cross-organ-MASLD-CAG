# Data sources

No expression data are redistributed in this repository. Every dataset used in the
study is publicly available from the accession listed below and must be downloaded
by the user. GEO series matrices are not mirrored here: several of them are tens of
gigabytes, and the GEO terms of use do not permit wholesale re-hosting.

Corresponding author for the manuscript: Yanzhu Chen, Binfang Zeng (Xinjiang Medical
University, School of Traditional Chinese Medicine).

---

## 1. Liver discovery cohorts (n = 6)

| Accession | First author / study | Platform | Role | Downloaded |
|---|---|---|---|---|
| GSE48452 | Ahrens et al. | Illumina HumanHT-12 V4.0 | Liver discovery | _(fill in)_ |
| GSE89632 | Arendt et al. | Affymetrix Human Genome U219 | Liver discovery | _(fill in)_ |
| GSE83452 | Lefebvre et al. | Illumina HumanHT-12 V4.0 | Liver discovery | _(fill in)_ |
| GSE126848 | Suppli et al. | Illumina HiSeq 2500 (RNA-seq) | Liver discovery | _(fill in)_ |
| GSE130970 | Hoang et al. | Illumina HiSeq 2500 (RNA-seq) | Liver discovery | _(fill in)_ |
| GSE162694 | Pantano et al. | Illumina NovaSeq 6000 (RNA-seq) | Liver discovery | _(fill in)_ |

## 2. Gastric discovery cohorts (n = 2)

| Accession | First author / study | Platform | Role | Downloaded |
|---|---|---|---|---|
| GSE116312 | — | — | CAG vs follicular gastritis | _(fill in)_ |
| GSE27411 | Nookaew et al. | Affymetrix Human Genome U133 Plus 2.0 | Atrophy vs non-infected stomach | _(fill in)_ |

## 3. External validation cohorts (n = 2, held out from discovery)

| Accession | Samples | Platform | Role | Downloaded |
|---|---|---|---|---|
| GSE135251 | n = 216 | Illumina (RNA-seq) | Liver validation; NAS and fibrosis staging available | _(fill in)_ |
| GSE153224 | 10 pooled libraries, 40 patients | Illumina (RNA-seq) | Gastric validation | _(fill in)_ |

These two cohorts were **held out of the discovery screen**: they were not used to
select the 14 candidate genes, and were interrogated once for replication only.

## 3b. Additional bulk cohort used in the deconvolution sensitivity analysis

| Accession | Samples | Role | Downloaded |
|---|---|---|---|
| GSE174478 | 94 | 4th bulk cohort for the MuSiC deconvolution sensitivity analysis (Tables S32–S36, Fig S2). **Not** part of the discovery set. Its endpoint is a numeric fibrosis stage (0–4), not case/control, so comparison arms 1 and 2 are not available for it. | _(fill in)_ |

⚠️ **GSE174478 is not listed in the manuscript's Data availability statement.** It is
analysed and reported in the supplementary tables and Figure S2, so it needs to be added
there. See `docs/manuscript_discrepancies.md` §19.


## 4. Single-cell / single-nucleus datasets (n = 3)

| Accession | Tissue | Samples | Role | Downloaded |
|---|---|---|---|---|
| GSE202379 | Liver | 47 | snRNA-seq cell-type localisation; fibrosis-stage analysis (Fig 5) | _(fill in)_ |
| GSE115469 | Liver | — | Cell-type confirmation / annotation reference | _(fill in)_ |
| GSE134520 | Stomach | 9 | Gastric snRNA-seq patient-level pseudobulk (Fig S1) | _(fill in)_ |

## 5. Non-GEO sources

| Source | Version | Used for | Accessed |
|---|---|---|---|
| eQTLGen consortium | https://www.eqtlgen.org | Cis-eQTL instruments for the exploratory MR | _(fill in)_ |
| GTEx Portal | v10 (https://gtexportal.org) | Additional eQTL instruments | _(fill in)_ |
| OpenGWAS / MR-Base | https://gwas.mrcieu.ac.uk | Outcome GWAS summary statistics | _(fill in)_ |
| Human Protein Atlas | version 25.0 (https://www.proteinatlas.org) | Normal-tissue protein and cell-type annotation check | 2026-09-05 |

---

## How to download

```r
# GEO series matrix (per accession)
#   https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE135251/matrix/
# Or from within R (requires GEOquery, not needed for the committed scripts):
#   GEOquery::getGEO("GSE135251", destdir = "data/GSE135251")
```

Place each accession in its own folder under `data/`, e.g. `data/GSE135251/`.
Every `R/*.R` script resolves paths relative to the repository root (see
`R/00_setup.R`) and will tell you if a dataset folder is missing.

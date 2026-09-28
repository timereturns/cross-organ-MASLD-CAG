# =============================================================================
# 08_singlecell_pseudobulk.R — cell-type localisation in three sc/snRNA datasets
# =============================================================================
# STATUS: RESTORATION (from the file the author called 跑 compact 版三件套.R).
#
# Three datasets, three jobs:
#
#   GSE202379  liver snRNA, 47 patients -- the primary localisation, and the source
#              of the fibrosis-stage analysis in Fig 5.
#   GSE115469  liver scRNA -- an independent liver localisation used to check that
#              the cell-type assignment is not an artefact of one dataset. It also
#              provides the lineage vocabulary GSE202379 is mapped onto.
#   GSE134520  gastric snRNA, 9 patients -- the gastric side, descriptive (Fig S1).
#
# Patient-level pseudobulk throughout: counts are summed per (patient x cell type),
# then CPM-normalised against that unit's library size. Nothing is analysed
# per-cell, because a per-cell test on 100,000 cells would report significance for
# differences far smaller than any biological effect.
#
# THE 20-CELL RULE. A patient x cell-type unit enters analysis only if it has at
# least MIN_CELLS_PER_PATIENT_CELLTYPE cells. The value is 20 everywhere it appears
# in the original, but it appears as an inline literal in several places rather
# than as one parameter. It is declared once in R/00_setup.R here.
#
# FIBROSIS STAGE. Patient stage is the MODE of the per-cell annotation. A rounded
# mean was computed as a comparator during the original work and rejected because
# it did not reproduce the three locked correlation values. Ties in the mode
# resolve to the lower stage, because which.max returns the first maximum and the
# table columns are ordered 0..4. That is a real edge case: 9 of 47 patients sit at
# F0 and the stage distribution is 5/9/12/12/9.
#
# Run:  Rscript R/08_singlecell_pseudobulk.R
# =============================================================================

source("R/00_setup.R")

# -----------------------------------------------------------------------------
# 1. GSE202379 — liver snRNA, primary localisation
# -----------------------------------------------------------------------------
# The compact object holds, for the 13 candidate genes present in this dataset,
# the expression matrix, the raw counts and the per-cell metadata. It is produced
# by the single-cell preprocessing that is not included in this repository (see
# docs/repository_scope.md); this script consumes it.

gse202379_localize <- function() {
  require_dataset("GSE202379")
  f <- file.path(cache_dir, "p1_gse202379_compact.rds")
  if (!file.exists(f)) {
    stop("p1_gse202379_compact.rds not found in cache/.\n",
         "  It is produced by the single-cell preprocessing, which is not part ",
         "of this repository.\n",
         "  See docs/repository_scope.md section 3.", call. = FALSE)
  }
  cc <- readRDS(f)
  md <- cc$meta
  genes <- cc$genes_present

  # NOTE: the original named this vector genes14 while it holds 13 entries,
  # because one candidate is absent from this dataset. The comment sits in the
  # original source. Reported, not silently renamed.
  message("[08] GSE202379: ", length(genes), " of ", length(candidates_14),
          " candidates present: ", paste(genes, collapse = ", "))
  message("[08]   absent: ",
          paste(setdiff(candidates_14, genes), collapse = ", "))

  ct <- factor(md$cell.annotation)

  # --- localisation: mean expression and fraction expressing, per cell type ----
  E <- cc$expr14
  loc <- do.call(rbind, lapply(genes, function(g) {
    x <- as.numeric(E[g, ])
    m <- tapply(x, ct, mean)
    p <- tapply(x > 0, ct, mean)
    data.frame(gene = g, celltype = names(m),
               mean_log2expr = round(m, 3), pct_expr = round(p, 3),
               stringsAsFactors = FALSE)
  }))
  write_table(loc, "p1_gse202379_localization.csv")

  # --- patient-level pseudobulk -----------------------------------------------
  D <- cc$counts14
  grp <- interaction(md$Patient.ID, ct, drop = TRUE)
  pb  <- rowsum(t(D), grp)
  tot <- rowsum(md$nCount_RNA, grp)[, 1]
  cpm <- pb / as.numeric(tot) * 1e6
  pb_df <- data.frame(
    group = rownames(pb),
    patient = sub("\\..*$", "", rownames(pb)),
    celltype = sub("^[^.]*\\.", "", rownames(pb)),
    n_cells = as.integer(table(grp)[rownames(pb)]),
    round(cpm, 3),
    stringsAsFactors = FALSE, check.names = FALSE)
  write_table(pb_df, "p1_gse202379_pseudobulk_cpm.csv")

  message("[08]   patient-celltype units: ", nrow(pb_df),
          " | patients: ", length(unique(pb_df$patient)))
  message("[08]   units with >= ", MIN_CELLS_PER_PATIENT_CELLTYPE, " cells: ",
          sum(pb_df$n_cells >= MIN_CELLS_PER_PATIENT_CELLTYPE))

  # --- fibrosis stage per patient: MODE ---------------------------------------
  fib_tab <- table(md$Patient.ID, md$Fibrosis.score..F0.4.)
  fib_num <- suppressWarnings(as.numeric(
    apply(fib_tab, 1, function(r) names(which.max(r)))))

  # Report the comparator too, so the choice between them is not invisible.
  fib_mean <- suppressWarnings(round(as.numeric(
    tapply(as.numeric(md$Fibrosis.score..F0.4.), md$Patient.ID, mean))))
  n_agree <- sum(fib_num == fib_mean, na.rm = TRUE)
  message("[08]   stage: mode and rounded mean agree for ", n_agree, "/",
          length(fib_num), " patients")
  message("[08]   stage distribution (mode): ",
          paste(names(table(fib_num)), table(fib_num), collapse = " / "))

  # --- association between candidate expression and stage ---------------------
  f <- fib_num[match(pb_df$patient, names(fib_num))]
  sel <- pb_df$n_cells >= MIN_CELLS_PER_PATIENT_CELLTYPE

  cor_res <- do.call(rbind, lapply(unique(pb_df$celltype), function(ctn) {
    s <- pb_df$celltype == ctn & sel
    if (sum(s) < 5) return(NULL)
    do.call(rbind, lapply(genes, function(g) {
      v <- pb_df[[g]][s]; ff <- f[s]
      ok <- !is.na(ff) & !is.na(v) & is.finite(v)
      if (sum(ok) < 5) return(NULL)
      r <- cor.test(v[ok], ff[ok], method = "spearman")
      data.frame(celltype = ctn, gene = g, n_patients = sum(ok),
                 rho = round(unname(r$estimate), 3),
                 p = signif(r$p.value, 3), stringsAsFactors = FALSE)
    }))
  }))
  write_table(cor_res, "p1_gse202379_fibrosis_cor.csv")

  # --- the three locked values -------------------------------------------------
  # Fig 5 rests on three gene x cell-type pairs. If they do not reproduce, the
  # stage definition or the pseudobulk pipeline has changed and the figure is no
  # longer the published one.
  locked <- list(c("IL32", "Hepatocytes", 0.639),
                 c("CADM2", "Stellate", -0.761),
                 c("ANXA4", "Cholangiocytes", -0.585))
  message("[08]   locked fibrosis correlations:")
  for (lk in locked) {
    hit <- cor_res[cor_res$gene == lk[1] & cor_res$celltype == lk[2], ]
    if (!nrow(hit)) {
      warning("[08] locked pair not found: ", lk[1], " / ", lk[2], call. = FALSE)
      next
    }
    target <- as.numeric(lk[3])
    ok <- abs(hit$rho[1] - target) < 0.01
    message(sprintf("[08]     %-20s n=%2d  rho=%+.3f  (locked %+.3f)  %s",
                    paste(lk[1], lk[2], sep = " / "), hit$n_patients[1],
                    hit$rho[1], target, if (ok) "OK" else "MISMATCH"))
    if (!ok) warning("[08] locked value not reproduced for ", lk[1], " / ",
                     lk[2], call. = FALSE)
  }

  list(loc = loc, pb = pb_df, cor = cor_res, stage = fib_num)
}

# -----------------------------------------------------------------------------
# 2. GSE115469 — liver scRNA, lineage confirmation
# -----------------------------------------------------------------------------
# Maps this dataset's finer cell-type labels onto the same lineage vocabulary used
# for GSE202379, so the two localisations can be compared. The mapping is the
# original's, unchanged.

liver_lineage_map <- function(ct) {
  out <- rep("Other", length(ct))
  out[grepl("^Hepatocyte", ct)]                  <- "Hepatocytes"
  out[ct == "Cholangiocytes"]                    <- "Cholangiocytes"
  out[grepl("Macrophage", ct)]                   <- "Macrophages"
  out[ct == "Hepatic_Stellate_Cells"]            <- "Stellate"
  out[grepl("LSEC|Endothelial", ct)]             <- "Endothelial"
  out[grepl("T_Cells|NK-like", ct)]              <- "Lymphocytes"
  out[grepl("B_Cells|Plasma", ct)]               <- "B_cells"
  out[grepl("Erythroid", ct)]                    <- "Erythroid"
  out
}

gse115469_localize <- function() {
  dir <- file.path(data_dir, "GSE115469")
  if (!dir.exists(dir)) { message("[08] GSE115469 absent; skipping"); return(NULL) }

  ann_f <- file.path(dir, "GSE115469_CellClusterType.txt.gz")
  dat_f <- file.path(dir, "GSE115469_Data.csv.gz")
  if (!file.exists(ann_f) || !file.exists(dat_f)) {
    message("[08] GSE115469 files not found; skipping"); return(NULL)
  }

  ann <- read.delim(gzfile(ann_f), sep = "\t", stringsAsFactors = FALSE)
  message("[08] GSE115469 annotation: ", nrow(ann), " cells, ",
          length(unique(ann$Sample)), " samples")
  message("[08]   cell types: ",
          paste(names(sort(table(ann$CellType), decreasing = TRUE))[1:min(12, length(table(ann$CellType)))],
                collapse = ", "))

  # 84k cells; the full matrix fits in memory but is read once and subset
  M <- read.csv(gzfile(dat_f), row.names = 1, check.names = FALSE)
  hit <- candidates_14[candidates_14 %in% rownames(M)]
  message("[08]   candidates present: ", length(hit), "/14",
          "  (absent: ", paste(setdiff(candidates_14, rownames(M)), collapse = ", "), ")")

  ct <- factor(liver_lineage_map(ann$CellType[match(colnames(M), ann$CellName)]))
  E  <- as.matrix(M[hit, , drop = FALSE])

  loc <- do.call(rbind, lapply(hit, function(g) {
    x <- as.numeric(E[g, ])
    m <- tapply(x, ct, mean)
    p <- tapply(x > 0, ct, mean)
    data.frame(gene = g, celltype = names(m),
               mean_expr = round(m, 3), pct_expr = round(p, 3),
               stringsAsFactors = FALSE)
  }))
  write_table(loc, "p1_gse115469_localization.csv")

  grp <- interaction(ann$Sample[match(colnames(M), ann$CellName)], ct, drop = TRUE)
  pb  <- rowsum(t(E), grp)
  tot <- rowsum(colSums(M), grp)[, 1]
  cpm <- pb / as.numeric(tot) * 1e6
  pb_df <- data.frame(
    group = rownames(pb),
    patient = sub("\\..*$", "", rownames(pb)),
    celltype = sub("^[^.]*\\.", "", rownames(pb)),
    n_cells = as.integer(table(grp)[rownames(pb)]),
    round(cpm, 3), stringsAsFactors = FALSE, check.names = FALSE)
  write_table(pb_df, "p1_gse115469_pseudobulk_cpm.csv")
  message("[08]   patient-celltype units: ", nrow(pb_df))

  rm(M, E); gc()
  list(loc = loc, pb = pb_df)
}

# -----------------------------------------------------------------------------
# 3. GSE134520 — gastric snRNA, descriptive
# -----------------------------------------------------------------------------
# 9 patients, 13 libraries, from pre-processed per-library matrices. Cell types are
# assigned here by marker score over 12 gastric lineages, since the dataset does not
# ship an annotation.

mk_stomach <- list(
  G_epithelium    = c("EPCAM", "KRT18", "KRT8", "TACSTD2", "CLDN18", "MUC13"),
  Pit             = c("GKN1", "GKN2", "MUC5AC", "TFF1"),
  Neck            = c("MUC6", "TFF2", "PGA3", "PGC"),
  Parietal        = c("ATP4A", "ATP4B"),
  Chief           = c("LIPF", "PGA4", "PGC"),
  Enteroendocrine = c("CHGA", "CHGB", "GHRL", "GAST"),
  Fibroblast      = c("COL1A1", "COL3A1", "DCN", "LUM", "ACTA2", "PDGFRB"),
  Endothelial     = c("PECAM1", "VWF", "CDH5", "ENG"),
  Myeloid         = c("CD68", "CD163", "LYZ", "C1QA", "C1QB", "FCER1G", "CD14"),
  T_NK            = c("CD3D", "CD3E", "CD2", "NKG7", "GZMA", "IL7R", "TRAC"),
  B_plasma        = c("MS4A1", "CD79A", "CD79B", "MZB1", "IGHG1", "IGKC"),
  Mast            = c("TPSAB1", "TPSB2", "CPA3"))

# Biopsy -> patient. P4, P7 and P8 each contribute two libraries; EGC is a single
# library with no subject id. This mapping is the original's and is not derivable
# from the data, so it is stated explicitly.
gse134520_patient_map <- c(
  "GSM3954946" = "P1", "GSM3954947" = "P2", "GSM3954948" = "P9",
  "GSM3954949" = "P3", "GSM3954950" = "P4", "GSM3954951" = "P4",
  "GSM3954952" = "P5", "GSM3954953" = "P6", "GSM3954954" = "P7",
  "GSM3954955" = "P7", "GSM3954956" = "P8", "GSM3954957" = "P8",
  "GSM3954958" = "EGC")

# Lesion class per patient, from the series metadata.
gse134520_lesion <- c(P1 = "NAG", P2 = "NAG", P9 = "NAG", P3 = "CAG", P4 = "CAG",
                      P5 = "IMW", P6 = "IMW", P7 = "IMS", P8 = "IMS", EGC = "EGC")

gse134520_localize <- function() {
  dir <- file.path(data_dir, "GSE134520")
  raw_dir <- file.path(dir, "GSE134520_raw")
  if (!dir.exists(raw_dir)) {
    message("[08] GSE134520_raw absent; skipping"); return(NULL)
  }
  fl <- list.files(raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  message("[08] GSE134520 matrices: ", length(fl))

  lst <- lapply(fl, function(f)
    read.delim(gzfile(f), sep = "\t", row.names = 1, check.names = FALSE))
  names(lst) <- basename(fl)

  # common genes across libraries, so the cbind is well defined
  cmn <- Reduce(intersect, lapply(lst, rownames))
  message("[08]   common genes: ", length(cmn))
  allM <- do.call(cbind, lapply(lst, function(x) x[cmn, , drop = FALSE]))

  cell_meta <- data.frame(
    cell = unlist(lapply(names(lst), function(nm) colnames(lst[[nm]]))),
    biopsy = rep(sub("_processed_.*$", "", names(lst)), vapply(lst, ncol, integer(1))),
    stringsAsFactors = FALSE)
  cell_meta$patient <- gse134520_patient_map[cell_meta$biopsy]
  cell_meta$lesion  <- unname(gse134520_lesion[cell_meta$patient])

  # --- cell-type assignment by marker score -----------------------------------
  # z-score per gene across cells, then take the lineage with the highest mean.
  allg <- intersect(unlist(mk_stomach), rownames(allM))
  message("[08]   markers found: ", length(allg), "/",
          length(unique(unlist(mk_stomach))))
  z <- t(scale(t(allM[allg, , drop = FALSE])))
  sc <- sapply(names(mk_stomach), function(ctn) {
    g <- mk_stomach[[ctn]][mk_stomach[[ctn]] %in% rownames(allM)]
    if (length(g) < 2) rep(NA_real_, ncol(allM))
    else colMeans(z[g, , drop = FALSE], na.rm = TRUE)
  })
  cell_meta$celltype <- apply(sc, 1, function(r) {
    if (all(is.na(r))) "Unassigned" else names(which.max(r))
  })
  print(table(cell_meta$celltype, cell_meta$lesion))

  # --- localisation (log1p of counts, this dataset is count-valued) -----------
  cand <- candidates_14[candidates_14 %in% rownames(allM)]
  E <- allM[cand, , drop = FALSE]
  E_log <- log1p(E)
  ct <- factor(cell_meta$celltype)
  loc <- do.call(rbind, lapply(cand, function(g) {
    x <- as.numeric(E_log[g, ])
    m <- tapply(x, ct, mean)
    p <- tapply(E[g, ] > 0, ct, mean)
    data.frame(gene = g, celltype = names(m),
               mean_log1p = round(m, 3), pct_expr = round(p, 3),
               stringsAsFactors = FALSE)
  }))
  write_table(loc, "p1_gse134520_localization.csv")

  # --- patient-level pseudobulk -----------------------------------------------
  grp <- interaction(cell_meta$patient, ct, drop = TRUE)
  pb  <- rowsum(t(E), grp)
  tot <- rowsum(colSums(allM), grp)[, 1]
  cpm <- pb / as.numeric(tot) * 1e6
  pb_df <- data.frame(
    group = rownames(pb),
    patient = sub("\\..*$", "", rownames(pb)),
    celltype = sub("^[^.]*\\.", "", rownames(pb)),
    n_cells = as.integer(table(grp)[rownames(pb)]),
    round(cpm, 3), stringsAsFactors = FALSE, check.names = FALSE)
  pb_df$lesion <- unname(gse134520_lesion[pb_df$patient])
  pb_df$hp     <- pb_df$patient %in% c("P5", "P7")
  write_table(pb_df, "p1_gse134520_pseudobulk_cpm.csv")

  # --- descriptive CAG vs NAG, in the main candidates only ------------------
  # Descriptive, not a test of the hypothesis: 3 NAG patients and 2 CAG patients
  # after the 20-cell filter. The original required >= 2 patients per side.
  main <- c("IL32", "CDHR2", "LGALS3", "ANXA4", "CADM2", "RPS6KA1")
  sub <- pb_df[pb_df$lesion %in% c("NAG", "CAG") &
                 pb_df$n_cells >= MIN_CELLS_PER_PATIENT_CELLTYPE, ]
  res <- do.call(rbind, lapply(main, function(g) {
    if (!g %in% names(sub)) return(NULL)
    do.call(rbind, lapply(unique(sub$celltype), function(ctn) {
      s <- sub[sub$celltype == ctn, ]
      v_c <- s[[g]][s$lesion == "CAG"]; v_n <- s[[g]][s$lesion == "NAG"]
      if (length(v_c) < 2 || length(v_n) < 2) return(NULL)
      p <- tryCatch(wilcox.test(v_c, v_n)$p.value, error = function(e) NA_real_)
      data.frame(gene = g, celltype = ctn,
                 n_CAG = length(v_c), n_NAG = length(v_n),
                 mean_CAG = round(mean(v_c), 1), mean_NAG = round(mean(v_n), 1),
                 log2FC = round(log2((mean(v_c) + 1) / (mean(v_n) + 1)), 2),
                 wilcox_p = signif(p, 3), stringsAsFactors = FALSE)
    }))
  }))
  write_table(res, "p1_gse134520_cag_vs_nag.csv")

  rm(lst, allM, z, sc, E, E_log, pb); gc()
  list(loc = loc, pb = pb_df)
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
sc202379 <- gse202379_localize()
sc115469 <- gse115469_localize()
sc134520 <- gse134520_localize()

message("[08] done: ",
        paste(c(if (!is.null(sc202379)) "GSE202379",
                if (!is.null(sc115469)) "GSE115469",
                if (!is.null(sc134520)) "GSE134520"), collapse = ", "))

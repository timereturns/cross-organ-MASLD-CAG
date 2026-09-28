# =============================================================================
# 13_figures.R — manuscript figures
# =============================================================================
# STATUS: RESTORATION of the FINAL version of each figure.
#
# The original working file was ~1,950 lines across 18 RStudio blocks. Most of it
# was iteration: four superseded revisions of the heatmaps, two of the schematic,
# two of the external-validation panels, and one abandoned block that loaded its
# inputs and then drew nothing. Only the version that produced the submitted file
# is kept here.
#
# ⚠️ FIGURE NUMBERS. The original's filenames do NOT match the manuscript numbers,
# and the renames were done by hand outside R -- no file.rename() call exists
# anywhere in the original. The mapping, and the evidence for it, is in
# figure_map.md:
#
#     original filename                   manuscript figure
#     fig1_localization_heatmap.png   ->  Figure 3
#     fig3_cross_organ_schematic.pdf  ->  Figure 4
#     fig2_fibrosis_boxplots.pdf      ->  Figure 5 (but superseded, see below)
#     Fig1_cross_cohort_audit.pdf     ->  Figure 1
#     Fig2_external_validation.pdf    ->  Figure 2
#     figS1_cag_nag_bars.pdf          ->  Figure S1
#     figS_B4_attenuation.pdf         ->  Figure S2
#
# This script names its outputs by the MANUSCRIPT number, which is the opposite of
# the original. That is a deliberate change: a script whose "fig1" is Figure 3 is a
# trap for every future reader.
#
# Figure 5 is NOT here. The submitted Fig 5 came from a later redraw script that
# uses the authoritative patient stage table rather than the modal stage computed
# in the single-cell pipeline; see R/14_fig5_redraw.R.
#
# DESIGN CHOICES WORTH KNOWING
#   * No clustering anywhere. Every heatmap sets cluster_rows/cols = FALSE and the
#     row and column orders are fixed by hand. There is no distance metric or
#     linkage method to report.
#   * Figure 3 is z-scored WITHIN dataset, across cell types, never across
#     datasets. A gene absent from a dataset is drawn as grey and is never imputed.
#   * Figure 1 clips its colour scale at the 98th percentile of |effect|, so the
#     legend shows +/- mx, not the full range. Values above are saturated.
#   * Attenuation in Figure S2 is 1 - adj/raw. Cells with att > 1 are sign flips,
#     not attenuations, and are dropped from the plot rather than drawn at >100%.
#
# Run:  Rscript R/13_figures.R   (after 05, 07, 08, 09)
# =============================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(pheatmap)
  library(ggplot2)
})

T <- function(f) file.path(tables_dir, f)

# --- shared palette -----------------------------------------------------------
# Colour carries gene identity across every figure, so the panels can be read
# against each other. One palette, defined once.
role_cols <- c(Anchor = "#E41A1C", Epithelial = "#1B9E77",
               `Immune-epithelial` = "#377EB8", `Stromal(liver)` = "#984EA3",
               Watch = "#E69F00", Demoted = "grey70")

role_of <- function(gene) {
  if (gene %in% c("IL32")) "Anchor"
  else if (gene %in% c("CDHR2", "ANXA4", "LGALS3")) "Epithelial"
  else if (gene == "CADM2") "Stromal(liver)"
  else if (gene == candidate_observation) "Watch"
  else "Demoted"
}

# Genes sharing a role form a visual block in the heatmaps.
gene_order <- c(candidates_main, candidate_observation, candidates_downgraded)
stopifnot(setequal(gene_order, candidates_14))

heat_ramp <- colorRampPalette(c("#4393C3", "#F7F7F7", "#D6604D"))(100)

save_pheatmap <- function(ph, file, width, height, dpi = 300) {
  out <- file.path(figures_dir, file)
  if (grepl("\\.png$", file)) {
    grDevices::png(out, width = width, height = height, units = "in", res = dpi)
  } else {
    grDevices::cairo_pdf(out, width = width, height = height)
  }
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  grDevices::dev.off()
  message("[13] [out] ", out)
  invisible(out)
}

save_gg <- function(plot, file, width, height, dpi = 300) {
  out <- file.path(figures_dir, file)
  ggplot2::ggsave(out, plot = plot, width = width, height = height,
                  units = "in", dpi = dpi,
                  device = if (grepl("\\.pdf$", file)) grDevices::cairo_pdf else NULL)
  message("[13] [out] ", out)
  invisible(out)
}

star_of <- function(p) {
  ifelse(is.na(p), "n.d.",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# =============================================================================
# Figure 1 — cross-cohort consistency of the 14 candidates across 10 cohorts
# =============================================================================
# 14 genes x 10 cohorts of per-cohort effect size, with BH-adjusted significance
# stars. The BH family is the whole flattened matrix: 14 x 10 = 140 tests. The
# cohort order is discovery liver, held-out liver, discovery gastric, held-out
# gastric, so the held-out columns sit adjacent to their own organ.

figure1 <- function() {
  eff <- read.csv(T("p0_cohort_effect_sizes.csv"), stringsAsFactors = FALSE)
  fin <- read.csv(T("candidate_genes_final.csv"), stringsAsFactors = FALSE)
  auc <- read.csv(T("shared_validation_auc.csv"), row.names = 1,
                  stringsAsFactors = FALSE)

  gsel <- gene_order
  liver <- c("Ahrens", "Arendt", "Lefebvre", "Suppli", "Hoang", "Pantano")
  gastric <- c("GSE116312", "GSE27411")

  # per-cohort effect matrix from the discovery cohorts
  E <- matrix(NA_real_, length(gsel), length(c(liver, gastric)),
              dimnames = list(gsel, c(liver, gastric)))
  P <- E
  for (g in gsel) {
    for (coh in colnames(E)) {
      row <- eff[eff$gene == g & eff$cohort == coh, ]
      if (nrow(row)) { E[g, coh] <- row$logFC[1]; P[g, coh] <- row$P[1] }
    }
  }

  # held-out columns. GSE135251 is a difference of group means on log2(CPM+1);
  # its p-value is the Wilcoxon from the validation table. GSE153224 is read from
  # the external validation results produced outside this repository.
  E$GSE135251 <- auc[gsel, "mean_case"] - auc[gsel, "mean_ctrl"]
  P$GSE135251 <- auc[gsel, "wilcox_p"]

  g153 <- T("GSE153224_candidate_effects.csv")
  if (file.exists(g153)) {
    e3 <- read.csv(g153, stringsAsFactors = FALSE)
    E$GSE153224 <- e3$log2FC[match(gsel, e3$gene)]
    P$GSE153224 <- e3$wilcox_p[match(gsel, e3$gene)]
  } else {
    message("[13] GSE153224 effects absent; Figure 1 will have one fewer column.\n",
            "      See docs/repository_scope.md - that script is not included.")
  }

  col_order <- c(liver[1:3], "GSE135251", liver[4:6], gastric, "GSE153224")
  col_order <- intersect(col_order, colnames(E))
  E <- E[, col_order, drop = FALSE]
  P <- P[, col_order, drop = FALSE]

  # --- BH over the whole matrix, 140 tests -----------------------------------
  Q <- matrix(p.adjust(as.vector(P), "BH"), nrow(E), ncol(E),
              dimnames = dimnames(P))
  message("[13] Fig 1: BH family = ", sum(!is.na(as.vector(P))), " tests (",
          nrow(E), " genes x ", ncol(E), " cohorts)")

  stars <- matrix(star_of(as.vector(Q)), nrow(Q), ncol(Q), dimnames = dimnames(Q))

  # --- annotations -----------------------------------------------------------
  exp_dir <- setNames(fin$direction, fin$gene)[gsel]
  dir_ok <- sign(E) == ifelse(exp_dir[rownames(E)] == "up", 1, -1)
  conc <- rowMeans(dir_ok, na.rm = TRUE)

  ann_row <- data.frame(
    role = factor(vapply(gsel, role_of, character(1)),
                  levels = c("Anchor", "Epithelial", "Immune-epithelial",
                             "Stromal(liver)", "Watch", "Demoted")),
    concordance = conc,
    row.names = gsel)

  ann_col <- data.frame(
    disease = factor(ifelse(colnames(E) %in% gastric, "CAG", "MASLD"),
                     levels = c("MASLD", "CAG")),
    stage = factor(ifelse(colnames(E) %in% validation_cohorts,
                          "held-out", "discovery"),
                   levels = c("discovery", "held-out")),
    row.names = colnames(E))

  ann_colors <- list(
    role = role_cols,
    concordance = colorRampPalette(c("#FFFFFF", "#08519C"))(100),
    disease = c(MASLD = "#2E7D32", CAG = "#E65100"),
    stage = c(discovery = "#525252", "held-out" = "#D9D9D9"))

  # 98th-percentile clipping: with 140 effects a few outliers would otherwise
  # flatten the whole colour range.
  mx <- ceiling(quantile(abs(E), 0.98, na.rm = TRUE) * 10) / 10
  breaks <- seq(-mx, mx, length.out = 101)
  message("[13] Fig 1 colour limits: +/-", mx, " (98th percentile clip)")

  # gaps so the role blocks and the organ groups are visually separated
  gaps_row <- cumsum(rle(as.character(ann_row$role))$lengths)
  gaps_row <- gaps_row[-length(gaps_row)]

  ph <- pheatmap(E,
                 color = heat_ramp,
                 breaks = breaks,
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 gaps_row = gaps_row,
                 gaps_col = c(3, 4, 6),   # liver | held-out liver | gastric ...
                 annotation_row = ann_row,
                 annotation_col = ann_col,
                 annotation_colors = ann_colors,
                 display_numbers = stars,
                 fontsize = 10, fontsize_row = 9.5, fontsize_col = 9,
                 fontsize_number = 7.5,
                 cellwidth = 18, cellheight = 16,
                 border_color = "white",
                 na_col = "grey90",
                 legend_breaks = seq(-mx, mx, length.out = 5),
                 legend_labels = round(seq(-mx, mx, length.out = 5), 1),
                 main = "Cross-cohort robustness audit: 14 candidates x 10 cohorts",
                 silent = TRUE)

  save_pheatmap(ph, "Fig1_cross_cohort_consistency.pdf", 10.5, 8)
  save_pheatmap(ph, "Fig1_cross_cohort_consistency.png", 10.5, 8)

  # the matrices behind the figure, so the numbers are traceable
  write_table(data.frame(gene = rownames(E), as.data.frame(E), check.names = FALSE),
              "Fig1_effect_matrix.csv")
  write_table(data.frame(gene = rownames(Q), as.data.frame(Q), check.names = FALSE),
              "Fig1_fdr_matrix.csv")
  write_table(data.frame(gene = rownames(stars), direction = exp_dir,
                         concordance = round(conc, 3),
                         as.data.frame(stars), check.names = FALSE),
              "Fig1_concordance.csv")
  invisible(ph)
}

# =============================================================================
# Figure 2 — replication in the two held-out cohorts
# =============================================================================
# Panel A GSE135251 (liver). Panel B GSE153224 (gastric). Horizontal bars of the
# five retained candidates, coloured by role, with arrows showing the expected
# direction.

figure2 <- function() {
  genes5 <- candidates_main
  exp_dir <- expected_direction[genes5]

  auc <- read.csv(T("shared_validation_auc.csv"), row.names = 1,
                  stringsAsFactors = FALSE)
  fcA <- auc[genes5, "mean_case"] - auc[genes5, "mean_ctrl"]
  pA  <- auc[genes5, "wilcox_p"]

  g153 <- T("GSE153224_candidate_effects.csv")
  if (!file.exists(g153)) {
    stop("Figure 2 panel B needs GSE153224_candidate_effects.csv, which comes ",
         "from the external validation script that is not included in this ",
         "repository. See docs/repository_scope.md.", call. = FALSE)
  }
  e3 <- read.csv(g153, stringsAsFactors = FALSE)
  fcB <- e3$log2FC[match(genes5, e3$gene)]
  pB  <- e3$wilcox_p[match(genes5, e3$gene)]

  panel <- function(fc, p, title, tag) {
    df <- data.frame(gene = factor(genes5, levels = rev(genes5)),
                     fc = fc, p = p,
                     role = vapply(genes5, role_of, character(1)),
                     dir = exp_dir, stringsAsFactors = FALSE)
    mx <- ceiling(max(abs(fc), na.rm = TRUE) * 10) / 10
    df$label <- ifelse(df$dir == "up", "\u2191 ", "\u2193 ")
    ggplot(df, aes(x = fc, y = gene, fill = role)) +
      geom_col(width = 0.62, linewidth = 0.5) +
      geom_text(aes(label = star_of(p), x = fc + sign(fc) * mx * 0.05),
                size = 3.3, vjust = 0.75) +
      scale_fill_manual(values = role_cols, guide = "none") +
      scale_y_discrete(labels = function(x) paste0(df$label[match(x, df$gene)], x)) +
      scale_x_continuous(expand = expansion(add = c(0.6, 1.15))) +
      labs(title = title, x = "log2 fold change", y = NULL, tag = tag) +
      theme_bw(base_size = 10.5) +
      theme(axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 11),
            plot.tag = element_text(size = 14, face = "bold"))
  }

  gA <- panel(fcA, pA, "GSE135251 - MASLD held-out (n = 216)", "A")
  gB <- panel(fcB, pB, "GSE153224 - CAG held-out (5 vs 5 libraries)", "B")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(gA, gB, nrow = 1)
  } else {
    message("[13] patchwork not installed; writing panels separately")
    combined <- gA
    save_gg(gB, "Fig2_panelB_external_validation.pdf", 5.5, 4.8)
  }

  save_gg(combined, "Fig2_external_validation.pdf", 11, 4.8)
  save_gg(combined, "Fig2_external_validation.png", 11, 4.8)

  write_table(data.frame(gene = genes5, role = vapply(genes5, role_of, character(1)),
                         GSE135251_log2FC = round(fcA, 3),
                         GSE135251_P = signif(pA, 3),
                         GSE153224_log2FC = round(fcB, 3),
                         GSE153224_P = signif(pB, 3),
                         both_direction_ok = sign(fcA) == sign(fcB)),
              "Fig2_external_validation_data.csv")
  invisible(combined)
}

# =============================================================================
# Figure 3 — cell-type localisation in three snRNA-seq datasets
# =============================================================================
# Two liver datasets (L1 GSE202379, L2 GSE115469) and one gastric (St GSE134520).
# Within-dataset z-scoring across cell types. Genes absent from a dataset are NA
# and drawn grey.

figure3 <- function() {
  gsel <- gene_order

  read_loc <- function(f, value_candidates) {
    if (!file.exists(f)) return(NULL)
    l <- read.csv(f, stringsAsFactors = FALSE)
    vcol <- intersect(value_candidates, names(l))[1]
    if (is.na(vcol)) vcol <- grep("^mean", names(l), value = TRUE)[1]
    l$v <- l[[vcol]]
    l
  }

  l1 <- read_loc(T("p1_gse202379_localization.csv"), "mean_log2expr")
  l2 <- read_loc(T("p1_gse115469_localization.csv"), "mean_expr")
  l3 <- read_loc(T("p1_gse134520_localization.csv"), "mean_log1p")
  sets <- list(L1 = l1, L2 = l2, St = l3)
  sets <- sets[!vapply(sets, is.null, logical(1))]
  if (!length(sets)) {
    stop("No localisation tables found. Run R/08 first.", call. = FALSE)
  }

  to_mat <- function(d, prefix) {
    m <- tapply(d$v, list(d$gene, d$celltype), function(x) x[1])
    # a cell type explicitly labelled unknown is not a cell type
    m <- m[, !grepl("unknown", colnames(m), ignore.case = TRUE), drop = FALSE]
    colnames(m) <- paste0(prefix, ".", colnames(m))
    m
  }

  mats <- lapply(names(sets), function(nm) to_mat(sets[[nm]], nm))
  big <- matrix(NA_real_, length(gsel), sum(vapply(mats, ncol, integer(1))),
                dimnames = list(gsel, unlist(lapply(mats, colnames))))
  for (m in mats) big[rownames(m), colnames(m)] <- m

  # z-score within dataset, across cell types. A constant row becomes 0 rather
  # than being dropped: dropping it would silently change the row count.
  zscore_rows <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  z <- matrix(NA_real_, nrow(big), ncol(big), dimnames = dimnames(big))
  for (nm in names(sets)) {
    cols <- grep(paste0("^", nm, "\\."), colnames(big))
    z[, cols] <- t(apply(big[, cols, drop = FALSE], 1, zscore_rows))
  }

  # dataset annotation, recovered from the column prefix
  ds <- sub("\\..*$", "", colnames(z))
  ann_col <- data.frame(
    dataset = factor(ds, levels = names(sets)),
    row.names = colnames(z))
  ann_colors <- list(dataset = c(L1 = "#2E7D32", L2 = "#81C784", St = "#E65100"))

  ann_row <- data.frame(
    role = factor(vapply(gsel, role_of, character(1)),
                  levels = c("Anchor", "Epithelial", "Immune-epithelial",
                             "Stromal(liver)", "Watch", "Demoted")),
    row.names = gsel)

  mx <- max(abs(z), na.rm = TRUE)
  lg <- ceiling(mx)

  ph <- pheatmap(z,
                 color = heat_ramp,
                 breaks = seq(-mx, mx, length.out = 101),
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 gaps_row = { g <- cumsum(rle(as.character(ann_row$role))$lengths); g[-length(g)] },
                 gaps_col = { n <- vapply(mats, ncol, integer(1)); cumsum(n)[-length(n)] },
                 annotation_row = ann_row,
                 annotation_col = ann_col,
                 annotation_colors = c(ann_colors,
                   list(role = role_cols)),
                 labels_col = sub("^(L1|L2|St)\\.", "", colnames(z)),
                 fontsize = 10, fontsize_row = 9.5, fontsize_col = 8.5,
                 angle_col = 45,
                 cellwidth = 15, cellheight = 17,
                 border_color = "white",
                 na_col = "grey90",
                 legend_breaks = c(-lg, -lg / 2, 0, lg / 2, lg),
                 main = "Candidate-gene cell-type localisation (within-dataset z-scored)",
                 silent = TRUE)

  save_pheatmap(ph, "Fig3_celltype_localization.pdf", 11, 6)
  save_pheatmap(ph, "Fig3_celltype_localization.png", 11, 6)

  message("[13] Fig 3: z-scored within dataset; ",
          sum(is.na(big)), " gene x cell-type cells are NA (gene absent from that ",
          "dataset) and are drawn grey, not imputed")
  invisible(ph)
}

# =============================================================================
# Figure 4 — cross-organ schematic
# =============================================================================
# Hand-composed. Every coordinate is a literal: this is a diagram, not a plot of
# data. The gastric CADM2 node carries the dashed treatment because its bulk
# signal could not be assigned to a cell type of origin.

figure4 <- function() {
  prog <- data.frame(
    program = c("Anchor", "Epithelial", "Immune-epithelial", "Stromal(liver)"),
    y = c(0.82, 0.62, 0.42, 0.22),
    fill = unname(role_cols[c("Anchor", "Epithelial", "Immune-epithelial",
                              "Stromal(liver)")]),
    stringsAsFactors = FALSE)
  prog$x <- 0.5

  gast <- data.frame(
    gene = c("IL32", "ANXA4", "LGALS3", "CDHR2", "CADM2"),
    level = c("High", "Medium", "High", "Medium", "Not detected"),
    resolved = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE)
  gast$x <- 0.15
  gast$y <- seq(0.8, 0.2, length.out = nrow(gast))

  liver <- data.frame(
    gene = c("IL32", "ANXA4", "CADM2", "LGALS3", "CDHR2"),
    level = c("Low", "Medium", "Not detected", "Low", "Medium"),
    resolved = TRUE, stringsAsFactors = FALSE)
  liver$x <- 0.85
  liver$y <- c(0.8, 0.65, 0.5, 0.35, 0.2)

  p <- ggplot() +
    annotate("rect", xmin = 0.02, xmax = 0.30, ymin = 0.10, ymax = 0.92,
             fill = "#E65100", alpha = 0.08) +
    annotate("rect", xmin = 0.70, xmax = 0.98, ymin = 0.10, ymax = 0.92,
             fill = "#2E7D32", alpha = 0.08) +
    annotate("text", x = 0.16, y = 0.96, label = "Gastric (CAG)",
             size = 4.2, fontface = "bold") +
    annotate("text", x = 0.84, y = 0.96, label = "Liver (MASLD)",
             size = 4.2, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.99, label = "SHARED CANDIDATE GENES",
             size = 3.4, fontface = "bold") +
    geom_rect(data = prog, aes(xmin = 0.44, xmax = 0.56,
                               ymin = y - 0.07, ymax = y + 0.07, fill = fill),
              alpha = 0.85, colour = NA) +
    geom_text(data = prog, aes(x = 0.5, y = y, label = program),
              size = 3.5, colour = "white", fontface = "bold")

  # curves from each shared program band to both organs
  for (side in c("gast", "liver")) {
    node <- if (side == "gast") gast else liver
    for (i in seq_len(nrow(prog))) {
      tgt <- node[i, ]
      if (is.na(tgt$x)) next
      x <- seq(if (side == "gast") 0.30 else 0.70,
               if (side == "gast") tgt$x + 0.09 else tgt$x - 0.09,
               length.out = 60)
      bend <- 0.55
      y0 <- prog$y[i]
      y1 <- min(tgt$y[1], y0)
      y <- y0 + (y1 - y0) * ((x - min(x)) / (max(x) - min(x)))^bend
      p <- p + annotate("path", x = x, y = y,
                        linewidth = 1.05, alpha = 0.8,
                        colour = prog$fill[i],
                        linetype = if (side == "gast" && !tgt$resolved[1]) 2 else 1)
    }
  }

  p <- p +
    geom_text(data = gast, aes(x = x, y = y, label = gene),
              size = 3.5, hjust = 1) +
    geom_text(data = liver, aes(x = x, y = y, label = gene),
              size = 3.5, hjust = 0) +
    annotate("text", x = 0.15, y = 0.13,
             label = "CADM2: bulk-detectable,\nlow abundance\n(cell of origin unresolved)",
             size = 2.4, lineheight = 0.95, hjust = 0.5) +
    labs(title = "Cross-organ molecular pattern",
         subtitle = "Solid curves: localisation reproduced across organs.  Dashed: signal present but unattributed.") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.07), clip = "off") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(size = 13.5, face = "bold"),
          plot.subtitle = element_text(size = 8.2),
          legend.position = "none")

  save_gg(p, "Fig4_cross_organ_schematic.pdf", 9, 6.8)
  save_gg(p, "Fig4_cross_organ_schematic.png", 9, 6.8)
  invisible(p)
}

# =============================================================================
# Figure S1 — GSE134520 gastric pseudobulk, descriptive
# =============================================================================
# Four candidates across five gastric cell types, CAG vs NAG. Descriptive: after
# the 20-cell filter there are 2 CAG and 3 NAG patients, so no test is reported.

figureS1 <- function() {
  f <- T("p1_gse134520_cag_vs_nag.csv")
  if (!file.exists(f)) {
    message("[13] Fig S1 input absent (run R/08); skipped"); return(invisible(NULL))
  }
  d <- read.csv(f, stringsAsFactors = FALSE)
  ct_order <- c("G_epithelium", "Pit", "Neck", "Chief", "Myeloid")
  d <- d[d$celltype %in% ct_order, , drop = FALSE]
  d$celltype <- factor(d$celltype, levels = ct_order)

  long <- rbind(
    data.frame(d[, c("gene", "celltype", "log2FC")],
               group = "CAG", expr = d$mean_CAG, stringsAsFactors = FALSE),
    data.frame(d[, c("gene", "celltype", "log2FC")],
               group = "NAG", expr = d$mean_NAG, stringsAsFactors = FALSE))

  p <- ggplot(long, aes(x = celltype, y = expr, fill = group)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    facet_wrap(~gene, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c(CAG = "#D6604D", NAG = "#4393C3")) +
    labs(title = "GSE134520 gastric pseudobulk by cell type",
         subtitle = "Patient-level pseudobulk CPM, group means; descriptive only",
         x = NULL, y = "CPM (group mean)", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(size = 12),
          axis.text.x = element_text(size = 8.5, angle = 45, hjust = 1),
          plot.title = element_text(size = 13),
          plot.subtitle = element_text(size = 9))

  save_gg(p, "FigS1_gse134520_cag_nag.pdf", 10, 4.6)
  save_gg(p, "FigS1_gse134520_cag_nag.png", 10, 4.6)
  invisible(p)
}

# =============================================================================
# Figure S2 — MuSiC-adjusted effect attenuation
# =============================================================================
# Per cohort, how much of the candidate's raw effect survives adjusting for the
# deconvolved cell-type proportions. The source table carries two curator
# annotation columns (verdict, note) that were appended by hand in the original;
# if they are absent the verdict is recomputed here from the computed columns.
#
# Cells with att > 1 are sign flips, not attenuations, and are dropped. This is why
# some gene-cohort cells have no bar.

figureS2 <- function() {
  f <- T("B_b4_compare.csv")
  if (!file.exists(f)) {
    message("[13] Fig S2 input absent (run R/09); skipped"); return(invisible(NULL))
  }
  b4 <- read.csv(f, stringsAsFactors = FALSE)

  # verdict: use the curator column if present, otherwise derive it
  if (!"verdict" %in% names(b4)) {
    message("[13] Fig S2: 'verdict' column absent; deriving it from the computed ",
            "columns. The published table has a curator-appended version.")
    b4$verdict <- with(b4, ifelse(abs(arm3_logFC_raw) < 0.1 & arm3_p_raw > 0.05,
                                  "no_signal",
                            ifelse(dir_keep_music & abs(att_music) < 0.5,
                                   "direction_kept",
                            ifelse(dir_keep_music,
                                   "direction_kept_attenuated",
                                   "direction_flipped"))))
  }

  adj_p_col <- if ("arm3_p_adj" %in% names(b4)) "arm3_p_adj" else
    grep("p.*adj", names(b4), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(adj_p_col)) stop("No adjusted-p column found in B_b4_compare.csv",
                             call. = FALSE)

  df <- b4
  df$att <- df$att_music
  # drop sign flips and out-of-range attenuations
  df$att[abs(df$att) > 2 | df$verdict == "direction_flipped" | df$att > 1] <- NA

  df$panel <- df$gene
  df$label <- ifelse(is.na(df$att), "",
                     ifelse(df$verdict == "no_signal", "no signal",
                            paste0(round(100 * df$att), "%")))

  p <- ggplot(df, aes(x = cohort, y = att)) +
    geom_col(fill = "#377EB8", width = 0.62) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
    annotate("text", x = 1, y = 0.52, label = "50% attenuation",
             size = 2.5, hjust = 0, colour = "grey40") +
    geom_text(aes(label = label), vjust = -0.4, size = 2.5) +
    facet_wrap(~panel, nrow = 1) +
    scale_x_discrete(labels = function(x)
      ifelse(x == "GSE174478", paste0(x, " \u2020"), x)) +
    scale_y_continuous(limits = c(0, 1.05),
                       labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "Effect attenuation after MuSiC cell-composition adjustment",
         subtitle = paste0("1 - adjusted/raw log2FC.  \u2020 = fibrosis-gradient arm ",
                           "(GSE174478).  Cells that flip sign are not shown."),
         x = NULL, y = "Attenuation") +
    theme_bw(base_size = 10.5) +
    theme(strip.text = element_text(size = 11),
          axis.text.x = element_text(size = 8.5, angle = 45, hjust = 1),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8.2))

  save_gg(p, "FigS2_music_attenuation.pdf", 9, 4.2)
  save_gg(p, "FigS2_music_attenuation.png", 9, 4.2)

  n_drop <- sum(is.na(df$att) & !is.na(df$att_music))
  message("[13] Fig S2: ", n_drop, " gene-cohort cells dropped as sign flips")
  invisible(p)
}

# =============================================================================
# Run
# =============================================================================
figs <- list(figure1 = figure1, figure2 = figure2, figure3 = figure3,
             figure4 = figure4, figureS1 = figureS1, figureS2 = figureS2)

ok <- character(0); failed <- character(0)
for (nm in names(figs)) {
  message("\n[13] ===== ", nm, " =====")
  r <- tryCatch({ figs[[nm]](); TRUE },
                error = function(e) { message("[13] FAILED: ", conditionMessage(e)); FALSE })
  if (r) ok <- c(ok, nm) else failed <- c(failed, nm)
}

message("\n[13] produced: ", if (length(ok)) paste(ok, collapse = ", ") else "none")
if (length(failed)) {
  message("[13] failed:   ", paste(failed, collapse = ", "))
  message("[13]   Figure 2 panel B needs the GSE153224 external validation results,",
          " which are not in this repository.")
  message("[13]   Figure 5 is produced by R/14_fig5_redraw.R, not here.")
}

# =============================================================================
# 14_fig5_redraw.R — Figure 5: candidate expression across fibrosis stages
# =============================================================================
# STATUS: RESTORATION (from 预审稿后图片重绘).
#
# WHY FIGURE 5 HAS ITS OWN SCRIPT, AND WHY THAT MATTERS
#
# The main figure script draws a fibrosis-stage panel too, but the submitted
# Figure 5 is NOT that panel. It was redrawn later, by this script, for three
# reasons that are worth stating rather than hiding:
#
#   1. The stage source changed. The redraw reads the authoritative patient-level
#      covariate table (GSE202379_患者级协变量与纤维化分期表.csv) instead of the
#      modal per-cell stage computed in the single-cell pipeline. The two were
#      compared and the authoritative table was adopted because it reproduced the
#      three locked correlation values. See manuscript_discrepancies.md.
#
#   2. The three headline correlations are asserted, not just computed. If
#      IL32/Hepatocytes, CADM2/Stellate and ANXA4/Cholangiocytes do not reproduce
#      their published rho values, the script stops. That is the strongest check in
#      the repository and it exists precisely because a figure was being redrawn
#      after the numbers had already been reported.
#
#   3. Fonts are embedded. The original switched to cairo_pdf here for submission
#      safety.
#
# THE THREE LOCKED VALUES
#   IL32  / Hepatocytes    rho = +0.639
#   CADM2 / Stellate       rho = -0.761
#   ANXA4 / Cholangiocytes rho = -0.585
# along with the stage distribution F0-F4 = 5 / 9 / 12 / 12 / 9 and the three
# per-pair patient counts (47 / 37 / 30).
#
# PATIENT-ID NORMALISATION. The pseudobulk table uses bare IDs (30, CL16, HL1);
# the authoritative table uses a "P" prefix (P30, PCL16, PHL1). The redraw strips
# the prefix on the reference side and matches on the stripped form. This is
# brittle but it is what the original did, and the assertions catch it if it
# breaks.
#
# Run:  Rscript R/14_fig5_redraw.R   (after 08)
# =============================================================================

source("R/00_setup.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# --- locked values ------------------------------------------------------------
LOCKED_RHO <- list(c("IL32", "Hepatocytes", 0.639),
                   c("CADM2", "Stellate", -0.761),
                   c("ANXA4", "Cholangiocytes", -0.585))
LOCKED_STAGE_COUNTS <- c(`0` = 5L, `1` = 9L, `2` = 12L, `3` = 12L, `4` = 9L)
LOCKED_N <- c("IL32/Hepatocytes" = 47L, "CADM2/Stellate" = 37L,
              "ANXA4/Cholangiocytes" = 30L)

PANELS <- list(c("IL32",  "Hepatocytes"),
               c("ANXA4", "Cholangiocytes"),
               c("CADM2", "Stellate"),
               c("IL32",  "Lymphocytes"),
               c("CDHR2", "Hepatocytes"),
               c("LGALS3", "Macrophages"))

# --- inputs -------------------------------------------------------------------
pb_path  <- file.path(tables_dir, "p1_gse202379_pseudobulk_cpm.csv")
cov_path <- file.path(data_dir, "GSE202379",
                      "GSE202379_患者级协变量与纤维化分期表.csv")

if (!file.exists(pb_path)) {
  stop("p1_gse202379_pseudobulk_cpm.csv missing. Run R/08 first.", call. = FALSE)
}
if (!file.exists(cov_path)) {
  stop("The authoritative patient stage table is missing:\n  ", cov_path, "\n",
       "  It is an input to this analysis and is not produced by any script in\n",
       "  this repository. See docs/repository_scope.md.", call. = FALSE)
}

pb  <- read.csv(pb_path, stringsAsFactors = FALSE)
ref <- read.csv(cov_path, stringsAsFactors = FALSE,
                colClasses = c(patient_id = "character"))

# --- stage per patient, from the authoritative table --------------------------
# Strip the "P" prefix so the two tables can be matched.
ref$pb_id <- sub("^P", "", ref$patient_id)

stage_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(ref$stage))))
names(stage_num) <- ref$pb_id
message("[14] authoritative stage table: ", nrow(ref), " patients, ",
        sum(!is.na(stage_num)), " with a numeric stage")
stopifnot(sum(!is.na(stage_num)) >= 37L)

# --- n_cells per patient x cell type -----------------------------------------
# Rebuilt from the pseudobulk table's own cell counts rather than re-reading the
# single-cell object.
cell_counts <- pb[, c("patient", "celltype", "n_cells")]
patients <- sort(unique(pb$patient))

stage <- stage_num[patients]
if (anyNA(stage)) {
  message("[14] patients without a numeric stage: ",
          paste(patients[is.na(stage)], collapse = ", "))
}

# --- assertion A: stage distribution -----------------------------------------
obs_stage <- table(stage, useNA = "ifany")
message("[14] stage distribution: ",
        paste(names(obs_stage), obs_stage, collapse = " / "))
exp_vec <- LOCKED_STAGE_COUNTS
got_vec <- as.integer(obs_stage[names(exp_vec)])
if (any(is.na(got_vec)) || any(got_vec != as.integer(exp_vec))) {
  stop("Stage distribution does not match the locked values F0-F4 = ",
       paste(exp_vec, collapse = "/"), "\n  observed: ",
       paste(names(obs_stage), obs_stage, collapse = " / "), call. = FALSE)
}
message("[14] stage distribution matches the locked values")

# --- correlation for one gene x cell type ------------------------------------
cor_pair <- function(gene, celltype) {
  key  <- paste(gene, celltype, sep = "_")
  if (!key %in% names(pb)) return(NULL)

  s <- pb[pb$celltype == celltype, , drop = FALSE]
  s <- s[s$n_cells >= MIN_CELLS_PER_PATIENT_CELLTYPE, , drop = FALSE]
  s$stage <- stage[s$patient]
  s <- s[!is.na(s$stage) & is.finite(s[[key]]), , drop = FALSE]
  if (nrow(s) < 5) return(NULL)

  # exact = FALSE: the t asymptotic approximation, matching the recomputed
  # correlation table. With these n, the exact permutation p-value is
  # computationally heavy and the two differ slightly.
  r <- suppressWarnings(cor.test(s[[key]], s$stage, method = "spearman",
                                 exact = FALSE))
  list(n = nrow(s), rho = unname(r$estimate), p = r$p.value, data = s, key = key)
}

# --- assertion B: the three locked correlations -------------------------------
message("[14] ==== locked correlations ====")
locked_ok <- TRUE
for (lk in LOCKED_RHO) {
  g <- lk[1]; ct <- lk[2]; target <- as.numeric(lk[3])
  r <- cor_pair(g, ct)
  if (is.null(r)) {
    message("[14]   ", g, "/", ct, ": not estimable (too few patients)")
    locked_ok <- FALSE; next
  }
  ok <- abs(r$rho - target) < 0.01
  message(sprintf("[14]   %-22s n=%2d  rho=%+.3f  (locked %+.3f)  p=%.2e  %s",
                  paste(g, ct, sep = " / "), r$n, r$rho, target, r$p,
                  if (ok) "OK" else "MISMATCH"))
  if (!ok) locked_ok <- FALSE

  # assertion C: patient count
  nm <- paste(g, ct, sep = "/")
  if (nm %in% names(LOCKED_N) && r$n != LOCKED_N[[nm]]) {
    message("[14]     note: n = ", r$n, ", published n = ", LOCKED_N[[nm]])
    locked_ok <- FALSE
  }
}

if (!locked_ok) {
  stop("The published correlation values were not reproduced.\n",
       "  The figure would not be the published figure. Investigate the stage\n",
       "  source or the pseudobulk construction before proceeding.", call. = FALSE)
}
message("[14] all three locked correlations reproduced")

# --- correlation table --------------------------------------------------------
cor_rows <- list()
for (p in PANELS) {
  r <- cor_pair(p[1], p[2])
  if (is.null(r)) next
  cor_rows[[length(cor_rows) + 1L]] <- data.frame(
    gene = p[1], celltype = p[2], n_patients = r$n,
    rho = round(r$rho, 3), p = signif(r$p, 3), stringsAsFactors = FALSE)
}
cor_tab <- do.call(rbind, cor_rows)
write_table(cor_tab, "p1_gse202379_fibrosis_cor_locked.csv")
print(cor_tab)

# --- the figure ---------------------------------------------------------------
# Panel order: the significant pairs first (up, then down), then the three
# non-significant ones. The order is fixed by hand so the figure reads as a
# narrative; it is not sorted by any statistic.
plot_df <- do.call(rbind, lapply(PANELS, function(p) {
  r <- cor_pair(p[1], p[2])
  if (is.null(r)) return(NULL)
  d <- r$data
  d$gene <- p[1]; d$celltype <- p[2]
  d$panel <- factor(paste(p[1], p[2], sep = " / "),
                    levels = vapply(PANELS, function(q) paste(q[1], q[2], sep = " / "),
                                    character(1)))
  d$rho <- r$rho; d$p <- r$p; d$expr <- d[[r$key]]
  d[, c("panel", "gene", "celltype", "patient", "stage", "expr", "rho", "p")]
}))
plot_df$stage_f <- factor(plot_df$stage, levels = 0:4)

labs <- unique(plot_df[, c("panel", "rho", "p")])
labs$label <- sprintf("rho = %+.3f\np = %.2g", labs$rho, labs$p)

p <- ggplot(plot_df, aes(x = stage_f, y = log2(expr + 1))) +
  geom_boxplot(outlier.shape = NA, width = 0.7, fill = "grey92") +
  geom_point(aes(colour = gene), position = position_jitter(width = 0.12,
                                                             height = 0, seed = JITTER_SEED),
             size = 1.6, alpha = 0.85) +
  geom_text(data = labs, aes(x = 1, y = Inf, label = label),
            hjust = 0, vjust = 1.2, size = 2.9, inherit.aes = FALSE) +
  facet_wrap(~panel, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = role_cols, guide = "none") +
  labs(x = "Fibrosis stage (F0-F4)", y = "log2(CPM + 1)",
       title = "Candidate-gene expression across fibrosis stages (GSE202379)",
       subtitle = paste0("Patient-level pseudobulk; patients with >= ",
                         MIN_CELLS_PER_PATIENT_CELLTYPE,
                         " cells of that type. Spearman rho with the asymptotic t approximation.")) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 9),
        plot.title = element_text(size = 12.5),
        plot.subtitle = element_text(size = 9))

out_pdf <- file.path(figures_dir, "Fig5_fibrosis_stages.pdf")
ggsave(out_pdf, p, width = 10, height = 6.2, units = "in",
       device = grDevices::cairo_pdf)   # cairo_pdf for embedded fonts
message("[14] [out] ", out_pdf)

out_png <- file.path(figures_dir, "Fig5_fibrosis_stages.png")
ggsave(out_png, p, width = 10, height = 6.2, units = "in", dpi = 300)
message("[14] [out] ", out_png)

# --- patients per stage for the analysis pairs (Table S45) --------------------
stage_counts <- function(gene, ct) {
  r <- cor_pair(gene, ct)
  if (is.null(r)) return(NULL)
  r$data |>
    dplyr::count(stage, name = "n_patients") |>
    dplyr::mutate(pair = paste(gene, ct, sep = " / "), .before = 1)
}
s45 <- dplyr::bind_rows(lapply(PANELS[1:3], function(p) stage_counts(p[1], p[2])))
write_table(s45, "p1_gse202379_fibrosis_stage_counts.csv")
print(s45)

message("[14] done")

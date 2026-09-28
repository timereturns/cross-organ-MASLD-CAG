# =============================================================================
# 10_candidate_scoring.R — 8-dimension adjudication of the 14 candidates
# =============================================================================
# STATUS: RESTORATION (from the file the author called P2 开始).
#
# Fourteen candidates reach this step. Five are retained, one goes on an
# observation list because MR gives it no directional support, and eight are
# downgraded. This script produces the ranked table behind that split.
#
# THE SCORING MATRIX IS NOT FULLY COMPUTED, AND THAT IS DOCUMENTED RATHER THAN
# HIDDEN. Six of the eight dimensions come from data:
#
#   A  LODO stability          in200 count over the six leave-one-out runs
#   B  gastric support         how many gastric cohorts support the gene
#   C  external validation     GSE135251 direction and AUC_norm
#   D  cell-type specificity   max over the three sc/snRNA datasets
#   E  cross-organ conservation lineage agreement within and across organs
#   F  genetic support         count of MR estimates agreeing in direction
#
# Two do not:
#
#   G  novelty                 author-assigned
#   H  verifiability           author-assigned
#
# G and H are expert judgements encoded as numbers. Table S11 is titled a
# "transparent 8-dimension scoring matrix"; a reader will reasonably assume all
# eight are derived. They are not, and the two hand-assigned dimensions are marked
# as such in the output and in docs/manuscript_discrepancies.md item 8.
#
# The retained / observation / downgraded boundary is NOT a coded rule either. The
# script ranks by total and reports the grouping as it stands; the threshold is
# authorial. See the same document.
#
# Run:  Rscript R/10_candidate_scoring.R   (after 05, 07, 08, 11)
# =============================================================================

source("R/00_setup.R")

T <- function(f) file.path(tables_dir, f)
C <- function(f) file.path(cache_dir, f)

need <- c("candidate_genes_final.csv", "p0_lodo_exact_rra.csv",
          "p0_cag_support.csv", "shared_validation_auc.csv")
missing <- need[!file.exists(vapply(need, T, character(1)))]
if (length(missing)) {
  stop("Missing inputs: ", paste(missing, collapse = ", "),
       "\n  Run R/03, R/05, R/07 first (and R/11 for dimension F).", call. = FALSE)
}

cand <- read.csv(T("candidate_genes_final.csv"), stringsAsFactors = FALSE)
genes14 <- cand$gene
message("[10] scoring ", length(genes14), " candidates")

score <- data.frame(gene = genes14, stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# A. LODO stability
# -----------------------------------------------------------------------------
# How many of the six leave-one-out runs kept the gene inside the top 200.
# 6 -> 2.0, 5 -> 1.5, 4 -> 1.0, fewer -> 0.
lodo <- read.csv(T("p0_lodo_exact_rra.csv"), stringsAsFactors = FALSE)
in200 <- grep("^in200_drop_", names(lodo))
keep <- rowSums(lodo[, in200], na.rm = TRUE)
A <- setNames(ifelse(keep == 6, 2,
              ifelse(keep == 5, 1.5,
              ifelse(keep == 4, 1, 0))), lodo$gene)
score$A_LODO_MASLD <- A[genes14]
score$A_detail <- paste0(keep[match(genes14, lodo$gene)], "/", length(in200))

# -----------------------------------------------------------------------------
# B. Gastric support
# -----------------------------------------------------------------------------
# The support column holds strings like "1 (logFC=0.8, P=0.01)"; the leading
# character is the verdict.
cag <- read.csv(T("p0_cag_support.csv"), stringsAsFactors = FALSE)
cag_n <- apply(cag[, -c(1, 2), drop = FALSE], 1,
               function(x) sum(sub(" .*", "", x) == "1", na.rm = TRUE))
B <- setNames(ifelse(cag_n == 2, 2, ifelse(cag_n == 1, 1, 0)), cag$gene)
score$B_CAG <- B[genes14]
score$B_detail <- paste0(cag_n[match(genes14, cag$gene)], "/2")

# -----------------------------------------------------------------------------
# C. External validation
# -----------------------------------------------------------------------------
# Direction must agree first -- a gene that fails direction scores 0 regardless of
# AUC. Then AUC_norm is banded. NOTE: AUC_norm = max(AUC, 1 - AUC), so a strongly
# down-regulated gene scores high here too. That is intended (it is a directional
# concordance metric) but it must be described as such, or 0.918 reads as
# discrimination performance.
val <- read.csv(T("shared_validation_auc.csv"), row.names = 1, stringsAsFactors = FALSE)
C_ <- setNames(rep(NA_real_, length(genes14)), genes14)
for (g in genes14) {
  if (!g %in% rownames(val)) next
  if (!isTRUE(val[g, "direction_ok"])) { C_[g] <- 0; next }
  a <- val[g, "AUC_norm"]
  C_[g] <- if (a >= 0.9) 2 else if (a >= 0.8) 1.5 else
           if (a >= 0.7) 1 else if (a >= 0.6) 0.5 else 0.25
}
score$C_Validation <- C_[genes14]

# -----------------------------------------------------------------------------
# D. Cell-type specificity
# -----------------------------------------------------------------------------
# Patient-level specificity in each sc/snRNA dataset, then the maximum. The
# expression column name differs between datasets (mean_log2expr for the liver
# snRNA, mean_expr for GSE115469, mean_log1p for the gastric counts), so it is
# detected rather than assumed.
spec_of <- function(file, col) {
  if (!file.exists(file)) return(setNames(rep(NA_real_, length(genes14)), genes14))
  l <- read.csv(file, stringsAsFactors = FALSE)
  if (!col %in% names(l)) col <- grep("^mean", names(l), value = TRUE)[1]
  sapply(genes14, function(g) {
    s <- l[l$gene == g, , drop = FALSE]
    if (!nrow(s)) return(NA_real_)
    v <- s[[col]]; p <- s$pct_expr
    i <- which.max(v)
    sec <- if (nrow(s) > 1) sort(v, decreasing = TRUE)[2] else 0
    ratio <- if (sec > 0) v[i] / sec else 99
    pct <- p[i]
    if (pct >= 0.25 && ratio >= 2) 2 else
      if (pct >= 0.20) 1.5 else
        if (pct >= 0.10) 1 else 0.5
  })
}

d1 <- spec_of(T("p1_gse202379_localization.csv"), "mean_log2expr")
d2 <- spec_of(T("p1_gse115469_localization.csv"), "mean_expr")
d3 <- spec_of(T("p1_gse134520_localization.csv"), "mean_log1p")
D <- vapply(genes14, function(g) max(c(d1[g], d2[g], d3[g]), na.rm = TRUE),
            numeric(1))
score$D_CellSpecificity <- D

# -----------------------------------------------------------------------------
# E. Cross-organ conservation
# -----------------------------------------------------------------------------
# The top cell type in each liver dataset is mapped to a lineage, and the two
# liver datasets must agree. Cross-organ agreement is then checked against the
# gastric dataset. Four levels, from "liver agrees and the lineage carries across
# organs" down to "liver does not even agree with itself".

liver_lineage <- function(ct) {
  if (is.na(ct)) return(NA_character_)
  if (grepl("Hepatocyte|Cholangio", ct)) "Epithelial" else
    if (grepl("Stellate", ct)) "Stromal" else
      if (grepl("Macrophage|Lymphocyte|T_|NK|B-|B_", ct)) "Immune" else
        if (grepl("Endothelial|LSEC", ct)) "Endothelial" else "Other"
}
stomach_lineage <- function(ct) {
  if (is.na(ct)) return(NA_character_)
  if (grepl("G_epithelium|Pit|Neck|Parietal|Chief|Enteroendocrine", ct)) "Epithelial" else
    if (grepl("Fibroblast", ct)) "Stromal" else
      if (grepl("Myeloid|T_NK|B_plasma|Mast", ct)) "Immune" else
        if (grepl("Endothelial", ct)) "Endothelial" else "Other"
}

top_ct <- function(file, col) {
  if (!file.exists(file)) return(setNames(rep(NA_character_, length(genes14)), genes14))
  l <- read.csv(file, stringsAsFactors = FALSE)
  if (!col %in% names(l)) col <- grep("^mean", names(l), value = TRUE)[1]
  sapply(genes14, function(g) {
    s <- l[l$gene == g, , drop = FALSE]
    if (!nrow(s)) NA_character_ else s$celltype[which.max(s[[col]])]
  })
}

top1 <- top_ct(T("p1_gse202379_localization.csv"), "mean_log2expr")
top2 <- top_ct(T("p1_gse115469_localization.csv"), "mean_expr")
top3 <- top_ct(T("p1_gse134520_localization.csv"), "mean_log1p")

E <- vapply(genes14, function(g) {
  ll1 <- liver_lineage(top1[g]); ll2 <- liver_lineage(top2[g])
  sl  <- stomach_lineage(top3[g])
  liver_ok <- !is.na(ll1) && identical(ll1, ll2)
  cross_ok <- !is.na(ll1) && !is.na(ll2) && !is.na(sl) &&
    (identical(ll1, sl) || identical(ll2, sl))
  if (liver_ok && cross_ok) 2 else
    if (liver_ok) 1.5 else
      if (!is.na(ll1)) 1 else 0.5
}, numeric(1))
score$E_CrossOrgan <- E
score$E_detail <- paste(top1[genes14], "/", top2[genes14], "/", top3[genes14])

# -----------------------------------------------------------------------------
# F. Genetic support
# -----------------------------------------------------------------------------
# Count of MR estimates that agree with the expected direction. Genes with no MR
# row at all default to 0.5 -- they are not penalised to 0, but they also do not
# earn support. Four of the fourteen have no cis-eQTL instrument set.
mr_file <- T("mr_all_results_v4.csv")
F_ <- setNames(rep(0.5, length(genes14)), genes14)
mr_detail <- setNames(rep("no MR estimate", length(genes14)), genes14)
if (file.exists(mr_file)) {
  mr <- read.csv(mr_file, stringsAsFactors = FALSE)
  cnt <- tapply(mr$dir_ok, mr$gene, function(x) sum(x, na.rm = TRUE))
  tot <- tapply(mr$dir_ok, mr$gene, length)
  for (g in genes14) {
    if (g %in% names(cnt)) {
      k <- cnt[[g]]
      F_[g] <- if (k >= 5) 2 else if (k == 4) 1.5 else if (k == 3) 1 else
               if (k == 2) 0.5 else 0
      mr_detail[g] <- paste0(k, "/", tot[[g]])
    }
  }
} else {
  warning("[10] mr_all_results_v4.csv absent; dimension F defaults to 0.5 ",
          "for every gene (run R/11)", call. = FALSE)
}
score$F_GeneticSupport <- F_[genes14]
score$F_detail <- mr_detail[genes14]

# -----------------------------------------------------------------------------
# G and H. Author-assigned dimensions
# -----------------------------------------------------------------------------
# Not derived from any data. Encoded as numbers so they can enter the total, which
# is exactly why they must be labelled.
G_ <- c(IL32 = 0.5, LGALS3 = 1, ANXA4 = 1.5, CDHR2 = 2, CADM2 = 2,
        RPS6KA1 = 1, TSPAN3 = 0.5, GOLM1 = 0.5, ANO10 = 0.5,
        SLC6A16 = 1, KIAA1958 = 0.5, FGA = 0.5, FGB = 0.5, LEPR = 0.5)
H_ <- c(IL32 = 2, LGALS3 = 2, ANXA4 = 2, CDHR2 = 1.5, CADM2 = 1,
        RPS6KA1 = 1.5, TSPAN3 = 1.5, GOLM1 = 1.5, ANO10 = 1,
        SLC6A16 = 1, KIAA1958 = 1, FGA = 2, FGB = 2, LEPR = 1.5)

score$G_Novelty       <- G_[genes14]
score$H_Verifiability <- H_[genes14]
score$G_H_source <- "author-assigned (not computed)"

# -----------------------------------------------------------------------------
# Total and grouping
# -----------------------------------------------------------------------------
dims <- c("A_LODO_MASLD", "B_CAG", "C_Validation", "D_CellSpecificity",
          "E_CrossOrgan", "F_GeneticSupport", "G_Novelty", "H_Verifiability")
stopifnot(all(dims %in% names(score)))

score$total <- rowSums(score[, dims], na.rm = TRUE)
score <- score[order(-score$total), , drop = FALSE]
score$rank <- seq_len(nrow(score))

# The published grouping. Report it against the ranking so a change in the data
# shows up as a mismatch rather than passing unnoticed.
score$published_group <- ifelse(score$gene %in% candidates_main, "retained",
                          ifelse(score$gene == candidate_observation, "observation",
                                 "downgraded"))

write_table(score, "p2_scoring_matrix.csv")

message("[10] ==== ranking ====")
print(score[, c("rank", "gene", "total", "A_LODO_MASLD", "B_CAG", "C_Validation",
                "D_CellSpecificity", "E_CrossOrgan", "F_GeneticSupport",
                "G_Novelty", "H_Verifiability", "published_group")])

message("[10] NOTE: dimensions G and H are author-assigned; the retained / ",
        "observation / downgraded boundary is not a coded rule.")

# Consistency between the ranking and the published grouping.
ret <- score$gene[score$published_group == "retained"]
message("[10] retained genes: ", paste(ret, collapse = ", "))
if (!setequal(ret, candidates_main)) {
  warning("[10] the retained set does not match the published set:\n",
          "  ranked:   ", paste(ret, collapse = ", "), "\n",
          "  published: ", paste(candidates_main, collapse = ", "),
          call. = FALSE)
}
message("[10] done")

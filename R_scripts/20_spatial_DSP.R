# ============================================================
# 20_spatial_DSP.R — Spatial co-localisation of the ferro-aging
#   score with cancer-associated fibroblasts (CAF) in intact CRC
#   tissue, and its compartment / tumour-vs-normal localisation.
#
#   Data source (downloaded, not simulated): GEO GSE281413
#   "Early onset and average onset colorectal cancer with matched
#   normal colon tissues" — GeoMx Digital Spatial Profiler (DSP),
#   Homo sapiens, 111 region-of-interest segments across 6 CRC
#   tumours and matched normal colon. Per-spot quantification
#   matrix  data/raw/GSE281413/GSE281413_processed_data.txt.gz
#   (14,364 targets x 111 spots); per-spot compartment / tumour
#   annotations parsed from the GEO series matrix into
#   data/raw/GSE281413/spot_annot.tsv (spot, segment, treatment,
#   cancer, roi). Each spot carries one WTA segment: PanCK
#   (epithelium), CD45 (immune) or Vimentin (stroma).
#
#   We test three questions, all pre-specified:
#     (1) co-localisation — is TIMP1 (the lead gene) spatially
#         correlated with a CAF marker program across spots and
#         specifically within the tumour stroma? (replicates the
#         bulk CAF-burden TIMP1~CAF rho=0.60 / scRNA 0.54 in intact tissue);
#     (2) compartment — where does the ferro-aging score sit
#         (stromal vs epithelial vs immune segments)?;
#     (3) context — is the score elevated in tumour vs matched
#         normal tissue?
#   Because GEO DSP normalisation units differ from bulk counts we
#   work on per-spot z-scores (scale within this dataset only) and
#   report associational, not causal, language. n = 111 spots from
#   6 donors, so estimates are noisy; we annotate n on every panel.
# ============================================================

# ---- 自动定位项目根目录（不写死本机路径）----
.script_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.script_file)) {
  .fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .script_file <- if (length(.fa)) normalizePath(sub("^--file=", "", .fa[1])) else file.path(getwd(), "x.R")
}
.ROOT <- dirname(dirname(.script_file))
while (!file.exists(file.path(.ROOT, "R_scripts", "00_config.R")) && !identical(.ROOT, dirname(.ROOT))) .ROOT <- dirname(.ROOT)
source(file.path(.ROOT, "R_scripts", "00_config.R"))
load_core_pkgs()
suppressPackageStartupMessages(library(ggpubr))

SPAT_DIR <- file.path(DATA_RAW, "GSE281413")

# ---- inputs --------------------------------------------------
mat <- read.delim(file.path(SPAT_DIR, "proc.txt"), row.names = "TargetName",
                  check.names = FALSE)                      # genes x spots
ann <- read.delim(file.path(SPAT_DIR, "spot_annot.tsv"), stringsAsFactors = FALSE)

# orient: keep only spots present in annotation, preserve order
ann <- ann[ann$spot %in% colnames(mat), ]
mat <- mat[, ann$spot, drop = FALSE]
message("Spots: ", ncol(mat), "  |  targets: ", nrow(mat))

# normalise annotation
ann$cancer_flag <- ifelse(tolower(ann$cancer) %in% c("yes"), "Tumour",
                    ifelse(tolower(ann$cancer) %in% c("no"), "Normal", NA))
ann$onset <- ifelse(grepl("^EO", ann$treatment), "Early-onset",
              ifelse(grepl("^AO", ann$treatment), "Average-onset", NA))
ann$segment <- factor(ann$segment, levels = c("PanCK", "CD45", "Vimentin"))

# ---- marker programs (only genes confirmed present in the panel) ----
CAF_MARKERS <- c("COL1A1","COL1A2","FAP","ACTA2","LUM","SPARC","PDGFRB","DCN","POSTN")
SIG5        <- c("CDKN2A","CHEK1","MMP1","MMP3","TIMP1")
stopifnot(all(CAF_MARKERS %in% rownames(mat)), all(SIG5 %in% rownames(mat)))

# per-gene z-score ACROSS SPOTS (dataset-internal scaling)
zmat <- t(scale(t(as.matrix(mat))))          # genes x spots, each gene mean0/sd1
zmat <- zmat[, !duplicated(colnames(zmat))]

sig_score   <- colMeans(zmat[SIG5, , drop = FALSE])
caf_score   <- colMeans(zmat[CAF_MARKERS, , drop = FALSE])
timpx       <- as.numeric(zmat["TIMP1", ])

spot <- data.frame(spot = colnames(zmat), ann,
                   score = as.numeric(sig_score), caf = as.numeric(caf_score),
                   TIMP1 = timpx, stringsAsFactors = FALSE)
spot$segment <- factor(spot$segment, levels = c("PanCK","CD45","Vimentin"))
message("Segment counts: ", paste(names(table(spot$segment)), table(spot$segment), sep="=", collapse="  "))
message("Tumour/normal:  ", paste(names(table(spot$cancer_flag)), table(spot$cancer_flag), sep="=", collapse="  "))

# ============================================================
# (1) Co-localisation: TIMP1 vs each CAF marker (Spearman)
# ============================================================
corr_df <- function(x, y, label, marker = NA_character_) {
  ok <- !is.na(x) & !is.na(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  data.frame(Scope = label, Marker = marker, rho = ct$estimate,
             p = ct$p.value, n = sum(ok), stringsAsFactors = FALSE)
}

# TIMP1 vs individual CAF markers, across all spots and within tumour
mk_rows <- function(subset_expr, label) {
  s <- spot[subset_expr, ]
  do.call(rbind, lapply(CAF_MARKERS, function(m){
    ct <- suppressWarnings(cor.test(s$TIMP1, as.numeric(zmat[m, match(s$spot, colnames(zmat))]),
                                    method = "spearman"))
    data.frame(Scope = label, Marker = m, rho = ct$estimate, p = ct$p.value,
               n = nrow(s), stringsAsFactors = FALSE)
  }))
}
tab14 <- rbind(mk_rows(rep(TRUE, nrow(spot)), "All spots"),
               mk_rows(!is.na(spot$cancer_flag) & spot$cancer_flag == "Tumour", "Tumour spots"),
               mk_rows(!is.na(spot$cancer_flag) & spot$cancer_flag == "Tumour" & spot$segment == "Vimentin", "Tumour stroma (Vimentin+)"))
# ---- Aggregate (programme-level) correlations -----------------
# These are the numbers quoted in Fig 14A / Results; they MUST be persisted in
# Table 14 so the manuscript is traceable (see AUDIT_BRIEF §provenance).
.tum   <- !is.na(spot$cancer_flag) & spot$cancer_flag == "Tumour"
.strom <- .tum & spot$segment == "Vimentin"
agg <- rbind(
  corr_df(spot$TIMP1,  spot$caf, "All spots",               "TIMP1 vs CAF score"),
  corr_df(spot$TIMP1[.tum],   spot$caf[.tum],   "Tumour spots",   "TIMP1 vs CAF score"),
  corr_df(spot$TIMP1[.strom], spot$caf[.strom], "Tumour stroma (Vimentin+)", "TIMP1 vs CAF score"),
  corr_df(spot$score,  spot$caf, "All spots",               "5-gene score vs CAF score"),
  corr_df(spot$score[.tum],  spot$caf[.tum],  "Tumour spots",   "5-gene score vs CAF score"),
  corr_df(spot$score[.strom], spot$caf[.strom], "Tumour stroma (Vimentin+)", "5-gene score vs CAF score")
)
tab14 <- rbind(tab14, agg)
save_table(tab14, "Table14_spatial_TIMP1_CAF_correlation")

# ============================================================
# (2) Compartment localisation  (Kruskal-Wallis + pairwise)
# ============================================================
kw_seg <- kruskal.test(score ~ segment, data = spot)
kw_seg_TIMP1 <- kruskal.test(TIMP1 ~ segment, data = spot)
pw <- pairwise.wilcox.test(spot$score, spot$segment, p.adjust.method = "BH")
seg_mean <- aggregate(cbind(score, TIMP1) ~ segment, data = spot, FUN = function(x) c(mean=mean(x), n=length(x)))
seg_summary <- data.frame(
  segment = seg_mean$segment,
  n = seg_mean$score[, "n"],
  score_mean = seg_mean$score[, "mean"],
  TIMP1_mean = seg_mean$TIMP1[, "mean"])
tab15 <- seg_summary
tab15$note <- ""
tab15$note[1] <- sprintf("score Kruskal-Wallis p=%.3g; TIMP1 p=%.3g", kw_seg$p.value, kw_seg_TIMP1$p.value)
save_table(tab15, "Table15_spatial_compartment_summary")

# ============================================================
# (3) Tumour vs matched normal
# ============================================================
tn <- spot[!is.na(spot$cancer_flag), ]
wt_score <- wilcox.test(score ~ cancer_flag, data = tn)
wt_timp  <- wilcox.test(TIMP1 ~ cancer_flag, data = tn)
tn_summary <- rbind(
  data.frame(Gene = "5-gene score", Tumour = mean(tn$score[tn$cancer_flag=="Tumour"]),
             Normal = mean(tn$score[tn$cancer_flag=="Normal"]),
             n_Tumour = sum(tn$cancer_flag=="Tumour"), n_Normal = sum(tn$cancer_flag=="Normal"),
             p = wt_score$p.value),
  data.frame(Gene = "TIMP1", Tumour = mean(tn$TIMP1[tn$cancer_flag=="Tumour"]),
             Normal = mean(tn$TIMP1[tn$cancer_flag=="Normal"]),
             n_Tumour = sum(tn$cancer_flag=="Tumour"), n_Normal = sum(tn$cancer_flag=="Normal"),
             p = wt_timp$p.value))
save_table(tn_summary, "Table15b_spatial_tumour_vs_normal")

message("\n=== KEY RESULTS ===")
message("TIMP1 vs CAF score, all spots: rho=", sprintf("%.2f", cor(spot$TIMP1, spot$caf, method="spearman")),
        "  p=", sprintf("%.2g", cor.test(spot$TIMP1, spot$caf, method="spearman")$p.value))
message("score ~ segment KW p=", sprintf("%.2g", kw_seg$p.value),
        "  |  Tumour vs Normal score Wilcoxon p=", sprintf("%.2g", wt_score$p.value))

# ============================================================
# FIGURES
# ============================================================
theme_set(theme_bw(base_size = 12))

# Fig14A — TIMP1 vs CAF score scatter, coloured by segment
pA <- ggplot(spot, aes(x = caf, y = TIMP1, colour = segment)) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6, size = 3.6,
           label = sprintf("Spearman rho = %.2f, p = %.2g\nn = %d spots",
                           cor.test(spot$TIMP1, spot$caf, method="spearman")$estimate,
                           cor.test(spot$TIMP1, spot$caf, method="spearman")$p.value, nrow(spot))) +
  scale_colour_manual(values = c(PanCK = "#6A8CAF", CD45 = "#C77FBF", Vimentin = "#E8843C"),
                      name = "Segment") +
  labs(title = "TIMP1 co-localises with CAF program",
       x = "CAF marker score (z)", y = "TIMP1 (z)") +
  theme(legend.position = "bottom")
save_fig(pA, "Fig14A_spatial_TIMP1_vs_CAF", w = 6.4, h = 5.4)

# Fig14B — TIMP1 by compartment segment
pB <- ggplot(spot, aes(x = segment, y = TIMP1, fill = segment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.6) +
  scale_fill_manual(values = c(PanCK = "#6A8CAF", CD45 = "#C77FBF", Vimentin = "#E8843C"), guide = "none") +
  stat_compare_means(method = "kruskal.test", label.y = max(spot$TIMP1, na.rm = TRUE) + 0.4) +
  labs(title = "TIMP1 by compartment", x = "Segment", y = "TIMP1 (z)")
save_fig(pB, "Fig14B_spatial_TIMP1_by_compartment", w = 5.6, h = 5.4)

# Fig14C — tumour vs matched normal (TIMP1)
pC <- ggplot(tn, aes(x = cancer_flag, y = TIMP1, fill = cancer_flag)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.55) +
  scale_fill_manual(values = c(Tumour = "#D7263D", Normal = "#3D5A80"), guide = "none") +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "TIMP1 in tumour vs matched normal", x = "Tissue", y = "TIMP1 (z)")
save_fig(pC, "Fig14C_spatial_tumour_vs_normal", w = 4.8, h = 5.4)

# Combined panel
suppressPackageStartupMessages(library(patchwork))
pALL <- pA + pB + pC + plot_annotation(title = "Spatial (GeoMx DSP) validation of the ferro-aging score in CRC",
                                       theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_fig(pALL, "Fig14_spatial_validation_panel", w = 16, h = 5.6)

# ---- persist -------------------------------------------------
save_obj(list(spot = spot, zmat = zmat, tab14 = tab14, tab15 = tab15,
              tn_summary = tn_summary, kw_seg = kw_seg, wt_score = wt_score),
         "spatial_DSP")

# ---- console summary for manuscript text ---------------------
allct  <- cor.test(spot$TIMP1, spot$caf, method = "spearman")
tumct  <- { s <- spot[spot$cancer_flag=="Tumour", ]; cor.test(s$TIMP1, s$caf, method="spearman") }
cat("\n----- MANUSCRIPT NUMBERS -----\n")
cat("All spots  TIMP1~CAF: rho=", sprintf("%.2f", allct$estimate), " p=", sprintf("%.3g", allct$p.value), " n=", nrow(spot), "\n")
cat("Tumour     TIMP1~CAF: rho=", sprintf("%.2f", tumct$estimate), " p=", sprintf("%.3g", tumct$p.value), " n=", sum(spot$cancer_flag=="Tumour"), "\n")
cat("Score~segment KW p=", sprintf("%.3g", kw_seg$p.value),
    " Vimentin mean=", sprintf("%.2f", seg_summary$score_mean[seg_summary$segment=="Vimentin"]),
    " PanCK=", sprintf("%.2f", seg_summary$score_mean[seg_summary$segment=="PanCK"]),
    " CD45=", sprintf("%.2f", seg_summary$score_mean[seg_summary$segment=="CD45"]), "\n")
cat("Tumour vs Normal score Wilcoxon p=", sprintf("%.3g", wt_score$p.value), "\n")
cat("TIMP1~segment KW p=", sprintf("%.3g", kw_seg_TIMP1$p.value),
    " Vimentin=", sprintf("%.2f", seg_summary$TIMP1_mean[seg_summary$segment=="Vimentin"]),
    " PanCK=", sprintf("%.2f", seg_summary$TIMP1_mean[seg_summary$segment=="PanCK"]),
    " CD45=", sprintf("%.2f", seg_summary$TIMP1_mean[seg_summary$segment=="CD45"]), "\n")
cat("TIMP1 Tumour vs Normal Wilcoxon p=", sprintf("%.3g", wt_timp$p.value),
    " (", sprintf("%.2f / %.2f", mean(tn$TIMP1[tn$cancer_flag=="Tumour"]), mean(tn$TIMP1[tn$cancer_flag=="Normal"])), ")\n")
cat("Tumour/Normal means:", sprintf("%.2f / %.2f", mean(tn$score[tn$cancer_flag=="Tumour"]), mean(tn$score[tn$cancer_flag=="Normal"])), "\n")
.score_all <- tab14[tab14$Scope == "All spots" &
                      tab14$Marker == "5-gene score vs CAF score", , drop = FALSE]
cat("5-gene~CAF all spots rho=", sprintf("%.2f", .score_all$rho[1]),
    " p=", sprintf("%.3g", .score_all$p[1]), "\n")

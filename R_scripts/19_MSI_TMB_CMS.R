# ============================================================
# 19_MSI_TMB_CMS.R — FARGs score vs molecular context
#   MSI (MANTIS), TMB (nonsynonymous mut/Mb) and the TCGA
#   PanCanAtlas molecular subtype (CIN vs hypermutated) are
#   candidate confounders for a prognostic score. This script
#   (i) tests whether the FARGs score is associated with any of
#   them, and (ii) checks whether the score stays an independent
#   prognostic factor after adjustment for MSI + TMB on top of
#   age/sex/stage.
#
#   Data source (fetched, not simulated): cBioPortal Public API,
#   study coadread_tcga_pan_can_atlas_2018, sample-level clinical
#   attributes MSI_SCORE_MANTIS / MSI_SENSOR_SCORE /
#   TMB_NONSYNONYMOUS / MUTATION_COUNT / SUBTYPE, pulled per
#   primary-tumour sample by scripts_fetch_cbio.py into
#   data/processed/cbio_msi_tmb_subtype.csv and joined on the
#   12-character TCGA patient barcode.
#
#   NOTE on CMS: cBioPortal exposes the PanCanAtlas consensus
#   subtype (COAD_CIN / COAD_hypermutated / READ_*), NOT the
#   Guinney CMS1-4 call. We therefore report the PanCanAtlas
#   subtype explicitly and do not relabel it as CMS.
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
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ggpubr))

# ---- inputs --------------------------------------------------
clin <- readRDS(file.path(RES_DIR, "clinical_table1_df.rds"))     # patient, score, os_*, age, sex, stage
cbio <- read.csv(file.path(DATA_PROC, "cbio_msi_tmb_subtype.csv"),
                 stringsAsFactors = FALSE, na.strings = c("", "NA"))
# SUBTYPE is a patient-level attribute (empty in the sample-level pull) -> merge from
# the patient-level table fetched by scripts_fetch_cbio_pat.py
sub_p <- read.csv(file.path(DATA_PROC, "cbio_subtype_patient.csv"),
                  stringsAsFactors = FALSE, na.strings = c("", "NA"))

num <- function(x) as.numeric(as.character(x))
cbio$MSI_SCORE_MANTIS <- num(cbio$MSI_SCORE_MANTIS)
cbio$MSI_SENSOR_SCORE <- num(cbio$MSI_SENSOR_SCORE)
cbio$TMB_NONSYNONYMOUS <- num(cbio$TMB_NONSYNONYMOUS)
cbio$MUTATION_COUNT <- num(cbio$MUTATION_COUNT)

# de-dup cbio to one row per patient (primary tumour), then attach patient-level SUBTYPE
cbio <- cbio[!duplicated(cbio$patientId), ]
cbio$SUBTYPE <- sub_p$SUBTYPE[match(cbio$patientId, sub_p$patientId)]

df <- merge(clin, cbio, by.x = "patient", by.y = "patientId", all.x = FALSE)
message("Modelled patients: ", nrow(clin), "  |  joined with molecular data: ", nrow(df))

# ---------- ESTIMATE-like stromal score (purity/stromal-confounding test) ----------
# Reviewer concern: a stroma-weighted index may simply proxy tumour purity / stromal
# content. We compute an ESTIMATE-like stromal score per patient (same marker set and
# per-gene z-scoring as script 09) and a TIMP1-free variant (dropping the hub gene to
# avoid circularity), then adjust the prognostic model for them below.
suppressPackageStartupMessages(library(data.table))
expr_dt19 <- data.table::fread(file.path(DATA_RAW, "TCGA_COADREAD_TPM_matrix.tsv"))
em19 <- as.matrix(expr_dt19[, -1, with = FALSE]); rownames(em19) <- expr_dt19[[1]]
colnames(em19) <- sub("^\\.\\d+", "", colnames(em19)); rm(expr_dt19)
sinfo19 <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_sample_types.csv"), stringsAsFactors = FALSE)
pt19 <- sinfo19 %>% filter(sample_type == "Primary Tumor", sample_id %in% colnames(em19)) %>%
  group_by(patient_barcode) %>% summarise(sample_id = min(sample_id), .groups = "drop") %>% as.data.frame()
ep19 <- log2(em19[, pt19$sample_id, drop = FALSE] + 1); rownames(ep19) <- rownames(em19)
z19 <- (ep19 - rowMeans(ep19)) / apply(ep19, 1, sd)          # per-gene z across patients
stromal_markers <- intersect(c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1",
                               "FN1","LUM","POSTN","SPARC","TIMP1","TIMP2","THBS1",
                               "MMP2","MMP14","ACTA2","TAGLN","PDGFRB","CSPG4","VIM"),
                             rownames(ep19))
stromal_markers_nt <- setdiff(stromal_markers, "TIMP1")       # TIMP1-free control
strom_pat <- colMeans(z19[stromal_markers, , drop = FALSE])
strom_nt_pat <- colMeans(z19[stromal_markers_nt, , drop = FALSE])
stromal_all <- as.numeric(scale(strom_pat)); names(stromal_all) <- pt19$patient_barcode
stromal_nt_all <- as.numeric(scale(strom_nt_pat)); names(stromal_nt_all) <- pt19$patient_barcode
df$stromal   <- stromal_all[df$patient]
df$stromal_nt<- stromal_nt_all[df$patient]
rho_strom <- suppressWarnings(cor.test(df$score, df$stromal, method = "spearman", use = "complete.obs"))
message(sprintf("Stromal score: n=%d | Spearman score~stromal rho=%.3f (p=%.3g) | markers=%d (TIMP1-free=%d)",
                sum(!is.na(df$stromal)), rho_strom$estimate, rho_strom$p.value,
                length(stromal_markers), length(stromal_markers_nt)))

# MSI status from MANTIS (thresholds per cBioPortal: >0.6 MSI, 0.4-0.6 indeterminate, <0.4 MSS)
df$MSI_status <- ifelse(df$MSI_SCORE_MANTIS > 0.6, "MSI-H",
                 ifelse(df$MSI_SCORE_MANTIS < 0.4, "MSS", "Indeterminate"))
# collapse PanCanAtlas subtype into CIN-like vs hypermutated
df$sub_cat <- ifelse(grepl("MSI|POLE", df$SUBTYPE), "Hypermutated",
              ifelse(grepl("CIN|GS", df$SUBTYPE), "CIN-like", NA))

# ---------- (a) score ~ MSI ----------
msi_two <- subset(df, MSI_status %in% c("MSS", "MSI-H"))
wil_msi <- wilcox.test(score ~ MSI_status, data = msi_two)
msi_summary <- msi_two %>% group_by(MSI_status) %>%
  summarise(n = n(), median_score = median(score), mean_score = mean(score),
            sd_score = sd(score), .groups = "drop")

# ---------- (b) score ~ TMB ----------
cor_tmb <- suppressMessages(cor.test(df$score, df$TMB_NONSYNONYMOUS,
                                     method = "spearman", use = "complete.obs"))
df$tmb_tertile <- cut(df$TMB_NONSYNONYMOUS,
                      breaks = quantile(df$TMB_NONSYNONYMOUS, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      labels = c("T1 low", "T2 mid", "T3 high"), include.lowest = TRUE)
krusk_tmb <- kruskal.test(score ~ tmb_tertile, data = df)

# ---------- (c) score ~ PanCan subtype ----------
wil_sub <- wilcox.test(score ~ sub_cat, data = df)
sub_summary <- df %>% group_by(sub_cat, SUBTYPE) %>%
  summarise(n = n(), median_score = median(score), mean_score = mean(score), .groups = "drop")
kw_sub <- kruskal.test(score ~ SUBTYPE, data = df[!is.na(df$SUBTYPE) & df$SUBTYPE != "", ])

# Label the residual bucket explicitly for the shipped table (df$sub_cat itself is
# deliberately left NA so that wilcox.test()/ggboxplot() stay strictly two-group).
sub_summary_out <- sub_summary
sub_summary_out$sub_cat[is.na(sub_summary_out$sub_cat)] <- "Unclassified"
sub_summary_out$SUBTYPE[is.na(sub_summary_out$SUBTYPE)] <- "SUBTYPE not called"

# ---------- Table 12: score vs molecular context ----------
tab12 <- data.frame(
  Variable = c(
    "MSI status (MANTIS)",
    "  MSS vs MSI-H",
    "TMB (nonsynonymous, per Mb)",
    "  Spearman rho",
    "  T1 vs T3 medians",
    "PanCanAtlas subtype",
    "  CIN vs Hypermutated"),
  Statistic = c(
    "Wilcoxon rank-sum", "",
    "Spearman correlation", sprintf("rho = %.3f", cor_tmb$estimate),
    sprintf("Kruskal-Wallis p = %.3g", krusk_tmb$p.value),
    "Wilcoxon rank-sum", ""),
  p_value = c(
    sprintf("%.3g", wil_msi$p.value), "",
    sprintf("%.3g", cor_tmb$p.value), "", "",
    sprintf("%.3g", wil_sub$p.value), ""),
  Detail = c(
    paste0("MSS median ", round(msi_summary$median_score[msi_summary$MSI_status == "MSS"], 3),
           " vs MSI-H median ", round(msi_summary$median_score[msi_summary$MSI_status == "MSI-H"], 3)),
    paste0("n = ", msi_summary$n[msi_summary$MSI_status == "MSS"], " (MSS) vs ",
           msi_summary$n[msi_summary$MSI_status == "MSI-H"], " (MSI-H)"),
    paste0("n = ", sum(!is.na(df$TMB_NONSYNONYMOUS))), "", "",
    paste0("CIN-like median ", round(sub_summary$median_score[sub_summary$sub_cat == "CIN-like"][1], 3),
           " vs Hypermut median ", round(sub_summary$median_score[sub_summary$sub_cat == "Hypermutated"][1], 3)),
    paste0("n = ", sum(df$sub_cat == "CIN-like", na.rm = TRUE), " vs ", sum(df$sub_cat == "Hypermutated", na.rm = TRUE))),
  stringsAsFactors = FALSE)
save_table(tab12, "Table12_score_vs_molecular_context")

# per-group detail tables
save_table(msi_summary, "Table12b_msi_groups")
save_table(sub_summary_out, "Table12c_subtype_groups")

# ---------- (d) independence: Cox adjusted for MSI + TMB + stromal content ----------
dsub <- df %>% filter(!is.na(MSI_SCORE_MANTIS), !is.na(TMB_NONSYNONYMOUS), !is.na(os_time))
dsub$os_status <- as.numeric(dsub$os_status)
dsub$stage_num <- ifelse(grepl("IV", dsub$stage_grp), 4,
                  ifelse(grepl("III", dsub$stage_grp), 3,
                  ifelse(grepl("II",  dsub$stage_grp), 2, 1)))
dsub$sex_m <- ifelse(dsub$sex == "male", 1, 0)
dsub$MSI_H <- ifelse(dsub$MSI_status == "MSI-H", 1, 0)

form_score_only <- coxph(Surv(os_time, os_status) ~ score, data = dsub)
form_clin       <- coxph(Surv(os_time, os_status) ~ score + age_years + sex_m + factor(stage_num), data = dsub)
form_clin_msi   <- coxph(Surv(os_time, os_status) ~ score + age_years + sex_m + factor(stage_num) + MSI_H + TMB_NONSYNONYMOUS, data = dsub)
# stromal / purity adjustment (reviewer concern that the index just proxies stroma)
dsub_st <- dsub %>% filter(!is.na(stromal), !is.na(stromal_nt))
form_clin_stromal    <- coxph(Surv(os_time, os_status) ~ score + age_years + sex_m + factor(stage_num) + stromal,    data = dsub_st)
form_clin_stromal_nt <- coxph(Surv(os_time, os_status) ~ score + age_years + sex_m + factor(stage_num) + stromal_nt, data = dsub_st)
form_stromal_only    <- coxph(Surv(os_time, os_status) ~ stromal, data = dsub_st)

# tidy HR extraction for any named term
row_of <- function(fit, term = "score") {
  s <- summary(fit); ci <- s$conf.int
  if (!term %in% rownames(s$coefficients)) return(NULL)
  data.frame(model = "", term = term, HR = ci[term, "exp(coef)"],
             lower = ci[term, 3], upper = ci[term, 4],
             p = s$coefficients[term, "Pr(>|z|)"],
             C_index = s$concordance[1],
             n_patients = fit$n, n_deaths = nobs(fit))
}
tab13 <- rbind(
  within(row_of(form_score_only),     { model <- "Univariable (score)" }),
  within(row_of(form_clin),           { model <- "+ age, sex, stage" }),
  within(row_of(form_clin_msi),       { model <- "+ age, sex, stage, MSI-H, TMB" }),
  within(row_of(form_clin_stromal),   { model <- "+ age, sex, stage, stromal score" }),
  within(row_of(form_clin_stromal_nt),{ model <- "+ age, sex, stage, stromal score (TIMP1-free)" }),
  within(row_of(form_stromal_only, "stromal"), { model <- "Stromal score alone (reference)" })
)
rownames(tab13) <- NULL
tab13$p <- as.numeric(tab13$p)
tab13 <- tab13[, c("model", "term", "HR", "lower", "upper", "p", "C_index", "n_patients", "n_deaths")]
save_table(tab13, "Table13_score_independence_molecular")

# full-cohort (all modelled patients) discrimination contrasts for prose / guards
fc <- clin; fc$os_status <- as.numeric(fc$os_status)
fc$stromal <- stromal_all[fc$patient]
fc <- fc %>% filter(!is.na(os_time), os_time > 0, !is.na(stromal))
c_score  <- survConcordance(Surv(os_time, os_status) ~ score, data = fc)$concordance
c_strom  <- survConcordance(Surv(os_time, os_status) ~ stromal, data = fc)$concordance
c_both   <- survConcordance(Surv(os_time, os_status) ~ predict(coxph(Surv(os_time, os_status) ~ score + stromal, data = fc)), data = fc)$concordance
rho_sc   <- suppressWarnings(cor.test(fc$score, fc$stromal, method = "spearman"))$estimate
strom_contrast <- data.frame(
  Quantity = c("Score-alone C-index", "Stromal-alone C-index", "Score + stromal C-index",
               "Spearman rho (score ~ stromal)", "n modelled (stroma-adj cohort)"),
  Value    = c(round(c_score, 4), round(c_strom, 4), round(c_both, 4), round(rho_sc, 3), nrow(fc)),
  stringsAsFactors = FALSE)
save_table(strom_contrast, "Table13c_stroma_purity_contrast")
message("Full-cohort stroma contrast: score C=", round(c_score,4),
        " stromal-alone C=", round(c_strom,4), " score+stroma C=", round(c_both,4),
        " rho=", round(rho_sc,3), " n=", nrow(fc))

# likelihood-ratio: does adding MSI+TMB improve the clinical model? (confounder relevance)
lr_msi <- anova(form_clin, form_clin_msi, test = "Chisq")

# full adjusted Cox coefficients (for supplementary transparency)
supp <- as.data.frame(summary(form_clin_msi)$coefficients)
supp$HR <- exp(supp$coef); supp$term <- rownames(supp)
save_table(supp[, c("term", "coef", "HR", "se(coef)", "z", "Pr(>|z|)")],
           "Table13b_full_adjusted_cox")

# ---------- figures ------------------------------------------
pA <- ggboxplot(msi_two, x = "MSI_status", y = "score", color = "MSI_status",
                add = "jitter", addParams = list(size = 1, alpha = 0.4), fill = "MSI_status") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
  labs(title = "Score by MSI status", x = "", y = "FARGs score",
       subtitle = sprintf("MSS n=%d | MSI-H n=%d",
                          sum(msi_two$MSI_status == "MSS"), sum(msi_two$MSI_status == "MSI-H"))) +
  scale_fill_manual(values = c(MSS = "#3D5A80", `MSI-H` = "#D7263D")) +
  scale_color_manual(values = c(MSS = "#3D5A80", `MSI-H` = "#D7263D")) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(size = 12, face = "bold"))

df$tmb_label <- factor(df$tmb_tertile, levels = c("T1 low", "T2 mid", "T3 high"))
pB <- ggplot(df[!is.na(df$TMB_NONSYNONYMOUS), ], aes(x = TMB_NONSYNONYMOUS, y = score)) +
  geom_point(size = 1.4, alpha = 0.55, colour = "#1B998B") +
  geom_smooth(method = "lm", se = TRUE, colour = "#D7263D", linewidth = 0.8) +
  annotate("text", x = min(df$TMB_NONSYNONYMOUS, na.rm = TRUE),
           y = max(df$score, na.rm = TRUE),
           label = sprintf("Spearman rho = %.2f (p = %.2g)", cor_tmb$estimate, cor_tmb$p.value),
           hjust = 0, vjust = 1.5, size = 3.8) +
  labs(title = "Score vs TMB", x = "TMB (nonsynonymous mut/Mb)",
       y = "FARGs score", subtitle = sprintf("n = %d", sum(!is.na(df$TMB_NONSYNONYMOUS)))) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(size = 12, face = "bold"))

pC <- ggboxplot(df[!is.na(df$sub_cat), ], x = "sub_cat", y = "score", color = "sub_cat",
                add = "box", width = 0.55,
                xlab = "", ylab = "FARGs score",
                title = "Score by consensus subtype",
                subtitle = sprintf("CIN-like n=%d | Hypermut n=%d",
                                   sum(df$sub_cat == "CIN-like", na.rm = TRUE),
                                   sum(df$sub_cat == "Hypermutated", na.rm = TRUE))) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5) +
  scale_color_manual(values = c(`CIN-like` = "#3D5A80", Hypermutated = "#E3A511")) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(size = 12, face = "bold"))

fig12 <- pA + pB + pC + patchwork::plot_layout(ncol = 3)
save_fig(fig12, "Fig13A-C_score_vs_MSI_TMB_subtype", w = 13, h = 4.6)

# KM by subtype (does the score's risk split differ across CIN/hypermut?) — supplementary figure
df_hi <- df %>% mutate(risk = ifelse(score > median(score, na.rm = TRUE), "High", "Low"))
fit_sub <- survfit(Surv(os_time, os_status) ~ sub_cat, data = df_hi)
pD <- ggsurvplot(fit_sub, data = df_hi, pval = TRUE, risk.table = TRUE, size = 1,
                 legend.title = "Subtype", legend.labs = c("CIN-like", "Hypermutated"),
                 palette = c("#3D5A80", "#E3A511"),
                 title = "Overall survival by PanCanAtlas subtype")
ggsave(file.path(FIG_DIR, "Fig13D_KM_by_subtype.png"), pD$plot, width = 6, height = 4.5, dpi = FIG_DPI, bg = "white")

# ---------- persist ------------------------------------------
save_obj(list(df = df, tab12 = tab12, tab13 = tab13,
              wil_msi = wil_msi, cor_tmb = cor_tmb, wil_sub = wil_sub,
              lr_msi = lr_msi, forms = list(uni = form_score_only, clin = form_clin,
                                            clin_msi = form_clin_msi)),
         "MSI_TMB_subtype")

cat("\n===== SUMMARY (MSI/TMB/subtype) =====\n")
cat("joined n:", nrow(df), "\n")
print(tab12); cat("\n"); print(tab13)
cat("LR test add MSI+TMB to clinical: chisq =", round(lr_msi[2, "Chisq"], 3),
    "p =", sprintf("%.4g", lr_msi[2, "Pr(>|Chi|)"]), "\n")

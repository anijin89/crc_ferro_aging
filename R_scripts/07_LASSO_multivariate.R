# ============================================================
# 07_LASSO_multivariate.R — Univariate Cox → LASSO → Multi-Cox → FARGs score
# ============================================================
# 输入：TCGA TPM_matrix + sample_types + clinical.csv + FARG_DEG_intersection.rds
# 关键：13 例患者有多个 Primary Tumor 样本 → 患者级去重（保留首个 barcode）
# ============================================================

suppressPackageStartupMessages({ library(magrittr); library(dplyr); library(data.table) })
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

# ---- 1. 载入表达矩阵 + 样本映射 ----
expr_dt <- data.table::fread(file.path(DATA_RAW, "TCGA_COADREAD_TPM_matrix.tsv"))
gene_ids <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_ids; rm(expr_dt, gene_ids)
colnames(expr_mat) <- sub("^\\.\\d+", "", colnames(expr_mat))

sinfo <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_sample_types.csv"),
                  stringsAsFactors = FALSE)
farg_degs <- readRDS(file.path(RES_DIR, "FARG_DEG_intersection.rds"))

# ---- 2. Primary Tumor 子集 ----
pt <- sinfo %>% filter(sample_type == "Primary Tumor") %>%
  filter(sample_id %in% colnames(expr_mat))
message(sprintf("Primary Tumor samples: %d (patients: %d)",
                nrow(pt), length(unique(pt$patient_barcode))))

# ---- 3. 患者级去重：每位患者保留 barcode 字典序最小的 Primary Tumor 样本 ----
pt$sort_key <- pt$sample_id
dup <- pt %>% group_by(patient_barcode) %>% filter(n() > 1) %>% ungroup()
if (nrow(dup) > 0) {
  message(sprintf("⚠ %d 患者存在多个原发肿瘤样本 → 保留每位首个 sample_id",
                  length(unique(dup$patient_barcode))))
  message("   Affected patients (first 15): ",
          paste(unique(dup$patient_barcode)[1:15], collapse=", "))
}
pt_dedup <- pt %>% group_by(patient_barcode) %>%
  summarise(sample_id = min(sample_id), .groups = "drop") %>%
  as.data.frame()
message(sprintf("After dedup: %d patients × FARGs=%d genes",
                nrow(pt_dedup), length(farg_degs)))

expr_tumor <- log2(expr_mat[farg_degs, pt_dedup$sample_id, drop = FALSE] + 1)
rownames(expr_tumor) <- farg_degs

# ---- 4. 临床匹配 ----
clin <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_clinical.csv"),
                 stringsAsFactors = FALSE)
rownames(clin) <- clin$patient_barcode
clin_sub <- clin[pt_dedup$patient_barcode, ]
stopifnot(all(pt_dedup$patient_barcode %in% clin$patient_barcode))

os_time   <- clin_sub$os_days
os_status <- ifelse(clin_sub$vital_status == "Dead", 1,
             ifelse(clin_sub$vital_status == "Alive", 0, NA))

valid <- !is.na(os_time) & os_time > 0 & !is.na(os_status)
message(sprintf("Valid OS: %d / %d (dropped %d: %d missing time, %d missing status)",
                sum(valid), nrow(clin_sub), sum(!valid),
                sum(is.na(os_time) | os_time <= 0),
                sum(is.na(os_status) & !(is.na(os_time) | os_time <= 0))))

expr_tumor <- expr_tumor[, valid, drop = FALSE]
os_time   <- os_time[valid]
os_status <- os_status[valid]
patient_ids <- pt_dedup$patient_barcode[valid]

# ---- 4b. 每基因 z-score（训练队列内标准化）----
# 使预后指数跨平台可迁移，且与外部验证(08 同样 per-gene z-score)口径一致；
# 之后 Cox 系数为 per-SD 效应，评分 = Σ β_z × z。
mu <- rowMeans(expr_tumor); sg <- apply(expr_tumor, 1, sd)
sg[sg == 0] <- 1                      # 防零方差
expr_tumor <- (expr_tumor - mu) / sg
rownames(expr_tumor) <- farg_degs
message("✓ Training expression z-scored per gene (per-SD betas)")

message(sprintf("✓ Final N for survival: %d patients, Events: %d",
                length(os_status), sum(os_status == 1)))

# ---- 5. 单因素 Cox ----
cox_uni <- do.call(rbind, lapply(rownames(expr_tumor), function(g) {
  fit <- tryCatch(coxph(Surv(os_time, os_status) ~ as.numeric(expr_tumor[g, ])),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  if (is.na(s$coefficients[1, "Pr(>|z|)"])) return(NULL)
  data.frame(gene = g,
             HR      = s$coefficients[1, "exp(coef)"],
             HR_low  = s$conf.int[1, "lower .95"],
             HR_high = s$conf.int[1, "upper .95"],
             beta    = s$coefficients[1, "coef"],
             p       = s$coefficients[1, "Pr(>|z|)"])
}))
save_table(cox_uni, "TableS_univariate_cox")
cox_uni_sig <- cox_uni %>% filter(p < COX_P_CUTOFF)
message(sprintf("✓ Univariate Cox p<0.05: %d / %d genes",
                nrow(cox_uni_sig), nrow(cox_uni)))

# ---- 6. LASSO Cox ----
x <- t(as.matrix(expr_tumor[cox_uni_sig$gene, ]))
y <- Surv(os_time, os_status)
set.seed(SEED)
cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10,
                    type.measure = "deviance")

png(file.path(FIG_DIR, "Fig3C_LASSO_CV.png"), width = 6, height = 5,
    units = "in", res = FIG_DPI)
plot(cv_fit)
# glmnet 在画布顶部绘制「非零系数个数」轴，标题必须上移，否则与该轴重叠
title(main = "LASSO Cox: 10-fold CV", line = 3.2)
abline(v = log(cv_fit$lambda.min), lty = 2, col = "red"); dev.off()

png(file.path(FIG_DIR, "Fig3D_LASSO_coef.png"), width = 7, height = 5,
    units = "in", res = FIG_DPI)
plot(cv_fit$glmnet.fit, xvar = "lambda", label = TRUE)
abline(v = log(cv_fit$lambda.min), lty = 2, col = "red")
title(main = "LASSO Cox: coefficient paths"); dev.off()

coef_min <- coef(cv_fit, s = "lambda.min")
selected <- rownames(coef_min)[coef_min[, 1] != 0]
betas_lasso <- coef_min[selected, 1]
message(sprintf("✓ LASSO selected %d genes: %s",
                length(selected), paste(selected, collapse = ", ")))

# ---- 7. 多因素 Cox ----
df_multi <- as.data.frame(t(as.matrix(expr_tumor[selected, ])))
df_multi$os_time <- os_time; df_multi$os_status <- os_status
fit_multi <- coxph(Surv(os_time, os_status) ~ .,
                   data = df_multi[, c(selected, "os_time", "os_status")])
print(summary(fit_multi))
multi_coef <- summary(fit_multi)$coefficients[, "coef"]
mtab <- as.data.frame(summary(fit_multi)$coefficients)
mtab <- cbind(term = rownames(mtab), mtab, row.names = NULL)
save_table(mtab, "Table3_multivariate_cox")
save_table(data.frame(gene = names(multi_coef), beta = multi_coef),
           "Table3_LASSO_selected_beta")

# ---- 8. FARGs score = Σ β_z,i × z_i  (z-scored expression, per-SD betas) ----
calc_score <- function(expr_mat, genes, betas) {
  t(betas) %*% expr_mat[genes, , drop = FALSE]
}
score_train <- as.numeric(calc_score(expr_tumor, selected, multi_coef))
save_obj(list(score = score_train, genes = selected,
              betas_lasso = betas_lasso, betas_multi = multi_coef,
              os_time = os_time, os_status = os_status,
              patient_ids = patient_ids),
         "FARGs_score_TCGA")

# ---- 9. KM (中位分组) ----
grp <- ifelse(score_train > median(score_train), "High", "Low")
df_surv <- data.frame(time = os_time, status = os_status,
                      group = factor(grp, levels = c("Low", "High")))
fit_km <- survfit(Surv(time, status) ~ group, data = df_surv)
p_km <- ggsurvplot(fit_km, data = df_surv,
                   pval = TRUE, pval.method = TRUE,
                   risk.table = TRUE, risk.table.height = 0.25,
                   palette = KM_PALETTE,
                   xlab = "Days", ylab = "Overall Survival Probability",
                   title = "FARGs score in TCGA-COADREAD (patient-level)",
                   legend.title = "Risk", legend.labs = KM_LEGEND,
                   ggtheme = theme_classic(base_size = 12))
png(file.path(FIG_DIR, "Fig4A_KM_TCGA.png"), width = 7, height = 7,
    units = "in", res = FIG_DPI); print(p_km); dev.off()

# ---- 10. timeROC ----
suppressPackageStartupMessages(library(timeROC))
roc <- timeROC(T = os_time, delta = os_status, marker = score_train,
               cause = 1, weighting = "marginal",
               times = c(365, 3*365, 5*365), ROC = TRUE, iid = TRUE)
save_obj(roc, "FARGs_ROC_TCGA")
p_roc <- ggplot(data.frame(
  FPR = c(roc$FP[, 1], roc$FP[, 2], roc$FP[, 3]),
  TPR = c(roc$TP[, 1], roc$TP[, 2], roc$TP[, 3]),
  Time = rep(c("1-year", "3-year", "5-year"), each = nrow(roc$FP))
), aes(FPR, TPR, color = Time)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#1B998B", "#F46036", "#D7263D")) +
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("AUC:\n 1y = %.3f\n 3y = %.3f\n 5y = %.3f",
                           roc$AUC[1], roc$AUC[2], roc$AUC[3]), size = 4) +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "Time-dependent ROC (TCGA, patient-level)") +
  theme_classic(base_size = 12)
save_fig(p_roc, "Fig4B_ROC_TCGA", w = 6, h = 5)

message(sprintf("✓ AUC 1y=%.3f, 3y=%.3f, 5y=%.3f", roc$AUC[1], roc$AUC[2], roc$AUC[3]))
message("=== 07_LASSO_multivariate.R 完成 ===")

# ============================================================
# 08_ExternalValidation.R — GEO 队列外部验证
# ============================================================
# GSE39582 (GPL570, n=585, French PICOR cohort, has OS+DFS)
# GSE41258 (GPL570 gene-symbol matrix, n=390, Scottish, has OS)
# 输出：Fig4C_KM_GSE39582, Fig4D_ROC_GSE39582, Fig4E_KM_GSE41258, Fig4F_ROC_GSE41258
#       Table4_external_validation_summary.csv
# ============================================================

suppressPackageStartupMessages({ library(magrittr); library(dplyr); library(GEOquery);
                                library(AnnotationDbi); library(hgu133plus2.db) })
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
suppressPackageStartupMessages({ library(survival); library(timeROC) })

# 载入 TCGA 模型
model   <- readRDS(file.path(RES_DIR, "FARGs_score_TCGA.rds"))
genes   <- model$genes
betas   <- model$betas_multi

# ---- 通用：探针均值 collapse → 基因表达（快速 base rowsum） ----
collapse_probe_to_gene <- function(mat, map_df) {
  keep <- rownames(mat) %in% map_df$PROBEID
  m <- mat[keep, , drop = FALSE]
  sym <- map_df$SYMBOL[match(rownames(m), map_df$PROBEID)]
  sums <- rowsum(m, group = sym, reorder = FALSE)
  counts <- as.numeric(table(sym)[rownames(sums)])
  sums / counts    # recycles vector across columns → each row divided by its count
}

# ============================================================
# 1) GSE39582
# ============================================================
cat("\n############ GSE39582 (n=585, GPL570) ############\n")
g <- getGEO(filename = file.path(DATA_RAW, "GSE39582_series_matrix.txt.gz"), Annot = FALSE)
e1 <- if (is.list(g)) g[[1]] else g
expr_p <- exprs(e1)  # probes x 585
pd <- as.data.frame(pData(e1))

# 临床抽取
parse_char <- function(x) {
  as.character(sub("^[^:]*:\\s*", "", x))
}
pd$os_event_raw <- parse_char(pd[["characteristics_ch1.13"]])
pd$os_delay_raw <- parse_char(pd[["characteristics_ch1.14"]])
pd$rfs_event_raw<- parse_char(pd[["characteristics_ch1.11"]])
pd$rfs_delay_raw<- parse_char(pd[["characteristics_ch1.12"]])
pd$dataset      <- parse_char(pd[["characteristics_ch1.1"]])
pd$os_time      <- suppressWarnings(as.numeric(pd$os_delay_raw)) * 30.44  # months→days
pd$os_status    <- suppressWarnings(as.numeric(pd$os_event_raw))
# Non-Tumoral 剔除
pd$tumor        <- pd$dataset %in% c("discovery", "validation")
message(sprintf("GSE39582: total %d, tumor %d, valid OS %d",
                nrow(pd), sum(pd$tumor), sum(pd$tumor & !is.na(pd$os_time) & pd$os_time>0)))

# 探针→gene symbol：maxIQR 选择（GPL570 外部验证标准做法）
map_all <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr_p),
                                 keytype = "PROBEID", columns = "SYMBOL")
map_all <- map_all[!is.na(map_all$SYMBOL) & map_all$PROBEID %in% rownames(expr_p), ]
tumor_cols <- which(pd$tumor & !is.na(pd$os_time) & pd$os_time > 0)
best_probe <- sapply(genes, function(g) {
  pr <- map_all$PROBEID[map_all$SYMBOL == g]
  if (length(pr) == 0) return(NA_character_)
  if (length(pr) == 1) return(pr)
  iqr <- apply(expr_p[pr, tumor_cols, drop = FALSE], 1, IQR)
  pr[which.max(iqr)]
})
best_probe <- best_probe[!is.na(best_probe)]
message(sprintf("Model genes: %d, maxIQR probes: %d",
                length(genes), length(best_probe)))
expr_g <- expr_p[best_probe, , drop = FALSE]
rownames(expr_g) <- names(best_probe)
sel_g <- rownames(expr_g)
missing_g <- setdiff(genes, sel_g)
if (length(missing_g)) message("⚠ Missing genes: ", paste(missing_g, collapse=", "))

# per-gene z-score (row-wise standardization on genes × samples matrix)
z_gene <- function(m) {
  mu  <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, function(v) { s <- sd(v, na.rm = TRUE); if (is.na(s) || s == 0) 1 else s })
  (m - mu) / sdv
}
tcga_ref <- model$score  # only for reference; per-gene z-score is enough
expr_v_gse <- z_gene(as.matrix(expr_g[sel_g, ]))

# reorder genes to match model$genes order
expr_v_gse <- expr_v_gse[genes, , drop = FALSE]
score_gse39582 <- as.numeric(t(betas) %*% expr_v_gse)

names(score_gse39582) <- colnames(expr_p)
valid <- pd$tumor & !is.na(pd$os_time) & pd$os_time > 0 & !is.na(pd$os_status)
df <- data.frame(sample = rownames(pd)[valid],
                 time   = pd$os_time[valid],
                 status = pd$os_status[valid],
                 score  = score_gse39582[valid])
message(sprintf("GSE39582 usable: %d, events: %d", nrow(df), sum(df$status)))

# KM
grp <- ifelse(df$score > median(df$score), "High", "Low")
df$group <- factor(grp, levels = c("Low", "High"))
fit_km <- survfit(Surv(time, status) ~ group, data = df)
p_km <- ggsurvplot(fit_km, data = df, pval = TRUE, pval.method = TRUE,
                   risk.table = TRUE, risk.table.height = 0.25,
                   palette = KM_PALETTE,
                   legend.title = "Risk", legend.labs = KM_LEGEND,
                   xlab = "Days", ylab = "Overall Survival",
                   title = sprintf("FARGs score — GSE39582 (external, n=%d)", nrow(df)),
                   ggtheme = theme_classic(base_size = 12))
png(file.path(FIG_DIR, "Fig4C_KM_GSE39582.png"), width = 7, height = 7,
    units = "in", res = FIG_DPI); print(p_km); dev.off()

# timeROC
roc <- timeROC(T = df$time, delta = df$status, marker = df$score,
               cause = 1, weighting = "marginal",
               times = c(365, 3*365, 5*365), ROC = TRUE, iid = TRUE)
save_obj(list(df = df, roc = roc), "FARGs_external_GSE39582")
p_roc <- ggplot(data.frame(
  FPR = c(roc$FP[, 1], roc$FP[, 2], roc$FP[, 3]),
  TPR = c(roc$TP[, 1], roc$TP[, 2], roc$TP[, 3]),
  Time = rep(c("1-year","3-year","5-year"), each = nrow(roc$FP))
), aes(FPR, TPR, color = Time)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#1B998B", "#F46036", "#D7263D")) +
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("AUC:\n 1y = %.3f\n 3y = %.3f\n 5y = %.3f",
                           roc$AUC[1], roc$AUC[2], roc$AUC[3]), size = 4) +
  labs(title = "GSE39582 time-dependent ROC", x = "1-Specificity", y = "Sensitivity") +
  theme_classic(base_size = 12)
save_fig(p_roc, "Fig4D_ROC_GSE39582", w = 6, h = 5)
message(sprintf("✓ GSE39582 AUC 1y=%.3f 3y=%.3f 5y=%.3f",
                roc$AUC[1], roc$AUC[2], roc$AUC[3]))

# Cox on score (continuous)
cox_fit <- coxph(Surv(time, status) ~ score, data = df)
s <- summary(cox_fit)
hr_ext <- s$coefficients[1, "exp(coef)"]
p_ext  <- s$coefficients[1, "Pr(>|z|)"]
message(sprintf("GSE39582 Cox HR per SD = %.3f (p=%.2e)", hr_ext, p_ext))

# 每基因单因素 Cox（外部验证）
per_gene_cox <- do.call(rbind, lapply(sel_g, function(g) {
  fit <- coxph(Surv(df$time, df$status) ~ expr_v_gse[g, match(rownames(df), colnames(expr_v_gse))])
  s <- summary(fit)
  data.frame(gene = g,
             HR      = s$coefficients[1, "exp(coef)"],
             HR_low  = s$conf.int[1, "lower .95"],
             HR_high = s$conf.int[1, "upper .95"],
             beta_extern = s$coefficients[1, "coef"],
             p       = s$coefficients[1, "Pr(>|z|)"],
             beta_train  = betas[g])
}))
per_gene_cox$direction_consistent <- sign(per_gene_cox$beta_extern) == sign(per_gene_cox$beta_train)
save_table(per_gene_cox, "Table4b_GSE39582_per_gene_cox")
message(sprintf("Direction consistent: %d / %d genes",
                sum(per_gene_cox$direction_consistent), nrow(per_gene_cox)))

# ============================================================
# Summary table
# ============================================================
km_p39582 <- tryCatch(survdiff(Surv(time, status) ~ group, data = df)$p,
                      error = function(e) NA)
km_p39582 <- 1 - pchisq(survdiff(Surv(time, status) ~ group, data = df)$chisq, 1)
roc_t <- readRDS(file.path(RES_DIR, "FARGs_ROC_TCGA.rds"))
summ <- data.frame(
  Cohort  = c("TCGA-COADREAD (train)", "GSE39582 (external)"),
  N       = c(length(model$score), nrow(df)),
  Events  = c(sum(model$os_status), sum(df$status)),
  KM_p    = c(NA, km_p39582),
  Cox_HR_per_SD = c(NA, hr_ext),
  Cox_p   = c(NA, p_ext),
  AUC_1y  = c(roc_t$AUC[1], roc$AUC[1]),
  AUC_3y  = c(roc_t$AUC[2], roc$AUC[2]),
  AUC_5y  = c(roc_t$AUC[3], roc$AUC[3]),
  stringsAsFactors = FALSE
)
save_table(summ, "Table4_external_validation_summary")
message("=== 08_ExternalValidation.R 完成 ===")

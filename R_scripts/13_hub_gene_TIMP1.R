# ============================================================
# 13_hub_gene_TIMP1.R — Hub gene TIMP1 深度刻画
# ============================================================
# 类比原文 DPEP1 的角色：TIMP1 是本研究的候选 hub gene
# 1) T vs N 表达  2) KM 单基因  3) timeROC  4) 多因素 Cox (adjust age/sex/stage)
# 5) 与 16 免疫细胞类型 Spearman  6) 与 24 检查点相关
# 7) 与 CAF 标志 (ACTA2/FAP/PDGFRB/COL1A1) 相关
# 8) 与临床分期关联  9) GSE39582 外部单基因验证
# ============================================================

suppressPackageStartupMessages({ library(magrittr); library(dplyr); library(data.table); library(GEOquery)
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

HUB <- "TIMP1"

# ---- 1. TCGA 数据 ----
expr_dt <- data.table::fread(file.path(DATA_RAW, "TCGA_COADREAD_TPM_matrix.tsv"))
gene_ids <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_ids; rm(expr_dt, gene_ids)
colnames(expr_mat) <- sub("^\\.\\d+", "", colnames(expr_mat))

sinfo <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_sample_types.csv"), stringsAsFactors = FALSE)
clin  <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_clinical.csv"), stringsAsFactors = FALSE)
rownames(clin) <- clin$patient_barcode

pt <- sinfo %>% filter(sample_type == "Primary Tumor") %>%
  filter(sample_id %in% colnames(expr_mat))
pt_dedup <- pt %>% group_by(patient_barcode) %>%
  summarise(sample_id = min(sample_id), .groups = "drop") %>% as.data.frame()
pt_dedup <- pt_dedup %>% filter(!is.na(patient_barcode))

expr_pt_log <- log2(expr_mat[, pt_dedup$sample_id, drop = FALSE] + 1)
rownames(expr_pt_log) <- rownames(expr_mat)
hub_expr <- expr_pt_log[HUB, ]
names(hub_expr) <- pt_dedup$patient_barcode

# normal 组（每患者一个）
nm <- sinfo %>% filter(sample_type == "Solid Tissue Normal") %>%
  filter(sample_id %in% colnames(expr_mat))
expr_nm_log <- log2(expr_mat[, nm$sample_id, drop = FALSE] + 1)
rownames(expr_nm_log) <- rownames(expr_mat)
hub_norm <- as.numeric(expr_nm_log[HUB, ])

# ---- 2. Fig 6A: Tumor vs Normal ----
df_tn <- data.frame(expression = c(hub_norm, as.numeric(hub_expr)),
                    group = factor(c(rep("Normal", length(hub_norm)),
                                     rep("Tumor", length(hub_expr))),
                                   levels = c("Normal","Tumor")))
p_tn <- ggplot(df_tn, aes(group, expression, fill = group)) +
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c("Tumor" = unname(COL_GROUP["Tumor"]),
                               "Normal" = unname(COL_GROUP["Normal"]))) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 5) +
  labs(x = NULL, y = expression ~ log[2] * "(TPM+1)",
       title = sprintf("%s: tumour vs normal", HUB),
       subtitle = "TCGA-COADREAD") +
  theme_classic(base_size = 12) + theme(legend.position = "none")
save_fig(p_tn, "Fig6A_TIMP1_TvsN", w = 4.5, h = 5)
wilcox_p <- wilcox.test(hub_norm, as.numeric(hub_expr))$p.value
message(sprintf("✓ %s T vs N: Wilcoxon p = %.2e", HUB, wilcox_p))

# ---- 3. 单基因 KM ----
clin_t <- clin[pt_dedup$patient_barcode, ]
os_time   <- clin_t$os_days
os_status <- ifelse(clin_t$vital_status == "Dead", 1,
             ifelse(clin_t$vital_status == "Alive", 0, NA))
valid <- !is.na(os_time) & os_time > 0 & !is.na(os_status)
h <- hub_expr[valid]; t <- os_time[valid]; e <- os_status[valid]
grp_h <- ifelse(h > median(h), "High", "Low")
df_km <- data.frame(time = t, status = e,
                    group = factor(grp_h, levels = c("Low","High")))
fit_km <- survfit(Surv(time, status) ~ group, data = df_km)
p_km <- ggsurvplot(fit_km, data = df_km, pval = TRUE, pval.method = TRUE,
                   risk.table = TRUE, risk.table.height = 0.25,
                   palette = KM_PALETTE,
                   legend.title = "Risk", legend.labs = KM_LEGEND,
                   xlab = "Days", ylab = "Overall Survival",
                   title = sprintf("%s (single gene) — TCGA-COADREAD", HUB),
                   ggtheme = theme_classic(base_size = 12))
png(file.path(FIG_DIR, "Fig6B_TIMP1_KM_TCGA.png"),
    width = 7, height = 7, units = "in", res = FIG_DPI); print(p_km); dev.off()
km_p <- 1 - pchisq(survdiff(Surv(time, status) ~ group, data = df_km)$chisq, 1)
message(sprintf("✓ %s KM log-rank p = %.2e", HUB, km_p))

# ---- 4. timeROC ----
roc <- timeROC(T = t, delta = e, marker = as.numeric(h),
               cause = 1, weighting = "marginal",
               times = c(365, 3*365, 5*365), ROC = TRUE)
p_roc <- ggplot(data.frame(
  FPR = c(roc$FP[,1], roc$FP[,2], roc$FP[,3]),
  TPR = c(roc$TP[,1], roc$TP[,2], roc$TP[,3]),
  Time = rep(c("1-year","3-year","5-year"), each = nrow(roc$FP))
), aes(FPR, TPR, color = Time)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#1B998B","#F46036","#D7263D")) +
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("AUC:\n 1y = %.3f\n 3y = %.3f\n 5y = %.3f",
                           roc$AUC[1], roc$AUC[2], roc$AUC[3]), size = 4) +
  labs(title = sprintf("Time-dependent ROC: %s (TCGA)", HUB),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_classic(base_size = 12)
save_fig(p_roc, "Fig6C_TIMP1_ROC_TCGA", w = 6, h = 5)
message(sprintf("✓ %s AUC 1y=%.3f 3y=%.3f 5y=%.3f", HUB,
                roc$AUC[1], roc$AUC[2], roc$AUC[3]))

# ---- 5. 多因素 Cox（调整 age/sex/stage）----
age_days <- clin_t$age_days[valid]
sex <- factor(ifelse(clin_t$gender[valid] == "male", "Male",
                     ifelse(clin_t$gender[valid] == "female", "Female", NA)))
stage_raw <- clin_t$stage[valid]
# 用 00_config.R 的统一归类（长模式优先），避免 "Stage IIIA" 被并入 I-II
stage_grp <- classify_stage(stage_raw)
df_adj <- data.frame(time = t, status = e,
                     TIMP1 = as.numeric(h),
                     age_days = age_days, sex = sex, stage = stage_grp)
# stage "Unknown" 不是真实临床类别，纳入模型会产生虚假显著 HR；
# 相对分期已知者属缺失数据，应从调整模型中剔除。
df_adj <- df_adj %>% filter(!is.na(age_days), !is.na(sex), stage != "Unknown")
df_adj$stage <- droplevels(df_adj$stage)
message(sprintf("Adjusted Cox set: N=%d, events=%d (excluded %d with unknown stage)",
                nrow(df_adj), sum(df_adj$status), sum(stage_grp == "Unknown")))
fit_adj <- coxph(Surv(time, status) ~ TIMP1 + age_days + sex + stage, data = df_adj)
s_adj <- summary(fit_adj)
cox_tab <- as.data.frame(s_adj$coefficients)
cox_tab <- cbind(term = rownames(cox_tab), cox_tab, row.names = NULL)
save_table(cox_tab, "Table6_TIMP1_multivariate_cox")
print(s_adj)

# ---- 6. 与免疫细胞浸润 / 检查点 关联（复用 09 结果）----
imm <- read.csv(file.path(TAB_DIR, "Table5_immune_infiltration.csv"),
                stringsAsFactors = FALSE)
rownames(imm) <- imm$patient
cell_cols <- setdiff(colnames(imm), c("patient","group","FARGs_score"))
hub_aligned <- imm[[1]]  # patient UUIDs from immune infiltration
# 需要 TIMP1 表达按 UUID 索引（immune_df 用 UUID）
uuid_to_pat <- setNames(pt_dedup$patient_barcode, pt_dedup$sample_id)
hub_by_uuid <- hub_expr[uuid_to_pat[rownames(imm)]]
hub_by_uuid <- as.numeric(hub_by_uuid)

imm_cor <- t(sapply(cell_cols, function(ct) {
  r <- suppressWarnings(cor.test(imm[[ct]], hub_by_uuid, method = "spearman"))
  c(rho = unname(r$estimate), p = r$p.value)
}))
imm_cor_df <- data.frame(cell_type = rownames(imm_cor), imm_cor) %>%
  arrange(desc(rho)) %>%
  mutate(sig = ifelse(p < 0.001,"***", ifelse(p < 0.01,"**",
                        ifelse(p < 0.05, "*", "ns"))))
save_table(imm_cor_df, "Table6b_TIMP1_immune_celltype")
p_imm <- ggplot(imm_cor_df, aes(reorder(cell_type, rho), rho, fill = rho)) +
  geom_col() + coord_flip() +
  geom_text(aes(label = sig), hjust = ifelse(imm_cor_df$rho > 0, -0.1, 1.1), size = 3) +
  scale_fill_gradient2(low = "#1B998B", mid = "white", high = "#D7263D", midpoint = 0) +
  labs(x = NULL, y = "Spearman rho with TIMP1",
       title = sprintf("%s vs immune cell infiltration", HUB)) +
  theme_classic(base_size = 11) + theme(legend.position = "none")
save_fig(p_imm, "Fig6D_TIMP1_immune", w = 6, h = 6)

# ---- 7. TIMP1 vs CAF markers ----
caf_markers <- c("ACTA2","FAP","PDGFRB","COL1A1","COL1A2","FN1","POSTN","SPARC","LUM","TAGLN")
caf_avail <- intersect(caf_markers, rownames(expr_mat))
caf_cor <- t(sapply(caf_avail, function(g) {
  v <- as.numeric(expr_mat[g, pt_dedup$sample_id])
  r <- suppressWarnings(cor.test(v, as.numeric(hub_expr), method = "spearman"))
  c(rho = unname(r$estimate), p = r$p.value)
}))
caf_df <- data.frame(marker = rownames(caf_cor), caf_cor) %>%
  arrange(desc(rho)) %>%
  mutate(sig = ifelse(p < 0.001,"***", ifelse(p < 0.01,"**",
                        ifelse(p < 0.05, "*", "ns"))))
save_table(caf_df, "Table6c_TIMP1_CAF_markers")
message(sprintf("✓ TIMP1 correlates with CAF markers: %d/%d significant (p<0.05)",
                sum(caf_df$p < 0.05), nrow(caf_df)))

p_caf <- ggplot(caf_df, aes(reorder(marker, rho), rho, fill = rho)) +
  geom_col() + coord_flip() +
  geom_text(aes(label = sig), hjust = -0.1, size = 3) +
  scale_fill_gradient2(low = "#1B998B", mid = "white", high = "#D7263D", midpoint = 0) +
  labs(x = NULL, y = "Spearman rho with TIMP1",
       title = sprintf("%s vs CAF/stromal markers", HUB)) +
  theme_classic(base_size = 11) + theme(legend.position = "none")
save_fig(p_caf, "Fig6E_TIMP1_CAF", w = 5, h = 5)

# ---- 7b. Composite CAF-burden score vs TIMP1 (single 'CAF burden' rho) ----
# z-score each available CAF marker across primary-tumour patients, average into
# one CAF-burden score per patient, then Spearman-correlate with patient-level TIMP1.
caf_z <- t(apply(expr_pt_log[caf_avail, , drop = FALSE], 1, function(v) (v - mean(v)) / sd(v)))  # markers x patients
caf_burden <- colMeans(caf_z)                                                                       # per patient_barcode (colnames = sample_id)
names(caf_burden) <- pt_dedup$patient_barcode[match(colnames(caf_z), pt_dedup$sample_id)]
common <- intersect(names(hub_expr), names(caf_burden))
caf_ct <- suppressWarnings(cor.test(as.numeric(hub_expr[common]), as.numeric(caf_burden[common]),
                                    method = "spearman"))
caf_burden_tbl <- data.frame(measure = "TIMP1 vs composite CAF-burden score (mean z of 10 markers)",
                             rho = unname(caf_ct$estimate), p = caf_ct$p.value, n = length(common))
save_table(caf_burden_tbl, "Table6c2_TIMP1_CAF_burden_composite")
message(sprintf("✓ Composite CAF-burden vs TIMP1: Spearman rho = %.3f, p = %.3g, n = %d",
                unname(caf_ct$estimate), caf_ct$p.value, length(common)))

# ---- 8. TIMP1 vs 临床分期 ----
# 用全部有表达的患者（与 hub_expr 等长），不受 OS 有效性限制
stage_full <- clin_t$stage
stage_grp_full <- classify_stage(stage_full)
df_stage <- data.frame(hub = as.numeric(hub_expr), stage = stage_grp_full) %>%
  filter(stage != "Unknown")
p_stage <- ggplot(df_stage, aes(stage, hub, fill = stage)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "YlOrRd") +
  stat_compare_means(method = "kruskal.test", label.y = max(df_stage$hub)*1.05) +
  labs(x = NULL, y = sprintf("%s log2(TPM+1)", HUB),
       title = sprintf("%s expression across AJCC stage", HUB)) +
  theme_classic(base_size = 11) + theme(legend.position = "none")
save_fig(p_stage, "Fig6F_TIMP1_stage", w = 5, h = 5)

# ---- 9. GSE39582 外部单基因验证 ----
cat("\n############ GSE39582 external TIMP1 validation ############\n")
g <- getGEO(filename = file.path(DATA_RAW, "GSE39582_series_matrix.txt.gz"), Annot = FALSE)
e1 <- if (is.list(g)) g[[1]] else g
expr_p <- exprs(e1); pd <- as.data.frame(pData(e1))
pc <- function(x) as.character(sub("^[^:]*:[[:space:]]*", "", x))
pd$os_time   <- as.numeric(pc(pd[["characteristics_ch1.14"]])) * 30.44
pd$os_status <- as.numeric(pc(pd[["characteristics_ch1.13"]]))
pd$dataset   <- pc(pd[["characteristics_ch1.1"]])
pd$tumor     <- pd$dataset %in% c("discovery","validation")
map_all <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr_p),
                                 keytype = "PROBEID", columns = "SYMBOL")
pr <- map_all$PROBEID[map_all$SYMBOL == HUB]
pr <- pr[pr %in% rownames(expr_p)]
message(sprintf("TIMP1 probes on GPL570: %d", length(pr)))
if (length(pr) > 1) {
  tv <- which(pd$tumor & !is.na(pd$os_time) & pd$os_time > 0)
  iqr <- apply(expr_p[pr, tv, drop=FALSE], 1, IQR)
  pr <- pr[which.max(iqr)]
}
hub_gse <- as.numeric(expr_p[pr, ])
v2 <- pd$tumor & !is.na(pd$os_time) & pd$os_time > 0 & !is.na(pd$os_status)
df2 <- data.frame(time = pd$os_time[v2], status = pd$os_status[v2],
                  TIMP1 = hub_gse[v2])
grp2 <- ifelse(df2$TIMP1 > median(df2$TIMP1), "High", "Low")
df2$group <- factor(grp2, levels = c("Low","High"))
fit_km2 <- survfit(Surv(time, status) ~ group, data = df2)
p_km2 <- ggsurvplot(fit_km2, data = df2, pval = TRUE, pval.method = TRUE,
                    risk.table = TRUE, risk.table.height = 0.25,
                    palette = KM_PALETTE,
                    legend.title = "Risk", legend.labs = KM_LEGEND,
                    xlab = "Days", ylab = "Overall Survival",
                    title = sprintf("%s (single gene) — GSE39582 external", HUB),
                    ggtheme = theme_classic(base_size = 12))
png(file.path(FIG_DIR, "Fig6G_TIMP1_KM_GSE39582.png"),
    width = 7, height = 7, units = "in", res = FIG_DPI); print(p_km2); dev.off()
km_p2 <- 1 - pchisq(survdiff(Surv(time, status) ~ group, data = df2)$chisq, 1)
cox_h <- summary(coxph(Surv(time, status) ~ TIMP1, data = df2))
message(sprintf("✓ GSE39582 TIMP1 KM log-rank p = %.2e, Cox HR = %.3f (p=%.3f)",
                km_p2,
                cox_h$coefficients[1, "exp(coef)"],
                cox_h$coefficients[1, "Pr(>|z|)"]))

# ---- 汇总表 ----
# NOTE: timeROC() with a *scalar* `times` returns AUC for c(t=0, t=times), so
# $AUC[1] is the degenerate t=0 value (always NA). Always pass the 3-element
# vector and index 1:3 — this bug previously produced all-NA external AUCs.
roc_ext <- timeROC(T = df2$time, delta = df2$status, marker = df2$TIMP1,
                   cause = 1, weighting = "marginal",
                   times = c(365, 3*365, 5*365))
stopifnot(length(roc_ext$AUC) == 3, length(roc$AUC) == 3)
summary_hub <- data.frame(
  Cohort  = c("TCGA-COADREAD","GSE39582 (external)"),
  N       = c(sum(valid), nrow(df2)),
  Events  = c(sum(e), sum(df2$status)),
  KM_p    = c(km_p, km_p2),
  AUC_1y  = c(roc$AUC[1],     roc_ext$AUC[1]),
  AUC_3y  = c(roc$AUC[2],     roc_ext$AUC[2]),
  AUC_5y  = c(roc$AUC[3],     roc_ext$AUC[3])
)
save_table(summary_hub, "Table6d_TIMP1_hub_summary")
message("=== 13_hub_gene_TIMP1.R 完成 ===")

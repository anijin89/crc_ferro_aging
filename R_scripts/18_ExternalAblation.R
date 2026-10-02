# ============================================================
# 18_ExternalAblation.R — 并集 vs 单臂 ablation 在外部队列 GSE39582 复跑
# ============================================================
# 目的：证明"ferro-aging 预后权重由衰老/SASP 臂承载、铁死亡臂无鉴别"跨队列成立。
# 做法：用 TCGA 训练好的各臂系数，套到 GSE39582 的 per-gene z-score 表达上：
#   - Union / Senescence：同一 5 基因模型 β（训练队列内监督选择结果，二者相同）
#   - Ferroptosis：训练期惩罚化 Cox(LASSO, 无 p 门) 的冻结系数（与监督臂同一估计量）
# 输出：Table11_external_ablation.csv, Fig12A_Ext_KM_3panels.png, Fig12B_Ext_ROC_3panels.png
# ============================================================

suppressPackageStartupMessages({ library(magrittr); library(dplyr); library(GEOquery)
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

# ---- 训练模型（并集=衰老臂 5 基因 β；铁死亡臂=训练期惩罚化 Cox 系数）----
model   <- readRDS(file.path(RES_DIR, "FARGs_score_TCGA.rds"))
union_genes <- model$genes          # CDKN2A, CHEK1, MMP1, MMP3, TIMP1
union_betas <- model$betas_multi
uva     <- readRDS(file.path(RES_DIR, "union_vs_arms.rds"))
ferro_model <- uva$models$Ferroptosis       # list(genes, betas) from LASSO-Cox (no gate)
ferroptosis_genes <- ferro_model$genes
ferroptosis_betas <- ferro_model$betas
message("Frozen ferroptosis-arm penalized Cox model: ",
        paste(sprintf("%s(%+.4f)", ferroptosis_genes, ferroptosis_betas), collapse=", "))
all_pool <- c(union_genes, ferroptosis_genes)   # 需外部探针映射的训练期模型基因

# ---- 载入 GSE39582 ----
g <- getGEO(filename = file.path(DATA_RAW, "GSE39582_series_matrix.txt.gz"), Annot = FALSE)
e1 <- if (is.list(g)) g[[1]] else g
expr_p <- exprs(e1)
pd <- as.data.frame(pData(e1))
parse_char <- function(x) as.character(sub("^[^:]*:\\s*", "", x))
pd$os_event_raw <- parse_char(pd[["characteristics_ch1.13"]])
pd$os_delay_raw <- parse_char(pd[["characteristics_ch1.14"]])
pd$dataset      <- parse_char(pd[["characteristics_ch1.1"]])
pd$os_time   <- suppressWarnings(as.numeric(pd$os_delay_raw)) * 30.44
pd$os_status <- suppressWarnings(as.numeric(pd$os_event_raw))
pd$tumor     <- pd$dataset %in% c("discovery", "validation")
valid <- pd$tumor & !is.na(pd$os_time) & pd$os_time > 0 & !is.na(pd$os_status)
message(sprintf("GSE39582: tumor+valid OS = %d, events = %d", sum(valid), sum(pd$os_status[valid])))

# ---- 探针→symbol：maxIQR（与 08 一致）----
map_all <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(expr_p),
                                 keytype = "PROBEID", columns = "SYMBOL")
map_all <- map_all[!is.na(map_all$SYMBOL) & map_all$PROBEID %in% rownames(expr_p), ]
tumor_cols <- which(valid)
best_probe <- sapply(all_pool, function(gn) {
  pr <- map_all$PROBEID[map_all$SYMBOL == gn]
  if (length(pr) == 0) return(NA_character_)
  if (length(pr) == 1) return(pr)
  iqr <- apply(expr_p[pr, tumor_cols, drop = FALSE], 1, IQR)
  pr[which.max(iqr)]
})
best_probe <- best_probe[!is.na(best_probe)]
message(sprintf("Pool genes %d → probes %d (missing: %s)", length(all_pool), length(best_probe),
                paste(setdiff(all_pool, names(best_probe)), collapse=", ")))
expr_g <- expr_p[best_probe, , drop = FALSE]
rownames(expr_g) <- names(best_probe)

# ---- per-gene z-score（外部队列内标准化，与 08 一致）----
z_gene <- function(m) {
  # m: genes × samples. 逐基因(逐行)标准化：mu/sdv 长度 = nrow(m)，
  # 直接 (m - mu)/sdv 按行对齐（勿用 t() 转置，否则按列回收会错）。
  mu  <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, function(v) { s <- sd(v, na.rm = TRUE); if (is.na(s) || s == 0) 1 else s })
  (m - mu) / sdv
}
ext_z <- z_gene(as.matrix(expr_g))
ext_time   <- pd$os_time[valid]
ext_status <- pd$os_status[valid]

# ---- 各臂评分（套训练系数；z 已按基因名索引，缺基因则跳过）----
score_arm_supervised <- function(genes, betas) {
  g <- intersect(genes, rownames(ext_z))
  if (length(g) < length(genes)) message("  ⚠ 外部缺基因: ", paste(setdiff(genes,g), collapse=", "))
  as.numeric(t(betas[g]) %*% ext_z[g, valid, drop = FALSE])
}
score_arm_equalweight <- function(genes) {
  g <- intersect(genes, rownames(ext_z))
  as.numeric(colSums(ext_z[g, valid, drop = FALSE]))
}

ARMS <- list(
  Union       = list(type="sup", genes=union_genes, betas=union_betas),
  Senescence  = list(type="sup", genes=union_genes, betas=union_betas),
  Ferroptosis = list(type="sup", genes=ferroptosis_genes, betas=ferroptosis_betas)
)

eval_arm <- function(score, label) {
  cidx <- survConcordance(Surv(ext_time, ext_status) ~ score)$concordance
  roc  <- timeROC(T=ext_time, delta=ext_status, marker=score, cause=1,
                  weighting="marginal", times=c(365,3*365,5*365), ROC=TRUE, iid=FALSE)
  grp  <- ifelse(score > median(score), "High","Low")
  dfs  <- data.frame(time=ext_time, status=ext_status, group=factor(grp, levels=c("Low","High")))
  sdf  <- survdiff(Surv(time,status)~group, data=dfs)
  kmp  <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
  cox  <- summary(coxph(Surv(time,status)~score, data=data.frame(time=ext_time,status=ext_status,score=score)))
  hr   <- cox$coefficients[1,"exp(coef)"]; hp <- cox$coefficients[1,"Pr(>|z|)"]
  fit_km <- survfit(Surv(time,status)~group, data=dfs)
  list(label=label, cindex=cidx, auc1=roc$AUC[1], auc3=roc$AUC[2], auc5=roc$AUC[3],
       km_p=kmp, hr=hr, hr_p=hp, km_fit=fit_km, km_df=dfs, roc=roc)
}

res <- lapply(names(ARMS), function(nm) {
  a <- ARMS[[nm]]
  sc <- if (a$type=="sup") score_arm_supervised(a$genes, a$betas) else score_arm_equalweight(a$genes)
  r <- eval_arm(sc, nm)
  message(sprintf("%-12s C=%.3f  AUC1/3/5=%.3f/%.3f/%.3f  KM p=%.3f  HR/SD=%.2f (p=%.3f)",
                  nm, r$cindex, r$auc1, r$auc3, r$auc5, r$km_p, r$hr, r$hr_p))
  r
})
names(res) <- names(ARMS)

# ---- 对比表 ----
cmp <- do.call(rbind, lapply(res, function(r) data.frame(
  Cohort="GSE39582 (external)", Signature=r$label,
  C_index=round(r$cindex,4), AUC_1y=round(r$auc1,4), AUC_3y=round(r$auc3,4),
  AUC_5y=round(r$auc5,4), KM_logrank_p=signif(r$km_p,3),
  Cox_HR_per_SD=round(r$hr,3), Cox_p=signif(r$hr_p,3), stringsAsFactors=FALSE)))
save_table(cmp, "Table11_external_ablation")
print(cmp)

# ---- 图 A：三面 KM（外部）----
km_plots <- lapply(res, function(r) ggsurvplot(
  r$km_fit, data=r$km_df, pval=TRUE, pval.method=TRUE, risk.table=FALSE, legend.title=NULL,
  palette=c(unname(COL_RISK["Low"]), unname(COL_RISK["High"])),
  xlab="Days", ylab="OS probability", title=sprintf("%s — GSE39582", r$label),
  ggtheme=theme_classic(base_size=11))$plot)
p_km <- patchwork::wrap_plots(km_plots, ncol=3) +
  patchwork::plot_annotation(title="External validation (GSE39582): union vs single-programme arms",
                             theme=theme(plot.title=element_text(face="bold",size=12)))
save_fig(p_km, "Fig12A_Ext_KM_3panels", w=12, h=4.3)

# ---- 图 B：三面 ROC（外部）----
roc_plots <- lapply(res, function(r) {
  df <- data.frame(FPR=c(r$roc$FP[,1],r$roc$FP[,2],r$roc$FP[,3]),
                   TPR=c(r$roc$TP[,1],r$roc$TP[,2],r$roc$TP[,3]),
                   Time=rep(c("1-y","3-y","5-y"), each=nrow(r$roc$FP)))
  ggplot(df, aes(FPR,TPR,color=Time)) + geom_line(linewidth=1) +
    geom_abline(slope=1,intercept=0,linetype="dashed",color="grey60") +
    scale_color_manual(values=c("#1B998B","#F46036","#D7263D")) +
    annotate("text",x=0.55,y=0.18,hjust=0,size=3.4,
             label=sprintf("AUC\n 1y=%.3f\n 3y=%.3f\n 5y=%.3f",r$auc1,r$auc3,r$auc5)) +
    coord_equal(xlim=c(0,1),ylim=c(0,1)) +
    labs(title=sprintf("%s (C=%.3f)",r$label,r$cindex),x="1 - Specificity",y="Sensitivity") +
    theme_classic(base_size=11)
})
p_roc <- patchwork::wrap_plots(roc_plots, ncol=3) +
  patchwork::plot_annotation(title="External time-dependent ROC (GSE39582)",
                             theme=theme(plot.title=element_text(face="bold",size=12)))
save_fig(p_roc, "Fig12B_Ext_ROC_3panels", w=12, h=4.6)

save_obj(list(external=cmp, models=lapply(res,function(r) r[c("label","cindex","auc1","auc3","auc5","km_p","hr","hr_p")])),
         "external_ablation")
message("=== 18_ExternalAblation.R 完成 ===")

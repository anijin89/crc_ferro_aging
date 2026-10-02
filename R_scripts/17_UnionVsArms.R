# ============================================================
# 17_UnionVsArms.R — 铁-衰老并集 vs 单臂（铁死亡 / 衰老）鉴别力对比
# ============================================================
# 背景：FARG 定义 = 铁死亡 ∪ 衰老（Cell Metab 2026 概念）。
# 审稿命门：并集签名相比任一单程序多提供了什么？
# 做法：同一套 univariate Cox → LASSO → multivariate Cox → z-score 评分流程，
#       仅更换候选基因池（union / ferroptosis-only / senescence-only），
#       比较 C-index、timeROC AUC(1/3/5y)、KM log-rank。
# 输出：Table10_union_vs_arms.csv, Fig11A_KM_3panels.png, Fig11B_ROC_3panels.png,
#       union_vs_arms.rds
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
suppressPackageStartupMessages({ library(survival); library(timeROC) })

# ---- 1. 表达矩阵 + 样本映射（与 07 完全一致）----
expr_dt <- data.table::fread(file.path(DATA_RAW, "TCGA_COADREAD_TPM_matrix.tsv"))
gene_ids <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_ids; rm(expr_dt, gene_ids)
colnames(expr_mat) <- sub("^\\.\\d+", "", colnames(expr_mat))

sinfo <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_sample_types.csv"), stringsAsFactors = FALSE)
clin  <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_clinical.csv"), stringsAsFactors = FALSE)
rownames(clin) <- clin$patient_barcode
farg_degs <- readRDS(file.path(RES_DIR, "FARG_DEG_intersection.rds"))

# ---- 2. Primary Tumor + 患者级去重 ----
pt <- sinfo %>% filter(sample_type == "Primary Tumor") %>%
  filter(sample_id %in% colnames(expr_mat))
pt_dedup <- pt %>% group_by(patient_barcode) %>%
  summarise(sample_id = min(sample_id), .groups = "drop") %>% as.data.frame()

# ---- 3. 临床 OS ----
clin_sub  <- clin[pt_dedup$patient_barcode, ]
os_time   <- clin_sub$os_days
os_status <- ifelse(clin_sub$vital_status == "Dead", 1,
             ifelse(clin_sub$vital_status == "Alive", 0, NA))
valid <- !is.na(os_time) & os_time > 0 & !is.na(os_status)
samples     <- pt_dedup$sample_id[valid]
os_time     <- os_time[valid]
os_status   <- os_status[valid]
patient_ids <- pt_dedup$patient_barcode[valid]
message(sprintf("✓ N=%d, events=%d", length(os_status), sum(os_status == 1)))

# ============================================================
# 4. 26 个 FARG-DEG 的铁死亡 / 衰老归属（FerrDb / CellAge canonical 注释）
#    脚本内固化，作为 curated assumption；两臂允许生物学重叠但此处按主注释互斥划分。
# ============================================================
ferroptosis_genes <- c("HMOX1","SLC11A2","SLC3A2","SLC7A11","NQO1","SQLE")          # 铁/谷胱甘肽/脂质过氧化核心
senescence_genes  <- c("CDKN1A","CDKN2A","CDKN2B","CCND1","CDK2","CDK4","E2F1",
                       "CHEK1","BMP4",                                                # 细胞周期停滞 / DDR
                       "CXCL8","GDF15","IL1A","IL1B","CSF2","SERPINE1",
                       "MMP1","MMP3","MMP9","TIMP1","VEGFA")                          # SASP

stopifnot(all(ferroptosis_genes %in% farg_degs))
stopifnot(all(senescence_genes  %in% farg_degs))
stopifnot(length(unique(c(ferroptosis_genes, senescence_genes))) == length(farg_degs))

# 归属表（供补表 / 回复信引用）
membership <- data.frame(
  gene        = farg_degs,
  ferroptosis = farg_degs %in% ferroptosis_genes,
  senescence  = farg_degs %in% senescence_genes,
  assignment  = ifelse(farg_degs %in% ferroptosis_genes, "ferroptosis", "senescence"),
  stringsAsFactors = FALSE
)
save_table(membership, "Table10b_gene_arm_membership")

POOLS <- list(
  Union        = farg_degs,
  Ferroptosis  = ferroptosis_genes,
  Senescence   = senescence_genes
)

# ============================================================
# 5. 对任一基因池跑完整建模流程，返回评分 + 鉴别力指标
# ============================================================
build_and_eval <- function(pool_genes, label) {
  # log2 + per-gene z-score（训练队列内，与 07 一致）
  e <- log2(expr_mat[pool_genes, samples, drop = FALSE] + 1)
  mu <- rowMeans(e); sg <- apply(e, 1, sd); sg[sg == 0] <- 1
  e  <- (e - mu) / sg; rownames(e) <- pool_genes

  # 记录每基因单因素 Cox p（仅作描述，不作筛选门槛）
  cox_uni <- do.call(rbind, lapply(rownames(e), function(g) {
    fit <- tryCatch(coxph(Surv(os_time, os_status) ~ as.numeric(e[g, ])), error = function(err) NULL)
    if (is.null(fit)) return(data.frame(gene = g, p = NA_real_, beta = NA_real_))
    s <- summary(fit)
    data.frame(gene = g, p = s$coefficients[1, "Pr(>|z|)"],
               beta = s$coefficients[1, "coef"])
  }))
  n_sig_uni <- sum(cox_uni$p < COX_P_CUTOFF, na.rm = TRUE)

  # ---- 与主模型 07 完全一致的流程：单因素门槛 → LASSO → 多因素 ----
  sig_genes <- cox_uni$gene[which(cox_uni$p < COX_P_CUTOFF)]
  supervised <- TRUE
  fallback_lambda <- FALSE
  method <- "UniCox+LASSO+Cox"
  if (length(sig_genes) >= 2) {
    x <- t(as.matrix(e[sig_genes, , drop = FALSE])); y <- Surv(os_time, os_status)
    set.seed(SEED)
    cv <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10, type.measure = "deviance")
    cf <- coef(cv, s = "lambda.min")
    sel <- rownames(cf)[cf[, 1] != 0]
    if (length(sel) == 0) sel <- sig_genes
    dfm <- as.data.frame(t(as.matrix(e[sel, , drop = FALSE])))
    dfm$os_time <- os_time; dfm$os_status <- os_status
    fitm <- coxph(Surv(os_time, os_status) ~ ., data = dfm[, c(sel, "os_time", "os_status")])
    betas <- summary(fitm)$coefficients[, "coef"]
    score <- as.numeric(t(betas) %*% e[sel, , drop = FALSE])
  } else {
    # 门槛后候选 <2（铁死亡臂）：为与监督臂保持同一估计量（惩罚化 Cox），
    # 不再退回「等权 z-score 和（非监督）」，而是把同一套 LASSO-Cox 直接施加到
    # 整个候选池（跳过单因素 p 门，因为该门对本臂为空）。惩罚本身完成变量选择。
    # 若 λ.min 将全部系数压为 0，则沿 λ 序列取「首个保留 ≥1 个非零系数」的
    # 最小惩罚 λ（least-penalty non-null），并对任何空门臂一视同仁地套用该规则。
    x <- t(as.matrix(e[pool_genes, , drop = FALSE])); y <- Surv(os_time, os_status)
    set.seed(SEED)
    cv <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10, type.measure = "deviance")
    cf <- coef(cv, s = "lambda.min")
    sel <- rownames(cf)[cf[, 1] != 0]; betas <- cf[sel, 1]
    if (length(sel) == 0) {
      fallback_lambda <- TRUE
      for (lam in cv$glmnet.fit$lambda) {
        cc <- coef(cv$glmnet.fit, s = lam); nz <- which(cc[, 1] != 0)
        if (length(nz) >= 1) { sel <- rownames(cc)[nz]; betas <- cc[sel, 1]; break }
      }
    }
    if (length(sel) == 0) {   # 极端情形：全 λ 序列皆零 → 无预后信号
      sel <- pool_genes; betas <- rep(0, length(sel))
    }
    names(betas) <- sel
    score <- as.numeric(t(betas) %*% e[sel, , drop = FALSE])
    method <- "LASSO-Cox (no gate)"
    message("  ⚠ ", label, ": 单因素门槛候选 ", length(sig_genes),
            " 个 → 改用惩罚化 Cox(LASSO, 无 p 门)",
            if (fallback_lambda) "；λ.min 全零，取最小非零惩罚 λ" else "",
            " → 保留 ", length(sel), " 基因 (", paste(sel, collapse = ", "), ")")
  }

  # C-index（survival 一致性）
  cidx <- survConcordance(Surv(os_time, os_status) ~ score)$concordance

  # timeROC 1/3/5y
  roc <- timeROC(T = os_time, delta = os_status, marker = score, cause = 1,
                 weighting = "marginal", times = c(365, 3*365, 5*365), ROC = TRUE, iid = FALSE)

  # KM 中位分组 log-rank
  grp <- ifelse(score > median(score), "High", "Low")
  dfs <- data.frame(time = os_time, status = os_status,
                    group = factor(grp, levels = c("Low","High")))
  sdf <- survdiff(Surv(time, status) ~ group, data = dfs)
  kmp <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  fit_km <- survfit(Surv(time, status) ~ group, data = dfs)

  list(label = label, n_pool = length(pool_genes), n_sig = n_sig_uni, n_sel = length(sel),
       supervised = supervised, method = method,
       genes = sel, betas = betas, score = score, patient_ids = patient_ids,
       cindex = cidx, auc1 = roc$AUC[1], auc3 = roc$AUC[2], auc5 = roc$AUC[3],
       km_p = kmp, km_fit = fit_km, km_df = dfs, roc = roc,
       uni = cox_uni)
}

res <- lapply(names(POOLS), function(nm) {
  message(sprintf("\n===== 建模: %s (pool n=%d) =====", nm, length(POOLS[[nm]])))
  r <- build_and_eval(POOLS[[nm]], nm)
  if (!is.null(r)) message(sprintf("  C=%.3f  AUC1=%.3f AUC3=%.3f AUC5=%.3f  KMrank p=%.3g  (sig %d → LASSO %d)",
                                   r$cindex, r$auc1, r$auc3, r$auc5, r$km_p, r$n_sig, r$n_sel))
  r
})
names(res) <- names(POOLS)
res <- res[!vapply(res, is.null, logical(1))]

# 铁死亡臂：6 个经典执行基因的单因素 Cox（证明其在本队列无独立 OS 关联）
ef <- log2(expr_mat[ferroptosis_genes, samples, drop = FALSE] + 1)
muf <- rowMeans(ef); sgf <- apply(ef, 1, sd); sgf[sgf == 0] <- 1
ef <- (ef - muf) / sgf; rownames(ef) <- ferroptosis_genes
ferro_uni <- do.call(rbind, lapply(ferroptosis_genes, function(g) {
  fit <- coxph(Surv(os_time, os_status) ~ as.numeric(ef[g, ]))
  s <- summary(fit)
  data.frame(gene = g,
             HR = round(s$coefficients[1,"exp(coef)"],3),
             CI = sprintf("%.2f-%.2f", s$conf.int[1,"lower .95"], s$conf.int[1,"upper .95"]),
             p = signif(s$coefficients[1,"Pr(>|z|)"],3))
}))
save_table(ferro_uni, "Table10d_ferroptosis_arm_univariate_cox")
message("\n--- Ferroptosis arm 单因素 Cox（6 基因）---"); print(ferro_uni)

# ============================================================
# 6. 对比表
# ============================================================
cmp <- do.call(rbind, lapply(res, function(r) data.frame(
  Signature      = r$label,
  Pool_genes     = r$n_pool,
  UniCox_sig     = r$n_sig,
  Model_genes    = r$n_sel,
  Selection      = r$method,
  Genes_in_model = paste(r$genes, collapse = "; "),
  C_index        = round(r$cindex, 4),
  AUC_1y         = round(r$auc1, 4),
  AUC_3y         = round(r$auc3, 4),
  AUC_5y         = round(r$auc5, 4),
  KM_logrank_p   = signif(r$km_p, 3),
  stringsAsFactors = FALSE
)))
save_table(cmp, "Table10_union_vs_arms")
print(cmp)

# ΔC 与 ΔAUC(并集 − 各单臂)
ref <- cmp[cmp$Signature == "Union", ]
delta <- cmp %>% filter(Signature != "Union") %>%
  mutate(Delta_C     = round(ref$C_index - C_index, 4),
         Delta_AUC3y = round(ref$AUC_3y - AUC_3y, 4)) %>%
  select(Signature, Delta_C, Delta_AUC3y)
save_table(delta, "Table10c_union_minus_arm")
message("\n--- Union − single arm ---"); print(delta)

# ============================================================
# 7. 图 A：三面 KM 并排
# ============================================================
km_plots <- lapply(res, function(r) {
  ggsurvplot(r$km_fit, data = r$km_df, pval = TRUE, pval.method = TRUE,
             risk.table = FALSE, legend.title = NULL,
             palette = c(unname(COL_RISK["Low"]), unname(COL_RISK["High"])),
             xlab = "Days", ylab = "OS probability",
             title = sprintf("%s signature", r$label),
             ggtheme = theme_classic(base_size = 11))$plot
})
p_km <- patchwork::wrap_plots(km_plots, ncol = 3) +
  patchwork::plot_annotation(title = "Ferro-aging union vs single-programme arms (TCGA-COADREAD)",
                             theme = theme(plot.title = element_text(face = "bold", size = 12)))
save_fig(p_km, "Fig11A_KM_3panels", w = 12, h = 4.3)

# ============================================================
# 8. 图 B：三面 timeROC 并排（1/3/5y）
# ============================================================
roc_plots <- lapply(res, function(r) {
  df <- data.frame(
    FPR = c(r$roc$FP[, 1], r$roc$FP[, 2], r$roc$FP[, 3]),
    TPR = c(r$roc$TP[, 1], r$roc$TP[, 2], r$roc$TP[, 3]),
    Time = rep(c("1-y","3-y","5-y"), each = nrow(r$roc$FP)))
  ggplot(df, aes(FPR, TPR, color = Time)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_manual(values = c("#1B998B","#F46036","#D7263D")) +
    annotate("text", x = 0.55, y = 0.18, hjust = 0, size = 3.4,
             label = sprintf("AUC\n 1y=%.3f\n 3y=%.3f\n 5y=%.3f", r$auc1, r$auc3, r$auc5)) +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    labs(title = sprintf("%s (C=%.3f)", r$label, r$cindex),
         x = "1 - Specificity", y = "Sensitivity") +
    theme_classic(base_size = 11)
})
p_roc <- patchwork::wrap_plots(roc_plots, ncol = 3) +
  patchwork::plot_annotation(title = "Time-dependent ROC: union vs ferroptosis-only vs senescence-only",
                             theme = theme(plot.title = element_text(face = "bold", size = 12)))
save_fig(p_roc, "Fig11B_ROC_3panels", w = 12, h = 4.6)

# 存对象
save_obj(list(comparison = cmp, delta = delta, membership = membership,
              models = lapply(res, function(r) r[c("label","genes","betas","score","patient_ids","cindex","auc1","auc3","auc5","km_p")])),
         "union_vs_arms")

message("\n=== 17_UnionVsArms.R 完成 ===")

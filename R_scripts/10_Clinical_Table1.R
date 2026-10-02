# ============================================================
# 10_Clinical_Table1.R — 按 FARGs 风险分组的临床基线特征
#   Table 1: 队列基线（High vs Low risk group），组间差异检验
#   Table 1b: 风险分组与临床病理特征的关联（列联 + chi-square）
#   Fig:     风险评分在不同临床亚组中的分布
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
suppressPackageStartupMessages({ library(dplyr); library(survival) })

# ---- 1. 载入评分对象 + 临床 ----
score_obj <- readRDS(file.path(RES_DIR, "FARGs_score_TCGA.rds"))
clin <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_clinical.csv"), stringsAsFactors = FALSE)
rownames(clin) <- clin$patient_barcode

df <- data.frame(
  patient   = score_obj$patient_ids,
  score     = score_obj$score,
  os_time   = score_obj$os_time,
  os_status = score_obj$os_status
)
cl <- clin[df$patient, ]
df$age_years <- cl$age_days / 365.25
df$sex      <- factor(cl$gender, levels = c("female","male"))
# stage 归组（用 00_config.R 的统一归类；Unknown 视作缺失，与原行为一致）
stage_full  <- cl$stage
df$stage_grp <- classify_stage(stage_full)
df$stage_grp[df$stage_grp == "Unknown"] <- NA
df$stage_grp <- droplevels(df$stage_grp)
# TNM 归组
df$T_grp <- factor(dplyr::case_when(grepl("^T4",cl$T)~"T4", grepl("^T3",cl$T)~"T3",
             grepl("^T2",cl$T)~"T2", grepl("^T1",cl$T)~"T1", TRUE~"Unknown"),
             levels=c("T1","T2","T3","T4"))
df$N_grp <- factor(dplyr::case_when(grepl("^N2",cl$N)~"N2", grepl("^N1",cl$N)~"N1",
             grepl("^N0",cl$N)~"N0", TRUE~"Unknown"), levels=c("N0","N1","N2"))
df$M_grp <- factor(dplyr::case_when(grepl("^M1",cl$M)~"M1", grepl("^M0",cl$M)~"M0",
             TRUE~"Unknown"), levels=c("M0","M1"))

# 中位切分风险组
df$group <- factor(ifelse(df$score > median(df$score), "High", "Low"),
                   levels = c("Low","High"))

# 分析集：stage 已知的样本（Table1 用全队列，关联分析仅可分类者）
cat(sprintf("Total patients with score: %d (High=%d, Low=%d)\n",
            nrow(df), sum(df$group=="High"), sum(df$group=="Low")))
save_obj(df, "clinical_table1_df")

# ---- 2. 通用统计函数 ----
summ_cont <- function(x, grp) {
  # mean(SD) per group + Wilcoxon rank-sum
  lv <- levels(grp); res <- list()
  for (g in lv) res[[g]] <- x[!is.na(x) & !is.na(grp) & grp == g]
  m <- sapply(res, function(v) if(length(v)) sprintf("%.1f (%.1f)", mean(v), sd(v)) else "NA")
  p_w  <- suppressWarnings(tryCatch(wilcox.test(res[[1]], res[[2]])$p.value, error=function(e) NA))
  p_str <- if (is.na(p_w)) "NA" else if (p_w < 0.001) "<0.001" else sprintf("%.3f", p_w)
  data.frame(Variable="Age (years)",
             Low=m[[1]], High=m[[2]],
             Overall=sprintf("%.1f (%.1f)", mean(x,na.rm=TRUE), sd(x,na.rm=TRUE)),
             Test="Wilcoxon", p.value = p_str,
             n1=length(res[[1]]), n2=length(res[[2]]),
             check.rows=FALSE, stringsAsFactors=FALSE)
}
summ_cat <- function(var, name, grp) {
  lv <- levels(factor(var)); res <- list(); pvals <- c()
  tab <- table(var, grp, useNA="no")
  # 去除全为 Unknown 的行
  tab <- tab[rowSums(tab) > 0, , drop=FALSE]
  pchisq <- tryCatch(chisq.test(tab)$p.value, warning=function(w) NA, error=function(e) NA)
  # 若期望频数过小用 Fisher（仅 2x2）
  if (all(dim(tab) == c(2,2))) {
    pf <- tryCatch(fisher.test(tab)$p.value, error=function(e) NA)
    if (!is.na(pf)) pchisq <- pf
  }
  out <- data.frame()
  p_str <- if (is.na(pchisq)) "NA" else if (pchisq < 0.001) "<0.001" else sprintf("%.3f", pchisq)
  for (r in rownames(tab)) {
    low_n <- tab[r,"Low"]; high_n <- tab[r,"High"]
    low_pct <- 100*low_n/sum(tab[,"Low"]); high_pct <- 100*high_n/sum(tab[,"High"])
    all_n <- low_n + high_n; all_pct <- 100*all_n/sum(tab)
    out <- rbind(out, data.frame(
      Variable = sprintf("%s: %s", name, r),
      Low = sprintf("%d (%.1f%%)", low_n, low_pct),
      High = sprintf("%d (%.1f%%)", high_n, high_pct),
      Overall = sprintf("%d (%.1f%%)", all_n, all_pct),
      Test = "Chi-square", p.value = p_str,
      n1 = sum(tab[,"Low"]), n2 = sum(tab[,"High"]), check.rows=FALSE, stringsAsFactors=FALSE))
  }
  out
}

# ---- 3. 组装 Table 1 ----
# 每个变量按 complete-case 同时取 var 与 group，避免错位
summ_cat_cc <- function(var, name) {
  ok <- !is.na(var) & !as.character(var) %in% c("Unknown","")
  summ_cat(factor(var[ok]), name, df$group[ok])
}
row_age <- summ_cont(df$age_years, df$group)
r_sex   <- summ_cat_cc(df$sex,   "Sex")
r_stage <- summ_cat_cc(df$stage_grp, "Stage")
r_T     <- summ_cat_cc(df$T_grp, "T stage")
r_N     <- summ_cat_cc(df$N_grp, "N stage")
r_M     <- summ_cat_cc(df$M_grp, "M stage")

tab1 <- rbind(row_age, r_sex, r_stage, r_T, r_N, r_M)
tab1 <- tab1 %>% select(Variable, Low, High, Overall, p.value)
save_table(tab1, "Table1_baseline_by_riskgroup")
print(tab1)

# ---- 4. 风险评分分布随 stage（Fig）----
suppressPackageStartupMessages(library(ggplot2))
p_dist <- ggplot(df %>% filter(stage_grp != "Unknown"),
                 aes(x=stage_grp, y=score, fill=stage_grp)) +
  geom_violin(alpha=0.6, na.rm=TRUE) + geom_boxplot(width=0.12, outlier.shape=NA, na.rm=TRUE) +
  scale_fill_brewer(palette="Set2") +
  labs(x="Stage", y="FARGs risk score",
       title="FARGs risk score across clinicopathological stages") +
  theme_classic(base_size=12) + theme(legend.position="none")
kw <- kruskal.test(score ~ stage_grp, data = df %>% filter(stage_grp!="Unknown"))
p_dist <- p_dist + annotate("text", x=Inf, y=Inf, hjust=1, vjust=1.5,
             label=sprintf("Kruskal-Wallis p = %.3g", kw$p.value))
save_fig(p_dist, "Fig9A_risk_by_stage", w=6, h=5)

# ---- 5. Table1b: 高风险占比随临床特征 ----
tab1b <- df %>% filter(stage_grp != "Unknown") %>%
  group_by(stage_grp) %>%
  summarise(N=n(), High_n=sum(group=="High"),
            High_pct = 100*mean(group=="High"), .groups="drop") %>%
  mutate(across(c(High_pct), ~round(.,1)))
save_table(as.data.frame(tab1b), "Table1b_highrisk_by_stage")
print(tab1b)

cat("\n=== Table 1 完成 ===\n")

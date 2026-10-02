# ============================================================
# 11_Nomogram_DCA.R — 预后列线图 + 校准 + DCA + C-index 增量
#   模型: Cox ~ FARGs风险评分 + 年龄 + 性别 + 分期
#   Fig9B 列线图   Fig9C 校准曲线(3y)   Fig9D DCA   Fig9E C-index
#   Table2: Cox 系数 / DCA / C-index 与 LR 检验
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
suppressPackageStartupMessages({
  library(survival); library(rms); library(dplyr); library(ggplot2)
})

# ---- 1. 数据 ----
df <- readRDS(file.path(RES_DIR, "clinical_table1_df.rds"))
d <- df %>%
  filter(!is.na(stage_grp), !is.na(age_years), !is.na(sex),
         !is.na(os_time), os_time > 0) %>%
  mutate(sex = relevel(factor(sex), ref="female"),
         stage_grp = factor(stage_grp),
         event = ifelse(os_status == 1, 1, 0))
cat(sprintf("Nomogram analysis set: N=%d, Events=%d\n", nrow(d), sum(d$event)))

dd <- datadist(d); options(datadist="dd")

# ---- 2. Cox 模型（rms 用于列线图，coxph 用于推断）----
fit_full <- cph(Surv(os_time, event) ~ score + age_years + sex + stage_grp,
                data=d, x=TRUE, y=TRUE, surv=TRUE)
m_cl_cox   <- coxph(Surv(os_time,event) ~ age_years + sex + stage_grp, data=d)
m_full_cox <- coxph(Surv(os_time,event) ~ score + age_years + sex + stage_grp, data=d)
cat("\n--- Full model coefficients (rms) ---\n"); print(summary(fit_full))

s_full <- summary(m_full_cox)
co_full <- as.data.frame(s_full$coefficients)
co_full$term <- rownames(co_full)
co_full$HR <- co_full$`exp(coef)`
co_full$HR_lower <- exp(co_full$coef - 1.96*co_full$`se(coef)`)
co_full$HR_upper <- exp(co_full$coef + 1.96*co_full$`se(coef)`)
save_table(co_full %>% select(term, coef, `se(coef)`, z, `Pr(>|z|)`, HR, HR_lower, HR_upper),
           "Table2_nomogram_cox")

# ---- 3. Fig9B 列线图（3 年 OS 概率标尺）----
S_full <- survfit(fit_full)          # 含 baseline
# 3 年基线生存
t_cal <- 3*365.25
lp0 <- 0
base3 <- {
  tm <- S_full$time; sv <- S_full$surv
  idx <- max(which(tm <= t_cal)); if (length(idx)==0 || idx==0) 1 else sv[idx]
}
fn3 <- function(x) 1 - base3^exp(x)   # x = linear predictor -> 3y 死亡概率
png(file.path(FIG_DIR, "Fig9B_nomogram.png"), width=10, height=7, units="in", res=300)
plot(nomogram(fit_full, fun=fn3, funlabel="3-year mortality probability"),
     xfrac=0.32, cex.vars=0.9, cex.axis=0.85)
title(main="FARGs-based prognostic nomogram (overall survival)", line=2.5)
dev.off()
message("✓ Fig9B saved")

# ---- 4. 个体 3 年预测死亡概率（用于校准与 DCA）----
lp <- predict(fit_full, d)                       # linear predictor
base <- survfit(fit_full)
S0_3 <- { idx <- max(which(base$time <= t_cal)); if (length(idx)==0||idx==0) 1 else base$surv[idx] }
d$pred3 <- 1 - S0_3^exp(lp)                       # 3 年死亡概率
cat(sprintf("Predicted 3-y mortality: mean=%.3f, range [%.3f, %.3f]\n",
            mean(d$pred3), min(d$pred3), max(d$pred3)))

# ---- 5. Fig9C 校准曲线（5 组，预测 vs KM 观测）----
cal <- d %>% mutate(bin = ntile(pred3, 5)) %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(pred3), n = n(),
            KM_obs = {
              f <- survfit(Surv(os_time, event) ~ 1, data=pick(everything()))
              idx <- max(which(f$time <= t_cal))
              obs <- if (length(idx)==0||idx==0) 0 else 1 - f$surv[idx]
              obs
            }) %>% ungroup()
save_table(cal, "Table2_calibration_3yr")

set.seed(SEED)
boot_cal <- replicate(200, suppressWarnings({
  db <- d[sample(nrow(d)), ]
  db %>% mutate(bin=ntile(pred3,5)) %>% group_by(bin) %>%
    summarise(km = { f<-survfit(Surv(os_time,event)~1,data=pick(everything()))
                     idx<-max(which(f$time<=t_cal)); if(length(idx)==0||idx==0) 0 else 1-f$surv[idx] },
              p  = mean(pred3), .groups="drop") %>% group_by(bin) %>%
    summarise(km=mean(km), lo=quantile(km,0.025), hi=quantile(km,0.975), p=mean(p), .groups="drop")
}), simplify = FALSE)
boot_df <- do.call(rbind, lapply(boot_cal, as.data.frame))
cal_ci <- boot_df %>% group_by(bin) %>%
  summarise(km=mean(km), lo=quantile(km,.025), hi=quantile(km,.975), p=mean(p), .groups="drop")

p_cal <- ggplot(cal_ci, aes(x=p, y=km)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey50") +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.02, color="grey40") +
  geom_point(size=3.2, color="#1B998B") +
  geom_smooth(aes(x=p, y=km), method="lm", se=FALSE, color="#D7263D", linetype="dotted") +
  labs(x="Predicted 3-year mortality", y="Observed (Kaplan-Meier) 3-year mortality",
       title="Calibration of FARGs nomogram (3-year OS)") +
  theme_classic(base_size=12) + xlim(0,1) + ylim(0,1)
ggsave(file.path(FIG_DIR, "Fig9C_calibration.png"), p_cal, width=6, height=5.5, dpi=300)
message("✓ Fig9C saved")

# ---- 6. Fig9D DCA（阈值概率 net benefit）----
tp_grid <- seq(0.02, 0.60, by=0.01)
# KM 基线在该时点总体事件概率
sf_all <- survfit(Surv(os_time, event) ~ 1, data=d)
idx_all <- max(which(sf_all$time <= t_cal))
prev <- if (length(idx_all)==0||idx_all==0) 0 else 1 - sf_all$surv[idx_all]
n <- nrow(d)
nb_model <- sapply(tp_grid, function(t) {
  flag <- d$pred3 > t
  tp <- sum(flag & d$event==1 & d$os_time <= t_cal)
  fp <- sum(flag) - tp
  tp/n - (fp/n) * (t/(1-t))
})
nb_all <- prev - (1-prev) * (tp_grid/(1-tp_grid))
dca <- data.frame(threshold=tp_grid, Model=nb_model, All=nb_all, None=0)
save_table(dca, "Table2b_DCA_3yr")
cat(sprintf("Overall 3-y event prevalence=%.3f; max model net benefit=%.4f at t=%.2f\n",
            prev, max(nb_model, na.rm=TRUE), tp_grid[which.max(nb_model)]))

p_dca <- ggplot(dca, aes(x=threshold)) +
  geom_line(aes(y=Model, color="Nomogram"), linewidth=1) +
  geom_line(aes(y=All,   color="Treat all"), linewidth=0.8, linetype="dashed") +
  geom_hline(yintercept=0, color="grey60") +
  scale_color_manual(values=c("Nomogram"="#D7263D","Treat all"="#3D5A80")) +
  labs(x="Threshold probability", y="Net benefit",
       title="Decision Curve Analysis (3-year OS)") +
  theme_classic(base_size=12) + ylim(-0.05, max(c(nb_model,nb_all), na.rm=TRUE)*1.2)
ggsave(file.path(FIG_DIR, "Fig9D_DCA.png"), p_dca, width=6.5, height=5, dpi=300)
message("✓ Fig9D saved")

# ---- 7. C-index 增量 + LR 检验 ----
c_cl   <- concordance(m_cl_cox)$concordance
c_full <- concordance(m_full_cox)$concordance
c_score<- concordance(coxph(Surv(os_time,event)~score, data=d))$concordance
lr_stat <- as.numeric(2*(logLik(m_full_cox) - logLik(m_cl_cox)))
lr_p <- pchisq(lr_stat, df=1, lower.tail=FALSE)
cat(sprintf("\nLogLik: clinical=%.3f full=%.3f | LR chi2=%.3f p=%.4g\n",
            logLik(m_cl_cox), logLik(m_full_cox), lr_stat, lr_p))

ci_tab <- data.frame(
  Model = c("FARGs score alone", "Clinical only (age+sex+stage)", "Clinical + FARGs score"),
  C_index = round(c(c_score, c_cl, c_full), 4),
  LR_stat = c(NA, NA, round(lr_stat,3)),
  LR_p = c(NA, NA, if (lr_p<0.001) "<0.001" else round(lr_p,4)))
save_table(ci_tab, "Table2c_cindex_comparison")
print(ci_tab)

p_ci <- ggplot(ci_tab, aes(reorder(Model,C_index), C_index, fill=Model)) +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.3f", C_index)), vjust=-0.3, size=4) +
  scale_fill_manual(values=c("#1B998B","#3D5A80","#D7263D")) +
  labs(x=NULL, y="Harrell's C-index",
       title="Discrimination: added value of FARGs score",
       subtitle=sprintf("Complete-case set, n = %d (score-alone C = %.3f)",
                        nrow(d), c_score)) +
  theme_classic(base_size=12) + theme(legend.position="none",
    axis.text.x = element_text(angle=12, hjust=1)) +
  expand_limits(y=0.85)
ggsave(file.path(FIG_DIR, "Fig9E_cindex.png"), p_ci, width=6.5, height=4.5, dpi=300)
message("✓ Fig9E saved")

cat("\n=== Nomogram / DCA / C-index 完成 ===\n")

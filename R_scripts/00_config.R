# ============================================================
# 00_config.R — 全局配置（所有脚本均 source 此文件）
# ============================================================

# ---- 路径（自动定位，不写死本机路径）----
# 优先级：环境变量 CRC_PROJ_ROOT > 调用方传入的 .ROOT > 由本文件位置向上查找
PROJ_ROOT <- Sys.getenv("CRC_PROJ_ROOT", unset = "")
if (!nzchar(PROJ_ROOT) || !dir.exists(PROJ_ROOT)) {
  if (exists(".ROOT", envir = globalenv()) && dir.exists(get(".ROOT", envir = globalenv()))) {
    PROJ_ROOT <- get(".ROOT", envir = globalenv())
  } else {
    .cfg_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
    if (is.na(.cfg_file)) {
      .fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      .cfg_file <- if (length(.fa)) normalizePath(sub("^--file=", "", .fa[1])) else getwd()
    }
    PROJ_ROOT <- dirname(dirname(.cfg_file))
    while (!file.exists(file.path(PROJ_ROOT, "R_scripts", "00_config.R")) &&
           !identical(PROJ_ROOT, dirname(PROJ_ROOT))) PROJ_ROOT <- dirname(PROJ_ROOT)
  }
}
message("PROJ_ROOT = ", PROJ_ROOT)

DATA_RAW   <- file.path(PROJ_ROOT, "data", "raw")
DATA_PROC  <- file.path(PROJ_ROOT, "data", "processed")
RES_DIR    <- file.path(PROJ_ROOT, "results")
FIG_DIR    <- file.path(PROJ_ROOT, "figures")
TAB_DIR    <- file.path(PROJ_ROOT, "tables")
for (d in c(DATA_RAW, DATA_PROC, RES_DIR, FIG_DIR, TAB_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- 全局参数 ----
SEED          <- 20260919
DEG_LOGFC     <- 1
DEG_FDR       <- 0.05
COX_P_CUTOFF  <- 0.05
FIG_DPI       <- 300
FIG_WIDTH     <- 8
FIG_HEIGHT    <- 6

set.seed(SEED)
options(stringsAsFactors = FALSE)
options(timeout = max(600, getOption("timeout")))  # GEO 大文件下载

# ---- 数据集 ----
TCGA_COHORT   <- "TCGA-COADREAD"
GEO_VALIDATION <- c("GSE39582", "GSE17536", "GSE41258")  # bulk 验证
GEO_DEG_VALIDATION <- c("GSE17538", "GSE33113")          # 表达验证（含 normal）
GSE_SCRNA     <- "GSE132465"                              # 单细胞
GSE_SPATIAL   <- "GSE281413"                              # GeoMx DSP 人CRC 肿瘤/正常 分区空间转录组

# ---- FARGs ----
FARG_FILE <- file.path(DATA_PROC, "FARG_95genes.csv")

# ---- 颜色 ----
COL_RISK <- c("High" = "#D7263D", "Low" = "#1B998B")
COL_GROUP <- c("Tumor" = "#D7263D", "Normal" = "#3D5A80")
# survminer 的 palette 只有在 strata 标签与键名完全一致时才按名匹配；
# 默认 strata 名为 "group=Low"/"group=High"，键名 "Low"/"High" 匹配失败会静默
# 回退到灰色调色板。因此统一改用「无名向量 + 显式 legend.labs」。
KM_PALETTE <- c(unname(COL_RISK["Low"]), unname(COL_RISK["High"]))  # Low, High 顺序
KM_LEGEND  <- c("Low", "High")

# ---- 通用函数 ----

# 统一 AJCC 分期归类（根治各脚本各写一份、顺序不同导致的分组漂移）
# 关键：必须「长模式优先」——"Stage IIIA" 也匹配 I{1,2}[ABC]?$，若先判 I-II
# 会把全部 III 期误并入 I-II（曾导致 181 例错分）。
classify_stage <- function(x, unknown = "Unknown") {
  out <- dplyr::case_when(
    is.na(x)                              ~ unknown,
    grepl("IV",   x, fixed = TRUE)        ~ "IV",
    grepl("III",  x, fixed = TRUE)        ~ "III",
    grepl("^I{1,2}$|^Stage I{1,2}[ABC]?$", x) ~ "I-II",
    grepl("I{1,2}[ABC]$", x)              ~ "I-II",
    TRUE                                  ~ unknown)
  factor(out, levels = c("I-II", "III", "IV", "Unknown"))
}

# p 值格式化：避免 %.2e 在 p < 0.005e-300 时下溢成 "0.00e+00"
fmt_p <- function(p, digits = 2) {
  if (length(p) != 1L || is.na(p)) return("NA")
  if (p == 0 || p < 1e-300) return("< 1e-300")
  if (p < 1e-4) return(sprintf(paste0("%.", digits, "e"), p))
  sprintf(paste0("%.", max(digits, 3), "g"), p)
}
save_fig <- function(p, name, w = FIG_WIDTH, h = FIG_HEIGHT) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = FIG_DPI, bg = "white")
  try(ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p,
             width = w, height = h), silent = TRUE)
  message("✓ Saved figure: ", name)
}

save_table <- function(df, name) {
  write.csv(df, file.path(TAB_DIR, paste0(name, ".csv")), row.names = FALSE)
  message("✓ Saved table: ", name)
}

save_obj <- function(x, name) {
  saveRDS(x, file.path(RES_DIR, paste0(name, ".rds")))
  message("✓ Saved object: ", name)
}

# ---- 包加载（按需） ----
load_core_pkgs <- function() {
  pkgs <- c("dplyr", "tidyr", "ggplot2", "ggpubr", "patchwork",
            "survival", "survminer", "timeROC", "glmnet",
            "limma", "pheatmap", "ComplexHeatmap", "ggrepel",
            "RColorBrewer", "viridis", "ggsci")
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
}

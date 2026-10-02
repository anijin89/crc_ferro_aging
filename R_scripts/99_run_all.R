# ============================================================
# 99_run_all.R — 一键全流程（按实际存在的脚本，分阶段执行）
# ============================================================
# 用法：
#   Rscript R_scripts/99_run_all.R                # 只跑 S0（数据准备）
#   STAGES=all Rscript R_scripts/99_run_all.R     # 跑全部阶段
#   STAGES=S1,S2 Rscript R_scripts/99_run_all.R   # 只跑指定阶段
#   STOP_ON_ERROR=0 ...                           # 单脚本失败后继续
#
# 阶段划分（含脚本间依赖）：
#   S0 数据准备   01_data_download → 03_FARG_input
#   S1 差异+富集  04_DEG_intersection → 05_GO_KEGG              （依赖 S0）
#   S2 建模+验证  07_LASSO_multivariate → 08_ExternalValidation （依赖 S1）
#   S3 微环境临床 09_immune_infiltration → 10_Clinical_Table1 →
#                 11_Nomogram_DCA → 12_Druggability_DGIdb →
#                 13_hub_gene_TIMP1                             （依赖 S2）
#   S4 单细胞扰动 14_scRNA_celltype → 16_scTenifoldKnk           （依赖 S2）
#   S5 消融分型   17_UnionVsArms → 18_ExternalAblation →
#                 19_MSI_TMB_CMS                                （依赖 S2）
#   S6 空间组学   20_spatial_DSP                                （依赖 S2）
#   S7 图元       21_study_design（Fig 1 研究设计图，无上游依赖）
# ============================================================

# ---- 自动定位项目根目录（不写死本机路径）----
.script_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.script_file)) {
  .fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .script_file <- if (length(.fa)) normalizePath(sub("^--file=", "", .fa[1])) else file.path(getwd(), "x.R")
}
.ROOT <- dirname(dirname(.script_file))
while (!file.exists(file.path(.ROOT, "R_scripts", "00_config.R")) && !identical(.ROOT, dirname(.ROOT))) .ROOT <- dirname(.ROOT)

STAGES <- list(
  S0 = c("01_data_download.R", "03_FARG_input.R"),
  S1 = c("04_DEG_intersection.R", "05_GO_KEGG.R"),
  S2 = c("07_LASSO_multivariate.R", "08_ExternalValidation.R"),
  S3 = c("09_immune_infiltration.R", "10_Clinical_Table1.R", "11_Nomogram_DCA.R",
         "12_Druggability_DGIdb.R", "13_hub_gene_TIMP1.R"),
  S4 = c("14_scRNA_celltype.R", "16_scTenifoldKnk.R"),
  S5 = c("17_UnionVsArms.R", "18_ExternalAblation.R", "19_MSI_TMB_CMS.R"),
  S6 = c("20_spatial_DSP.R"),
  S7 = c("21_study_design.R")
)

# ---- 选择要跑的阶段 ----
sel <- Sys.getenv("STAGES", unset = "S0")
if (identical(tolower(sel), "all")) {
  to_run_stages <- names(STAGES)
} else {
  to_run_stages <- trimws(strsplit(sel, ",")[[1]])
}
bad <- setdiff(to_run_stages, names(STAGES))
if (length(bad)) stop("未知阶段: ", paste(bad, collapse = ", "),
                      "\n可用阶段: ", paste(names(STAGES), collapse = ", "), " 或 all")

to_run <- unlist(STAGES[to_run_stages], use.names = FALSE)

# ---- 存在性检查（先报错，不要跑到一半才发现）----
missing <- to_run[!file.exists(file.path(.ROOT, "R_scripts", to_run))]
if (length(missing)) stop("以下脚本不存在: ", paste(missing, collapse = ", "))

cat("================================================================\n")
cat("  RUN ALL — 项目根目录: ", .ROOT, "\n", sep = "")
cat("  阶段: ", paste(to_run_stages, collapse = ", "), "\n", sep = "")
cat("  脚本: ", paste(to_run, collapse = " -> "), "\n", sep = "")
cat("================================================================\n")

# ---- 顺序执行 ----
timing <- data.frame(script = character(), minutes = numeric(), status = character(),
                     stringsAsFactors = FALSE)
for (s in to_run) {
  message(sprintf("\n############ Running %s ############", s))
  t0 <- Sys.time()
  st <- tryCatch({
    source(file.path(.ROOT, "R_scripts", s))
    "OK"
  }, error = function(e) {
    message("x ", s, " 失败: ", conditionMessage(e))
    "FAIL"
  })
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  timing <- rbind(timing, data.frame(script = s, minutes = round(mins, 2), status = st))
  message(sprintf("%s %s done in %.1f min", if (st == "OK") "OK" else "FAIL", s, mins))
  if (st == "FAIL" && identical(Sys.getenv("STOP_ON_ERROR", "1"), "1")) {
    message("遇到错误，已中止（设 STOP_ON_ERROR=0 可继续）")
    break
  }
}

cat("\n=== 执行汇总 ===\n")
print(timing, row.names = FALSE)
if (dir.exists(file.path(.ROOT, "results"))) {
  saveRDS(timing, file.path(.ROOT, "results", "run_all_timing.rds"))
}
message("\n=== 完成 ===")

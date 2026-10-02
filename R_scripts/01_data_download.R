# ============================================================
# 01_data_download.R — TCGA + GEO 数据下载
# ============================================================
# 依赖：TCGAbiolinks, GEOquery, SummarizedExperiment
# 输出：data/raw/*.rds
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
  library(TCGAbiolinks)
  library(GEOquery)
  library(SummarizedExperiment)
})

# ------------------------------------------------------------
# 1. TCGA-COADREAD bulk RNA-seq + 临床
# ------------------------------------------------------------
download_tcga <- function() {
  message("=== 1. 下载 TCGA-COADREAD ===")

  # 1.1 RNA-seq (HTSeq-FPKM)
  if (!file.exists(file.path(DATA_RAW, "TCGA_COADREAD_fpkm.rds"))) {
    q <- GDCquery(project = "TCGA-COADREAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")
    GDCdownload(q, method = "api", files.per.chunk = 50)
    fpkm <- GDCprepare(q, save = TRUE,
                       save.filename = file.path(DATA_RAW, "TCGA_COADREAD_SE.rda"))
    saveRDS(fpkm, file.path(DATA_RAW, "TCGA_COADREAD_fpkm.rds"))
    message("✓ TCGA RNA-seq saved")
  }

  # 1.2 临床（TCGA-CDR via UCSC Xena or cBioPortal）
  if (!file.exists(file.path(DATA_RAW, "TCGA_COADREAD_clinical.csv"))) {
    clin <- GDCquery_clinic(project = "TCGA-COADREAD", type = "clinical")
    write.csv(clin, file.path(DATA_RAW, "TCGA_COADREAD_clinical.csv"), row.names = FALSE)
    message("✓ TCGA clinical saved")
  }

  # 1.3 突变 (MC3 MAF)
  if (!file.exists(file.path(DATA_RAW, "TCGA_COADREAD_maf.rds"))) {
    qmaf <- GDCquery(project = "TCGA-COADREAD",
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     access = "open")
    GDCdownload(qmaf)
    maf <- GDCprepare(qmaf)
    saveRDS(maf, file.path(DATA_RAW, "TCGA_COADREAD_maf.rds"))
    message("✓ TCGA MAF saved")
  }
}

# ------------------------------------------------------------
# 2. GEO bulk 验证集
# ------------------------------------------------------------
download_geo_bulk <- function(gses) {
  message("=== 2. 下载 GEO bulk 数据集 ===")
  for (gse in gses) {
    f <- file.path(DATA_RAW, paste0(gse, "_series_matrix.rds"))
    if (file.exists(f)) { message(sprintf("✓ %s already exists", gse)); next }
    tryCatch({
      g <- getGEO(GEO = gse, destdir = DATA_RAW, getGPL = TRUE)
      saveRDS(g, f)
      message(sprintf("✓ %s saved (n=%d)", gse, length(g[[1]]@dataDescription)))
    }, error = function(e) message(sprintf("✗ %s failed: %s", gse, e$message)))
  }
}

# ------------------------------------------------------------
# 3. GEO scRNA-seq (GSE132465)
# ------------------------------------------------------------
download_geo_scrna <- function(gse = GSE_SCRNA) {
  message("=== 3. 下载 scRNA-seq ===")
  # GSE132465 提供 supplementary files（10x 三个文件 per sample）
  # 用 GEOquery::getGEOSuppFiles 下载 supplementary
  out_dir <- file.path(DATA_RAW, gse)
  if (!dir.exists(out_dir) || length(list.files(out_dir)) < 5) {
    getGEOSuppFiles(GEO = gse, baseDir = DATA_RAW, fetchFiles = TRUE)
    message(sprintf("✓ %s supplementary downloaded", gse))
  } else {
    message(sprintf("✓ %s already exists", gse))
  }
}

# ------------------------------------------------------------
# 4. GEO 空间转录组（可选）
# ------------------------------------------------------------
download_geo_spatial <- function(gse = GSE_SPATIAL) {
  message("=== 4. 下载空间转录组（可选）===")
  out_dir <- file.path(DATA_RAW, gse)
  if (!dir.exists(out_dir)) {
    tryCatch({
      getGEOSuppFiles(GEO = gse, baseDir = DATA_RAW, fetchFiles = TRUE)
      message(sprintf("✓ %s downloaded", gse))
    }, error = function(e) message(sprintf("✗ %s failed: %s", gse, e$message)))
  }
}

# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------
download_tcga()
download_geo_bulk(c(GEO_VALIDATION, GEO_DEG_VALIDATION))
download_geo_scrna()
# download_geo_spatial()  # 取消注释启用

message("\n=== 全部下载完成，检查 data/raw/ ===")
list.files(DATA_RAW)

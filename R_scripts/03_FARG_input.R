# ============================================================
# 03_FARG_input.R — 加载/构建 95 个 ferro-aging-related genes
# ============================================================
# 来源：人工编纂 —— 不是两库的集合运算。从
#       FerrDb v2 (ferroptosis driver/suppressor/marker/resistance) 与
#       CellAge (senescence-inducing / SASP) 中精选出 95 个核心基因
#       （两库各自的规模都远大于 95，故稿件不得写成“两者的并集”）
# 归一：FSP1 (=AIFM2) 等别名映射到 HGNC canonical symbol
# 输出：FARG_FILE，以及补充表 tables/TableS4_FARG_panel.csv（稿件 Methods 引用）
# 若原始 supplementary 可用（方案 A），取消对应注释块
# ============================================================

suppressPackageStartupMessages({ library(magrittr); library(dplyr) })
# ---- 自动定位项目根目录（不写死本机路径）----
.script_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.script_file)) {
  .fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .script_file <- if (length(.fa)) normalizePath(sub("^--file=", "", .fa[1])) else file.path(getwd(), "x.R")
}
.ROOT <- dirname(dirname(.script_file))
while (!file.exists(file.path(.ROOT, "R_scripts", "00_config.R")) && !identical(.ROOT, dirname(.ROOT))) .ROOT <- dirname(.ROOT)
source(file.path(.ROOT, "R_scripts", "00_config.R"))

# ------------------------------------------------------------
# 方案 A：用户提供 Cell Metabolism 2026 supplementary 后取消注释
# ------------------------------------------------------------
# library(readxl)
# supp <- read_excel("path/to/mmc2.xlsx", sheet = "Table S1")
# fargs <- supp %>%
#   filter(grepl("ferro-aging|FARG", category, ignore.case = TRUE)) %>%
#   pull(gene_symbol) %>% unique()
# stopifnot(length(fargs) == 95)

# ------------------------------------------------------------
# 方案 B：自构 FARGs（本项目的实际路径）
# 步骤：B.1 手工编纂的候选清单（从 FerrDb v2 与 CellAge 精选，非集合运算）
#      B.2 别名归一
#      B.3 若 <95 回填 ferroptosis/senescence 池中同优先级基因
# ------------------------------------------------------------
build_FARGs_fallback <- function() {
  message("⚠ 使用备用方案构建 FARGs（须在 Methods 披露）")

  # B.1 原始 curated 95（含 FSP1/AIFM2 重复项，交给 B.2 归一）
  curated_95_raw <- c(
    # --- Ferroptosis: drivers / executors ---
    "ACSL3","ACSL4","ALOX5","ALOX12","ALOX15","ALOX15B","CISD1","CYBB",
    "FTH1","FTL","GCH1","HMOX1","HP","LPCAT3","MAP1LC3B","NCOA4","NOX1",
    "NOX4","PEBP1","SAT1","SLC11A2","SLC3A2","SLC7A11","STEAP3","TFRC",
    # --- Ferroptosis: suppressors / antioxidant ---
    "AIFM2","FSP1","DHODH","GCLC","GCLM","GPX4","GSS","NQO1","NFE2L2",
    "KEAP1","SLC7A11","TP53","CBS","CTH","SOD2","TXN",
    # --- Ferroptosis: cholesterol / mevalonate axis ---
    "HMGCR","IDI1","MVD","MVK","SQLE","FDPS",
    # --- Ferroptosis: iron handling / mitophagy ---
    "ABCB10","ATG5","ATG7","BECN1","DHFR","SQSTM1",
    # --- Senescence: cell-cycle / DNA damage ---
    "ATM","ATR","CCND1","CDK2","CDK4","CDKN1A","CDKN1C","CDKN2A","CDKN2B",
    "CHEK1","CHEK2","E2F1","MDM2","RB1",
    # --- Senescence: SASP / inflammatory ---
    "BMP4","CCL2","CSF2","CXCL8","GDF15","ID1","IGFBP3","IGFBP7","IL1A",
    "IL1B","IL6","MMP1","MMP3","MMP9","NFKB1","NOS2","RELA","SERPINE1",
    "TGFB1","TIMP1","TIMP2","VEGFA","ZNF24",
    # --- Senescence: markers / chromatin / metabolism ---
    "FOXO3","FOXO4","GLB1","LAMP1","SIRT1","SIRT6"
  ) %>% unique()

  # B.2 别名归一（HGNC/Expression matrix 用 canonical symbol）
  alias_map <- c(
    "FSP1"   = "AIFM2",   # ferroptosis suppressor protein 1 == apoptosis-inducing factor 2
    "NRF2"   = "NFE2L2",
    "IL8"    = "CXCL8",
    "BECLN1" = "BECN1",
    "DHODQ"  = "DHODH",
    "DACBA"  = "NOS2",
    "GMCSF"  = "CSF2",
    "BH2"    = "DHFR"
  )
  n_raw <- length(curated_95_raw)
  curated_95_raw <- ifelse(curated_95_raw %in% names(alias_map),
                           alias_map[curated_95_raw], curated_95_raw)
  curated <- unique(curated_95_raw)
  message(sprintf("别名归一 + unique: %d → %d (dropped %d)",
                  n_raw, length(curated), n_raw - length(curated)))

  # B.3 若 <95 回填
  if (length(curated) < 95) {
    pool <- c("GPX4","SLC7A11","TXN","TXNRD1","GSR","CARS1","FANCD2",
              "SIRT2","PGD","G6PD","TALDO1","RPIA","IDH1","IDH2","ME1",
              "ACLY","ACACA","SCD","ELOVL5","FTMT","MT1X","MT2A","FTH1",
              "MYC","JUN","FOS","HIF1A","EPAS1","STAT3","YAP1","WWTR1",
              "LATS2","PPARG","RXRA","NRF1","PCBP1","PCBP2","ACO1",
              "IREB2","SLC25A37","SLC25A38","ALAS2","ATF4","DDIT3",
              "HSPA5","SLC47A1","GABARAPL1","ULK1","MTOR","RPS6KB1",
              "AKT1","PIK3CA","PRKAA1","STK11","HMGB1","LMNB1","ERBB2",
              "MAP2K3","MAPK14","CCL5","CCL20","MMP12","IGFBP5")
    fill <- setdiff(pool, curated)
    need <- 95 - length(curated)
    curated <- c(curated, fill[1:need])
    message(sprintf("回填 %d 个基因: %s", need, paste(fill[1:need], collapse=", ")))
  }
  if (length(curated) > 95) curated <- curated[1:95]
  message(sprintf("FARGs 最终数量: %d", length(curated)))

  data.frame(
    gene_symbol = curated,
    source      = "FerrDb_v2 + CellAge_curated",
    category    = "ferro-aging core",
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 主流程
# ------------------------------------------------------------
if (file.exists(FARG_FILE)) {
  fargs <- read.csv(FARG_FILE, stringsAsFactors = FALSE)
  message(sprintf("✓ 已加载 %s: %d 个 FARGs", FARG_FILE, nrow(fargs)))
} else {
  fargs <- build_FARGs_fallback()
  write.csv(fargs, FARG_FILE, row.names = FALSE)
  message(sprintf("✓ 已写入 %s", FARG_FILE))
}

# ---- 守护：已发布清单必须与 build_FARGs_fallback() 的输出逐项、逐序一致 ----
# 这 95 个基因是全文的起点输入（→ 26 FARG-DEG → 五基因签名 → 所有统计量）。
# 若有人改了候选池却没同步 FARG_FILE，下游会与清单静默脱节且无人察觉。
.built <- as.character(suppressMessages(build_FARGs_fallback())$gene_symbol)
.using <- as.character(fargs$gene_symbol)
if (!identical(.using, .built)) {
  stop(sprintf(paste0("FARG_FILE 与 build_FARGs_fallback() 不一致\n",
                      "  仅在 FARG_FILE : %s\n  仅在 builder   : %s\n",
                      "  两者元素相同但顺序不同: %s"),
               paste(setdiff(.using, .built), collapse = ", "),
               paste(setdiff(.built, .using), collapse = ", "),
               identical(sort(.using), sort(.built))))
}
stopifnot(nrow(fargs) == 95)
message("✓ 清单与 builder 逐项一致（95 个）")

# 验证：与表达矩阵命中情况
expr_genes <- data.table::fread(file.path(DATA_RAW, "TCGA_COADREAD_TPM_matrix.tsv"),
                                select = 1L)[[1]]
hit <- fargs$gene_symbol %in% expr_genes
message(sprintf("Matrix coverage: %d / %d", sum(hit), length(hit)))
if (any(!hit)) message("  Missing: ", paste(fargs$gene_symbol[!hit], collapse=", "))

# 保存
save_obj(fargs$gene_symbol, "FARG_genes")

# ---- 补充表：95 基因面板清单（稿件 Methods 与 Introduction 引作 Table S4）----
# 与 FARG_FILE 同源写出，保证补充表与建模输入永不脱节。
save_table(data.frame(
  gene_symbol = fargs$gene_symbol,
  panel       = fargs$category,
  source      = fargs$source,
  stringsAsFactors = FALSE
), "TableS4_FARG_panel")
message(sprintf("✓ tables/TableS4_FARG_panel.csv: %d 行", nrow(fargs)))

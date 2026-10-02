# ============================================================
# 21_study_design.R — Fig 1：研究设计与样本筛选示意图
# ============================================================
#   这是一个「图元」脚本：不做任何统计，只把论文的数据结构与样本
#   筛选链条画成一张矢量图（PDF + 300 dpi PNG）。
#
#   Panel A — 四个数据集与研究流程的映射：
#       TCGA-COADREAD（bulk RNA-seq，发现/训练）
#       GSE39582（microarray，外部验证）
#       GSE132465（scRNA-seq，hub 基因的细胞定位）
#       GSE281413（GeoMx DSP，原位空间共定位）
#   Panel B — TCGA-COADREAD 的样本筛选漏斗。
#
#   依赖：无上游脚本（只读 data/raw 的样本/临床注释 + 若干已产出的
#   results / tables）。所有画在图上的数字都由脚本现算，并在
#   stopifnot() 中断言；原始数据若变化会直接报错而不是静默画错。
#
#   设备说明：本机 R 构建的 cairo 不可用（缺 X11 依赖），PNG 用
#   quartz 设备；图内文字全部使用 ASCII + Latin-1 字符，使 pdf()
#   与 quartz 两个设备输出一致，不触发字符替换。
#
#   输出：figures/Fig1_study_design.png / .pdf
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

suppressPackageStartupMessages(library(grid))

# ---- 引出线走线样式 --------------------------------------------------------
#   diag     : 第 4 框中心斜向连到总线中点（**当前正式稿**，2026-09-20 作者选定）
#   elbow    : 向下 -> 弯折向左 -> 收到中线 -> 接入总线（弯折指向中间，备选）
#   straight : 第 4 框中心直接垂直落到总线右端（改版前的旧样式，备选）
# 为什么不用 straight：第 4 框中心 x=0.8445 与第 3 列色条中心 x=0.840 几乎重合，
#   垂直落下的线看着只接了第 3 列，旁边那条总线倒像另一根独立的线。
# FIG1_HEAD：在引出线的**末端**补一个向下的箭头，落在总线上。
#   为什么需要：走线收到中线后向下接入总线那一段，与中列自己的箭头
#   （总线 -> 中列色条）完全共线，两条线在总线处连成一条，看不出
#   「引出线到此收尾、总线再把流量分给三列」。补一个箭头落在总线上，
#   这个收尾就有了明确的落点。diag 与 elbow 默认补，straight 强制不补。
# 出预览/复现旧稿：FIG1_FEED=<样式> FIG1_HEAD=<0|1> FIG1_TAG=<后缀>
#   （带 FIG1_TAG 时只写 _preview/<名>_<后缀>.png，不碰 figures/ 的正式产物）
FEED_STYLE  <- Sys.getenv("FIG1_FEED", "elbow")
FEED_HEAD   <- identical(Sys.getenv("FIG1_HEAD", "1"), "1")
PREVIEW_TAG <- Sys.getenv("FIG1_TAG", "")
stopifnot(FEED_STYLE %in% c("straight", "elbow", "diag"))
# 直落式本身就落在总线右端，没有「弯折到中间」的收尾可言，故强制不补箭头 ——
# 这样单给一个 FIG1_FEED=straight 就能复现改版前的旧图，不会被默认的 HEAD 拖住
if (identical(FEED_STYLE, "straight")) FEED_HEAD <- FALSE

# ============================================================
# 1. 现算所有要标注的数字
# ============================================================
message("--- 核算样本筛选漏斗 ---")
stypes <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_sample_types.csv"), stringsAsFactors = FALSE)
clin   <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_clinical.csv"),   stringsAsFactors = FALSE)

# 表达矩阵中的样本列（只读首行，不加载 223 MB 矩阵）
.con <- file(file.path(DATA_RAW, "TCGA_COADREAD_TPM_matrix.tsv"), "r")
.hdr <- readLines(.con, n = 1); close(.con)
mat_samples <- strsplit(.hdr, "\t")[[1]]
mat_samples <- mat_samples[mat_samples %in% stypes$sample_id]

st <- stypes[stypes$sample_id %in% mat_samples, ]
tt <- table(st$sample_type)
n_total   <- nrow(st)
n_primary <- as.integer(tt[["Primary Tumor"]])
n_normal  <- as.integer(tt[["Solid Tissue Normal"]])
n_rec     <- as.integer(tt[["Recurrent Tumor"]])
n_met     <- as.integer(tt[["Metastatic"]])

pt        <- st[st$sample_type == "Primary Tumor", ]
n_pat_all <- length(unique(pt$patient_barcode))
dup_tab   <- table(pt$patient_barcode)
n_dup_pat <- sum(dup_tab > 1)
n_surplus <- nrow(pt) - n_pat_all

# 与 07_LASSO_multivariate.R 相同的去重规则：每位患者保留 sample_id 最小者
pt_d <- pt[order(pt$sample_id), ]
pt_d <- pt_d[!duplicated(pt_d$patient_barcode), ]
clin_sub  <- clin[match(pt_d$patient_barcode, clin$patient_barcode), ]
os_time   <- clin_sub$os_days
os_status <- ifelse(clin_sub$vital_status == "Dead", 1,
             ifelse(clin_sub$vital_status == "Alive", 0, NA))
valid     <- !is.na(os_time) & os_time > 0 & !is.na(os_status)
n_drop_os <- sum(!valid)
n_final   <- sum(valid)
n_events  <- sum(os_status[valid] == 1, na.rm = TRUE)

cat(sprintf("  %d 样本 -> %d 原发（%d 患者）-> 去重后 %d -> 有 OS %d（%d 例死亡）\n",
            n_total, n_primary, n_pat_all, nrow(pt_d), n_final, n_events))
cat(sprintf("  扣除：非原发 %d；多余原发 %d（%d 位患者）；无可用 OS %d\n",
            n_total - n_primary, n_surplus, n_dup_pat, n_drop_os))

message("--- 读取其他产物中的数字 ---")
n_farg     <- length(readRDS(file.path(RES_DIR, "FARG_genes.rds")))              # 95
farg_degs  <- readRDS(file.path(RES_DIR, "FARG_DEG_intersection.rds"))           # 26
n_farg_deg <- length(farg_degs)

deg_tab <- read.csv(file.path(TAB_DIR, "TableS1_DEG_full.csv"), stringsAsFactors = FALSE)
n_up    <- sum(deg_tab$sig == "Up")
n_down  <- sum(deg_tab$sig == "Down")

sig_genes <- read.csv(file.path(TAB_DIR, "Table3_LASSO_selected_beta.csv"), stringsAsFactors = FALSE)$gene
sig_lab   <- paste(sig_genes, collapse = " \u00b7 ")

ext  <- read.csv(file.path(TAB_DIR, "Table4_external_validation_summary.csv"), stringsAsFactors = FALSE)
gse  <- ext[grepl("GSE39582", ext$Cohort, fixed = TRUE), ]
n_gse  <- gse$N[1]
ev_gse <- gse$Events[1]
hr_gse <- gse$Cox_HR_per_SD[1]
p_gse  <- gse$Cox_p[1]

spat   <- readRDS(file.path(RES_DIR, "spatial_DSP.rds"))$spot
n_roi  <- nrow(spat)
.base  <- sub("-[0-9]+$", "", trimws(spat$treatment))          # 去掉载玻片序号后缀
.donor <- sub("^[A-Z]+\\s*[A-Z]{2}\\s*", "", .base)            # 去掉 EO/AO + CA/NL 前缀
n_donor <- length(unique(.donor[spat$cancer_flag == "Tumour"]))

N_CELLS <- 63689L   # GSE132465 单细胞数；由 14_scRNA_celltype.R 运行日志得到，
                    # results/ 未保存汇总对象。00_verify_data.R 第 16 节会与正文对账。

# ---- 断言：图上每个数字都必须来自数据，不能被默默改掉 ----
stopifnot(
  n_total == 701L, n_primary == 647L, n_normal == 51L, n_rec == 2L, n_met == 1L,
  n_pat_all == 624L, n_dup_pat == 13L, n_surplus == 23L,
  n_drop_os == 34L, n_final == 590L, n_events == 121L,
  n_farg == 95L, n_farg_deg == 26L, n_up == 1715L, n_down == 1641L,
  length(sig_genes) == 5L,
  n_roi == 111L, n_donor == 6L
)
message("OK 所有标注数字与项目产物一致")

# ============================================================
# 2. 绘图工具箱
# ============================================================
COL <- list(
  disc  = list(line = "#3D5A80", fill = "#E4EBF3", head = "#3D5A80"),
  val   = list(line = "#1B998B", fill = "#E2F1EE", head = "#1B998B"),
  mech  = list(line = "#7B5EA7", fill = "#EDE9F5", head = "#7B5EA7"),
  clin  = list(line = "#C08552", fill = "#FAF1E6", head = "#C08552"),
  fun   = list(line = "#565F6E", fill = "#EEF1F5", head = "#565F6E"),
  hit   = list(line = "#1B998B", fill = "#D8EDE8", head = "#1B998B"),
  arrow = "#5A6472", text = "#1A1A1A", sub = "#5A6472", minus = "#B03A48",
  rule  = "#C9D1DA"
)

box_grob <- function(x0, x1, y0, y1, label, col, fontsize = 6.8, lwd = 0.9,
                     r = 0.02, textcol = COL$text, bold = FALSE) {
  grid.roundrect(x = unit((x0 + x1) / 2, "npc"), y = unit((y0 + y1) / 2, "npc"),
                 width = unit(x1 - x0, "npc"), height = unit(y1 - y0, "npc"),
                 r = unit(r, "npc"), gp = gpar(fill = col$fill, col = col$line, lwd = lwd))
  grid.text(label, x = unit((x0 + x1) / 2, "npc"), y = unit((y0 + y1) / 2, "npc"),
            gp = gpar(fontsize = fontsize, col = textcol, lineheight = 1.22,
                      fontface = if (bold) 2 else 1))
}
h_arrow <- function(x0, x1, y, lwd = 1.0) {
  grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y, y), "npc"),
             arrow = arrow(length = unit(1.7, "mm"), type = "closed"),
             gp = gpar(lwd = lwd, col = COL$arrow, lineend = "butt"))
}
v_arrow <- function(x, y0, y1, lwd = 1.0) {
  grid.lines(x = unit(c(x, x), "npc"), y = unit(c(y0, y1), "npc"),
             arrow = arrow(length = unit(1.7, "mm"), type = "closed"),
             gp = gpar(lwd = lwd, col = COL$arrow, lineend = "butt"))
}
hline <- function(x0, x1, y, col = COL$rule, lwd = 0.5, lty = 1) {
  grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y, y), "npc"),
             gp = gpar(lwd = lwd, col = col, lty = lty))
}

# Panel A 视口的物理尺寸（draw_all() 里的 0.91 x 0.510），用于把「英寸」换算成
# npc。面板宽高比约 2.3:1，同一个物理半径在两个方向的 npc 增量并不相等，
# 画圆角弯折必须分开换算，否则圆角会被拉成椭圆。
PANEL_A_W_IN <- 0.91 * 7.2
PANEL_A_H_IN <- 0.51 * 5.6

# 第 4 框的引出线：向下 -> 弯折向左 -> 收到中线 -> 向下接入总线。
# 两个拐点各用一个四分之一椭圆弧代替直角。注意 npc 的 y 越小越靠下，圆角弧
# 必须切在「拐点之前」的两条边上（竖直边在拐点上方、水平边在拐点左侧），
# 圆心因此在转弯的外侧，否则弧会绕到拐点后面去。半径按英寸给定、两个方向
# 分别换算成 npc，保证圆角在任何输出尺寸下都是正圆角而不是椭圆。
elbow_path <- function(x4, x_mid, y0, y_feed, y_bus, r_in = 0.045, n = 12,
                       head = FALSE) {
  rx <- r_in / PANEL_A_W_IN; ry <- r_in / PANEL_A_H_IN
  k1 <- min(1, (y0 - y_feed) / (2 * ry), (x4 - x_mid) / (2 * rx))
  k2 <- min(1, (y_feed - y_bus) / (2 * ry), (x4 - x_mid) / (2 * rx))
  e1x <- rx * k1; e1y <- ry * k1
  e2x <- rx * k2; e2y <- ry * k2
  c1x <- x4 - e1x;    c1y <- y_feed + e1y     # 右上角：竖直段在拐点上方切出
  c2x <- x_mid + e2x; c2y <- y_feed - e2y     # 左下角：竖直段在拐点下方切出
  a1 <- seq(0, -pi / 2, length.out = n)       # 从竖直（角度 0）转到水平（-pi/2）
  a2 <- seq(pi / 2, pi, length.out = n)       # 从水平（pi/2）转到竖直（pi）
  xs <- c(x4, x4, c1x + e1x * cos(a1), x4 - e1x, x_mid + e2x,
          c2x + e2x * cos(a2), x_mid)
  ys <- c(y0, y_feed + e1y, c1y + e1y * sin(a1), y_feed, y_feed,
          c2y + e2y * sin(a2), y_feed - e2y)
  # head = TRUE 时最后一段竖直留给调用方用带箭头的线画，这里不重复描一遍，
  # 否则圆角末端会被无箭头的直线盖住，接缝处会看到一小段毛刺。
  if (!head) {
    xs <- c(xs, x_mid); ys <- c(ys, y_bus)
  }
  stopifnot(
            length(xs) == length(ys),
            all(diff(ys) <= 1e-9),                # 走线自始至终只向下 / 向左，无回头
            all(diff(xs) <= 1e-9))
  list(x = xs, y = ys, y_tail = y_feed - e2y)
}

# ============================================================
# 3. Panel A — 研究设计与分析流程
# ============================================================
panel_A <- function() {
  grid.text("A", x = unit(0.001, "npc"), y = unit(0.997, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 10.5, fontface = 2, col = COL$text))

  # ---- 版式：两种走线共用同一套面板坐标，只有第 4 框向下的引出线画法不同 ----
  # 走线要在三列色条上方占用竖直净空：elbow 要「两段竖直 + 一条水平」，
  # 补箭头时末端那根箭头杆也要挤得出来，所以把三列色条与下面两排子框整体
  # 下移，把这段净空让给走线。四种配置里向下的三根箭头长度都保持一致
  # （总线到色条上沿一律 0.076）。
  if (FEED_HEAD) {
    # 末端箭头头部 1.7 mm，杆长还要留够（否则画出来只有一个三角形）。
    # 第 4 框下沿是 0.715，所以把水平段抬到 0.700、总线压到 0.628，
    # 中间 0.072（约 5.2 mm）扣掉 1.14 mm 的圆角，杆长还有约 2.4 mm。
    y_feed <- 0.700; y_bus <- 0.628
    hd <- c(0.492, 0.552)
    a0 <- 0.294; a1 <- 0.464
    b0 <- 0.094; b1 <- 0.264
  } else if (identical(FEED_STYLE, "elbow")) {
    y_feed <- 0.684                      # 弯折的水平段
    y_bus  <- 0.646                      # 总线
    hd     <- c(0.510, 0.570)            # 三列色条
    a0 <- 0.312; a1 <- 0.482             # 上层子框
    b0 <- 0.112; b1 <- 0.282             # 下层子框
  } else {
    y_feed <- NA_real_
    y_bus  <- 0.672
    hd     <- c(0.535, 0.595)
    a0 <- 0.335; a1 <- 0.505
    b0 <- 0.135; b1 <- 0.305
  }

  # ---- 第一行：发现与建模（TCGA）----
  # 框间留 0.058 的净空（原先 0.037），箭头杆长才够 —— 杆长若与箭头头部
  # （1.7 mm）相当，画出来只剩一个三角形，看着像没画箭头。
  y0 <- 0.715; y1 <- 0.935
  yc <- (y0 + y1) / 2
  r1 <- list(c(0.005, 0.200), c(0.258, 0.398), c(0.456, 0.636), c(0.694, 0.995))
  box_grob(r1[[1]][1], r1[[1]][2], y0, y1,
           sprintf("TCGA-COADREAD\nbulk RNA-seq\n%d primary tumours\n%d solid tissue normals", n_primary, n_normal),
           COL$disc, 6.7)
  box_grob(r1[[2]][1], r1[[2]][2], y0, y1,
           sprintf("Differential\nexpression\n%d DEGs\n%d up /\n%d down", n_up + n_down, n_up, n_down),
           COL$fun, 6.4)
  box_grob(r1[[3]][1], r1[[3]][2], y0, y1,
           sprintf("Intersection with the\n%d-gene ferro-aging\npanel -> %d FARG-DEGs", n_farg, n_farg_deg),
           COL$fun, 6.4)
  box_grob(r1[[4]][1], r1[[4]][2], y0, y1,
           sprintf("Univariate Cox screening + LASSO\n-> five-gene prognostic score\n%s", sig_lab),
           COL$disc, 6.8)
  h_arrow(r1[[1]][2] + 0.008, r1[[2]][1] - 0.008, yc)
  h_arrow(r1[[2]][2] + 0.008, r1[[3]][1] - 0.008, yc)
  h_arrow(r1[[3]][2] + 0.008, r1[[4]][1] - 0.008, yc)

  # ---- 分流总线：score 框下沿 -> 三列 ----
  # v_arrow 的杆长 = y_bus - 表头下沿，必须明显大于箭头头部，否则三根向下的
  # 箭头看起来只是贴在色条上的三角形。表头色条下移以让出这段。
  cols <- list(c(0.005, 0.310), c(0.345, 0.650), c(0.685, 0.995))
  .xs  <- vapply(cols, function(z) mean(z), numeric(1))
  .x_from <- mean(r1[[4]])
  if (identical(FEED_STYLE, "elbow")) {
    # 向下 -> 弯折向左 -> 收到中线 -> 向下，等于把三列收敛的那一点挪到正中，
    # 三根向下箭头由中间一点对称发散，而不是全部挂在总线的右端。
    pth <- elbow_path(.x_from, .xs[2], y0, y_feed, y_bus, head = FEED_HEAD)
    grid.lines(x = unit(pth$x, "npc"), y = unit(pth$y, "npc"),
               gp = gpar(lwd = 1.0, col = COL$arrow, lineend = "butt",
                         linejoin = "round"))
    if (FEED_HEAD) v_arrow(.xs[2], pth$y_tail, y_bus)
  } else if (identical(FEED_STYLE, "diag")) {
    # 斜段本身没有竖直段可收尾：补箭头时斜段收到最高处（y_land = y_feed），
    # 末端再补一小段竖直箭头落到总线上，杆长才能明显大于箭头头部。
    y_land <- if (FEED_HEAD) y_feed else y_bus
    grid.lines(x = unit(c(.x_from, .xs[2]), "npc"), y = unit(c(y0, y_land), "npc"),
               gp = gpar(lwd = 1.0, col = COL$arrow))
    if (FEED_HEAD) v_arrow(.xs[2], y_land, y_bus)
  } else {
    grid.lines(x = unit(c(.x_from, .x_from), "npc"), y = unit(c(y0, y_bus), "npc"),
               gp = gpar(lwd = 1.0, col = COL$arrow))
  }
  hline(min(.xs), max(.xs), y_bus, col = COL$arrow, lwd = 1.0)
  for (i in 1:3) v_arrow(.xs[i], y_bus, hd[2] + 0.003)

  # ---- 三个方向 ----
  heads <- c("External validation", "Mechanistic dissection of TIMP1",
             "Clinical translation and robustness")
  sch   <- list(COL$val, COL$mech, COL$clin)
  for (i in 1:3) {
    box_grob(cols[[i]][1], cols[[i]][2], hd[1], hd[2], heads[i],
             list(line = sch[[i]]$head, fill = sch[[i]]$head),
             6.9, textcol = "white", bold = TRUE, r = 0.012)
  }

  # 子框的纵向位置在上面的版式块里给定（a0/a1/b0/b1）

  # 外部队列验证
  box_grob(cols[[1]][1], cols[[1]][2], a0, a1,
           sprintf("GSE39582 (microarray)\nn = %d \u00b7 %d deaths", n_gse, ev_gse), COL$val, 6.7)
  box_grob(cols[[1]][1], cols[[1]][2], b0, b1,
           sprintf("Kaplan-Meier \u00b7 time-dependent ROC\nper-SD HR %.2f (p = %.3f)", hr_gse, p_gse), COL$val, 6.7)

  # 机制层（TIMP1）
  box_grob(cols[[2]][1], cols[[2]][2], a0, a1,
           sprintf("GSE132465 (scRNA-seq)\n%s cells \u00b7 6 major cell types",
                   format(N_CELLS, big.mark = ",")), COL$mech, 6.7)
  box_grob(cols[[2]][1], cols[[2]][2], b0, b1,
           sprintf("GSE281413 (GeoMx DSP)\n%d ROIs \u00b7 %d patients\nspatial TIMP1-CAF co-localisation", n_roi, n_donor), COL$mech, 6.7)

  # 临床转化与稳健性
  box_grob(cols[[3]][1], cols[[3]][2], a0, a1,
           "Nomogram \u00b7 calibration \u00b7 DCA\nDGIdb druggability annotation", COL$clin, 6.7)
  box_grob(cols[[3]][1], cols[[3]][2], b0, b1,
           "Compositional ablation\n(union vs single-programme arms)\nIndependence of MSI, TMB, subtype", COL$clin, 6.7)

  grid.text("All panels are reported in Fig 2-14 (Results sections 1-9).",
            x = unit(0.5, "npc"), y = unit(0.045, "npc"),
            gp = gpar(fontsize = 6.1, col = COL$sub, fontface = 3))
}

# ============================================================
# 4. Panel B — TCGA-COADREAD 样本筛选漏斗
# ============================================================
panel_B <- function() {
  grid.text("B", x = unit(0.001, "npc"), y = unit(0.997, "npc"), just = c("left", "top"),
            gp = gpar(fontsize = 10.5, fontface = 2, col = COL$text))

  y0 <- 0.580; y1 <- 0.920; yc <- (y0 + y1) / 2
  # 框间净空 0.082（原先 0.070）、徽标收窄到 0.030：箭头杆在徽标两侧各留出
  # 0.020 的可见长度，否则徽标几乎把整根箭头吃掉。
  fx <- list(c(0.005, 0.190), c(0.272, 0.417), c(0.499, 0.644), c(0.726, 0.995))

  box_grob(fx[[1]][1], fx[[1]][2], y0, y1,
           sprintf("TCGA-COADREAD\n%d samples profiled\n%d primary tumours\n%d solid tissue normals\n%d recurrent, %d metastatic",
                   n_total, n_primary, n_normal, n_rec, n_met), COL$fun, 6.3)
  box_grob(fx[[2]][1], fx[[2]][2], y0, y1,
           sprintf("%d primary\ntumours from\n%d patients", n_primary, n_pat_all), COL$fun, 6.6)
  box_grob(fx[[3]][1], fx[[3]][2], y0, y1,
           sprintf("%d patients\nwith expression\ndata", nrow(pt_d)), COL$fun, 6.6)
  box_grob(fx[[4]][1], fx[[4]][2], y0, y1,
           sprintf("%d patients \u00b7 %d deaths\nsurvival modelling cohort", n_final, n_events),
           COL$hit, 6.9)

  # 箭头 + 扣除量徽标
  rem <- c(n_total - n_primary, n_surplus, n_drop_os)
  for (i in 1:3) {
    gx0 <- fx[[i]][2] + 0.006; gx1 <- fx[[i + 1]][1] - 0.006
    h_arrow(gx0, gx1, yc)
    gcx <- (gx0 + gx1) / 2
    grid.roundrect(x = unit(gcx, "npc"), y = unit(yc, "npc"),
                   width = unit(0.030, "npc"), height = unit(0.105, "npc"),
                   r = unit(0.012, "npc"),
                   gp = gpar(fill = "white", col = COL$minus, lwd = 0.8))
    grid.text(sprintf("-%d", rem[i]), x = unit(gcx, "npc"), y = unit(yc, "npc"),
              gp = gpar(fontsize = 6.2, col = COL$minus, fontface = 2))
  }

  # 扣除内容说明
  hline(0.005, 0.995, 0.520)
  desc <- c(sprintf("non-primary samples (%d solid tissue normals, %d recurrent, %d metastatic tumour)",
                    n_normal, n_rec, n_met),
            sprintf("surplus primaries from %d patients contributing more than one primary tumour (%d samples -> %d)",
                    n_dup_pat, n_dup_pat + n_surplus, n_dup_pat),
            sprintf("patients without a usable overall-survival time"))
  yy <- c(0.400, 0.280, 0.160)
  for (i in 1:3) {
    grid.text(sprintf("-%d", rem[i]), x = unit(0.052, "npc"), y = unit(yy[i], "npc"),
              just = "right", gp = gpar(fontsize = 6.4, col = COL$minus, fontface = 2))
    grid.text(desc[i], x = unit(0.062, "npc"), y = unit(yy[i], "npc"),
              just = "left", gp = gpar(fontsize = 6.4, col = COL$sub))
  }
  grid.text(sprintf("Differential expression (Fig 2A-B) used all %d primary tumours versus the %d solid tissue normals.",
                    n_primary, n_normal),
            x = unit(0.005, "npc"), y = unit(0.030, "npc"), just = "left",
            gp = gpar(fontsize = 6.1, col = COL$sub, fontface = 3))
  grid.text("Retained for downstream analysis",
            x = unit(0.995, "npc"), y = unit(0.030, "npc"), just = "right",
            gp = gpar(fontsize = 6.1, col = COL$hit$line, fontface = 3))
}

# ============================================================
# 5. 输出
# ============================================================
W <- 7.2; H <- 5.6

draw_all <- function() {
  grid.newpage()
  grid.text("Study design and sample selection",
            x = unit(0.5, "npc"), y = unit(0.985, "npc"), just = c("center", "top"),
            gp = gpar(fontsize = 10.5, fontface = 2, col = COL$text))
  pushViewport(viewport(x = 0.045, y = 0.415, width = 0.91, height = 0.510, just = c("left", "bottom")))
  panel_A(); popViewport()
  pushViewport(viewport(x = 0.045, y = 0.040, width = 0.91, height = 0.320, just = c("left", "bottom")))
  panel_B(); popViewport()
}

.fig <- file.path(FIG_DIR, "Fig1_study_design")
.ex  <- c("png", "pdf")
if (nzchar(PREVIEW_TAG)) {                    # 预览模式：只出 PNG，不碰正式产物
  .dir <- file.path(.ROOT, "_preview")
  dir.create(.dir, showWarnings = FALSE)
  .fig <- file.path(.dir, sprintf("Fig1_study_design_%s", PREVIEW_TAG))
  .ex  <- "png"
}

if ("png" %in% .ex) {
  png(paste0(.fig, ".png"), width = W, height = H, units = "in",
      res = FIG_DPI, type = "quartz", bg = "white")
  draw_all(); dev.off()
}
if ("pdf" %in% .ex) {
  pdf(paste0(.fig, ".pdf"), width = W, height = H, bg = "white")
  draw_all(); dev.off()
}
message("OK Saved figure: ", basename(.fig), "  [feed = ", FEED_STYLE, "]")

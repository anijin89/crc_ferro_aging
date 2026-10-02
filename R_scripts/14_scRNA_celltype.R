# ============================================================
# 14_scRNA_celltype.R — GSE132465 单细胞验证 TIMP1+CAF 轴
# ============================================================
# 数据：Lee et al. 2020 (Gustafson) CRC scRNA-seq, GEO processed
#   - natural_log_TPM_matrix (ln(TPM+1)), 30 genes × 63,689 cells (已抽取)
#   - cell_annotation: Cell_type(6) / Cell_subtype / Patient / Sample
# 目标：
#   1) TIMP1 是否富集于 Myofibroblast(CAF)
#   2) 5 个模型基因 (TIMP1/MMP1/MMP3/CDKN2A/CHEK1) 的细胞定位
#   3) CAF 标志与 TIMP1 的单细胞共表达
#   4) ferroptosis/senescence 基因在 CAF 中的程序活性
# 输出：Fig7A-D, Table7a-d
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

MODEL_GENES <- c("TIMP1","MMP1","MMP3","CDKN2A","CHEK1")
CAF_MARKERS <- c("ACTA2","FAP","COL1A1","COL1A2","POSTN","SPARC","LUM","TAGLN","PDGFRB","FN1")
FERRO_GENES <- c("ACSL4","GPX4","SLC7A11","TFRC","FTH1","HMOX1","SLC11A2","AIFM2")
SENS_GENES  <- c("CDKN1A","CDKN2A","IL6","CXCL8","SERPINE1","GDF15","TP53","LMNB1")

# ---- 1. 读取表达子集 (genes × cells) ----
sub <- data.table::fread(file.path(DATA_PROC, "GSE132465_genes_subset.tsv"))
gene_ids <- sub[[1]]
expr <- as.matrix(sub[, -1, with = FALSE])
rownames(expr) <- gene_ids; rm(sub)
# ln(TPM+1) → 线性 TPM（用于 %positive 与均值更直观）
expr_lin <- expr2lin <- expr
expr_lin[] <- (exp(expr) - 1)          # 还原 TPM
storage.mode(expr_lin) <- "double"

# ---- 2. 读取注释并对齐 ----
ann <- data.table::fread(file.path(DATA_RAW, "GSE132465",
             "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"))
setnames(ann, tolower(names(ann)))
ann <- as.data.frame(ann)
rownames(ann) <- ann$index
ann <- ann[colnames(expr_lin), ]       # 与表达列严格对齐
stopifnot(all(ann$index == colnames(expr_lin)))
message(sprintf("✓ aligned cells: %d", ncol(expr_lin)))

# 细胞类型顺序（按丰度）
ct_levels <- names(sort(table(ann$cell_type), decreasing = TRUE))
ann$cell_type <- factor(ann$cell_type, levels = ct_levels)

# ---- 3. 每细胞类型平均表达 & 阳性率 ----
agg_by_ct <- function(genes, group_vec) {
  grps <- levels(factor(group_vec))
  mean_mat <- pct_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(grps),
                                dimnames = list(genes, grps))
  for (g in grps) {
    idx <- which(group_vec == g)
    m <- expr_lin[genes, idx, drop = FALSE]
    mean_mat[, g] <- rowMeans(m)
    pct_mat[, g]  <- rowSums(m > 0.1) / length(idx) * 100
  }
  list(mean = mean_mat, pct = pct_mat)
}

all_genes <- unique(c(MODEL_GENES, CAF_MARKERS, FERRO_GENES, SENS_GENES))
all_genes <- intersect(all_genes, rownames(expr_lin))
ct_res <- agg_by_ct(all_genes, ann$cell_type)
save_table(as.data.frame(ct_res$mean) %>% mutate(gene = rownames(ct_res$mean)),
           "Table7a_scRNA_mean_by_celltype")
save_table(as.data.frame(ct_res$pct) %>% mutate(gene = rownames(ct_res$pct)),
           "Table7b_scRNA_pctpos_by_celltype")

# ---- 4. Fig 7A: 5 模型基因 + CAF 标志 热图（按 Cell_type 平均 z-score）----
heat_genes <- intersect(c(MODEL_GENES, CAF_MARKERS), rownames(ct_res$mean))
hm <- ct_res$mean[heat_genes, , drop = FALSE]
hm_z <- t(scale(t(log2(hm + 1))))
hm_z[hm_z > 2] <- 2; hm_z[hm_z < -2] <- -2
png(file.path(FIG_DIR, "Fig7A_scRNA_modelCAF_heatmap.png"),
    width = 7, height = 7, units = "in", res = FIG_DPI)
pheatmap::pheatmap(hm_z,
                   cluster_rows = TRUE, cluster_cols = FALSE,
                   color = colorRampPalette(c("#1B998B","white","#D7263D"))(100),
                   display_numbers = TRUE, number_format = "%.1f", fontsize_number = 6,
                   main = "Model + CAF genes: mean expr by cell type (z-score)")
dev.off()

# ---- 5. Fig 7B: TIMP1 在各细胞类型 violin（%positive 标注）----
tim_df <- data.frame(
  TPM = as.numeric(expr_lin["TIMP1", ]),
  cell_type = ann$cell_type
)
p_tim <- ggplot(tim_df, aes(cell_type, log2(TPM + 1), fill = cell_type)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  scale_fill_brewer(palette = "Set2") +
  stat_summary(aes(y = log2(TPM + 1)), fun = mean, geom = "point",
               shape = 23, fill = "black", size = 2) +
  labs(x = NULL, y = expression("log2(TPM+1)"),
       title = "TIMP1 expression across cell types (GSE132465)") +
  theme_classic(base_size = 11) + theme(legend.position = "none",
                                        axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_tim, "Fig7B_TIMP1_violin_celltype", w = 7, h = 5)

# TIMP1 Myofibroblast 特异性统计
myo <- which(ann$cell_subtype == "Myofibroblasts")
nonmyo_stromal <- which(ann$cell_type == "Stromal cells" & ann$cell_subtype != "Myofibroblasts")
epi <- which(ann$cell_type == "Epithelial cells")
tim_all <- as.numeric(expr_lin["TIMP1", ])
wil_myo <- suppressWarnings(wilcox.test(tim_all[myo], tim_all[epi])$p.value)
message(sprintf("✓ TIMP1 Myofibroblast vs Epithelial: Wilcoxon p = %.2e (mean TPM %.1f vs %.1f)",
                wil_myo, mean(tim_all[myo]), mean(tim_all[epi])))

# ---- 6. Fig 7C: 5 模型基因 dotplot (mean TPM 大小 + %positive 颜色) ----
mod_df <- do.call(rbind, lapply(MODEL_GENES, function(g) {
  data.frame(gene = g, cell_type = rep(ct_levels, each = 1),
             mean = ct_res$mean[g, ct_levels],
             pct  = ct_res$pct[g, ct_levels])
}))
mod_df$cell_type <- factor(mod_df$cell_type, levels = rev(ct_levels))
p_dot <- ggplot(mod_df, aes(cell_type, reorder(gene, -match(gene, MODEL_GENES)))) +
  geom_point(aes(size = pct, color = mean)) +
  scale_color_gradient2(low = "#1B998B", mid = "white", high = "#D7263D", midpoint = 0) +
  scale_size_area(max_size = 8) +
  labs(x = NULL, y = NULL, size = "% positive", color = "Mean TPM",
       title = "5-gene FARGs model: cellular localization") +
  theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_dot, "Fig7C_model_genes_dotplot", w = 7, h = 4.5)

# ---- 7. Fig 7D: 基质细胞内 TIMP1 vs CAF 标志 单细胞相关 ----
strom <- which(ann$cell_type == "Stromal cells")
cor_genes <- intersect(CAF_MARKERS, rownames(expr_lin))
sc_cor <- sapply(cor_genes, function(g) {
  suppressWarnings(cor(tim_all[strom], as.numeric(expr_lin[g, strom]), method = "spearman"))
})
sc_df <- data.frame(marker = names(sc_cor), rho = as.numeric(sc_cor)) %>% arrange(desc(rho))
save_table(sc_df, "Table7c_scRNA_TIMP1_CAF_cor_stromal")
p_sc <- ggplot(sc_df, aes(reorder(marker, rho), rho, fill = rho)) +
  geom_col() + coord_flip() +
  scale_fill_gradient2(low = "#1B998B", mid = "white", high = "#D7263D", midpoint = 0) +
  labs(x = NULL, y = "Spearman rho (stromal cells)",
       title = "TIMP1 co-expression with CAF markers (single-cell)") +
  theme_classic(base_size = 11) + theme(legend.position = "none")
save_fig(p_sc, "Fig7D_TIMP1_CAF_scRNA", w = 5, h = 5)

# ---- 8. ferroptosis / senescence 程序在 CAF(Myofibroblast) 中的活性 ----
prog_genes <- intersect(c(FERRO_GENES, SENS_GENES), rownames(expr_lin))
prog_ct <- agg_by_ct(prog_genes, ann$cell_type)
caf_vs_other <- data.frame(
  gene = prog_genes,
  Myofibroblast_mean = prog_ct$mean[prog_genes, "Stromal cells"],  # 用 subtype 更精确
  stringsAsFactors = FALSE
)
# 用 Myofibroblast subtype 精确均值
myo_mean <- rowMeans(expr_lin[prog_genes, myo, drop = FALSE])
strom_mean <- rowMeans(expr_lin[prog_genes, which(ann$cell_type=="Stromal cells"), drop = FALSE])
epi_mean <- rowMeans(expr_lin[prog_genes, epi, drop = FALSE])
prog_tab <- data.frame(gene = prog_genes,
                       Myofibroblast = myo_mean,
                       AllStromal = strom_mean,
                       Epithelial = epi_mean) %>%
  mutate(Myofib_vs_Epi_logFC = log2((Myofibroblast + 1e-3)/(Epithelial + 1e-3)))
save_table(prog_tab, "Table7d_ferro_senescence_program_CAF")

message("=== 14_scRNA_celltype.R 完成 ===")

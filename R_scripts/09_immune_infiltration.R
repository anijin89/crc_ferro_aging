# ============================================================
# 09_immune_infiltration.R — ssGSEA 免疫浸润 + 免疫检查点 + ESTIMATE-like
# ============================================================
# GSVA 包已装 (ssGSEA)；ESTIMATE/immunedeconv 因 R 4.6 兼容性无法安装
# 采用双策略：
#   1) ssGSEA + 手工 immune cell marker gene sets (16 类，覆盖 T/B/NK/髓系/基质)
#   2) 免疫检查点表达 (24 个 checkpoint genes)
#   3) ESTIMATE 打分法：手工 STROMAL_SCORE / IMMUNE_SCORE marker 平均 z
# 输出：Fig5A (boxplot High/Low × 16 cell types)
#       Fig5B (heat correlation score × cell types)
#       Fig5C (checkpoint barplot)
#       Table5_immune_infiltration.csv, Table5b_checkpoint_correlation.csv
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
suppressPackageStartupMessages(library(GSVA))

# ---- 载入 log2(TPM+1) + FARGs score ----
expr_dt <- data.table::fread(file.path(DATA_RAW, "TCGA_COADREAD_TPM_matrix.tsv"))
gene_ids <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_ids; rm(expr_dt, gene_ids)
colnames(expr_mat) <- sub("^\\.\\d+", "", colnames(expr_mat))

sinfo <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_sample_types.csv"),
                  stringsAsFactors = FALSE)
pt <- sinfo %>% filter(sample_type == "Primary Tumor") %>%
  filter(sample_id %in% colnames(expr_mat))
pt_dedup <- pt %>% group_by(patient_barcode) %>%
  summarise(sample_id = min(sample_id), .groups = "drop") %>% as.data.frame()

expr_pt <- log2(expr_mat[, pt_dedup$sample_id, drop = FALSE] + 1)
rownames(expr_pt) <- rownames(expr_mat)
message(sprintf("Primary tumor matrix (patient-level): %d genes × %d patients",
                nrow(expr_pt), ncol(expr_pt)))

model <- readRDS(file.path(RES_DIR, "FARGs_score_TCGA.rds"))
# 对齐到有效样本
common_patients <- intersect(pt_dedup$patient_barcode, model$patient_ids)
message(sprintf("Aligned samples with FARGs score: %d / %d",
                length(common_patients), ncol(expr_pt)))

score_map <- setNames(model$score, model$patient_ids)
sample_map <- setNames(pt_dedup$sample_id, pt_dedup$patient_barcode)[common_patients]
expr_pt <- expr_pt[, sample_map, drop = FALSE]
score   <- unname(score_map[common_patients])
names(score) <- unname(sample_map)                  # UUID keys
grp     <- ifelse(score > median(score), "High", "Low")
names(grp) <- names(score)

# ---- 1. ssGSEA 免疫浸润评分（16 cell types）----
immune_gene_sets <- list(
  CD8_T         = c("CD8A","CD8B","GZMB","PRF1","IFNG","CXBR3","KLRB1","LAG3"),
  CD4_T         = c("CD4","IL4","IL5","IL13","CD28","ICOS","IL2","TNFRSF4"),
  Treg          = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","TNFRSF4","CCR8"),
  Th17          = c("RORC","IL17A","IL17F","IL22","CCR6","IL23R"),
  NK            = c("NCAM1","KLRD1","NKG7","GNLY","KLRF1","PRF1","GZMA","GZMB"),
  B_cell        = c("CD19","MS4A1","CD79A","CD79B","PAX5","IGHM","BLK"),
  Plasma        = c("MZB1","PRDM1","SDC1","XBP1","JCHAIN","IGHG1"),
  Monocyte      = c("CD14","FCGR3A","FCGR1A","VCAN","S100A8","S100A9"),
  Macrophage_M1 = c("CD68","iNOS","NOS2","IL1B","TNF","CXCL9","CXCL10","CD80","CD86"),
  Macrophage_M2 = c("CD163","MRC1","MSR1","CCL18","IL10","CD209","STAB1"),
  DC            = c("ITGAX","CD1C","FCER1A","CLEC9A","CLEC10A","HLA-DRA"),
  Neutrophil    = c("CSF3R","FCGR3B","CEACAM8","CXCR2","ELANE","MPO"),
  Mast          = c("KIT","CPA3","TPSAB1","TPSB2","MS4A2","HDC"),
  CAF           = c("ACTA2","FAP","THY1","PDGFRB","COL1A1","COL1A2","FN1","VIM"),
  Endothelial   = c("PECAM1","CDH5","VWF","KDR","TEK","CLDN5","ESAM"),
  Tumor_Epith   = c("EPCAM","KRT8","KRT18","KRT19","CDH1","CLDN4")
)
# 过滤基因不在矩阵的
immune_gene_sets <- lapply(immune_gene_sets, function(g) intersect(g, rownames(expr_pt)))
immune_gene_sets <- immune_gene_sets[sapply(immune_gene_sets, length) >= 3]
message(sprintf("Using %d immune gene sets after filter", length(immune_gene_sets)))

set.seed(SEED)
ssp <- GSVA::ssgseaParam(exprData = as.matrix(expr_pt),
                         geneSets = immune_gene_sets,
                         normalize = TRUE)
ssgsea_res <- GSVA::gsva(ssp)
message(sprintf("ssGSEA result: %d cell types × %d samples",
                nrow(ssgsea_res), ncol(ssgsea_res)))

immune_df <- as.data.frame(t(ssgsea_res))
cell_types <- colnames(immune_df)   # 16 cell type names
immune_df$patient <- rownames(immune_df)
immune_df$group <- factor(grp[rownames(immune_df)], levels = c("Low","High"))
immune_df$FARGs_score <- score[rownames(immune_df)]
save_table(immune_df, "Table5_immune_infiltration")

# ---- 2. High vs Low 差异检验（Wilcoxon）----
test_results <- sapply(cell_types, function(ct) {
  suppressWarnings(wilcox.test(immune_df[[ct]] ~ immune_df$group)$p.value)
})
test_df <- data.frame(cell_type = names(test_results), p_value = as.numeric(test_results),
                      stringsAsFactors = FALSE) %>%
  arrange(p_value) %>%
  mutate(sig = ifelse(p_value < 0.001, "***",
              ifelse(p_value < 0.01,  "**",
              ifelse(p_value < 0.05,  "*", "ns"))))
save_table(test_df, "Table5a_immune_celltype_diff")
message("Immune cell type significance (High vs Low):")
print(head(test_df, 15))

# ---- 3. Fig 5A: boxplot top 6 cell types ----
top_ct <- test_df$cell_type[1:min(6, nrow(test_df))]
box_df <- immune_df %>% select(all_of(c(top_ct, "group"))) %>%
  tidyr::pivot_longer(-group, names_to = "cell_type", values_to = "score")
box_df$cell_type <- factor(box_df$cell_type, levels = rev(top_ct))
p_box <- ggplot(box_df, aes(group, score, fill = group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Low" = unname(COL_RISK["Low"]),
                               "High" = unname(COL_RISK["High"]))) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3) +
  labs(x = NULL, y = "ssGSEA score",
       title = "Immune infiltration: FARGs score High vs Low") +
  theme_classic(base_size = 10) + theme(legend.position = "none")
save_fig(p_box, "Fig5A_immune_boxplot", w = 10, h = 7)

# ---- 4. Fig 5B: 相关性 heat (FARGs score vs 各细胞)----
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(list(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  list(rho = unname(ct$estimate), p = ct$p.value)
}
cor_list <- lapply(cell_types, function(ct)
  safe_cor(immune_df[[ct]], immune_df$FARGs_score))
names(cor_list) <- cell_types
cor_vec <- sapply(cor_list, `[[`, "rho")
p_cor_vec <- sapply(cor_list, `[[`, "p")
cor_df <- data.frame(cell_type = cell_types,
                     rho = as.numeric(cor_vec),
                     p = as.numeric(p_cor_vec),
                     stringsAsFactors = FALSE) %>%
  filter(is.finite(rho)) %>% arrange(desc(rho))
cor_df$sig <- ifelse(cor_df$p < 0.001, "***",
             ifelse(cor_df$p < 0.01,  "**",
             ifelse(cor_df$p < 0.05,  "*", "")))
save_table(cor_df, "Table5b_FARGs_immune_correlation")
p_cor_bar <- ggplot(cor_df, aes(reorder(cell_type, rho), rho, fill = rho)) +
  geom_col() + coord_flip() +
  geom_text(aes(label = sig), hjust = -0.1, size = 3.5) +
  scale_fill_gradient2(low = "#1B998B", mid = "white", high = "#D7263D", midpoint = 0) +
  labs(x = NULL, y = "Spearman rho (with FARGs score)",
       title = "FARGs score vs immune cell type correlation") +
  theme_classic(base_size = 11) + theme(legend.position = "none")
save_fig(p_cor_bar, "Fig5B_immune_correlation", w = 6, h = 6)

# ---- 5. Fig 5C: 24 免疫检查点基因相关 ----
checkpoints <- c("PD-1"="PDCD1","PD-L1"="CD274","PD-L2"="PDCD1LG2",
                 "CTLA-4"="CTLA4","LAG3"="LAG3","TIM-3"="HAVCR2",
                 "TIGIT"="TIGIT","GITR"="TNFRSF18","4-1BB"="TNFRSF9",
                 "ICOS"="ICOS","CD80"="CD80","CD86"="CD86",
                 "CD27"="CD27","BTLA"="BTLA","VISTA"="VSIR",
                 "CD40"="CD40","CD40L"="CD40LG","B7-H3"="CD276",
                 "SIGLEC15"="SIGLEC15","TGFBR1"="TGFBR1",
                 "TNFRSF4"="TNFRSF4","TNFRSF14"="TNFRSF14",
                 "IDO1"="IDO1","CD70"="CD27")
cp_avail <- intersect(checkpoints, rownames(expr_pt))
cp_rows <- lapply(cp_avail, function(g) {
  v <- as.numeric(expr_pt[g, rownames(immune_df)])
  r <- safe_cor(immune_df$FARGs_score, v)
  data.frame(gene = g, rho = r$rho, p = r$p, stringsAsFactors = FALSE)
})
cp_df <- do.call(rbind, cp_rows)
cp_df$sig <- ifelse(cp_df$p < 0.001, "***",
             ifelse(cp_df$p < 0.01,  "**",
             ifelse(cp_df$p < 0.05,  "*", "")))
cp_df <- cp_df %>% arrange(desc(rho))
save_table(cp_df, "Table5c_checkpoint_correlation")
p_cp <- ggplot(cp_df, aes(reorder(gene, rho), rho, fill = rho)) +
  geom_col() + coord_flip() +
  geom_text(aes(label = sig), hjust = ifelse(cp_df$rho > 0, -0.1, 1.1), size = 3) +
  scale_fill_gradient2(low = "#1B998B", mid = "white", high = "#D7263D", midpoint = 0) +
  labs(x = NULL, y = "Spearman rho (with FARGs score)",
       title = "Immune checkpoint genes vs FARGs score") +
  theme_classic(base_size = 10) + theme(legend.position = "none")
save_fig(p_cp, "Fig5C_checkpoint_correlation", w = 6, h = 8)

# ---- 6. ESTIMATE-like 分数（手工 marker 集平均 z）----
stromal_markers <- intersect(c("COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1",
                                "FN1","LUM","POSTN","SPARC","TIMP1","TIMP2","THBS1",
                                "MMP2","MMP14","ACTA2","TAGLN","PDGFRB","CSPG4","VIM"),
                              rownames(expr_pt))
immune_markers <- unlist(immune_gene_sets, use.names = FALSE) |> unique() |>
  (\(x) intersect(x, rownames(expr_pt)))()
z_row <- function(m) (m - rowMeans(m))/apply(m, 1, sd)
z_expr <- z_row(expr_pt)
stromal_score <- colMeans(z_expr[stromal_markers, , drop = FALSE])
immune_score  <- colMeans(z_expr[immune_markers, , drop = FALSE])
est_df <- data.frame(patient = colnames(expr_pt),
                     stromal = stromal_score, immune = immune_score,
                     tumor_purity = stromal_score + immune_score,
                     FARGs = score[colnames(expr_pt)],
                     group = factor(grp[colnames(expr_pt)], levels=c("Low","High")))
save_table(est_df, "Table5d_estimate_like_scores")
p_est <- ggplot(est_df, aes(stromal, immune, color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.6) +
  scale_color_manual(values = c("Low" = unname(COL_RISK["Low"]),
                                "High" = unname(COL_RISK["High"]))) +
  annotate("text", x = -1.5, y = 2.5,
           label = {
             .ct <- suppressWarnings(cor.test(est_df$stromal, est_df$immune,
                                              method = "spearman"))
             # p 值极小（~1e-90）时 %.2e 会下溢成 "0.00e+00"，改用安全格式化
             sprintf("rho = %.2f, p %s", .ct$estimate,
                     if (is.na(.ct$p.value)) "= NA" else if (.ct$p.value < 1e-300)
                       "< 1e-300" else sprintf("= %.2g", .ct$p.value))
           },
           size = 3.5) +
  labs(x = "Stromal score (ESTIMATE-like)", y = "Immune score (ESTIMATE-like)",
       title = "Tumor microenvironment scores", color = "Risk") +
  theme_classic(base_size = 11)
save_fig(p_est, "Fig5D_estimate_scatter", w = 6, h = 5)

message("=== 09_immune_infiltration.R 完成 ===")

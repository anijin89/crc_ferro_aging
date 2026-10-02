# ============================================================
# 04_DEG_intersection.R — DEG 鉴定 + FARG 交集 + 相关性热图
# ============================================================
# 输入：TPM_matrix.tsv + sample_types.csv (仅 Primary Tumor + Solid Tissue Normal)
# 输出：Fig 2A/1B/1C/1D；DEG.rds；FARG_DEG_intersection.rds；TPM_log2.rds
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
load_core_pkgs()
suppressPackageStartupMessages({ library(ggvenn); library(data.table) })

# ---- 1. 载入表达矩阵与样本注释 ----
expr_dt <- data.table::fread(file.path(DATA_RAW, "TCGA_COADREAD_TPM_matrix.tsv"))
sinfo   <- read.csv(file.path(DATA_RAW, "TCGA_COADREAD_sample_types.csv"),
                    stringsAsFactors = FALSE)
gene_ids  <- expr_dt[[1]]
expr_mat  <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_ids
rm(expr_dt, gene_ids)
colnames(expr_mat) <- sub("^\\.\\d+", "", colnames(expr_mat))  # safety
message(sprintf("Raw TPM: %d genes × %d samples", nrow(expr_mat), ncol(expr_mat)))

# ---- 2. 严格过滤样本类型：只保留 Primary Tumor + Solid Tissue Normal ----
keep_types <- c("Primary Tumor", "Solid Tissue Normal")
sinfo_keep <- sinfo %>% filter(sample_type %in% keep_types)
excluded   <- sinfo %>% filter(!sample_type %in% keep_types)
message(sprintf("Filtered to %d samples (dropped %d Met/Recurrent/Other)",
                nrow(sinfo_keep), nrow(excluded)))
if (nrow(excluded) > 0) {
  message("  Excluded sample types:")
  print(table(excluded$sample_type))
}

# 只保留表达矩阵中真实存在的样本
avail <- intersect(sinfo_keep$sample_id, colnames(expr_mat))
sinfo_keep <- sinfo_keep %>% filter(sample_id %in% avail)
stopifnot(nrow(sinfo_keep) == length(avail))

expr <- expr_mat[, avail, drop = FALSE]
storage.mode(expr) <- "numeric"
expr <- log2(expr + 1)  # log2(TPM+1)

# 组别
samp_type <- sinfo_keep$sample_type[match(colnames(expr), sinfo_keep$sample_id)]
group <- factor(ifelse(samp_type == "Primary Tumor", "Tumor", "Normal"),
                levels = c("Normal", "Tumor"))
message(sprintf("Tumor=%d, Normal=%d", sum(group=="Tumor"), sum(group=="Normal")))

# ---- 3. 保存样本-患者映射（不保存大矩阵，下游脚本可复用同一过滤规则重算）----
sample_patient_map <- data.frame(
  sample_id       = colnames(expr),
  sample_type     = samp_type,
  patient_barcode = sinfo_keep$patient_barcode[match(avail, sinfo_keep$sample_id)],
  stringsAsFactors = FALSE
)
save_obj(sample_patient_map, "TCGA_COADREAD_sample_patient_map")

# ---- 4. 低表达过滤（至少 10 个样本 TPM>1） ----
tpm_sub <- expr_mat[, avail, drop = FALSE]
keep_g  <- rowSums(tpm_sub > 1) >= 10
expr    <- expr[keep_g, ]
rm(tpm_sub)
message(sprintf("After low-expression filter: %d genes", nrow(expr)))

# ---- 5. limma DEG ----
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr, design)
cont.matrix  <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg  <- topTable(fit2, number = Inf, sort.by = "P")
deg$gene <- rownames(deg)
deg$sig  <- ifelse(deg$adj.P.Val < DEG_FDR & abs(deg$logFC) > DEG_LOGFC,
                   ifelse(deg$logFC > 0, "Up", "Down"), "NS")
message("DEG summary:"); print(table(deg$sig))
save_obj(deg, "TCGA_COADREAD_DEG")
save_table(deg, "TableS1_DEG_full")

# ---- 6. 火山图 (Fig 2A) ----
p_volcano <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "#D7263D", "Down" = "#1B998B", "NS" = "grey70")) +
  geom_vline(xintercept = c(-DEG_LOGFC, DEG_LOGFC), linetype = "dashed") +
  geom_hline(yintercept = -log10(DEG_FDR), linetype = "dashed") +
  geom_text_repel(data = subset(deg, sig != "NS" & abs(logFC) > 2)[1:20, ],
                  aes(label = gene), size = 3, max.overlaps = 30) +
  labs(x = "log2 Fold Change (Tumor vs Normal)",
       y = "-log10 (FDR)",
       title = sprintf("TCGA-COADREAD DEGs (PT=%d vs STN=%d)",
                       sum(group=="Tumor"), sum(group=="Normal"))) +
  theme_classic(base_size = 12) + theme(legend.position = "right")
save_fig(p_volcano, "Fig2A_volcano", w = 8, h = 6)

# ---- 7. 热图：Top 50 DEG (Fig 2B) ----
top50 <- deg %>% filter(sig != "NS") %>%
  arrange(desc(abs(logFC))) %>% head(50) %>% pull(gene)
mat <- expr[top50, ]; mat <- t(scale(t(mat)))
mat[mat > 2] <- 2; mat[mat < -2] <- -2
png(file.path(FIG_DIR, "Fig2B_FARGDEG_heatmap.png"),
    width = 10, height = 10, units = "in", res = FIG_DPI)
pheatmap(mat, annotation_col = data.frame(Group = group[colnames(mat)]),
         show_colnames = FALSE,
         color = colorRampPalette(c("#1B998B", "white", "#D7263D"))(100),
         main = "Top 50 DEGs")
dev.off()

# ---- 8. FARG ∩ DEG (Fig 2C Venn) ----
fargs <- readRDS(file.path(RES_DIR, "FARG_genes.rds"))
deg_sig <- deg %>% filter(sig != "NS") %>% pull(gene)
fargs_in_deg <- intersect(fargs, deg_sig)
message(sprintf("✓ FARG-DEGs: %d (|FARG|=%d, |DEG_sig|=%d)",
                length(fargs_in_deg), length(fargs), length(deg_sig)))

lst <- list("FARGs" = fargs, "DEGs" = deg_sig)
p_venn <- ggvenn(lst, fill_color = c("#1B998B", "#D7263D"),
                 stroke_size = 0.5, set_name_size = 4, text_size = 4) +
  ggtitle("Ferro-aging-related DEGs")
save_fig(p_venn, "Fig2C_venn", w = 6, h = 5)

save_obj(fargs_in_deg, "FARG_DEG_intersection")
farg_deg_meta <- deg[deg$gene %in% fargs_in_deg, ]
save_table(farg_deg_meta, "TableS2_FARG_DEG")

# ---- 9. 相关性热图 (Fig 2D) ----
expr_farg <- expr[fargs_in_deg, ]
cor_mat <- cor(t(expr_farg), method = "pearson")
png(file.path(FIG_DIR, "Fig2D_correlation.png"),
    width = 9, height = 8, units = "in", res = FIG_DPI)
pheatmap(cor_mat,
         color = colorRampPalette(c("#1B998B", "white", "#D7263D"))(100),
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 6, clustering_method = "ward.D2",
         main = "Pearson correlation among FARG-DEGs")
dev.off()

message("=== 04_DEG_intersection.R 完成 ===")

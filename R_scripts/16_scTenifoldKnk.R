# ============================================================
# 16_scTenifoldKnk.R — TIMP1 虚拟敲除（核心创新分析）
# ============================================================
# 数据：GSE132465 Myofibroblast(CAF) 细胞 raw UMI counts
#       (33,695 genes × 1,155 cells, 已抽取)
# 理由：TIMP1 主要来源于 CAF/Myofibroblast（Fig7 证实），
#       故在 CAF 中做虚拟敲除最贴合机制叙事。
# 方法：scTenifoldKnk (CRAN, cailab-tamu) —
#       构建 scGRN (PC regression + adaptive LASSO) → 网络扰动
#       → dRegulation 给出下游差异调控基因
# 输出：Fig8A volcano, Fig8B GO, Fig8C KEGG, Table8a/b
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
suppressPackageStartupMessages(library(scTenifoldKnk))

HUB <- "TIMP1"
CELLTYPE <- "Myofibroblasts (CAF)"

# ---- 1. 读取 CAF UMI 矩阵 ----
dt <- data.table::fread(file.path(DATA_PROC, "GSE132465_myofibroblast_UMI.tsv"))
genes <- dt[[1]]
cm <- as.matrix(dt[, -1, with = FALSE])
rownames(cm) <- genes; rm(dt, genes)
storage.mode(cm) <- "integer"
message(sprintf("✓ raw matrix: %d genes × %d cells", nrow(cm), ncol(cm)))

# ---- 2. 基因过滤（去低表达 + 去 MT/rRNA + 取 top 800 高变基因加速）----
detected <- rowSums(cm > 0)
keep <- detected >= 0.05 * ncol(cm) & rowSums(cm) >= 10
cm <- cm[keep, ]
cm <- cm[!grepl("^MT-|^RP[SL]|^RNU", rownames(cm)), ]
# log-CPM 方差选 top 800 高变基因（大幅加速 GRN 构建）
libsize <- colSums(cm)
logcpm <- sweep(log2(cm + 1), 2, log2(libsize + 1), FUN = "-")
gvar <- rowMeans(logcpm^2) - rowMeans(logcpm)^2
top <- order(gvar, decreasing = TRUE)[1:min(800, length(gvar))]
top <- unique(c(top, which(rownames(cm) == HUB)))
cm <- cm[sort(top), ]
stopifnot(HUB %in% rownames(cm))
message(sprintf("✓ filtered matrix: %d genes × %d cells (TIMP1 kept)",
                nrow(cm), ncol(cm)))

# ---- 3. 虚拟敲除（结果持久化，避免重复计算）----
knk_rds <- file.path(RES_DIR, "knk_raw_TIMP1.rds")
if (file.exists(knk_rds)) {
  message("✓ 发现已缓存 knk_raw_TIMP1.rds，直接读取")
  knk <- readRDS(knk_rds)
} else {
  message("⏳ 构建 scGRN + 网络扰动（nc_nNet=3, nc_nCells=400）...")
  set.seed(SEED)
  t0 <- Sys.time()
  knk <- scTenifoldKnk(countMatrix = cm, gKO = HUB,
                       transcriptomeWide = FALSE,
                       nc_nNet = 3, nc_nCells = 400, nCores = 4)
  message(sprintf("✓ scTenifoldKnk done in %.1f min",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(knk, knk_rds)
}

# ---- 4. 差异调控基因（scTenifoldKnk 已返回 knk$diffRegulation）----
deg <- as.data.frame(knk$diffRegulation)
# Z = 带符号 z-score（方向/效应量）; FC = 卡方统计量（幅值）; p.adj = FDR
deg <- deg %>%
  rename(effect = Z, FDR = p.adj) %>%
  mutate(effect = as.numeric(effect), FDR = as.numeric(FDR),
         direction = ifelse(FDR < 0.05 & effect >  1, "Up",
                     ifelse(FDR < 0.05 & effect < -1, "Down", "NS")))
message("deg columns: ", paste(colnames(deg), collapse = ", "))
save_table(deg, "Table8a_scTenifoldKnk_TIMP1_full")

sig <- deg %>% filter(direction != "NS") %>% arrange(FDR)
message(sprintf("✓ KO-driven DEGs (FDR<0.05, |Z|>1): %d (Up=%d, Down=%d)",
                nrow(sig), sum(sig$direction=="Up"), sum(sig$direction=="Down")))
save_table(sig, "Table8b_scTenifoldKnk_TIMP1_significant")

# ---- 5. Fig 8A: 火山图 (x = signed Z, y = -log10 FDR) ----
p_lab <- sig %>% arrange(FDR) %>% head(25)
p_vol <- ggplot(deg, aes(effect, -log10(FDR + 1e-300), color = direction)) +
  geom_point(alpha = 0.65, size = 1.5) +
  scale_color_manual(values = c("Up"="#D7263D","Down"="#1B998B","NS"="grey70"),
                     name = "Direction",
                     guide = guide_legend(override.aes = list(alpha = 1, size = 2.6))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # show.legend = FALSE 是必须的：ggrepel 的图例键画的是一个小写 "a"，
  # 而它默认继承上一层 aes(color = direction)，于是红色的 "a" 正好压在
  # 红色 "Up" 图例点上，看起来像红点被一个字母盖住。
  ggrepel::geom_text_repel(data = p_lab, aes(label = gene), size = 3,
                           max.overlaps = 40, show.legend = FALSE) +
  labs(x = sprintf("Regulation effect (signed Z, %s KO vs WT)", HUB),
       y = "-log10 (FDR)",
       title = sprintf("Virtual knockout of %s in %s", HUB, CELLTYPE)) +
  theme_classic(base_size = 12)
save_fig(p_vol, "Fig8A_TIMP1_KO_volcano", w = 8, h = 6)

# ---- 6. 下游扰动程序（hypothesis-generating：按 perturbation distance 排名）----
# 说明：scTenifoldKnk 的 FDR 检验中 KO 基因自身主导流形，其余基因 FDR 均不显著；
#       这是该方法的固有特性。因此下游分析采用"网络扰动幅度排名"的假设生成策略，
#       取除 TIMP1 外 perturbation distance 最大的前 150 个基因作为 TIMP1 调控程序。
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
})
topN <- deg %>% filter(gene != HUB) %>%
  arrange(desc(distance)) %>% head(150)
save_table(topN, "Table8c_TIMP1_KO_top_perturbed")
message(sprintf("✓ top-perturbed downstream genes: %d (distance %.1e ~ %.1e)",
                nrow(topN), min(topN$distance), max(topN$distance)))

# Fig 8B: top 25 扰动基因排名条形图
top25 <- head(topN, 25)
p_rank <- ggplot(top25 %>% mutate(gene = factor(gene, levels = rev(gene))),
                 aes(gene, distance, fill = effect)) +
  geom_col() + coord_flip() +
  scale_fill_gradient2(low = "#1B998B", mid = "grey85", high = "#D7263D", midpoint = 0) +
  labs(x = NULL, y = "Network perturbation distance",
       fill = "signed effect",
       title = sprintf("Top 25 downstream genes perturbed by %s KO", HUB)) +
  theme_classic(base_size = 11)
save_fig(p_rank, "Fig8B_TIMP1_KO_top_perturbed", w = 7, h = 7)

# GO / KEGG on top-perturbed set
eg <- bitr(topN$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
bg_eg <- bitr(rownames(cm), fromType="SYMBOL", toType="ENTREZID",
              OrgDb=org.Hs.eg.db, drop=TRUE)
message(sprintf("✓ %d/%d top genes mapped to ENTREZID", nrow(eg), nrow(topN)))
ego <- enrichGO(gene=eg$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
                readable=TRUE, universe=bg_eg$ENTREZID)
ekegg <- tryCatch(enrichKEGG(gene=eg$ENTREZID, organism="hsa",
                pvalueCutoff=0.05, qvalueCutoff=0.2), error=function(e) NULL)

if (!is.null(ego) && nrow(ego) > 0) {
  p_go <- barplot(ego, showCategory=15,
                  title="GO BP: programs perturbed by TIMP1 KO (top-ranked)") +
    theme(axis.text.y = element_text(size=9))
  save_fig(p_go, "Fig8C_TIMP1_KO_GO", w=9, h=7)
  save_table(as.data.frame(ego), "Table8d_TIMP1_KO_GO")
} else message("  GO: no enriched term")
if (!is.null(ekegg) && nrow(ekegg) > 0) {
  p_kegg <- barplot(ekegg, showCategory=15,
                    title="KEGG: pathways perturbed by TIMP1 KO") +
    theme(axis.text.y = element_text(size=9))
  save_fig(p_kegg, "Fig8D_TIMP1_KO_KEGG", w=9, h=7)
  save_table(as.data.frame(ekegg), "Table8e_TIMP1_KO_KEGG")
} else message("  KEGG: no enriched pathway")

# ---- 7. 保存对象 ----
save_obj(list(deg=deg, sig=sig, topN=topN, ego=ego, ekegg=ekegg,
              n_cells=ncol(cm), n_genes=nrow(cm)),
         "scTenifoldKnk_TIMP1")

message("=== 16_scTenifoldKnk.R 完成 ===")

# ============================================================
# 05_GO_KEGG.R — FARG-DEG 的 GO-BP / KEGG 富集（Fig 3A–2B + Table S3）
# ============================================================
# 输入：results/FARG_DEG_intersection.rds  （26 个 FARG-DEG，由 04 产出）
# 输出：figures/Fig3A_GO_BP.png/pdf, Fig3B_KEGG.png/pdf
#       tables/TableS3_GO_BP.csv, TableS3_KEGG.csv
# 说明：输入基因集严格取自 04 的产出（26 个），
#       不重新定义基因集，保证 TableS2 / TableS3 / Fig2C 完全同源。
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
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
})

TOP_N <- 15          # 图中展示的通路条数
PSIG  <- 0.05        # p / q 阈值（与 clusterProfiler 默认一致）

# ---- 1. 输入基因集（严格来自 04 的产出）----
farg_deg <- readRDS(file.path(RES_DIR, "FARG_DEG_intersection.rds"))
farg_deg <- unique(as.character(farg_deg))
message(sprintf("输入 FARG-DEG: %d 个", length(farg_deg)))

# 与 Table S2 交叉核对，防止 04 被改动后本脚本静默不一致
tab_s2 <- read.csv(file.path(TAB_DIR, "TableS2_FARG_DEG.csv"), stringsAsFactors = FALSE)
if (nrow(tab_s2) != length(farg_deg) || !setequal(tab_s2$gene, farg_deg)) {
  warning(sprintf("!! FARG_DEG_intersection.rds (%d) 与 TableS2 (%d) 不一致，请先重跑 04",
                  length(farg_deg), nrow(tab_s2)))
} else {
  message("✓ 与 TableS2_FARG_DEG.csv 一致（", nrow(tab_s2), " 个基因）")
}

# ---- 2. SYMBOL → ENTREZID ----
id_map <- suppressWarnings(
  bitr(farg_deg, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
)
unmapped <- setdiff(farg_deg, id_map$SYMBOL)
if (length(unmapped)) message("  未映射到 Entrez 的基因: ", paste(unmapped, collapse = ", "))
message(sprintf("映射成功: %d / %d", nrow(id_map), length(farg_deg)))

# ---- 3. GO Biological Process ----
ego <- enrichGO(gene          = farg_deg,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = PSIG,
                qvalueCutoff  = PSIG,
                readable      = TRUE)
ego_df <- as.data.frame(ego)
if (nrow(ego_df) == 0) stop("GO-BP 无显著结果，请检查输入基因集")
# 去掉 ONTOLOGY 列，保持与既有 TableS3 列结构一致
ego_out <- ego_df[, c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
                      "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue",
                      "geneID", "Count")]
save_table(ego_out, "TableS3_GO_BP")
message(sprintf("GO-BP 显著条目: %d（展示前 %d）", nrow(ego_out), min(TOP_N, nrow(ego_out))))

p_go <- ego_out %>% slice_max(Count, n = TOP_N, with_ties = FALSE) %>%
  ggplot(aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  labs(x = NULL, y = "Count", title = "GO Biological Process", fill = "p.adjust") +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 13))
save_fig(p_go, "Fig3A_GO_BP", w = 8, h = 6)

# ---- 4. KEGG pathway ----
ekg <- enrichKEGG(gene          = id_map$ENTREZID,
                  organism      = "hsa",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = PSIG,
                  qvalueCutoff  = PSIG)
if (is.null(ekg) || nrow(as.data.frame(ekg)) == 0) stop("KEGG 无显著结果，请检查输入基因集")
ekg_out <- as.data.frame(ekg)
save_table(ekg_out, "TableS3_KEGG")
message(sprintf("KEGG 显著通路: %d（展示前 %d）", nrow(ekg_out), min(TOP_N, nrow(ekg_out))))

p_kegg <- ekg_out %>% slice_max(Count, n = TOP_N, with_ties = FALSE) %>%
  ggplot(aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  labs(x = NULL, y = "Count", title = "KEGG Pathway", fill = "p.adjust") +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 13))
save_fig(p_kegg, "Fig3B_KEGG", w = 8, h = 6)

# ---- 5. 控制台摘要 ----
cat("\n=== GO-BP top 10 ===\n"); print(head(ego_out[, c("Description", "Count", "p.adjust")], 10))
cat("\n=== KEGG top 10 ===\n");  print(head(ekg_out[, c("Description", "Count", "p.adjust")], 10))

message("=== 05_GO_KEGG.R 完成 ===")

# ============================================================
# 12_Druggability_DGIdb.R — 模型基因可成药性 / 已知药物互作
#   动机: oncoPredict 所需的 GDSC 细胞系参考矩阵在本机网络不可达
#         (Sanger oncoPredict 主机公网 404/无响应, GitHub 镜像缺失)。
#         改用可访问的 DGIdb GraphQL API, 对 5 个模型基因 + hub 基因
#         TIMP1 做"可成药性 / 已知治疗药物"注释, 作为转化意义补充。
#   输出: Table9_druggability.csv  Fig10A_药物互作计数  Fig10B_抗癌药分布
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
  library(httr); library(jsonlite); library(dplyr); library(ggplot2)
})

# 分析基因 = 5 模型基因 + hub TIMP1 (去重)
score_obj <- readRDS(file.path(RES_DIR, "FARGs_score_TCGA.rds"))
model_genes <- score_obj$genes
QUERY_GENES <- unique(c(model_genes, "TIMP1"))
message(sprintf("查询基因: %s", paste(QUERY_GENES, collapse=", ")))

API <- "https://dgidb.org/api/graphql"
# DGIdb schema: genes(names:[..]) { edges { node { name
#   interactions { drug{name approved antiNeoplastic} interactionTypes{type} evidenceScore } } } }
q <- sprintf(
  'query { genes(names:[%s]) { edges { node { name interactions { drug { name approved antiNeoplastic } interactionTypes { type } evidenceScore } } } } }',
  paste0('"', QUERY_GENES, '"', collapse=","))

# ---- 稳健请求（重试 + 长超时）----
fetch_dgidb <- function(q, tries = 4) {
  for (i in seq_len(tries)) {
    res <- tryCatch(
      POST(API, body = toJSON(list(query = q), auto_unbox = TRUE),
           add_headers("Content-Type" = "application/json"), timeout(120)),
      error = function(e) { message(sprintf("  尝试 %d 网络错误: %s", i, conditionMessage(e))); NULL })
    if (!is.null(res) && status_code(res) == 200) {
      txt <- content(res, as = "text", encoding = "UTF-8")
      js  <- tryCatch(fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(js) && is.null(js$errors)) return(js$data$genes$edges)
      if (!is.null(js$errors)) message(sprintf("  尝试 %d API 报错: %s", i, js$errors$message[1]))
    } else if (!is.null(res)) {
      message(sprintf("  尝试 %d HTTP %d", i, status_code(res)))
    }
    Sys.sleep(3 * i)
  }
  stop("DGIdb 查询失败")
}

edges <- fetch_dgidb(q)

# ---- 展平 ----
rows <- list(); drug_rows <- list()
for (e in edges) {
  gene <- e$node$name
  ints <- e$node$interactions
  if (is.null(ints) || length(ints) == 0) {
    rows[[gene]] <- data.frame(gene=gene, n_int=0, n_approved=0, n_antineoplastic=0,
                               top_approved="", stringsAsFactors=FALSE)
    next
  }
  # ints 为 list<interaction>（simplifyVector=FALSE）
  drug_name <- sapply(ints, function(x) x$drug$name %||% NA_character_)
  approved  <- sapply(ints, function(x) isTRUE(as.logical(x$drug$approved)))
  anti_neo  <- sapply(ints, function(x) isTRUE(as.logical(x$drug$antiNeoplastic)))
  iScore    <- sapply(ints, function(x) { s <- x$interactionTypes; if (is.null(s)||length(s)==0) NA_character_ else paste(unlist(lapply(s, function(t) t$type)), collapse="|") })
  ev        <- sapply(ints, function(x) as.numeric(x$evidenceScore))
  ev[is.na(ev)] <- 0

  rows[[gene]] <- data.frame(
    gene = gene,
    n_int = length(ints),
    n_approved = sum(approved),
    n_antineoplastic = sum(approved & anti_neo),
    top_approved = {
      idx <- which(approved & anti_neo)
      if (length(idx) == 0) "" else {
        ord <- idx[order(-ev[idx])][1:min(3, length(idx))]
        paste(na.omit(drug_name[ord]), collapse="; ")
      }
    },
    stringsAsFactors = FALSE)

  drug_rows[[length(drug_rows)+1]] <- data.frame(
    gene=gene, drug=drug_name, approved=approved, antineoplastic=anti_neo,
    itype=iScore, evidence=ev, stringsAsFactors=FALSE)
}
drugg <- do.call(rbind, rows)
drugg$gene <- factor(drugg$gene, levels = QUERY_GENES)
drugg <- drugg %>% arrange(gene)
save_table(drugg, "Table9_druggability")
cat("\n=== Table 9: 模型基因可成药性 (DGIdb) ===\n"); print(drugg)

all_drugs <- do.call(rbind, drug_rows)
all_drugs$gene <- factor(all_drugs$gene, levels = QUERY_GENES)
save_table(all_drugs %>% filter(approved) %>% arrange(gene, desc(evidence)),
           "Table9b_approved_interactions")

# ---- Fig 10A: 每基因互作计数 (总 / 已获批 / 获批抗癌) ----
p_long <- drugg %>%
  select(gene, Total=n_int, Approved=n_approved, Approved_Antineoplastic=n_antineoplastic) %>%
  tidyr::pivot_longer(c(Total, Approved, Approved_Antineoplastic), names_to="metric", values_to="n")
p_long$metric <- factor(p_long$metric, levels=c("Total","Approved","Approved_Antineoplastic"))
pA <- ggplot(p_long, aes(gene, n, fill=metric)) +
  geom_col(position=position_dodge(0.8), width=0.72) +
  scale_fill_manual(values=c(Total="#3D5A80", Approved="#1B998B", Approved_Antineoplastic="#D7263D"),
                    name=NULL) +
  labs(x=NULL, y="Number of drug interactions",
       title="Druggability of ferro-aging model genes (DGIdb)",
       subtitle=sprintf("Model genes: %s", paste(model_genes, collapse=", "))) +
  theme_classic(base_size=12) + theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="top")
save_fig(pA, "Fig10A_druggability_counts", w=8, h=5)

# ---- Fig 10B: 已获批抗癌药 evidence 排序 (top 12) ----
top_drugs <- all_drugs %>% filter(approved & antineoplastic) %>%
  group_by(gene, drug) %>% summarise(evidence=max(evidence), .groups="drop") %>%
  arrange(desc(evidence)) %>% slice_max(order_by=evidence, n=12)
top_drugs$gene <- factor(top_drugs$gene, levels=QUERY_GENES)
if (nrow(top_drugs) > 0) {
  pB <- ggplot(top_drugs, aes(reorder(drug, evidence), evidence, color=gene)) +
    geom_point(size=3) + geom_segment(aes(x=reorder(drug, evidence), xend=reorder(drug, evidence), y=0, yend=evidence), color="grey80") +
    coord_flip() +
    labs(x=NULL, y="Evidence score", color="Gene",
         title="Approved anti-neoplastic agents linked to model genes") +
    theme_classic(base_size=12) + scale_color_brewer(palette="Dark2")
  save_fig(pB, "Fig10B_top_approved_drugs", w=8, h=6)
} else {
  message("  (无获批抗癌药互作, 跳过 Fig10B — 与 TIMP1/MMP 类直接靶向药匮乏一致)")
}

save_obj(list(drugg=drugg, all_drugs=all_drugs), "druggability_DGIdb")
cat("\n=== 脚本12 完成 ===\n")

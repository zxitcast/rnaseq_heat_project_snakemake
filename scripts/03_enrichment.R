#!/usr/bin/env Rscript
# =============================================
# 03_enrichment.R
# 功能: GO 富集分析（仅 Biological Process）
#       输入: 按 logFC 绝对值排序的前 100 个上调/下调基因
#       绘图: X轴 = Fold Enrichment, 点大小 = Count, 排序 = Fold Enrichment
#       过滤: 
#          - 只保留 BP
#          - 剔除宽泛/错误术语 (broad_terms)
#          - 只保留 Count >= 3 的 GO term
#          - 按富集倍数排序，取前 15 个
# =============================================

.libPaths(c("~/R/library", .libPaths()))
setwd("~/biotrainee/snakemake_rnaseq")

library(org.At.tair.db)
library(GO.db)
library(ggplot2)

# ---------- 1. 读取差异分析结果 ----------
res <- read.csv("results/diff_expr/DEG_all.csv", row.names = 1)

# 按 logFC 排序，取前 100 个上调和前 100 个下调基因
res_up <- res[res$logFC > 0, ]
res_down <- res[res$logFC < 0, ]

top_n <- 100
up_genes <- rownames(head(res_up[order(-res_up$logFC), ], top_n))
down_genes <- rownames(head(res_down[order(res_down$logFC), ], top_n))
bg_genes <- rownames(res)

cat("上调基因数（Top", top_n, "）:", length(up_genes), "\n")
cat("下调基因数（Top", top_n, "）:", length(down_genes), "\n")
cat("背景基因数:", length(bg_genes), "\n")

# ---------- 2. 构建基因→GO 映射表（只保留 Biological Process）---------
gene_go_all <- AnnotationDbi::select(org.At.tair.db,
                                     keys = bg_genes,
                                     columns = c("GO", "ONTOLOGY"),
                                     keytype = "TAIR")
gene_go_all <- unique(gene_go_all[!is.na(gene_go_all$GO), ])
gene_go_bp <- gene_go_all[gene_go_all$ONTOLOGY == "BP", c("TAIR", "GO")]

# ---------- 3. 超几何分布检验函数（计算 Fold Enrichment）----------
hyperGO <- function(de_genes, bg_genes, gene_go_map) {
    annotated_bg <- unique(gene_go_map$TAIR)
    de_annotated <- intersect(de_genes, annotated_bg)
    total_genes <- length(annotated_bg)
    total_de <- length(de_annotated)
    
    go_terms <- unique(gene_go_map$GO)
    results <- lapply(go_terms, function(go_term) {
        term_genes <- gene_go_map$TAIR[gene_go_map$GO == go_term]
        term_bg_count <- length(intersect(annotated_bg, term_genes))
        term_de_count <- length(intersect(de_annotated, term_genes))
        if (term_de_count < 3) return(NULL)   # 提前过滤，不计算 Count<3 的 GO term
        p_value <- phyper(term_de_count - 1, term_bg_count, total_genes - term_bg_count, total_de, lower.tail = FALSE)
        fold_enrichment <- (term_de_count / term_bg_count) / (total_de / total_genes)
        data.frame(
            GO = go_term,
            Term = Term(go_term),
            Count = term_de_count,
            Background = term_bg_count,
            FoldEnrichment = fold_enrichment,
            P.Value = p_value,
            stringsAsFactors = FALSE
        )
    })
    results <- do.call(rbind, results)
    if (is.null(results)) return(data.frame())
    results$FDR <- p.adjust(results$P.Value, "BH")
    # 按富集倍数降序排序
    results <- results[order(-results$FoldEnrichment), ]
    return(results)
}

# ---------- 4. 执行富集分析 ----------
go_up <- hyperGO(up_genes, bg_genes, gene_go_bp)
if (nrow(go_up) > 0) write.csv(go_up, "results/diff_expr/GO_up.csv", row.names = FALSE)

go_down <- hyperGO(down_genes, bg_genes, gene_go_bp)
if (nrow(go_down) > 0) write.csv(go_down, "results/diff_expr/GO_down.csv", row.names = FALSE)

# ---------- 5. 可视化气泡图（过滤宽泛/错误术语，且 Count >=3）----------
broad_terms <- c(
    "nucleus", "cytoplasm", "chloroplast", "mitochondrion",
    "extracellular region", "protein binding", "biological_process",
    "cellular_component", "molecular_function", "gene expression",
    "response to stress", "defense response",
    "regulation of transcription, DNA-templated",
    "negative regulation of egg-laying behavior"   # 动物特有错误术语
)

plot_bubble <- function(data, title, color, max_terms = 15) {
    if (nrow(data) == 0) return(NULL)
    # 剔除宽泛术语
    data <- data[!data$Term %in% broad_terms, ]
    # 已经在前面的 hyperGO 中过滤了 Count<3，但保险再加一次
    data <- data[data$Count >= 3, ]
    # 取 FoldEnrichment 最高的前 max_terms 个
    plot_df <- head(data[order(-data$FoldEnrichment), ], max_terms)
    if (nrow(plot_df) == 0) return(NULL)
    plot_df$Term <- factor(plot_df$Term, levels = rev(plot_df$Term))
    
    p <- ggplot(plot_df, aes(x = FoldEnrichment, y = Term, size = Count)) +
        geom_point(color = color, alpha = 0.8) +
        labs(title = paste(title, "(Biological Process, Top by Fold Enrichment, Count >=3)"),
             x = "Fold Enrichment (higher = more enriched)",
             y = "",
             caption = "Note: Due to lack of biological replicates, only fold enrichment and gene count are shown for biological trend. Terms with <3 genes are excluded.") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 8))
    return(p)
}

p1 <- plot_bubble(go_up, "GO Enrichment (Up-regulated)", "red")
p2 <- plot_bubble(go_down, "GO Enrichment (Down-regulated)", "blue")

if (!is.null(p1)) ggsave("results/diff_expr/GO_up_bubble.png", p1, width = 12, height = 8)
if (!is.null(p2)) ggsave("results/diff_expr/GO_down_bubble.png", p2, width = 12, height = 8)

cat("===== GO 富集分析完成! 结果保存在 results/diff_expr/ 下 =====\n")

#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
setwd("~/biotrainee/rnaseq_heat_project")

library(org.At.tair.db)
library(ggplot2)

sig_genes <- read.csv("results/diff_expr/DEG_significant.csv", row.names = 1)
res <- read.csv("results/diff_expr/DEG_all.csv", row.names = 1)
up_genes <- rownames(subset(sig_genes, logFC > 1))
down_genes <- rownames(subset(sig_genes, logFC < -1))
bg_genes <- rownames(res)

gene_kegg <- AnnotationDbi::select(org.At.tair.db, keys = bg_genes, columns = "PATH", keytype = "TAIR")
gene_kegg <- unique(gene_kegg[!is.na(gene_kegg$PATH), c("TAIR", "PATH")])
colnames(gene_kegg) <- c("TAIR", "KEGG")

# 读取映射表
kegg_names <- read.csv("ath_kegg_pathways.csv", stringsAsFactors = FALSE)

# 统一清洗ID
clean_id <- function(x) { x <- as.character(x); x <- gsub("^ath", "", x, ignore.case = TRUE); x <- gsub("^0+", "", x); x[x == ""] <- "0"; x }
kegg_names$clean_id <- clean_id(kegg_names$KEGG)
gene_kegg$clean_id <- clean_id(gene_kegg$KEGG)

# 超几何检验
hyperKEGG <- function(de_genes, bg_genes, gene_kegg_map) {
    annotated_bg <- unique(gene_kegg_map$TAIR)
    de_annotated <- intersect(de_genes, annotated_bg)
    total_genes <- length(annotated_bg)
    total_de <- length(de_annotated)
    kegg_terms <- unique(gene_kegg_map$clean_id)
    results <- lapply(kegg_terms, function(kegg_term) {
        term_genes <- gene_kegg_map$TAIR[gene_kegg_map$clean_id == kegg_term]
        term_bg_count <- length(intersect(annotated_bg, term_genes))
        term_de_count <- length(intersect(de_annotated, term_genes))
        p_value <- phyper(term_de_count - 1, term_bg_count, total_genes - term_bg_count, total_de, lower.tail = FALSE)
        orig_id <- gene_kegg_map$KEGG[gene_kegg_map$clean_id == kegg_term][1]
        data.frame(KEGG = orig_id, clean_id = kegg_term, Count = term_de_count, Background = term_bg_count, P.Value = p_value, stringsAsFactors = FALSE)
    })
    results <- do.call(rbind, results)
    results <- results[results$Count > 0, ]; results$FDR <- p.adjust(results$P.Value, "BH"); results
}

kegg_up <- hyperKEGG(up_genes, bg_genes, gene_kegg)
kegg_down <- hyperKEGG(down_genes, bg_genes, gene_kegg)

# 添加名称
kegg_up$Name <- kegg_names$Description[match(kegg_up$clean_id, kegg_names$clean_id)]
kegg_down$Name <- kegg_names$Description[match(kegg_down$clean_id, kegg_names$clean_id)]
kegg_up$Name[is.na(kegg_up$Name)] <- kegg_up$KEGG[is.na(kegg_up$Name)]
kegg_down$Name[is.na(kegg_down$Name)] <- kegg_down$KEGG[is.na(kegg_down$Name)]

cat("上调前5条通路（已匹配名称）:\n")
print(head(kegg_up[, c("KEGG", "Name")], 5))
cat("下调前5条通路:\n")
print(head(kegg_down[, c("KEGG", "Name")], 5))

kegg_up <- kegg_up[order(kegg_up$FDR), ]; kegg_down <- kegg_down[order(kegg_down$FDR), ]
write.csv(kegg_up, "results/diff_expr/KEGG_up.csv", row.names = FALSE)
write.csv(kegg_down, "results/diff_expr/KEGG_down.csv", row.names = FALSE)

# 气泡图
plot_bubble <- function(data, title, color, file) {
    if (nrow(data) == 0) return()
    sig <- subset(data, FDR < 0.05); if (nrow(sig) == 0) sig <- head(data, 15)
    plot_df <- head(sig[order(sig$FDR), ], 15)
    plot_df$Name <- factor(plot_df$Name, levels = rev(plot_df$Name))
    p <- ggplot(plot_df, aes(x = -log10(FDR), y = Name, size = Count)) + geom_point(color = color) +
        labs(title = title, x = "-log10(FDR)", y = "") + theme_minimal()
    ggsave(file, p, width = 14, height = 8); cat("已保存:", file, "\n")
}
plot_bubble(kegg_up, "KEGG Enrichment (Up-regulated)", "red", "results/diff_expr/KEGG_up_bubble.png")
plot_bubble(kegg_down, "KEGG Enrichment (Down-regulated)", "blue", "results/diff_expr/KEGG_down_bubble.png")
cat("===== 全部完成！ =====\n")

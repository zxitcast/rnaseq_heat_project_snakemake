#!/usr/bin/env Rscript
library(tximport)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(readr)
library(stringr)

# ---------------------------- 读取命令行参数 ----------------------------
args <- commandArgs(trailingOnly = TRUE)
quant_dir <- args[1]          # "results/quant"
samples <- strsplit(args[2], " ")[[1]]   # 样本列表
ref_fasta <- args[3]          # cDNA fasta 路径
bcv <- as.numeric(args[4])    # 固定 BCV，如 0.1
logFC_thr <- as.numeric(args[5])
fdr_thr <- as.numeric(args[6])
out_de_all <- args[7]
out_de_sig <- args[8]
out_volcano <- args[9]
out_heatmap <- args[10]

# ---------------------------- 构建 quant.sf 文件路径 ----------------------------
files <- file.path(quant_dir, samples, "quant.sf")
names(files) <- samples

# ---------------------------- 构建 tx2gene（适用于拟南芥 TAIR10）--------------------
fasta_lines <- read_lines(ref_fasta)
tx_lines <- fasta_lines[grepl("^>", fasta_lines)]
tx_ids_full <- str_extract(tx_lines, "(?<=>)[^ ]+")
tx_ids <- str_replace(tx_ids_full, "\\.[0-9]+$", "")   # 去掉版本号，得到基因 ID
gene_ids <- tx_ids
tx2gene <- data.frame(tx = tx_ids, gene = gene_ids, stringsAsFactors = FALSE)
tx2gene <- tx2gene[!is.na(tx2gene$gene), ]
tx2gene <- unique(tx2gene)

# ---------------------------- 导入 Salmon 定量 ----------------------------
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# ---------------------------- 构建 DGEList，过滤低表达基因 ----------------------------
# 注意：默认样本顺序为 config.yaml 中 samples 的顺序
# 本例假设第一个样本为 Control，第二个为 Treatment
group <- factor(c("Control", "Treatment"))   
y <- DGEList(counts = txi$counts, group = group)
keep <- rowSums(cpm(y) >= 1) >= 1
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# ---------------------------- edgeR exactTest（无重复，固定 BCV）------------------------
et <- exactTest(y, dispersion = bcv^2)
res <- topTags(et, n = nrow(y))$table
res <- res[order(res$PValue), ]

write.csv(res, out_de_all, row.names = TRUE)

# 筛选显著差异基因（仅用于展示，无统计效力）
sig <- subset(res, FDR < fdr_thr & abs(logFC) > logFC_thr)
write.csv(sig, out_de_sig, row.names = TRUE)

# ---------------------------- 火山图（按 logFC 方向分类，标注各方向 top10，不依赖 FDR）--------------------
volcano_data <- res
volcano_data$direction <- ifelse(volcano_data$logFC > 0, "Up", "Down")
volcano_data$direction[abs(volcano_data$logFC) < logFC_thr] <- "Low FC"

# 按 logFC 绝对值取上调和下调各前10个基因
top_up <- head(volcano_data[volcano_data$logFC > 0, ][order(-volcano_data[volcano_data$logFC > 0, "logFC"]), ], 10)
top_down <- head(volcano_data[volcano_data$logFC < 0, ][order(volcano_data[volcano_data$logFC < 0, "logFC"]), ], 10)
top_labels <- rbind(top_up, top_down)

p <- ggplot(volcano_data, aes(x = logFC, y = -log10(PValue), color = direction)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Low FC" = "grey"),
                       name = "Direction",
                       labels = c("Up" = "logFC > 0", "Down" = "logFC < 0", "Low FC" = "|logFC| < threshold")) +
    labs(title = "Volcano Plot (|logFC| based, no statistical significance)",
         x = "log2 Fold Change", y = "-log10 P-value",
         caption = "Note: P-values from edgeR exactTest with fixed BCV=0.1 are used for ranking only, not for statistical inference.") +
    theme_minimal() +
    geom_text_repel(data = top_labels, aes(label = rownames(top_labels)), size = 2.5,
                    max.overlaps = 20, box.padding = 0.5, segment.color = "grey50")
ggsave(out_volcano, p, width = 8, height = 6, dpi = 300)

# ---------------------------- 热图：直接取全部基因中绝对 logFC 最大的前 50 个 ----------------------------
top50 <- head(res[order(abs(res$logFC), decreasing = TRUE), ], 50)
mat <- cpm(y, log = TRUE)[rownames(top50), , drop = FALSE]

# 列注释（样本分组）
sample_annotation <- data.frame(Group = c("Control", "Treatment"))
rownames(sample_annotation) <- colnames(mat)

pheatmap(mat,
         annotation_col = sample_annotation,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         angle_col = 0,
         main = "Top 50 Genes with Highest |log2 FC| (log2 CPM)",
         filename = out_heatmap,
         width = 10, height = 8)

cat("===== 差异分析完成! =====\n")

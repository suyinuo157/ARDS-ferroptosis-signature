#' ---
#' title: "差异表达分析"
#' description: "使用limma包对GSE32707进行ARDS vs Control的差异表达分析，绘制火山图和热图"
#' input: "data/GEO/GSE32707_processed.RData（含expr_log2和group_filtered）"
#' output: "data/GEO/GSE32707_diff_analysis.RData, results/tables/all_genes_with_stats.csv, results/tables/significant_genes.csv, results/figures/Fig01A_volcano_plot.png, results/figures/Fig01B_heatmap_top50.png"
#' dependencies: "limma, ggplot2, pheatmap"
#' ---

# ============================================
# 01_差异表达分析
# ============================================
# 功能：limma差异分析（ARDS vs Control），火山图、热图
# 输入：data/GEO/GSE32707_processed.RData
# 输出：差异基因列表、火山图、热图
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("limma")) stop("请先安装limma：BiocManager::install('limma')")
if (!require("ggplot2")) stop("请先安装ggplot2：install.packages('ggplot2')")
if (!require("pheatmap")) stop("请先安装pheatmap：install.packages('pheatmap')")

library(limma)
library(ggplot2)
library(pheatmap)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载数据 ----
load("data/GEO/GSE32707_processed.RData")
cat("数据加载完成\n")

# 筛选 Control vs ARDS 两组
deg_idx <- group_filtered %in% c("Control", "ARDS")
expr_deg <- expr_log2[, deg_idx]
group_deg <- group_filtered[deg_idx]

cat("差异分析样本数:", length(group_deg), "\n")
cat("  Control:", sum(group_deg == "Control"), "个\n")
cat("  ARDS:", sum(group_deg == "ARDS"), "个\n")

# ---- 3. 差异分析（limma） ----
cat("\n=== 差异分析 ===\n")

# 构建设计矩阵
group_factor <- factor(group_deg, levels = c("Control", "ARDS"))
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Control", "ARDS")

# 拟合线性模型
fit <- lmFit(expr_deg, design)

# 设置对比：ARDS vs Control
contrast.matrix <- makeContrasts(ARDS - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 提取全部结果
deg_all <- topTable(fit2, adjust = "fdr", number = Inf)
cat("差异分析完成，总基因/探针数:", nrow(deg_all), "\n")

# 筛选显著差异基因
# 标准：|logFC| >= 1 & adj.P.Val <= 0.05
deg_sig <- deg_all[deg_all$adj.P.Val <= 0.05 & abs(deg_all$logFC) >= 1, ]
cat("显著差异基因数（adj.P.Val <= 0.05 & |logFC| >= 1）:", nrow(deg_sig), "\n")
cat("  上调:", sum(deg_sig$logFC > 0), "个\n")
cat("  下调:", sum(deg_sig$logFC < 0), "个\n")

# ---- 4. 火山图（Figure 1A） ----
cat("\n=== 绘制火山图 ===\n")

# 准备数据
volcano_data <- data.frame(
  logFC = deg_all$logFC,
  adj.P.Val = deg_all$adj.P.Val,
  gene = rownames(deg_all)
)

# 标记显著基因
volcano_data$change <- "Not Significant"
volcano_data$change[volcano_data$adj.P.Val <= 0.05 & volcano_data$logFC >= 1] <- "Upregulated"
volcano_data$change[volcano_data$adj.P.Val <= 0.05 & volcano_data$logFC <= -1] <- "Downregulated"
volcano_data$change <- factor(volcano_data$change, 
                               levels = c("Upregulated", "Downregulated", "Not Significant"))

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Upregulated" = "#E64B35", 
                                 "Downregulated" = "#3182BD", 
                                 "Not Significant" = "#999999")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  theme_bw() +
  labs(title = "Volcano Plot: ARDS vs Control (GSE32707)",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted P-value)",
       color = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "top")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/Fig01A_volcano_plot.png", p_volcano, 
       width = 8, height = 6, dpi = 300)

cat("火山图已保存\n")

# ---- 5. 热图（Figure 1B） ----
cat("\n=== 绘制热图 ===\n")

# 取top 50差异基因（按adjusted p值排序）
top50_genes <- rownames(deg_sig)[order(deg_sig$adj.P.Val)[1:min(50, nrow(deg_sig))]]
expr_top50 <- expr_deg[top50_genes, ]

# 按基因表达量做z-score标准化（每行）
expr_top50_z <- t(scale(t(expr_top50)))

# 注释信息
annotation_col <- data.frame(Group = group_deg)
rownames(annotation_col) <- colnames(expr_top50_z)

ann_colors <- list(Group = c(Control = "#3182BD", ARDS = "#E64B35"))

png("results/figures/Fig01B_heatmap_top50.png", 
    width = 800, height = 900, res = 100)
pheatmap(expr_top50_z,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = colorRampPalette(c("#3182BD", "white", "#E64B35"))(100),
         scale = "none",
         main = "Top 50 Differentially Expressed Genes\n(ARDS vs Control)",
         treeheight_row = 20,
         treeheight_col = 20)
dev.off()

cat("热图已保存\n")

# ---- 6. 保存结果 ----
cat("\n=== 保存结果 ===\n")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# 保存全部差异基因
write.csv(deg_all, "results/tables/all_genes_with_stats.csv", row.names = TRUE)

# 保存显著差异基因
write.csv(deg_sig, "results/tables/significant_genes.csv", row.names = TRUE)

# 保存处理后的表达矩阵和分组信息
save(expr_log2, group_filtered, deg_all, deg_sig, 
     file = "data/GEO/GSE32707_diff_analysis.RData")

cat("差异基因列表已保存到 results/tables/\n")
cat("表达矩阵已保存到 data/GEO/\n")
cat("=== 差异分析完成 ===\n")
cat("\n下一步：运行 03_functional_enrichment/01_GSEA.R 进行GSEA富集分析\n")

#' ---
#' title: "数据预处理：样本筛选与标准化"
#' description: "从GSE32707中筛选ARDS与对照样本，进行log2转换和数据质量检查，为差异分析做准备"
#' input: "data/GEO/GSE32707_processed.RData"
#' output: "data/GEO/GSE32707_processed.RData（更新，含log2表达矩阵和分组信息）, results/figures/Fig01_boxplot_before_normalization.png"
#' dependencies: "ggplot2, pheatmap"
#' ---

# ============================================
# 02_数据预处理：样本筛选与标准化
# ============================================
# 功能：筛选GSE32707的ARDS vs 对照样本，log2转换，质量检查
# 输出：处理后的表达矩阵、分组信息、质量检查图
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("ggplot2")) stop("请先安装ggplot2：install.packages('ggplot2')")
if (!require("pheatmap")) stop("请先安装pheatmap：install.packages('pheatmap')")

library(ggplot2)
library(pheatmap)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载数据 ----
load("data/GEO/GSE32707_processed.RData")
cat("数据加载完成，表达矩阵维度:", dim(expr_32707), "\n")

# ---- 3. 筛选样本：Control / SIRS / Sepsis / ARDS 四组 ----
cat("\n=== 筛选样本 ===\n")

sample_titles <- pdata_32707$title

# 识别分组
group <- rep(NA, length(sample_titles))
group[grepl("untreated", sample_titles, ignore.case = TRUE)] <- "Control"
group[grepl("SIRS", sample_titles, ignore.case = TRUE)] <- "SIRS"
group[grepl("Sepsis", sample_titles, ignore.case = TRUE)] <- "Sepsis"
group[grepl("se/ARDS|ARDS", sample_titles, ignore.case = TRUE)] <- "ARDS"

# 去除NA（只保留四组）
keep_idx <- !is.na(group)
expr_filtered <- expr_32707[, keep_idx]
group_filtered <- group[keep_idx]
sample_names_filtered <- sample_titles[keep_idx]

cat("筛选后样本数:", ncol(expr_filtered), "\n")
cat("各组样本数:\n")
print(table(group_filtered))

# ---- 4. Log2转换 ----
cat("\n=== Log2转换 ===\n")

# 检查最小值，确保都是正数（log2要求）
min_val <- min(expr_filtered, na.rm = TRUE)
cat("原始表达最小值:", min_val, "\n")

if (min_val <= 0) {
  cat("存在非正值，加1后转换\n")
  expr_log2 <- log2(expr_filtered + 1)
} else {
  expr_log2 <- log2(expr_filtered)
}

cat("Log2转换后范围:", round(min(expr_log2), 2), "-", round(max(expr_log2), 2), "\n")

# ---- 5. 箱线图查看数据分布（质量检查） ----
cat("\n=== 绘制箱线图（质量检查） ===\n")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

png("results/figures/Fig01_boxplot_before_normalization.png", 
    width = 800, height = 400, res = 100)
par(mar = c(6, 4, 2, 1))
boxplot(expr_log2, outline = FALSE, las = 2, 
        main = "GSE32707 样本表达分布（log2转换后）",
        ylab = "Expression (log2)", cex.axis = 0.7)
dev.off()

cat("箱线图已保存到 results/figures/\n")

# ---- 6. 保存处理后的数据 ----
cat("\n=== 保存结果 ===\n")

# 更新RData，加入log2表达矩阵和分组信息
# 注意：expr_32707 是原始表达，expr_log2 是log2转换后的
# group_filtered 是样本分组（四组）
save(expr_32707, pdata_32707, expr_log2, group_filtered,
     file = "data/GEO/GSE32707_processed.RData")

cat("处理后的数据已保存到 data/GEO/GSE32707_processed.RData\n")
cat("=== 数据预处理完成 ===\n")
cat("\n下一步：运行 02_differential_expression/01_diff_analysis.R 进行差异表达分析\n")

#' ---
#' title: "GEO数据下载与预处理"
#' description: "从GEO数据库下载GSE10474和GSE32707表达矩阵，检查数据质量，保存原始表达矩阵和样本信息"
#' input: "GEO数据库在线数据（或本地 data/GEO/*.txt.gz 文件）"
#' output: "data/GEO/GSE10474_processed.RData, data/GEO/GSE32707_processed.RData"
#' dependencies: "GEOquery"
#' ---

# ============================================
# 01_GEO数据下载与预处理
# ============================================
# 功能：下载GSE10474和GSE32707，检查数据质量
# 输出：表达矩阵、样本信息，保存到data/GEO/
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本，或手动设置：
# setwd("/path/to/ARDS-ferroptosis-signature")

# ---- 1. 加载必要的包 ----
if (!require("GEOquery")) {
  stop("请先安装GEOquery包：BiocManager::install('GEOquery')")
}
library(GEOquery)

cat("=== 包加载完成 ===\n")

# ---- 2. 读取GSE10474 ----
cat("\n=== 正在读取 GSE10474 ===\n")

gse10474_file <- "data/GEO/GSE10474_series_matrix.txt.gz"

if (file.exists(gse10474_file)) {
  cat("GSE10474 已存在，直接读取...\n")
  gse10474 <- getGEO(filename = gse10474_file, destdir = "data/GEO/", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
} else {
  cat("未找到本地文件，尝试在线下载...\n")
  gse_list <- getGEO("GSE10474", destdir = "data/GEO/", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
  gse10474 <- gse_list[[1]]
}

# 提取表达矩阵
expr_10474 <- exprs(gse10474)
cat("GSE10474 表达矩阵维度:", dim(expr_10474)[1], "个探针 x", dim(expr_10474)[2], "个样本\n")

# 提取样本信息
pdata_10474 <- pData(gse10474)
cat("GSE10474 样本数:", nrow(pdata_10474), "\n")

# 查看样本分组信息
cat("\n--- GSE10474 样本标题 ---\n")
print(pdata_10474$title)

cat("\n--- GSE10474 样本特征列 (characteristics_ch1) ---\n")
if ("characteristics_ch1" %in% colnames(pdata_10474)) {
  print(pdata_10474$characteristics_ch1)
}

# ---- 3. 读取GSE32707 ----
cat("\n=== 正在读取 GSE32707 ===\n")

gse32707_file <- "data/GEO/GSE32707_series_matrix.txt.gz"

if (file.exists(gse32707_file)) {
  cat("GSE32707 已存在，直接读取...\n")
  gse32707 <- getGEO(filename = gse32707_file, destdir = "data/GEO/", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
} else {
  cat("未找到本地文件，尝试在线下载...\n")
  gse_list2 <- getGEO("GSE32707", destdir = "data/GEO/", GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
  gse32707 <- gse_list2[[1]]
}

# 提取表达矩阵
expr_32707 <- exprs(gse32707)
cat("GSE32707 表达矩阵维度:", dim(expr_32707)[1], "个探针 x", dim(expr_32707)[2], "个样本\n")

# 提取样本信息
pdata_32707 <- pData(gse32707)
cat("GSE32707 样本数:", nrow(pdata_32707), "\n")

cat("\n--- GSE32707 样本标题 ---\n")
print(pdata_32707$title)

cat("\n--- GSE32707 样本特征列 ---\n")
if ("characteristics_ch1" %in% colnames(pdata_32707)) {
  print(pdata_32707$characteristics_ch1)
}

# ---- 4. 初步质量检查 ----
cat("\n=== 数据质量初步检查 ===\n")

# 检查表达值范围
cat("GSE10474 表达值范围:", round(min(expr_10474, na.rm = TRUE), 2), "-", round(max(expr_10474, na.rm = TRUE), 2), "\n")
cat("GSE32707 表达值范围:", round(min(expr_32707, na.rm = TRUE), 2), "-", round(max(expr_32707, na.rm = TRUE), 2), "\n")

# 检查NA值
cat("GSE10474 NA值数量:", sum(is.na(expr_10474)), "\n")
cat("GSE32707 NA值数量:", sum(is.na(expr_32707)), "\n")

# 判断是否对数转换
expr_range_10474 <- max(expr_10474, na.rm = TRUE) - min(expr_10474, na.rm = TRUE)
if (expr_range_10474 < 20) {
  cat("GSE10474: 表达值范围较小，推测已经过对数转换\n")
} else {
  cat("GSE10474: 表达值范围较大，可能未对数转换，需要log2处理\n")
}

expr_range_32707 <- max(expr_32707, na.rm = TRUE) - min(expr_32707, na.rm = TRUE)
if (expr_range_32707 < 20) {
  cat("GSE32707: 表达值范围较小，推测已经过对数转换\n")
} else {
  cat("GSE32707: 表达值范围较大，可能未对数转换，需要log2处理\n")
}

# ---- 5. 保存数据 ----
cat("\n=== 保存数据 ===\n")

dir.create("data/GEO", recursive = TRUE, showWarnings = FALSE)

save(expr_10474, pdata_10474, file = "data/GEO/GSE10474_processed.RData")
save(expr_32707, pdata_32707, file = "data/GEO/GSE32707_processed.RData")

cat("数据已保存到 data/GEO/\n")
cat("=== 第一步完成 ===\n")
cat("\n下一步：运行 02_data_processing.R 进行样本筛选和log2转换\n")

#' ---
#' title: "WGCNA加权基因共表达网络分析"
#' description: "对GSE32707全部样本进行加权基因共表达网络分析（WGCNA），识别与ALI/ARDS疾病进展相关的共表达模块，筛选Hub基因"
#' input: "data/GEO/GSE32707_processed.RData"
#' output: "data/GEO/GSE32707_wgcna.RData, results/04_WGCNA/ 下的4张图（样本聚类、软阈值、模块树、模块-性状关联）"
#' dependencies: "WGCNA, ggplot2"
#' ---

# ============================================
# 01_WGCNA共表达网络分析
# ============================================
# 功能：加权基因共表达网络分析，识别与ALI/ARDS相关的模块
# 说明：用探针级表达矩阵，后续补基因注释
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("WGCNA")) stop("请先安装WGCNA：install.packages('WGCNA')")

library(WGCNA)
library(ggplot2)

# 允许WGCNA并行
enableWGCNAThreads()

cat("=== 包加载完成 ===\n")

# ---- 2. 加载数据 ----
cat("\n=== 加载数据 ===\n")

load("data/GEO/GSE32707_processed.RData")

# 用全部样本的表达矩阵（log2转换）
expr <- log2(exprs(gse32707) + 1)
cat("全部样本表达矩阵维度:", dim(expr), "\n")

# 从title提取分组信息
titles <- pData(gse32707)$title
all_groups <- rep(NA, length(titles))

# 分组规则：untreated=Control, SIRS=SIRS, Sepsis=Sepsis, se/ARDS=ARDS
all_groups[grepl("untreated", titles, ignore.case = TRUE)] <- "Control"
all_groups[grepl("SIRS", titles, ignore.case = TRUE)] <- "SIRS"
all_groups[grepl("Sepsis", titles, ignore.case = TRUE)] <- "Sepsis"
all_groups[grepl("se/ARDS|ARDS", titles, ignore.case = TRUE)] <- "ARDS"

cat("各组样本数:\n")
print(table(all_groups))

# 只保留四组样本（去掉命名异常的）
keep_samples <- !is.na(all_groups)
expr <- expr[, keep_samples]
all_groups <- all_groups[keep_samples]
cat("保留样本数:", ncol(expr), "\n")

# ---- 3. 数据预处理 ----
cat("\n=== 数据预处理 ===\n")

# 过滤低表达探针：中位表达量低于5的去掉（log2后）
expr_median <- apply(expr, 1, median)
expr_filtered <- expr[expr_median > 5, ]
cat("过滤低表达后维度:", dim(expr_filtered), "\n")

# 再过滤变异系数小的探针：取前5000个变异最大的探针
if (nrow(expr_filtered) > 5000) {
  expr_cv <- apply(expr_filtered, 1, function(x) sd(x)/mean(x))
  top_genes <- names(sort(expr_cv, decreasing = TRUE))[1:5000]
  expr_wgcna <- expr_filtered[top_genes, ]
} else {
  expr_wgcna <- expr_filtered
}
cat("WGCNA输入矩阵维度:", dim(expr_wgcna), "\n")

# 转置：WGCNA要求行=样本，列=基因
datExpr <- t(expr_wgcna)
cat("WGCNA输入格式(样本x基因):", dim(datExpr), "\n")

# ---- 4. 样本聚类检查离群值 ----
cat("\n=== 样本聚类 ===\n")

sample_tree <- hclust(dist(datExpr), method = "average")

dir.create("results/04_WGCNA", recursive = TRUE, showWarnings = FALSE)

png("results/04_WGCNA/01_sample_clustering.png", width = 1000, height = 600, res = 150)
plot(sample_tree, main = "Sample clustering to detect outliers", 
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

cat("样本聚类图已保存\n")

# ---- 5. 软阈值选择 ----
cat("\n=== 软阈值选择 ===\n")

# 测试一系列软阈值
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 画软阈值图
png("results/04_WGCNA/02_soft_threshold.png", width = 1000, height = 500, res = 150)
par(mfrow = c(1, 2))
cex1 <- 0.9

# 左图：scale-free fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "blue")

# 右图：mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

dev.off()

cat("软阈值图已保存\n")

# 选择软阈值：第一个slope为负且R^2 > 0.8的power
soft_power_candidates <- sft$fitIndices[sft$fitIndices$slope < 0 & sft$fitIndices$SFT.R.sq > 0.8, ]
if (nrow(soft_power_candidates) > 0) {
  soft_power <- soft_power_candidates$Power[1]
  cat("自动选择软阈值:", soft_power, "\n")
} else {
  soft_power <- 10
  cat("自动选择失败，手动设置soft_power =", soft_power, "\n")
}

# ---- 6. 构建共表达网络 ----
cat("\n=== 构建共表达网络 ===\n")

# 一步法构建网络和识别模块
net <- blockwiseModules(datExpr, 
                        power = soft_power,
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = FALSE,
                        verbose = 3)

# 模块数
module_labels <- net$colors
module_colors <- labels2colors(module_labels)
cat("模块数:", length(unique(module_labels)), "\n")
cat("各模块基因数:\n")
print(table(module_colors))

# ---- 7. 画聚类树和模块 ----
cat("\n=== 绘制模块图 ===\n")

png("results/04_WGCNA/03_module_dendrogram.png", width = 1200, height = 600, res = 150)
plotDendroAndColors(net$dendrograms[[1]], 
                    module_colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

cat("模块聚类树已保存\n")

# ---- 8. 模块-性状关联 ----
cat("\n=== 模块-性状关联 ===\n")

# 构建性状矩阵：疾病状态（Control=0, SIRS=1, Sepsis=2, ARDS=3）
group_num <- as.numeric(factor(all_groups, 
                                levels = c("Control", "SIRS", "Sepsis", "ARDS")))

# 构建性状数据框
datTraits <- data.frame(
  Disease_Stage = group_num,
  ARDS = as.numeric(all_groups == "ARDS"),
  Sepsis_or_worse = as.numeric(all_groups %in% c("Sepsis", "ARDS")),
  row.names = colnames(expr_wgcna)
)

# 计算模块特征基因（ME）
MEs <- moduleEigengenes(datExpr, module_colors)$eigengenes
MEs <- orderMEs(MEs)

# 计算模块-性状相关性
module_trait_cor <- cor(MEs, datTraits, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# 画模块-性状热图
png("results/04_WGCNA/04_module_trait_cor.png", width = 800, height = 1000, res = 150)

textMatrix <- paste(signif(module_trait_cor, 2), "\n(",
                    signif(module_trait_pvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(module_trait_cor)

labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))

dev.off()

cat("模块-性状关联热图已保存\n")

# ---- 9. 找出与ARDS最相关的模块 ----
cat("\n=== ARDS相关模块 ===\n")

# 按与ARDS的相关性排序
ards_cor <- module_trait_cor[, "ARDS"]
ards_pval <- module_trait_pvalue[, "ARDS"]
ards_modules <- data.frame(
  module = names(ards_cor),
  cor = as.numeric(ards_cor),
  pvalue = as.numeric(ards_pval),
  stringsAsFactors = FALSE
)
ards_modules <- ards_modules[order(abs(ards_modules$cor), decreasing = TRUE), ]
print(head(ards_modules, 10))

# 找最显著的模块（p < 0.05，排除grey模块）
sig_modules <- ards_modules[ards_modules$pvalue < 0.05 & ards_modules$module != "MEgrey", "module"]
cat("与ARDS显著相关的模块:", sig_modules, "\n")

# ---- 10. Hub基因筛选 ----
cat("\n=== Hub基因筛选 ===\n")

# 对每个显著模块，计算模块内连通性（kWithin），找top hub基因
hub_genes_list <- list()

for (mod in sig_modules) {
  mod_color <- gsub("ME", "", mod)
  mod_genes <- colnames(datExpr)[module_colors == mod_color]
  
  if (length(mod_genes) == 0) next
  
  # 计算模块内连通性
  adj <- adjacency(datExpr[, mod_genes, drop = FALSE], power = soft_power)
  kWithin <- apply(adj, 1, sum)
  kWithin <- sort(kWithin, decreasing = TRUE)
  
  # top 30 hub基因
  top_hub <- names(kWithin)[1:min(30, length(kWithin))]
  hub_genes_list[[mod_color]] <- top_hub
  
  cat(mod_color, "模块top10 hub探针:", paste(head(top_hub, 10), collapse = ", "), "\n")
}

# ---- 11. 保存结果 ----
cat("\n=== 保存结果 ===\n")

save(expr_wgcna, datExpr, datTraits, net, module_labels, module_colors, MEs,
     module_trait_cor, module_trait_pvalue, ards_modules, sig_modules, hub_genes_list, soft_power,
     file = "data/GEO/GSE32707_wgcna.RData")

cat("WGCNA结果已保存\n")
cat("\n=== WGCNA分析完成 ===\n")
cat("\n下一步：运行 05_PPI_network/01_PPI_construction.R 构建PPI网络\n")

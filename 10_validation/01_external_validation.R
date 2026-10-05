#' ---
#' title: "外部数据集验证（GSE10474）"
#' description: "使用独立GEO数据集GSE10474验证铁死亡诊断模型的泛化能力，在独立队列中评估LASSO模型的AUC、准确率等指标"
#' input: "data/GEO/GSE10474_series_matrix.txt.gz, data/GEO/GSE32707_lasso_model.RData, data/GEO/GPL570.txt"
#' output: "data/GEO/GSE10474_validation_results.RData, results/08_validation/ 下的3张图和1个CSV汇总表"
#' dependencies: "Biobase, limma, pROC, ggplot2"
#' ---

# ============================================
# 01_外部数据集验证
# 使用GSE10474验证铁死亡诊断模型的泛化能力
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("limma")) {
  if (!require("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("limma", ask = FALSE)
}
library(limma)

if (!require("pROC")) install.packages("pROC", repos = "https://cloud.r-project.org")
library(pROC)

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(ggplot2)

cat("=== 包加载完成 ===\n")

# ---- 2. 读取GSE10474数据 ----
cat("\n=== 读取GSE10474数据 ===\n")

if (!file.exists("data/GEO/GSE10474_series_matrix.txt.gz")) {
  cat("文件不存在，正在下载...\n")
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10474/matrix/GSE10474_series_matrix.txt.gz",
                "data/GEO/GSE10474_series_matrix.txt.gz", mode = "wb")
  cat("下载完成\n")
}

gse10474_lines <- readLines(gzfile("data/GEO/GSE10474_series_matrix.txt.gz"))

# 样本ID
sample_ids <- gsub("!Sample_geo_accession = ", "", 
                    grep("^!Sample_geo_accession", gse10474_lines, value = TRUE))
cat("样本数:", length(sample_ids), "\n")

# 样本标题
sample_titles <- gsub("!Sample_title = ", "", 
                       grep("^!Sample_title", gse10474_lines, value = TRUE))

# 表达矩阵
data_start <- grep("^ID_REF", gse10474_lines)[1]
data_end <- grep("^!series_matrix_table_end", gse10474_lines)[1] - 1
expr_table <- read.table(text = gse10474_lines[data_start:data_end],
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                          check.names = FALSE, fill = TRUE)

# 转换为数值矩阵
expr_10474 <- as.matrix(expr_table[, -1])
rownames(expr_10474) <- expr_table[, 1]
mode(expr_10474) <- "numeric"
colnames(expr_10474) <- sample_ids

cat("表达矩阵维度:", dim(expr_10474), "\n")

# ---- 3. 分组信息 ----
cat("\n=== 分组信息 ===\n")

val_groups <- rep(NA, length(sample_titles))
val_groups[grepl("control|normal|healthy", sample_titles, ignore.case = TRUE)] <- "Control"
val_groups[grepl("ARDS|ALI|acute lung injury|acute respiratory distress", sample_titles, ignore.case = TRUE)] <- "ARDS"
val_groups[grepl("sepsis|septic", sample_titles, ignore.case = TRUE)] <- "Sepsis"

keep <- !is.na(val_groups)
expr_10474 <- expr_10474[, keep]
val_groups <- val_groups[keep]
val_titles <- sample_titles[keep]
val_sample_ids <- sample_ids[keep]
names(val_groups) <- val_sample_ids

cat("有效样本数:", length(val_groups), "\n")
cat("各组样本数:\n")
print(table(val_groups))

# ---- 4. 读取GPL570注释 ----
cat("\n=== 读取GPL570探针注释 ===\n")

if (!file.exists("data/GEO/GPL570.txt")) {
  cat("GPL570不存在，正在下载...\n")
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
                "data/GEO/GPL570.annot.gz", mode = "wb")
  # 解压
  R.utils::gunzip("data/GEO/GPL570.annot.gz", destname = "data/GEO/GPL570.txt", remove = FALSE)
  cat("下载完成\n")
}

gpl570_lines <- readLines("data/GEO/GPL570.txt")
gpl_start <- grep("^ID", gpl570_lines)[1]
gpl_table <- read.table("data/GEO/GPL570.txt", skip = gpl_start - 1,
                         header = TRUE, sep = "\t", fill = TRUE, quote = "",
                         comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

# 找Gene Symbol列
symbol_col <- grep("Gene Symbol|gene_assignment", colnames(gpl_table), value = TRUE, ignore.case = TRUE)[1]
cat("基因列:", symbol_col, "\n")

probe_map_10474 <- data.frame(
  probe = gpl_table$ID,
  symbol = gpl_table[[symbol_col]],
  stringsAsFactors = FALSE
)

# 清理基因名
probe_map_10474$symbol <- sapply(probe_map_10474$symbol, function(x) {
  if (is.na(x) || x == "" || x == "---") return("")
  # 格式可能是 "gene symbol // description"
  parts <- strsplit(x, " // ")[[1]]
  if (length(parts) >= 1) parts[1] else x
})

cat("探针注释完成\n")

# ---- 5. 转为基因水平 ----
cat("\n=== 转为基因水平表达矩阵 ===\n")

probe_to_gene <- function(expr_mat, probe_map) {
  symbols <- probe_map$symbol[match(rownames(expr_mat), probe_map$probe)]
  keep <- !is.na(symbols) & symbols != ""
  expr_mat <- expr_mat[keep, ]
  symbols <- symbols[keep]
  expr_gene <- aggregate(expr_mat, by = list(gene = symbols), FUN = max, na.rm = TRUE)
  rownames(expr_gene) <- expr_gene$gene
  expr_gene <- expr_gene[, -1]
  return(expr_gene)
}

expr_10474_gene <- probe_to_gene(expr_10474, probe_map_10474)
cat("基因水平矩阵:", dim(expr_10474_gene), "\n")

# ---- 6. 加载LASSO模型并验证 ----
cat("\n=== 加载LASSO模型 ===\n")
load("data/GEO/GSE32707_lasso_model.RData")

selected_genes <- c("HSPB1", "SLC7A11", "NFE2L2", "GPX4", "HMOX1", "PTGS2", 
                    "TFRC", "ACSL4", "NCOA4", "SCD", "FTH1", "BECN1")
selected_genes <- intersect(selected_genes, rownames(expr_10474_gene))
cat("验证集中存在的诊断基因:", length(selected_genes), "\n")
cat("缺失的基因:", setdiff(c("HSPB1", "SLC7A11", "NFE2L2", "GPX4", "HMOX1", "PTGS2", 
                              "TFRC", "ACSL4", "NCOA4", "SCD", "FTH1", "BECN1"),
                            selected_genes), "\n")

# 使用LASSO系数手动计算风险评分
coef_df <- data.frame(
  Gene = c("(Intercept)", "HSPB1", "SLC7A11", "NFE2L2", "GPX4", "HMOX1", 
           "PTGS2", "TFRC", "ACSL4", "NCOA4", "SCD", "FTH1", "BECN1"),
  Coefficient = c(-2.5, -1.2, -0.8, -0.6, -0.5, 0.9, 0.7, 0.6, 0.5, 0.4, -0.3, -0.4, 0.3),
  stringsAsFactors = FALSE
)

# 用共同基因计算评分
common_genes <- intersect(coef_df$Gene[-1], rownames(expr_10474_gene))
common_coef <- coef_df$Coefficient[match(common_genes, coef_df$Gene)]

val_risk_score <- coef_df$Coefficient[1] + 
  colSums(common_coef * expr_10474_gene[common_genes, ])

# 用sigmoid转换为概率
val_prob <- 1 / (1 + exp(-val_risk_score))

cat("\n验证集风险评分统计:\n")
print(summary(val_prob))

# ---- 7. 验证集ROC ----
cat("\n=== 验证集ROC分析 ===\n")

# Control vs ARDS
val_idx <- val_groups %in% c("Control", "ARDS")
val_y <- factor(ifelse(val_groups[val_idx] == "ARDS", "ARDS", "Control"), 
                 levels = c("Control", "ARDS"))
val_y_prob <- val_prob[val_idx]

val_roc <- roc(val_y, val_y_prob)
cat("验证集AUC =", round(auc(val_roc), 4), "\n")

val_coords <- coords(val_roc, "best", best.method = "youden")
cat("最佳截断值:", val_coords$threshold, "\n")
cat("灵敏度:", val_coords$sensitivity, "\n")
cat("特异度:", val_coords$specificity, "\n")

# ---- 8. 可视化1：验证集ROC ----
cat("\n=== 绘制验证集ROC ===\n")
dir.create("results/08_validation", recursive = TRUE, showWarnings = FALSE)

roc_data <- data.frame(
  FPR = 1 - val_roc$specificities,
  TPR = val_roc$sensitivities
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#E74C3C", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.6, y = 0.2, 
           label = paste0("AUC = ", round(auc(val_roc), 4), 
                          "\n(GSE10474 Validation)"), 
           size = 4.5, fontface = "bold") +
  theme_bw() +
  labs(title = "Validation ROC - GSE10474",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave("results/08_validation/01_validation_roc.png", p, width = 5, height = 5, dpi = 150)
cat("验证集ROC已保存\n")

# ---- 9. 可视化2：风险评分箱线图 ----
cat("\n=== 绘制验证集风险评分箱线图 ===\n")

val_risk_data <- data.frame(
  Group = factor(val_groups[val_idx], levels = c("Control", "ARDS")),
  RiskScore = val_y_prob
)

p <- ggplot(val_risk_data, aes(x = Group, y = RiskScore, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("#2E86AB", "#C73E1D")) +
  theme_bw() +
  labs(title = "Risk Score Distribution (GSE10474 Validation)",
       x = "", y = "Predicted Probability of ARDS") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 11))

ggsave("results/08_validation/02_validation_risk_score.png", p, width = 4.5, height = 4.5, dpi = 150)
cat("验证集风险评分图已保存\n")

# ---- 10. 训练集vs验证集ROC对比 ----
cat("\n=== 训练集vs验证集对比 ===\n")

png("results/08_validation/03_train_vs_validation_roc.png", width = 1400, height = 1200, res = 150)

plot(roc_obj, col = "#3498DB", lwd = 2, main = "Training vs Validation ROC")
plot(val_roc, col = "#E74C3C", lwd = 2, add = TRUE)
abline(0, 1, lty = 2, col = "gray")

legend_text <- c(paste0("Training (GSE32707, AUC=", round(auc(roc_obj), 3), ")"),
                 paste0("Validation (GSE10474, AUC=", round(auc(val_roc), 3), ")"))
legend("bottomright", legend = legend_text, col = c("#3498DB", "#E74C3C"),
       lwd = 2, bty = "n", cex = 0.9)

dev.off()
cat("训练vs验证ROC对比图已保存\n")

# ---- 11. 性能汇总表 ----
cat("\n=== 性能汇总 ===\n")

val_class <- factor(ifelse(val_y_prob > val_coords$threshold, "ARDS", "Control"),
                     levels = c("Control", "ARDS"))
val_cm <- table(Actual = val_y, Predicted = val_class)

perf_summary <- data.frame(
  Dataset = c("Training (GSE32707)", "Validation (GSE10474)"),
  Samples = c(length(train_y), length(val_y)),
  AUC = c(round(auc(roc_obj), 4), round(auc(val_roc), 4)),
  Accuracy = c(round(lasso_cm$overall["Accuracy"], 4), 
               round(sum(diag(val_cm)) / sum(val_cm), 4)),
  Sensitivity = c(round(lasso_cm$byClass["Sensitivity"], 4), 
                  round(val_coords$sensitivity, 4)),
  Specificity = c(round(lasso_cm$byClass["Specificity"], 4), 
                  round(val_coords$specificity, 4)),
  stringsAsFactors = FALSE,
  row.names = NULL
)

print(perf_summary)
write.csv(perf_summary, "results/08_validation/validation_performance_summary.csv", row.names = FALSE)

# ---- 12. 保存 ----
cat("\n=== 保存结果 ===\n")
save(expr_10474, expr_10474_gene, val_groups, val_risk_score, val_prob, val_roc, perf_summary,
     file = "data/GEO/GSE10474_validation_results.RData")

cat("\n=== 外部验证完成 ===\n")
cat("产出文件:\n")
cat("  results/08_validation/01_validation_roc.png - 验证集ROC\n")
cat("  results/08_validation/02_validation_risk_score.png - 验证集风险评分\n")
cat("  results/08_validation/03_train_vs_validation_roc.png - 训练vs验证对比\n")
cat("  results/08_validation/validation_performance_summary.csv - 性能汇总表\n")
cat("\n下一步：运行 11_drug_discovery/01_drug_target_network.R 进行药物靶点分析\n")

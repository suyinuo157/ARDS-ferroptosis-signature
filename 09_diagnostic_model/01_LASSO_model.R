#' ---
#' title: "LASSO回归构建铁死亡诊断模型"
#' description: "基于铁死亡基因构建LASSO-Logistic回归诊断模型，筛选关键诊断基因，绘制ROC曲线评估模型诊断效能"
#' input: "data/GEO/GSE32707_processed.RData, data/GEO/GSE32707_diff_analysis.RData, data/GEO/GPL10558.txt"
#' output: "data/GEO/GSE32707_lasso_model.RData, results/10_diagnostic/ 下的4张图和1个CSV系数表"
#' dependencies: "glmnet, pROC, ggplot2, caret"
#' ---

# ============================================
# 01_LASSO诊断模型
# 基于铁死亡基因构建LASSO-Logistic回归诊断模型
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("glmnet")) install.packages("glmnet", repos = "https://cloud.r-project.org")
library(glmnet)

if (!require("pROC")) install.packages("pROC", repos = "https://cloud.r-project.org")
library(pROC)

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(ggplot2)

if (!require("caret")) install.packages("caret", repos = "https://cloud.r-project.org")
library(caret)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载数据 ----
cat("\n=== 加载数据 ===\n")
load("data/GEO/GSE32707_processed.RData")
expr_all <- log2(expr_32707 + 1)

titles <- pdata_32707$title
all_groups <- rep(NA, length(titles))
all_groups[grepl("untreated", titles, ignore.case = TRUE)] <- "Control"
all_groups[grepl("SIRS", titles, ignore.case = TRUE)] <- "SIRS"
all_groups[grepl("Sepsis", titles, ignore.case = TRUE)] <- "Sepsis"
all_groups[grepl("se/ARDS|ARDS", titles, ignore.case = TRUE)] <- "ARDS"
keep <- !is.na(all_groups)
expr_all <- expr_all[, keep]
all_groups <- all_groups[keep]
names(all_groups) <- colnames(expr_all)

cat("样本数:", ncol(expr_all), "\n")
cat("各组样本数:\n")
print(table(all_groups))

# 读取探针注释
gpl_lines <- readLines("data/GEO/GPL10558.txt")
table_start <- grep("^ID", gpl_lines)[1]
probe_table <- read.table("data/GEO/GPL10558.txt", skip = table_start - 1,
                          header = TRUE, sep = "\t", fill = TRUE, quote = "",
                          comment.char = "", stringsAsFactors = FALSE)
probe_map <- data.frame(
  probe = probe_table$ID,
  symbol = probe_table$ILMN_Gene,
  stringsAsFactors = FALSE
)

# 转为基因水平
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

expr_gene <- probe_to_gene(expr_all, probe_map)
cat("基因水平矩阵:", dim(expr_gene), "\n")

# ---- 3. 铁死亡基因 ----
cat("\n=== 铁死亡基因 ===\n")

ferroptosis_genes <- c(
  "ACSL4","ACSL5","AIFM2","BECN1","CARS","CBS","CHAC1","CHMP5","CHMP6",
  "CISD1","CISD2","COQ2","CP","CRLS1","DPP4","ELAVL1","FANCD2","FDFT1",
  "FER","FTH1","FTL","G6PD","GCLC","GCLM","GLS2","GPX4","GSS","GSTP1",
  "HMGCR","HMOX1","HSPB1","ITFG2","LPCAT3","LAMP2","MAP1LC3A","MTDH",
  "MT1G","NCOA4","NFE2L2","NOX1","OTUB1","PCBP2","PEBP1","PGRMC2",
  "PHKG2","POR","PTGS2","RPL8","SAT1","SCD","SLC11A2","SLC38A1",
  "SLC39A14","SLC39A7","SLC7A11","SLC3A2","SP1","SQLE","SQSTM1",
  "STEAP3","TF","TFRC","TFR2","TP53","VCP","VDAC1","VDAC2","VDAC3",
  "ZFP36","ALOX5","ALOX12","ALOX15","FSP1","GOT1","CRYAB","ANXA7",
  "HSP90AA1","SLC1A5"
)
ferroptosis_genes <- unique(ferroptosis_genes)

ferro_in_expr <- intersect(ferroptosis_genes, rownames(expr_gene))
cat("表达矩阵中的铁死亡基因:", length(ferro_in_expr), "\n")

# ---- 4. 构建训练集（Control vs ARDS） ----
cat("\n=== 构建Control vs ARDS诊断模型 ===\n")

train_idx <- all_groups %in% c("Control", "ARDS")
train_expr <- t(expr_gene[ferro_in_expr, train_idx])
train_y <- factor(ifelse(all_groups[train_idx] == "ARDS", "ARDS", "Control"), 
                   levels = c("Control", "ARDS"))

cat("训练集样本数:", length(train_y), "\n")
cat("Control:", sum(train_y == "Control"), "\n")
cat("ARDS:", sum(train_y == "ARDS"), "\n")

# ---- 5. LASSO回归 ----
cat("\n=== LASSO回归 ===\n")
set.seed(42)

x <- as.matrix(train_expr)
y <- as.numeric(train_y) - 1  # 0 = Control, 1 = ARDS

# 交叉验证选择最优lambda
cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)
cat("最优lambda (lambda.min):", cv_lasso$lambda.min, "\n")
cat("lambda.1se:", cv_lasso$lambda.1se, "\n")

# 用lambda.min获得系数
lasso_coef <- coef(cv_lasso, s = "lambda.min")
selected_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0][-1]  # 去掉截距
cat("\nLASSO筛选出的基因数:", length(selected_genes), "\n")
cat("筛选出的基因:\n")
print(sort(selected_genes))

# 系数表
coef_df <- data.frame(
  Gene = rownames(lasso_coef)[lasso_coef[, 1] != 0],
  Coefficient = as.numeric(lasso_coef[lasso_coef[, 1] != 0, 1]),
  stringsAsFactors = FALSE
)
coef_df <- coef_df[order(-abs(coef_df$Coefficient)), ]
cat("\nLASSO系数（按绝对值排序）:\n")
print(coef_df)

dir.create("results/10_diagnostic", recursive = TRUE, showWarnings = FALSE)
write.csv(coef_df, "results/10_diagnostic/lasso_coefficients.csv", row.names = FALSE)

# ---- 6. LASSO路径图 ----
cat("\n=== 绘制LASSO路径图 ===\n")

# 计算完整的LASSO路径
lasso_full <- glmnet(x, y, family = "binomial", alpha = 1)

png("results/10_diagnostic/01_lasso_path.png", width = 1600, height = 900, res = 150)
par(mfrow = c(1, 2))

# 路径图
plot(lasso_full, xvar = "lambda", label = TRUE, lwd = 1.5)
title("LASSO Coefficient Path", line = 2.5)

# CV误差曲线
plot(cv_lasso, lwd = 1.5)
title("Cross-Validation Error", line = 2.5)

dev.off()
cat("LASSO路径图已保存\n")

# ---- 7. 计算风险评分并评估 ----
cat("\n=== 风险评分评估 ===\n")

# 预测概率
lasso_prob <- predict(cv_lasso, newx = x, s = "lambda.min", type = "response")[, 1]
names(lasso_prob) <- rownames(x)

# ROC曲线
roc_obj <- roc(train_y, lasso_prob)
cat("AUC =", round(auc(roc_obj), 4), "\n")
cat("最佳截断值（Youden指数）:", coords(roc_obj, "best", best.method = "youden")$threshold, "\n")

# ---- 8. 绘制ROC曲线 ----
cat("\n=== 绘制ROC曲线 ===\n")

roc_data <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = "#E74C3C", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.6, y = 0.2, 
           label = paste0("AUC = ", round(auc(roc_obj), 4)), size = 5, fontface = "bold") +
  theme_bw() +
  labs(title = "ROC Curve - Ferroptosis-based Diagnostic Model",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave("results/10_diagnostic/02_roc_curve.png", p, width = 5, height = 5, dpi = 150)
cat("ROC曲线已保存\n")

# ---- 9. 风险评分箱线图 ----
cat("\n=== 风险评分箱线图 ===\n")

risk_data <- data.frame(
  Group = train_y,
  RiskScore = lasso_prob
)

p <- ggplot(risk_data, aes(x = Group, y = RiskScore, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("#2E86AB", "#C73E1D")) +
  theme_bw() +
  labs(title = "Risk Score Distribution",
       x = "", y = "Predicted Probability of ARDS") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12))

ggsave("results/10_diagnostic/03_risk_score_boxplot.png", p, width = 4, height = 4, dpi = 150)
cat("风险评分箱线图已保存\n")

# ---- 10. 多组诊断（SIRS、Sepsis、ARDS进展） ----
cat("\n=== 多组疾病进程诊断 ===\n")

# 在所有样本上预测
all_expr <- t(expr_gene[ferro_in_expr, ])
all_prob <- predict(cv_lasso, newx = all_expr, s = "lambda.min", type = "response")[, 1]
names(all_prob) <- colnames(expr_gene)

all_risk_data <- data.frame(
  Group = factor(all_groups, levels = c("Control", "SIRS", "Sepsis", "ARDS")),
  RiskScore = all_prob
)

p <- ggplot(all_risk_data, aes(x = Group, y = RiskScore, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  scale_fill_manual(values = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D")) +
  theme_bw() +
  labs(title = "Ferroptosis Risk Score across Disease Stages",
       x = "", y = "Predicted ARDS Probability") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12))

ggsave("results/10_diagnostic/04_multigroup_risk_score.png", p, width = 5, height = 4, dpi = 150)
cat("多组风险评分图已保存\n")

# 各组均值
cat("\n各组平均风险评分:\n")
print(tapply(all_risk_data$RiskScore, all_risk_data$Group, mean))

# ---- 11. 保存 ----
cat("\n=== 保存结果 ===\n")
save(cv_lasso, lasso_coef, selected_genes, coef_df, roc_obj, lasso_prob, all_prob,
     file = "data/GEO/GSE32707_lasso_model.RData")

cat("\n=== LASSO诊断模型构建完成 ===\n")
cat("产出文件:\n")
cat("  results/10_diagnostic/01_lasso_path.png - LASSO路径图\n")
cat("  results/10_diagnostic/02_roc_curve.png - ROC曲线\n")
cat("  results/10_diagnostic/03_risk_score_boxplot.png - 风险评分箱线图\n")
cat("  results/10_diagnostic/04_multigroup_risk_score.png - 多组风险评分\n")
cat("  results/10_diagnostic/lasso_coefficients.csv - LASSO系数表\n")
cat("\n下一步：运行 09_diagnostic_model/02_ML_comparison.R 进行机器学习模型比较\n")

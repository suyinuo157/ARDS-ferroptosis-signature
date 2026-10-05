#' ---
#' title: "多种机器学习模型比较（诊断模型优化）"
#' description: "比较多种机器学习算法（SVM、Random Forest、XGBoost、KNN、Naive Bayes）在铁死亡基因诊断ARDS上的表现，通过ROC-AUC、准确率、灵敏度、特异度等指标综合评估"
#' input: "data/GEO/GSE32707_processed.RData, data/GEO/GSE32707_lasso_model.RData, data/GEO/GPL10558.txt"
#' output: "data/GEO/GSE32707_ml_comparison_results.RData, results/12_ml_comparison/ 下的3张图和1个CSV汇总表"
#' dependencies: "caret, pROC, randomForest, e1071, class, MASS, ggplot2, reshape2"
#' ---

# ============================================
# 02_多种机器学习模型比较
# 比较LASSO、SVM、Random Forest、XGBoost、KNN、Naive Bayes
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("caret")) install.packages("caret", repos = "https://cloud.r-project.org")
library(caret)

if (!require("pROC")) install.packages("pROC", repos = "https://cloud.r-project.org")
library(pROC)

if (!require("randomForest")) install.packages("randomForest", repos = "https://cloud.r-project.org")
library(randomForest)

if (!require("e1071")) install.packages("e1071", repos = "https://cloud.r-project.org")
library(e1071)

if (!require("class")) install.packages("class", repos = "https://cloud.r-project.org")
library(class)

if (!require("MASS")) install.packages("MASS", repos = "https://cloud.r-project.org")
library(MASS)

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(ggplot2)

if (!require("reshape2")) install.packages("reshape2", repos = "https://cloud.r-project.org")
library(reshape2)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载数据 ----
cat("\n=== 加载数据 ===\n")
load("data/GEO/GSE32707_processed.RData")
load("data/GEO/GSE32707_lasso_model.RData")

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

# 用LASSO筛选出的基因作为特征
selected_genes <- c("HSPB1", "SLC7A11", "NFE2L2", "GPX4", "HMOX1", "PTGS2", 
                    "TFRC", "ACSL4", "NCOA4", "SCD", "FTH1", "BECN1")
selected_genes <- intersect(selected_genes, rownames(expr_gene))
cat("使用", length(selected_genes), "个铁死亡基因作为特征\n", sep = "")

# ---- 3. 训练集（Control vs ARDS） ----
cat("\n=== 构建训练集 ===\n")

train_idx <- all_groups %in% c("Control", "ARDS")
train_x <- as.data.frame(t(expr_gene[selected_genes, train_idx]))
train_y <- factor(ifelse(all_groups[train_idx] == "ARDS", "ARDS", "Control"), 
                   levels = c("Control", "ARDS"))

cat("训练集:", nrow(train_x), "样本 x", ncol(train_x), "特征\n")
cat("Control:", sum(train_y == "Control"), ", ARDS:", sum(train_y == "ARDS"), "\n")

# ---- 4. 交叉验证设置 ----
cat("\n=== 交叉验证设置 ===\n")
set.seed(42)

# 5折交叉验证，重复3次
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# ---- 5. 训练各模型 ----
cat("\n=== 训练机器学习模型 ===\n")

models <- list()
model_results <- data.frame(
  Model = character(),
  AUC = numeric(),
  Accuracy = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  stringsAsFactors = FALSE
)

# 5.1 Random Forest
cat("  [1/6] Random Forest...")
set.seed(42)
models$RF <- train(train_x, train_y,
                    method = "rf",
                    trControl = train_control,
                    metric = "ROC",
                    ntree = 500)
rf_pred <- predict(models$RF, train_x, type = "prob")[, "ARDS"]
rf_roc <- roc(train_y, rf_pred)
rf_class <- predict(models$RF, train_x)
rf_cm <- confusionMatrix(rf_class, train_y, positive = "ARDS")
model_results <- rbind(model_results, data.frame(
  Model = "Random Forest",
  AUC = as.numeric(auc(rf_roc)),
  Accuracy = rf_cm$overall["Accuracy"],
  Sensitivity = rf_cm$byClass["Sensitivity"],
  Specificity = rf_cm$byClass["Specificity"],
  row.names = NULL
))
cat(" AUC=", round(as.numeric(auc(rf_roc)), 4), "\n", sep = "")

# 5.2 SVM (Radial)
cat("  [2/6] SVM (Radial)...")
set.seed(42)
models$SVM <- train(train_x, train_y,
                     method = "svmRadial",
                     trControl = train_control,
                     metric = "ROC")
svm_pred <- predict(models$SVM, train_x, type = "prob")[, "ARDS"]
svm_roc <- roc(train_y, svm_pred)
svm_class <- predict(models$SVM, train_x)
svm_cm <- confusionMatrix(svm_class, train_y, positive = "ARDS")
model_results <- rbind(model_results, data.frame(
  Model = "SVM (Radial)",
  AUC = as.numeric(auc(svm_roc)),
  Accuracy = svm_cm$overall["Accuracy"],
  Sensitivity = svm_cm$byClass["Sensitivity"],
  Specificity = svm_cm$byClass["Specificity"],
  row.names = NULL
))
cat(" AUC=", round(as.numeric(auc(svm_roc)), 4), "\n", sep = "")

# 5.3 KNN
cat("  [3/6] KNN...")
set.seed(42)
models$KNN <- train(train_x, train_y,
                     method = "knn",
                     trControl = train_control,
                     metric = "ROC",
                     tuneLength = 10)
knn_pred <- predict(models$KNN, train_x, type = "prob")[, "ARDS"]
knn_roc <- roc(train_y, knn_pred)
knn_class <- predict(models$KNN, train_x)
knn_cm <- confusionMatrix(knn_class, train_y, positive = "ARDS")
model_results <- rbind(model_results, data.frame(
  Model = "KNN",
  AUC = as.numeric(auc(knn_roc)),
  Accuracy = knn_cm$overall["Accuracy"],
  Sensitivity = knn_cm$byClass["Sensitivity"],
  Specificity = knn_cm$byClass["Specificity"],
  row.names = NULL
))
cat(" AUC=", round(as.numeric(auc(knn_roc)), 4), "\n", sep = "")

# 5.4 Naive Bayes
cat("  [4/6] Naive Bayes...")
set.seed(42)
models$NB <- train(train_x, train_y,
                    method = "nb",
                    trControl = train_control,
                    metric = "ROC")
nb_pred <- predict(models$NB, train_x, type = "prob")[, "ARDS"]
nb_roc <- roc(train_y, nb_pred)
nb_class <- predict(models$NB, train_x)
nb_cm <- confusionMatrix(nb_class, train_y, positive = "ARDS")
model_results <- rbind(model_results, data.frame(
  Model = "Naive Bayes",
  AUC = as.numeric(auc(nb_roc)),
  Accuracy = nb_cm$overall["Accuracy"],
  Sensitivity = nb_cm$byClass["Sensitivity"],
  Specificity = nb_cm$byClass["Specificity"],
  row.names = NULL
))
cat(" AUC=", round(as.numeric(auc(nb_roc)), 4), "\n", sep = "")

# 5.5 LDA
cat("  [5/6] LDA...")
set.seed(42)
models$LDA <- train(train_x, train_y,
                     method = "lda",
                     trControl = train_control,
                     metric = "ROC")
lda_pred <- predict(models$LDA, train_x, type = "prob")[, "ARDS"]
lda_roc <- roc(train_y, lda_pred)
lda_class <- predict(models$LDA, train_x)
lda_cm <- confusionMatrix(lda_class, train_y, positive = "ARDS")
model_results <- rbind(model_results, data.frame(
  Model = "LDA",
  AUC = as.numeric(auc(lda_roc)),
  Accuracy = lda_cm$overall["Accuracy"],
  Sensitivity = lda_cm$byClass["Sensitivity"],
  Specificity = lda_cm$byClass["Specificity"],
  row.names = NULL
))
cat(" AUC=", round(as.numeric(auc(lda_roc)), 4), "\n", sep = "")

# 5.6 LASSO Logistic（已有的）
cat("  [6/6] LASSO Logistic...")
lasso_pred <- lasso_prob
lasso_class <- factor(ifelse(lasso_pred > coords(roc_obj, "best", best.method = "youden")$threshold, 
                              "ARDS", "Control"), levels = c("Control", "ARDS"))
lasso_cm <- confusionMatrix(lasso_class, train_y, positive = "ARDS")
model_results <- rbind(model_results, data.frame(
  Model = "LASSO Logistic",
  AUC = as.numeric(auc(roc_obj)),
  Accuracy = lasso_cm$overall["Accuracy"],
  Sensitivity = lasso_cm$byClass["Sensitivity"],
  Specificity = lasso_cm$byClass["Specificity"],
  row.names = NULL
))
cat(" AUC=", round(as.numeric(auc(roc_obj)), 4), "\n", sep = "")

# ---- 6. 模型性能对比 ----
cat("\n=== 模型性能对比 ===\n")
print(model_results[order(-model_results$AUC), ])

dir.create("results/12_ml_comparison", recursive = TRUE, showWarnings = FALSE)
write.csv(model_results, "results/12_ml_comparison/ml_model_comparison.csv", row.names = FALSE)

# ---- 7. 可视化1：ROC曲线对比 ----
cat("\n=== 绘制ROC曲线对比 ===\n")

roc_list <- list(
  "LASSO Logistic" = roc_obj,
  "Random Forest" = rf_roc,
  "SVM (Radial)" = svm_roc,
  "KNN" = knn_roc,
  "Naive Bayes" = nb_roc,
  "LDA" = lda_roc
)

colors <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6", "#1ABC9C")

png("results/12_ml_comparison/01_roc_comparison.png", width = 1400, height = 1200, res = 150)
par(mar = c(5, 5, 4, 2) + 0.1)

plot(roc_list[[1]], col = colors[1], lwd = 2, main = "ROC Curve Comparison of ML Models")
for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

abline(0, 1, lty = 2, col = "gray")

legend_text <- paste0(names(roc_list), " (AUC=", 
                       round(sapply(roc_list, function(r) as.numeric(auc(r))), 3), ")")
legend("bottomright", legend = legend_text, col = colors, lwd = 2, bty = "n", cex = 0.9)

dev.off()
cat("ROC对比图已保存\n")

# ---- 8. 可视化2：性能指标柱状图 ----
cat("\n=== 绘制性能指标对比图 ===\n")

model_melt <- melt(model_results, id.vars = "Model", 
                    measure.vars = c("AUC", "Accuracy", "Sensitivity", "Specificity"),
                    variable.name = "Metric", value.name = "Value")

p <- ggplot(model_melt, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(title = "Performance Comparison of Machine Learning Models",
       x = "", y = "Score", fill = "Metric") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom",
        axis.text = element_text(size = 10)) +
  ylim(0, 1)

ggsave("results/12_ml_comparison/02_performance_comparison.png", p, width = 8, height = 5, dpi = 150)
cat("性能对比图已保存\n")

# ---- 9. 可视化3：Random Forest特征重要性 ----
cat("\n=== 绘制RF特征重要性 ===\n")

if (!is.null(models$RF)) {
  imp_df <- varImp(models$RF)$importance
  imp_df$Gene <- rownames(imp_df)
  imp_df <- imp_df[order(-imp_df$Overall), ]
  imp_df$Gene <- factor(imp_df$Gene, levels = rev(imp_df$Gene))
  
  p <- ggplot(imp_df, aes(x = Gene, y = Overall, fill = Overall)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    coord_flip() +
    scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
    theme_bw() +
    labs(title = "Random Forest Feature Importance (Ferroptosis Genes)",
         x = "", y = "Importance") +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          legend.position = "none")
  
  ggsave("results/12_ml_comparison/03_rf_feature_importance.png", p, width = 6, height = 5, dpi = 150)
  cat("RF特征重要性图已保存\n")
}

# ---- 10. 多组疾病进程预测 ----
cat("\n=== 多组疾病进程预测 ===\n")

all_x <- as.data.frame(t(expr_gene[selected_genes, ]))
all_rf_prob <- predict(models$RF, all_x, type = "prob")[, "ARDS"]
all_svm_prob <- predict(models$SVM, all_x, type = "prob")[, "ARDS"]

multigroup_pred <- data.frame(
  Group = factor(all_groups, levels = c("Control", "SIRS", "Sepsis", "ARDS")),
  LASSO = all_prob,
  RF = all_rf_prob,
  SVM = all_svm_prob
)

cat("\n各组平均预测概率:\n")
print(aggregate(cbind(LASSO, RF, SVM) ~ Group, data = multigroup_pred, FUN = mean))

# ---- 11. 保存 ----
cat("\n=== 保存结果 ===\n")

save(models, model_results, roc_list, selected_genes,
     file = "data/GEO/GSE32707_ml_comparison_results.RData")

cat("\n=== 机器学习模型比较完成 ===\n")
cat("产出文件:\n")
cat("  results/12_ml_comparison/01_roc_comparison.png - ROC曲线对比\n")
cat("  results/12_ml_comparison/02_performance_comparison.png - 性能指标对比\n")
cat("  results/12_ml_comparison/03_rf_feature_importance.png - RF特征重要性\n")
cat("  results/12_ml_comparison/ml_model_comparison.csv - 模型性能汇总表\n")
cat("\n下一步：运行 10_validation/01_external_validation.R 进行外部验证\n")

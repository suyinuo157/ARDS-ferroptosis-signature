#' ---
#' title: "免疫浸润分析 + 铁死亡-免疫相关性分析"
#' description: "使用ssGSEA评估22种免疫细胞在ALI/ARDS疾病进程中的浸润比例，分析铁死亡评分和Nrf2评分与免疫细胞浸润的相关性"
#' input: "data/GEO/GSE32707_processed.RData, data/GEO/GSE32707_gsea_results.RData, data/GEO/GPL10558.txt"
#' output: "data/GEO/GSE32707_immune_results.RData, results/07_immune/ 下的4张图和2个CSV表"
#' dependencies: "Biobase, ggplot2, pheatmap, reshape2"
#' ---

# ============================================
# 01_免疫浸润分析 + 铁死亡-免疫相关性分析
# 方法：ssGSEA评估22种免疫细胞比例
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("Biobase")) {
  if (!require("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("Biobase", ask = FALSE)
}
library(Biobase)

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(ggplot2)

if (!require("pheatmap")) install.packages("pheatmap", repos = "https://cloud.r-project.org")
library(pheatmap)

if (!require("reshape2")) install.packages("reshape2", repos = "https://cloud.r-project.org")
library(reshape2)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载数据 ----
cat("\n=== 加载数据 ===\n")
load("data/GEO/GSE32707_processed.RData")
load("data/GEO/GSE32707_gsea_results.RData")

expr_log2 <- log2(expr_32707 + 1)

titles <- pdata_32707$title
all_groups <- rep(NA, length(titles))
all_groups[grepl("untreated", titles, ignore.case = TRUE)] <- "Control"
all_groups[grepl("SIRS", titles, ignore.case = TRUE)] <- "SIRS"
all_groups[grepl("Sepsis", titles, ignore.case = TRUE)] <- "Sepsis"
all_groups[grepl("se/ARDS|ARDS", titles, ignore.case = TRUE)] <- "ARDS"
keep <- !is.na(all_groups)
expr_log2 <- expr_log2[, keep]
all_groups <- all_groups[keep]
names(all_groups) <- colnames(expr_log2)

cat("样本数:", ncol(expr_log2), "\n")
cat("各组样本数:", table(all_groups), "\n")

# ---- 3. 读取探针注释 ----
cat("\n=== 读取GPL10558探针注释 ===\n")
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
cat("探针注释完成\n")

# 基因水平表达矩阵（取最大探针）
expr_gene <- aggregate(expr_log2, 
                        by = list(gene = probe_map$symbol[match(rownames(expr_log2), probe_map$probe)]),
                        FUN = max, na.rm = TRUE)
expr_gene <- expr_gene[!is.na(expr_gene$gene) & expr_gene$gene != "", ]
rownames(expr_gene) <- expr_gene$gene
expr_gene <- expr_gene[, -1]

cat("基因水平矩阵维度:", dim(expr_gene), "\n")

# ---- 4. ssGSEA函数（向量化加速版） ----
cat("\n=== 计算ssGSEA免疫评分 ===\n")

ssgsea_score <- function(expr_mat, gene_set) {
  n_samples <- ncol(expr_mat)
  n_genes <- nrow(expr_mat)
  
  common_genes <- intersect(rownames(expr_mat), gene_set)
  cat("  基因集中匹配到的基因数:", length(common_genes), "/", length(gene_set), sep = "")
  
  if (length(common_genes) < 5) return(rep(NA, n_samples))
  
  scores <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    expr_vec <- sort(expr_mat[, i], decreasing = TRUE)
    gene_ranks <- match(common_genes, names(expr_vec))
    gene_ranks <- sort(gene_ranks[!is.na(gene_ranks)])
    
    k <- length(gene_ranks)
    n <- length(expr_vec)
    
    if (k < 5) {
      scores[i] <- NA
      next
    }
    
    total <- 0
    for (t in 0:k) {
      if (t == 0) {
        start <- 1
        end <- gene_ranks[1] - 1
      } else if (t == k) {
        start <- gene_ranks[k]
        end <- n
      } else {
        start <- gene_ranks[t]
        end <- gene_ranks[t+1] - 1
      }
      
      if (start > end) next
      
      len <- end - start + 1
      sum_j <- (start + end) * len / 2
      total <- total + t/k * len - sum_j / (n - k)
    }
    
    scores[i] <- total
  }
  
  names(scores) <- colnames(expr_mat)
  return(scores)
}

# ---- 5. 22种免疫细胞标志基因 ----
cat("\n=== 加载免疫细胞标志基因 ===\n")

immune_gene_sets <- list(
  `B cells naive` = c("IGHD", "TCL1A", "FCER2", "CR2", "CD22", "MS4A1", "CD19", "BANK1"),
  `B cells memory` = c("CD27", "AIM2", "FCRL1", "FCRL2", "FCRL3", "FCRL5"),
  `Plasma cells` = c("SDC1", "CD38", "MZB1", "DERL3", "XBP1", "PRDM1", "IRF4"),
  `T cells CD8` = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7", "KLRD1", "CD6"),
  `T cells CD4 naive` = c("CCR7", "SELL", "LEF1", "TCF7", "S1PR1"),
  `T cells CD4 memory resting` = c("CD4", "IL7R", "CCR7", "SELL", "S1PR1"),
  `T cells CD4 memory activated` = c("CD4", "IL7R", "CD44", "CD69", "HLA-DRA", "HLA-DPA1"),
  `T cells follicular helper` = c("CXCR5", "ICOS", "CD40LG", "SAP", "IL21"),
  `T cells regulatory (Tregs)` = c("FOXP3", "CTLA4", "IL2RA", "TIGIT", "IKZF2", "TNFRSF18"),
  `T cells gamma delta` = c("TRDC", "TRGC1", "TRGC2", "KLRB1", "NKG2D"),
  `NK cells resting` = c("KLRB1", "KLRC1", "KLRD1", "NKG7", "CD160", "FCGR3A"),
  `NK cells activated` = c("GZMB", "PRF1", "IFNG", "CD69", "NKG7", "KLRK1"),
  `Monocytes` = c("CD14", "S100A8", "S100A9", "LYZ", "VCAN", "FCN1", "CD68", "CSF1R"),
  `Macrophages M0` = c("CD68", "CD163", "MSR1", "MRC1", "CSF1R", "CD14"),
  `Macrophages M1` = c("NOS2", "IL1B", "TNF", "IL6", "CD80", "CD86", "IL12B"),
  `Macrophages M2` = c("CD163", "MRC1", "MSR1", "IL10", "TGFB1", "ARG1", "VEGFA"),
  `Dendritic cells resting` = c("FLT3", "ITGAX", "ITGAM", "CD1C", "FCER1A", "CLEC10A"),
  `Dendritic cells activated` = c("CD80", "CD86", "HLA-DRA", "HLA-DPA1", "CD40", "CCR7"),
  `Mast cells resting` = c("CPA3", "TPSAB1", "TPSB2", "KIT", "HDC", "FCER1A"),
  `Mast cells activated` = c("TPSAB1", "TPSB2", "IL13", "IL4", "TNF"),
  `Eosinophils` = c("PRG2", "RNASE2", "RNASE3", "CCR3", "EPX", "SIGLEC8"),
  `Neutrophils` = c("ELANE", "MPO", "PRTN3", "CXCR1", "CXCR2", "FCGR3B", "S100A8", "S100A9", "CD177")
)

cat("免疫细胞基因集数:", length(immune_gene_sets), "\n")

# ---- 6. 计算ssGSEA评分 ----
cat("\n=== 计算免疫细胞ssGSEA评分 ===\n")

immune_scores <- matrix(NA, nrow = length(immune_gene_sets), ncol = ncol(expr_gene))
rownames(immune_scores) <- names(immune_gene_sets)
colnames(immune_scores) <- colnames(expr_gene)

for (i in seq_along(immune_gene_sets)) {
  cell_type <- names(immune_gene_sets)[i]
  cat("  [", i, "/", length(immune_gene_sets), "] ", cell_type, sep = "")
  score <- ssgsea_score(as.matrix(expr_gene), immune_gene_sets[[i]])
  immune_scores[i, ] <- score
  na_count <- sum(is.na(score))
  if (na_count > 0) {
    cat(" (NA:", na_count, ")", sep = "")
  }
  cat(" 完成\n")
}

# 过滤掉全是NA的细胞类型
valid_rows <- apply(immune_scores, 1, function(x) sum(!is.na(x)) > 2)
cat("\n有效免疫细胞类型数:", sum(valid_rows), "/", length(valid_rows), "\n")
if (sum(valid_rows) < length(valid_rows)) {
  cat("被过滤的细胞类型:", names(valid_rows)[!valid_rows], "\n")
}
immune_scores <- immune_scores[valid_rows, ]

# 标准化（按行做min-max归一化）
immune_scores_norm <- t(apply(immune_scores, 1, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] - rng[1] == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}))

cat("免疫浸润评分计算完成\n")

# ---- 7. 可视化1：各组免疫细胞组成热图 ----
cat("\n=== 绘制免疫浸润热图 ===\n")
dir.create("results/07_immune", recursive = TRUE, showWarnings = FALSE)

sample_order <- unlist(lapply(c("Control", "SIRS", "Sepsis", "ARDS"), 
                                function(g) names(all_groups[all_groups == g])))
immune_heatmap <- immune_scores_norm[, sample_order]

# 替换NA为行均值
immune_heatmap <- t(apply(immune_heatmap, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

sample_anno <- data.frame(
  Group = factor(all_groups[sample_order], levels = c("Control", "SIRS", "Sepsis", "ARDS")),
  row.names = sample_order
)

anno_colors <- list(
  Group = c(Control = "#2E86AB", SIRS = "#A23B72", Sepsis = "#F18F01", ARDS = "#C73E1D")
)

png("results/07_immune/01_immune_infiltration_heatmap.png", 
    width = 1400, height = 1000, res = 150)
pheatmap(immune_heatmap,
         annotation_col = sample_anno,
         annotation_colors = anno_colors,
         show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "none",
         color = colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"))(100),
         main = "Immune Cell Infiltration across Disease Stages",
         fontsize_row = 9,
         cellwidth = 3,
         cellheight = 14)
dev.off()
cat("免疫浸润热图已保存\n")

# ---- 8. 可视化2：关键免疫细胞箱线图 ----
cat("\n=== 绘制关键免疫细胞箱线图 ===\n")

key_cells <- c("Neutrophils", "Macrophages M1", "Macrophages M2", 
               "T cells CD8", "T cells regulatory (Tregs)", "NK cells activated")

key_cells <- intersect(key_cells, rownames(immune_scores))

plot_data_list <- list()
for (cell in key_cells) {
  plot_data_list[[cell]] <- data.frame(
    Group = factor(all_groups, levels = c("Control", "SIRS", "Sepsis", "ARDS")),
    Score = as.numeric(immune_scores_norm[cell, ]),
    CellType = cell
  )
}

plot_df <- do.call(rbind, plot_data_list)

p <- ggplot(plot_df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D")) +
  theme_bw() +
  labs(title = "Key Immune Cell Populations across Disease Stages",
       x = "", y = "ssGSEA Score (normalized)") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 9, face = "bold"))

ggsave("results/07_immune/02_key_immune_cell_boxplot.png", p, width = 9, height = 6, dpi = 150)
cat("免疫细胞箱线图已保存\n")

# ---- 9. 铁死亡评分与免疫细胞相关性 ----
cat("\n=== 铁死亡-免疫相关性分析 ===\n")

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

ferro_score <- ssgsea_score(as.matrix(expr_gene), ferroptosis_genes)

# 计算相关性
ferro_immune_cor <- data.frame(
  CellType = rownames(immune_scores_norm),
  Correlation = apply(immune_scores_norm, 1, function(x) cor(x, ferro_score, use = "pairwise.complete.obs")),
  stringsAsFactors = FALSE
)
ferro_immune_cor$PValue <- apply(immune_scores_norm, 1, function(x) {
  cor.test(x, ferro_score, use = "pairwise.complete.obs")$p.value
})
ferro_immune_cor <- ferro_immune_cor[order(-abs(ferro_immune_cor$Correlation)), ]

cat("铁死亡评分与免疫细胞相关性Top 10:\n")
print(head(ferro_immune_cor, 10))

write.csv(ferro_immune_cor, "results/07_immune/ferroptosis_immune_correlation.csv", row.names = FALSE)

# 相关性柱状图
ferro_immune_cor$CellType <- factor(ferro_immune_cor$CellType, 
                                      levels = rev(ferro_immune_cor$CellType))
ferro_immune_cor$Significant <- ifelse(ferro_immune_cor$PValue < 0.05, "Significant", "Not significant")
ferro_immune_cor$Direction <- ifelse(ferro_immune_cor$Correlation > 0, "Positive", "Negative")

p <- ggplot(ferro_immune_cor, aes(x = CellType, y = Correlation, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Positive" = "#C0392B", "Negative" = "#2E86AB")) +
  theme_bw() +
  labs(title = "Correlation between Ferroptosis Score and Immune Cell Infiltration",
       x = "", y = "Pearson Correlation Coefficient",
       fill = "Direction") +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.text = element_text(size = 9),
        legend.position = "bottom")

ggsave("results/07_immune/03_ferroptosis_immune_correlation.png", p, width = 7, height = 6, dpi = 150)
cat("铁死亡-免疫相关性图已保存\n")

# ---- 10. Nrf2评分与免疫的相关性 ----
cat("\n=== Nrf2与免疫相关性 ===\n")

nrf2_target_genes <- c("HMOX1", "NQO1", "GPX4", "SLC7A11", "GCLC", "GCLM", "GSS", 
                       "FTH1", "FTL", "SOD1", "CAT", "PRDX1", "PRDX6")
nrf2_score <- ssgsea_score(as.matrix(expr_gene), nrf2_target_genes)

nrf2_immune_cor <- data.frame(
  CellType = rownames(immune_scores_norm),
  Correlation = apply(immune_scores_norm, 1, function(x) cor(x, nrf2_score, use = "pairwise.complete.obs")),
  PValue = apply(immune_scores_norm, 1, function(x) {
    cor.test(x, nrf2_score, use = "pairwise.complete.obs")$p.value
  }),
  stringsAsFactors = FALSE
)
nrf2_immune_cor <- nrf2_immune_cor[order(-abs(nrf2_immune_cor$Correlation)), ]

cat("Nrf2评分与免疫细胞相关性Top 10:\n")
print(head(nrf2_immune_cor, 10))

write.csv(nrf2_immune_cor, "results/07_immune/Nrf2_immune_correlation.csv", row.names = FALSE)

# ---- 11. 散点图：铁死亡 vs 中性粒细胞/M1巨噬 ----
cat("\n=== 绘制铁死亡-免疫散点图 ===\n")

scatter_cells <- c("Neutrophils", "Macrophages M1")
scatter_list <- list()
for (cell in scatter_cells) {
  scatter_list[[cell]] <- data.frame(
    Ferroptosis = ferro_score,
    ImmuneScore = as.numeric(immune_scores_norm[cell, ]),
    Group = all_groups,
    CellType = cell
  )
}
scatter_df <- do.call(rbind, scatter_list)

p <- ggplot(scatter_df, aes(x = Ferroptosis, y = ImmuneScore, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  facet_wrap(~ CellType, scales = "free", ncol = 2) +
  scale_color_manual(values = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D")) +
  theme_bw() +
  labs(title = "Ferroptosis Score vs Immune Cell Infiltration",
       x = "Ferroptosis ssGSEA Score", y = "Immune Cell Score (normalized)",
       color = "Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom")

ggsave("results/07_immune/04_ferroptosis_immune_scatter.png", p, width = 8, height = 4, dpi = 150)
cat("铁死亡-免疫散点图已保存\n")

# ---- 12. 保存 ----
cat("\n=== 保存结果 ===\n")
save(immune_scores, immune_scores_norm, ferro_score, nrf2_score,
     ferro_immune_cor, nrf2_immune_cor, immune_gene_sets,
     file = "data/GEO/GSE32707_immune_results.RData")

cat("\n=== 免疫浸润分析全部完成 ===\n")
cat("产出文件:\n")
cat("  results/07_immune/01_immune_infiltration_heatmap.png - 免疫浸润热图\n")
cat("  results/07_immune/02_key_immune_cell_boxplot.png - 关键免疫细胞箱线图\n")
cat("  results/07_immune/03_ferroptosis_immune_correlation.png - 铁死亡-免疫相关性\n")
cat("  results/07_immune/04_ferroptosis_immune_scatter.png - 铁死亡-免疫散点图\n")
cat("  results/07_immune/ferroptosis_immune_correlation.csv - 相关性表\n")
cat("  results/07_immune/Nrf2_immune_correlation.csv - Nrf2-免疫相关性表\n")
cat("\n下一步：运行 08_ceRNA_network/01_ceRNA_network.R 进行ceRNA调控网络分析\n")

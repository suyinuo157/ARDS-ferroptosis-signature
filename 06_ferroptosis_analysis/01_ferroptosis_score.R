#' ---
#' title: "铁死亡关键基因表达分析"
#' description: "展示铁死亡关键基因（Driver和Suppressor）在Control/SIRS/Sepsis/ARDS四组疾病进程中的表达变化，绘制热图和关键基因箱线图"
#' input: "data/GEO/GSE32707_processed.RData, data/GEO/GPL10558.txt"
#' output: "results/06_ferroptosis/ferroptosis_gene_heatmap.png, results/06_ferroptosis/*_boxplot.png（6个关键基因）"
#' dependencies: "pheatmap, ggplot2"
#' ---

# ============================================
# 01_铁死亡关键基因表达热图
# ============================================
# 功能：展示铁死亡关键基因在Control/SIRS/Sepsis/ARDS四组中的表达变化
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
library(pheatmap)
library(ggplot2)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载数据 ----
cat("\n=== 加载数据 ===\n")

load("data/GEO/GSE32707_processed.RData")
expr_all <- log2(expr_32707 + 1)

# 从title提取分组
titles <- pdata_32707$title
all_groups <- rep(NA, length(titles))
all_groups[grepl("untreated", titles, ignore.case = TRUE)] <- "Control"
all_groups[grepl("SIRS", titles, ignore.case = TRUE)] <- "SIRS"
all_groups[grepl("Sepsis", titles, ignore.case = TRUE)] <- "Sepsis"
all_groups[grepl("se/ARDS|ARDS", titles, ignore.case = TRUE)] <- "ARDS"
keep <- !is.na(all_groups)
expr_all <- expr_all[, keep]
all_groups <- all_groups[keep]

cat("样本数:", ncol(expr_all), "\n")
print(table(all_groups))

# ---- 3. 铁死亡关键基因探针列表 ----
ferro_probe_map <- data.frame(
  probe = c(
    # 铁死亡抑制基因（Suppressor）
    "ILMN_1734353",  # GPX4
    "ILMN_1655229",  # SLC7A11 (xCT)
    "ILMN_1790909",  # NFE2L2 (Nrf2)
    "ILMN_1674236",  # HSPB1
    "ILMN_1657634",  # FANCD2
    "ILMN_2122952",  # CISD1
    "ILMN_1796397",  # CISD2
    "ILMN_2175601",  # VDAC1
    "ILMN_1744822",  # BECN1
    "ILMN_1697559",  # G6PD
    "ILMN_1689329",  # SCD
    "ILMN_1753342",  # SAT1
    "ILMN_1712639",  # AIFM2 (FSP1)
    "ILMN_1764873",  # ELAVL1
    "ILMN_1679041",  # SLC3A2
    "ILMN_1805225",  # LPCAT3
    # 铁死亡促进基因（Driver）
    "ILMN_1683598",  # ACSL4
    "ILMN_1674243",  # TFRC
    "ILMN_1800512",  # HMOX1
    "ILMN_1773906",  # NCOA4
    "ILMN_1677511",  # PTGS2
    "ILMN_1739241",  # CHAC1
    "ILMN_1651466",  # ALOX5
    "ILMN_1713731",  # ALOX12
    "ILMN_1676042",  # ALOX15
    "ILMN_1677768",  # POR
    "ILMN_1765355",  # NOX1
    "ILMN_1761577",  # STEAP3
    "ILMN_1745034",  # SLC11A2 (DMT1)
    "ILMN_1764629",  # SLC39A14
    "ILMN_1696066",  # CARS
    "ILMN_1692535"   # DPP4
  ),
  gene = c(
    "GPX4", "SLC7A11", "NFE2L2", "HSPB1", "FANCD2",
    "CISD1", "CISD2", "VDAC1", "BECN1", "G6PD",
    "SCD", "SAT1", "AIFM2", "ELAVL1", "SLC3A2", "LPCAT3",
    "ACSL4", "TFRC", "HMOX1", "NCOA4", "PTGS2",
    "CHAC1", "ALOX5", "ALOX12", "ALOX15", "POR",
    "NOX1", "STEAP3", "SLC11A2", "SLC39A14", "CARS", "DPP4"
  ),
  role = c(
    rep("Suppressor", 16),
    rep("Driver", 16)
  ),
  stringsAsFactors = FALSE
)

# 过滤出在表达矩阵中的探针
ferro_probe_map <- ferro_probe_map[ferro_probe_map$probe %in% rownames(expr_all), ]
cat("铁死亡基因数:", nrow(ferro_probe_map), "\n")

# ---- 4. 提取表达矩阵并标准化 ----
cat("\n=== 生成热图数据 ===\n")

ferro_expr <- expr_all[ferro_probe_map$probe, ]
rownames(ferro_expr) <- ferro_probe_map$gene

# 给all_groups加上名字
names(all_groups) <- colnames(expr_all)

# 按行（基因）做Z-score标准化
ferro_expr_scaled <- t(scale(t(ferro_expr)))

# ---- 5. 样本按分组排序 ----
group_order <- c("Control", "SIRS", "Sepsis", "ARDS")
sample_order <- unlist(lapply(group_order, function(g) {
  names(all_groups[all_groups == g])
}))
ferro_expr_scaled <- ferro_expr_scaled[, sample_order]

# ---- 6. 注释信息 ----
# 样本分组注释
sample_anno <- data.frame(
  Group = factor(all_groups[sample_order], levels = group_order),
  row.names = sample_order
)

# 基因角色注释
gene_anno <- data.frame(
  Role = factor(ferro_probe_map$role, levels = c("Driver", "Suppressor")),
  row.names = ferro_probe_map$gene
)

# 颜色
anno_colors <- list(
  Group = c(Control = "#2E86AB", SIRS = "#A23B72", Sepsis = "#F18F01", ARDS = "#C73E1D"),
  Role = c(Driver = "#E63946", Suppressor = "#2A9D8F")
)

# ---- 7. 画热图 ----
cat("\n=== 绘制热图 ===\n")

dir.create("results/06_ferroptosis", recursive = TRUE, showWarnings = FALSE)

png("results/06_ferroptosis/ferroptosis_gene_heatmap.png", 
    width = 1000, height = 800, res = 150)

pheatmap(ferro_expr_scaled,
         annotation_col = sample_anno,
         annotation_row = gene_anno,
         annotation_colors = anno_colors,
         show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "none",
         color = colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"))(100),
         main = "Ferroptosis-related gene expression in ALI/ARDS progression",
         fontsize_row = 9,
         cellwidth = 2,
         cellheight = 12,
         gaps_col = cumsum(c(sum(all_groups=="Control"), 
                              sum(all_groups=="SIRS"), 
                              sum(all_groups=="Sepsis"))))

dev.off()

cat("铁死亡基因热图已保存\n")

# ---- 8. 画箱线图（几个关键基因） ----
cat("\n=== 绘制关键基因箱线图 ===\n")

key_genes <- c("HSPB1", "SLC7A11", "NFE2L2", "GPX4", "PTGS2", "TFRC")

for (gene in key_genes) {
  probe <- ferro_probe_map$probe[ferro_probe_map$gene == gene]
  if (length(probe) == 0) next
  
  plot_data <- data.frame(
    Group = factor(all_groups, levels = group_order),
    Expression = as.numeric(expr_all[probe, ])
  )
  
  p <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    scale_fill_manual(values = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D")) +
    theme_bw() +
    labs(title = paste(gene, "expression"), x = "", y = "Expression (log2)") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12),
          axis.text = element_text(size = 10))
  
  ggsave(paste0("results/06_ferroptosis/", gene, "_boxplot.png"), 
         p, width = 4, height = 4, dpi = 150)
}

cat("关键基因箱线图已保存\n")

cat("\n=== 铁死亡基因可视化完成 ===\n")
cat("\n下一步：运行 06_ferroptosis_analysis/02_Nrf2_regulation.R 进行Nrf2调控分析\n")

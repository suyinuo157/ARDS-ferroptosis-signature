#' ---
#' title: "GSEA富集分析 + ssGSEA铁死亡评分"
#' description: "基于GPL10558探针注释进行KEGG和GO的GSEA富集分析，同时使用ssGSEA计算铁死亡评分并验证其在疾病进程中的变化"
#' input: "data/GEO/GSE32707_diff_analysis.RData, data/GEO/GSE32707_processed.RData, data/GEO/GPL10558.txt"
#' output: "data/GEO/GSE32707_gsea_results.RData, results/03_GSEA/KEGG_GSEA_results.csv, results/03_GSEA/01_KEGG_GSEA_dotplot.png, results/03_GSEA/02_GO_GSEA_dotplot.png, results/03_GSEA/03_ferroptosis_ssGSEA_boxplot.png, results/03_GSEA/04_ferroptosis_gene_heatmap.png"
#' dependencies: "clusterProfiler, org.Hs.eg.db, ggplot2, ggpubr, pheatmap"
#' ---

# ============================================
# 01_GSEA富集分析 + ssGSEA铁死亡评分
# ============================================
# 功能：GSEA富集分析（KEGG+GO）+ ssGSEA铁死亡评分
# 说明：从GPL10558.txt提取探针注释
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(pheatmap)

cat("=== 包加载完成 ===\n")

# ---- 手动实现ssGSEA函数 ----
ssgsea_score <- function(expr_matrix, gene_set) {
  n_genes <- nrow(expr_matrix)
  n_samples <- ncol(expr_matrix)
  gene_set <- intersect(gene_set, rownames(expr_matrix))
  n_set <- length(gene_set)
  
  if (n_set < 5) {
    warning("基因集太小，结果不可靠")
    return(rep(NA, n_samples))
  }
  
  scores <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    ranked <- sort(expr_matrix[, i], decreasing = TRUE)
    ranks <- 1:n_genes
    in_set <- rownames(expr_matrix) %in% gene_set
    sum_in <- sum(ranks[in_set[match(names(ranked), rownames(expr_matrix))]])
    expected <- n_set * (n_genes + 1) / 2
    scores[i] <- (sum_in - expected) / (n_genes - n_set)
  }
  
  names(scores) <- colnames(expr_matrix)
  return(scores)
}

# ---- 2. 读取GPL10558注释 ----
cat("\n=== 读取GPL10558探针注释 ===\n")

# 读取SOFT格式文件，找到数据表格开始的位置
gpl_lines <- readLines("data/GEO/GPL10558.txt")

# 找表头行（包含ID和Gene Symbol的行）
table_start <- grep("^ID", gpl_lines)[1]
cat("表格起始行:", table_start, "\n")

# 读取数据表格
probe_table <- read.table("data/GEO/GPL10558.txt",
                          skip = table_start - 1,
                          header = TRUE,
                          sep = "\t",
                          fill = TRUE,
                          quote = "",
                          comment.char = "",
                          stringsAsFactors = FALSE)

cat("探针总数:", nrow(probe_table), "\n")

# 构建探针-Symbol映射
# ILMN_GENE列是基因Symbol
if ("ILMN_GENE" %in% colnames(probe_table)) {
  probe_map <- data.frame(
    probe = probe_table$ID,
    symbol = probe_table$ILMN_GENE,
    stringsAsFactors = FALSE
  )
} else if ("Gene.Symbol" %in% colnames(probe_table)) {
  probe_map <- data.frame(
    probe = probe_table$ID,
    symbol = probe_table$Gene.Symbol,
    stringsAsFactors = FALSE
  )
} else {
  gene_col <- grep("gene|symbol", colnames(probe_table), ignore.case = TRUE, value = TRUE)[1]
  cat("使用基因列:", gene_col, "\n")
  probe_map <- data.frame(
    probe = probe_table$ID,
    symbol = probe_table[[gene_col]],
    stringsAsFactors = FALSE
  )
}

# 去除空值和NA
probe_map <- probe_map[!is.na(probe_map$symbol) & probe_map$symbol != "" & probe_map$symbol != "---", ]
cat("有基因注释的探针数:", nrow(probe_map), "\n")

# ---- 3. 加载差异分析结果 ----
cat("\n=== 加载差异分析结果 ===\n")

load("data/GEO/GSE32707_diff_analysis.RData")

# 给deg_all加symbol
deg_df <- as.data.frame(deg_all)
deg_df$probe <- rownames(deg_df)
deg_df$symbol <- probe_map$symbol[match(deg_df$probe, probe_map$probe)]

# 去除NA
deg_mapped <- deg_df[!is.na(deg_df$symbol) & deg_df$symbol != "", ]
cat("映射后基因数:", nrow(deg_mapped), "\n")

# 多探针对应同一基因时，取logFC绝对值最大的
deg_mapped$abs_logFC <- abs(deg_mapped$logFC)
deg_mapped <- deg_mapped[order(-deg_mapped$abs_logFC), ]
deg_mapped <- deg_mapped[!duplicated(deg_mapped$symbol), ]
rownames(deg_mapped) <- deg_mapped$symbol

cat("去重后基因数:", nrow(deg_mapped), "\n")

# ---- 4. GSEA分析 ----
cat("\n=== GSEA富集分析 ===\n")

# 准备基因列表（命名的logFC向量，按从大到小排序）
gene_list <- deg_mapped$logFC
names(gene_list) <- rownames(deg_mapped)
gene_list <- sort(gene_list, decreasing = TRUE)

cat("GSEA基因列表长度:", length(gene_list), "\n")

# Symbol转Entrez ID
gene_df <- bitr(names(gene_list), 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)
cat("成功转换Entrez ID的基因数:", nrow(gene_df), "\n")

# 构建Entrez ID的gene list
gene_list_entrez <- gene_list[gene_df$SYMBOL]
names(gene_list_entrez) <- gene_df$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# ---- 4a. KEGG GSEA ----
cat("\n=== KEGG GSEA ===\n")

set.seed(123)
gsea_kegg <- gseKEGG(geneList = gene_list_entrez,
                     organism = "hsa",
                     pvalueCutoff = 0.5,
                     pAdjustMethod = "BH",
                     verbose = FALSE)

cat("KEGG富集通路数:", nrow(gsea_kegg), "\n")
if (nrow(gsea_kegg) > 0) {
  print(head(gsea_kegg[, c("ID", "Description", "NES", "pvalue", "p.adjust")], 15))
}

# 保存结果
dir.create("results/03_GSEA", recursive = TRUE, showWarnings = FALSE)

if (nrow(gsea_kegg) > 0) {
  write.csv(as.data.frame(gsea_kegg), "results/03_GSEA/KEGG_GSEA_results.csv", row.names = FALSE)
  
  # 气泡图（前20个）
  top_kegg <- gsea_kegg[order(gsea_kegg$pvalue), ][1:min(20, nrow(gsea_kegg)), ]
  
  png("results/03_GSEA/01_KEGG_GSEA_dotplot.png", width = 1000, height = 800, res = 150)
  print(dotplot(gsea_kegg, showCategory = 20, title = "KEGG GSEA - ARDS vs Control"))
  dev.off()
  
  cat("KEGG GSEA图已保存\n")
}

# ---- 4b. GO GSEA ----
cat("\n=== GO GSEA (BP) ===\n")

set.seed(123)
gsea_go <- gseGO(geneList = gene_list_entrez,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pvalueCutoff = 0.5,
                 verbose = FALSE)

cat("GO富集通路数:", nrow(gsea_go), "\n")
if (nrow(gsea_go) > 0) {
  print(head(gsea_go[, c("ID", "Description", "NES", "pvalue", "p.adjust")], 15))
}

# 保存结果
if (nrow(gsea_go) > 0) {
  write.csv(as.data.frame(gsea_go), "results/03_GSEA/GO_GSEA_results.csv", row.names = FALSE)
  
  png("results/03_GSEA/02_GO_GSEA_dotplot.png", width = 1000, height = 800, res = 150)
  print(dotplot(gsea_go, showCategory = 20, title = "GO BP GSEA - ARDS vs Control"))
  dev.off()
  
  cat("GO GSEA图已保存\n")
}

# ---- 5. 铁死亡基因集验证 ----
cat("\n=== 铁死亡基因集验证 ===\n")

# 铁死亡基因列表
ferroptosis_genes <- c(
  "ACSL4", "ACSL5", "ACSL6", "AIFM2", "ATP5MC3", "ATP5G3",
  "BECN1", "CARS", "CBS", "CHAC1", "CHMP5", "CHMP6",
  "COQ2", "CP", "CRLS1", "CISD1", "CISD2", "DPP4",
  "ELAVL1", "EMC2", "FANCD2", "FDFT1", "FER", "FTH1",
  "FTL", "G6PD", "GCLC", "GCLM", "GLS2", "GPX4",
  "GSS", "GSTP1", "HMGCR", "HMOX1", "HSPB1", "ITFG2",
  "LPCAT3", "LAMP2", "MAP1LC3A", "MTDH", "MT1G",
  "NCOA4", "NFE2L2", "NOX1", "OTUB1", "PCBP2",
  "PEBP1", "PGRMC2", "PHKG2", "POR", "PTGS2",
  "RPL8", "SAT1", "SCD", "SLC11A2", "SLC38A1",
  "SLC39A14", "SLC39A7", "SLC7A11", "SLC3A2",
  "SP1", "SQLE", "SQSTM1", "STEAP3", "TF", "TFRC",
  "TFR2", "TP53", "VCP", "VDAC1", "VDAC2", "VDAC3",
  "ZFP36", "ALOX5", "ALOX12", "ALOX15", "FSP1",
  "SLC1A5", "GOT1", "CRYAB", "ANXA7", "HSP90AA1"
)
ferroptosis_genes <- unique(ferroptosis_genes)

# 差异基因中的铁死亡基因
ferro_in_deg <- intersect(ferroptosis_genes, rownames(deg_mapped))
cat("差异表达的铁死亡基因数:", length(ferro_in_deg), "\n")

if (length(ferro_in_deg) > 0) {
  cat("铁死亡基因列表:\n")
  print(deg_mapped[ferro_in_deg, c("logFC", "P.Value", "adj.P.Val")])
  
  # Fisher精确检验
  all_genes <- rownames(deg_mapped)
  deg_genes <- rownames(deg_mapped)[deg_mapped$adj.P.Val < 0.05]
  
  ferro_deg <- intersect(ferroptosis_genes, deg_genes)
  ferro_not_deg <- setdiff(intersect(ferroptosis_genes, all_genes), deg_genes)
  nonferro_deg <- setdiff(deg_genes, ferroptosis_genes)
  nonferro_not_deg <- length(all_genes) - length(deg_genes) - length(ferro_not_deg)
  
  fisher_matrix <- matrix(c(
    length(ferro_deg), length(ferro_not_deg),
    length(nonferro_deg), nonferro_not_deg
  ), nrow = 2, byrow = TRUE)
  rownames(fisher_matrix) <- c("Ferroptosis", "Non-ferroptosis")
  colnames(fisher_matrix) <- c("DEG", "Non-DEG")
  
  cat("\nFisher精确检验列联表:\n")
  print(fisher_matrix)
  cat("Fisher检验p值:", fisher.test(fisher_matrix)$p.value, "\n")
}

# ---- 6. ssGSEA铁死亡评分 ----
cat("\n=== ssGSEA铁死亡评分 ===\n")

# 重新加载全部样本数据（四组）
load("data/GEO/GSE32707_processed.RData")
expr_all_gsea <- log2(exprs(gse32707) + 1)

# 从title提取分组
titles <- pData(gse32707)$title
all_groups_gsea <- rep(NA, length(titles))
all_groups_gsea[grepl("untreated", titles, ignore.case = TRUE)] <- "Control"
all_groups_gsea[grepl("SIRS", titles, ignore.case = TRUE)] <- "SIRS"
all_groups_gsea[grepl("Sepsis", titles, ignore.case = TRUE)] <- "Sepsis"
all_groups_gsea[grepl("se/ARDS|ARDS", titles, ignore.case = TRUE)] <- "ARDS"
keep_gsea <- !is.na(all_groups_gsea)
expr_all_gsea <- expr_all_gsea[, keep_gsea]
all_groups_gsea <- all_groups_gsea[keep_gsea]

# 找铁死亡基因对应的探针
ferro_probes <- probe_map$probe[probe_map$symbol %in% ferroptosis_genes]
ferro_probes <- intersect(ferro_probes, rownames(expr_all_gsea))
cat("铁死亡基因对应探针数:", length(ferro_probes), "\n")

cat("运行ssGSEA...\n")
ferroptosis_scores <- ssgsea_score(expr_all_gsea, ferro_probes)

# 各组评分
cat("各组铁死亡ssGSEA评分均值:\n")
print(tapply(ferroptosis_scores, all_groups_gsea, mean))

# Control vs ARDS t检验
control_idx <- all_groups_gsea == "Control"
ards_idx <- all_groups_gsea == "ARDS"
pval_ssgsea <- t.test(ferroptosis_scores[control_idx], ferroptosis_scores[ards_idx])$p.value
cat("Control vs ARDS p值:", pval_ssgsea, "\n")

# 画箱线图
plot_data <- data.frame(
  Group = factor(all_groups_gsea, levels = c("Control", "SIRS", "Sepsis", "ARDS")),
  Score = as.numeric(ferroptosis_scores)
)

p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  scale_fill_manual(values = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D")) +
  theme_bw() +
  labs(title = "Ferroptosis ssGSEA Score across Disease Stages",
       x = "", y = "ssGSEA Score") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 10))

ggsave("results/03_GSEA/03_ferroptosis_ssGSEA_boxplot.png", p, width = 5, height = 4, dpi = 150)
cat("ssGSEA箱线图已保存\n")

# ---- 7. 铁死亡基因热图 ----
cat("\n=== 铁死亡基因热图 ===\n")

# 取显著差异的铁死亡基因
ferro_sig <- deg_mapped[ferro_in_deg, ]
ferro_sig <- ferro_sig[order(ferro_sig$logFC, decreasing = TRUE), ]

if (nrow(ferro_sig) > 5) {
  # 取这些基因的表达
  ferro_probes_sig <- probe_map$probe[match(rownames(ferro_sig), probe_map$symbol)]
  ferro_expr_sig <- expr_all_gsea[ferro_probes_sig, ]
  rownames(ferro_expr_sig) <- rownames(ferro_sig)
  
  # Z-score标准化
  ferro_scaled <- t(scale(t(ferro_expr_sig)))
  
  # 样本排序
  sample_order <- unlist(lapply(c("Control", "SIRS", "Sepsis", "ARDS"), 
                                  function(g) names(all_groups_gsea[all_groups_gsea == g])))
  ferro_scaled <- ferro_scaled[, sample_order]
  
  # 注释
  sample_anno <- data.frame(
    Group = factor(all_groups_gsea[sample_order], levels = c("Control", "SIRS", "Sepsis", "ARDS")),
    row.names = sample_order
  )
  
  anno_colors <- list(
    Group = c(Control = "#2E86AB", SIRS = "#A23B72", Sepsis = "#F18F01", ARDS = "#C73E1D")
  )
  
  png("results/03_GSEA/04_ferroptosis_gene_heatmap.png", 
      width = 1000, height = 800, res = 150)
  pheatmap(ferro_scaled,
           annotation_col = sample_anno,
           annotation_colors = anno_colors,
           show_colnames = FALSE,
           show_rownames = TRUE,
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           scale = "none",
           color = colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"))(100),
           main = "Ferroptosis-related gene expression",
           fontsize_row = 9,
           cellwidth = 2,
           cellheight = 10)
  dev.off()
  
  cat("铁死亡基因热图已保存\n")
}

# ---- 保存结果 ----
save(gsea_kegg, gsea_go, deg_mapped, probe_map, ferroptosis_scores,
     file = "data/GEO/GSE32707_gsea_results.RData")

cat("\n=== GSEA分析完成 ===\n")
cat("\n下一步：运行 03_functional_enrichment/02_extended_enrichment.R 进行扩展富集分析\n")

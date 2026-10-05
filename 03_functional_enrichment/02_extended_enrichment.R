#' ---
#' title: "扩展富集分析：铁死亡多维度 + 免疫相关 + 自定义基因集GSEA"
#' description: "对铁死亡相关的多个功能维度（铁死亡调控、铁代谢、脂质代谢、谷胱甘肽系统、Nrf2通路、铁自噬、炎症反应等）进行GSEA富集分析和ssGSEA评分，并绘制火山图标注铁死亡基因"
#' input: "data/GEO/GSE32707_diff_analysis.RData, data/GEO/GSE32707_processed.RData, data/GEO/GPL10558.txt"
#' output: "data/GEO/GSE32707_extended_enrichment.RData, results/13_extended_enrichment/ 下的6张图和结果CSV"
#' dependencies: "clusterProfiler, org.Hs.eg.db, ggplot2, pheatmap, enrichplot, ggrepel"
#' ---

# ============================================================================
# 02. 扩展富集分析（铁死亡多维度 + 免疫相关 + 自定义基因集GSEA）
# ============================================================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

cat("=== 扩展富集分析 ===\n\n")

# 加载包
cat("加载包...\n")
required_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "ggplot2", "pheatmap", "enrichplot", "ggrepel")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}
cat("包加载完成\n\n")

# 创建输出目录
dir.create("results/13_extended_enrichment", recursive = TRUE, showWarnings = FALSE)

# 加载数据
cat("加载数据...\n")
load("data/GEO/GSE32707_diff_analysis.RData")

# 检查deg_all的行名是探针还是基因
head_rownames <- head(rownames(deg_all))
cat("前6个行名:", paste(head_rownames, collapse=", "), "\n")

# 如果是探针ID，用GPL注释转成基因名
if (grepl("^ILMN_", head_rownames[1])) {
  cat("检测到Illumina探针ID，加载GPL注释...\n")
  
  gpl_file <- "data/GEO/GPL10558.txt"
  if (file.exists(gpl_file)) {
    lines <- readLines(gpl_file, warn = FALSE)
    start_idx <- grep("!platform_table_begin", lines)
    end_idx <- grep("!platform_table_end", lines)
    
    if (length(start_idx) > 0 && length(end_idx) > 0) {
      gpl <- read.delim(text = lines[(start_idx + 1):(end_idx - 1)],
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                        comment.char = "", check.names = FALSE)
      
      symbol_col <- grep("Symbol|gene_symbol|Gene.Symbol", colnames(gpl), value = TRUE, ignore.case = TRUE)[1]
      probe_col <- grep("ID|Probe_ID", colnames(gpl), value = TRUE, ignore.case = TRUE)[1]
      
      if (length(symbol_col) > 0 && length(probe_col) > 0) {
        probe_map <- data.frame(
          probe = gpl[, probe_col],
          symbol = gpl[, symbol_col],
          stringsAsFactors = FALSE
        )
        
        deg_df <- as.data.frame(deg_all)
        deg_df$probe <- rownames(deg_df)
        deg_df$symbol <- probe_map$symbol[match(deg_df$probe, probe_map$probe)]
        deg_df <- deg_df[!is.na(deg_df$symbol) & deg_df$symbol != "", ]
        
        deg_df <- do.call(rbind, by(deg_df, deg_df$symbol, function(x) {
          x[which.min(x$adj.P.Val), , drop = FALSE]
        }))
        rownames(deg_df) <- deg_df$symbol
        deg_all_gene <- deg_df[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
        cat("基因水平DEG数:", nrow(deg_all_gene), "\n")
      }
    }
  }
} else {
  deg_all_gene <- deg_all
}

# 准备GSEA的geneList
gsea_gene <- deg_all_gene$logFC
names(gsea_gene) <- rownames(deg_all_gene)
gsea_gene <- sort(gsea_gene, decreasing = TRUE)
gsea_gene <- gsea_gene[!is.na(gsea_gene) & !is.infinite(gsea_gene)]
gsea_gene <- gsea_gene[!duplicated(names(gsea_gene))]
cat("GSEA基因数:", length(gsea_gene), "\n\n")

# === 1. 铁死亡相关多维度基因集GSEA ===
cat("=== 1. 铁死亡多维度基因集GSEA ===\n")

ferroptosis_gene_sets <- list(
  "Ferroptosis Regulators" = c("ACSL4", "GPX4", "SLC7A11", "ACSL3", "LPCAT3", "POR", 
                                "LOX", "ALOX15", "ALOX5", "NOX1", "NOX2", "NOX4",
                                "FSP1", "DHODH", "GCH1", "CHAC1", "PTGS2", "CHAC2",
                                "SLC3A2", "GCLC", "GCLM", "GSS", "GSR"),
  "Iron Metabolism" = c("FTH1", "FTL", "TF", "TFRC", "SLC40A1", "HMOX1", 
                        "STEAP3", "IREB2", "NCOA4", "DMT1", "SLC11A2",
                        "HAMP", "TFR2", "HJV", "BMP6"),
  "Lipid Metabolism (PUFA)" = c("ACSL4", "LPCAT3", "PLA2G6", "ELOVL6", "FADS2", 
                                 "SCD", "ACACA", "FASN", "DGAT1", "DGAT2",
                                 "SOAT1", "PTGS2", "ALOX5", "ALOX12", "ALOX15"),
  "Glutathione System" = c("GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX2", "GPX3",
                           "GPX4", "SLC7A11", "SLC3A2", "CHAC1", "CHAC2",
                           "GLS", "GLS2", "SLC1A5"),
  "Nrf2 Antioxidant Pathway" = c("NFE2L2", "KEAP1", "MAFG", "MAFK", "NQO1", "HMOX1",
                                 "GCLC", "GCLM", "TXNRD1", "TXN", "SLC7A11", 
                                 "GPX4", "FTH1", "FTL", "ABCC1", "GSTP1"),
  "Ferritinophagy" = c("ATG5", "ATG7", "BECN1", "MAP1LC3A", "MAP1LC3B",
                       "SQSTM1", "NCOA4", "FTH1", "FTL", "TFRC", "ATG16L1"),
  "Inflammatory Response" = c("NFKB1", "RELA", "TNF", "IL1B", "IL6", "CXCL8",
                              "PTGS2", "HMOX1", "SOD2", "CAT", "NOS2",
                              "STAT1", "STAT3", "IL10", "TGFB1"),
  "Cell Death Pathways" = c("BAX", "BAK1", "BCL2", "BCL2L1", "CASP3", "CASP8",
                            "CASP9", "APAF1", "CYCS", "RIPK1", "RIPK3", "MLKL",
                            "GSDMD", "GSDME", "PARP1", "AIFM1"),
  "Immune Cell Markers" = c("CD4", "CD8A", "CD19", "CD14", "CD68", "CD163",
                            "MS4A1", "NCAM1", "NKG7", "GZMB", "PRF1",
                            "FOXP3", "IL2RA", "CSF3R", "S100A8", "S100A9"),
  "Cytokine-Cytokine Receptor" = c("TNF", "IL1B", "IL6", "IL8", "CXCL10", "CCL2",
                                   "CCL3", "CCL4", "CCL5", "IFNG", "IL10", "TGFB1",
                                   "IL1R1", "IL6R", "TNFR1", "TNFRSF1A"),
  "Oxidative Phosphorylation" = c("ND1", "ND2", "ND4", "CYTB", "COX1", "COX2",
                                   "ATP5A1", "ATP5B", "SDHA", "SDHB",
                                   "UQCRC2", "UQCRC1"),
  "Autophagy" = c("ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1", "BECN1",
                  "MAP1LC3A", "MAP1LC3B", "SQSTM1", "ULK1", "ATG13",
                  "PINK1", "PRKN", "BNIP3", "BNIP3L"),
  "p53 Signaling" = c("TP53", "MDM2", "CDKN1A", "BAX", "PUMA", "BBC3",
                       "NOXA", "PMAIP1", "GADD45A", "SESN1", "SESN2",
                       "SLC7A11", "GPX4", "DPP4"),
  "Hypoxia Response" = c("HIF1A", "EPAS1", "VEGFA", "EPO", "HK2", "LDHA",
                         "PDK1", "GLUT1", "SLC2A1", "NDRG1", "CA9",
                         "BNIP3", "BNIP3L", "HMOX1"),
  "Mitochondrial Function" = c("PPARGC1A", "TFAM", "POLG", "DRP1", "DNM1L",
                               "MFN1", "MFN2", "OPA1", "PINK1", "PRKN",
                               "CYCS", "APAF1", "BAX", "BAK1")
)

# 转换为TERM2GENE格式
term2gene_ferro <- data.frame(
  term = rep(names(ferroptosis_gene_sets), sapply(ferroptosis_gene_sets, length)),
  gene = unlist(ferroptosis_gene_sets),
  stringsAsFactors = FALSE
)

# GSEA
set.seed(42)
gsea_ferro <- GSEA(gsea_gene, TERM2GENE = term2gene_ferro,
                   pvalueCutoff = 0.5, minGSSize = 5, maxGSSize = 200,
                   eps = 0, by = "fgsea")

if (!is.null(gsea_ferro) && nrow(gsea_ferro) > 0) {
  cat("铁死亡相关GSEA结果:\n")
  print(as.data.frame(gsea_ferro)[, c("ID", "NES", "pvalue", "p.adjust")])
  
  # 图1：GSEA点图
  p1 <- dotplot(gsea_ferro, showCategory = 15, split = ".sign",
                title = "Ferroptosis-Related Gene Set Enrichment Analysis") +
    facet_grid(.~.sign) +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(size = 12, face = "bold"))
  ggsave("results/13_extended_enrichment/01_ferroptosis_gsea_dotplot.png", p1, width = 14, height = 10, dpi = 150)
  cat("铁死亡GSEA点图已保存\n")
  
  # 图2：GSEA曲线（Top 6）
  p2 <- gseaplot2(gsea_ferro, geneSetID = 1:min(6, nrow(gsea_ferro)),
                  title = "Top Ferroptosis-Related Gene Sets",
                  pvalue_table = TRUE)
  ggsave("results/13_extended_enrichment/02_ferroptosis_gsea_curves.png", p2, width = 12, height = 10, dpi = 150)
  cat("铁死亡GSEA曲线已保存\n")
  
  # 保存结果
  write.csv(as.data.frame(gsea_ferro),
            "results/13_extended_enrichment/ferroptosis_gsea_results.csv", row.names = FALSE)
}

# === 2. 铁死亡评分分组差异表达热图 ===
cat("\n=== 2. 铁死亡关键通路基因热图 ===\n")

load("data/GEO/GSE32707_processed.RData")

# 用expr_32707和pdata_32707
if (exists("expr_gene") && exists("pdata_32707")) {
  expr_mat <- expr_gene
} else if (exists("expr_32707") && exists("pdata_32707")) {
  cat("从探针水平转基因水平...\n")
  if (exists("probe_map")) {
    expr_df <- as.data.frame(expr_32707)
    expr_df$symbol <- probe_map$symbol[match(rownames(expr_df), probe_map$probe)]
    expr_df <- expr_df[!is.na(expr_df$symbol) & expr_df$symbol != "", ]
    
    expr_df$mean_expr <- rowMeans(expr_df[, 1:ncol(expr_32707)], na.rm = TRUE)
    expr_df <- do.call(rbind, by(expr_df, expr_df$symbol, function(x) {
      x[which.max(x$mean_expr), , drop = FALSE]
    }))
    rownames(expr_df) <- expr_df$symbol
    expr_mat <- as.matrix(expr_df[, 1:ncol(expr_32707)])
    cat("基因水平矩阵:", dim(expr_mat), "\n")
  }
}

if (exists("expr_mat")) {
  groups <- rep("Other", ncol(expr_mat))
  
  if ("source_name_ch1" %in% colnames(pdata_32707)) {
    src <- tolower(pdata_32707$source_name_ch1)
    groups[grepl("untreated|control|normal|healthy", src)] <- "Control"
    groups[grepl("ards", src)] <- "ARDS"
    groups[grepl("sepsis", src) & !grepl("ards", src)] <- "Sepsis"
    groups[grepl("sirs", src)] <- "SIRS"
  }
  
  if (sum(groups == "Other") > length(groups) / 2 && "description" %in% colnames(pdata_32707)) {
    desc <- tolower(pdata_32707$description)
    groups[grepl("control|normal|healthy", desc)] <- "Control"
    groups[grepl("ards", desc)] <- "ARDS"
    groups[grepl("sepsis", desc)] <- "Sepsis"
    groups[grepl("sirs", desc)] <- "SIRS"
  }
  
  names(groups) <- colnames(expr_mat)
  cat("样本分组统计:\n")
  print(table(groups))
  
  key_genes <- unique(unlist(ferroptosis_gene_sets))
  key_genes_found <- intersect(key_genes, rownames(expr_mat))
  cat("找到关键基因:", length(key_genes_found), "/", length(key_genes), "\n")
  
  if (length(key_genes_found) > 20) {
    group_order <- order(factor(groups, levels = c("Control", "SIRS", "Sepsis", "ARDS")))
    expr_sorted <- expr_mat[key_genes_found, group_order]
    groups_sorted <- groups[group_order]
    
    expr_scaled <- t(scale(t(expr_sorted)))
    
    ann_col <- data.frame(Group = groups_sorted)
    rownames(ann_col) <- colnames(expr_scaled)
    
    key_genes_unique <- unique(key_genes_found)
    gene_pathway <- data.frame(Pathway = rep("Other", length(key_genes_unique)),
                               row.names = key_genes_unique,
                               stringsAsFactors = FALSE)
    for (pathway_name in names(ferroptosis_gene_sets)) {
      for (g in ferroptosis_gene_sets[[pathway_name]]) {
        if (g %in% rownames(gene_pathway) && gene_pathway[g, "Pathway"] == "Other") {
          gene_pathway[g, "Pathway"] <- pathway_name
        }
      }
    }
    
    gene_var <- apply(expr_scaled[key_genes_unique, ], 1, var, na.rm = TRUE)
    top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(60, length(gene_var))]
    
    pheatmap(expr_scaled[top_genes, ],
             annotation_col = ann_col,
             annotation_row = gene_pathway[top_genes, , drop = FALSE],
             show_colnames = FALSE,
             show_rownames = TRUE,
             fontsize_row = 8,
             scale = "none",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "euclidean",
             clustering_method = "ward.D2",
             main = "Ferroptosis-Related Gene Expression Profile",
             filename = "results/13_extended_enrichment/03_ferroptosis_gene_heatmap.png",
             width = 12, height = 14,
             res = 150)
    cat("铁死亡基因热图已保存\n")
  }
}

# === 3. ssGSEA评分箱线图（多基因集） ===
cat("\n=== 3. 多基因集ssGSEA评分 ===\n")

# 手动ssGSEA函数
ssgsea_score <- function(expr_mat, gene_set) {
  common_genes <- intersect(gene_set, rownames(expr_mat))
  if (length(common_genes) < 3) return(rep(NA, ncol(expr_mat)))
  
  ranks <- apply(expr_mat, 2, function(x) rank(x, na.last = "keep"))
  rownames(ranks) <- rownames(expr_mat)
  
  n_genes <- nrow(expr_mat)
  n_set <- length(common_genes)
  
  scores <- numeric(ncol(expr_mat))
  for (j in 1:ncol(expr_mat)) {
    rank_order <- order(ranks[, j], decreasing = TRUE, na.last = NA)
    rank_order <- rank_order[!is.na(rank_order)]
    
    hit <- numeric(n_genes)
    hit[match(common_genes, rownames(expr_mat)[rank_order])] <- 1 / n_set
    miss <- (1 - hit) / (n_genes - n_set)
    
    es <- cumsum(hit - miss)
    scores[j] <- max(es) - min(es)
  }
  
  return(scores)
}

if (exists("expr_mat") && exists("groups")) {
  
  ssgsea_scores <- list()
  pathway_names <- c("Ferroptosis Regulators", "Iron Metabolism", "Glutathione System",
                     "Nrf2 Antioxidant Pathway", "Inflammatory Response", "Ferritinophagy",
                     "Autophagy", "Hypoxia Response")
  
  for (pn in pathway_names) {
    scores <- ssgsea_score(expr_mat, ferroptosis_gene_sets[[pn]])
    if (!all(is.na(scores))) {
      ssgsea_scores[[pn]] <- scores
    }
  }
  
  if (length(ssgsea_scores) > 0) {
    score_df <- data.frame(
      Sample = rep(colnames(expr_mat), length(ssgsea_scores)),
      Pathway = rep(names(ssgsea_scores), each = ncol(expr_mat)),
      Score = unlist(ssgsea_scores),
      Group = rep(groups, length(ssgsea_scores)),
      stringsAsFactors = FALSE
    )
    
    p3 <- ggplot(score_df, aes(x = Group, y = Score, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
      facet_wrap(~Pathway, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = c("Control" = "#2ECC71", "SIRS" = "#F1C40F",
                                    "Sepsis" = "#E67E22", "ARDS" = "#E74C3C")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8),
            strip.text = element_text(size = 9, face = "bold"),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(x = "", y = "ssGSEA Score", title = "ssGSEA Scores of Ferroptosis-Related Pathways")
    
    ggsave("results/13_extended_enrichment/04_ssgsea_pathways.png", p3, width = 16, height = 12, dpi = 150)
    cat("多通路ssGSEA箱线图已保存\n")
    
    cat("\n各组评分比较（Control vs ARDS）:\n")
    for (pn in names(ssgsea_scores)) {
      ctrl_scores <- ssgsea_scores[[pn]][groups == "Control"]
      ards_scores <- ssgsea_scores[[pn]][groups == "ARDS"]
      p_val <- wilcox.test(ctrl_scores, ards_scores)$p.value
      cat(sprintf("  %s: p=%.4f\n", pn, p_val))
    }
  }
}

# === 4. 铁死亡评分与免疫评分相关性 ===
cat("\n=== 4. 铁死亡与免疫通路相关性 ===\n")

if (length(ssgsea_scores) > 2) {
  score_matrix <- do.call(cbind, ssgsea_scores)
  cor_matrix <- cor(score_matrix, method = "spearman", use = "pairwise.complete.obs")
  
  pheatmap(cor_matrix,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 8,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           main = "Correlation of Ferroptosis-Related Pathway Scores",
           filename = "results/13_extended_enrichment/05_pathway_correlation.png",
           width = 10, height = 9,
           res = 150)
  cat("通路相关性热图已保存\n")
}

# === 5. DEG火山图（标注铁死亡基因） ===
cat("\n=== 5. 铁死亡基因火山图标注 ===\n")

deg_df_plot <- as.data.frame(deg_all_gene)
deg_df_plot$Gene <- rownames(deg_df_plot)
deg_df_plot$logP <- -log10(deg_df_plot$adj.P.Val)

ferro_genes_all <- unique(unlist(ferroptosis_gene_sets))
deg_df_plot$IsFerroptosis <- ifelse(deg_df_plot$Gene %in% ferro_genes_all, "Ferroptosis", "Other")
deg_df_plot$Significance <- ifelse(deg_df_plot$adj.P.Val < 0.05 & abs(deg_df_plot$logFC) > 1,
                                    ifelse(deg_df_plot$logFC > 0, "Up", "Down"), "Not Sig")

key_ferro_genes <- c("GPX4", "SLC7A11", "ACSL4", "NFE2L2", "PTGS2", "FTH1", 
                     "HMOX1", "NQO1", "GCLC", "CHAC1", "BECN1", "MAP1LC3A")
deg_df_plot$Label <- ifelse(deg_df_plot$Gene %in% key_ferro_genes & 
                              deg_df_plot$Significance != "Not Sig", 
                            deg_df_plot$Gene, "")

p5 <- ggplot(deg_df_plot, aes(x = logFC, y = logP)) +
  geom_point(aes(color = IsFerroptosis, size = IsFerroptosis), alpha = 0.6) +
  scale_color_manual(values = c("Ferroptosis" = "#E74C3C", "Other" = "#95A5A6")) +
  scale_size_manual(values = c("Ferroptosis" = 2.5, "Other" = 1)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20,
                   color = "#2C3E50", fontface = "bold") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(x = "log2(Fold Change)", y = "-log10(adj.P.Val)",
       title = "Volcano Plot with Ferroptosis Genes Highlighted",
       color = "Gene Type", size = "Gene Type")

ggsave("results/13_extended_enrichment/06_volcano_ferroptosis.png", p5, width = 10, height = 8, dpi = 150)
cat("铁死亡火山图已保存\n")

# === 保存 ===
cat("\n=== 保存结果 ===\n")
save(gsea_ferro, ferroptosis_gene_sets, ssgsea_scores,
     file = "data/GEO/GSE32707_extended_enrichment.RData")

cat("\n=== 扩展富集分析完成 ===\n")
cat("\n下一步：运行 04_WGCNA/01_WGCNA.R 进行加权基因共表达网络分析\n")

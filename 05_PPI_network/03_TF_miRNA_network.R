#' ---
#' title: "转录因子调控网络分析"
#' description: "从PPI网络中识别转录因子（TF），构建TF-靶基因调控网络，重点分析Nrf2（NFE2L2）对铁死亡基因的调控作用，进行Fisher精确检验验证富集"
#' input: "data/GEO/GSE32707_ppi_results.RData"
#' output: "data/GEO/GSE32707_tf_results.RData, results/05_PPI/05_TF_regulatory_network.png, results/05_PPI/06_TF_ferroptosis_barplot.png, results/05_PPI/transcription_factor_summary.csv"
#' dependencies: "igraph, ggplot2, httr, jsonlite"
#' ---

# ============================================
# 03_TF调控网络分析
# 转录因子调控网络预测 + 铁死亡基因上游调控分析
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("igraph")) stop("需要igraph包")
library(igraph)

if (!require("ggplot2")) stop("需要ggplot2包")
library(ggplot2)

if (!require("httr")) install.packages("httr", repos = "https://cloud.r-project.org")
library(httr)

if (!require("jsonlite")) install.packages("jsonlite", repos = "https://cloud.r-project.org")
library(jsonlite)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载PPI结果 ----
cat("\n=== 加载PPI结果 ===\n")
load("data/GEO/GSE32707_ppi_results.RData")
cat("加载完成\n")

# 铁死亡基因列表
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

# ---- 3. 从已有PPI网络中识别转录因子 ----
cat("\n=== 从PPI网络中识别转录因子 ===\n")

ppi_genes_all <- V(g)$name
cat("PPI网络总基因数:", length(ppi_genes_all), "\n")

# ---- 4. 识别转录因子 ----
cat("\n=== 识别转录因子 ===\n")

# 转录因子列表（从TRRUST/Homo sapiens TFs整理的核心TF）
tf_list <- c(
  "NFE2L2","NFE2L1","NFE2L3","NRF1","ATF4","ATF6","XBP1","HIF1A","EPAS1",
  "NFKB1","RELA","RELB","NFKB2","STAT1","STAT2","STAT3","STAT5A","STAT6",
  "JUN","JUNB","JUND","FOS","FOSL1","FOSL2","SP1","SP3","KLF4","KLF2",
  "TP53","TP63","TP73","E2F1","E2F4","MYC","MYCN","MAX","MAD","USF1","USF2",
  "PPARA","PPARG","PPARD","RXRA","RXRB","RXRG","THRA","THRB",
  "CREB1","CREBBP","EP300","HDAC1","HDAC2","SIRT1","SIRT3",
  "FOXO1","FOXO3","FOXO4","FOXM1","GATA1","GATA2","GATA3","GATA4",
  "IRF1","IRF3","IRF5","IRF7","CEBPA","CEBPB","CEBPD",
  "MAFF","MAFG","MAFK","BACH1","BACH2"
)
tf_list <- unique(tf_list)

# 在PPI网络中找TF
tf_in_network <- intersect(ppi_genes_all, tf_list)
cat("PPI网络中的转录因子数:", length(tf_in_network), "\n")
if (length(tf_in_network) > 0) {
  cat("转录因子:", paste(sort(tf_in_network), collapse = ", "), "\n")
}

# 从igraph对象提取边
edges_df <- as_data_frame(g, what = "edges")
cat("PPI网络边数:", nrow(edges_df), "\n")

# 构建TF-靶基因边（TF在from列，target在to列）
tf_edges_from <- edges_df[edges_df$from %in% tf_in_network, ]
tf_edges_to <- edges_df[edges_df$to %in% tf_in_network, ]
# 交换to的列
tf_edges_to_swapped <- data.frame(
  from = tf_edges_to$to,
  to = tf_edges_to$from,
  score = tf_edges_to$score,
  stringsAsFactors = FALSE
)

tf_target_edges <- rbind(
  data.frame(TF = tf_edges_from$from, Target = tf_edges_from$to, 
             score = if("score" %in% colnames(tf_edges_from)) tf_edges_from$score else rep(NA, nrow(tf_edges_from)),
             stringsAsFactors = FALSE),
  data.frame(TF = tf_edges_to_swapped$from, Target = tf_edges_to_swapped$to,
             score = if("score" %in% colnames(tf_edges_to_swapped)) tf_edges_to_swapped$score else rep(NA, nrow(tf_edges_to_swapped)),
             stringsAsFactors = FALSE)
)
# 去重
tf_target_edges <- tf_target_edges[!duplicated(tf_target_edges[, c("TF", "Target")]), ]

cat("TF-靶基因互作数:", nrow(tf_target_edges), "\n")
cat("调控铁死亡基因的TF数:", length(unique(tf_target_edges$TF[tf_target_edges$Target %in% ferroptosis_genes])), "\n")

# 统计每个TF调控的铁死亡基因数
tf_ferro_count <- aggregate(Target ~ TF, 
                             data = tf_target_edges[tf_target_edges$Target %in% ferroptosis_genes, ],
                             FUN = length)
colnames(tf_ferro_count) <- c("TF", "Ferroptosis_Targets")
tf_ferro_count <- tf_ferro_count[order(-tf_ferro_count$Ferroptosis_Targets), ]
cat("\n调控铁死亡基因最多的TF:\n")
print(tf_ferro_count)

# 统计每个TF调控的总靶基因数
tf_total_count <- aggregate(Target ~ TF, data = tf_target_edges, FUN = length)
colnames(tf_total_count) <- c("TF", "Total_Targets")

# 合并
tf_summary <- merge(tf_total_count, tf_ferro_count, by = "TF", all.x = TRUE)
tf_summary$Ferroptosis_Targets[is.na(tf_summary$Ferroptosis_Targets)] <- 0
tf_summary <- tf_summary[order(-tf_summary$Ferroptosis_Targets, -tf_summary$Total_Targets), ]

dir.create("results/05_PPI", recursive = TRUE, showWarnings = FALSE)
write.csv(tf_summary, "results/05_PPI/transcription_factor_summary.csv", row.names = FALSE)
cat("\nTF汇总表已保存\n")

# ---- 5. TF-靶基因网络图 ----
cat("\n=== 绘制TF-靶基因网络图 ===\n")

if (nrow(tf_target_edges) > 0 && nrow(tf_ferro_count) > 0) {
  # 取Top 10 TF（按调控铁死亡基因数排序）
  top_tfs <- head(tf_ferro_count$TF, 10)
  
  # 筛选这些TF的边
  top_edges <- tf_target_edges[tf_target_edges$TF %in% top_tfs, ]
  cat("Top 10 TF的靶基因互作数:", nrow(top_edges), "\n")
  
  # 构建igraph
  g_tf <- graph_from_data_frame(top_edges, directed = TRUE)
  
  # 节点类型：TF vs 靶基因 vs 铁死亡基因
  V(g_tf)$type <- ifelse(V(g_tf)$name %in% tf_list, "TF", "Target")
  V(g_tf)$is_ferroptosis <- V(g_tf)$name %in% ferroptosis_genes
  
  # 节点颜色
  V(g_tf)$color <- ifelse(V(g_tf)$type == "TF", "#8E44AD",
                            ifelse(V(g_tf)$is_ferroptosis, "#C0392B", "#5DADE2"))
  
  # 节点大小：TF按调控的铁死亡基因数
  tf_size_map <- setNames(tf_ferro_count$Ferroptosis_Targets, tf_ferro_count$TF)
  V(g_tf)$size <- ifelse(V(g_tf)$type == "TF",
                          8 + tf_size_map[V(g_tf)$name] * 1.5,
                          5)
  
  # 标签：TF全部显示，铁死亡靶基因显示，其他不显示
  V(g_tf)$label <- ifelse(V(g_tf)$type == "TF" | V(g_tf)$is_ferroptosis,
                           V(g_tf)$name, NA)
  
  # 箭头方向：TF -> target
  E(g_tf)$arrow.size <- 0.3
  E(g_tf)$color <- rgb(0.6, 0.6, 0.6, 0.3)
  E(g_tf)$width <- 0.5
  
  png("results/05_PPI/05_TF_regulatory_network.png", width = 2200, height = 1800, res = 150)
  set.seed(42)
  plot(g_tf,
       layout = layout_with_fr(g_tf),
       vertex.size = V(g_tf)$size,
       vertex.color = V(g_tf)$color,
       vertex.label = V(g_tf)$label,
       vertex.label.cex = 0.7,
       vertex.label.font = 2,
       vertex.label.color = "black",
       vertex.frame.color = "white",
       edge.arrow.size = 0.3,
       edge.color = E(g_tf)$color,
       edge.width = E(g_tf)$width,
       main = "Transcription Factor Regulatory Network of Ferroptosis Genes")
  legend("topright",
         legend = c("Transcription Factor", "Ferroptosis Target", "Other Target"),
         fill = c("#8E44AD", "#C0392B", "#5DADE2"),
         cex = 0.9,
         bty = "n")
  dev.off()
  cat("TF调控网络图已保存\n")
}

# ---- 6. Nrf2靶基因富集验证 ----
cat("\n=== Nrf2靶基因验证 ===\n")

# 找NFE2L2（Nrf2）的靶基因
nrf2_targets <- tf_target_edges$Target[tf_target_edges$TF == "NFE2L2"]
cat("NFE2L2（Nrf2）靶基因数:", length(nrf2_targets), "\n")
cat("Nrf2靶基因:", paste(sort(nrf2_targets), collapse = ", "), "\n")

nrf2_ferro <- intersect(nrf2_targets, ferroptosis_genes)
cat("其中铁死亡基因数:", length(nrf2_ferro), "\n")
cat("Nrf2调控的铁死亡基因:", paste(sort(nrf2_ferro), collapse = ", "), "\n")

# Fisher精确检验：Nrf2靶基因中铁死亡基因是否富集
all_genes_network <- unique(c(tf_target_edges$TF, tf_target_edges$Target))
ferro_in_net <- intersect(all_genes_network, ferroptosis_genes)
nrf2_in_net <- intersect(nrf2_targets, all_genes_network)
non_nrf2_in_net <- setdiff(all_genes_network, nrf2_in_net)

fisher_mat <- matrix(c(
  length(intersect(nrf2_in_net, ferroptosis_genes)),
  length(setdiff(nrf2_in_net, ferroptosis_genes)),
  length(intersect(non_nrf2_in_net, ferroptosis_genes)),
  length(setdiff(non_nrf2_in_net, ferroptosis_genes))
), nrow = 2, byrow = TRUE)
rownames(fisher_mat) <- c("Nrf2_target", "Non_Nrf2_target")
colnames(fisher_mat) <- c("Ferroptosis", "Non_ferroptosis")
cat("\nFisher精确检验列联表:\n")
print(fisher_mat)
cat("Fisher检验p值:", fisher.test(fisher_mat)$p.value, "\n")

# ---- 7. TF富集柱状图 ----
cat("\n=== 绘制TF富集图 ===\n")

if (nrow(tf_ferro_count) >= 3) {
  top_n <- min(15, nrow(tf_ferro_count))
  top15_tf <- head(tf_ferro_count, top_n)
  top15_tf$TF <- factor(top15_tf$TF, levels = rev(top15_tf$TF))
  top15_tf$is_Nrf2 <- top15_tf$TF %in% c("NFE2L2")
  
  p <- ggplot(top15_tf, aes(x = TF, y = Ferroptosis_Targets, fill = is_Nrf2)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("#5DADE2", "#8E44AD"),
                      labels = c("Other TF", "NFE2L2 (Nrf2)"),
                      name = "") +
    theme_bw() +
    labs(title = "Top 15 TFs Regulating Ferroptosis Genes",
         x = "", y = "Number of Ferroptosis Target Genes") +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.text = element_text(size = 9),
          legend.position = "bottom")
  
  ggsave("results/05_PPI/06_TF_ferroptosis_barplot.png", p, width = 6, height = 5, dpi = 150)
  cat("TF富集柱状图已保存\n")
}

# ---- 8. 保存 ----
cat("\n=== 保存结果 ===\n")
save(tf_target_edges, tf_summary, tf_ferro_count, nrf2_targets, nrf2_ferro,
     file = "data/GEO/GSE32707_tf_results.RData")

cat("\n=== 转录因子分析全部完成 ===\n")
cat("产出文件:\n")
cat("  results/05_PPI/05_TF_regulatory_network.png - TF调控网络图\n")
cat("  results/05_PPI/06_TF_ferroptosis_barplot.png - 调控铁死亡的Top TF柱状图\n")
cat("  results/05_PPI/transcription_factor_summary.csv - TF汇总表\n")
cat("\n下一步：运行 06_ferroptosis_analysis/01_ferroptosis_score.R 进行铁死亡分析\n")

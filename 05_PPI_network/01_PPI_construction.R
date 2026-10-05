#' ---
#' title: "PPI蛋白互作网络构建与Hub基因筛选"
#' description: "基于STRING数据库构建turquoise模块中铁死亡相关基因的PPI网络，计算网络中心性（Degree/Betweenness/MCC）识别Hub基因，绘制PPI网络图和MM-GS散点图"
#' input: "data/GEO/GSE32707_wgcna.RData, data/GEO/GSE32707_gsea_results.RData"
#' output: "data/GEO/GSE32707_ppi_results.RData, results/05_PPI/ 下的4张图和Hub基因汇总表"
#' dependencies: "igraph, ggplot2, WGCNA, httr, jsonlite"
#' ---

# ============================================
# 01_PPI网络构建 + Hub基因筛选 + MM-GS散点图
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
if (!require("igraph")) install.packages("igraph", repos = "https://cloud.r-project.org")
library(igraph)

if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cloud.r-project.org")
library(ggplot2)

if (!require("WGCNA")) stop("需要WGCNA包")
library(WGCNA)

if (!require("httr")) install.packages("httr", repos = "https://cloud.r-project.org")
library(httr)

if (!require("jsonlite")) install.packages("jsonlite", repos = "https://cloud.r-project.org")
library(jsonlite)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载WGCNA和GSEA结果 ----
cat("\n=== 加载WGCNA结果 ===\n")
load("data/GEO/GSE32707_wgcna.RData")
cat("WGCNA变量加载完成\n")
cat("模块数:", length(unique(module_colors)), "\n")
cat("各模块基因数:\n")
print(table(module_colors))

cat("\n=== 加载GSEA结果（获取probe_map）===\n")
load("data/GEO/GSE32707_gsea_results.RData")
cat("probe_map维度:", dim(probe_map), "\n")

# ---- 3. 提取turquoise模块基因并映射到symbol ----
cat("\n=== 提取turquoise模块基因 ===\n")

turquoise_probes <- colnames(datExpr)[module_colors == "turquoise"]
cat("turquoise模块探针数:", length(turquoise_probes), "\n")

# 映射到gene symbol
turquoise_mapped <- probe_map[probe_map$probe %in% turquoise_probes, ]
turquoise_mapped <- turquoise_mapped[turquoise_mapped$symbol != "" & !is.na(turquoise_mapped$symbol), ]
turquoise_mapped <- turquoise_mapped[!duplicated(turquoise_mapped$symbol), ]
cat("映射到symbol后基因数:", nrow(turquoise_mapped), "\n")

# 取交集：turquoise模块 ∩ 差异基因
turquoise_genes <- turquoise_mapped$symbol
deg_genes <- rownames(deg_mapped)
turquoise_deg <- intersect(turquoise_genes, deg_genes)
cat("turquoise ∩ DEG基因数:", length(turquoise_deg), "\n")

# 铁死亡基因交集
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
turquoise_ferro <- intersect(turquoise_genes, ferroptosis_genes)
cat("turquoise ∩ 铁死亡基因数:", length(turquoise_ferro), "\n")
cat("铁死亡基因:", paste(turquoise_ferro, collapse = ", "), "\n")

# ---- 4. 筛选PPI核心基因 ----
cat("\n=== 筛选PPI核心基因 ===\n")

# 计算所有turquoise基因的MM和GS
turquoise_probes_all <- colnames(datExpr)[module_colors == "turquoise"]
ME_turquoise <- MEs[, "MEturquoise"]
MM_vals <- abs(as.numeric(cor(datExpr[, turquoise_probes_all], ME_turquoise, use = "p")))
GS_vals <- abs(as.numeric(cor(datExpr[, turquoise_probes_all], datTraits$ARDS, use = "p")))
MMGS <- MM_vals * GS_vals
names(MMGS) <- turquoise_probes_all
MMGS <- sort(MMGS, decreasing = TRUE)

# 映射到symbol
mmgs_df <- data.frame(
  probe = names(MMGS),
  MMGS = MMGS,
  stringsAsFactors = FALSE
)
mmgs_df$symbol <- probe_map$symbol[match(mmgs_df$probe, probe_map$probe)]
mmgs_df <- mmgs_df[!is.na(mmgs_df$symbol) & mmgs_df$symbol != "", ]
mmgs_df <- mmgs_df[!duplicated(mmgs_df$symbol), ]

# 取top 80个非铁死亡基因
top_mmgs_genes <- mmgs_df$symbol[!mmgs_df$symbol %in% turquoise_ferro][1:min(80, sum(!mmgs_df$symbol %in% turquoise_ferro))]

# PPI基因 = 铁死亡基因 + top MM*GS基因
ppi_genes <- unique(c(turquoise_ferro, top_mmgs_genes))
cat("用于PPI构建的基因数:", length(ppi_genes), "\n")
cat("  其中铁死亡基因:", length(intersect(ppi_genes, turquoise_ferro)), "\n")

# ---- 5. 从STRING获取PPI数据 ----
cat("\n=== 从STRING获取PPI数据 ===\n")

cat("查询STRING数据库（POST方法）...\n")
string_data <- NULL

tryCatch({
  batch_size <- 100
  all_edges <- list()
  
  for (batch_start in seq(1, length(ppi_genes), by = batch_size)) {
    batch_end <- min(batch_start + batch_size - 1, length(ppi_genes))
    batch_genes <- ppi_genes[batch_start:batch_end]
    
    cat("  批次", ceiling(batch_start/batch_size), ":", length(batch_genes), "个基因\n")
    
    body <- list(
      identifiers = paste(batch_genes, collapse = "\n"),
      species = "9606",
      required_score = 400,
      caller_identity = "ALI_research"
    )
    
    response <- POST("https://string-db.org/api/json/network",
                     body = body,
                     encode = "form",
                     timeout(120))
    
    if (!http_error(response)) {
      batch_data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
      if (nrow(batch_data) > 0) {
        all_edges[[length(all_edges) + 1]] <- batch_data
      }
    }
  }
  
  if (length(all_edges) > 0) {
    string_data <- do.call(rbind, all_edges)
    string_data <- string_data[!duplicated(string_data[, c("preferredName_A", "preferredName_B")]), ]
    cat("STRING返回交互数:", nrow(string_data), "\n")
    cat("STRING网络基因数:", length(unique(c(string_data$preferredName_A, string_data$preferredName_B))), "\n")
  } else {
    cat("STRING未返回数据\n")
    string_data <- NULL
  }
}, error = function(e) {
  cat("STRING API错误:", e$message, "\n")
  string_data <- NULL
})

# 如果STRING API失败，使用MM*GS方法找hub基因
if (is.null(string_data) || nrow(string_data) == 0) {
  cat("STRING数据不足，使用MM×GS方法识别Hub基因\n")
  stop("STRING API不可用，请检查网络连接")
} else {
  # 用STRING数据构建网络
  edges <- data.frame(
    from = string_data$preferredName_A,
    to = string_data$preferredName_B,
    score = string_data$score,
    stringsAsFactors = FALSE
  )
  # 去掉自环
  edges <- edges[edges$from != edges$to, ]
  cat("有效交互边数:", nrow(edges), "\n")
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 简化网络：去掉孤立节点
  g <- delete_vertices(g, degree(g) == 0)
  cat("PPI网络节点数:", vcount(g), "\n")
  cat("PPI网络边数:", ecount(g), "\n")
}

# ---- 6. 计算网络中心性指标 ----
cat("\n=== 计算网络中心性 ===\n")

# Degree
deg_cent <- degree(g)
# Betweenness
btw_cent <- betweenness(g, normalized = TRUE)
# Closeness
close_cent <- closeness(g, normalized = TRUE)

# MCC (Maximum Clique Centrality)
cat("计算MCC...\n")
mcc_scores <- sapply(V(g), function(v) {
  neighbors_v <- neighbors(g, v)
  if (length(neighbors_v) < 2) return(1)
  subg <- induced_subgraph(g, c(v, neighbors_v))
  if (vcount(subg) < 3) return(length(neighbors_v))
  cliques_v <- max_cliques(subg, min = 2)
  if (length(cliques_v) == 0) return(length(neighbors_v))
  max(lengths(cliques_v) - 1)
})

# 汇总
hub_df <- data.frame(
  Gene = V(g)$name,
  Degree = deg_cent,
  Betweenness = btw_cent,
  Closeness = close_cent,
  MCC = mcc_scores,
  stringsAsFactors = FALSE
)
hub_df <- hub_df[order(-hub_df$Degree), ]
cat("\nTop 20 Hub基因（按Degree排序）:\n")
print(head(hub_df, 20))

dir.create("results/05_PPI", recursive = TRUE, showWarnings = FALSE)
write.csv(hub_df, "results/05_PPI/hub_genes_all.csv", row.names = FALSE)

# 取Top 10
top10_hub <- head(hub_df$Gene, 10)
cat("\nTop 10 Hub基因:", paste(top10_hub, collapse = ", "), "\n")

# 标注铁死亡基因
hub_df$is_ferroptosis <- hub_df$Gene %in% ferroptosis_genes
cat("Top 20中铁死亡基因数:", sum(head(hub_df, 20)$is_ferroptosis), "\n")

# ---- 7. PPI网络可视化 ----
cat("\n=== PPI网络可视化 ===\n")

# 只取Top 100节点可视化
if (vcount(g) > 100) {
  top_nodes <- head(hub_df$Gene, 100)
  g_sub <- induced_subgraph(g, V(g)[name %in% top_nodes])
} else {
  g_sub <- g
}

# 设置节点属性
deg_sub <- degree(g_sub)
V(g_sub)$size <- 6 + (deg_sub / max(deg_sub)) * 14
V(g_sub)$color <- ifelse(V(g_sub)$name %in% ferroptosis_genes, "#C0392B", "#5DADE2")

# 标注：铁死亡基因 + Top 10 hub
V(g_sub)$label <- ifelse(V(g_sub)$name %in% c(top10_hub, turquoise_ferro), V(g_sub)$name, NA)

png("results/05_PPI/01_PPI_network.png", width = 2000, height = 2000, res = 150)
set.seed(42)
plot(g_sub,
     layout = layout_with_fr(g_sub),
     vertex.size = V(g_sub)$size,
     vertex.color = V(g_sub)$color,
     vertex.label = V(g_sub)$label,
     vertex.label.cex = 0.75,
     vertex.label.color = "black",
     vertex.label.font = 2,
     vertex.frame.color = "white",
     edge.color = rgb(0.6, 0.6, 0.6, 0.4),
     edge.width = 0.6,
     main = "PPI Network of Ferroptosis-related Genes in Turquoise Module")
legend("topright",
       legend = c("Ferroptosis gene", "Other gene"),
       fill = c("#C0392B", "#5DADE2"),
       cex = 0.9,
       bty = "n")
dev.off()
cat("PPI网络图已保存\n")

# ---- 8. Hub基因排序柱状图 ----
top20_df <- head(hub_df, 20)
top20_df$Gene <- factor(top20_df$Gene, levels = rev(top20_df$Gene))

p <- ggplot(top20_df, aes(x = Gene, y = Degree, fill = is_ferroptosis)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("#5DADE2", "#C0392B"),
                    labels = c("Non-ferroptosis", "Ferroptosis"),
                    name = "") +
  theme_bw() +
  labs(title = "Top 20 Hub Genes by Degree",
       x = "", y = "Degree") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 9))

ggsave("results/05_PPI/02_hub_genes_barplot.png", p, width = 6, height = 5, dpi = 150)
cat("Hub基因柱状图已保存\n")

# ---- 9. 三种算法比较图 ----
hub_df$Degree_scaled <- scale(hub_df$Degree)[,1]
hub_df$Betweenness_scaled <- scale(hub_df$Betweenness)[,1]
hub_df$MCC_scaled <- scale(hub_df$MCC)[,1]

top30_df <- head(hub_df, 30)

p2 <- ggplot(top30_df, aes(x = Degree, y = Betweenness, color = MCC, size = MCC)) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = ifelse(Gene %in% top10_hub, Gene, "")),
            hjust = -0.1, vjust = 0.5, size = 3) +
  scale_color_gradient(low = "#5DADE2", high = "#C0392B") +
  theme_bw() +
  labs(title = "Hub Gene Comparison: Degree vs Betweenness vs MCC",
       x = "Degree", y = "Betweenness (normalized)") +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

ggsave("results/05_PPI/03_hub_comparison.png", p2, width = 7, height = 5, dpi = 150)
cat("Hub基因比较图已保存\n")

# ---- 10. MM-GS散点图 ----
cat("\n=== MM-GS散点图 ===\n")

turquoise_probes_all <- colnames(datExpr)[module_colors == "turquoise"]

if (length(turquoise_probes_all) > 0) {
  ME_turquoise <- MEs[, "MEturquoise"]
  
  MM <- cor(datExpr[, turquoise_probes_all], ME_turquoise, use = "p")
  GS_ards <- cor(datExpr[, turquoise_probes_all], datTraits$ARDS, use = "p")
  
  mm_gs_df <- data.frame(
    probe = turquoise_probes_all,
    MM = as.numeric(MM),
    GS = as.numeric(GS_ards),
    stringsAsFactors = FALSE
  )
  mm_gs_df$symbol <- probe_map$symbol[match(mm_gs_df$probe, probe_map$probe)]
  mm_gs_df$is_ferroptosis <- mm_gs_df$symbol %in% ferroptosis_genes
  
  mm_gs_df <- mm_gs_df[!is.na(mm_gs_df$symbol) & mm_gs_df$symbol != "", ]
  
  p3 <- ggplot(mm_gs_df, aes(x = MM, y = GS, color = is_ferroptosis)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "grey50", linetype = "dashed") +
    geom_text(data = mm_gs_df[mm_gs_df$is_ferroptosis, ],
              aes(label = symbol), vjust = -0.8, size = 2.5, color = "#C0392B") +
    scale_color_manual(values = c("#5DADE2", "#C0392B"),
                       labels = c("Other", "Ferroptosis"),
                       name = "") +
    theme_bw() +
    labs(title = "Module Membership vs Gene Significance (Turquoise)",
         x = "Module Membership (MM)", y = "Gene Significance (GS, ARDS)") +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          legend.position = "bottom")
  
  ggsave("results/05_PPI/04_MM_GS_scatter.png", p3, width = 6, height = 5, dpi = 150)
  cat("MM-GS散点图已保存\n")
  
  cor_mm_gs <- cor(mm_gs_df$MM, mm_gs_df$GS, use = "p")
  cat("MM-GS相关性:", round(cor_mm_gs, 3), "\n")
}

# ---- 11. 保存所有结果 ----
cat("\n=== 保存结果 ===\n")
save(g, hub_df, top10_hub, ppi_genes, turquoise_genes, turquoise_ferro,
     file = "data/GEO/GSE32707_ppi_results.RData")

write.csv(hub_df, "results/05_PPI/hub_genes_summary.csv", row.names = FALSE)

cat("\n=== PPI + Hub基因分析全部完成 ===\n")
cat("产出文件:\n")
cat("  results/05_PPI/01_PPI_network.png - PPI网络图\n")
cat("  results/05_PPI/02_hub_genes_barplot.png - Top20 Hub基因柱状图\n")
cat("  results/05_PPI/03_hub_comparison.png - Hub基因三算法比较\n")
cat("  results/05_PPI/04_MM_GS_scatter.png - MM-GS散点图\n")
cat("  results/05_PPI/hub_genes_summary.csv - Hub基因汇总表\n")
cat("\n下一步：运行 05_PPI_network/02_MCODE_modules.R 进行模块分析\n")

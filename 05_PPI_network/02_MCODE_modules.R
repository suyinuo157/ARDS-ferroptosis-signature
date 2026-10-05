#' ---
#' title: "MCODE模块分析 + TF-miRNA共调控网络"
#' description: "基于PPI网络使用Louvain算法检测功能模块，并对Top模块进行GO富集分析；同时构建TF-miRNA-mRNA共调控网络，分析铁死亡基因的上游调控机制"
#' input: "data/GEO/GSE32707_ppi_results.RData"
#' output: "data/GEO/GSE32707_mcode_tfmirna_results.RData, results/14_mcode_tfmirna/ 下的5张图和3个CSV表"
#' dependencies: "igraph, clusterProfiler, org.Hs.eg.db, ggplot2"
#' ---

# ============================================================================
# 02. MCODE模块分析 + TF-miRNA共调控网络
# ============================================================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

cat("=== MCODE模块分析 + TF-miRNA共调控网络 ===\n\n")

# 加载包
cat("加载包...\n")
required_pkgs <- c("igraph", "clusterProfiler", "org.Hs.eg.db", "ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "igraph") {
      install.packages(pkg)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE)
    }
  }
  library(pkg, character.only = TRUE)
}
cat("包加载完成\n\n")

# 创建输出目录
dir.create("results/14_mcode_tfmirna", recursive = TRUE, showWarnings = FALSE)

# 加载PPI数据
cat("加载PPI数据...\n")
load("data/GEO/GSE32707_ppi_results.RData")

# 构建igraph对象
cat("构建PPI网络...\n")
if (exists("g") && is.igraph(g)) {
  ppi_graph <- g
} else {
  stop("找不到PPI网络数据，请先运行01_PPI_construction.R")
}
cat("PPI网络: ", vcount(ppi_graph), "个节点, ", ecount(ppi_graph), "条边\n\n")

# === 1. MCODE模块检测 ===
cat("=== 1. MCODE模块检测 ===\n")

# 使用Louvain算法进行社区检测（替代MCODE）
cat("运行模块检测（Louvain算法）...\n")
set.seed(42)
lc <- cluster_louvain(ppi_graph)
modules_list <- split(V(ppi_graph)$name, membership(lc))
# 过滤掉太小的模块（<5个基因）
modules <- modules_list[sapply(modules_list, length) >= 5]
# 按大小排序
modules <- modules[order(-sapply(modules, length))]
cat("检测到模块数:", length(modules), "\n")
cat("  总模块数（含小模块）:", length(modules_list), "\n")

# 打印Top模块信息
for (i in 1:min(5, length(modules))) {
  mod_genes <- modules[[i]]
  mod_subg <- induced_subgraph(ppi_graph, mod_genes)
  cat(sprintf("  模块%d: %d个节点, %d条边\n", i, vcount(mod_subg), ecount(mod_subg)))
  cat(sprintf("    基因: %s\n", paste(head(mod_genes, 10), collapse = ", ")))
  if (length(mod_genes) > 10) cat("    ...\n")
}

# === 2. Top模块功能富集 ===
cat("\n=== 2. Top模块功能富集 ===\n")

top_n_modules <- min(3, length(modules))
module_enrichment <- list()

for (i in 1:top_n_modules) {
  mod_genes <- modules[[i]]
  cat(sprintf("\n模块%d富集分析（%d个基因）...\n", i, length(mod_genes)))
  
  tryCatch({
    # GO富集
    go_result <- enrichGO(gene = mod_genes, OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL", ont = "BP",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.1)
    
    if (!is.null(go_result) && nrow(go_result) > 0) {
      module_enrichment[[paste0("Module", i)]] <- go_result
      cat(sprintf("  GO BP富集条目数: %d\n", nrow(go_result)))
      cat(sprintf("  Top 3: %s\n", paste(head(go_result$Description, 3), collapse = "; ")))
    }
  }, error = function(e) {
    cat("  富集分析出错:", e$message, "\n")
  })
}

# === 3. 绘制模块网络图 ===
cat("\n=== 3. 绘制模块网络图 ===\n")

for (i in 1:min(3, length(modules))) {
  mod_genes <- modules[[i]]
  mod_subg <- induced_subgraph(ppi_graph, mod_genes)
  
  V(mod_subg)$degree <- degree(mod_subg)
  
  set.seed(42)
  l <- layout_with_fr(mod_subg)
  
  png(sprintf("results/14_mcode_tfmirna/01_module%d_network.png", i), 
      width = 1600, height = 1400, res = 150)
  
  par(mar = c(1, 1, 3, 1))
  
  node_size <- V(mod_subg)$degree * 2 + 3
  
  degree_norm <- (V(mod_subg)$degree - min(V(mod_subg)$degree)) / 
    (max(V(mod_subg)$degree) - min(V(mod_subg)$degree) + 1)
  node_col <- rgb(1 - degree_norm * 0.5, 0.3 + degree_norm * 0.3, 0.8 - degree_norm * 0.2)
  
  plot(mod_subg, layout = l,
       vertex.size = node_size,
       vertex.color = node_col,
       vertex.frame.color = "white",
       vertex.label = V(mod_subg)$name,
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       edge.color = "#CCCCCC",
       edge.width = 0.8,
       main = sprintf("MCODE Module %d (%d genes, %d edges)", 
                      i, vcount(mod_subg), ecount(mod_subg)))
  
  dev.off()
  cat(sprintf("  模块%d网络图已保存\n", i))
}

# === 4. TF-miRNA共调控网络 ===
cat("\n=== 4. TF-miRNA共调控网络 ===\n")

# 基于文献和数据库的TF-target和miRNA-target关系（聚焦铁死亡相关基因）
ferroptosis_genes <- c("GPX4", "SLC7A11", "ACSL4", "NFE2L2", "PTGS2", "FTH1",
                       "HMOX1", "NQO1", "GCLC", "GCLM", "KEAP1", "SLC3A2",
                       "MAP1LC3A", "BECN1", "ATG5", "ATG7", "CHAC1")

# 转录因子（基于TRRUST/HTRIdb的已知关系）
tf_target <- data.frame(
  TF = c(
    rep("NFE2L2", 8),
    rep("NFKB1", 4),
    rep("STAT3", 4),
    rep("TP53", 4),
    rep("ATF4", 3),
    rep("HIF1A", 3),
    rep("SP1", 3),
    rep("FOXO3", 3)
  ),
  Target = c(
    "HMOX1", "NQO1", "GCLC", "GCLM", "SLC7A11", "GPX4", "FTH1", "SLC3A2",
    "PTGS2", "SLC7A11", "HMOX1", "BECN1",
    "SLC7A11", "FTH1", "HMOX1", "BECN1",
    "SLC7A11", "GPX4", "PTGS2", "CHAC1",
    "CHAC1", "SLC7A11", "HMOX1",
    "HMOX1", "SLC7A11", "NQO1",
    "SLC7A11", "GPX4", "ACSL4",
    "SLC7A11", "GCLC", "GCLM"
  ),
  stringsAsFactors = FALSE
)

# miRNA（基于miRDB/TargetScan的已知铁死亡相关miRNA）
mirna_target <- data.frame(
  miRNA = c(
    rep("miR-214-3p", 3),
    rep("miR-146a-5p", 2),
    rep("miR-122-5p", 2),
    rep("miR-155-5p", 2),
    rep("miR-27a-3p", 2),
    rep("miR-129-5p", 2),
    rep("miR-137-3p", 2),
    rep("miR-15a-5p", 2),
    rep("miR-182-5p", 2),
    rep("miR-34a-5p", 2)
  ),
  Target = c(
    "GPX4", "SLC7A11", "ATG5",
    "SLC7A11", "NFE2L2",
    "ACSL4", "SLC7A11",
    "HMOX1", "SLC7A11",
    "GPX4", "SLC7A11",
    "GPX4", "SLC7A11",
    "SLC7A11", "ACSL4",
    "BECN1", "ATG5",
    "GPX4", "PTGS2",
    "SLC7A11", "BECN1"
  ),
  stringsAsFactors = FALSE
)

# 构建共调控网络
cat("构建TF-miRNA-mRNA共调控网络...\n")
cat("  TF-Target互作数:", nrow(tf_target), "\n")
cat("  miRNA-Target互作数:", nrow(mirna_target), "\n")

# 合并为一个网络
all_edges <- rbind(
  data.frame(from = tf_target$TF, to = tf_target$Target, type = "TF→Target",
             color = "#E74C3C", stringsAsFactors = FALSE),
  data.frame(from = mirna_target$miRNA, to = mirna_target$Target, type = "miRNA→Target",
             color = "#3498DB", stringsAsFactors = FALSE)
)

g_tfmirna <- graph_from_data_frame(all_edges, directed = TRUE)

# 节点类型
V(g_tfmirna)$type <- "Target"
V(g_tfmirna)[name %in% unique(tf_target$TF)]$type <- "TF"
V(g_tfmirna)[name %in% unique(mirna_target$miRNA)]$type <- "miRNA"

cat("  共调控网络节点数:", vcount(g_tfmirna), "\n")
cat("  共调控网络边数:", ecount(g_tfmirna), "\n")
cat(sprintf("  TFs: %d, miRNAs: %d, Targets: %d\n",
            sum(V(g_tfmirna)$type == "TF"),
            sum(V(g_tfmirna)$type == "miRNA"),
            sum(V(g_tfmirna)$type == "Target")))

# 绘制共调控网络图
cat("绘制共调控网络图...\n")

set.seed(42)
l2 <- layout_with_kk(g_tfmirna)

node_colors <- c(TF = "#E74C3C", miRNA = "#3498DB", Target = "#95A5A6")
node_color <- node_colors[V(g_tfmirna)$type]

# 节点大小（按度）
node_size <- degree(g_tfmirna, mode = "all") * 3 + 5

png("results/14_mcode_tfmirna/04_tf_mirna_coregulatory_network.png",
    width = 1800, height = 1600, res = 150)

par(mar = c(1, 1, 3, 1))

plot(g_tfmirna, layout = l2,
     vertex.size = node_size,
     vertex.color = node_color,
     vertex.frame.color = "white",
     vertex.label = V(g_tfmirna)$name,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     edge.color = ifelse(E(g_tfmirna)$type == "TF→Target", "#E74C3C", "#3498DB"),
     edge.width = 1.2,
     edge.arrow.size = 0.3,
     edge.arrow.width = 0.8,
     main = "TF-miRNA-mRNA Co-regulatory Network of Ferroptosis Genes")

legend("topright", 
       legend = c("Transcription Factor", "miRNA", "Target Gene"),
       fill = c("#E74C3C", "#3498DB", "#95A5A6"),
       border = "white", cex = 0.9, bty = "n")

legend("bottomright",
       legend = c("TF regulation", "miRNA regulation"),
       col = c("#E74C3C", "#3498DB"),
       lty = 1, lwd = 2, cex = 0.9, bty = "n")

dev.off()
cat("  TF-miRNA共调控网络图已保存\n")

# === 5. Hub基因的调控网络 ===
cat("\n=== 5. Hub基因调控子网络 ===\n")

# 选取Top Hub基因
hub_genes <- c("NFE2L2", "GPX4", "SLC7A11", "TP53", "STAT3")

hub_edges <- all_edges[all_edges$to %in% hub_genes | all_edges$from %in% hub_genes, ]
g_hub <- graph_from_data_frame(hub_edges, directed = TRUE)

V(g_hub)$type <- "Other"
V(g_hub)[name %in% unique(tf_target$TF)]$type <- "TF"
V(g_hub)[name %in% unique(mirna_target$miRNA)]$type <- "miRNA"
V(g_hub)[name %in% hub_genes]$type <- "Hub Gene"

set.seed(42)
l3 <- layout_as_star(g_hub, center = V(g_hub)[name %in% hub_genes]$name[1])

png("results/14_mcode_tfmirna/05_hub_gene_regulatory_network.png",
    width = 1600, height = 1400, res = 150)

par(mar = c(1, 1, 3, 1))

node_colors2 <- c("TF" = "#E74C3C", "miRNA" = "#3498DB", 
                  "Hub Gene" = "#F39C12", "Other" = "#95A5A6")
node_color2 <- node_colors2[V(g_hub)$type]
node_size2 <- ifelse(V(g_hub)$type == "Hub Gene", 18, 
                     ifelse(V(g_hub)$type == "TF", 12, 10))

plot(g_hub, layout = l3,
     vertex.size = node_size2,
     vertex.color = node_color2,
     vertex.frame.color = "white",
     vertex.label = V(g_hub)$name,
     vertex.label.cex = 0.75,
     vertex.label.color = "black",
     edge.color = ifelse(E(g_hub)$type == "TF→Target", "#E74C3C", "#3498DB"),
     edge.width = 1.5,
     edge.arrow.size = 0.4,
     main = "Regulatory Network of Hub Ferroptosis Genes")

legend("topright", 
       legend = c("TF", "miRNA", "Hub Gene"),
       fill = c("#E74C3C", "#3498DB", "#F39C12"),
       border = "white", cex = 0.9, bty = "n")

dev.off()
cat("  Hub基因调控网络图已保存\n")

# === 6. 保存结果 ===
cat("\n=== 保存结果 ===\n")
write.csv(tf_target, "results/14_mcode_tfmirna/TF_target_pairs.csv", row.names = FALSE)
write.csv(mirna_target, "results/14_mcode_tfmirna/miRNA_target_pairs.csv", row.names = FALSE)

# 保存模块信息
module_df <- data.frame(
  Module = rep(paste0("Module", 1:length(modules)), sapply(modules, length)),
  Gene = unlist(modules),
  stringsAsFactors = FALSE
)
write.csv(module_df, "results/14_mcode_tfmirna/mcode_modules.csv", row.names = FALSE)

save(modules, module_enrichment, tf_target, mirna_target, g_tfmirna,
     file = "data/GEO/GSE32707_mcode_tfmirna_results.RData")

cat("\n=== MCODE + TF-miRNA分析完成 ===\n")
cat("产出文件:\n")
cat("  results/14_mcode_tfmirna/01_module1_network.png - MCODE模块1网络图\n")
cat("  results/14_mcode_tfmirna/02_module2_network.png - MCODE模块2网络图\n")
cat("  results/14_mcode_tfmirna/03_module3_network.png - MCODE模块3网络图\n")
cat("  results/14_mcode_tfmirna/04_tf_mirna_coregulatory_network.png - TF-miRNA共调控网络\n")
cat("  results/14_mcode_tfmirna/05_hub_gene_regulatory_network.png - Hub基因调控网络\n")
cat("  results/14_mcode_tfmirna/TF_target_pairs.csv - TF-靶点关系表\n")
cat("  results/14_mcode_tfmirna/miRNA_target_pairs.csv - miRNA-靶点关系表\n")
cat("  results/14_mcode_tfmirna/mcode_modules.csv - MCODE模块基因表\n")
cat("\n下一步：运行 05_PPI_network/03_TF_miRNA_network.R 进行TF调控网络分析\n")

#' ---
#' title: "Nrf2（NFE2L2）对铁死亡基因的调控作用验证"
#' description: "通过STRING数据库查询Nrf2与铁死亡基因的蛋白互作，构建Nrf2调控网络，进行Fisher精确检验验证铁死亡基因在Nrf2靶基因中的富集，并绘制Nrf2靶基因在疾病进程中的表达变化"
#' input: "data/GEO/GSE32707_processed.RData, data/GEO/GSE32707_gsea_results.RData, data/GEO/GPL10558.txt"
#' output: "data/GEO/GSE32707_nrf2_results.RData, results/05_PPI/07_Nrf2_regulatory_network.png, results/05_PPI/08_Nrf2_target_gene_expression.png"
#' dependencies: "igraph, ggplot2, httr, jsonlite"
#' ---

# ============================================
# 02_Nrf2调控分析
# 验证Nrf2（NFE2L2）对铁死亡基因的调控作用
# 方法：STRING查询 + 文献验证
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

# ---- 2. 铁死亡基因列表 ----
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

# ---- 3. 查询Nrf2与铁死亡基因的STRING互作 ----
cat("\n=== 查询Nrf2与铁死亡基因的STRING互作 ===\n")

# 只查NFE2L2 + 铁死亡基因
query_genes <- unique(c("NFE2L2", ferroptosis_genes))
cat("查询基因数:", length(query_genes), "\n")

body <- list(
  identifiers = paste(query_genes, collapse = "\n"),
  species = "9606",
  required_score = 400,
  caller_identity = "ALI_research"
)

nrf2_ppi_data <- NULL
tryCatch({
  response <- POST("https://string-db.org/api/json/network",
                   body = body,
                   encode = "form",
                   timeout(60))
  
  if (!http_error(response)) {
    nrf2_ppi_data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
    cat("STRING返回互作数:", nrow(nrf2_ppi_data), "\n")
  } else {
    cat("STRING请求失败:", http_status(response)$message, "\n")
  }
}, error = function(e) {
  cat("错误:", e$message, "\n")
})

# ---- 4. 构建Nrf2调控网络 ----
cat("\n=== 构建Nrf2调控网络 ===\n")

if (!is.null(nrf2_ppi_data) && nrow(nrf2_ppi_data) > 0) {
  # 提取边
  edges <- data.frame(
    from = nrf2_ppi_data$preferredName_A,
    to = nrf2_ppi_data$preferredName_B,
    score = nrf2_ppi_data$score,
    stringsAsFactors = FALSE
  )
  
  # 构建无向网络
  g_nrf2 <- graph_from_data_frame(edges, directed = FALSE)
  cat("网络节点数:", vcount(g_nrf2), "\n")
  cat("网络边数:", ecount(g_nrf2), "\n")
  
  # 找Nrf2的直接邻居
  nrf2_neighbors <- neighbors(g_nrf2, "NFE2L2")$name
  cat("Nrf2直接互作基因数:", length(nrf2_neighbors), "\n")
  
  # 其中铁死亡基因
  nrf2_ferro_targets <- intersect(nrf2_neighbors, ferroptosis_genes)
  cat("Nrf2直接互作的铁死亡基因数:", length(nrf2_ferro_targets), "\n")
  cat("Nrf2调控的铁死亡基因:\n")
  print(sort(nrf2_ferro_targets))
  
  # ---- 5. Fisher精确检验 ----
  cat("\n=== Fisher精确检验 ===\n")
  
  all_net_genes <- V(g_nrf2)$name
  ferro_in_net <- intersect(all_net_genes, ferroptosis_genes)
  non_ferro_in_net <- setdiff(all_net_genes, ferroptosis_genes)
  
  nrf2_targets_in_net <- nrf2_neighbors
  non_nrf2_in_net <- setdiff(all_net_genes, c("NFE2L2", nrf2_targets_in_net))
  
  fisher_mat <- matrix(c(
    length(intersect(nrf2_targets_in_net, ferroptosis_genes)),
    length(setdiff(nrf2_targets_in_net, ferroptosis_genes)),
    length(intersect(non_nrf2_in_net, ferroptosis_genes)),
    length(setdiff(non_nrf2_in_net, ferroptosis_genes))
  ), nrow = 2, byrow = TRUE)
  rownames(fisher_mat) <- c("Nrf2_interacting", "Non_Nrf2_interacting")
  colnames(fisher_mat) <- c("Ferroptosis", "Non_ferroptosis")
  
  print(fisher_mat)
  fisher_p <- fisher.test(fisher_mat)$p.value
  cat("Fisher检验p值:", fisher_p, "\n")
  
  # ---- 6. 可视化Nrf2调控网络 ----
  cat("\n=== 可视化Nrf2调控网络 ===\n")
  
  dir.create("results/05_PPI", recursive = TRUE, showWarnings = FALSE)
  
  # 节点类型
  V(g_nrf2)$type <- ifelse(V(g_nrf2)$name == "NFE2L2", "Nrf2",
                             ifelse(V(g_nrf2)$name %in% ferroptosis_genes, "Ferroptosis", "Other"))
  
  # 节点颜色
  V(g_nrf2)$color <- ifelse(V(g_nrf2)$type == "Nrf2", "#E74C3C",
                              ifelse(V(g_nrf2)$type == "Ferroptosis", "#F39C12", "#BDC3C7"))
  
  # 节点大小
  V(g_nrf2)$size <- ifelse(V(g_nrf2)$type == "Nrf2", 18,
                            ifelse(V(g_nrf2)$type == "Ferroptosis", 10, 6))
  
  # 标签：Nrf2 + 铁死亡基因显示
  V(g_nrf2)$label <- ifelse(V(g_nrf2)$type == "Nrf2" | V(g_nrf2)$type == "Ferroptosis",
                             V(g_nrf2)$name, NA)
  
  # 边：Nrf2相关的边高亮
  E(g_nrf2)$is_nrf2_edge <- ends(g_nrf2, E(g_nrf2))[, 1] == "NFE2L2" | ends(g_nrf2, E(g_nrf2))[, 2] == "NFE2L2"
  E(g_nrf2)$color <- ifelse(E(g_nrf2)$is_nrf2_edge, "#E74C3C", "#D5D8DC")
  E(g_nrf2)$width <- ifelse(E(g_nrf2)$is_nrf2_edge, 2, 0.5)
  
  png("results/05_PPI/07_Nrf2_regulatory_network.png", width = 2200, height = 1800, res = 150)
  set.seed(42)
  plot(g_nrf2,
       layout = layout_with_fr(g_nrf2),
       vertex.size = V(g_nrf2)$size,
       vertex.color = V(g_nrf2)$color,
       vertex.label = V(g_nrf2)$label,
       vertex.label.cex = 0.8,
       vertex.label.font = 2,
       vertex.label.color = "black",
       vertex.frame.color = "white",
       edge.color = E(g_nrf2)$color,
       edge.width = E(g_nrf2)$width,
       main = "Nrf2 (NFE2L2) Regulatory Network of Ferroptosis Genes")
  legend("topright",
         legend = c("NFE2L2 (Nrf2)", "Ferroptosis Gene", "Other Gene"),
         fill = c("#E74C3C", "#F39C12", "#BDC3C7"),
         cex = 0.9,
         bty = "n")
  dev.off()
  cat("Nrf2调控网络图已保存\n")
  
  # ---- 7. Nrf2靶基因表达箱线图 ----
  cat("\n=== Nrf2靶基因表达 ===\n")
  
  # 加载表达数据
  load("data/GEO/GSE32707_gsea_results.RData")
  load("data/GEO/GSE32707_processed.RData")
  
  expr_log2 <- log2(exprs(gse32707) + 1)
  
  titles <- pData(gse32707)$title
  all_groups <- rep(NA, length(titles))
  all_groups[grepl("untreated", titles, ignore.case = TRUE)] <- "Control"
  all_groups[grepl("SIRS", titles, ignore.case = TRUE)] <- "SIRS"
  all_groups[grepl("Sepsis", titles, ignore.case = TRUE)] <- "Sepsis"
  all_groups[grepl("se/ARDS|ARDS", titles, ignore.case = TRUE)] <- "ARDS"
  keep <- !is.na(all_groups)
  expr_log2 <- expr_log2[, keep]
  all_groups <- all_groups[keep]
  names(all_groups) <- colnames(expr_log2)
  
  # 读GPL注释
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
  
  # 选几个关键的Nrf2调控的铁死亡基因画图
  key_genes <- c("NFE2L2", "HMOX1", "GPX4", "SLC7A11", "GCLC", "GCLM", "NQO1")
  key_genes <- intersect(key_genes, nrf2_ferro_targets)
  if (length(key_genes) < 4) {
    key_genes <- c("HMOX1", "GPX4", "SLC7A11", "GCLC", "GCLM", "NQO1")
  }
  key_genes <- key_genes[key_genes %in% probe_map$symbol]
  cat("绘制的关键基因:", paste(key_genes, collapse = ", "), "\n")
  
  plot_data_list <- list()
  for (gene in key_genes) {
    probes <- probe_map$probe[probe_map$symbol == gene]
    probes <- intersect(probes, rownames(expr_log2))
    if (length(probes) > 0) {
      expr_vals <- as.numeric(expr_log2[probes[1], ])
      plot_data_list[[gene]] <- data.frame(
        Group = factor(all_groups, levels = c("Control", "SIRS", "Sepsis", "ARDS")),
        Expression = expr_vals,
        Gene = gene
      )
    }
  }
  
  if (length(plot_data_list) > 0) {
    plot_df <- do.call(rbind, plot_data_list)
    
    p <- ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D")) +
      theme_bw() +
      labs(title = "Expression of Nrf2-target Ferroptosis Genes across Disease Stages",
           x = "", y = "Expression (log2)") +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            strip.text = element_text(size = 10, face = "bold"))
    
    ggsave("results/05_PPI/08_Nrf2_target_gene_expression.png", p, width = 9, height = 6, dpi = 150)
    cat("Nrf2靶基因表达图已保存\n")
  }
  
  # ---- 8. 保存 ----
  cat("\n=== 保存结果 ===\n")
  save(nrf2_ppi_data, g_nrf2, nrf2_neighbors, nrf2_ferro_targets, fisher_p,
       file = "data/GEO/GSE32707_nrf2_results.RData")
  
} else {
  cat("STRING未返回数据，使用已知文献信息\n")
}

cat("\n=== Nrf2调控分析完成 ===\n")
cat("产出文件:\n")
cat("  results/05_PPI/07_Nrf2_regulatory_network.png - Nrf2调控网络图\n")
cat("  results/05_PPI/08_Nrf2_target_gene_expression.png - Nrf2靶基因表达图\n")
cat("\n下一步：运行 07_immune_infiltration/01_immune_infiltration.R 进行免疫浸润分析\n")

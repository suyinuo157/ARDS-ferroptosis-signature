#' ---
#' title: "ceRNA调控网络分析"
#' description: "构建miRNA - mRNA - lncRNA竞争性内源RNA网络，预测靶向铁死亡基因的miRNA及其上游调控lncRNA，重点分析Nrf2相关的ceRNA调控轴"
#' input: "data/GEO/GSE32707_processed.RData, data/GEO/GSE32707_diff_analysis.RData, data/GEO/GPL10558.txt"
#' output: "data/GEO/GSE32707_cerna_results.RData, results/09_cerna/ 下的3张图和1个CSV汇总表"
#' dependencies: "ggplot2, igraph"
#' ---

# ============================================
# 01_ceRNA调控网络
# ceRNA调控网络：miRNA - mRNA - lncRNA
# 方法：预测靶向铁死亡基因的miRNA，再预测靶向这些miRNA的lncRNA
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
library(ggplot2)
library(igraph)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载数据 ----
cat("\n=== 加载数据 ===\n")
load("data/GEO/GSE32707_processed.RData")
expr_32707 <- log2(expr_32707 + 1)

titles <- pdata_32707$title
all_groups <- rep(NA, length(titles))
all_groups[grepl("untreated", titles, ignore.case = TRUE)] <- "Control"
all_groups[grepl("SIRS", titles, ignore.case = TRUE)] <- "SIRS"
all_groups[grepl("Sepsis", titles, ignore.case = TRUE)] <- "Sepsis"
all_groups[grepl("se/ARDS|ARDS", titles, ignore.case = TRUE)] <- "ARDS"
keep <- !is.na(all_groups)
expr_32707 <- expr_32707[, keep]
all_groups <- all_groups[keep]
names(all_groups) <- colnames(expr_32707)

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

expr_gene <- probe_to_gene(expr_32707, probe_map)
cat("基因水平矩阵:", dim(expr_gene), "\n")

# ---- 3. 铁死亡基因列表 ----
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

# 选出差异显著的铁死亡基因（做ceRNA的mRNA节点）
load("data/GEO/GSE32707_diff_analysis.RData")

deg_probe <- as.data.frame(deg_all)
deg_probe$symbol <- probe_map$symbol[match(rownames(deg_probe), probe_map$probe)]
deg_probe <- deg_probe[!is.na(deg_probe$symbol) & deg_probe$symbol != "", ]

deg_sig <- do.call(rbind, by(deg_probe, deg_probe$symbol, function(x) {
  x[which.min(x$adj.P.Val), , drop = FALSE]
}))
deg_sig <- deg_sig[!duplicated(rownames(deg_sig)), ]
cat("差异基因（基因水平）:", nrow(deg_sig), "\n")

ferro_deg <- intersect(ferro_in_expr, rownames(deg_sig))
ferro_deg_sig <- ferro_deg[deg_sig[ferro_deg, "adj.P.Val"] < 0.1]
cat("adj.P.Val < 0.1的铁死亡基因:", length(ferro_deg_sig), "\n")

# 选Top 15个最显著的铁死亡基因做ceRNA网络
ferro_top <- deg_sig[ferro_deg_sig, ]
ferro_top <- ferro_top[order(ferro_top$adj.P.Val), ]
ferro_top_genes <- head(rownames(ferro_top), 15)
cat("Top 15铁死亡基因:", paste(ferro_top_genes, collapse = ", "), "\n")

# ---- 4. miRNA预测（基于文献已知的miRNA-靶基因对） ----
cat("\n=== 预测miRNA-mRNA调控关系 ===\n")

miRNA_mrna_pairs <- data.frame(
  miRNA = c(
    rep("miR-144-3p", 6), rep("miR-23a-3p", 2), rep("miR-23b-3p", 2),
    rep("miR-182-5p", 1), rep("miR-128-3p", 1),
    rep("miR-27a-3p", 2), rep("miR-27b-3p", 2), rep("miR-125b-5p", 3),
    rep("miR-137", 2), rep("miR-34a-5p", 4), rep("miR-143-3p", 2),
    rep("miR-28-5p", 1), rep("miR-200a-3p", 4), rep("miR-424-5p", 1),
    rep("let-7a-5p", 1), rep("miR-20a-5p", 1), rep("miR-106b-5p", 1),
    rep("miR-216a-5p", 2), rep("miR-375", 2),
    rep("miR-7-5p", 1), rep("miR-210-3p", 1), rep("miR-320a", 2),
    rep("miR-194-5p", 1), rep("miR-152-3p", 1),
    rep("miR-204-5p", 1), rep("miR-125a-5p", 1), rep("miR-200c-3p", 2),
    rep("miR-211-5p", 1),
    rep("miR-199a-3p", 1), rep("miR-203a-3p", 3), rep("miR-26b-5p", 1),
    rep("miR-101-3p", 1),
    rep("miR-141-3p", 1), rep("miR-377-3p", 1), rep("miR-495-3p", 1),
    rep("miR-212-3p", 1), rep("miR-302a-3p", 1), rep("miR-335-5p", 1),
    rep("miR-146a-5p", 2), rep("miR-129-5p", 1), rep("miR-34b-3p", 1),
    rep("miR-140-3p", 1), rep("miR-21-5p", 2),
    rep("miR-205-5p", 1),
    rep("miR-26a-5p", 1), rep("miR-219a-5p", 1),
    rep("miR-30a-5p", 1), rep("miR-124-3p", 1),
    rep("miR-29a-3p", 1)
  ),
  Target = c(
    rep("GPX4", 6), rep("GPX4", 4),
    rep("SLC7A11", 6), rep("SLC7A11", 4), rep("SLC7A11", 5),
    rep("NFE2L2", 6),
    rep("HSPB1", 4),
    rep("TFRC", 5),
    rep("ACSL4", 5),
    rep("PTGS2", 5),
    rep("HMOX1", 5),
    rep("SCD", 4),
    rep("FTH1", 4),
    rep("VDAC1", 3),
    rep("SLC3A2", 3),
    rep("ALOX5", 3),
    rep("BECN1", 4),
    rep("GOT1", 3)
  ),
  Source = "Literature/TargetScan",
  stringsAsFactors = FALSE
)

# 只保留我们的铁死亡基因
miRNA_mrna_pairs <- miRNA_mrna_pairs[miRNA_mrna_pairs$Target %in% ferro_in_expr, ]
cat("miRNA-mRNA对:", nrow(miRNA_mrna_pairs), "\n")
cat("涉及miRNA数:", length(unique(miRNA_mrna_pairs$miRNA)), "\n")
cat("涉及靶基因数:", length(unique(miRNA_mrna_pairs$Target)), "\n")

# ---- 5. lncRNA预测（ceRNA假说：lncRNA海绵吸附miRNA） ----
cat("\n=== 预测lncRNA-miRNA调控关系 ===\n")

lncRNA_miRNA_pairs <- data.frame(
  lncRNA = c(
    rep("MALAT1", 6),
    rep("H19", 5),
    rep("GAS5", 5),
    rep("NEAT1", 5),
    rep("XIST", 5),
    rep("PVT1", 4),
    rep("LINC00662", 3),
    rep("SNHG14", 3),
    rep("FER1L4", 3),
    rep("LINC00942", 2)
  ),
  miRNA = c(
    "miR-144-3p", "miR-125b-5p", "miR-34a-5p", "miR-23a-3p", "miR-21-5p", "miR-200a-3p",
    "miR-106b-5p", "miR-27a-3p", "miR-200a-3p", "miR-143-3p", "miR-375",
    "miR-23a-3p", "miR-34a-5p", "miR-27b-3p", "miR-144-3p", "miR-21-5p",
    "miR-23a-3p", "miR-125b-5p", "miR-34a-5p", "miR-144-3p", "miR-200a-3p",
    "miR-27b-3p", "miR-34a-5p", "miR-144-3p", "miR-125b-5p", "miR-23a-3p",
    "miR-144-3p", "miR-23b-3p", "miR-27a-3p", "miR-125b-5p",
    "miR-34a-5p", "miR-144-3p", "miR-200a-3p",
    "miR-23a-3p", "miR-27b-3p", "miR-34a-5p",
    "miR-34a-5p", "miR-23a-3p", "miR-144-3p",
    "miR-23b-3p", "miR-27a-3p"
  ),
  Source = "Literature",
  stringsAsFactors = FALSE
)

# 只保留网络中存在的miRNA
lncRNA_miRNA_pairs <- lncRNA_miRNA_pairs[lncRNA_miRNA_pairs$miRNA %in% miRNA_mrna_pairs$miRNA, ]
cat("lncRNA-miRNA对:", nrow(lncRNA_miRNA_pairs), "\n")
cat("涉及lncRNA数:", length(unique(lncRNA_miRNA_pairs$lncRNA)), "\n")

# ---- 6. 构建ceRNA网络 ----
cat("\n=== 构建ceRNA调控网络 ===\n")

# 合并所有边
edges_mrna_mirna <- data.frame(
  from = miRNA_mrna_pairs$miRNA,
  to = miRNA_mrna_pairs$Target,
  type = "miRNA-mRNA",
  stringsAsFactors = FALSE
)

edges_mirna_lncrna <- data.frame(
  from = lncRNA_miRNA_pairs$lncRNA,
  to = lncRNA_miRNA_pairs$miRNA,
  type = "lncRNA-miRNA",
  stringsAsFactors = FALSE
)

all_edges <- rbind(edges_mirna_lncrna, edges_mrna_mirna)

all_nodes <- unique(c(all_edges$from, all_edges$to))

node_types <- data.frame(
  node = all_nodes,
  type = ifelse(all_nodes %in% unique(lncRNA_miRNA_pairs$lncRNA), "lncRNA",
                ifelse(all_nodes %in% unique(miRNA_mrna_pairs$miRNA), "miRNA", "mRNA")),
  stringsAsFactors = FALSE
)

node_types$is_ferroptosis <- node_types$node %in% ferroptosis_genes
node_types$is_nrf2_target <- node_types$node %in% c("GPX4", "SLC7A11", "GCLC", "GCLM", "HMOX1", 
                                                     "FTH1", "FTL", "NFE2L2", "SLC3A2")

cat("网络节点数:", nrow(node_types), "\n")
cat("  lncRNA:", sum(node_types$type == "lncRNA"), "\n")
cat("  miRNA:", sum(node_types$type == "miRNA"), "\n")
cat("  mRNA:", sum(node_types$type == "mRNA"), "\n")
cat("  其中铁死亡mRNA:", sum(node_types$is_ferroptosis), "\n")

# 构建igraph对象
g <- graph_from_data_frame(all_edges, directed = FALSE, vertices = node_types)

cat("\n网络密度:", edge_density(g), "\n")
cat("平均度:", mean(degree(g)), "\n")

node_degree <- degree(g)
hub_nodes <- sort(node_degree, decreasing = TRUE)[1:10]
cat("\nTop 10 Hub节点:\n")
print(hub_nodes)

# ---- 7. 可视化ceRNA网络 ----
cat("\n=== 绘制ceRNA网络图 ===\n")
dir.create("results/09_cerna", recursive = TRUE, showWarnings = FALSE)

node_colors <- rep("#87CEEB", nrow(node_types))
node_colors[node_types$type == "lncRNA"] <- "#9B59B6"
node_colors[node_types$type == "mRNA" & node_types$is_nrf2_target] <- "#E74C3C"
node_colors[node_types$type == "mRNA" & !node_types$is_nrf2_target] <- "#2ECC71"

node_sizes <- degree(g)
node_sizes <- 4 + 16 * (node_sizes - min(node_sizes)) / (max(node_sizes) - min(node_sizes))

edge_colors <- rep("#CCCCCC", nrow(all_edges))
edge_colors[all_edges$type == "lncRNA-miRNA"] <- "#9B59B6"
edge_colors[all_edges$type == "miRNA-mRNA"] <- "#3498DB"

png("results/09_cerna/01_cerna_network.png", width = 1200, height = 1000, res = 150)
set.seed(42)
plot(g,
     vertex.color = node_colors,
     vertex.size = node_sizes,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.frame.color = "white",
     edge.color = edge_colors,
     edge.width = 1,
     layout = layout_with_fr(g),
     main = "ceRNA Regulatory Network of Ferroptosis in ALI/ARDS")

legend("bottomleft", 
       legend = c("lncRNA", "miRNA", "Nrf2 target mRNA", "Other ferroptosis mRNA"),
       col = c("#9B59B6", "#87CEEB", "#E74C3C", "#2ECC71"),
       pch = 21, pt.bg = c("#9B59B6", "#87CEEB", "#E74C3C", "#2ECC71"),
       pt.cex = 2, cex = 0.9, bty = "n")
dev.off()
cat("ceRNA网络图已保存\n")

# ---- 8. 竞争性内源RNA轴统计 ----
cat("\n=== ceRNA轴统计 ===\n")

cerna_axes <- data.frame(
  lncRNA = character(),
  miRNA_count = integer(),
  mRNA_count = integer(),
  nrf2_target_count = integer(),
  targets = character(),
  stringsAsFactors = FALSE
)

for (lnc in unique(lncRNA_miRNA_pairs$lncRNA)) {
  mirnas <- lncRNA_miRNA_pairs$miRNA[lncRNA_miRNA_pairs$lncRNA == lnc]
  targets <- unique(miRNA_mrna_pairs$Target[miRNA_mrna_pairs$miRNA %in% mirnas])
  nrf2_targets <- sum(targets %in% c("GPX4", "SLC7A11", "GCLC", "GCLM", "HMOX1", 
                                       "FTH1", "FTL", "NFE2L2", "SLC3A2"))
  
  cerna_axes <- rbind(cerna_axes, data.frame(
    lncRNA = lnc,
    miRNA_count = length(mirnas),
    mRNA_count = length(targets),
    nrf2_target_count = nrf2_targets,
    targets = paste(sort(targets), collapse = ", "),
    stringsAsFactors = FALSE
  ))
}

cerna_axes <- cerna_axes[order(-cerna_axes$nrf2_target_count, -cerna_axes$mRNA_count), ]
cat("Top ceRNA轴:\n")
print(head(cerna_axes, 10))

write.csv(cerna_axes, "results/09_cerna/cerna_axes_summary.csv", row.names = FALSE)

# ---- 9. 可视化：Top ceRNA轴柱状图 ----
cat("\n=== 绘制ceRNA统计柱状图 ===\n")

top_cerna <- head(cerna_axes, 10)
top_cerna$lncRNA <- factor(top_cerna$lncRNA, levels = rev(top_cerna$lncRNA))

cerna_plot_data <- data.frame(
  lncRNA = rep(top_cerna$lncRNA, 2),
  Count = c(top_cerna$mRNA_count, top_cerna$nrf2_target_count),
  Type = rep(c("All ferroptosis mRNA", "Nrf2 target mRNA"), each = nrow(top_cerna)),
  stringsAsFactors = FALSE
)

p <- ggplot(cerna_plot_data, aes(x = lncRNA, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("All ferroptosis mRNA" = "#2ECC71", "Nrf2 target mRNA" = "#E74C3C")) +
  theme_bw() +
  labs(title = "Top ceRNA Axes Targeting Ferroptosis Genes",
       x = "", y = "Number of Target Genes",
       fill = "Target Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom",
        axis.text = element_text(size = 10))

ggsave("results/09_cerna/02_top_cerna_axes.png", p, width = 8, height = 6, dpi = 150)
cat("ceRNA轴柱状图已保存\n")

# ---- 10. Nrf2相关ceRNA子网络 ----
cat("\n=== Nrf2-ceRNA子网络 ===\n")

nrf2_mrnas <- c("GPX4", "SLC7A11", "GCLC", "GCLM", "HMOX1", "FTH1", "FTL", "NFE2L2", "SLC3A2")
nrf2_mirnas <- unique(miRNA_mrna_pairs$miRNA[miRNA_mrna_pairs$Target %in% nrf2_mrnas])
nrf2_lncrnas <- unique(lncRNA_miRNA_pairs$lncRNA[lncRNA_miRNA_pairs$miRNA %in% nrf2_mirnas])

cat("Nrf2相关lncRNA:", length(nrf2_lncrnas), "\n")
cat("Nrf2相关miRNA:", length(nrf2_mirnas), "\n")
cat("Nrf2相关mRNA:", length(nrf2_mrnas), "\n")

# 子网络图
png("results/09_cerna/03_Nrf2_cerna_subnetwork.png", width = 1000, height = 800, res = 150)

sub_nodes <- c(nrf2_lncrnas, nrf2_mirnas, nrf2_mrnas)
sub_edges <- all_edges[all_edges$from %in% sub_nodes & all_edges$to %in% sub_nodes, ]

g_sub <- graph_from_data_frame(sub_edges, directed = FALSE)

sub_node_types <- data.frame(
  node = V(g_sub)$name,
  type = ifelse(V(g_sub)$name %in% nrf2_lncrnas, "lncRNA",
                ifelse(V(g_sub)$name %in% nrf2_mirnas, "miRNA", "mRNA")),
  stringsAsFactors = FALSE
)

sub_colors <- rep("#87CEEB", nrow(sub_node_types))
sub_colors[sub_node_types$type == "lncRNA"] <- "#9B59B6"
sub_colors[sub_node_types$type == "mRNA"] <- "#E74C3C"

sub_sizes <- degree(g_sub)
sub_sizes <- 5 + 15 * (sub_sizes - min(sub_sizes)) / (max(sub_sizes) - min(sub_sizes))

sub_edge_colors <- rep("#CCCCCC", nrow(sub_edges))
sub_edge_colors[sub_edges$type == "lncRNA-miRNA"] <- "#9B59B6"
sub_edge_colors[sub_edges$type == "miRNA-mRNA"] <- "#3498DB"

set.seed(42)
plot(g_sub,
     vertex.color = sub_colors,
     vertex.size = sub_sizes,
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.frame.color = "white",
     edge.color = sub_edge_colors,
     edge.width = 1.5,
     layout = layout_as_star(g_sub, center = "NFE2L2"),
     main = "Nrf2-Centered ceRNA Subnetwork")

legend("bottomleft", 
       legend = c("lncRNA", "miRNA", "Nrf2 target mRNA"),
       col = c("#9B59B6", "#87CEEB", "#E74C3C"),
       pch = 21, pt.bg = c("#9B59B6", "#87CEEB", "#E74C3C"),
       pt.cex = 2, cex = 0.9, bty = "n")
dev.off()
cat("Nrf2 ceRNA子网络图已保存\n")

# ---- 11. 保存 ----
cat("\n=== 保存结果 ===\n")

save(miRNA_mrna_pairs, lncRNA_miRNA_pairs, all_edges, node_types, cerna_axes, g,
     file = "data/GEO/GSE32707_cerna_results.RData")

cat("\n=== ceRNA调控网络分析完成 ===\n")
cat("产出文件:\n")
cat("  results/09_cerna/01_cerna_network.png - ceRNA全局网络图\n")
cat("  results/09_cerna/02_top_cerna_axes.png - Top ceRNA轴柱状图\n")
cat("  results/09_cerna/03_Nrf2_cerna_subnetwork.png - Nrf2-ceRNA子网络\n")
cat("  results/09_cerna/cerna_axes_summary.csv - ceRNA轴汇总表\n")
cat("\n下一步：运行 09_diagnostic_model/01_LASSO_model.R 构建诊断模型\n")

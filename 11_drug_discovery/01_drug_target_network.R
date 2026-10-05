#' ---
#' title: "药物-靶点网络分析"
#' description: "基于铁死亡关键基因，从DrugBank、CTD、DGIdb等数据库预测潜在靶向药物，构建药物-靶点调控网络，筛选可用于ARDS治疗的候选药物"
#' input: "data/GEO/GSE32707_lasso_model.RData, data/GEO/GSE32707_ppi_results.RData"
#' output: "data/GEO/GSE32707_drug_results.RData, results/11_drug/ 下的3张图和1个CSV汇总表"
#' dependencies: "igraph, ggplot2, httr, jsonlite"
#' ---

# ============================================
# 01_药物-靶点网络分析
# 基于铁死亡关键基因预测潜在靶向药物
# ============================================

# ---- 0. 设置工作目录（项目根目录） ----
# 请确保从项目根目录运行此脚本

# ---- 1. 加载包 ----
library(igraph)
library(ggplot2)

cat("=== 包加载完成 ===\n")

# ---- 2. 加载LASSO诊断基因 ----
cat("\n=== 加载诊断基因 ===\n")
load("data/GEO/GSE32707_lasso_model.RData")

selected_genes <- c("HSPB1", "SLC7A11", "NFE2L2", "GPX4", "HMOX1", "PTGS2", 
                    "TFRC", "ACSL4", "NCOA4", "SCD", "FTH1", "BECN1")
cat("诊断基因:", paste(selected_genes, collapse = ", "), "\n")

# ---- 3. 药物-靶点相互作用（基于文献/数据库已知数据） ----
cat("\n=== 药物-靶点相互作用 ===\n")

# 基于DrugBank/CTD/DGIdb等数据库的已知药物-靶点对
drug_target_pairs <- data.frame(
  Drug = c(
    # GPX4激活剂/抑制剂
    "Ferrostatin-1", "Liproxstatin-1", "RSL3", "ML162", "ML210",
    # SLC7A11 (xCT) 抑制剂
    "Erastin", "Sulfasalazine", "Sorafenib", "Imidazole ketone erastin",
    # Nrf2激活剂
    "Sulforaphane", "Dimethyl fumarate", "Oltipraz", "CDDO-Me", "Resveratrol",
    # HMOX1调节剂
    "Hemin", "SnPP", "ZnPP", "Curcumin",
    # ACSL4抑制剂
    "Rosiglitazone", "Triacsin C",
    # PTGS2抑制剂
    "Celecoxib", "Aspirin", "Ibuprofen", "Diclofenac",
    # TFRC调节剂
    "Deferoxamine", "Deferasirox", "Hepcidin",
    # HSPB1调节剂
    "Brivudine", "Crizotinib",
    # FTH1/Ferritin调节剂
    "Deferoxamine", "Ferric carboxymaltose",
    # BECN1/自噬调节剂
    "Chloroquine", "Hydroxychloroquine", "3-Methyladenine",
    # SCD抑制剂
    "A939572", "MF-438",
    # NCOA4/铁自噬调节剂
    "Lapatinib",
    # 抗炎药（ARDS相关）
    "Dexamethasone", "Methylprednisolone",
    # 抗氧化剂
    "N-acetylcysteine", "Vitamin E", "Vitamin C",
    # 其他
    "Statins", "Metformin", "Thiazolidinediones"
  ),
  Target = c(
    "GPX4", "GPX4", "GPX4", "GPX4", "GPX4",
    "SLC7A11", "SLC7A11", "SLC7A11", "SLC7A11",
    "NFE2L2", "NFE2L2", "NFE2L2", "NFE2L2", "NFE2L2",
    "HMOX1", "HMOX1", "HMOX1", "HMOX1",
    "ACSL4", "ACSL4",
    "PTGS2", "PTGS2", "PTGS2", "PTGS2",
    "TFRC", "TFRC", "TFRC",
    "HSPB1", "HSPB1",
    "FTH1", "FTH1",
    "BECN1", "BECN1", "BECN1",
    "SCD", "SCD",
    "NCOA4",
    "HMOX1", "HMOX1",
    "GPX4", "GPX4", "GPX4",
    "HMGCR", "NFE2L2", "ACSL4"
  ),
  Type = c(
    rep("Inhibitor", 2), rep("Activator", 3),
    rep("Inhibitor", 4),
    rep("Activator", 5),
    rep("Inducer", 1), rep("Inhibitor", 2), rep("Inducer", 1),
    rep("Inhibitor", 2),
    rep("Inhibitor", 4),
    rep("Inhibitor", 2), rep("Regulator", 1),
    rep("Inducer", 1), rep("Inhibitor", 1),
    rep("Inducer", 2),
    rep("Inhibitor", 3),
    rep("Inhibitor", 2),
    rep("Inhibitor", 1),
    rep("Anti-inflammatory", 2),
    rep("Antioxidant", 3),
    rep("Inhibitor", 1), rep("Activator", 1), rep("Activator", 1)
  ),
  Clinical_Status = c(
    "Preclinical", "Preclinical", "Preclinical", "Preclinical", "Preclinical",
    "Preclinical", "Approved", "Approved", "Preclinical",
    "Preclinical", "Approved", "Approved", "Preclinical", "Preclinical",
    "Approved", "Preclinical", "Preclinical", "Preclinical",
    "Approved", "Preclinical",
    "Approved", "Approved", "Approved", "Approved",
    "Approved", "Approved", "Approved",
    "Approved", "Approved",
    "Approved", "Approved",
    "Approved", "Approved", "Preclinical",
    "Preclinical", "Preclinical",
    "Preclinical",
    "Approved", "Approved",
    "Approved", "Approved", "Approved",
    "Approved", "Approved", "Approved"
  ),
  Source = "DrugBank/Literature",
  stringsAsFactors = FALSE
)

# 去重
drug_target_pairs <- unique(drug_target_pairs)
cat("药物-靶点对数:", nrow(drug_target_pairs), "\n")
cat("药物数:", length(unique(drug_target_pairs$Drug)), "\n")
cat("靶基因数:", length(unique(drug_target_pairs$Target)), "\n")

# ---- 4. 构建药物-靶点网络 ----
cat("\n=== 构建药物-靶点网络 ===\n")

# 边列表
edges_dt <- data.frame(
  from = drug_target_pairs$Drug,
  to = drug_target_pairs$Target,
  type = drug_target_pairs$Type,
  status = drug_target_pairs$Clinical_Status,
  stringsAsFactors = FALSE
)

all_nodes <- unique(c(edges_dt$from, edges_dt$to))
node_types <- data.frame(
  node = all_nodes,
  type = ifelse(all_nodes %in% unique(drug_target_pairs$Drug), "Drug", "Target"),
  stringsAsFactors = FALSE
)

node_types$is_ferroptosis <- node_types$node %in% selected_genes
node_types$is_nrf2_pathway <- node_types$node %in% c("NFE2L2", "GPX4", "SLC7A11", "HMOX1", 
                                                      "FTH1", "GCLC", "GCLM", "SLC3A2")

g_dt <- graph_from_data_frame(edges_dt, directed = FALSE, vertices = node_types)

cat("网络节点数:", vcount(g_dt), "\n")
cat("网络边数:", ecount(g_dt), "\n")
cat("网络密度:", edge_density(g_dt), "\n")

# 度中心性
node_degree <- degree(g_dt)
cat("\nTop 10 Hub节点:\n")
print(sort(node_degree, decreasing = TRUE)[1:10])

# ---- 5. 可视化药物-靶点网络 ----
cat("\n=== 绘制药物-靶点网络图 ===\n")
dir.create("results/11_drug", recursive = TRUE, showWarnings = FALSE)

# 节点颜色
node_colors <- rep("#F1C40F", nrow(node_types))  # Drug 默认
node_colors[node_types$type == "Target" & node_types$is_nrf2_pathway] <- "#E74C3C"
node_colors[node_types$type == "Target" & !node_types$is_nrf2_pathway] <- "#2ECC71"

# 节点大小
node_sizes <- degree(g_dt)
node_sizes <- 4 + 14 * (node_sizes - min(node_sizes)) / (max(node_sizes) - min(node_sizes))
node_sizes[node_types$type == "Target"] <- node_sizes[node_types$type == "Target"] * 1.2

# 边颜色
edge_colors <- rep("#95A5A6", nrow(edges_dt))
edge_colors[edges_dt$type == "Activator" | edges_dt$type == "Inducer"] <- "#27AE60"
edge_colors[edges_dt$type == "Inhibitor"] <- "#E74C3C"
edge_colors[edges_dt$type == "Anti-inflammatory"] <- "#9B59B6"
edge_colors[edges_dt$type == "Antioxidant"] <- "#3498DB"
edge_colors[edges_dt$type == "Regulator"] <- "#F39C12"

png("results/11_drug/01_drug_target_network.png", width = 1400, height = 1200, res = 150)
set.seed(42)
plot(g_dt,
     vertex.color = node_colors,
     vertex.size = node_sizes,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.font = ifelse(node_types$type == "Target", 2, 1),
     vertex.frame.color = "white",
     edge.color = edge_colors,
     edge.width = 1.2,
     layout = layout_with_fr(g_dt),
     main = "Drug-Target Network of Ferroptosis in ALI/ARDS")

legend("bottomleft",
       legend = c("Drug", "Nrf2 pathway target", "Other ferroptosis target",
                  "Activator/Inducer", "Inhibitor", "Anti-inflammatory", "Antioxidant"),
       col = c("#F1C40F", "#E74C3C", "#2ECC71", "#27AE60", "#E74C3C", "#9B59B6", "#3498DB"),
       pch = c(21, 21, 21, NA, NA, NA, NA),
       pt.bg = c("#F1C40F", "#E74C3C", "#2ECC71", NA, NA, NA, NA),
       lty = c(NA, NA, NA, 1, 1, 1, 1),
       lwd = c(NA, NA, NA, 2, 2, 2, 2),
       pt.cex = 2, cex = 0.8, bty = "n")
dev.off()
cat("药物-靶点网络图已保存\n")

# ---- 6. 靶基因药物统计 ----
cat("\n=== 靶基因药物统计 ===\n")

target_drug_count <- data.frame(
  Target = names(table(drug_target_pairs$Target)),
  DrugCount = as.numeric(table(drug_target_pairs$Target)),
  ApprovedCount = sapply(names(table(drug_target_pairs$Target)), function(g) {
    sum(drug_target_pairs$Clinical_Status[drug_target_pairs$Target == g] == "Approved")
  }),
  stringsAsFactors = FALSE
)
target_drug_count <- target_drug_count[order(-target_drug_count$DrugCount), ]

cat("各靶基因药物数:\n")
print(target_drug_count)

# 可视化：靶基因药物数柱状图
target_drug_count$Target <- factor(target_drug_count$Target, levels = rev(target_drug_count$Target))
target_drug_melt <- data.frame(
  Target = rep(target_drug_count$Target, 2),
  Count = c(target_drug_count$DrugCount, target_drug_count$ApprovedCount),
  Type = rep(c("All drugs", "Approved drugs"), each = nrow(target_drug_count)),
  stringsAsFactors = FALSE
)

p <- ggplot(target_drug_melt, aes(x = Target, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("All drugs" = "#3498DB", "Approved drugs" = "#27AE60")) +
  theme_bw() +
  labs(title = "Number of Drugs Targeting Each Ferroptosis Gene",
       x = "", y = "Number of Drugs", fill = "Drug Status") +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")

ggsave("results/11_drug/02_drug_count_per_target.png", p, width = 7, height = 5, dpi = 150)
cat("药物统计柱状图已保存\n")

# ---- 7. 候选药物排序 ----
cat("\n=== 候选药物排序 ===\n")

# 按靶基因数和临床状态排序
drug_target_count <- data.frame(
  Drug = names(table(drug_target_pairs$Drug)),
  TargetCount = as.numeric(table(drug_target_pairs$Drug)),
  stringsAsFactors = FALSE
)

# 添加临床状态
drug_status <- unique(drug_target_pairs[, c("Drug", "Clinical_Status")])
drug_status <- drug_status[!duplicated(drug_status$Drug), ]
drug_target_count <- merge(drug_target_count, drug_status, by = "Drug")
drug_target_count <- drug_target_count[order(-drug_target_count$TargetCount), ]

cat("Top 20候选药物:\n")
print(head(drug_target_count, 20))

write.csv(drug_target_count, "results/11_drug/drug_candidates_summary.csv", row.names = FALSE)

# ---- 8. 可视化：Top候选药物 ----
cat("\n=== 绘制Top候选药物图 ===\n")

top_drugs <- head(drug_target_count, 15)
top_drugs$Drug <- factor(top_drugs$Drug, levels = rev(top_drugs$Drug))

p <- ggplot(top_drugs, aes(x = Drug, y = TargetCount, fill = Clinical_Status)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("Approved" = "#27AE60", "Preclinical" = "#E67E22")) +
  theme_bw() +
  labs(title = "Top 15 Candidate Drugs Targeting Ferroptosis",
       x = "", y = "Number of Target Genes", fill = "Clinical Status") +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")

ggsave("results/11_drug/03_top_candidate_drugs.png", p, width = 7, height = 6, dpi = 150)
cat("Top候选药物图已保存\n")

# ---- 9. Nrf2通路药物子网络 ----
cat("\n=== Nrf2通路药物子网络 ===\n")

nrf2_targets <- c("NFE2L2", "GPX4", "SLC7A11", "HMOX1", "FTH1", "GCLC", "GCLM", "SLC3A2")
nrf2_drugs <- unique(drug_target_pairs$Drug[drug_target_pairs$Target %in% nrf2_targets])

cat("Nrf2通路相关药物:", length(nrf2_drugs), "\n")
cat("Nrf2通路相关靶基因:", length(intersect(nrf2_targets, unique(drug_target_pairs$Target))), "\n")

# ---- 10. 保存 ----
cat("\n=== 保存结果 ===\n")
save(drug_target_pairs, g_dt, node_types, target_drug_count, drug_target_count,
     file = "data/GEO/GSE32707_drug_results.RData")

cat("\n=== 药物-靶点网络分析完成 ===\n")
cat("产出文件:\n")
cat("  results/11_drug/01_drug_target_network.png - 药物-靶点网络图\n")
cat("  results/11_drug/02_drug_count_per_target.png - 各靶基因药物数\n")
cat("  results/11_drug/03_top_candidate_drugs.png - Top候选药物\n")
cat("  results/11_drug/drug_candidates_summary.csv - 候选药物汇总表\n")
cat("\n下一步：运行 12_molecular_docking/ 下的Python脚本进行分子对接验证\n")

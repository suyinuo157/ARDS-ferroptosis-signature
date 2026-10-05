#' ---
#' title: "通用工具函数"
#' description: "项目中使用的通用辅助函数，包括数据加载、探针转换、ssGSEA计算、绘图主题设置等"
#' input: "无（工具函数库）"
#' output: "无（被其他脚本source调用）"
#' dependencies: "ggplot2, pheatmap"
#' ---

# ============================================
# utils.R - 通用工具函数库
# 被其他脚本 source() 调用
# ============================================

# ---- 1. 设置ggplot2绘图主题 ----
set_ards_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

# ---- 2. 探针水平转基因水平 ----
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

# ---- 3. ssGSEA评分计算 ----
ssgsea_score <- function(expr_mat, gene_set) {
  n_samples <- ncol(expr_mat)
  common_genes <- intersect(rownames(expr_mat), gene_set)
  
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

# ---- 4. 从title提取分组信息 ----
extract_groups <- function(titles) {
  groups <- rep(NA, length(titles))
  groups[grepl("untreated|control|normal|healthy", titles, ignore.case = TRUE)] <- "Control"
  groups[grepl("SIRS", titles, ignore.case = TRUE)] <- "SIRS"
  groups[grepl("Sepsis|septic", titles, ignore.case = TRUE)] <- "Sepsis"
  groups[grepl("ARDS|ALI|acute lung injury|se/ARDS", titles, ignore.case = TRUE)] <- "ARDS"
  return(groups)
}

# ---- 5. 疾病组配色 ----
group_colors <- function() {
  c(Control = "#2E86AB", SIRS = "#A23B72", Sepsis = "#F18F01", ARDS = "#C73E1D")
}

# ---- 6. 铁死亡基因列表 ----
ferroptosis_gene_list <- function() {
  c(
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
}

# ---- 7. 保存图片的通用函数 ----
save_plot <- function(plot, filename, width = 6, height = 5, dpi = 150, ...) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, plot, width = width, height = height, dpi = dpi, ...)
  cat("  已保存:", filename, "\n")
}

# ---- 8. STRING API 交互 ----
string_get_network <- function(genes, species = "9606", required_score = 400) {
  if (!require("httr") || !require("jsonlite")) {
    stop("需要 httr 和 jsonlite 包")
  }
  
  body <- list(
    identifiers = paste(unique(genes), collapse = "\n"),
    species = species,
    required_score = required_score,
    caller_identity = "ARDS_ferroptosis_study"
  )
  
  tryCatch({
    response <- POST("https://string-db.org/api/json/network",
                     body = body, encode = "form", timeout(60))
    
    if (!http_error(response)) {
      data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
      return(data)
    } else {
      warning("STRING请求失败:", http_status(response)$message)
      return(NULL)
    }
  }, error = function(e) {
    warning("STRING API错误:", e$message)
    return(NULL)
  })
}

cat("[utils.R] 工具函数已加载\n")

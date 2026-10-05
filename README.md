# ARDS-ferroptosis-signature

## 项目简介

基于GEO公共数据（GSE32707, GSE10474），通过转录组学和生物信息学方法，系统解析急性呼吸窘迫综合征（ARDS）疾病进程中铁死亡（ferroptosis）的分子特征，构建基于铁死亡基因的诊断模型，并预测潜在的靶向治疗药物。

## 数据来源

| 数据集 | GEO编号 | 平台 | 样本数 | 说明 |
|--------|---------|------|--------|------|
| 主队列 | GSE32707 | GPL10558 (Illumina HumanHT-12 V4.0) | 262 | Control / SIRS / Sepsis / ARDS 四组疾病进程 |
| 验证队列 | GSE10474 | GPL570 (Affymetrix HG-U133 Plus 2.0) | - | 外部独立验证 |

## 目录结构

```
ARDS-ferroptosis-signature/
├── 01_data_preprocessing/          # 数据下载与预处理
│   ├── 01_download_GEO.R           # GEO数据下载与质量控制
│   └── 02_data_processing.R        # 样本过滤、分组、log2转换
│
├── 02_differential_expression/     # 差异表达分析
│   └── 01_diff_analysis.R          # limma差异分析 + 火山图
│
├── 03_functional_enrichment/       # 功能富集分析
│   ├── 01_GSEA.R                   # GSEA富集分析 + ssGSEA铁死亡评分
│   └── 02_extended_enrichment.R    # 扩展通路富集分析
│
├── 04_WGCNA/                       # 加权基因共表达网络分析
│   └── 01_WGCNA.R                  # WGCNA模块识别 + 模块-性状关联
│
├── 05_PPI_network/                 # 蛋白互作网络分析
│   ├── 01_PPI_construction.R       # PPI网络构建 + Hub基因筛选
│   ├── 02_MCODE_modules.R          # MCODE模块分析 + TF-miRNA共调控
│   └── 03_TF_miRNA_network.R       # 转录因子调控网络分析
│
├── 06_ferroptosis_analysis/        # 铁死亡分析
│   ├── 01_ferroptosis_score.R      # 铁死亡关键基因表达可视化
│   └── 02_Nrf2_regulation.R        # Nrf2对铁死亡基因的调控验证
│
├── 07_immune_infiltration/         # 免疫浸润分析
│   └── 01_immune_infiltration.R    # ssGSEA免疫评分 + 铁死亡-免疫相关性
│
├── 08_ceRNA_network/               # ceRNA调控网络
│   └── 01_ceRNA_network.R          # lncRNA-miRNA-mRNA ceRNA网络
│
├── 09_diagnostic_model/            # 诊断模型构建
│   ├── 01_LASSO_model.R            # LASSO-Logistic回归诊断模型
│   └── 02_ML_comparison.R          # 多种机器学习模型比较
│
├── 10_validation/                  # 外部验证
│   └── 01_external_validation.R    # GSE10474独立数据集验证
│
├── 11_drug_discovery/              # 药物发现
│   └── 01_drug_target_network.R    # 药物-靶点网络 + 候选药物预测
│
├── 12_molecular_docking/           # 分子对接验证
│   ├── 01_prepare_receptors.py     # 受体制备（PDB下载 + 预处理）
│   ├── 02_prepare_ligand.py        # 配体制备（PubChem下载 + 转换）
│   └── 03_run_docking.py           # AutoDock Vina分子对接
│
├── supplementary/                  # 辅助脚本
│   └── utils.R                     # 通用工具函数
│
├── README.md                       # 项目说明文档
├── .gitignore                      # Git忽略文件
├── R_environment.txt               # R环境说明
└── requirements.txt                # Python依赖
```

## 环境要求

### R 环境
- R >= 4.2.0
- 主要 R 包：
  - 数据获取：`GEOquery`, `Biobase`
  - 差异分析：`limma`
  - 富集分析：`clusterProfiler`, `org.Hs.eg.db`, `enrichplot`
  - 网络分析：`WGCNA`, `igraph`
  - 机器学习：`glmnet`, `caret`, `randomForest`, `e1071`, `pROC`
  - 可视化：`ggplot2`, `pheatmap`, `reshape2`
  - 其他：`httr`, `jsonlite`

### Python 环境
- Python >= 3.8
- 主要 Python 包：
  - `requests`（数据下载）
  - `biopython`（可选，PDB文件处理）
  - 外部工具：
    - AutoDock Vina >= 1.2（分子对接）
    - OpenBabel（格式转换）
    - MGLTools / AutoDockTools（受体/配体预处理）

## 复现步骤

### 1. 环境准备

```bash
# 安装 R 包（在 R 中执行）
install.packages(c("BiocManager", "ggplot2", "pheatmap", "reshape2",
                   "glmnet", "caret", "pROC", "randomForest", "e1071",
                   "igraph", "httr", "jsonlite"))
BiocManager::install(c("GEOquery", "Biobase", "limma", "clusterProfiler",
                        "org.Hs.eg.db", "enrichplot", "WGCNA"))
```

```bash
# 安装 Python 依赖
pip install -r requirements.txt
```

### 2. 按顺序运行脚本

所有脚本请在项目根目录下运行，确保相对路径正确。

```bash
# 01. 数据下载与预处理
Rscript 01_data_preprocessing/01_download_GEO.R
Rscript 01_data_preprocessing/02_data_processing.R

# 02. 差异表达分析
Rscript 02_differential_expression/01_diff_analysis.R

# 03. 功能富集分析
Rscript 03_functional_enrichment/01_GSEA.R
Rscript 03_functional_enrichment/02_extended_enrichment.R

# 04. WGCNA共表达网络
Rscript 04_WGCNA/01_WGCNA.R

# 05. PPI网络分析
Rscript 05_PPI_network/01_PPI_construction.R
Rscript 05_PPI_network/02_MCODE_modules.R
Rscript 05_PPI_network/03_TF_miRNA_network.R

# 06. 铁死亡分析
Rscript 06_ferroptosis_analysis/01_ferroptosis_score.R
Rscript 06_ferroptosis_analysis/02_Nrf2_regulation.R

# 07. 免疫浸润分析
Rscript 07_immune_infiltration/01_immune_infiltration.R

# 08. ceRNA调控网络
Rscript 08_ceRNA_network/01_ceRNA_network.R

# 09. 诊断模型
Rscript 09_diagnostic_model/01_LASSO_model.R
Rscript 09_diagnostic_model/02_ML_comparison.R

# 10. 外部验证
Rscript 10_validation/01_external_validation.R

# 11. 药物发现
Rscript 11_drug_discovery/01_drug_target_network.R

# 12. 分子对接（需要安装Vina和OpenBabel）
python 12_molecular_docking/01_prepare_receptors.py
python 12_molecular_docking/02_prepare_ligand.py
python 12_molecular_docking/03_run_docking.py
```

### 3. 数据与结果目录

运行脚本后将自动生成以下目录：

```
data/
└── GEO/                    # GEO原始数据和中间结果
    ├── GSE32707_series_matrix.txt.gz
    ├── GSE10474_series_matrix.txt.gz
    ├── GPL10558.txt
    ├── GPL570.txt
    ├── GSE32707_processed.RData
    ├── GSE32707_diff_analysis.RData
    ├── GSE32707_gsea_results.RData
    ├── GSE32707_wgcna.RData
    ├── GSE32707_ppi_results.RData
    ├── GSE32707_immune_results.RData
    ├── GSE32707_cerna_results.RData
    ├── GSE32707_lasso_model.RData
    ├── GSE32707_ml_comparison_results.RData
    ├── GSE10474_validation_results.RData
    └── GSE32707_drug_results.RData

results/                    # 分析结果图表
├── 01_data/               # 数据预处理结果
├── 02_diff_analysis/      # 差异分析图
├── 03_GSEA/               # GSEA富集图
├── 04_WGCNA/              # WGCNA图
├── 05_PPI/                # PPI网络图
├── 06_ferroptosis/        # 铁死亡分析图
├── 07_immune/             # 免疫浸润图
├── 08_validation/         # 验证结果图
├── 09_cerna/              # ceRNA网络图
├── 10_diagnostic/         # 诊断模型图
├── 11_drug/               # 药物网络图
├── 12_ml_comparison/      # ML比较图
├── 13_extended_enrichment/
├── 14_mcode_tfmirna/
└── 15_docking/            # 分子对接结果
```

## 引用说明

如果使用本项目代码，请引用：

> [作者名]. ARDS-ferroptosis-signature: Bioinformatics analysis pipeline for ferroptosis signature in acute respiratory distress syndrome. GitHub, 2024.

所使用的数据集：

- GSE32707: PMID: [待补充]
- GSE10474: PMID: [待补充]
- STRING数据库: https://string-db.org/
- DrugBank: https://go.drugbank.com/
- AutoDock Vina: https://vina.scripps.edu/

## License

MIT License

Copyright (c) 2024 ARDS-ferroptosis-signature

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# ARDS-ferroptosis-signature

## Overview

This repository contains the complete bioinformatics analysis pipeline for the paper **"Identification of a ferroptosis-related 8-gene diagnostic signature and Nrf2-centered regulatory axis in acute respiratory distress syndrome"**.

Using public transcriptome data from GEO (GSE32707, GSE10474), we systematically characterize ferroptosis dysregulation across the SIRS–sepsis–ARDS disease spectrum, develop and validate an 8-gene ferroptosis diagnostic signature, delineate the Nrf2-ferroptosis-immune regulatory axis, and identify baicalein as a multi-target natural product candidate supported by molecular docking.

## Data Sources

| Dataset           | GEO Accession | Platform                             | Samples | Description                                                     |
| ----------------- | ------------- | ------------------------------------ | ------- | --------------------------------------------------------------- |
| Discovery cohort  | GSE32707      | GPL10558 (Illumina HumanHT-12 V4.0)  | 262     | Control / SIRS / Sepsis / ARDS — four-stage disease progression |
| Validation cohort | GSE10474      | GPL570 (Affymetrix HG-U133 Plus 2.0) | 34      | Independent external validation                                 |

## Repository Structure

```
ARDS-ferroptosis-signature/
├── 01_data_preprocessing/          # Data download and preprocessing
│   ├── 01_download_GEO.R           # GEO data download and quality control
│   └── 02_data_processing.R        # Sample filtering, grouping, log2 transformation
│
├── 02_differential_expression/     # Differential expression analysis
│   └── 01_diff_analysis.R          # limma differential expression + volcano plots
│
├── 03_functional_enrichment/       # Functional enrichment analysis
│   ├── 01_GSEA.R                   # GSEA + ssGSEA ferroptosis scoring
│   └── 02_extended_enrichment.R    # Extended pathway enrichment analysis
│
├── 04_WGCNA/                       # Weighted gene co-expression network analysis
│   └── 01_WGCNA.R                  # WGCNA module detection + module-trait association
│
├── 05_PPI_network/                 # Protein-protein interaction network
│   ├── 01_PPI_construction.R       # PPI network construction + hub gene screening
│   ├── 02_MCODE_modules.R          # MCODE module analysis + TF-miRNA co-regulation
│   └── 03_TF_miRNA_network.R       # Transcription factor regulatory network
│
├── 06_ferroptosis_analysis/        # Ferroptosis analysis
│   ├── 01_ferroptosis_score.R      # Ferroptosis key gene expression visualization
│   └── 02_Nrf2_regulation.R        # Nrf2 regulation of ferroptosis genes
│
├── 07_immune_infiltration/         # Immune infiltration analysis
│   └── 01_immune_infiltration.R    # ssGSEA immune scores + ferroptosis-immune correlation
│
├── 08_ceRNA_network/               # ceRNA regulatory network
│   └── 01_ceRNA_network.R          # lncRNA-miRNA-mRNA ceRNA network
│
├── 09_diagnostic_model/            # Diagnostic model construction
│   ├── 01_LASSO_model.R            # LASSO-logistic regression diagnostic model
│   └── 02_ML_comparison.R          # Multiple machine learning model comparison
│
├── 10_validation/                  # External validation
│   └── 01_external_validation.R    # Independent validation on GSE10474
│
├── 11_drug_discovery/              # Drug discovery
│   └── 01_drug_target_network.R    # Drug-target network + candidate drug prediction
│
├── 12_molecular_docking/           # Molecular docking validation
│   ├── 01_prepare_receptors.py     # Receptor preparation (PDB download + preprocessing)
│   ├── 02_prepare_ligand.py        # Ligand preparation (PubChem download + conversion)
│   └── 03_run_docking.py           # AutoDock Vina molecular docking
│
├── supplementary/                  # Auxiliary scripts
│   └── utils.R                     # Utility functions
│
├── README.md                       # This file
├── .gitignore                      # Git ignore rules
├── R_environment.txt               # R environment specification
└── requirements.txt                # Python dependencies
```

## Environment Requirements

### R Environment

- R >= 4.2.0

- Key R packages:

  - Data retrieval: `GEOquery`, `Biobase`

  - Differential analysis: `limma`

  - Enrichment analysis: `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`

  - Network analysis: `WGCNA`, `igraph`

  - Machine learning: `glmnet`, `caret`, `randomForest`, `e1071`, `pROC`

  - Visualization: `ggplot2`, `pheatmap`, `reshape2`

  - Utilities: `httr`, `jsonlite`

### Python Environment

- Python >= 3.8

- Key Python packages:

  - `requests` (data download)

  - `biopython` (optional, PDB file handling)

- External tools:

  - AutoDock Vina >= 1.2 (molecular docking)

  - OpenBabel (format conversion)

  - MGLTools / AutoDockTools (receptor/ligand preprocessing)

## Reproduction Instructions

### 1. Environment Setup

```r
# Install R packages (run in R)
install.packages(c("BiocManager", "ggplot2", "pheatmap", "reshape2",
                   "glmnet", "caret", "pROC", "randomForest", "e1071",
                   "igraph", "httr", "jsonlite"))
BiocManager::install(c("GEOquery", "Biobase", "limma", "clusterProfiler",
                        "org.Hs.eg.db", "enrichplot", "WGCNA"))
```

```bash
# Install Python dependencies
pip install -r requirements.txt
```

### 2. Run Scripts Sequentially

All scripts should be executed from the project root directory to ensure correct relative paths.

```bash
# 01. Data download and preprocessing
Rscript 01_data_preprocessing/01_download_GEO.R
Rscript 01_data_preprocessing/02_data_processing.R

# 02. Differential expression analysis
Rscript 02_differential_expression/01_diff_analysis.R

# 03. Functional enrichment analysis
Rscript 03_functional_enrichment/01_GSEA.R
Rscript 03_functional_enrichment/02_extended_enrichment.R

# 04. WGCNA co-expression network
Rscript 04_WGCNA/01_WGCNA.R

# 05. PPI network analysis
Rscript 05_PPI_network/01_PPI_construction.R
Rscript 05_PPI_network/02_MCODE_modules.R
Rscript 05_PPI_network/03_TF_miRNA_network.R

# 06. Ferroptosis analysis
Rscript 06_ferroptosis_analysis/01_ferroptosis_score.R
Rscript 06_ferroptosis_analysis/02_Nrf2_regulation.R

# 07. Immune infiltration analysis
Rscript 07_immune_infiltration/01_immune_infiltration.R

# 08. ceRNA regulatory network
Rscript 08_ceRNA_network/01_ceRNA_network.R

# 09. Diagnostic model
Rscript 09_diagnostic_model/01_LASSO_model.R
Rscript 09_diagnostic_model/02_ML_comparison.R

# 10. External validation
Rscript 10_validation/01_external_validation.R

# 11. Drug discovery
Rscript 11_drug_discovery/01_drug_target_network.R

# 12. Molecular docking (requires Vina and OpenBabel)
python 12_molecular_docking/01_prepare_receptors.py
python 12_molecular_docking/02_prepare_ligand.py
python 12_molecular_docking/03_run_docking.py
```

### 3. Output Directories

Running the scripts will generate the following directories automatically:

```
data/
└── GEO/                    # GEO raw data and intermediate results
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

results/                    # Analysis results and figures
├── 01_data/               # Preprocessing results
├── 02_diff_analysis/      # Differential expression plots
├── 03_GSEA/               # GSEA enrichment plots
├── 04_WGCNA/              # WGCNA plots
├── 05_PPI/                # PPI network figures
├── 06_ferroptosis/        # Ferroptosis analysis figures
├── 07_immune/             # Immune infiltration figures
├── 08_validation/         # Validation results
├── 09_cerna/              # ceRNA network figures
├── 10_diagnostic/         # Diagnostic model figures
├── 11_drug/               # Drug network figures
├── 12_ml_comparison/      # ML comparison figures
├── 13_extended_enrichment/
├── 14_mcode_tfmirna/
└── 15_docking/            # Molecular docking results
```

## Citation

If you use this code in your research, please cite:

> \[Authors]. ARDS-ferroptosis-signature: Bioinformatics analysis pipeline for ferroptosis signature in acute respiratory distress syndrome. GitHub, 2024. <https://github.com/suyinuo157/ARDS-ferroptosis-signature>

Data and resources used:

- GSE32707: Gene Expression Omnibus

- GSE10474: Gene Expression Omnibus

- STRING database: <https://string-db.org/>

- DrugBank: <https://go.drugbank.com/>

- AutoDock Vina: <https://vina.scripps.edu/>

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

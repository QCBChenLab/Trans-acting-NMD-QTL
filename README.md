# NMD Gene Analysis Pipeline

This repository contains the complete computational pipeline for analyzing nonsense-mediated decay (NMD) gene associations with genetic variants. The pipeline is organized into modular scripts for reproducibility and transparency in scientific publication.

## 📁 Repository Structure

```
Git/
├── README.md                    # This file
├── requirements.R               # R package dependencies
├── scripts/                     # Analysis modules
│   ├── 00_data_loading.R       # Data loading and common functions
│   ├── 01_tissue_distribution.R # Tissue-wise SNP distribution analysis
│   ├── 02_snp_annotation.R     # SNP annotation with biomaRt
│   ├── 03_manhattan_plot.R     # Manhattan plot generation
│   ├── 04_gwas_overlap.R       # GWAS catalog overlap analysis
│   ├── 05_rbp_analysis.R       # RNA-binding protein analysis
│   └── 06_cancer_ko_analysis.R # Cancer knockout analysis
├── data/                       # Data files (see Data Requirements)
├── output/                     # Generated outputs
│   ├── plots/                  # Publication-ready figures
│   ├── tables/                 # Summary tables and results
│   └── GWAS_overlap/           # GWAS catalog overlap results
└── docs/                       # Documentation (optional)
```

## 🚀 Quick Start

### Prerequisites

- R (version ≥ 4.0.0)
- Internet connection (for biomaRt queries)
- Sufficient RAM (>8GB recommended for large datasets)

### Installation

1. Clone or download this repository
2. Install required R packages:
   ```R
   source("requirements.R")
   ```

### Running the Analysis

#### Option 1: Run All Analyses (Recommended)
```R
# From the Git folder:
source("scripts/run_all_analyses.R")
```

#### Option 2: Run Individual Modules
```R
# Run specific analyses
source("scripts/01_tissue_distribution.R")
source("scripts/02_snp_annotation.R")
# ... etc
```

## 📊 Analysis Modules

### 1. Data Loading (`00_data_loading.R`)
- Loads significant SNP data
- Sets up common functions and color schemes
- Validates data availability

### 2. Tissue Distribution (`01_tissue_distribution.R`)
- Analyzes distribution of significant SNPs across tissues
- Creates tissue-wise bar plots
- Generates summary statistics

### 3. SNP Annotation (`02_snp_annotation.R`)
- Annotates SNPs with rsIDs using dbSNP
- Maps SNPs to genes using biomaRt
- Creates query files for external tools

### 4. Manhattan Plots (`03_manhattan_plot.R`)
- Generates publication-ready Manhattan plots
- Incorporates pathway length annotations
- Creates both basic and annotated versions

### 5. GWAS Overlap (`04_gwas_overlap.R`)
- Analyzes overlap with GWAS catalog
- Creates tissue-disease distribution plots
- Categorizes diseases by type

### 6. RBP Analysis (`05_rbp_analysis.R`)
- Maps SNPs to RNA-binding proteins
- Analyzes expression data from knockout experiments
- Performs statistical testing

### 7. Cancer Knockout (`06_cancer_ko_analysis.R`)
- Analyzes cancer cell line knockout data
- Creates differential expression boxplots
- Performs statistical significance testing

## 📋 Data Requirements

The following data files are expected in the project root (one level up from Git folder):

### Required Files
- `significant_SNP.RData` - Significant SNP results
- `NMD_gene_list.RData` - NMD gene annotations
- `plot_man_anno.csv` - Manhattan plot annotation data
- `sig_snp_gene_id_mapped_res.csv` - SNP-gene mapping results

### Optional Files (for full functionality)
- `gwas_catalog_v1.0.2-associations_e112_r2024-09-22.tsv` - GWAS catalog data
- `pathway/pathways_info.csv` - Pathway information
- `RBP_map/RBP2GO_dataset.txt` - RBP dataset
- `Cancer_KO/level5_beta_trt_xpr_n142901x12328_nmd_genes.RData` - Cancer KO data

## 📈 Output Description

### Plots (`output/plots/`)
- `tissue_distribution.pdf/jpeg` - Tissue distribution bar plot
- `manhattan_plot.pdf/jpeg` - Annotated Manhattan plot
- `manhattan_plot_basic.pdf/jpeg` - Basic Manhattan plot
- `tissue_gwas_diseases.pdf/jpeg` - GWAS-tissue distribution
- `oxr1_expression_boxplot_*.pdf/jpeg` - Expression boxplots
- `cancer_ko_boxplot_*.pdf/jpeg` - Cancer knockout boxplots

### Tables (`output/tables/`)
- `tissue_distribution_summary.csv` - Tissue distribution statistics
- `manhattan_plot_data.csv` - Manhattan plot coordinates
- `annotation_summary.csv` - SNP annotation summary
- `gwas_overlap_summary.csv` - GWAS overlap statistics

### Analysis Results
- `output/GWAS_overlap/` - GWAS catalog overlap results

## 🔧 Configuration Options

### Modifying Analysis Parameters

Edit parameters at the top of individual scripts:

```R
# In 05_gwas_overlap.R
window_size <- 1e5  # Change window size around SNPs

# In 04_go_analysis.R
pvalue_cutoff <- 0.05  # Change GO analysis p-value threshold

# In 00_data_loading.R
# Modify color schemes and plotting parameters
```

### Adding New Analyses

1. Create new script in `scripts/` folder
2. Follow naming convention: `XX_analysis_name.R`
3. Add to `run_all_analyses.R` pipeline

## 🛠 Troubleshooting

### Common Issues

1. **biomaRt connection errors**
   ```R
   # Try different mirrors
   mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", 
                      mirror="useast")
   ```

2. **Memory issues with large datasets**
   - Increase R memory limit
   - Process data in chunks
   - Use data.table for large files

3. **Missing data files**
   - Check file paths in error messages
   - Verify data file locations
   - Some analyses will skip gracefully if data unavailable

### Error Reporting

When reporting issues, please include:
- R version and platform
- Error message and stack trace
- Which analysis script failed
- Available data files

## 📚 Dependencies

### R Packages

**CRAN Packages:**
- tidyverse, ggplot2, ggrepel
- viridis, RColorBrewer, gridExtra
- data.table, reshape2, broom
- igraph, coloc

**Bioconductor Packages:**
- biomaRt, GenomicRanges
- clusterProfiler, org.Hs.eg.db
- AnnotationDbi, enrichplot

### External Resources
- Ensembl biomaRt (for gene annotation)
- GWAS Catalog (for disease associations)
- dbSNP (for SNP annotations)

## 📄 Citation

If you use this pipeline in your research, please cite:

```

```

## 🤝 Contributing

This pipeline was developed for a specific research project. For modifications or extensions:

1. Fork the repository
2. Create analysis branches
3. Submit pull requests with detailed descriptions
4. Include test data when possible

## 📞 Contact

For questions or issues:
- Create an issue in this repository
- Contact: 

## 📜 License

This code is provided for research purposes. Please see the license file for details.

---

**Last Updated:** January 2025  
**Version:** 1.0.0  
**Compatibility:** R ≥ 4.0.0 
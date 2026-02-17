# Trans-acting NMD Modulator Mapping Pipeline

This repository contains the computational pipeline for paper ***Identification of tissue-specific NMD modulators via novel trans-acting molecular QTL mapping***. The pipeline is organized into modular scripts for reproducibility and transparency.

## Repository Structure

```
Trans-acting-NMD-QTL/
├── README.md                        # This file
├── requirements.R                   # R package dependencies
├── setup_data.R                     # Data setup helper
├── derive_deviation_scores.R        # Deviation score calculation
├── scripts/                         # Analysis modules
│   ├── 00_data_loading.R            # Data loading and common functions
│   ├── 01_tissue_distribution.R     # Tissue-wise SNP distribution analysis
│   ├── 02_snp_annotation.R          # SNP annotation with biomaRt
│   ├── 03_manhattan_plot.R          # Manhattan plot generation
│   ├── 04_gwas_overlap.R           # GWAS catalog overlap analysis
│   ├── 05_rbp_analysis.R           # RNA-binding protein analysis
│   ├── 06_cancer_ko_analysis.R     # Cancer knockout analysis
│   ├── 07_generate_deviation_data.R # NMD deviation data generation
│   ├── 08_outlier_deviation_analysis.R # Outlier deviation analysis
│   ├── 09_oxr1_knockdown_analysis.R   # OXR1 knockdown analysis (Fig 5A)
│   └── 10_oxr1_deviation_regression.R # OXR1 deviation regression (Fig 5B)
├── data/                            # Data files (see Data Requirements)
└── output/                          # Generated outputs
    ├── plots/                       # Publication-ready figures
    └── tables/                      # Summary tables and results
```

## Quick Start

### Prerequisites

- R (version >= 4.0.0)
- Internet connection (for biomaRt queries)
- Sufficient RAM (>8GB recommended for large datasets)

### Installation

1. Clone or download this repository
2. Install required R packages:
   ```R
   source("requirements.R")
   ```

### Running the Analysis

Scripts are designed to be run individually from the `scripts/` directory:

```R
setwd("scripts/")
source("01_tissue_distribution.R")
source("02_snp_annotation.R")
# ... etc
```

Scripts `07`-`10` require external data (GTEx, Cancer KO) not included in this repository. Set the `NMD_PROJECT_DIR` environment variable to point to the project root containing these data directories, or adjust the `BASE_DIR` variable at the top of each script.

```R
Sys.setenv(NMD_PROJECT_DIR = "/path/to/your/project/root")
```

## Analysis Modules

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

### 8. Deviation Data Generation (`07_generate_deviation_data.R`)
- Generates deviation data for all 81 overlapped genes between NMD-QTL mapping and GTEx
- Computes RINT-based deviation scores for NMD transcripts
- Creates per-gene deviation files used by downstream analyses
- Requires: GTEx expression data, Manuscript/Tables.xlsx

### 9. Outlier Deviation Analysis (`08_outlier_deviation_analysis.R`)
- Analyzes relationship between NMD factor expression and outlier deviations
- Creates outlier proportion line plots across expression deciles
- Generates scatter plots and correlation bar plots
- Produces heatmaps of outlier proportions for significant genes

### 10. OXR1 Knockdown Analysis (`09_oxr1_knockdown_analysis.R`) - Figure 5A
- Analyzes NMD transcript expression changes upon OXR1 knockdown
- Compares NMD-targeted vs non-NMD-targeted transcripts
- Uses signed log fold change to preserve directionality
- Generates boxplots, violin plots, and density plots

### 11. OXR1 Deviation Regression (`10_oxr1_deviation_regression.R`) - Figure 5B
- Analyzes OXR1 expression vs average absolute deviation in Brain Cerebellum
- Performs linear regression analysis
- Generates scatter plots with regression lines

## Data Requirements

### Included in Repository (`data/`)
- `significant_SNP.RData` - Significant SNP results
- `NMD_gene_list.RData` - NMD gene annotations
- `plot_man_anno.csv` - Manhattan plot annotation data
- `sig_snp_gene_id_mapped_res.csv` - SNP-gene mapping results
- `sig_snp_gwas_catalog.csv` - GWAS catalog overlap results
- `oxr1_tx_count_data_nmd_anno.csv` - OXR1 knockdown transcript data
- `near_gens_NMD_SNP.txt` - Nearest genes to NMD SNPs

### External Data (not included)
- `RBP2GO_dataset.txt` - RBP dataset (required for script 05; download from [RBP2GO](https://rbp2go.dkfz.de/))
- GTEx gene expression data v8 TPM (required for script 07)
- GTEx NMD transcript expression data (required for script 07)
- GTEx sample attributes (required for script 07)
- `Manuscript/Tables.xlsx` - Supplementary tables (required for script 07)
- `Cancer_KO/` - Cancer knockout deviation data (required for scripts 08, 10)

## Output Description

### Plots (`output/plots/`)
- `tissue_distribution.pdf/jpeg` - Tissue distribution bar plot
- `manhattan_plot.pdf/jpeg` - Annotated Manhattan plot
- `tissue_gwas_diseases.pdf/jpeg` - GWAS-tissue distribution
- `Fig5a_oxr1_*.pdf/jpeg` - OXR1 knockdown expression plots
- `Fig5b_OXR1_vs_deviation_*.pdf/jpeg` - OXR1 deviation regression plots

### Tables (`output/tables/`)
- `tissue_distribution_summary.csv` - Tissue distribution statistics
- `manhattan_plot_data.csv` - Manhattan plot coordinates
- `annotation_summary.csv` - SNP annotation summary
- `gwas_overlap_summary.csv` - GWAS overlap statistics

## Dependencies

### R Packages

**CRAN Packages:**
- tidyverse, ggplot2, ggrepel, ggpubr
- viridis, RColorBrewer, gridExtra, patchwork
- data.table, reshape2, broom, readxl, scales

**Bioconductor Packages:**
- biomaRt, GenomicRanges
- clusterProfiler, org.Hs.eg.db
- AnnotationDbi, enrichplot

### External Resources
- Ensembl biomaRt (for gene annotation)
- GWAS Catalog (for disease associations)
- dbSNP (for SNP annotations)
- GTEx Portal (for expression data)

## Citation

If you use this pipeline in your research, please cite:

```
[Citation to be added upon publication]
```

## License

This code is provided for research purposes. Please see the license file for details.

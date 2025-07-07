# R Package Requirements for NMD Gene Analysis
# Install and load required packages

# CRAN packages
cran_packages <- c(
  "tidyverse", "dplyr", "readr", "stringr", "ggplot2", "ggrepel",
  "viridis", "gridExtra", "RColorBrewer", "ComplexUpset", "ggsci",
  "grid", "patchwork", "data.table", "reshape2", "broom",
  "qqman", "GenomicRanges", "tximport", "DESeq2", "RNOmni",
  "igraph", "coloc"
)

# Bioconductor packages
bioc_packages <- c(
  "biomaRt", "AnnotationDbi", "org.Hs.eg.db", "EnsDb.Hsapiens.v86",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38", "clusterProfiler", "enrichplot",
  "cmapR", "rtracklayer"
)

# Function to install packages if not already installed
install_if_missing <- function(packages, type = "CRAN") {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  if(length(new_packages) > 0) {
    if(type == "CRAN") {
      install.packages(new_packages, dependencies = TRUE)
    } else if(type == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(new_packages)
    }
  }
}

# Install packages
install_if_missing(cran_packages, "CRAN")
install_if_missing(bioc_packages, "Bioconductor")

# Load all packages
library_packages <- function(packages) {
  for(pkg in packages) {
    library(pkg, character.only = TRUE)
  }
}

library_packages(cran_packages)
library_packages(bioc_packages)

cat("All required packages loaded successfully!\n") 
# ==============================================================================
# Data Loading and Setup for NMD Gene Analysis
# This script loads all common data files used across multiple analyses
# ==============================================================================

# Load required packages
source("../requirements.R")

# Set working directory (adjust as needed)
# setwd("path/to/your/project")

# ==============================================================================
# Load Core Data Files
# ==============================================================================

cat("Loading core data files...\n")

# Load significant SNPs data
if(file.exists("../data/significant_SNP.RData")) {
  load("../data/significant_SNP.RData")
  sig_SNP <- R
  rm(R)
  cat("✓ Significant SNPs loaded\n")
} else {
  stop("significant_SNP.RData file not found!")
}

# Load NMD gene list
if(file.exists("../data/NMD_gene_list.RData")) {
  load("../data/NMD_gene_list.RData")
  cat("✓ NMD gene list loaded\n")
} else {
  cat("⚠ NMD_gene_list.RData not found\n")
}

# Load Manhattan plot annotation data
if(file.exists("../data/plot_man_anno.csv")) {
  plot_man_anno <- read.csv("../data/plot_man_anno.csv", header = TRUE)[, -1]
  cat("✓ Manhattan plot annotations loaded\n")
} else {
  cat("⚠ plot_man_anno.csv not found\n")
}

# Load significant SNP gene mapping results
if(file.exists("../data/sig_snp_gene_id_mapped_res.csv")) {
  sig_snp_gene_id_mapped_res <- read.csv("../data/sig_snp_gene_id_mapped_res.csv", header = TRUE)
  cat("✓ SNP gene mapping results loaded\n")
} else {
  cat("⚠ sig_snp_gene_id_mapped_res.csv not found\n")
}

# Load NMD related genes and proteins
if(file.exists("../data/NMD_Related_Genes_and_Proteins.csv")) {
  nmd_related_genes <- read.csv("../data/NMD_Related_Genes_and_Proteins.csv", header = TRUE)
  cat("✓ NMD related genes and proteins loaded\n")
} else {
  cat("⚠ NMD_Related_Genes_and_Proteins.csv not found\n")
  nmd_related_genes <- NULL
}

# Load GWAS catalog overlap data
if(file.exists("../../GWAS_overlap/GWAS_Catalog/sig_snp_gwas_catalog.csv")) {
  sig_snp_gwas_catalog <- read.csv("../../GWAS_overlap/GWAS_Catalog/sig_snp_gwas_catalog.csv", header = TRUE)
  cat("✓ GWAS catalog overlap data loaded\n")
} else {
  cat("⚠ sig_snp_gwas_catalog.csv not found\n")
  sig_snp_gwas_catalog <- NULL
}

# ==============================================================================
# Common Functions
# ==============================================================================

# Function to clean and prepare SNP data
prepare_snp_data <- function(snp_data) {
  snp_data %>%
    separate(pos, c("CHR", "pos"), sep = ":") %>%
    mutate(CHR = gsub("chr", "", CHR), 
           pos = as.numeric(pos)) %>%
    mutate(SNP = variant_id, 
           start = pos, 
           end = pos, 
           variant_pos = paste0(CHR, ":", pos)) %>%
    dplyr::select(SNP, CHR, start, end, variant_pos, p_val, n, tissue)
}

# Function to create custom color palette
get_custom_colors <- function(n = 16) {
  c("#1b9e77", "#b2df8a", "#7570b3", "#e7298a", "#66a61e", 
    "#e6ab02", "#a6761d", "#fdbf6f", "#1f78b4", "#d95f02",
    "#fb9a99", "#cab2d6", "#666666", "#ffff99", "#a6cee3", 
    "#b15928")[1:n]
}

# Function to save plots with consistent formatting
save_plot <- function(plot, filename, width = 9, height = 6, dpi = 450) {
  output_path <- file.path("../output/plots", filename)
  
  # Save as both PDF and JPEG
  ggsave(paste0(tools::file_path_sans_ext(output_path), ".pdf"), 
         plot = plot, width = width, height = height, dpi = 300)
  ggsave(paste0(tools::file_path_sans_ext(output_path), ".jpeg"), 
         plot = plot, width = width, height = height, dpi = dpi)
  
  cat("Plot saved:", filename, "\n")
}

# Function to clean tissue names
clean_tissue_names <- function(tissue_names) {
  tissue_names %>%
    gsub(" - | ", "_", .) %>%
    gsub("[()]", "", .)
}

cat("✓ Data loading and setup complete!\n")
cat("Available objects: sig_SNP, plot_man_anno, sig_snp_gene_id_mapped_res, nmd_related_genes, sig_snp_gwas_catalog\n") 
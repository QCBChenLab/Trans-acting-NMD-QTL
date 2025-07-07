# ==============================================================================
# RNA-Binding Protein (RBP) Analysis - Focus on Downregulated Transcripts
# This script analyzes RBP binding and OXR1 downregulated transcript expression
# 
# Key features:
# - Maps significant SNPs to RBP genes
# - Analyzes expression of downregulated transcripts only
# - Performs statistical tests on NMD vs Non-NMD transcripts
# - Follows Manuscript.Rmd approach focusing on downregulation
# ==============================================================================

# Load data and functions
source("00_data_loading.R")

cat("Starting RBP binding analysis...\n")

# ==============================================================================
# Load RBP and Related Data
# ==============================================================================

# Function to load RBP dataset
load_rbp_data <- function() {
  # Try multiple potential file locations
  potential_files <- c(
    "../../RBP_map/RBP2GO_dataset.txt",
    "../data/RBP2GO_dataset.txt",
    "../RBP_map/RBP2GO_dataset.txt"
  )
  
  for(rbp_file in potential_files) {
    if(file.exists(rbp_file)) {
      rbp_data <- read.table(rbp_file, header = TRUE)
      cat("✓ RBP2GO dataset loaded from:", rbp_file, "\n")
      cat("  Total entries:", nrow(rbp_data), "\n")
      return(rbp_data)
    }
  }
  
  cat("⚠ RBP2GO dataset not found in any expected location\n")
  return(NULL)
}

# Function to load SNP gene mapping data
load_snp_gene_mapping <- function() {
  # Try multiple potential file locations
  potential_files <- c(
    "../../protein_map/snp_nexus_results/near_gens_NMD_SNP.txt",
    "../data/near_gens_NMD_SNP.txt",
    "../protein_map/snp_nexus_results/near_gens_NMD_SNP.txt"
  )
  
  for(mapping_file in potential_files) {
    if(file.exists(mapping_file)) {
      sig_SNP_anno_res <- read.table(mapping_file, header = TRUE, sep = "\t")
      cat("✓ SNP gene mapping loaded from:", mapping_file, "\n")
      cat("  Total entries:", nrow(sig_SNP_anno_res), "\n")
      return(sig_SNP_anno_res)
    }
  }
  
  cat("⚠ SNP gene mapping file not found in any expected location\n")
  return(NULL)
}

# Load data
rbp_data <- load_rbp_data()
sig_SNP_anno_res <- load_snp_gene_mapping()

# ==============================================================================
# Map SNPs to RBP Genes
# ==============================================================================

# Function to extract genes from SNP mapping
extract_snp_genes <- function(snp_mapping_data) {
  if(is.null(snp_mapping_data)) {
    return(list(overlapped = character(), upstream = character(), downstream = character()))
  }
  
  # Extract different types of genes
  sig_SNP_ovlp_gene <- snp_mapping_data %>% 
    select(Variation.ID, Overlapped.Gene) %>%
    filter(Overlapped.Gene != "None") %>%
    pull(Overlapped.Gene)
  
  sig_SNP_up_gene <- snp_mapping_data %>% 
    select(Variation.ID, Nearest.Upstream.Gene) %>%
    filter(Nearest.Upstream.Gene != "None") %>%
    pull(Nearest.Upstream.Gene)
  
  sig_SNP_dw_gene <- snp_mapping_data %>% 
    select(Variation.ID, Nearest.Downstream.Gene) %>%
    filter(Nearest.Downstream.Gene != "None") %>%
    pull(Nearest.Downstream.Gene)
  
  cat("✓ Extracted genes: Overlapped =", length(sig_SNP_ovlp_gene), 
      "Upstream =", length(sig_SNP_up_gene), 
      "Downstream =", length(sig_SNP_dw_gene), "\n")
  
  return(list(
    overlapped = sig_SNP_ovlp_gene,
    upstream = sig_SNP_up_gene,
    downstream = sig_SNP_dw_gene
  ))
}

# Function to map SNP genes to RBPs
map_snp_genes_to_rbp <- function(snp_genes_list, rbp_data) {
  if(is.null(rbp_data)) {
    cat("⚠ No RBP data available\n")
    return(data.frame())
  }
  
  # Map each gene type to RBPs
  sig_SNP_rbp_ovlp <- rbp_data %>% 
    filter(Gene_Name %in% snp_genes_list$overlapped) %>%
    mutate(Gene_Type = "Overlapped")
  
  sig_SNP_rbp_up <- rbp_data %>%
    filter(Gene_Name %in% snp_genes_list$upstream) %>%
    mutate(Gene_Type = "Upstream")
  
  sig_SNP_rbp_dw <- rbp_data %>%
    filter(Gene_Name %in% snp_genes_list$downstream) %>%
    mutate(Gene_Type = "Downstream")
  
  # Combine all RBP mappings
  sig_SNP_rbp <- rbind(sig_SNP_rbp_ovlp, sig_SNP_rbp_up, sig_SNP_rbp_dw)
  
  cat("✓ Mapped", nrow(sig_SNP_rbp), "SNP genes to RBPs\n")
  return(sig_SNP_rbp)
}

# Extract SNP genes and map to RBPs
snp_genes <- extract_snp_genes(sig_SNP_anno_res)
sig_SNP_rbp <- map_snp_genes_to_rbp(snp_genes, rbp_data)

# Save RBP mapping results
if(nrow(sig_SNP_rbp) > 0) {
  write.csv(sig_SNP_rbp, "../output/tables/sig_SNP_rbp_mapping.csv", row.names = FALSE)
  cat("✓ RBP mapping results saved\n")
}

# ==============================================================================
# Load Expression Data for OXR1 - Focus on Downregulated
# ==============================================================================

# Function to load OXR1 expression data
load_oxr1_expression <- function() {
  # Try multiple potential file locations
  potential_files <- c(
    "../../RBP_map/Gene_Expression/oxr1_tx_count_data_nmd_anno.csv",
    "../data/oxr1_tx_count_data_nmd_anno.csv",
    "../RBP_map/Gene_Expression/oxr1_tx_count_data_nmd_anno.csv"
  )
  
  for(exp_file in potential_files) {
    if(file.exists(exp_file)) {
      oxr1_exp_data <- read.csv(exp_file, header = TRUE)[, -1]
      cat("✓ OXR1 expression data loaded from:", exp_file, "\n")
      cat("  Total transcripts:", nrow(oxr1_exp_data), "\n")
      return(oxr1_exp_data)
    }
  }
  
  cat("⚠ OXR1 expression data not found in any expected location\n")
  return(NULL)
}

# Function to load downregulated transcripts only
load_downregulated_transcripts <- function() {
  # Try multiple potential file locations
  potential_files <- c(
    "../../RBP_map/Gene_Expression/GTEx_mutation/down_tx.txt",
    "../data/down_tx.txt",
    "../RBP_map/Gene_Expression/down_tx.txt"
  )
  
  for(down_file in potential_files) {
    if(file.exists(down_file)) {
      down_tx <- read.table(down_file, header = FALSE)$V1
      cat("✓ Downregulated transcripts loaded from:", down_file, "\n")
      cat("  Number of downregulated transcripts:", length(down_tx), "\n")
      return(down_tx)
    }
  }
  
  cat("⚠ Downregulated transcripts file not found\n")
  return(character())
}

# Load expression data
oxr1_tx_count_data_nmd_anno <- load_oxr1_expression()
down_tx <- load_downregulated_transcripts()
de_transcripts <- list(down = down_tx)  # Keep structure for compatibility

# ==============================================================================
# Expression Analysis - Focus on Downregulated Transcripts
# ==============================================================================

# Function to create expression boxplot for downregulated transcripts
create_downregulated_expression_boxplot <- function(expression_data, down_tx) {
  if(is.null(expression_data) || length(down_tx) == 0) {
    cat("⚠ No expression data or downregulated transcripts available\n")
    return(NULL)
  }
  
  # Filter for downregulated transcripts only
  plot_data <- expression_data %>%
    filter(Transcript_ID %in% down_tx)
  
  cat("✓ Analyzing", nrow(plot_data), "downregulated transcripts\n")
  
  # Convert to long format
  plot_data_long <- plot_data %>%
    pivot_longer(cols = starts_with("si"),  
                 names_to = "Condition", 
                 values_to = "Expression") %>% 
    mutate(Condition = factor(Condition, 
                              levels = c("siCon_NT", "siCon_R0h", 
                                         "siOXR1_NT", "siOXR1_R0h")), 
           Expression = as.numeric(Expression))
  
  # Create boxplot following Manuscript.Rmd style
  boxplot_p <- ggplot(plot_data_long, aes(x = Condition, y = Expression, fill = NMD)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    scale_fill_manual(values = c("Non_NMD" = "#1f77b4", "NMD" = "#ff7f0e")) +
    theme_bw() +
    labs(title = "Downregulated Transcripts Expression Level", 
         x = "Condition", 
         y = "Transcripts Expression Level", 
         fill = "NMD Status") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14), 
          axis.title.y = element_text(size = 14), 
          plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
  
  return(boxplot_p)
}

# Generate expression boxplot for downregulated transcripts only
if(!is.null(oxr1_tx_count_data_nmd_anno) && length(de_transcripts$down) > 0) {
  exp_boxplot_down <- create_downregulated_expression_boxplot(
    oxr1_tx_count_data_nmd_anno, 
    de_transcripts$down
  )
  
  if(!is.null(exp_boxplot_down)) {
    save_plot(exp_boxplot_down, "oxr1_tx_exp_box_down", width = 9, height = 6)
    cat("✓ Downregulated transcripts expression plot saved\n")
  }
}

# ==============================================================================
# Statistical Analysis for Downregulated Transcripts
# ==============================================================================

# Function to perform statistical tests on downregulated expression data
perform_downregulated_stats <- function(expression_data, down_tx) {
  if(is.null(expression_data) || length(down_tx) == 0) {
    cat("⚠ No data available for statistical analysis\n")
    return(NULL)
  }
  
  # Prepare data for statistical testing - downregulated only
  stats_data <- expression_data %>%
    filter(Transcript_ID %in% down_tx) %>%
    pivot_longer(cols = starts_with("si"), 
                 names_to = "Condition", 
                 values_to = "Expression") %>%
    mutate(Group = case_when(
             grepl("siCon", Condition) ~ "Control",
             grepl("siOXR1", Condition) ~ "OXR1_depleted",
             TRUE ~ NA_character_
           )) %>%
    filter(!is.na(Group), NMD %in% c("NMD", "Non_NMD"))
  
  # Perform Wilcoxon tests comparing NMD vs Non-NMD
  wilcox_results <- stats_data %>%
    group_by(Condition) %>%
    group_modify(~ {
      if(length(unique(.x$NMD)) < 2) {
        return(data.frame(statistic = NA, p.value = NA))
      }
      wt <- wilcox.test(Expression ~ NMD, data = .x, exact = FALSE)
      data.frame(statistic = wt$statistic, p.value = wt$p.value)
    }) %>%
    ungroup() %>%
    mutate(p.value = round(p.value, 5))
  
  cat("✓ Statistical tests on downregulated transcripts completed\n")
  return(wilcox_results)
}

# Perform statistical analysis on downregulated transcripts
if(!is.null(oxr1_tx_count_data_nmd_anno) && length(down_tx) > 0) {
  expression_stats <- perform_downregulated_stats(oxr1_tx_count_data_nmd_anno, down_tx)
  
  if(!is.null(expression_stats)) {
    write.csv(expression_stats, "../output/tables/downregulated_expression_stats.csv", row.names = FALSE)
    cat("\nWilcoxon test results (NMD vs Non-NMD in downregulated transcripts):\n")
    print(expression_stats)
  }
}

# ==============================================================================
# RBP Analysis Summary - Focus on Downregulated
# ==============================================================================

# Create RBP analysis summary
create_rbp_summary <- function() {
  # Count NMD vs Non-NMD in downregulated transcripts
  nmd_counts <- if(!is.null(oxr1_tx_count_data_nmd_anno) && length(down_tx) > 0) {
    oxr1_tx_count_data_nmd_anno %>%
      filter(Transcript_ID %in% down_tx) %>%
      group_by(NMD) %>%
      summarise(n = n()) %>%
      deframe()
  } else {
    list(NMD = 0, Non_NMD = 0)
  }
  
  summary_data <- data.frame(
    Metric = c(
      "Total RBP entries",
      "SNP-associated RBP genes",
      "Overlapped gene RBPs",
      "Upstream gene RBPs", 
      "Downstream gene RBPs",
      "Total expression data transcripts",
      "Downregulated transcripts analyzed",
      "Downregulated NMD transcripts",
      "Downregulated Non-NMD transcripts"
    ),
    Count = c(
      ifelse(!is.null(rbp_data), nrow(rbp_data), 0),
      ifelse(exists("sig_SNP_rbp"), nrow(sig_SNP_rbp), 0),
      ifelse(exists("sig_SNP_rbp"), sum(sig_SNP_rbp$Gene_Type == "Overlapped"), 0),
      ifelse(exists("sig_SNP_rbp"), sum(sig_SNP_rbp$Gene_Type == "Upstream"), 0),
      ifelse(exists("sig_SNP_rbp"), sum(sig_SNP_rbp$Gene_Type == "Downstream"), 0),
      ifelse(!is.null(oxr1_tx_count_data_nmd_anno), nrow(oxr1_tx_count_data_nmd_anno), 0),
      length(down_tx),
      ifelse("NMD" %in% names(nmd_counts), nmd_counts["NMD"], 0),
      ifelse("Non_NMD" %in% names(nmd_counts), nmd_counts["Non_NMD"], 0)
    )
  )
  
  write.csv(summary_data, "../output/tables/rbp_analysis_summary.csv", row.names = FALSE)
  return(summary_data)
}

# Generate and display summary
rbp_summary <- create_rbp_summary()
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("RBP Analysis Summary:\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
print(rbp_summary)
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("\n✓ RBP binding analysis complete!\n")
cat("✓ Analysis focused on downregulated transcripts only\n")
cat("✓ Results saved to output/tables/\n")
cat("✓ Plots saved to output/plots/\n") 
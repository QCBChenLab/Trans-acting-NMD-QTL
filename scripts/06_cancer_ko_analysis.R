# ==============================================================================
# Cancer Knockout Analysis
# This script analyzes cancer cell line knockout data and creates related plots
# ==============================================================================

# Load data and functions
source("00_data_loading.R")

cat("Starting cancer knockout analysis...\n")

# ==============================================================================
# Load Cancer KO Data
# ==============================================================================

# Function to load cancer knockout data
load_cancer_ko_data <- function() {
  ko_file <- "../Cancer_KO/level5_beta_trt_xpr_n142901x12328_nmd_genes.RData"
  sig_info_file <- "../Cancer_KO/siginfo_beta.txt"
  
  # Try to load pre-processed data
  if(file.exists(ko_file)) {
    load(ko_file)
    cat("✓ Pre-processed cancer KO data loaded\n")
    
    # Load signature information
    if(file.exists(sig_info_file)) {
      sig_info <- read_tsv(sig_info_file) %>%
        select(sig_id, cell_iname, cmap_name)
      cat("✓ Signature information loaded\n")
    } else {
      cat("⚠ Signature information not found\n")
      sig_info <- NULL
    }
    
    return(list(data_matrix = data_matrix, sig_info = sig_info))
    
  } else {
    cat("⚠ Cancer KO data file not found\n")
    return(NULL)
  }
}

# Function to load RBP genes from SNP analysis
load_rbp_genes <- function() {
  rbp_file <- "../protein_map/sig_SNP_anno_res_genes.csv"
  
  if(file.exists(rbp_file)) {
    rbp_genes_df <- read_csv(rbp_file)
    rbp_genes <- unique(rbp_genes_df$Gene)
    cat("✓ RBP genes loaded:", length(rbp_genes), "genes\n")
    return(rbp_genes)
  } else {
    cat("⚠ RBP genes file not found, using dummy genes\n")
    return(c("ACVR1B", "VDR", "IL1RAP", "SDC1", "NRG1", "WNK3"))
  }
}

# Load data
cancer_ko_data <- load_cancer_ko_data()
rbp_genes <- load_rbp_genes()

# ==============================================================================
# Cancer KO Analysis
# ==============================================================================

# Function to analyze cancer knockout data
analyze_cancer_ko <- function(ko_data, rbp_genes) {
  if(is.null(ko_data) || is.null(ko_data$data_matrix)) {
    cat("⚠ No cancer KO data available\n")
    return(NULL)
  }
  
  data_matrix <- ko_data$data_matrix
  sig_info <- ko_data$sig_info
  
  if(is.null(sig_info)) {
    cat("⚠ No signature information available\n")
    return(NULL)
  }
  
  # Create sample information dataframe
  sample_info <- data.frame(sig_id = rownames(data_matrix)) %>%
    left_join(sig_info, by = "sig_id") %>%
    column_to_rownames("sig_id") %>% 
    filter(!is.na(cmap_name))
  
  # Find overlap between RBP genes and available treatments
  rbp_gene_overlap <- intersect(sample_info$cmap_name, rbp_genes)
  cat("✓ Found", length(rbp_gene_overlap), "RBP genes in cancer KO data\n")
  
  if(length(rbp_gene_overlap) == 0) {
    cat("⚠ No overlap between RBP genes and cancer KO data\n")
    return(NULL)
  }
  
  # Compute mean NMD expression per RBP
  boxplot_data <- list()
  
  for(rbp in rbp_gene_overlap) {
    ids <- rownames(sample_info)[sample_info$cmap_name == rbp]
    if(length(ids) > 0) {
      boxplot_data[[rbp]] <- colMeans(data_matrix[ids, , drop = FALSE], na.rm = TRUE)
    }
  }
  
  # Remove empty entries
  boxplot_data <- Filter(Negate(is.null), boxplot_data)
  
  return(boxplot_data)
}

# Perform cancer KO analysis
ko_analysis_results <- analyze_cancer_ko(cancer_ko_data, rbp_genes)

# ==============================================================================
# Statistical Testing
# ==============================================================================

# Function to perform t-tests on knockout data
perform_ko_statistics <- function(boxplot_data) {
  if(is.null(boxplot_data) || length(boxplot_data) == 0) {
    return(NULL)
  }
  
  # Perform t-tests
  raw_pvals <- sapply(boxplot_data, function(values) {
    if(length(values) > 1) {
      t.test(values, mu = 0)$p.value
    } else {
      NA
    }
  })
  
  # FDR correction
  fdr_corrected_pvals <- p.adjust(raw_pvals, method = "fdr")
  
  # Create significance labels
  signif_label <- function(p) {
    if(is.na(p)) return("")
    if(p < 1e-35) return("****")
    if(p < 1e-25) return("***")
    if(p < 1e-15) return("**")
    if(p < 1e-5) return("*")
    return("")
  }
  
  signif_labels <- sapply(fdr_corrected_pvals, signif_label)
  
  # Create results dataframe
  stats_results <- data.frame(
    Gene = names(raw_pvals),
    Raw_pvalue = raw_pvals,
    FDR_pvalue = fdr_corrected_pvals,
    Significance = signif_labels,
    stringsAsFactors = FALSE
  )
  
  return(stats_results)
}

# Perform statistical testing
if(!is.null(ko_analysis_results)) {
  ko_stats <- perform_ko_statistics(ko_analysis_results)
  
  if(!is.null(ko_stats)) {
    write.csv(ko_stats, "../output/tables/cancer_ko_statistics.csv", row.names = FALSE)
    cat("✓ Statistical results saved\n")
  }
}

# ==============================================================================
# Create Boxplots
# ==============================================================================

# Function to create cancer KO boxplots
create_ko_boxplots <- function(boxplot_data, stats_results = NULL) {
  if(is.null(boxplot_data) || length(boxplot_data) == 0) {
    return(NULL)
  }
  
  # Convert to plotting format
  df_plot <- stack(boxplot_data)
  colnames(df_plot) <- c("expression", "rbp_gene")
  
  # Add significance labels if available
  if(!is.null(stats_results)) {
    # Calculate label positions
    label_positions <- tapply(df_plot$expression, df_plot$rbp_gene, max, na.rm = TRUE) + 0.1
    
    signif_df <- data.frame(
      rbp_gene = names(label_positions),
      label = stats_results$Significance[match(names(label_positions), stats_results$Gene)],
      y_pos = label_positions
    )
  } else {
    signif_df <- data.frame(rbp_gene = character(), label = character(), y_pos = numeric())
  }
  
  # Create full boxplot
  box_p_full <- ggplot(df_plot, aes(x = rbp_gene, y = expression)) +
    geom_boxplot(fill = "lightblue", outlier.shape = 1, outlier.size = 1.5) +
    geom_text(data = signif_df, 
              aes(x = rbp_gene, y = y_pos, label = label), 
              hjust = 0.3, vjust = 0.8, size = 5, angle = 90, color = "red") +
    geom_hline(yintercept = 0, color = "orange", linetype = "dashed") +
    coord_cartesian(ylim = c(-1.2, 1.6)) + 
    theme_bw() +
    labs(title = "NMD Gene Differential Expression Across RBP Knockouts",
         x = "RBP Gene",
         y = "Mean Differential Expression of NMD Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create selected genes plot
  selected_genes <- c("ACVR1B", "VDR", "IL1RAP", "SDC1", "NRG1", "WNK3")
  available_selected <- intersect(selected_genes, unique(df_plot$rbp_gene))
  
  if(length(available_selected) > 0) {
    df_plot_sel <- df_plot %>% filter(rbp_gene %in% available_selected)
    signif_df_sel <- signif_df %>% filter(rbp_gene %in% available_selected)
    
    box_p_selected <- ggplot(df_plot_sel, aes(x = rbp_gene, y = expression)) +
      geom_boxplot(width = 0.3, fill = "lightblue", outlier.shape = 1, outlier.size = 1.5) +
      geom_text(data = signif_df_sel, 
                aes(x = rbp_gene, y = y_pos, label = label), 
                hjust = 0.3, vjust = 0.8, size = 5, angle = 90, color = "red") +
      geom_hline(yintercept = 0, color = "orange", linetype = "dashed") +
      coord_cartesian(ylim = c(-1.2, 1.6)) + 
      theme_bw() +
      labs(title = "NMD Gene Differential Expression (Selected RBP Genes)",
           x = "RBP Gene",
           y = "Mean Differential Expression of NMD Genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    box_p_selected <- NULL
  }
  
  return(list(full = box_p_full, selected = box_p_selected))
}

# Create boxplots
if(!is.null(ko_analysis_results)) {
  ko_plots <- create_ko_boxplots(ko_analysis_results, ko_stats)
  
  if(!is.null(ko_plots$full)) {
    save_plot(ko_plots$full, "cancer_ko_boxplot_full", width = 12, height = 6)
  }
  
  if(!is.null(ko_plots$selected)) {
    save_plot(ko_plots$selected, "cancer_ko_boxplot_selected", width = 9, height = 6)
  }
}

# ==============================================================================
# Create Summary
# ==============================================================================

# Function to create cancer KO summary
create_ko_summary <- function(boxplot_data, stats_results) {
  summary_data <- data.frame(
    Metric = c(
      "Total RBP genes tested",
      "Significantly downregulated",
      "Significantly upregulated", 
      "No significant effect"
    ),
    Count = c(
      ifelse(!is.null(boxplot_data), length(boxplot_data), 0),
      ifelse(!is.null(stats_results), sum(stats_results$FDR_pvalue < 0.05 & 
                                          sapply(boxplot_data[stats_results$Gene], mean) < 0, na.rm = TRUE), 0),
      ifelse(!is.null(stats_results), sum(stats_results$FDR_pvalue < 0.05 & 
                                          sapply(boxplot_data[stats_results$Gene], mean) > 0, na.rm = TRUE), 0),
      ifelse(!is.null(stats_results), sum(stats_results$FDR_pvalue >= 0.05, na.rm = TRUE), 0)
    )
  )
  
  write.csv(summary_data, "../output/tables/cancer_ko_summary.csv", row.names = FALSE)
  return(summary_data)
}

# Generate summary
if(!is.null(ko_analysis_results)) {
  ko_summary <- create_ko_summary(ko_analysis_results, ko_stats)
  print(ko_summary)
}

cat("✓ Cancer knockout analysis complete!\n")
cat("✓ Results saved to output/tables/\n")
cat("✓ Plots saved to output/plots/\n") 
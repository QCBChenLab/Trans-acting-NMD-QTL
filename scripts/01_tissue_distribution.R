# ==============================================================================
# Tissue Distribution Analysis
# This script analyzes and visualizes the distribution of significant SNPs across tissues
# ==============================================================================

# Load data and functions
source("00_data_loading.R")

cat("Starting tissue distribution analysis...\n")

# ==============================================================================
# Significant SNPs Tissue-wise Distribution
# ==============================================================================

# Create tissue distribution plot
create_tissue_distribution_plot <- function(sig_SNP_data) {
  # Count SNPs per tissue
  tissue_counts <- table(sig_SNP_data$tissue)
  
  # Create plot
  tissue_dist_p <- ggplot(sig_SNP_data, aes(x = tissue, fill = tissue)) +
    geom_bar() +
    coord_cartesian(ylim = c(0, 100)) + 
    geom_text(stat = "count", 
              aes(label = after_stat(count), 
                  color = ifelse(after_stat(count) > 20, "highlight", "normal"), 
                  fontface = ifelse(after_stat(count) > 20, "bold", "plain")),
              vjust = -0.5, size = 3) +
    labs(title = "Distribution of Significant SNPs Across Tissues",
         x = "Tissue",
         y = "Number of Significant SNPs") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14),  
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none") +
    scale_fill_manual(values = rep("grey", length(unique(sig_SNP_data$tissue)))) +
    scale_color_manual(values = c("highlight" = "#D73027", "normal" = "darkgrey"))
  
  return(tissue_dist_p)
}

# Generate tissue distribution plot
tissue_plot <- create_tissue_distribution_plot(sig_SNP)

# Save plot
save_plot(tissue_plot, "tissue_distribution", width = 8, height = 4)

# ==============================================================================
# Ordered SNPs Analysis
# ==============================================================================

# Order SNPs by chromosome and position
sig_SNP_ordered <- sig_SNP %>%
  separate(pos, c("chr", "pos"), sep = ":") %>%
  mutate(chr = as.numeric(gsub("chr", "", gsub("chrX", "23", chr))), 
         pos = as.numeric(pos)) %>%
  arrange(chr, pos) %>% 
  select(chr, pos, p_val, tissue)

# Save ordered SNPs
write.csv(sig_SNP_ordered, "../output/tables/sig_SNP_ordered.csv", row.names = FALSE)

cat("✓ Tissue distribution analysis complete!\n")
cat("✓ Plot saved to output/plots/\n")
cat("✓ Summary tables saved to output/tables/\n") 
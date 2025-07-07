# ==============================================================================
# Manhattan Plot Generation
# This script creates Manhattan plots with pathway annotations
# ==============================================================================

# Load data and functions
source("00_data_loading.R")

cat("Starting Manhattan plot generation...\n")

# ==============================================================================
# Load and Prepare Data
# ==============================================================================

# Load plot_man_anno.csv which contains all pathway information
load_plot_annotation_data <- function() {
  tryCatch({
    if(file.exists("../../plot_man_anno.csv")) {
      plot_man_anno <- read.csv("../../plot_man_anno.csv", header = TRUE)[, -1]
      cat("✓ Plot annotation data loaded (", nrow(plot_man_anno), " entries)\n")
      return(plot_man_anno)
    } else {
      cat("⚠ plot_man_anno.csv not found\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Error loading plot annotation data:", e$message, "\n")
    return(NULL)
  })
}

# Load sig_snp_gene_id_mapped_res.csv for closest gene filtering
load_closest_gene_mapping <- function() {
  tryCatch({
    if(file.exists("../../sig_snp_gene_id_mapped_res.csv")) {
      closest_genes <- read.csv("../../sig_snp_gene_id_mapped_res.csv", header = TRUE)
      cat("✓ Closest gene mapping loaded (", nrow(closest_genes), " entries)\n")
      return(closest_genes)
    } else {
      cat("⚠ sig_snp_gene_id_mapped_res.csv not found\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Error loading closest gene mapping:", e$message, "\n")
    return(NULL)
  })
}

plot_man_anno <- load_plot_annotation_data()
closest_genes <- load_closest_gene_mapping()

if(is.null(plot_man_anno)) {
  stop("Cannot proceed without plot annotation data")
}

if(is.null(closest_genes)) {
  stop("Cannot proceed without closest gene mapping data")
}

# ==============================================================================
# Filter Data Using Closest Gene Mapping
# ==============================================================================

# Filter plot_man_anno to use only the closest gene for each SNP
filter_by_closest_gene <- function(plot_man_anno, closest_genes) {
  
  # Create a mapping of Variation.ID to closest gene
  closest_gene_map <- closest_genes %>%
    select(Variation.ID, Gene) %>%
    filter(!is.na(Gene))
  
  cat("✓ Filtering annotations to use only closest genes...\n")
  
  # Join and filter to keep only rows matching the closest gene
  filtered_data <- plot_man_anno %>%
    left_join(closest_gene_map, by = "Variation.ID", suffix = c("", "_closest")) %>%
    filter(Gene == Gene_closest | is.na(Gene_closest)) %>%
    select(-Gene_closest)
  
  # Check for SNPs that had multiple gene options but now have only one
  multi_gene_snps <- plot_man_anno %>%
    group_by(Variation.ID) %>%
    summarise(n_genes = n_distinct(Gene)) %>%
    filter(n_genes > 1)
  
  filtered_multi_gene_snps <- filtered_data %>%
    filter(Variation.ID %in% multi_gene_snps$Variation.ID) %>%
    group_by(Variation.ID) %>%
    summarise(n_genes_after = n_distinct(Gene))
  
  cat("✓ SNPs with multiple gene options before filtering:", nrow(multi_gene_snps), "\n")
  cat("✓ SNPs correctly reduced to single gene:", sum(filtered_multi_gene_snps$n_genes_after == 1), "\n")
  cat("✓ Total annotations after filtering:", nrow(filtered_data), "\n")
  
  return(filtered_data)
}

plot_man_anno_filtered <- filter_by_closest_gene(plot_man_anno, closest_genes)

# ==============================================================================
# Process Data
# ==============================================================================

# Group by variant_id and select the shortest path for each SNP
process_manhattan_data <- function(plot_man_anno) {
  
  # Calculate chromosome midpoints for plotting
  chr_midpoints <- plot_man_anno %>% 
    group_by(CHR) %>%
    summarise(chr_start = min(BPcum), chr_end = max(BPcum)) %>%
    ungroup() %>%
    mutate(chr_start = lag(chr_end, default = 0)) %>%
    mutate(mid = (chr_start + chr_end) / 2)
  
  # Process annotation data - now with filtered data
  plot_man_anno_wide <- plot_man_anno %>%
    group_by(variant_id) %>%
    mutate(nmd_gene_name = ifelse(all(is.na(path_len)), 
                                  NA, 
                                  nmd_gene_name[which.min(path_len)])) %>%
    summarise(Gene_Type = paste(unique(Gene.Type), collapse = ", "),
              shortest_path_len = ifelse(all(is.na(path_len)), NA, min(path_len, na.rm = TRUE)), 
              nmd_gene_name = unique(nmd_gene_name), 
              Gene_name = paste(unique(Gene), collapse = ", ")) %>%
    ungroup() %>% 
    left_join(plot_man_anno %>% 
                select(variant_id, P, BPcum) %>%
                distinct(), by = "variant_id") %>% 
    distinct()
  
  cat("✓ Processed", nrow(plot_man_anno_wide), "unique SNPs for plotting\n")
  cat("✓ SNPs with pathway annotations:", sum(!is.na(plot_man_anno_wide$shortest_path_len)), "\n")
  cat("✓ SNPs with gene annotations:", sum(!is.na(plot_man_anno_wide$Gene_name) & plot_man_anno_wide$Gene_name != ""), "\n")
  
  return(list(
    plot_data = plot_man_anno_wide,
    chr_midpoints = chr_midpoints
  ))
}

# Process the filtered data
processed_data <- process_manhattan_data(plot_man_anno_filtered)
plot_man_anno_wide <- processed_data$plot_data
chr_midpoints <- processed_data$chr_midpoints

# ==============================================================================
# Create Manhattan Plot
# ==============================================================================

create_manhattan_plot <- function(plot_data_wide, chr_midpoints) {
  
  manhattan_p <- ggplot(plot_data_wide, 
                       aes(x = BPcum, y = -log10(P), 
                           color = shortest_path_len, 
                           shape = Gene_Type)) +
    geom_point(size = 3, alpha = 0.6) + 
    labs(x = "Chromosome", y = "-log10(p-value)", title = "Manhattan Plot") +
    theme_bw() +
    geom_rect(data = chr_midpoints, 
              aes(xmin = chr_start, xmax = chr_end, ymin = -Inf, ymax = Inf, fill = factor(CHR)), 
              alpha = 0.2, inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = rep(c("gray90", "white"), 
                                   length.out = length(unique(chr_midpoints$CHR)))) +
    scale_x_continuous(breaks = unique(chr_midpoints$mid), labels = unique(chr_midpoints$CHR)) + 
    scale_color_gradient(name = "Path Length", 
                         low = "red", high = "orange", 
                         na.value = "gray0") + 
    scale_shape_manual(
      name = "Gene Mapping Type",
      values = c("overlapped" = 16, "upstream" = 17, 
                 "downstream" = 15, "upstream, downstream" = 18)) +
    geom_text_repel(data = plot_data_wide, 
                    aes(label = ifelse(is.na(nmd_gene_name), 
                                      Gene_name, 
                                      paste0(Gene_name, " → ", nmd_gene_name))), 
                    fontface = ifelse(is.na(plot_data_wide$nmd_gene_name), 
                                      "plain", "bold"),
                    max.overlaps = 10, 
                    size = 2) +
    theme(axis.text.x = element_text(hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right", 
          legend.title = element_text(size = 10),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  return(manhattan_p)
}

# Generate Manhattan plot
manhattan_plot <- create_manhattan_plot(plot_man_anno_wide, chr_midpoints)

# Save the plot
save_plot(manhattan_plot, "manhattan_plot", width = 9, height = 6)

cat("✓ Manhattan plot saved to output/plots/\n")

# ==============================================================================
# Export Data
# ==============================================================================

# Save only the essential processed Manhattan plot data
write.csv(plot_man_anno_wide, "../output/tables/manhattan_plot_annotated.csv", row.names = FALSE)

# Generate and display summary statistics
annotation_summary <- data.frame(
  metric = c("Total SNPs in plot", 
             "SNPs with gene annotations", 
             "SNPs with pathway annotations",
             "Unique genes annotated",
             "Overlapped genes",
             "Upstream genes", 
             "Downstream genes",
             "SNPs with shortest path length 3",
             "SNPs with shortest path length 4", 
             "SNPs with shortest path length 5",
             "SNPs with shortest path length 6+"),
  count = c(nrow(plot_man_anno_wide),
            sum(!is.na(plot_man_anno_wide$Gene_name) & plot_man_anno_wide$Gene_name != ""),
            sum(!is.na(plot_man_anno_wide$shortest_path_len)),
            length(unique(plot_man_anno_wide$Gene_name[!is.na(plot_man_anno_wide$Gene_name) & plot_man_anno_wide$Gene_name != ""])),
            sum(grepl("overlapped", plot_man_anno_wide$Gene_Type), na.rm = TRUE),
            sum(grepl("upstream", plot_man_anno_wide$Gene_Type), na.rm = TRUE),
            sum(grepl("downstream", plot_man_anno_wide$Gene_Type), na.rm = TRUE),
            sum(plot_man_anno_wide$shortest_path_len == 3, na.rm = TRUE),
            sum(plot_man_anno_wide$shortest_path_len == 4, na.rm = TRUE),
            sum(plot_man_anno_wide$shortest_path_len == 5, na.rm = TRUE),
            sum(plot_man_anno_wide$shortest_path_len >= 6, na.rm = TRUE))
)

# Print summary to console
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Manhattan Plot Annotation Summary:\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
print(annotation_summary)
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("\n✓ Manhattan plot generation complete!\n")
cat("✓ Plot saved to output/plots/manhattan_plot.*\n")
cat("✓ Data saved to output/tables/manhattan_plot_annotated.csv\n")
cat("✓ Used closest gene mapping to filter annotations\n")
cat("✓ Each SNP now has only its closest gene annotation\n") 
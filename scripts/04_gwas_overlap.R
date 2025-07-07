# ==============================================================================
# GWAS Catalog Overlap Analysis
# This script analyzes overlap between significant SNPs and GWAS catalog data
# ==============================================================================

# Load data and functions
source("00_data_loading.R")

cat("Starting GWAS catalog overlap analysis...\n")

# ==============================================================================
# Process GWAS Catalog Data
# ==============================================================================

# Check if GWAS catalog data is loaded
if(is.null(sig_snp_gwas_catalog)) {
  cat("⚠ No GWAS catalog data available. Skipping analysis.\n")
  stop("sig_snp_gwas_catalog.csv must be available in GWAS_overlap/GWAS_Catalog/")
}

cat("✓ Using pre-loaded GWAS catalog overlap data with", nrow(sig_snp_gwas_catalog), "entries\n")

# Define disease categories following Manuscript.Rmd
diseases_category <- c("Cancer", "Cardiovascular disease", 
                      "Neurological disorder", "Metabolic disorder", 
                      "Immune system disorder", "Digestive system disorder", 
                      "Other disease")

# Filter for disease categories only
sig_snp_gwas_catalog_diseases <- sig_snp_gwas_catalog %>% 
  filter(Category %in% diseases_category)

cat("✓ Filtered to", nrow(sig_snp_gwas_catalog_diseases), "disease-related entries\n")

# ==============================================================================
# Tissue-wise GWAS Disease Distribution
# ==============================================================================

# Create tissue-wise disease count following Manuscript.Rmd approach
# First ensure tissue is a character column
sig_snp_gwas_catalog_diseases <- sig_snp_gwas_catalog_diseases %>%
  mutate(tissue = as.character(tissue))

tissue_gwas_data <- sig_snp_gwas_catalog_diseases %>% 
  select(tissue, DISEASE.TRAIT) %>%
  distinct() %>%
  group_by(tissue) %>%
  summarise(n = n()) %>%
  ungroup()

# Add any missing tissues (like Spleen if not present)
if(!"Spleen" %in% tissue_gwas_data$tissue) {
  tissue_gwas_data <- bind_rows(tissue_gwas_data, 
                                data.frame(tissue = "Spleen", n = 0))
}

# Create bar plot for tissue distribution
tissue_gwas_p <- ggplot(tissue_gwas_data, aes(x = tissue, 
                                              y = n, 
                                              fill = tissue)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 150)) + 
  geom_text(aes(label = n, 
                color = ifelse(n > 50, "highlight", "normal"), 
                fontface = ifelse(n > 50, "bold", "plain")),
            vjust = -0.5, size = 3) +
  labs(title = "Distribution of Matched GWAS Catalog Diseases Across Tissues",
       x = "Tissue",
       y = "Number of Matched Traits") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14),  
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = rep("grey", length(unique(tissue_gwas_data$tissue)))) +
  scale_color_manual(values = c("highlight" = "#D73027", "normal" = "darkgrey"))

# Save tissue distribution plot
save_plot(tissue_gwas_p, "tissue_gwas_diseases", width = 7.2, height = 3.6)

# ==============================================================================
# Disease Category Distribution by Tissue
# ==============================================================================

# Prepare data for pie charts
tissue_gwas_cat_data <- sig_snp_gwas_catalog_diseases %>% 
  mutate(tissue = as.character(tissue),
         Category = as.character(Category)) %>%
  select(tissue, DISEASE.TRAIT, Category) %>%
  distinct() %>%
  group_by(tissue, Category) %>%
  summarise(n = n(), .groups = "drop")

# Add missing combinations including Spleen
all_tissues_cat <- unique(c(tissue_gwas_cat_data$tissue, "Spleen"))
all_categories <- unique(tissue_gwas_cat_data$Category)

# Create complete grid
complete_grid <- expand.grid(tissue = all_tissues_cat, 
                            Category = all_categories,
                            stringsAsFactors = FALSE)

# Merge with actual data
tissue_gwas_cat_data <- complete_grid %>%
  left_join(tissue_gwas_cat_data, by = c("tissue", "Category")) %>%
  mutate(n = ifelse(is.na(n), 0, n))

# Create pie chart faceted by tissue
tissue_pie_chart_p <- ggplot(tissue_gwas_cat_data, 
                            aes(x = "", y = n, fill = Category)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar(theta = "y") + 
  facet_wrap(~tissue, scales = "free", nrow = 1, 
             strip.position = "bottom") + 
  labs(title = "Distribution of Disease Categories Across Tissues", 
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom", 
        strip.text = element_text(size = 6, angle = 90), 
        panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA),   
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)) +
  scale_fill_brewer(palette = "Set3")

# Save pie chart
save_plot(tissue_pie_chart_p, "tissue_gwas_cat_pie_diseases", width = 8, height = 4)

# ==============================================================================
# Export Summary Tables
# ==============================================================================

# Save tissue distribution data
write.csv(tissue_gwas_data, "../output/tables/tissue_gwas_distribution.csv", row.names = FALSE)

# Save category distribution data
write.csv(tissue_gwas_cat_data, "../output/tables/tissue_gwas_category_distribution.csv", row.names = FALSE)

# ==============================================================================
# Summary Statistics
# ==============================================================================

# Create comprehensive summary
gwas_summary <- data.frame(
  metric = c(
    "Total GWAS catalog entries",
    "Disease-related entries",
    "Unique tissues with GWAS overlaps",
    "Unique diseases/traits",
    "Unique disease categories",
    "Most common disease category",
    "Tissue with most GWAS overlaps"
  ),
  value = c(
    nrow(sig_snp_gwas_catalog),
    nrow(sig_snp_gwas_catalog_diseases),
    length(unique(sig_snp_gwas_catalog_diseases$tissue)),
    length(unique(sig_snp_gwas_catalog_diseases$DISEASE.TRAIT)),
    length(unique(sig_snp_gwas_catalog_diseases$Category)),
    names(sort(table(sig_snp_gwas_catalog_diseases$Category), decreasing = TRUE)[1]),
    tissue_gwas_data$tissue[which.max(tissue_gwas_data$n)]
  )
)

# Save summary
write.csv(gwas_summary, "../output/tables/gwas_overlap_summary.csv", row.names = FALSE)

# Print summary
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("GWAS Catalog Overlap Summary:\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
print(gwas_summary)
cat(paste(rep("=", 50), collapse = ""), "\n")

# ==============================================================================
# Disease-specific Analysis
# ==============================================================================

# Identify top diseases by frequency
top_diseases <- sig_snp_gwas_catalog_diseases %>%
  mutate(DISEASE.TRAIT = as.character(DISEASE.TRAIT)) %>%
  group_by(DISEASE.TRAIT) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  head(20)

# Create plot for top diseases
top_diseases_p <- ggplot(top_diseases, aes(x = reorder(DISEASE.TRAIT, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 GWAS Catalog Diseases/Traits",
       x = "Disease/Trait",
       y = "Number of SNP Overlaps") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 14))

# Save top diseases plot
save_plot(top_diseases_p, "top_gwas_diseases", width = 8, height = 6)

# ==============================================================================
# Tissue-specific Disease Profiles
# ==============================================================================

# Get top tissues by GWAS overlap count
top_tissues <- tissue_gwas_data %>%
  filter(n > 0) %>%
  arrange(desc(n)) %>%
  head(5)

# For each top tissue, show its disease profile
for(i in 1:nrow(top_tissues)) {
  tissue_name <- top_tissues$tissue[i]
  
  tissue_diseases <- sig_snp_gwas_catalog_diseases %>%
    filter(tissue == tissue_name) %>%
    mutate(Category = as.character(Category)) %>%
    group_by(Category) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(desc(n))
  
  tissue_disease_p <- ggplot(tissue_diseases, aes(x = reorder(Category, n), y = n)) +
    geom_bar(stat = "identity", fill = "coral") +
    coord_flip() +
    labs(title = paste("Disease Categories in", tissue_name),
         x = "Disease Category",
         y = "Number of Overlaps") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 12))
  
  # Save individual tissue plots
  save_plot(tissue_disease_p, 
            paste0("gwas_profile_", gsub(" ", "_", tissue_name)), 
            width = 6, height = 4)
}

cat("\n✓ GWAS catalog overlap analysis complete!\n")
cat("✓ Plots saved to output/plots/\n")
cat("✓ Summary tables saved to output/tables/\n")
cat("✓ Analysis based on pre-processed GWAS catalog overlap data\n") 
#!/usr/bin/env Rscript
# ==============================================================================
# Outlier Deviation Analysis for NMD Factors
# ==============================================================================
# Analyzes relationship between NMD factor expression and outlier deviations
# Optimized for large datasets - processes files individually
# ==============================================================================

library(tidyverse)
library(data.table)
library(patchwork)
library(RColorBrewer)

# ==============================================================================
# Configuration
# ==============================================================================
# BASE_DIR: path to the main project directory containing Cancer_KO/, Manuscript/ etc.
BASE_DIR <- Sys.getenv("NMD_PROJECT_DIR", unset = normalizePath("../../..", mustWork = FALSE))

data_dir <- file.path(BASE_DIR, "Cancer_KO/Boxplots/NMD_Deviation_by_ExprGroup/")
output_dir <- "../output/plots/deviation_analysis/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Auto-detect all available genes from data files
data_files <- list.files(data_dir, pattern = "^plot_data_.*\\.csv$", full.names = FALSE)
genes <- gsub("^plot_data_|\\.csv$", "", data_files)
cat(sprintf("Found %d genes with data files\n", length(genes)))

thresholds <- c(1.5, 2.0, 2.5)

# ==============================================================================
# Process each gene file individually (memory efficient)
# ==============================================================================
compute_gene_stats <- function(gene) {
  file_path <- file.path(data_dir, paste0("plot_data_", gene, ".csv"))
  if (!file.exists(file_path)) return(NULL)
  
  cat(sprintf("Processing %s...\n", gene))
  dt <- fread(file_path, select = c("SAMPID", "Gene_log", "Tissue", "Deviation"))
  dt[, Gene_Symbol := gene]
  
  # Decile assignment
  dt[, Decile := ntile(Gene_log, 10)]
  
  # Outlier proportions by decile and threshold
  outlier_props <- rbindlist(lapply(thresholds, function(thr) {
    dt[, .(
      N_total = .N,
      N_outliers = sum(abs(Deviation) >= thr, na.rm = TRUE),
      Prop_outliers = sum(abs(Deviation) >= thr, na.rm = TRUE) / .N,
      Threshold = thr
    ), by = .(Gene_Symbol, Decile)]
  }))
  
  # Per-sample statistics
  sample_stats <- dt[, .(
    Gene_log = first(Gene_log),
    Tissue = first(Tissue),
    Mean_absD = mean(abs(Deviation), na.rm = TRUE),
    N_transcripts = .N
  ), by = .(Gene_Symbol, SAMPID)]
  
  # Correlation
  cor_result <- sample_stats[, {
    ct <- cor.test(Gene_log, Mean_absD)
    .(Correlation = ct$estimate, P_value = ct$p.value, 
      N_samples = .N, Tissue = first(Tissue))
  }, by = Gene_Symbol]
  
  list(outlier_props = outlier_props, 
       sample_stats = sample_stats, 
       correlation = cor_result)
}

# Process all genes
results <- lapply(genes, compute_gene_stats)
results <- results[!sapply(results, is.null)]

# Combine results
outlier_props <- rbindlist(lapply(results, `[[`, "outlier_props"))
sample_means <- rbindlist(lapply(results, `[[`, "sample_stats"))
correlations <- rbindlist(lapply(results, `[[`, "correlation")) %>%
  arrange(desc(Correlation))

cat(sprintf("\nProcessed %d genes successfully\n", length(results)))

# ==============================================================================
# Plot 1: Line plots - Outlier proportions by expression decile
# ==============================================================================
plot_outlier_trends <- function(data, threshold) {
  df <- as.data.frame(data[Threshold == threshold])
  
  ggplot(df, aes(x = Decile, y = Prop_outliers, color = Gene_Symbol)) +
    geom_line(linewidth = 0.9, alpha = 0.85) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_x_continuous(breaks = 1:10, labels = paste0("D", 1:10)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(genes))) +
    labs(
      title = sprintf("Proportion of Outlier Deviations (|Dev| >= %.1f)", threshold),
      subtitle = "By NMD factor expression decile",
      x = "Expression Decile (Low -> High)",
      y = "Proportion of Outliers",
      color = "NMD Factor"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold", size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.3, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.spacing.x = unit(0.05, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
      legend.box.margin = margin(0, 0, 0, 0)
    ) +
    guides(color = guide_legend(nrow = 7, byrow = TRUE, 
                                override.aes = list(size = 1.5, linewidth = 0.8)))
}

# Save individual threshold plots
for (thr in thresholds) {
  p <- plot_outlier_trends(outlier_props, thr)
  ggsave(file.path(output_dir, sprintf("outlier_prop_threshold_%.1f.pdf", thr)),
         plot = p, width = 14, height = 8, dpi = 300)
  ggsave(file.path(output_dir, sprintf("outlier_prop_threshold_%.1f.jpeg", thr)),
         plot = p, width = 14, height = 8, dpi = 450)
}
cat("Saved outlier proportion line plots\n")

# Combined faceted plot
p_combined <- as.data.frame(outlier_props) %>%
  mutate(Threshold_label = sprintf("|Dev| >= %.1f", Threshold)) %>%
  ggplot(aes(x = Decile, y = Prop_outliers, color = Gene_Symbol)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.8, alpha = 0.85) +
  facet_wrap(~ Threshold_label, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(genes))) +
  labs(
    title = "Proportion of Outlier NMD Transcript Deviations",
    subtitle = "Higher NMD factor expression is associated with more extreme deviations",
    x = "Expression Decile (Low -> High)",
    y = "Proportion of Outliers",
    color = "NMD Factor"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.35, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.spacing.y = unit(0.03, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(0, 0, 0, 0),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE,
                              override.aes = list(size = 1.5, linewidth = 0.8)))

ggsave(file.path(output_dir, "outlier_prop_combined.pdf"),
       plot = p_combined, width = 15, height = 7, dpi = 300)
ggsave(file.path(output_dir, "outlier_prop_combined.jpeg"),
       plot = p_combined, width = 15, height = 7, dpi = 450)
cat("Saved combined outlier proportion plot\n")

# ==============================================================================
# Plot 2: Scatter plots - Gene expression vs Mean |Deviation|
# ==============================================================================
scatter_dir <- file.path(output_dir, "scatter_individual")
dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)

# Function for grid plots (smaller)
plot_scatter_small <- function(gene, sample_means, correlations) {
  df <- as.data.frame(sample_means[Gene_Symbol == gene])
  cor_info <- as.data.frame(correlations[Gene_Symbol == gene])
  
  ggplot(df, aes(x = Gene_log, y = Mean_absD)) +
    geom_point(alpha = 0.4, color = "steelblue", size = 1.8) +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE, linewidth = 1) +
    labs(
      title = gene,
      subtitle = sprintf("r = %.3f, p = %.2e", 
                         cor_info$Correlation, cor_info$P_value),
      x = "log(TPM + 1)",
      y = "Mean |Deviation|"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40")
    )
}

# Function for individual plots (larger, more detailed)
plot_scatter_individual <- function(gene, sample_means, correlations, output_dir) {
  df <- as.data.frame(sample_means[Gene_Symbol == gene])
  cor_info <- as.data.frame(correlations[Gene_Symbol == gene])
  
  p <- ggplot(df, aes(x = Gene_log, y = Mean_absD)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 2.5) +
    geom_smooth(method = "lm", color = "firebrick", se = TRUE, linewidth = 1.2) +
    labs(
      title = sprintf("%s Expression vs Mean Absolute Deviation", gene),
      subtitle = sprintf("Pearson r = %.3f, p = %.2e, n = %d samples", 
                         cor_info$Correlation, cor_info$P_value, cor_info$N_samples),
      x = sprintf("%s Expression [log(TPM + 1)]", gene),
      y = "Mean |Deviation| of NMD Transcripts"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
  
  ggsave(file.path(output_dir, sprintf("scatter_%s.pdf", gene)),
         plot = p, width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, sprintf("scatter_%s.jpeg", gene)),
         plot = p, width = 8, height = 6, dpi = 450)
  
  return(p)
}

# Generate individual scatter plots for ALL genes
cat("Saving individual scatter plots...\n")
for (gene in genes) {
  plot_scatter_individual(gene, sample_means, correlations, scatter_dir)
}
cat(sprintf("Saved %d individual scatter plots to %s\n", length(genes), scatter_dir))

# Generate grid plot for SIGNIFICANT genes only (p < 0.05)
sig_correlations <- as.data.frame(correlations) %>%
  filter(!is.na(Correlation), P_value < 0.05) %>%
  arrange(desc(Correlation))
sig_genes_for_grid <- sig_correlations$Gene_Symbol

cat(sprintf("Creating grid plot for %d significant genes (p < 0.05)\n", length(sig_genes_for_grid)))

scatter_plots <- lapply(sig_genes_for_grid, plot_scatter_small, sample_means, correlations)

# Calculate grid dimensions
n_sig <- length(sig_genes_for_grid)
ncol_grid <- 5
nrow_grid <- ceiling(n_sig / ncol_grid)
grid_height <- max(8, nrow_grid * 2.5)

p_grid <- wrap_plots(scatter_plots, ncol = ncol_grid) +
  plot_annotation(
    title = "NMD Factor Expression vs Mean Absolute Deviation",
    subtitle = sprintf("Significant correlations only (p < 0.05, n = %d genes)", n_sig),
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 12)
    )
  )

ggsave(file.path(output_dir, "scatter_mean_absD_grid.pdf"),
       plot = p_grid, width = 16, height = grid_height, dpi = 300)
ggsave(file.path(output_dir, "scatter_mean_absD_grid.jpeg"),
       plot = p_grid, width = 16, height = grid_height, dpi = 450)
cat("Saved scatter plot grid (significant genes only)\n")

# ==============================================================================
# Plot 3: Correlation bar plots
# ==============================================================================
cor_plot_data <- as.data.frame(correlations) %>%
  filter(!is.na(Correlation)) %>%
  mutate(
    Gene_Symbol = fct_reorder(Gene_Symbol, Correlation),
    # More granular significance labels
    Sig_label = case_when(
      P_value < 1e-50 ~ "****",
      P_value < 1e-20 ~ "***",
      P_value < 1e-5 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Color-coded significance category
    Sig_category = case_when(
      P_value < 1e-50 ~ "p < 1e-50",
      P_value < 1e-20 ~ "p < 1e-20",
      P_value < 1e-5 ~ "p < 1e-5",
      P_value < 0.05 ~ "p < 0.05",
      TRUE ~ "n.s."
    ),
    Sig_category = factor(Sig_category, levels = c("p < 1e-50", "p < 1e-20", "p < 1e-5", "p < 0.05", "n.s.")),
    Is_significant = P_value < 0.05
  )

# Plot 3a: Horizontal bar plot (ALL genes)
p_cor_bar_all <- cor_plot_data %>%
  ggplot(aes(x = Gene_Symbol, y = Correlation, fill = Correlation)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_text(aes(label = Sig_label), hjust = -0.2, size = 2.5, fontface = "bold") +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-0.15, 0.7)) +
  labs(
    title = "Correlation: NMD Factor Expression vs Mean |Deviation|",
    subtitle = "All genes (* p<0.05, ** p<1e-5, *** p<1e-20, **** p<1e-50)",
    x = NULL,
    y = "Pearson Correlation"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10),
    axis.text.y = element_text(size = 7),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  coord_flip(ylim = c(-0.15, 0.75))

ggsave(file.path(output_dir, "correlation_barplot_all.pdf"),
       plot = p_cor_bar_all, width = 8, height = 14, dpi = 300)
ggsave(file.path(output_dir, "correlation_barplot_all.jpeg"),
       plot = p_cor_bar_all, width = 8, height = 14, dpi = 450)
cat("Saved correlation bar plot (all genes, horizontal)\n")

# Plot 3b: Vertical bar plot (SIGNIFICANT genes only, on x-axis)
sig_cor_data <- cor_plot_data %>%
  filter(P_value < 0.05) %>%
  mutate(
    Gene_Symbol = fct_reorder(Gene_Symbol, Correlation),
    # Dynamic positioning: above for positive, below for negative correlations
    label_vjust = ifelse(Correlation >= 0, -0.3, 1.3)
  )

p_cor_bar_sig <- sig_cor_data %>%
  ggplot(aes(x = Gene_Symbol, y = Correlation, fill = Sig_category)) +
  geom_col(width = 0.75, alpha = 0.9) +
  geom_text(aes(label = Sig_label, vjust = label_vjust), size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c("p < 1e-50" = "#67000d", "p < 1e-20" = "#a50f15", 
               "p < 1e-5" = "#ef3b2c", "p < 0.05" = "#fc9272"),
    name = "Significance"
  ) +
  labs(
    title = "Correlation: NMD Factor Expression vs Mean |Deviation|",
    subtitle = sprintf("Significant genes only (p < 0.05, n = %d)", nrow(sig_cor_data)),
    x = "Gene",
    y = "Pearson Correlation"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8, face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  ) +
  coord_cartesian(ylim = c(min(sig_cor_data$Correlation) * 1.3, max(sig_cor_data$Correlation) * 1.15))

ggsave(file.path(output_dir, "correlation_barplot.pdf"),
       plot = p_cor_bar_sig, width = 14, height = 7, dpi = 300)
ggsave(file.path(output_dir, "correlation_barplot.jpeg"),
       plot = p_cor_bar_sig, width = 14, height = 7, dpi = 450)
cat(sprintf("Saved correlation bar plot (significant genes, n=%d)\n", nrow(sig_cor_data)))

# ==============================================================================
# Plot 4: Heatmap of outlier proportions (significant genes only)
# ==============================================================================
sig_genes <- cor_plot_data %>% filter(P_value < 0.05) %>% pull(Gene_Symbol) %>% as.character()

heatmap_df <- as.data.frame(outlier_props[Threshold == 2.0]) %>%
  filter(Gene_Symbol %in% sig_genes) %>%
  mutate(
    Decile_label = factor(paste0("D", Decile), levels = paste0("D", 1:10)),
    Gene_Symbol = factor(Gene_Symbol, levels = sig_genes[order(
      cor_plot_data$Correlation[match(sig_genes, cor_plot_data$Gene_Symbol)], decreasing = TRUE
    )])
  )

n_sig_genes <- length(sig_genes)
heatmap_height <- max(8, n_sig_genes * 0.3)

p_heatmap <- ggplot(heatmap_df, aes(x = Decile_label, y = Gene_Symbol, fill = Prop_outliers)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", Prop_outliers * 100)), 
            size = 2.2, color = "black") +
  scale_fill_gradient(low = "white", high = "firebrick",
                      labels = scales::percent_format(accuracy = 1),
                      name = "Outlier %") +
  labs(
    title = sprintf("Outlier Proportion (|Dev| >= 2.0) by Expression Decile (n=%d significant genes)", n_sig_genes),
    x = "Expression Decile (Low -> High)",
    y = "NMD Factor"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "outlier_heatmap.pdf"),
       plot = p_heatmap, width = 12, height = heatmap_height, dpi = 300)
ggsave(file.path(output_dir, "outlier_heatmap.jpeg"),
       plot = p_heatmap, width = 12, height = heatmap_height, dpi = 450)
cat(sprintf("Saved heatmap with %d significant genes\n", n_sig_genes))

# ==============================================================================
# Save summary tables
# ==============================================================================
write_csv(as.data.frame(correlations), file.path(output_dir, "correlation_summary.csv"))

decile_fold <- as.data.frame(outlier_props[Threshold == 2.0]) %>%
  group_by(Gene_Symbol) %>%
  summarise(
    D1 = Prop_outliers[Decile == 1],
    D10 = Prop_outliers[Decile == 10],
    Fold_increase = D10 / max(D1, 0.001)
  ) %>%
  arrange(desc(Fold_increase))
write_csv(decile_fold, file.path(output_dir, "decile_fold_increase.csv"))

cat("\n=== Correlation Summary ===\n")
print(as.data.frame(correlations))

cat("\n=== Fold Increase (D10/D1) at threshold 2.0 ===\n")
print(decile_fold)

cat(sprintf("\n=== All outputs saved to: %s ===\n", output_dir))

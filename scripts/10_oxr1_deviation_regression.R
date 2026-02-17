#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5B: OXR1 Expression vs Average Absolute Deviation
# Brain Cerebellum Samples
# ==============================================================================

library(tidyverse)
library(ggpubr)

# ==============================================================================
# Configuration
# ==============================================================================
# BASE_DIR: path to the main project directory containing Cancer_KO/ etc.
BASE_DIR <- Sys.getenv("NMD_PROJECT_DIR", unset = normalizePath("../../..", mustWork = FALSE))

data_path <- file.path(BASE_DIR, "Cancer_KO/Boxplots/NMD_Deviation_by_ExprGroup/plot_data_OXR1.csv")
output_dir <- "../output/plots/"

# ==============================================================================
# Load and Process Data
# ==============================================================================
cat("Loading OXR1 plot data...\n")

# Read just the header first to understand structure
con <- file(data_path, "r")
header <- readLines(con, n = 1)
close(con)
cat("Header:", header, "\n")

# Read with data.table for large files
library(data.table)
df <- fread(data_path)

cat(sprintf("Loaded %d rows\n", nrow(df)))
cat(sprintf("Columns: %s\n", paste(names(df), collapse = ", ")))

# Check tissue
cat(sprintf("Tissue: %s\n", unique(df$Tissue)))

# Calculate per-sample mean absolute deviation
sample_stats <- df %>%
  group_by(SAMPID) %>%
  summarise(
    Gene_log = first(Gene_log),
    Tissue = first(Tissue),
    Mean_absD = mean(abs(Deviation), na.rm = TRUE),
    N_transcripts = n(),
    .groups = "drop"
  )

cat(sprintf("Unique samples: %d\n", nrow(sample_stats)))

# ==============================================================================
# Linear Regression Analysis
# ==============================================================================
lm_model <- lm(Mean_absD ~ Gene_log, data = sample_stats)
lm_summary <- summary(lm_model)

# Extract regression statistics
slope <- coef(lm_model)["Gene_log"]
intercept <- coef(lm_model)["(Intercept)"]
r_squared <- lm_summary$r.squared
adj_r_squared <- lm_summary$adj.r.squared
f_stat <- lm_summary$fstatistic[1]
p_value <- pf(lm_summary$fstatistic[1], lm_summary$fstatistic[2], lm_summary$fstatistic[3], lower.tail = FALSE)
slope_se <- lm_summary$coefficients["Gene_log", "Std. Error"]
slope_pval <- lm_summary$coefficients["Gene_log", "Pr(>|t|)"]

cat(sprintf("\nLinear Regression Analysis:\n"))
cat(sprintf("  Slope (β) = %.4f\n", slope))
cat(sprintf("  Slope SE = %.4f\n", slope_se))
cat(sprintf("  Slope p-value = %.4e\n", slope_pval))
cat(sprintf("  Intercept = %.4f\n", intercept))
cat(sprintf("  R² = %.4f\n", r_squared))
cat(sprintf("  Adjusted R² = %.4f\n", adj_r_squared))
cat(sprintf("  F-statistic = %.4f\n", f_stat))
cat(sprintf("  Model p-value = %.4e\n", p_value))
cat(sprintf("  n = %d samples\n", nrow(sample_stats)))

# ==============================================================================
# Figure 5B: Scatter plot - OXR1 Expression vs Mean |Deviation|
# ==============================================================================
cat("\nGenerating Figure 5B...\n")

# Format p-value
format_pval <- function(p) {
  if (p < 0.0001) {
    return(sprintf("%.2e", p))
  } else if (p < 0.01) {
    return(sprintf("%.4f", p))
  } else {
    return(sprintf("%.2f", p))
  }
}

# Subtitle with regression results
subtitle_text <- sprintf("Brain - Cerebellum (n = %d samples)\nSlope = %.3f, R-squared = %.3f, p = %s",
                         nrow(sample_stats), slope, r_squared, format_pval(slope_pval))

p_fig5b <- ggplot(sample_stats, aes(x = Gene_log, y = Mean_absD)) +
  geom_point(alpha = 0.6, color = "#4292C6", size = 2.5) +
  geom_smooth(method = "lm", color = "#D62728", se = TRUE, 
              linewidth = 1.2, fill = "#FFCCCC", alpha = 0.3) +
  labs(
    title = "OXR1 Expression vs Mean Absolute Deviation",
    subtitle = subtitle_text,
    x = "OXR1 Expression [log(TPM + 1)]",
    y = "Mean |Deviation| of NMD Transcripts"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white")
  )

ggsave(file.path(output_dir, "Fig5b_OXR1_vs_deviation_cerebellum.pdf"),
       plot = p_fig5b, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "Fig5b_OXR1_vs_deviation_cerebellum.jpeg"),
       plot = p_fig5b, width = 6, height = 5, dpi = 450)

cat("Saved: Fig5b_OXR1_vs_deviation_cerebellum.pdf/jpeg\n")

# ==============================================================================
# Alternative version with annotation inside plot
# ==============================================================================
annotation_text <- sprintf("Slope = %.3f\nR-squared = %.3f\np = %s\nn = %d", 
                           slope, r_squared, format_pval(slope_pval), nrow(sample_stats))

p_fig5b_v2 <- ggplot(sample_stats, aes(x = Gene_log, y = Mean_absD)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE, 
              linewidth = 1, fill = "pink", alpha = 0.3) +
  annotate("text", x = min(sample_stats$Gene_log) + 0.3, 
           y = max(sample_stats$Mean_absD) - 0.02,
           label = annotation_text, hjust = 0, vjust = 1, size = 3.5) +
  labs(
    title = "OXR1 Expression vs NMD Transcript Deviation",
    subtitle = "Brain - Cerebellum samples",
    x = "OXR1 Expression [log(TPM + 1)]",
    y = "Mean |Deviation| of NMD Transcripts"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(output_dir, "Fig5b_OXR1_vs_deviation_v2.pdf"),
       plot = p_fig5b_v2, width = 5.5, height = 4.5, dpi = 300)
ggsave(file.path(output_dir, "Fig5b_OXR1_vs_deviation_v2.jpeg"),
       plot = p_fig5b_v2, width = 5.5, height = 4.5, dpi = 450)

cat("Saved: Fig5b_OXR1_vs_deviation_v2.pdf/jpeg\n")

cat("\n=== Figure 5B Complete ===\n")

#!/usr/bin/env Rscript
# ==============================================================================
# OXR1 Knockdown Analysis - Figure 5a
# ==============================================================================
# Analyzes NMD transcript expression changes upon OXR1 knockdown
# Uses SIGNED log fold change to preserve direction (downregulation)
# ==============================================================================

library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# ==============================================================================
# Configuration
# ==============================================================================
# Uses data file included in repo data/ directory
data_path <- "../data/oxr1_tx_count_data_nmd_anno.csv"
output_dir <- "../output/plots/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Load and process data
# ==============================================================================
cat("Loading OXR1 knockdown data...\n")
df <- read.csv(data_path, row.names = 1)

cat(sprintf("Total transcripts: %d\n", nrow(df)))
cat(sprintf("NMD transcripts: %d\n", sum(df$NMD == "NMD")))
cat(sprintf("Non-NMD transcripts: %d\n", sum(df$NMD == "Non_NMD")))

# ==============================================================================
# Calculate SIGNED log fold change (preserving direction)
# ==============================================================================
# Data is already log-transformed, so logFC = siOXR1 - siCon
# Positive = upregulated upon OXR1 KO
# Negative = downregulated upon OXR1 KO

df <- df %>%
  mutate(
    # Log fold change for No Treatment condition
    logFC_NT = siOXR1_NT - siCon_NT,
    # Log fold change for Recovery 0h condition  
    logFC_R0h = siOXR1_R0h - siCon_R0h,
    # Average log fold change across conditions
    logFC_mean = (logFC_NT + logFC_R0h) / 2,
    # NMD status as factor
    NMD_status = factor(NMD, levels = c("Non_NMD", "NMD"))
  )

# ==============================================================================
# Summary statistics
# ==============================================================================
cat("\n=== Summary Statistics (Signed Log Fold Change) ===\n")

summary_stats <- df %>%
  group_by(NMD_status) %>%
  summarise(
    N = n(),
    Mean_logFC = mean(logFC_mean, na.rm = TRUE),
    Median_logFC = median(logFC_mean, na.rm = TRUE),
    SD_logFC = sd(logFC_mean, na.rm = TRUE),
    Pct_downregulated = mean(logFC_mean < 0, na.rm = TRUE) * 100,
    Pct_upregulated = mean(logFC_mean > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

print(summary_stats)

# Statistical test: Compare NMD vs Non-NMD
wilcox_test <- wilcox.test(
  logFC_mean ~ NMD_status, 
  data = df,
  alternative = "two.sided"
)
cat(sprintf("\nWilcoxon test (NMD vs Non-NMD): p = %.2e\n", wilcox_test$p.value))

# One-sample test: Is NMD logFC significantly different from 0?
nmd_data <- df %>% filter(NMD_status == "NMD") %>% pull(logFC_mean)
nmd_ttest <- t.test(nmd_data, mu = 0)
cat(sprintf("One-sample t-test (NMD logFC != 0): p = %.2e, mean = %.4f\n", 
            nmd_ttest$p.value, nmd_ttest$estimate))

non_nmd_data <- df %>% filter(NMD_status == "Non_NMD") %>% pull(logFC_mean)
non_nmd_ttest <- t.test(non_nmd_data, mu = 0)
cat(sprintf("One-sample t-test (Non-NMD logFC != 0): p = %.2e, mean = %.4f\n", 
            non_nmd_ttest$p.value, non_nmd_ttest$estimate))

# ==============================================================================
# Plot 1: Boxplot - Signed Log Fold Change by NMD Status
# ==============================================================================
cat("\nGenerating Figure 5a - Signed Log Fold Change boxplot...\n")

# Color palette
colors <- c("Non_NMD" = "#3498db", "NMD" = "#e74c3c")

# Create significance label
if (wilcox_test$p.value < 1e-100) {
  sig_label <- "p < 1e-100"
} else if (wilcox_test$p.value < 1e-50) {
  sig_label <- sprintf("p = %.0e", wilcox_test$p.value)
} else {
  sig_label <- sprintf("p = %.2e", wilcox_test$p.value)
}

p1 <- ggplot(df, aes(x = NMD_status, y = logFC_mean, fill = NMD_status)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_manual(values = colors) +
  labs(
    title = "NMD Transcript Expression Changes Upon OXR1 Knockdown",
    subtitle = sprintf("Wilcoxon test: %s", sig_label),
    x = "Transcript Type",
    y = "Log2 Fold Change (siOXR1 / siControl)",
    fill = "Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(labels = c("Non_NMD" = "Non-NMD\nTranscripts", 
                               "NMD" = "NMD\nTranscripts")) +
  annotate("text", x = 1.5, y = max(df$logFC_mean, na.rm = TRUE) * 0.9,
           label = sig_label, size = 4, fontface = "italic")

ggsave(file.path(output_dir, "Fig5a_oxr1_logFC_signed.pdf"),
       plot = p1, width = 6, height = 6, dpi = 300)
ggsave(file.path(output_dir, "Fig5a_oxr1_logFC_signed.jpeg"),
       plot = p1, width = 6, height = 6, dpi = 450)

cat("Saved: Fig5a_oxr1_logFC_signed.pdf/jpeg\n")

# ==============================================================================
# Plot 2: Violin plot with individual points (more detail)
# ==============================================================================
p2 <- ggplot(df, aes(x = NMD_status, y = logFC_mean, fill = NMD_status)) +
  geom_violin(alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_manual(values = colors) +
  labs(
    title = "NMD Transcript Expression Changes Upon OXR1 Knockdown",
    subtitle = sprintf("Wilcoxon test: %s | NMD mean logFC = %.3f", 
                       sig_label, mean(nmd_data, na.rm = TRUE)),
    x = "Transcript Type",
    y = "Log2 Fold Change (siOXR1 / siControl)",
    fill = "Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(labels = c("Non_NMD" = "Non-NMD\nTranscripts", 
                               "NMD" = "NMD\nTranscripts"))

ggsave(file.path(output_dir, "Fig5a_oxr1_logFC_violin.pdf"),
       plot = p2, width = 6, height = 6, dpi = 300)
ggsave(file.path(output_dir, "Fig5a_oxr1_logFC_violin.jpeg"),
       plot = p2, width = 6, height = 6, dpi = 450)

cat("Saved: Fig5a_oxr1_logFC_violin.pdf/jpeg\n")

# ==============================================================================
# Plot 3: DOWNREGULATED TRANSCRIPTS ONLY - Highlighting downregulation
# ==============================================================================
cat("\n=== Analyzing DOWNREGULATED transcripts specifically ===\n")

df_long <- df %>%
  select(Transcript_ID, NMD_status, logFC_NT, logFC_R0h) %>%
  pivot_longer(
    cols = c(logFC_NT, logFC_R0h),
    names_to = "Condition",
    values_to = "logFC"
  ) %>%
  mutate(
    Condition = factor(Condition, 
                       levels = c("logFC_NT", "logFC_R0h"),
                       labels = c("NT", "R0h")),
    NMD_status = factor(NMD_status, 
                        levels = c("NMD", "Non_NMD"),
                        labels = c("NMD-targeted", "non-NMD-targeted")),
    # For downregulated transcripts, use absolute value of logFC (magnitude of downregulation)
    abs_logFC = abs(logFC),
    is_downregulated = logFC < 0
  )

# Filter to DOWNREGULATED transcripts only
df_down <- df_long %>% filter(is_downregulated)

cat(sprintf("Downregulated transcripts in NT: NMD=%d, Non-NMD=%d\n",
            sum(df_down$Condition == "NT" & df_down$NMD_status == "NMD-targeted"),
            sum(df_down$Condition == "NT" & df_down$NMD_status == "non-NMD-targeted")))
cat(sprintf("Downregulated transcripts in R0h: NMD=%d, Non-NMD=%d\n",
            sum(df_down$Condition == "R0h" & df_down$NMD_status == "NMD-targeted"),
            sum(df_down$Condition == "R0h" & df_down$NMD_status == "non-NMD-targeted")))

# Calculate p-values for downregulated transcripts (magnitude comparison)
pval_NT_down <- wilcox.test(abs_logFC ~ NMD_status, 
                            data = df_down %>% filter(Condition == "NT"))$p.value
pval_R0h_down <- wilcox.test(abs_logFC ~ NMD_status, 
                             data = df_down %>% filter(Condition == "R0h"))$p.value

cat(sprintf("Wilcoxon p-value (downregulation magnitude) NT: %.4g\n", pval_NT_down))
cat(sprintf("Wilcoxon p-value (downregulation magnitude) R0h: %.4g\n", pval_R0h_down))

# Summary of downregulation magnitude
down_summary <- df_down %>%
  group_by(Condition, NMD_status) %>%
  summarise(
    N = n(),
    Mean_magnitude = mean(abs_logFC),
    Median_magnitude = median(abs_logFC),
    .groups = "drop"
  )
print(down_summary)

# Create annotation dataframe for p-values
pval_labels_down <- data.frame(
  Condition = factor(c("NT", "R0h"), levels = c("NT", "R0h")),
  label = c(sprintf("p-value: %.2g", pval_NT_down), sprintf("p-value: %.4f", pval_R0h_down)),
  x = 1.5,
  y = quantile(df_down$abs_logFC, 0.98, na.rm = TRUE)
)

# Colors: orange for NMD-targeted, blue for non-NMD-targeted
colors_ref <- c("NMD-targeted" = "darkorange", "non-NMD-targeted" = "steelblue")

# ==============================================================================
# Figure 5A: Expression levels of downregulated transcripts upon OXR1 depletion
# ==============================================================================
# Comparing NMD-targeted (orange) vs non-NMD-targeted (blue)
# Two-sided Wilcoxon rank-sum tests

# Colors: orange for NMD-targeted, light blue for non-NMD-targeted (matching reference)
colors_fig5a <- c("NMD-targeted" = "#E69F00", "non-NMD-targeted" = "#56B4E9")

# Create p-value annotations
pval_labels_fig5a <- data.frame(
  Condition = factor(c("NT", "R0h"), levels = c("NT", "R0h")),
  label = c(sprintf("p-value: %.2f", pval_NT_down), sprintf("p-value: %.4f", pval_R0h_down)),
  x = 1.5,
  y = quantile(df_down$abs_logFC, 0.97, na.rm = TRUE)
)

p_fig5a <- ggplot(df_down, aes(x = NMD_status, y = abs_logFC, fill = NMD_status)) +
  geom_violin(alpha = 0.6, width = 0.85, trim = FALSE, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", linewidth = 0.4,
               fill = "white", alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 1.5, color = "black") +
  facet_wrap(~ Condition, nrow = 1) +
  geom_text(data = pval_labels_fig5a, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.8, fontface = "plain") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_fill_manual(values = colors_fig5a) +
  scale_x_discrete(labels = c("NMD-targeted" = "NMD-targeted", 
                               "non-NMD-targeted" = "non-NMD-targeted")) +
  labs(
    title = "Distribution of transcript-level expression differences upon OXR1 depletion",
    x = NULL,
    y = "Transcripts expression level difference"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain", hjust = 0.5, size = 12),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = 11),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    panel.spacing = unit(1, "lines")
  ) +
  coord_cartesian(ylim = c(0, quantile(df_down$abs_logFC, 0.995, na.rm = TRUE) * 1.1))

ggsave(file.path(output_dir, "Fig5a_oxr1_expression_downregulated.pdf"),
       plot = p_fig5a, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(output_dir, "Fig5a_oxr1_expression_downregulated.jpeg"),
       plot = p_fig5a, width = 7, height = 4.5, dpi = 450)
cat("Saved: Fig5a_oxr1_expression_downregulated.pdf/jpeg\n")

# Also save with original filename for backward compatibility
ggsave(file.path(output_dir, "Fig5a_oxr1_downregulation_magnitude.pdf"),
       plot = p_fig5a, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(output_dir, "Fig5a_oxr1_downregulation_magnitude.jpeg"),
       plot = p_fig5a, width = 7, height = 4.5, dpi = 450)
cat("Saved: Fig5a_oxr1_downregulation_magnitude.pdf/jpeg\n")

# ==============================================================================
# Figure 5A-2: Same but with NEGATIVE values (showing actual downregulation)
# ==============================================================================
# This shows the actual negative log fold change values (no absolute value)

# Calculate p-values for the signed (negative) values
pval_NT_signed <- wilcox.test(logFC ~ NMD_status, 
                              data = df_down %>% filter(Condition == "NT"))$p.value
pval_R0h_signed <- wilcox.test(logFC ~ NMD_status, 
                               data = df_down %>% filter(Condition == "R0h"))$p.value

cat(sprintf("\nWilcoxon p-value (signed logFC, downregulated only):\n"))
cat(sprintf("  NT: p = %.4g\n", pval_NT_signed))
cat(sprintf("  R0h: p = %.4g\n", pval_R0h_signed))

# Summary stats for signed values
signed_summary <- df_down %>%
  group_by(Condition, NMD_status) %>%
  summarise(
    N = n(),
    Mean_logFC = mean(logFC),
    Median_logFC = median(logFC),
    .groups = "drop"
  )
cat("\nSigned logFC summary (downregulated transcripts only):\n")
print(signed_summary)

# Create p-value annotations for signed plot
pval_labels_signed <- data.frame(
  Condition = factor(c("NT", "R0h"), levels = c("NT", "R0h")),
  label = c(sprintf("p-value: %.2f", pval_NT_signed), sprintf("p-value: %.4f", pval_R0h_signed)),
  x = 1.5,
  y = min(df_down$logFC, na.rm = TRUE) * 0.1
)

p_fig5a_signed <- ggplot(df_down, aes(x = NMD_status, y = logFC, fill = NMD_status)) +
  geom_violin(alpha = 0.6, width = 0.85, trim = FALSE, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", linewidth = 0.4,
               fill = "white", alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 16, size = 1.5, color = "black") +
  facet_wrap(~ Condition, nrow = 1) +
  geom_text(data = pval_labels_signed, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.8, fontface = "plain") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_fill_manual(values = colors_fig5a) +
  scale_x_discrete(labels = c("NMD-targeted" = "NMD-targeted", 
                               "non-NMD-targeted" = "non-NMD-targeted")) +
  labs(
    title = "Distribution of transcript-level expression differences upon OXR1 depletion",
    subtitle = "Downregulated transcripts (negative values indicate downregulation)",
    x = NULL,
    y = "Log2 Fold Change (siOXR1 / siControl)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "plain", hjust = 0.5, size = 12),
    plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 9, color = "gray40"),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = 11),
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    panel.spacing = unit(1, "lines")
  ) +
  coord_cartesian(ylim = c(min(df_down$logFC, na.rm = TRUE) * 1.05, 0.05))

ggsave(file.path(output_dir, "Fig5a_oxr1_downregulated_signed.pdf"),
       plot = p_fig5a_signed, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(output_dir, "Fig5a_oxr1_downregulated_signed.jpeg"),
       plot = p_fig5a_signed, width = 7, height = 4.5, dpi = 450)
cat("Saved: Fig5a_oxr1_downregulated_signed.pdf/jpeg\n")

# ==============================================================================
# Plot 3b: Proportion of DOWNREGULATED transcripts (bar plot)
# ==============================================================================
prop_down <- df_long %>%
  group_by(Condition, NMD_status) %>%
  summarise(
    Total = n(),
    N_down = sum(is_downregulated),
    Pct_down = N_down / Total * 100,
    .groups = "drop"
  )

# Chi-square test for proportion differences
chisq_NT <- chisq.test(matrix(c(
  prop_down$N_down[prop_down$Condition == "NT" & prop_down$NMD_status == "NMD-targeted"],
  prop_down$Total[prop_down$Condition == "NT" & prop_down$NMD_status == "NMD-targeted"] - 
    prop_down$N_down[prop_down$Condition == "NT" & prop_down$NMD_status == "NMD-targeted"],
  prop_down$N_down[prop_down$Condition == "NT" & prop_down$NMD_status == "non-NMD-targeted"],
  prop_down$Total[prop_down$Condition == "NT" & prop_down$NMD_status == "non-NMD-targeted"] - 
    prop_down$N_down[prop_down$Condition == "NT" & prop_down$NMD_status == "non-NMD-targeted"]
), nrow = 2))$p.value

chisq_R0h <- chisq.test(matrix(c(
  prop_down$N_down[prop_down$Condition == "R0h" & prop_down$NMD_status == "NMD-targeted"],
  prop_down$Total[prop_down$Condition == "R0h" & prop_down$NMD_status == "NMD-targeted"] - 
    prop_down$N_down[prop_down$Condition == "R0h" & prop_down$NMD_status == "NMD-targeted"],
  prop_down$N_down[prop_down$Condition == "R0h" & prop_down$NMD_status == "non-NMD-targeted"],
  prop_down$Total[prop_down$Condition == "R0h" & prop_down$NMD_status == "non-NMD-targeted"] - 
    prop_down$N_down[prop_down$Condition == "R0h" & prop_down$NMD_status == "non-NMD-targeted"]
), nrow = 2))$p.value

cat(sprintf("\nProportion downregulated - Chi-square p-values:\n"))
cat(sprintf("  NT: p = %.4g\n", chisq_NT))
cat(sprintf("  R0h: p = %.4g\n", chisq_R0h))
print(prop_down)

pval_labels_prop <- data.frame(
  Condition = factor(c("NT", "R0h"), levels = c("NT", "R0h")),
  label = c(sprintf("p = %.2g", chisq_NT), sprintf("p = %.2g", chisq_R0h)),
  x = 1.5,
  y = max(prop_down$Pct_down) + 3
)

p3b <- ggplot(prop_down, aes(x = NMD_status, y = Pct_down, fill = NMD_status)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", Pct_down)), vjust = -0.5, size = 3.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  facet_wrap(~ Condition, nrow = 1) +
  geom_text(data = pval_labels_prop, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "italic") +
  scale_fill_manual(values = colors_ref) +
  labs(
    title = "Proportion of transcripts downregulated upon OXR1 depletion",
    x = NULL,
    y = "Percent downregulated (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white")
  ) +
  scale_y_continuous(limits = c(0, max(prop_down$Pct_down) + 8), expand = c(0, 0))

ggsave(file.path(output_dir, "Fig5a_oxr1_proportion_downregulated.pdf"),
       plot = p3b, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_dir, "Fig5a_oxr1_proportion_downregulated.jpeg"),
       plot = p3b, width = 8, height = 5, dpi = 450)
cat("Saved: Fig5a_oxr1_proportion_downregulated.pdf/jpeg\n")

# ==============================================================================
# Plot 3c: Combined view - all transcripts but colored by direction
# ==============================================================================
df_long_dir <- df_long %>%
  mutate(Direction = ifelse(logFC < 0, "Downregulated", "Upregulated"))

# Mean logFC by group
mean_logFC <- df_long_dir %>%
  group_by(Condition, NMD_status) %>%
  summarise(Mean_logFC = mean(logFC), .groups = "drop")

p3c <- ggplot(df_long_dir, aes(x = NMD_status, y = logFC, fill = NMD_status)) +
  geom_violin(alpha = 0.7, width = 0.9, trim = FALSE, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  geom_text(data = mean_logFC, aes(x = NMD_status, y = -2.5, 
                                    label = sprintf("mean=%.3f", Mean_logFC)),
            size = 3, color = "red", fontface = "bold") +
  facet_wrap(~ Condition, nrow = 1) +
  scale_fill_manual(values = colors_ref) +
  labs(
    title = "Transcript expression changes upon OXR1 depletion",
    subtitle = "Red diamond = mean; values below 0 = downregulation",
    x = NULL,
    y = "Log2 Fold Change (siOXR1 / siControl)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white")
  )

ggsave(file.path(output_dir, "Fig5a_oxr1_logFC_by_condition.pdf"),
       plot = p3c, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_dir, "Fig5a_oxr1_logFC_by_condition.jpeg"),
       plot = p3c, width = 8, height = 5, dpi = 450)
cat("Saved: Fig5a_oxr1_logFC_by_condition.pdf/jpeg\n")

# ==============================================================================
# Plot 4: Distribution histogram
# ==============================================================================
p4 <- ggplot(df, aes(x = logFC_mean, fill = NMD_status)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_vline(data = summary_stats, aes(xintercept = Mean_logFC, color = NMD_status),
             linetype = "solid", linewidth = 1) +
  scale_fill_manual(values = colors, labels = c("Non-NMD", "NMD")) +
  scale_color_manual(values = colors, guide = "none") +
  labs(
    title = "Distribution of Log Fold Changes Upon OXR1 Knockdown",
    subtitle = sprintf("NMD mean = %.3f, Non-NMD mean = %.3f", 
                       summary_stats$Mean_logFC[2], summary_stats$Mean_logFC[1]),
    x = "Log2 Fold Change (siOXR1 / siControl)",
    y = "Density",
    fill = "Transcript Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "Fig5a_oxr1_logFC_density.pdf"),
       plot = p4, width = 8, height = 5, dpi = 300)
ggsave(file.path(output_dir, "Fig5a_oxr1_logFC_density.jpeg"),
       plot = p4, width = 8, height = 5, dpi = 450)

cat("Saved: Fig5a_oxr1_logFC_density.pdf/jpeg\n")

# ==============================================================================
# Save summary data
# ==============================================================================
write.csv(summary_stats, file.path(output_dir, "Fig5a_oxr1_summary_stats.csv"), 
          row.names = FALSE)

# Save full data with logFC
df_export <- df %>%
  select(Transcript_ID, gene_name, transcript_type, NMD, 
         siCon_NT, siCon_R0h, siOXR1_NT, siOXR1_R0h,
         logFC_NT, logFC_R0h, logFC_mean)
write.csv(df_export, file.path(output_dir, "Fig5a_oxr1_logFC_data.csv"), 
          row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Output saved to: %s\n", output_dir))
cat("\nKey findings:\n")
cat(sprintf("- NMD transcripts: mean logFC = %.4f (direction: %s)\n", 
            summary_stats$Mean_logFC[2],
            ifelse(summary_stats$Mean_logFC[2] < 0, "DOWNREGULATED", "UPREGULATED")))
cat(sprintf("- Non-NMD transcripts: mean logFC = %.4f (direction: %s)\n",
            summary_stats$Mean_logFC[1],
            ifelse(summary_stats$Mean_logFC[1] < 0, "DOWNREGULATED", "UPREGULATED")))
cat(sprintf("- Wilcoxon test p-value: %.2e\n", wilcox_test$p.value))

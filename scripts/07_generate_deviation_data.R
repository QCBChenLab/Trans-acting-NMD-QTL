#!/usr/bin/env Rscript
# ==============================================================================
# Generate NMD Deviation Data for All Overlapped Genes
# ==============================================================================
# This script generates plot_data_*.csv files for all genes that overlap
# between NMD-QTL mapping and GTEx expression data.
# ==============================================================================

library(tidyverse)
library(readxl)
library(data.table)

# ==============================================================================
# Configuration
# ==============================================================================
# BASE_DIR: path to the main project directory containing Cancer_KO/, RBP_map/, Manuscript/ etc.
# Adjust this path according to your local setup.
BASE_DIR <- Sys.getenv("NMD_PROJECT_DIR", unset = normalizePath("../../..", mustWork = FALSE))

output_dir <- file.path(BASE_DIR, "Cancer_KO/Boxplots/NMD_Deviation_by_ExprGroup/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# All 81 overlapped genes (NMD-mapped genes present in GTEx)
all_genes <- c(
  "AC004012.1", "AC104088.1", "ACVR1B", "ADGRV1", "AKR1B15", "ANGPTL5", 
  "ANK2", "ASB4", "BCHE", "BPGM", "C1QBPP1", "CCDC162P", "CCDC26", 
  "CERS6", "CFAP54", "CIART", "CNIH3", "COLCA2", "CXCR5", "DGKG", 
  "DHRSX", "DIAPH2", "EDN2", "EEF1A1P17", "EVC2", "FAM110C", "FAM155A", 
  "FARS2", "FOXN3", "GEMIN8", "GLP1R", "HNRNPCP10", "IL1RAP", "ITPK1", 
  "KLF8", "KSR1P1", "LINC00276", "LINC00393", "LINC02020", "LINC02054", 
  "MAOA", "MEMO1P4", "MIR2054", "MIR4491", "MPPE1P1", "MRPS21", 
  "MTHFD2P1", "MTND4P17", "MYLK", "NANS", "NDUFB4P12", "NRG1", "OXR1", 
  "PABPC5-AS1", "PRKCE", "PTCSC2", "PUS7L", "RAB22A", "RN7SKP82", 
  "RN7SL193P", "RN7SL691P", "RNA5SP45", "RNF44", "RNU6-102P", 
  "RNU6-215P", "RNU6-918P", "RPL7AP30", "RRM2P4", "SDC1", "SLC7A15P", 
  "SRGAP1", "STPG1", "THOC2", "TMEM117", "TMEM179", "TSR2", "UBQLN2", 
  "VDR", "WNK3", "YEATS2", "ZNF890P"
)

cat(sprintf("Processing %d genes...\n", length(all_genes)))

# ==============================================================================
# 1. Load mapping tables to get gene-tissue relationships
# ==============================================================================
cat("Loading mapping tables...\n")
table_s1 <- read_excel(file.path(BASE_DIR, "Manuscript/Tables.xlsx"), sheet = "Table S1")
table_s2 <- read_excel(file.path(BASE_DIR, "Manuscript/Tables.xlsx"), sheet = "Table S2")[-1, ]

# Create gene-SNP mapping for ALL genes (not just cancer_ko_genes)
gene_snp_map <- table_s2 %>%
  mutate(
    Gene_Overlapped = `Overlapped Gene`,
    Gene_Upstream = `Nearest Upstream Gene`,
    Gene_Downstream = `Nearest Downstream Gene`
  ) %>%
  pivot_longer(
    cols = c(Gene_Overlapped, Gene_Upstream, Gene_Downstream),
    names_to = "Gene_Type",
    values_to = "Gene"
  ) %>%
  filter(Gene %in% all_genes) %>%
  transmute(Gene, SNP_ID = `Variation ID`) %>%
  distinct()

# Format Table S1 SNP IDs to match
table_s1_formatted <- table_s1 %>%
  mutate(
    SNP_ID = paste0(
      "chr", CHR, ":", Position, ":",
      sub(".*_([ACGT])_([ACGT])_.*", "\\1/\\2", SNP),
      ":1"
    )
  )

# Create gene-tissue mapping
gene_tissue_map <- gene_snp_map %>%
  inner_join(
    table_s1_formatted %>% select(SNP_ID, Tissue = tissue),
    by = "SNP_ID"
  ) %>%
  select(Gene, Tissue) %>%
  distinct()

cat(sprintf("Found tissue mappings for %d genes\n", n_distinct(gene_tissue_map$Gene)))

# ==============================================================================
# 2. Load GTEx gene expression data
# ==============================================================================
cat("Loading GTEx gene expression data...\n")
gtex_exp_genes <- read_tsv(
  file.path(BASE_DIR, "RBP_map/Gene_Expression/GTEx_mutation/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"),
  skip = 2,
  show_col_types = FALSE
) %>%
  filter(Description %in% all_genes)

cat(sprintf("Found expression data for %d genes in GTEx\n", nrow(gtex_exp_genes)))

# ==============================================================================
# 3. Load GTEx NMD transcript expression data
# ==============================================================================
cat("Loading GTEx NMD transcript data (this may take a while)...\n")
gtex_nmd_tx_exp <- fread(
file.path(BASE_DIR, "RBP_map/Gene_Expression/GTEx_mutation/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_nmd_transcript_tpm.gct"),
  header = FALSE, sep = "\t"
)

# ==============================================================================
# 4. Load sample attributes
# ==============================================================================
sample_attr <- read.delim(
  file.path(BASE_DIR, "RBP_map/Gene_Expression/GTEx_mutation/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
  header = TRUE, sep = "\t"
) %>%
  select(SAMPID, Tissue = SMTSD)

# ==============================================================================
# 5. Reshape gene expression to long format
# ==============================================================================
cat("Reshaping gene expression data...\n")
gene_long <- gtex_exp_genes %>%
  pivot_longer(
    cols = -(Name:Description),
    names_to = "SAMPID",
    values_to = "GeneExp"
  ) %>%
  transmute(
    Gene = Description,
    SAMPID,
    Gene_log = log1p(GeneExp)
  )

# ==============================================================================
# 6. Reshape NMD data & compute RINT deviation
# ==============================================================================
cat("Computing RINT deviations for NMD transcripts...\n")

# Assign sample IDs to NMD matrix
nmd_sample_cols <- colnames(gtex_nmd_tx_exp)[-c(1, 2)]
sample_ids <- unique(gene_long$SAMPID)
colnames(gtex_nmd_tx_exp) <- c(
  "Transcript_ID", "Gene_ID", sample_ids[seq_along(nmd_sample_cols)]
)

# Convert to long format
nmd_long <- gtex_nmd_tx_exp %>%
  select(-Gene_ID) %>%
  pivot_longer(
    cols = -Transcript_ID,
    names_to = "SAMPID",
    values_to = "NMD_Exp"
  ) %>%
  inner_join(sample_attr, by = "SAMPID")

# Compute per-transcript per-tissue RINT + deviation
cat("Computing RINT (this is the slowest step)...\n")
nmd_rint <- nmd_long %>%
  group_by(Transcript_ID, Tissue) %>%
  mutate(
    N = sum(!is.na(NMD_Exp)),
    Rank = rank(NMD_Exp, ties.method = "average", na.last = "keep"),
    RINT = qnorm((Rank - 0.375) / (N + 0.25)),
    Mean_RINT = mean(RINT, na.rm = TRUE),
    Deviation = RINT - Mean_RINT
  ) %>%
  ungroup()

rm(nmd_long, gtex_nmd_tx_exp)  # Free memory
gc()

# ==============================================================================
# 7. Generate plot_data files for each gene
# ==============================================================================
cat("\nGenerating plot_data files for each gene...\n")

genes_with_tissue <- unique(gene_tissue_map$Gene)
genes_with_expr <- unique(gene_long$Gene)
genes_to_process <- intersect(genes_with_tissue, genes_with_expr)

cat(sprintf("Will process %d genes with both tissue mapping and expression data\n", 
            length(genes_to_process)))

results_summary <- data.frame(
  Gene = character(),
  N_samples = integer(),
  N_observations = integer(),
  Tissues = character(),
  Status = character(),
  stringsAsFactors = FALSE
)

for (gene_symbol in genes_to_process) {
  
  # Get tissues for this gene
  tissues <- gene_tissue_map %>%
    filter(Gene == gene_symbol) %>%
    pull(Tissue) %>%
    unique()
  
  if (length(tissues) == 0) {
    results_summary <- rbind(results_summary, data.frame(
      Gene = gene_symbol, N_samples = 0, N_observations = 0,
      Tissues = "", Status = "No tissue mapping"
    ))
    next
  }
  
  # Get gene expression data
  gene_data <- gene_long %>% filter(Gene == gene_symbol)
  if (nrow(gene_data) == 0) {
    results_summary <- rbind(results_summary, data.frame(
      Gene = gene_symbol, N_samples = 0, N_observations = 0,
      Tissues = paste(tissues, collapse = "; "), Status = "No expression data"
    ))
    next
  }
  
  # Join with NMD deviation data
  dat <- gene_data %>%
    inner_join(nmd_rint %>% filter(Tissue %in% tissues), by = "SAMPID") %>%
    filter(!is.na(Gene_log), !is.na(Deviation))
  
  if (nrow(dat) == 0) {
    results_summary <- rbind(results_summary, data.frame(
      Gene = gene_symbol, N_samples = 0, N_observations = 0,
      Tissues = paste(tissues, collapse = "; "), Status = "No overlapping data"
    ))
    next
  }
  
  # Group into tertiles per tissue
  dat <- dat %>%
    group_by(Tissue, Transcript_ID) %>%
    mutate(
      ExprGroup = ntile(Gene_log, 3),
      ExprGroup = factor(ExprGroup, levels = 1:3,
                         labels = c("Low", "Medium", "High")),
      NumGroup = as.numeric(ExprGroup)
    ) %>%
    ungroup()
  
  # Export CSV
  csv_path <- file.path(output_dir, paste0("plot_data_", gene_symbol, ".csv"))
  write.csv(dat, csv_path, row.names = FALSE)
  
  # Compute trend statistics
  trend_stats <- dat %>%
    group_by(Tissue) %>%
    summarise(
      Slope = coef(lm(Deviation ~ NumGroup))[2],
      P_value = summary(lm(Deviation ~ NumGroup))$coefficients[2, 4],
      N = n(),
      .groups = "drop"
    )
  
  # Save trend results
  txt_path <- file.path(output_dir, paste0("trend_results_", gene_symbol, ".txt"))
  write.table(trend_stats, txt_path, row.names = FALSE, quote = FALSE, sep = "\t")
  
  results_summary <- rbind(results_summary, data.frame(
    Gene = gene_symbol,
    N_samples = n_distinct(dat$SAMPID),
    N_observations = nrow(dat),
    Tissues = paste(tissues, collapse = "; "),
    Status = "Success"
  ))
  
  cat(sprintf("  %s: %d samples, %d observations\n", 
              gene_symbol, n_distinct(dat$SAMPID), nrow(dat)))
}

# ==============================================================================
# 8. Summary
# ==============================================================================
cat("\n=== Summary ===\n")
cat(sprintf("Total genes attempted: %d\n", length(genes_to_process)))
cat(sprintf("Successfully processed: %d\n", sum(results_summary$Status == "Success")))
cat(sprintf("Failed: %d\n", sum(results_summary$Status != "Success")))

# Save summary
write.csv(results_summary, file.path(output_dir, "processing_summary.csv"), row.names = FALSE)

# List successful genes
successful_genes <- results_summary %>%
  filter(Status == "Success") %>%
  pull(Gene)

cat("\nSuccessful genes:\n")
cat(paste(successful_genes, collapse = ", "), "\n")

cat(sprintf("\n=== Output saved to: %s ===\n", output_dir))

# ==============================================================================
# SNP Annotation Analysis
# This script annotates significant SNPs using biomaRt and other databases
# ==============================================================================

# Load data and functions
source("00_data_loading.R")

cat("Starting SNP annotation analysis...\n")

# ==============================================================================
# SNP Data Preparation
# ==============================================================================

# Clean up significant SNPs data
sig_SNP_cln <- sig_SNP %>% 
  separate(pos, c("CHR", "pos"), sep = ":") %>% 
  mutate(CHR = gsub("chr", "", CHR), 
         pos = as.numeric(pos)) %>% 
  mutate(SNP = variant_id, 
         start = pos, 
         end = pos, 
         variant_pos = paste0(CHR, ":", pos)) %>% 
  dplyr::select(SNP, CHR, start, end, variant_pos, p_val, n, tissue)

# Create GRanges object
sig_SNP_cln_gr <- makeGRangesFromDataFrame(sig_SNP_cln, keep.extra.columns = TRUE)

cat("✓ SNP data prepared\n")

# ==============================================================================
# Reference SNP Mapping
# ==============================================================================

# Load reference SNPs and map to significant SNPs
annotate_with_rsid <- function(snp_gr) {
  tryCatch({
    ref_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
    known_snps <- snpsByOverlaps(ref_snps, snp_gr)
    
    # Create mapping table
    sig_SNP_rsid <- sig_SNP_cln %>% 
      left_join(as.data.frame(known_snps) %>% 
                  mutate(variant_pos = paste0(seqnames, ":", pos), 
                         rsid = RefSNP_id) %>% 
                  dplyr::select(variant_pos, rsid), 
                by = "variant_pos")
    
    return(sig_SNP_rsid)
  }, error = function(e) {
    cat("Warning: Could not load reference SNPs:", e$message, "\n")
    return(sig_SNP_cln)
  })
}

sig_SNP_rsid <- annotate_with_rsid(sig_SNP_cln_gr)
cat("✓ Reference SNP mapping completed\n")

# ==============================================================================
# Gene Annotation with biomaRt
# ==============================================================================

# Function to annotate SNPs with genes using biomaRt
annotate_snps_with_genes <- function(rsid_list) {
  tryCatch({
    # Set up biomaRt connections
    mart_snp <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
    mart_gene <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    
    cat("Connecting to Ensembl databases...\n")
    
    # Get SNP-gene mappings
    sig_snp_gene_mapping <- getBM(
      attributes = c('refsnp_id', 'chr_name', 
                     'chrom_start', 'chrom_end', 
                     'ensembl_gene_stable_id', 'ensembl_gene_name'),
      filters = 'snp_filter',
      values = rsid_list,
      mart = mart_snp
    )
    
    # Get gene symbols
    sig_gene_symbol_mapping <- getBM(
      attributes = c('ensembl_gene_id', 'external_gene_name'),
      filters = 'ensembl_gene_id',
      values = sig_snp_gene_mapping %>% 
        filter(ensembl_gene_stable_id != "") %>%
        pull(ensembl_gene_stable_id),
      mart = mart_gene
    )
    
    # Combine annotations
    sig_SNP_anno_res <- sig_SNP_rsid %>%
      left_join(sig_snp_gene_mapping %>% 
                  dplyr::filter(!str_starts(chr_name, "HSC")) %>%
                  dplyr::select(refsnp_id, ensembl_gene_stable_id),
                by = c("rsid" = "refsnp_id")) %>%
      left_join(sig_gene_symbol_mapping, 
                by = c("ensembl_gene_stable_id" = "ensembl_gene_id")) %>% 
      distinct()
    
    return(sig_SNP_anno_res)
    
  }, error = function(e) {
    cat("Error in biomaRt annotation:", e$message, "\n")
    cat("Returning basic annotation...\n")
    return(sig_SNP_rsid)
  })
}

# Perform gene annotation
if("rsid" %in% colnames(sig_SNP_rsid) && any(!is.na(sig_SNP_rsid$rsid))) {
  valid_rsids <- sig_SNP_rsid$rsid[!is.na(sig_SNP_rsid$rsid)]
  sig_SNP_anno_res <- annotate_snps_with_genes(valid_rsids)
  cat("✓ Gene annotation completed\n")
} else {
  sig_SNP_anno_res <- sig_SNP_rsid
  cat("⚠ No valid RSIDs found for annotation\n")
}

# ==============================================================================
# Export Results
# ==============================================================================

# Save annotation results
write.csv(sig_SNP_anno_res, "../output/tables/sig_SNP_annotation_biomart.csv", row.names = FALSE)

# ==============================================================================
# Create Query File for External Tools
# ==============================================================================

# Prepare file for SNP-Nexus annotation
query_sig_snp <- sig_SNP %>% 
  separate(pos, c("CHR", "pos"), sep = ":") %>% 
  mutate(CHR = gsub("chr", "", CHR), 
         pos = as.numeric(pos), 
         chr = "chromosome", 
         strand = 1) %>% 
  dplyr::select(chr, CHR, pos, A1, A2, strand)

write.table(query_sig_snp, "../output/tables/query_sig_snp_for_snpnexus.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("✓ SNP annotation analysis complete!\n")
cat("✓ Annotation results saved to output/tables/\n")
cat("✓ Query file for SNP-Nexus created\n") 
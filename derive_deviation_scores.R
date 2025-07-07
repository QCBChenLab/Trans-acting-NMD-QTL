library(CePa)
library(edgeR)
library(RNOmni)


tx_read_tpm <- read.gct("~/GTEx_Analysis_transcript_tpm.gct")

sample_attributes <- read.delim("~/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                                header = TRUE)

load('~/data/NMD_gene_list.RData')

print('Finished loading files.')


read_tpm_sample_id <- colnames(tx_read_tpm)
read_tpm_sample_id_mod <- unlist(lapply(read_tpm_sample_id, function(x) {
  paste0(unlist(strsplit(x,'[.]')), collapse = '-')
}))
colnames(tx_read_tpm) <- read_tpm_sample_id_mod

read_tx_id <- rownames(tx_read_tpm)
# load("C:/Users/zheyuli/research/NMD/data/GTEx_read_tpm_tx_id.RData")
read_tx_id_mod <- unlist(lapply(read_tx_id, function(x) {
  unlist(strsplit(x,'[.]'))[1]
}))
rownames(tx_read_tpm) <- read_tx_id_mod

redundant_tx_id <- read_tx_id_mod[duplicated(read_tx_id_mod)]
clean_tx_id <- read_tx_id_mod[!read_tx_id_mod %in% redundant_tx_id]
complete_NMD_gene_ind <- c()
for (i in 1:length(NMD_gene_tx_list_clean)) {
  tx_within_gene <- unlist(NMD_gene_tx_list_clean[[i]])
  if (sum(tx_within_gene %in% clean_tx_id) == length(tx_within_gene)) {
    complete_NMD_gene_ind <- c(complete_NMD_gene_ind, i)
  }
}
NMD_gene_tx_list_clean_2 <- NMD_gene_tx_list_clean[complete_NMD_gene_ind]

tx_read_tpm_clean <- tx_read_tpm[read_tx_id_mod %in% clean_tx_id,]
read_tpm_sample_id_clean <- rownames(tx_read_tpm_clean)


tissues <- unique(sample_attributes$SMTSD)

output_dir <- '~/cmb0/NMD/data/deviation_score'
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

output_dir2 <- '~/data/expression_by_tissue'
if (!dir.exists(output_dir2)) {
  dir.create(output_dir2)
}

for (i in 1:length(tissues)) { 
  # extract expression values per tissue
  tissue <- tissues[i]
  print(tissue)
  GTEx_sample_id <- sample_attributes$SAMPID[sample_attributes$SMTSD==tissue]
  # tissue_ct <- t(gene_read_ct[,read_ct_sample_id%in%GTEx_sample_id])
  tissue_tpm <- t(tx_read_tpm_clean[,read_tpm_sample_id_mod %in% GTEx_sample_id])
  if (nrow(tissue_tpm)  == 0) {
    next
  }
  
  # # find cols with not enough ct
  # col_ind_valid_ct <- c()
  # for (j in 1:ncol(tissue_ct)) {
  #   valid_ct_num <- sum(tissue_ct[,j]>=6)
  #   if (valid_ct_num > nrow(tissue_ct)*.2) {
  #     col_ind_valid_ct <- c(col_ind_valid_ct,j)
  #   }
  # }
  
  # # find cols with not enough tpm
  # col_ind_valid_tpm <- c()
  # for (j in 1:ncol(tissue_tpm)) {
  #   valid_ct_num <- sum(tissue_tpm[,j]>.1)
  #   if (valid_ct_num > nrow(tissue_tpm)*.2) {
  #     col_ind_valid_tpm <- c(col_ind_valid_tpm,j)
  #   }
  # }
  # 
  # # remove genes w/o enough data
  # col_ind_valid <- union(col_ind_valid_ct,col_ind_valid_tpm)
  # tissue_tpm_clean <- tissue_tpm[,col_ind_valid]
  tissue_tpm_clean <- tissue_tpm[,colSums(tissue_tpm) != 0]
  
  
  # normalization
  # TMM
  d <- DGEList(counts=tissue_tpm_clean)
  TMM <- calcNormFactors(d, method="TMM")
  tissue_tpm_TMM <- cpm(TMM)
  
  donor_id <- rownames(tissue_tpm_TMM)
  donor_id_mod <- unlist(lapply(donor_id, function(x) {
    unlist(strsplit(x,'-'))[2] 
  }))
  
  tissue_tx_id <- colnames(tissue_tpm_TMM)
  # load("C:/Users/zheyuli/research/NMD/data/GTEx_read_tpm_tx_id.RData")
  # tissue_tx_id_mod <- unlist(lapply(tissue_tx_id, function(x) {
  #   unlist(strsplit(x,'[.]'))[1]
  # }))
  # rownames(tx_read_tpm) <- read_tpm_tx_id_mod
  
  NMD_expression <- data.frame(matrix(nrow = length(donor_id_mod),
                                      ncol = length(NMD_gene_tx_list_clean_2)))
  rownames(NMD_expression) <- donor_id_mod
  colnames(NMD_expression) <- names(NMD_gene_tx_list_clean_2)
  non_NMD_expression <- NMD_expression
  
  bad_gene_ind <- c()
  for (j in 1:length(NMD_gene_tx_list_clean_2)) {
    gene_id <- names(NMD_gene_tx_list_clean_2)[j]
    tx_list <- NMD_gene_tx_list_clean_2[[gene_id]]
    
    NMD_tpm_TMM <- c()
    tx_id_NMD <- tx_list$tx_id_NMD
    hit_tx_id_NMD <- c()
    for (k in 1:length(tx_id_NMD)) {
      tx_id <- tx_id_NMD[k]
      tx_col_ind <- match(tx_id, tissue_tx_id)
      if (!is.na(tx_col_ind)) {
        NMD_tpm_TMM <- cbind(NMD_tpm_TMM, tissue_tpm_TMM[, tx_col_ind])
        hit_tx_id_NMD <- c(hit_tx_id_NMD, tx_id)
      }
    }

    if (is.null(NMD_tpm_TMM)) {
      bad_gene_ind <- c(bad_gene_ind, j)
      next
    }
    colnames(NMD_tpm_TMM) <- hit_tx_id_NMD
    
    # donor_id_2 <- unlist(lapply(rownames(NMD_tpm_TMM), function(x) {
    #   unlist(strsplit(x,'-'))[2]
    # }))
    NMD_expression[,j] <- rowSums(NMD_tpm_TMM)
    
    
    non_NMD_tpm_TMM <- c()
    tx_id_non_NMD <- tx_list$tx_id_non_NMD
    hit_tx_id_non_NMD <- c()
    for (k in 1:length(tx_id_non_NMD)) {
      tx_id <- tx_id_non_NMD[k]
      tx_col_ind <- match(tx_id, tissue_tx_id)
      if (!is.na(tx_col_ind)) {
        non_NMD_tpm_TMM <- cbind(non_NMD_tpm_TMM, tissue_tpm_TMM[, tx_col_ind])
        hit_tx_id_non_NMD <- c(hit_tx_id_non_NMD, tx_id)
      }
    }

    if (is.null(non_NMD_tpm_TMM)) {
      bad_gene_ind <- c(bad_gene_ind, j)
      next
    }
    colnames(non_NMD_tpm_TMM) <- hit_tx_id_non_NMD
    
    non_NMD_expression[,j] <- rowSums(non_NMD_tpm_TMM)
    
  }
  
  file_path2 <- file.path(output_dir2, paste0(tissue, '.RData'))
  save(list = c('NMD_expression', 'non_NMD_expression', 'bad_gene_ind'), 
       file = file_path2)
  
  NMD_expression_clean <- NMD_expression[, -bad_gene_ind]
  non_NMD_expression_clean <- non_NMD_expression[, -bad_gene_ind]
  
  NMD_expression_INT <- c()
  # circ_gene_id <- colnames(circadian_gene_tpm)
  for (j in 1:ncol(NMD_expression_clean)) {
    new_col <- RankNorm(NMD_expression_clean[,j], k = 0.375)
    # weight <- circ_gene_by_tissue_clean$log_p[circ_gene_by_tissue_clean$gene_id==circ_gene_id[j]]
    deviation_col <- abs(new_col - mean(new_col))
    NMD_expression_INT <- cbind(NMD_expression_INT,deviation_col)
  }
  rownames(NMD_expression_INT) <- rownames(NMD_expression_clean)
  colnames(NMD_expression_INT) <- colnames(NMD_expression_clean)
  
  NMD_deviation <- rowSums(NMD_expression_INT) / ncol(NMD_expression_INT)
  
  non_NMD_expression_INT <- c()
  # circ_gene_id <- colnames(circadian_gene_tpm)
  for (j in 1:ncol(non_NMD_expression_clean)) {
    new_col <- RankNorm(non_NMD_expression_clean[,j], k = 0.375)
    # weight <- circ_gene_by_tissue_clean$log_p[circ_gene_by_tissue_clean$gene_id==circ_gene_id[j]]
    deviation_col <- abs(new_col - mean(new_col))
    non_NMD_expression_INT <- cbind(non_NMD_expression_INT,deviation_col)
  }
  rownames(non_NMD_expression_INT) <- rownames(non_NMD_expression_clean)
  colnames(non_NMD_expression_INT) <- colnames(non_NMD_expression_clean)
  
  non_NMD_deviation <- rowSums(non_NMD_expression_INT) / ncol(non_NMD_expression_INT)
  
  deviation_score_tissue <- list(NMD_deviation = NMD_deviation,
                                 non_NMD_deviation = non_NMD_deviation)
  file_path <- file.path(output_dir, paste0(tissue, '.RData'))
  save('deviation_score_tissue', file = file_path)
}
  
  
  

  
  
  
  

  
  
  












#!/usr/bin/env Rscript

# Load libraries
library(coloc)
library(dplyr)

# Parse command line arguments
args <- commandArgs(trailingOnly=TRUE)
gene <- args[1]
merged_data <- args[2]
gwasN <- as.integer(args[3])
gwasType <- args[4]
cellType <- args[5]
gwas <- args[6]
eqtlN <- as.integer(args[7])

# Load merged data
merged_df <- read.csv(merged_data, sep = "\t")
print(str(merged_df))

# Prepare GWAS data
gwas_data <- list(
  beta = merged_df$beta_gwas,
  varbeta = merged_df$varbeta_gwas,
  snp = merged_df$snp_id_eqtl,
  position = merged_df$position_gwas,
  type = gwasType,
  N = as.integer(gwasN)
)

# If GWAS type is quant, add MAF column to GWAS data object
if (gwasType == 'quant') {
  gwas_data$MAF <- merged_df$maf_gwas
}
print(str(gwas_data))

# Prepare eQTL data
eqtl_data <- list(
    beta = merged_df$beta_eqtl,
    varbeta = merged_df$varbeta_eqtl,
    snp = merged_df$snp_id_eqtl,
    position = merged_df$snp_position_eqtl,
    type = "quant",
    N = eqtlN,
    MAF = merged_df$maf_eqtl
)
print(str(eqtl_data))

# Perform colocalization
res <- coloc.abf(dataset1=eqtl_data, dataset2=gwas_data)

# Create summary dataframe
summary_df <- data.frame(
  gene = gene,
  PPH0 = as.numeric(res$summary[2]),
  PPH1 = as.numeric(res$summary[3]),
  PPH2 = as.numeric(res$summary[4]),
  PPH3 = as.numeric(res$summary[5]),
  PPH4 = as.numeric(res$summary[6]),
  cell_type = cellType,
  gwas = gwas,
  stringsAsFactors = FALSE
)

# Write results to file
write.table(
  summary_df,
  file = paste0(gwas, "_", cellType, "_", gene, "_sv_coloc_out.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

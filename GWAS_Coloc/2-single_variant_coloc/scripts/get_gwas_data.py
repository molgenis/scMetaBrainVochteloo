#!/usr/bin/env python3

import pandas as pd 
import tabix
import gzip
import sys

# Parse command line arguments
gene = sys.argv[1]
chromosome = sys.argv[2]
gwas = sys.argv[3]
gwas_path = sys.argv[4]
pval_col = sys.argv[5]
beta_col = sys.argv[6]
se_col = sys.argv[7]
pos_col = sys.argv[8]
chrom_col = sys.argv[9]
ea_col = sys.argv[10]
oa_col = sys.argv[11]
snp_col = sys.argv[12]
cell_type = sys.argv[13]
annotation_path = sys.argv[14]

# Load annotation data
annotation_df = pd.read_csv(annotation_path, sep='\t', index_col=0)
print(annotation_df.head())

# Define locus start and end positions
locus_start = annotation_df.loc[gene]['ChrStart'] - 500000
locus_end = annotation_df.loc[gene]['ChrStart'] + 500000

tb = tabix.open(gwas_path)

# Retrieve variants from sum stats
records = tb.query(str(chromosome), locus_start, locus_end)
gwas_df = pd.DataFrame(list(records))
print(gwas_df)

# Get columns names
with gzip.open(gwas_path, 'rt') as gwas_file:
    for line in gwas_file:
        line = line.strip()
        columns = line.split('\t')
        break

gwas_df.columns = columns

# Set data type to numeric
gwas_df[beta_col] = pd.to_numeric(gwas_df[beta_col], errors='coerce')
gwas_df[se_col] = pd.to_numeric(gwas_df[se_col], errors='coerce')
gwas_df[pval_col] = pd.to_numeric(gwas_df[pval_col], errors='coerce')
gwas_df[pos_col] = pd.to_numeric(gwas_df[pos_col], errors='coerce')

# Drop rows with missing values in required colummns
gwas_df = gwas_df.dropna(subset=[
    pval_col,
    beta_col,
    se_col,
    pos_col,
    chrom_col,
    ea_col,
    oa_col,
    snp_col,
])

# Remove rows with a standard error of exactly 0; these will cause problems
gwas_df = gwas_df[gwas_df[se_col] != 0]

# Make sure alleles are uppercase
gwas_df[ea_col] = gwas_df[ea_col].str.upper()
gwas_df[oa_col] = gwas_df[oa_col].str.upper()

def format_snp_id_gwas(row):
    '''Function that creates a variant ID with alphabetically sorted alleles with format 1:12345:A,C'''
    
    chrom = row[chrom_col]
    pos = row[pos_col]
    ea = row[ea_col]
    oa = row[oa_col]
    alleles = sorted([ea, oa])

    return f'{chrom}:{pos}:{alleles[0]},{alleles[1]}'

# Calculate varbeta
gwas_df['varbeta'] = gwas_df[se_col] ** 2

gwas_df['snp_formatted_sorted'] = gwas_df.apply(format_snp_id_gwas, axis=1)
print(gwas_df.head())
print(gwas_df.shape)

# Rename columns
gwas_df.rename(columns={
    pval_col: 'p_value',
    beta_col: 'beta',
    se_col: 'standard_error',
    pos_col: 'position',
    chrom_col: 'chromosome',
    ea_col: 'effect_allele',
    oa_col: 'other_allele',
    snp_col: 'snp',
}, inplace=True)


# Write to file
gwas_df.to_csv(f'{gwas}_{cell_type}_{gene}_gwas.txt', sep='\t', index=False)
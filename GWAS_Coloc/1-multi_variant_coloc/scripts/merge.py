#!/usr/bin/env python3

import pandas as pd 
import sys

# Parse command line arguments
gene = sys.argv[1]
eqtl_data_path = sys.argv[2]
gwas_data_path = sys.argv[3]
cell_type = sys.argv[4]
gwas = sys.argv[5]

# Load eQTL data
eqtl_df = pd.read_csv(eqtl_data_path, sep='\t')

# Calculate varbeta for coloc_abf()
eqtl_df['varbeta'] = eqtl_df['beta_se'] ** 2

def format_snp_id_eqtl(row):
    '''Function that creates a variant ID with alphabetically sorted alleles with format 1:12345:A,C'''
    
    chrom = row['feature_chromosome']
    pos = row['snp_position']
    ea = row['assessed_allele']
    oa = row['other_allele']
    alleles = sorted([ea, oa])
    return f'{chrom}:{pos}:{alleles[0]},{alleles[1]}'

# Create formatted snp column with alleles alphabetically sorted for merging
eqtl_df['snp_formatted_sorted'] = eqtl_df.apply(format_snp_id_eqtl, axis=1)

# Load finemapped GWAS data
gwas_df = pd.read_csv(gwas_data_path, sep='\t')

# Give dataframe columns _eqtl and _gwas to be able to distinguish between column origin after merging
eqtl_df.rename(columns={c: f'{c}_eqtl' for c in eqtl_df.columns if c != 'snp_formatted_sorted'}, inplace=True)
gwas_df.rename(columns={c: f'{c}_gwas' for c in gwas_df.columns if c != 'snp_formatted_sorted' and not c.endswith('_ref')}, inplace=True)

# Merge data
merged = pd.merge(gwas_df, eqtl_df, on='snp_formatted_sorted')

# At this point betas are flipped towards reference panel
# So flip if eQTL tested allele does not match reference panel allele
def flip_beta(row):
    '''Function that flips GWAS betas towards eQTLs if eQTL assessed allele does not match reference panel allele'''

    if row['assessed_allele_eqtl'] != row['a1_ref']:
        return row['beta_gwas'] * -1
    return row['beta_gwas']

# Flip betas towards reference panel
merged['beta_gwas'] = merged.apply(flip_beta, axis=1)
print(merged.head())

# Save to file
merged.to_csv(f'{gwas}_{cell_type}_{gene}_merged.txt', sep='\t', index=False)


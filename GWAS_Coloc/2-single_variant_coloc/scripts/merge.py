#!/usr/bin/env python3

import pandas as pd 
import tabix 
import gzip
import sys

# Parse command line arguments
gene = sys.argv[1]
chromosome = sys.argv[2]
cell_type = sys.argv[3]
eqtl_path = sys.argv[4]
gwas = sys.argv[5]
gwas_path = sys.argv[6]
gwas_eaf_col = sys.argv[7]
annotation_path = sys.argv[8]

# Load annotation data
annotation_df = pd.read_csv(annotation_path, sep='\t', index_col=0)
print(annotation_df.head())

# Define locus start and end postitions
locus_start = annotation_df.loc[gene]['ChrStart'] - 500000
locus_end = annotation_df.loc[gene]['ChrStart'] + 500000

# Load eQTL data
eqtl_df = pd.read_csv(eqtl_path, sep='\t')

# Create eQTL dataframe with reduced window size of 1mb total (500kb upstream and 500kb downstream)
eqtl_df_reduced = eqtl_df[
    (eqtl_df['snp_position'] >= locus_start) &
    (eqtl_df['snp_position'] <= locus_end)
]

eqtl_df_reduced.to_csv(f'{gwas}_{cell_type}_{gene}_reduced.txt', sep='\t', index=False)

# Calculate varbeta for eQTLs
eqtl_df['varbeta'] = eqtl_df['beta_se'] ** 2
print(eqtl_df.head())

def format_snp_id(row):
    '''Function that creates a variant ID with alphabetically sorted alleles with format 1:12345:A,C'''

    chrom = row['feature_chromosome']
    pos = row['snp_position']
    ea = row['assessed_allele']
    oa = row['other_allele']
    alleles = sorted([ea, oa])
    
    return f'{chrom}:{pos}:{alleles[0]},{alleles[1]}'

# Create column with format "chrom:pos:a1,a2" where a1 and a2 are sorted alphabetically for merging purposes
eqtl_df['snp_formatted_sorted'] = eqtl_df.apply(format_snp_id, axis=1)
print(eqtl_df.head())

# Load GWAS data
gwas_df = pd.read_csv(gwas_path, sep='\t')
print(gwas_df.head())

# Add maf column for GWAS if effect allele frequency column is given
if gwas_eaf_col != 'none':
    gwas_df[gwas_eaf_col] = pd.to_numeric(gwas_df[gwas_eaf_col], errors='coerce')
    gwas_df.dropna(subset=[gwas_eaf_col], inplace=True)
    gwas_df['maf'] = gwas_df[gwas_eaf_col]
    gwas_df.loc[gwas_df['maf'] > 0.5, 'maf'] = 1 - gwas_df['maf']
    gwas_df.rename(columns={gwas_eaf_col: 'ea_frequency'}, inplace=True)

# Rename columns with _eqtl and _gwas underscores
eqtl_df.rename(columns={c: f'{c}_eqtl' for c in eqtl_df.columns if c != 'snp_formatted_sorted'}, inplace=True)
gwas_df.rename(columns={c: f'{c}_gwas' for c in gwas_df.columns if c != 'snp_formatted_sorted'}, inplace=True)
print(len(gwas_df.head()))

# Merge GWAS and eQTL data
merged = pd.merge(gwas_df, eqtl_df, on='snp_formatted_sorted')

# Drop variants with duplicate id
merged.drop_duplicates('snp_gwas', keep='first', inplace=True)
merged.drop_duplicates('snp_id_eqtl', keep='first', inplace=True)

# Flip GWAS betas towards eQTLs if the alleles are flipped
merged.loc[merged[f'effect_allele_gwas'] != merged['assessed_allele_eqtl'], 'beta_gwas'] = merged['beta_gwas'] * -1
print(merged.head())

# Write to file
merged.to_csv(f'{gwas}_{cell_type}_{gene}_merged.txt', sep='\t', index=False)

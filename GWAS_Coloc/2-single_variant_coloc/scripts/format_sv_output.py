#!/usr/bin/env python3

import pandas as pd 
import sys

# Parse command line arguments
gene = sys.argv[1]
cell_type = sys.argv[2]
gwas_data_path = sys.argv[3]
eqtl_data_path = sys.argv[4]
coloc_data_path = sys.argv[5]
merged_data_path = sys.argv[6]
sig_eqtl_path = sys.argv[7]
gwas = sys.argv[8]
out_dir = sys.argv[9]
annotation_path = sys.argv[10]

# Load GWAS data
gwas_df = pd.read_csv(gwas_data_path, sep='\t')
num_gwas_variants = len(gwas_df)

# Load eQTL data
eqtl_df = pd.read_csv(eqtl_data_path, sep='\t')
num_eqtl_variants = len(eqtl_df)

# Load coloc results
coloc_df = pd.read_csv(coloc_data_path, sep='\t')
print(coloc_df.head())

# Load merged eQTL/GWAS data
merged_df = pd.read_csv(merged_data_path, sep='\t')
print(merged_df.head())
num_overlapping_variants = len(merged_df)

# Add number of GWAS/eQTL/overlapping variants to output
coloc_df['Number of GWAS variants'] = num_gwas_variants
coloc_df['Number of eQTL variants'] = num_eqtl_variants
coloc_df['Number of overlapping variants'] = num_overlapping_variants

# Get GWAS info
gwas_top_row = merged_df.nsmallest(1, 'p_value_gwas').iloc[0]
coloc_df['GWAS variant'] = gwas_top_row['snp_gwas']
coloc_df['GWAS variant chromosome'] = gwas_top_row['chromosome_gwas']
coloc_df['GWAS variant position'] = gwas_top_row['position_gwas']
coloc_df['GWAS variant p-value'] = gwas_top_row['p_value_gwas']
coloc_df['GWAS variant beta'] = gwas_top_row['beta_gwas']
coloc_df['GWAS variant standard error'] = gwas_top_row['standard_error_gwas']

# Get eQTL info
eqtl_top_row = merged_df.nsmallest(1, 'p_value_eqtl').iloc[0]
coloc_df['eQTL variant'] = eqtl_top_row['snp_id_eqtl']
coloc_df['eQTL variant chromosome'] = eqtl_top_row['feature_chromosome_eqtl']
coloc_df['eQTL variant position'] = eqtl_top_row['snp_position_eqtl']
coloc_df['eQTL variant p-value'] = eqtl_top_row['p_value_eqtl']
coloc_df['eQTL variant beta'] = eqtl_top_row['beta_eqtl']
coloc_df['eQTL variant standard error'] = eqtl_top_row['beta_se_eqtl']

# Calculate absolute difference between top GWAS and top eQTL variant
coloc_df['Distance between variants'] = abs(coloc_df['GWAS variant position'] - coloc_df['eQTL variant position'])

# Load annotation data
annotation_df = pd.read_csv(annotation_path, sep='\t')
print(annotation_df.head())

def get_nearest_gene(row):
    '''Function that retrieves the gene closest to a genomic position (chrom:pos) from the annotation'''

    gene = row['gene']
    chromosome = row['GWAS variant chromosome']
    position = row['GWAS variant position']

    # Filter annotation on chromosome
    annotation_filtered = annotation_df[annotation_df['Chr'] == chromosome]

    # For each gene on chromosome, calculate absolute distance between gene TSS and input position
    annotation_filtered['distance_to_variant'] = abs(annotation_filtered['ChrStart'] - position)

    nearest_gene_row = annotation_filtered.nsmallest(1, 'distance_to_variant').iloc[0]
    nearest_gene = nearest_gene_row['Gene']
    nearest_gene_distance = nearest_gene_row['distance_to_variant']

    return nearest_gene, nearest_gene_distance

coloc_df[['Nearest gene to GWAS variant', 'Distance between GWAS variant and nearest gene']] = coloc_df.apply(get_nearest_gene, axis=1, result_type='expand')
coloc_df['Nearest gene is coloc gene'] = coloc_df['gene'] == coloc_df['Nearest gene to GWAS variant']

# Indicate whether a top GWAS variant is significant
gwas_sig_threshold = 0.00000005
coloc_df['GWAS variant is significant'] = coloc_df['GWAS variant p-value'] < gwas_sig_threshold

# Load eQTL top effects file
sig_eqtl_df = pd.read_csv(sig_eqtl_path, sep='\t')
sig_eqtl_df = sig_eqtl_df[sig_eqtl_df['qval'] < 0.05]
sig_eqtl_df.set_index('feature_id', inplace=True)

def check_eqtl_significance(row):
    '''Function that checks if an eQTL variant is significant'''

    # Get the nominal p-value threshold from the top effects file for a gene
    gene = row['gene']
    eqtl_p = row['eQTL variant p-value']
    pval_threshold = sig_eqtl_df.loc[gene]['PvalueNominalThreshold']

    return eqtl_p < pval_threshold

# Indicate whether a top eQTL variant is significant
coloc_df['eQTL variant is significant'] = coloc_df.apply(check_eqtl_significance, axis=1)

coloc_df['Analysis type'] = 'coloc.abf()'
coloc_df['GWAS lambda'] = '-'
coloc_df['eQTL lambda'] = '-'
coloc_df['GWAS credible set'] = '-'
coloc_df['eQTL credible set'] = '-'
coloc_df['GWAS'] = gwas
print(coloc_df.head())

# Rename columns
coloc_df.rename(columns={
    'gene': 'Gene',
    'cell_type': 'Cell type',
}, inplace=True)

# Filter and reorder relevant columns
coloc_df = coloc_df[[
    'Gene',
    'Number of GWAS variants',
    'Number of eQTL variants',
    'Number of overlapping variants',
    'GWAS variant',
    'GWAS variant chromosome',
    'GWAS variant position',
    'GWAS variant p-value',
    'GWAS variant beta',
    'GWAS variant standard error',
    'eQTL variant',
    'eQTL variant chromosome',
    'eQTL variant position',
    'eQTL variant p-value',
    'eQTL variant beta',
    'eQTL variant standard error',
    'Distance between variants',
    'PPH0',
    'PPH1',
    'PPH2',
    'PPH3',
    'PPH4',
    'GWAS credible set',
    'eQTL credible set',
    'GWAS lambda',
    'eQTL lambda',
    'Cell type',
    'Nearest gene to GWAS variant',
    'Distance between GWAS variant and nearest gene',
    'Nearest gene is coloc gene',
    'GWAS variant is significant',
    'eQTL variant is significant',
    'Analysis type',
    'GWAS',
]]

# Save to file
coloc_df.to_csv(f'{out_dir}/{gene}_sv_coloc_out_formatted.txt', sep='\t', index=False)
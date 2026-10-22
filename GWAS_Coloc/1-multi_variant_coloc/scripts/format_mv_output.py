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
gwas_summary_path = sys.argv[9]
annotation_path = sys.argv[10]

# Load GWAS summary file with credible set count per gene
gwas_summary_df = pd.read_csv(gwas_summary_path, sep='\t')
gwas_summary_df['gwas_gene'] = gwas_summary_df['gwas'] + '_' + gwas_summary_df['gene']
gwas_summary_df.set_index('gwas_gene', inplace=True)
num_gwas_variants = gwas_summary_df.loc[f'{gwas}_{gene}']['num_variants']
num_gwas_cs = gwas_summary_df.loc[f'{gwas}_{gene}']['num_cs']

# Load eQTL data
eqtl_df = pd.read_csv(eqtl_data_path, sep='\t')
num_eqtl_variants = len(eqtl_df)

# Get number of eQTL credible sets
eqtl_cs_df = eqtl_df[~eqtl_df['SusieRss_CS'].isna()]
if len(eqtl_cs_df) == 0:
    num_eqtl_cs = 0
else:
    num_eqtl_cs = eqtl_cs_df['SusieRss_CS'].nunique()

# Load coloc output
coloc_df = pd.read_csv(coloc_data_path, sep='\t')
print(coloc_df.head())

# Load the merged eQTL/GWAS data
merged_df = pd.read_csv(merged_data_path, sep='\t')
print(merged_df.head())

# Count variants in data (number of overlapping variants)
num_overlapping_variants_filtered = len(merged_df)

# Get GWAS info
merged_df_gwas = merged_df[['snp_formatted_sorted', 'snp_gwas', 'chromosome_gwas', 'position_gwas', 'p_value_gwas', 'beta_gwas', 'standard_error_gwas']]

# Merge GWAS info with coloc output
coloc_df = pd.merge(coloc_df, merged_df_gwas, left_on='sorted_snp_gwas', right_on='snp_formatted_sorted')

# Get eQTL info
merged_df_eqtl = merged_df[['snp_formatted_sorted', 'snp_id_eqtl', 'feature_chromosome_eqtl', 'snp_position_eqtl', 'p_value_eqtl', 'beta_eqtl', 'beta_se_eqtl']]

# Merge eQTL info with coloc output
coloc_df = pd.merge(coloc_df, merged_df_eqtl, left_on='sorted_snp_eqtl', right_on='snp_formatted_sorted')

# Calculate distance between top eQTL and top GWAS variant
coloc_df['var_distance'] = abs(coloc_df['position_gwas'] - coloc_df['snp_position_eqtl'])

# Add number of GWAS/eQTL/overlapping variants
coloc_df['num_gwas_variants'] = num_gwas_variants
coloc_df['num_eqtl_variants'] = num_eqtl_variants
coloc_df['num_overlapping_variants'] = num_overlapping_variants_filtered
coloc_df['gene'] = gene

# Set credible set number to 0 if not credible set is found for the GWAS and eQTL
if num_gwas_cs == 0:
    coloc_df['idx1'] = 0

if num_eqtl_cs == 0:
    coloc_df['idx2'] = 0

# Rename columns
coloc_df.rename(columns={
    'nsnps': 'Number of overlapping variants',
    'PP.H0.abf': 'PPH0',
    'PP.H1.abf': 'PPH1',
    'PP.H2.abf': 'PPH2',
    'PP.H3.abf': 'PPH3',
    'PP.H4.abf': 'PPH4',
    'idx1': 'GWAS credible set',
    'idx2': 'eQTL credible set',
    'gwas_lambda': 'GWAS lambda',
    'eqtl_lambda': 'eQTL lambda',
    'cell_type': 'Cell type',
    'snp_gwas': 'GWAS variant',
    'chromosome_gwas': 'GWAS variant chromosome',
    'position_gwas': 'GWAS variant position',
    'p_value_gwas': 'GWAS variant p-value',
    'beta_gwas': 'GWAS variant beta',
    'standard_error_gwas': 'GWAS variant standard error',
    'snp_id_eqtl': 'eQTL variant',
    'feature_chromosome_eqtl': 'eQTL variant chromosome',
    'snp_position_eqtl': 'eQTL variant position',
    'p_value_eqtl': 'eQTL variant p-value',
    'beta_eqtl': 'eQTL variant beta',
    'beta_se_eqtl': 'eQTL variant standard error',
    'var_distance': 'Distance between variants',
    'num_gwas_variants': 'Number of GWAS variants',
    'num_eqtl_variants': 'Number of eQTL variants',
    'gene': 'Gene',
}, inplace=True)

# Load annotation data
annotation_df = pd.read_csv(annotation_path, sep='\t')
print(annotation_df.head())

def get_nearest_gene(row):
    '''Function that retrieves the gene closest to a genomic position (chrom:pos) from the annotation'''

    gene = row['Gene']
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
coloc_df['Nearest gene is coloc gene'] = coloc_df['Gene'] == coloc_df['Nearest gene to GWAS variant']

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
    gene = row['Gene']
    eqtl_p = row['eQTL variant p-value']
    pval_threshold = sig_eqtl_df.loc[gene]['PvalueNominalThreshold']

    return eqtl_p < pval_threshold

# Indicate whether a top eQTL variant is significant
coloc_df['eQTL variant is significant'] = coloc_df.apply(check_eqtl_significance, axis=1)

coloc_df['Analysis type'] = 'SuSiE + coloc.bf_bf()'

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
]]

# Save to file
coloc_df.to_csv(f'{gwas}_{cell_type}_{gene}_mv_coloc_out_formatted.txt', sep='\t', index=False)
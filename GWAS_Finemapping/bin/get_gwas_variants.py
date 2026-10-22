#!/usr/bin/env python3

import pandas as pd 
import tabix
import gzip
import sys

# Parse command line arguments
gene = sys.argv[1]
pval_col = sys.argv[2]
snp_col = sys.argv[3]
beta_col = sys.argv[4]
se_col = sys.argv[5]
chrom_col = sys.argv[6]
ea_col = sys.argv[7]
oa_col = sys.argv[8]
pos_col = sys.argv[9]
eqtl_variant_file = sys.argv[10]
gwas_path = sys.argv[11]
ref_snp_list_path = sys.argv[12]
eaf_col = sys.argv[13]
chromosome = sys.argv[14]

# Load variants tested in eQTL analyses
eqtl_variant_df = pd.read_csv(eqtl_variant_file, sep='\t', names=['snp_formatted_sorted', 'position'])
print(eqtl_variant_df.head())

# Extract locus start and end from the dataframe
locus_start = eqtl_variant_df['position'].min()
locus_end = eqtl_variant_df['position'].max()

tb = tabix.open(gwas_path)

# Retrieve variants from sum stats
records = tb.query(str(chromosome), locus_start, locus_end)
gwas_df = pd.DataFrame(list(records))

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
    '''Function that creates a variant ID with alphabetically sorted alleles for GWAS data'''
    chrom = row[chrom_col]
    pos = row[pos_col]
    ea = row[ea_col]
    oa = row[oa_col]
    alleles = sorted([ea, oa])
    return f'{chrom}:{pos}:{alleles[0]},{alleles[1]}'

gwas_df['snp_formatted_sorted'] = gwas_df.apply(format_snp_id_gwas, axis=1)

# Filter variants that are tested in eQTLs
gwas_df = gwas_df[gwas_df['snp_formatted_sorted'].isin(eqtl_variant_df['snp_formatted_sorted'].tolist())]
print(gwas_df.head())
print(gwas_df.shape)


def format_snp_id_ref(row):
    '''Function that creates a variant ID with alphabetically sorted alleles for reference panel'''
    snp_id = row['ref_variant']
    chrom, pos, alleles = snp_id.split(':')
    a1, a2 = alleles.split(',')
    alleles = sorted([a1, a2])
    return f'{chrom}:{pos}:{alleles[0]},{alleles[1]}', a1, a2

# Load list of reference variants
ref_snp_df = pd.read_csv(ref_snp_list_path, sep='\t', names=['ref_variant'])

# Create a sorted variant ID for reference snps, along with the alleles
ref_snp_df[['snp_formatted_sorted', 'a1_ref', 'a2_ref']] = ref_snp_df.apply(format_snp_id_ref, axis=1, result_type='expand')
print(ref_snp_df.head())

# Merge GWAS data with reference data
gwas_df = pd.merge(gwas_df, ref_snp_df, on='snp_formatted_sorted')

# Drop duplicate SNP IDs if any exist
gwas_df.drop_duplicates('snp_formatted_sorted', inplace=True)

def flip_beta(row):
    '''Function that flips GWAS beta if effect allele does not match reference effect allele'''

    beta = float(row[beta_col])
    if row[ea_col] != row['a1_ref']:
        return beta * -1 
    else:
        return beta

gwas_df[beta_col] = gwas_df.apply(flip_beta, axis=1)
gwas_df.sort_values(pos_col, inplace=True)

# Rename columns
gwas_df.rename(columns={
    pval_col: 'p_value',
    snp_col: 'snp',
    beta_col: 'beta',
    se_col: 'standard_error',
    chrom_col: 'chromosome',
    ea_col: 'effect_allele',
    oa_col: 'other_allele',
    pos_col: 'position',
}, inplace=True)

def get_maf(eaf):
    '''Function that calculates MAF from AF'''

    if eaf > 0.5:
        return 1 - eaf
    return eaf

# If effect allele frequency column is not none, rename, clean and calculate MAF
if eaf_col != 'none':
    gwas_df.rename(columns={eaf_col: 'effect_allele_freq'}, inplace=True)
    gwas_df['effect_allele_freq'] = pd.to_numeric(gwas_df['effect_allele_freq'], errors='coerce')
    gwas_df.dropna(subset='effect_allele_freq', inplace=True)
    gwas_df['maf'] = gwas_df['effect_allele_freq'].apply(get_maf)

# Calculate varbeta
gwas_df['varbeta'] = gwas_df['standard_error'].astype(float) ** 2

# Save to file
gwas_df.to_csv(f'{gene}_gwas.txt', sep='\t', index=False)
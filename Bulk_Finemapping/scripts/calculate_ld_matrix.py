#!/usr/bin/env python3

import tabix
import numpy as np 
import pandas as pd 
import sys, gzip

gene = sys.argv[1]
chromosome = sys.argv[2]
eqtl_path = sys.argv[3]
cell_type = sys.argv[4]
vcf_path = sys.argv[5]
gte_path = sys.argv[6]

# Load GTE file to get a list of included samples
gte_df = pd.read_csv(gte_path, sep='\t', names=['DNA', 'RNA', 'dataset'])
included_samples = gte_df['DNA'].tolist()

# Load eQTL data
eqtl_df = pd.read_csv(eqtl_path, sep='\t')
print(eqtl_df.head())

# Extract start and end position for locus from eQTL data
start = eqtl_df['SNPPos'].min()
end = eqtl_df['SNPPos'].max()

def format_variant_id(row):
    rsid = row['SNP']
    alleles = row['SNPAlleles'].split('/')
    alt_allele = row['SNPEffectAllele']
    ref_allele = [a for a  in alleles if a != alt_allele][0]
    return f'{rsid}_{ref_allele}_{alt_allele}'

eqtl_df['variant_id_formatted'] = eqtl_df.apply(format_variant_id, axis=1)
eqtl_variants = eqtl_df['variant_id_formatted'].tolist()
print(eqtl_df.head())

# Extract variants in locus from VCF
tb = tabix.open(vcf_path)
records = tb.query(str(chromosome), start - 1, end + 1)

def count_genotypes(gt_string):
    '''Function that converts a genotype count string to a numerical value'''

    if gt_string.isnumeric():
        return gt_string
    
    if gt_string in ['0|0', '0/0']:
        return 0
    if gt_string in ['0|1', '1|0', '1/0', '0/1']:
        return 1
    if gt_string in ['1|1', '1/1']:
        return 2
    
    return np.nan

genotype_dict = {}

line_ctr = 0

vcf_variants = []

# Get a list of included sample indices
with gzip.open(vcf_path, 'rt') as vcf_file:
    for line in vcf_file:
        if line.startswith('#CHROM'):
            samples = [s.strip() for s in line.split('\t')[9:]]
            included_sample_indices = [i for i in range(len(samples)) if samples[i] in included_samples]
            break

# Loop through variants from VCF
for elems in records:

    # Get indices of format field and GT column
    format_field = elems[8]
    gt_col_index = format_field.split(':').index('GT')

    # It could be that not all records have a DS field
    ds_col_index = -1
    if 'DS' in format_field:
        ds_col_index = format_field.split(':').index('DS')

    gt_col_index = format_field.split(':').index('GT')

    # Create a formatted id field to match
    id_formatted = f'{elems[2]}_{elems[3]}_{elems[4]}'
    
    # Skip if variant is not in eQTL data
    if not id_formatted in eqtl_variants:
        continue

    vcf_variants.append(id_formatted)

    # Get genotype field
    genotypes = elems[9:]

    # Filter samples
    genotypes = [genotypes[i] for i in included_sample_indices]

    # Get dosages and genotype counts
    ds_values = [g.split(':')[ds_col_index] for g in genotypes]
    gt_counts = [count_genotypes(gt.split(':')[gt_col_index]) for gt in genotypes]

    # Replace missing dosage values by nan; otherwise the correlation function crashes
    ds_values = [ds_values[i] if ds_values[i] != '.' else np.nan for i in range(len(ds_values))]

    # If no dosages are available, use hard genotype counts
    if ds_col_index != -1:
        genotype_dict[snp_id] = ds_values
    else:
        genotype_dict[snp_id] = gt_counts


    line_ctr += 1

    if line_ctr % 100000 == 0:
        print(f'{line_ctr} variants processed')

    if len(genotype_dict) == len(eqtl_variants):
        print('All variants retrieved')
        break # Break if all variants are retrieved

res_df = pd.DataFrame(genotype_dict)

print(res_df.shape) 
print(res_df.head())

print('Correlating...')
corr_df = res_df.corr()
print(corr_df)

# Write to file
corr_df.to_csv(f'{gene}_ld.txt', sep='\t')
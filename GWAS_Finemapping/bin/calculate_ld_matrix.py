#!/usr/bin/env python3

import tabix
import numpy as np 
import pandas as pd 
import sys, gzip
import pyreadr

locus = sys.argv[1]
locus_file = sys.argv[2]
vcf_path = sys.argv[3]

# load sum stats
locus_df = pd.read_csv(locus_file, sep='\t')
print(locus_df.head())

# Get gene chromosome from first row
chromosome = locus_df.iloc[0]['chromosome']

# Extract start and end position for locus from eQTL data
start = locus_df['position'].min()
end = locus_df['position'].max()

locus_variants = locus_df['ref_variant'].tolist()

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

def flip_variant_id(variant_id):
    '''Function that flips alleles in a variant identifier'''

    chrom, pos, alleles = variant_id.split(':')
    a1, a2 = alleles.split(',')
    return f'{chrom}:{pos}:{a2},{a1}'

genotype_dict = {}

line_ctr = 0

vcf_variants = []

# Loop through variants from VCF
for elems in records:
    # exit()

    # Get indices of format field and GT column
    format_field = elems[8]
    gt_col_index = format_field.split(':').index('GT')

    # Create a formatted id field to match
    variant_id = elems[2]
    
    # Skip if variant is not in eQTL data
    if not variant_id in locus_variants:
        continue

    vcf_variants.append(variant_id)

    # Get genotype field
    genotypes = elems[9:]

    # Get genotype counts
    gt_counts = [count_genotypes(gt.split(':')[gt_col_index]) for gt in genotypes]

    # Add genotype counts to output dict
    genotype_dict[variant_id] = gt_counts

    line_ctr += 1

    if line_ctr % 100000 == 0:
        print(f'{line_ctr} variants processed')

    if len(genotype_dict) == len(locus_variants):
        print('All variants retrieved')
        break # Break if all variants are retrieved

res_df = pd.DataFrame(genotype_dict)


if len(locus_variants) != len(res_df):
    print([v for v in locus_variants if v not in genotype_dict.keys()])
    # exit()

print(res_df.shape) 
print(res_df.head())

print('Correlating...')
corr_mat = np.corrcoef(res_df.values, rowvar=False)
corr_df = pd.DataFrame(corr_mat, index=res_df.columns, columns=res_df.columns)
print(corr_df)

# Save a numpy file to reduce storage
np.save(f'{locus}_ld.npy', corr_df)

#!/usr/bin/env python3 

import pandas as pd 
import gzip
import tabix
import sys

gene = sys.argv[1]
chromosome = sys.argv[2]
variant_list_path = sys.argv[3]
eqtl_vcf_path = sys.argv[4]

variant_df = pd.read_csv(variant_list_path, sep='\t')
variant_df = variant_df[variant_df['gene'] == gene]
print(variant_df.head())

def separate_ids(row):
    '''Function that splits all properties in a SNP ID like rs1234567:1:23512422:T into separate fields'''

    variant = row['variant']
    elems = variant.split(':')
    ea = elems[-1]
    pos = elems[-2]
    chrom = elems[-3]
    snp_id = ':'.join(elems[0:-3])
    return chrom, int(pos), ea, snp_id

variant_df[['chrom', 'pos', 'ea', 'snp_id']] = variant_df.apply(separate_ids, axis=1, result_type='expand')
variant_df.sort_values('pos', inplace=True)
print(variant_df.head())

# Create list of eQTL variants
eqtl_variants = variant_df['snp_id'].tolist()

# Extract locus start and end
locus_start = variant_df['pos'].min()
locus_end = variant_df['pos'].max()

# Extract variants in locus from VCF
tb = tabix.open(eqtl_vcf_path)
records = tb.query(str(chromosome), locus_start - 1, locus_end + 1)

vcf_variants = []

# Create a list of variants with sorted alleles
for elems in records:
    vcf_snp_id = elems[2]

    if vcf_snp_id not in eqtl_variants:
        continue

    vcf_chrom = elems[0]
    vcf_pos = elems[1]
    vcf_ref = elems[3]
    vcf_alt = elems[4]
    vcf_alleles = sorted([vcf_ref, vcf_alt])

    sorted_variant_id = f'{vcf_chrom}:{vcf_pos}:{vcf_alleles[0]},{vcf_alleles[1]}'
    vcf_variants.append(f'{sorted_variant_id}\t{vcf_pos}')

# Write to file
with open(f'{gene}_eqtl_variants_alleles_alpha_sorted.txt', 'w') as out_file:
    out_file.write('\n'.join(vcf_variants))
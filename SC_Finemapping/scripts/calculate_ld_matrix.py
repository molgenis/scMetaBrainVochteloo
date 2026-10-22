import tabix
import numpy as np 
import pandas as pd 
import sys, gzip

gene = sys.argv[1]
chromosome = sys.argv[2]
eqtl_path = sys.argv[3]
cell_type = sys.argv[4]
gte_path = sys.argv[5]
vcf_path = sys.argv[6]

# Get a list of included samples
gte_df = pd.read_csv(gte_path, sep='\t', names=['DNA', 'RNA', 'dataset'])
included_samples = gte_df['DNA'].tolist()

# Load eQTL data
eqtl_df = pd.read_csv(eqtl_path, sep='\t')
print(eqtl_df.head())

# Extract start and end position for locus from eQTL data
start = int(eqtl_df['snp_position'].min())
end = int(eqtl_df['snp_position'].max())

eqtl_variants = eqtl_df['snp_id'].tolist()
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
both_variants = []

# Get a list of included sample indices from the VCF
with gzip.open(url, 'rt') as vcf_file:
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

    # Get snp id
    snp_id = elems[2]
    
    # Skip if variant is not in eQTL data
    if not snp_id in eqtl_variants:
        continue

    vcf_variants.append(snp_id)
    alleles = f'{elems[3]}/{elems[4]}'
    both_variants.append(alleles)

    # Get genotype fields
    genotypes = elems[9:]

    # Get dosages
    ds_values = [g.split(':')[ds_col_index] for g in genotypes]
    ds_values = [ds_values[i] for i in included_sample_indices]

    # Get genotype counts
    gt_counts = [count_genotypes(g.split(':')[gt_col_index]) for g in genotypes]
    gt_counts = [gt_counts[i] for i in included_sample_indices]
    
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

# Convert genotype dict to dataframe
res_df = pd.DataFrame(genotype_dict)

print('Correlating...')
corr_df = res_df.corr()
print(corr_df)

eqtl_df = eqtl_df[eqtl_df['snp_id'].isin(vcf_variants)]

# Add other allele column to output
allele_df = pd.DataFrame({'snp_id': vcf_variants, 'alleles': both_variants})
eqtl_df = pd.merge(eqtl_df, allele_df, on='snp_id')
eqtl_df.to_csv(f'{gene}_eqtls.txt', sep='\t', index=False)

# Save LD matrix to file
corr_df.to_csv(f'{gene}_ld.txt', sep='\t')
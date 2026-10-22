import pandas as pd 
import sys

ds1 = sys.argv[1]
ds1_path = sys.argv[2]
ds2 = sys.argv[3]
ds2_path = sys.argv[4]
gene = sys.argv[5]
annotation_path = sys.argv[6]

df1 = pd.read_csv(ds1_path, sep='\t')
df2 = pd.read_csv(ds2_path, sep='\t')

# Give describing prefixes to LBF columns
df1.rename(columns={c: f'{c}_ds1' for c in df1.columns if c.startswith('lbf')}, inplace=True)
df2.rename(columns={c: f'{c}_ds2' for c in df2.columns if c.startswith('lbf')}, inplace=True)

# Get locus start and end position from annotation
annotation_df = pd.read_csv(annotation_path, sep='\t')
tss_dict = annotation_df.set_index('Gene')['ChrStart'].to_dict()
lower_bound = tss_dict[gene] - 1000000
upper_bound = tss_dict[gene] + 1000000

def format_bulk_col(row):
    '''Function that creates a variant ID with alleles alphabetically sorted like 1:12345:A,C for bulk rows'''
    chrom = row['GeneChr']
    pos = row['SNPPos']
    alleles = sorted(row['SNPAlleles'].split('/'))
    return f'{chrom}:{pos}:{alleles[0]},{alleles[1]}'

def format_sc_col(row):
    '''Function that creates a variant ID with alleles alphabetically sorted like 1:12345:A,C for single cell rows'''
    chrom = row['feature_chromosome']
    pos = row['snp_position']
    ea = row['assessed_allele']
    oa = row['other_allele']
    alleles = sorted([ea, oa])
    return f'{chrom}:{pos}:{alleles[0]},{alleles[1]}'


def get_maf(af):
    if af > 0.5:
        return 1 - af
    return af

bulk_datasets = ['basalganglia-EUR', 'cerebellum-EUR', 'cortex-EUR', 'hippocampus-EUR', 'spinalcord-EUR']

# Format data for the first dataset
if ds1 in bulk_datasets:
    df1['snp_formatted'] = df1.apply(format_bulk_col, axis=1)
    df1['maf'] = df1['SNPEffectAlleleFreq'].apply(get_maf)
    df1['varbeta'] = df1['MetaSE'] ** 2
    df1.rename(columns={'MetaBeta': 'beta', 'SNPChr': 'chromosome', 'MetaSE': 'standard_error', 'SNPPos': 'position', 'MetaP': 'p_value', 'SNP': 'snp'}, inplace=True)
else:
    df1['snp_formatted'] = df1.apply(format_sc_col, axis=1)
    df1['varbeta'] = df1['beta_se'] ** 2
    df1.rename(columns={'snp_position': 'position', 'beta_se': 'standard_error', 'feature_chromosome': 'chromosome', 'snp_id': 'snp'}, inplace=True)


# Format data for the second dataset
if ds2 in bulk_datasets:
    df2['snp_formatted'] = df2.apply(format_bulk_col, axis=1)
    df2['maf'] = df2['SNPEffectAlleleFreq'].apply(get_maf)
    df2['varbeta'] = df2['MetaSE'] ** 2
    df2.rename(columns={'MetaBeta': 'beta', 'SNPChr': 'chromosome', 'MetaSE': 'standard_error', 'SNPPos': 'position', 'MetaP': 'p_value', 'SNP': 'snp'}, inplace=True)
else:
    df2['snp_formatted'] = df2.apply(format_sc_col, axis=1)
    df2['varbeta'] = df2['beta_se'] ** 2
    df2.rename(columns={'snp_position': 'position', 'beta_se': 'standard_error', 'feature_chromosome': 'chromosome', 'snp_id': 'snp'}, inplace=True)

merged = pd.merge(df1, df2, on='snp_formatted', suffixes=('_ds1', '_ds2'))

# Save to file
merged.to_csv(f'{gene}_{ds1}_{ds2}.txt', sep='\t', index=False)

# Reduce window to 1mb for single variant coloc
df1 = df1[(df1['position'] >= lower_bound) & (df1['position'] <= upper_bound)]
df2 = df2[(df2['position'] >= lower_bound) & (df2['position'] <= upper_bound)]
merged = merged[(merged['position_ds1'] >= lower_bound) & (merged['position_ds1'] <= upper_bound)]

# Write reduced dataframes to file
df1.to_csv(f'{gene}_{ds1}_reduced_window.txt', sep='\t', index=False)
df2.to_csv(f'{gene}_{ds2}_reduced_window.txt', sep='\t', index=False)
merged.to_csv(f'{gene}_{ds1}_{ds2}_reduced_window.txt', sep='\t', index=False)

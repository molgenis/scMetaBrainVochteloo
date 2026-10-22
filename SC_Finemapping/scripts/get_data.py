import pandas as pd 
import tabix
import gzip, sys

gene = sys.argv[1]
chromosome = sys.argv[2]
eqtl_path = sys.argv[3]
sig_eqtl_path = sys.argv[4]
sample_size = int(sys.argv[5])

# Extract cell type from path
cell_type = eqtl_path.split('/')[-2]

sig_eqtl_df = pd.read_csv(sig_eqtl_path, sep='\t')
sig_eqtl_df.set_index('feature_id', inplace=True)

# Extract start and end coordinates for gene
feature_start = int(sig_eqtl_df.loc[gene]['feature_start'])
feature_end = int(sig_eqtl_df.loc[gene]['feature_end'])

tb = tabix.open(eqtl_path)

# Retrieve variants from sum stats
records = tb.query(str(chromosome), feature_start - 1000001, feature_end + 1000001)
res_df = pd.DataFrame(list(records))
print(res_df.head())

# Get columns names
with gzip.open(eqtl_path, 'rt') as eqtl_file:
    for line in eqtl_file:
        line = line.strip()
        columns = line.split('\t')
        break

res_df.columns = columns
res_df = res_df[res_df['feature_id'] == gene]
print(res_df.head())

# Filter HWE P
res_df['hwe_p'] = res_df['hwe_p'].astype(float)
res_df = res_df[res_df['hwe_p'] > 1e-4]

# Add z-score column
res_df['beta'] = res_df['beta'].astype(float)
res_df['beta_se'] = res_df['beta_se'].astype(float)
res_df['z_score'] = res_df['beta'] / res_df['beta_se']
print(res_df.head())

# Save file
res_df.to_csv(f'{gene}_eqtls.txt', sep='\t', index=False)

with open(f'{gene}_snps.txt', 'w') as snp_file:
    snp_file.write('\n'.join(res_df['snp_id'].to_list()))
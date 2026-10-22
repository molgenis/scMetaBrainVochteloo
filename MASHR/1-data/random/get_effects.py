import sys
import pandas as pd 

cell_type = sys.argv[1]
chromosome = int(sys.argv[2])
base_cell_type = sys.argv[3]
random_effects_dir = sys.argv[4]
eqtl_base_dir = sys.argv[5]
out_dir = sys.argv[6]

# Load effects from base cell type
base_path = f'{random_effects_dir}/{base_cell_type}_chr{chromosome}_betas.txt'
base_df = pd.read_csv(base_path, sep='\t')
base_eqtls = base_df['snp_gene'].tolist()

path = f'{eqtl_base_dir}/{cell_type}/chr{chromosome}_qtl_results_all_sortedbysnp.bgz.txt.gz'

dfs = []

# Loop through eQTL sum stats and extract base eQTLs
for df in pd.read_csv(path, sep='\t', chunksize=10000):

    df['snp_gene'] = df['snp_id'] + '_' + df['feature_id']
    df = df[df['snp_gene'].isin(base_eqtls)]

    if len(df) > 0:
        dfs.append(df)
        print(len(dfs))

# Concat final dataframe
res_df = pd.concat(dfs)

# Save betas to file
res_beta_df = res_df[['snp_gene', 'beta']].rename(columns={'beta': cell_type})
res_beta_df.to_csv(f'{out_dir}/{cell_type}_chr{chromosome}_betas.txt', sep='\t', index=False)

# Save standard errors to file
res_se_df = res_df[['snp_gene', 'beta_se']].rename(columns={'beta_se': cell_type})
res_se_df.to_csv(f'{out_dir}/{cell_type}_chr{chromosome}_standard_errors.txt', sep='\t', index=False)
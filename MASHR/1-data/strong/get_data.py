import sys
import pandas as pd 
import os

cell_type = sys.argv[1]
chromosome = int(sys.argv[2])
eqtl_base_dir = sys.argv[3]
out_dir = sys.argv[4]

# Load top effects file and extract list of top effects for specified chromosome
top_effects_df = pd.read_csv('top_effects_all_cell_types.txt', sep='\t')
top_effects_df = top_effects_df[top_effects_df['feature_chromosome'] == chromosome]
top_effects_df.drop_duplicates('snp_gene', keep='first', inplace=True)
top_effect_eqtls = top_effects_df['snp_gene'].tolist()
print(top_effect_eqtls)

eqtl_path = f'{eqtl_base_dir}/{cell_type}/chr{chromosome}_qtl_results_all_sortedbysnp.bgz.txt.gz'

dfs = []

# Loop through eQTL sum stats and extract top effects
for df in pd.read_csv(eqtl_path, sep='\t', chunksize=10000):

    df['snp_gene'] = df['snp_id'] + '_' + df['feature_id']
    df = df[df['snp_gene'].isin(top_effect_eqtls)]

    if len(df) > 0:
        dfs.append(df)
        print(len(dfs))

# Merge dataframes
res_df = pd.concat(dfs)
print(res_df.head())

# Filter necessary columns and write to output file 
res_df = res_df[['snp_gene', 'beta', 'beta_se']]
res_df.to_csv(f'{out_dir}/{cell_type}_chr{chromosome}.txt', sep='\t', index=False)

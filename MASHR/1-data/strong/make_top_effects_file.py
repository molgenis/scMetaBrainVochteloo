import pandas as pd 
import sys 

eqtl_base_dir = sys.argv[1]

cell_types = ['AST', 'END', 'EX', 'IN', 'MIC', 'OLI', 'OPC', 'PER']

dfs = []

# For every cell type, load significant top effects
for cell_type in cell_types:

    top_eqtl_path = f'{eqtl_base_dir}/{cell_type}/top_qtl_results_all_MultiTest.txt'
    top_df = pd.read_csv(top_eqtl_path, sep='\t')
    top_df = top_df[top_df['qval'] < 0.05]

    top_df['snp_gene'] = top_df['snp_id'] + '_' + top_df['feature_id']
    top_df['cell_type'] = cell_type

    top_df = top_df[['snp_gene', 'feature_chromosome', 'cell_type']]
    print(top_df.head())

    dfs.append(top_df)

# Concat top effects for all cell types
res_df = pd.concat(dfs)
print(res_df.head())

# Save to file
res_df.to_csv('top_effects_all_cell_types.txt', sep='\t', index=False)

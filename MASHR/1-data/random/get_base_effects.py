import sys, os
import pandas as pd 

cell_type = sys.argv[1]
chromosome = int(sys.argv[2])
eqtl_base_dir = sys.argv[3]
out_dir = sys.argv[4]

# Get list of unique eGenes for base cell type
top_effects_path = '../strong/top_effects_all_cell_types.txt'
top_effect_df = pd.read_csv(top_effects_path, sep='\t')
top_effect_df = top_effect_df[top_effect_df['feature_chromosome'] == chromosome]
top_effects = top_effect_df['snp_gene'].unique().tolist()
genes = [e.split('_')[1] for e in top_effects]

eqtl_path = f'{eqtl_base_dir}/{cell_type}/chr{chromosome}_qtl_results_all_sortedbysnp.bgz.txt.gz'

# Load, filter and format eQTL sum stats
df = pd.read_csv(eqtl_path, sep='\t')
df['snp_gene'] = df['snp_id'] + '_' + df['feature_id']
df = df[~df['snp_gene'].isin(top_effects)]
df = df[['snp_gene', 'beta', 'beta_se', 'feature_id']]

dfs = []

# For each eGene, randomly extract variants
for i, gene in enumerate(genes):
    print(f'{i}/{len(genes)} processed')

    filtered = df[df['feature_id'] == gene]

    # If possible, extract 10 random variants
    # If less than 10 variants are available, extract all variants
    if len(filtered) > 0:
        sample_size = 10

        if len(filtered) < 10:
            sample_size = len(filtered)

        random_rows = filtered.sample(sample_size)
        dfs.append(random_rows)

# Concat final dataframe
res_df = pd.concat(dfs)

# Save betas to file
res_beta_df = res_df[['snp_gene', 'beta']].rename(columns={'beta': cell_type})
res_beta_df.to_csv(f'{out_dir}/{cell_type}_chr{chromosome}_betas.txt', sep='\t', index=False)

# Save standard errors to file
res_se_df = res_df[['snp_gene', 'beta_se']].rename(columns={'beta_se': cell_type})
res_se_df.to_csv(f'{out_dir}/{cell_type}_chr{chromosome}_standard_errors.txt', sep='\t', index=False)
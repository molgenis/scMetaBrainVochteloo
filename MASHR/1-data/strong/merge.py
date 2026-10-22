import pandas as pd 
import sys 

strong_effects_dir = sys.argv[1]

cell_types = [
    'AST',
    'END',
    'EX',
    'IN',
    'MIC',
    'OLI',
    'OPC',
    'PER',
]

chromosomes = [i for i in range(1,23)]

res_beta_dfs = []
res_se_dfs = []

# Loop through cell types
for cell_type in cell_types:
    print(cell_type)

    # Create empty lists for betas and standard errors
    beta_dfs = []
    se_dfs = []

    # Loop through chromosomes
    for chromosome in chromosomes:

        path = f'{strong_effects_dir}/{cell_type}_chr{chromosome}.txt'

        # Load data
        df = pd.read_csv(path, sep='\t')
        print(len(df))

        # Create beta dataframe
        beta_df = df[['snp_gene', 'beta']].rename(columns={'beta': cell_type})
        beta_dfs.append(beta_df)

        # Create standard error dataframe
        se_df = df[['snp_gene', 'beta_se']].rename(columns={'beta_se': cell_type})
        se_dfs.append(se_df)

    cell_type_beta_df = pd.concat(beta_dfs)
    res_beta_dfs.append(cell_type_beta_df)
    cell_type_se_df = pd.concat(se_dfs)
    res_se_dfs.append(cell_type_se_df)


# Save strong betas to file
res_beta_df = res_beta_dfs[0]

for df in res_beta_dfs[1:]:
    res_beta_df = pd.merge(res_beta_df, df, on='snp_gene', how='outer')

res_beta_df.to_csv('strong_betas.txt', sep='\t', index=False)

# Save strong standard errors to file
res_se_df = res_se_dfs[0]

for df in res_se_dfs[1:]:
    res_se_df = pd.merge(res_se_df, df, on='snp_gene', how='outer')

res_se_df.to_csv('strong_standard_errors.txt', sep='\t', index=False)
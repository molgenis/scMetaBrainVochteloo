import pandas as pd 
import glob 
import sys 

random_effects_dir = sys.argv[1]

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


# BETAS
dfs = []

for cell_type in cell_types:
    print(cell_type)

    files = sorted(glob.glob(f'{random_effects_dir}/{cell_type}_chr*_betas.txt'))
    cell_type_dfs = []

    for file in files:
        df = pd.read_csv(file, sep='\t')
        cell_type_dfs.append(df)

    cell_type_df = pd.concat(cell_type_dfs)
    dfs.append(cell_type_df)

res_df = dfs[0]

for df in dfs[1:]:
    res_df = pd.merge(res_df, df, on='snp_gene', how='outer')

res_df.drop_duplicates('snp_gene', keep='first', inplace=True)
res_df.to_csv('random_betas.txt', sep='\t', index=False)


# STANDARD ERRORS
dfs = []

for cell_type in cell_types:
    print(cell_type)

    files = sorted(glob.glob(f'{random_effects_dir}/{cell_type}_chr*_standard_errors.txt'))
    cell_type_dfs = []

    for file in files:
        df = pd.read_csv(file, sep='\t')
        cell_type_dfs.append(df)

    cell_type_df = pd.concat(cell_type_dfs)
    dfs.append(cell_type_df)

res_df = dfs[0]

for df in dfs[1:]:
    res_df = pd.merge(res_df, df, on='snp_gene', how='outer')

res_df.drop_duplicates('snp_gene', keep='first', inplace=True)
res_df.to_csv('random_standard_errors.txt', sep='\t', index=False)
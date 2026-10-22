import pandas as pd 
import numpy as np 

datasets = [
    'AST',
    'END',
    'EX',
    'IN',
    'MIC',
    'OLI',
    'OPC',
    'PER',
]

dfs = []

# Load correct row of every cell type output
for dataset in datasets:
    file = f'output/{dataset}.txt'

    df = pd.read_csv(file, sep='\t', index_col=0)
    dfs.append(df.loc[[dataset]])

res_df = pd.concat(dfs)
print(res_df)

# Create an ordered dataframe
res_df_ordered = pd.DataFrame(index=datasets, columns=datasets)

for dataset1 in datasets:
    for dataset2 in datasets:
        res_df_ordered.loc[dataset1, dataset2] = res_df.loc[dataset1, dataset2]

# Transform to numeric and transpose
res_df_ordered = res_df_ordered.apply(pd.to_numeric)
res_df_ordered = res_df_ordered.transpose()

# Save to file
res_df_ordered.to_csv('mashr_sharing_all_cell_types.txt', sep='\t')

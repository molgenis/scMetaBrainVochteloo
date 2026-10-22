#!/usr/bin/env python3

import sys 
import pandas as pd

gene = sys.argv[1]
chromosome = sys.argv[2]
eqtl_path = sys.argv[3]
sample_size = int(sys.argv[4])
tissue = sys.argv[5]

# Extract dataset name from path
dataset = eqtl_path.split('/')[-4]

res_df = None

# Loop through eQTL file and extract all rows for this gene
for df in pd.read_csv(eqtl_path, sep='\t', chunksize=10000):
    df = df[df['Gene'] == gene]

    if len(df) > 0:
        print(df.head())
        if res_df is None:
            res_df = df
        else:
            res_df = pd.concat([res_df, df])

# Sort by position
res_df = res_df.sort_values(by='SNPPos')

# Calculate call-rate
res_df['call_rate'] = res_df['MetaPN'] / sample_size

# Filter call rate
res_df = res_df[res_df['call_rate'] > 0.98]

# Save to file
res_df.to_csv(f'{gene}_eqtls.txt', sep='\t', index=False)

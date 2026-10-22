#!/usr/bin/env python3 

import sys
import pandas as pd 

# Parse command line arguments
gwas_path = sys.argv[1]
chromosome_column = sys.argv[2]
position_column = sys.argv[3]
p_column = sys.argv[4]
sig_cutoff = float(sys.argv[5])
window_size = int(sys.argv[6])

sig_df = None

# Loop through GWAS sum stats and extract all significant variants
for df in pd.read_csv(gwas_path, sep='\t', index_col=None, chunksize=10000, engine='python'):
    df[p_column] = pd.to_numeric(df[p_column], errors ='coerce')
    filtered = df[df[p_column] < sig_cutoff]

    if len(filtered) > 0:
        
        if sig_df is None:
            sig_df = filtered
        else:
            sig_df = pd.concat([sig_df, filtered])

        
# Write dataframe with significant variants to file
sig_df.to_csv(f"{'.'.join(gwas_path.split('/')[-1].replace('.gz','').split('.')[:-1])}_significant.txt", sep='\t', index=False)

# Create columns with intervals around the positions of significant variants
sig_df['lower_bound'] = sig_df[position_column].astype(int) - (window_size / 2)
sig_df['upper_bound'] = sig_df[position_column].astype(int) + (window_size / 2)
interval_df = sig_df[[chromosome_column, 'lower_bound', 'upper_bound']]
interval_df['lower_bound'] = interval_df['lower_bound'].astype(int)
interval_df['upper_bound'] = interval_df['upper_bound'].astype(int)

# Write interval data to file
interval_df.to_csv('intervals.txt', sep='\t', index=False, header=None)
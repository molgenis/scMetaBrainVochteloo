import pandas as pd 
import glob 

dfs = []

# Get multi variant results
mv_coloc_out_dir = 'path/to/multi_variant/output'
mv_files = glob.glob(f'{mv_coloc_out_dir}/*/*/*_mv_coloc_out_formatted.txt')

for file in mv_files:
    gwas = file.split('/')[-3]
    df = pd.read_csv(file, sep='\t')
    df['GWAS'] = gwas
    dfs.append(df)

# Get single variant results
sv_coloc_out_dir = 'path/to_single_variant/output'
sv_files = glob.glob(f'{sv_coloc_out_dir}/*/*/*_sv_coloc_out_formatted')

for file in sv_files:
    gwas = file.split('/')[-3]
    df = pd.read_csv(file, sep='\t')
    df['GWAS'] = gwas
    dfs.append(df)


# Merge all results
res_df = pd.concat(dfs)

# Write to file
res_df.to_csv('coloc_results_all_cell_types_all_traits.txt', sep='\t', index=False)
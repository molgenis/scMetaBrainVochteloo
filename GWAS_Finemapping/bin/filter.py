#!/usr/bin/env python3

import pandas as pd 
import glob, os 
import numpy as np
import sys
import tabix

locus = sys.argv[1]
susie_file = sys.argv[2]
vcf_path = sys.argv[3]
chromosome = sys.argv[4]
min_p_threshold = float(sys.argv[5])

def remove_insignificant_p_cs(df):
    '''Function that removes credible sets that have no signficicant variant in it'''
    min_p_threshold = 0.00000005

    no_cs_df = df[df['SusieRss_CS'].isna()]
    cs_df = df[~df['SusieRss_CS'].isna()]

    cs_dfs = []

    # Loop through credible sets and set credible set col to nan if there is no significant variant
    for cs in cs_df['SusieRss_CS'].unique():
        cs_filtered = cs_df[cs_df['SusieRss_CS'] == cs]
        
        if cs_filtered['p_value'].min() > min_p_threshold:
            cs_filtered['SusieRss_CS'] = np.nan

        cs_dfs.append(cs_filtered)

    no_cs_df['SusieRss_CS'] = np.nan
    res_df = pd.concat([no_cs_df] + cs_dfs)
    return res_df


def remove_low_lbf_cs(df):
    '''Function that removes credible sets that have no variant with LBF >= 2 in it'''

    no_cs_df = df[df['SusieRss_CS'].isna()].copy()
    cs_df = df[~df['SusieRss_CS'].isna()].copy()

    cs_dfs = []

    # Loop through credbile sets and set credible set col to nan if there is no variant with LBF >= 2
    for cs in cs_df['SusieRss_CS'].unique():
        filtered = cs_df[cs_df['SusieRss_CS'] == cs]
        lbs_cols = filtered[[col for col in filtered.columns if 'lbf_' in col]]
        max_lbf = lbs_cols.max().max()

        if max_lbf < 2:
            filtered['SusieRss_CS'] = np.nan

        cs_dfs.append(filtered)

    res_df = pd.concat(cs_dfs + [no_cs_df], ignore_index=True)

    return res_df


def count_genotypes(gt_string):
    '''Function that converts a genotype count string to a numerical value'''

    if gt_string.isnumeric():
        return gt_string
    
    if gt_string in ['0|0', '0/0']:
        return 0
    if gt_string in ['0|1', '1|0', '1/0', '0/1']:
        return 1
    if gt_string in ['1|1', '1/1']:
        return 2
    
    return np.nan


def remove_ld_cs(df):
    '''Function that removes the least significant credible set(s) if multiple credible sets are in LD with each other'''

    ld_value_threshold = 0.5

    # Split data frame into variants with and without a credible set
    no_cs_df = df[df['SusieRss_CS'].isna()].copy()
    cs_df = df[~df['SusieRss_CS'].isna()].copy()

    top_variant_dict = {}
    top_variants = []

    # Get top variant for each credible set
    for cs in cs_df['SusieRss_CS'].unique():
        filtered = cs_df[cs_df['SusieRss_CS'] == cs]
        top_variant_row = filtered.nsmallest(1, 'p_value').iloc[0]
        top_variants.append(top_variant_row['ref_variant'])
        top_variant_dict[cs] = {
            'variant_id': top_variant_row['ref_variant'],
            'position': top_variant_row['position'],
            'p_value': top_variant_row['p_value'],
        }


    # Sort by p-value
    top_variant_dict = dict(
        sorted(
            top_variant_dict.items(),
            key=lambda x: x[1]['p_value'],
            reverse=True
        )
    )

    tb = tabix.open(vcf_path)

    genotype_dict = {}
    credible_sets = []

    # For each variant, extract genotypes from VCF
    for cs, variant in top_variant_dict.items():
        records = tb.query(str(chromosome), variant['position'] - 1, variant['position'] + 1)

        for elems in records:
            record_variant = elems[2]
            
            if not record_variant in top_variants:
                continue
            
            credible_sets.append(cs)
            format_field = elems[8]
            gt_col_index = format_field.split(':').index('GT')
            genotypes = elems[9:]
            gt_counts = [count_genotypes(gt.split(':')[gt_col_index]) for gt in genotypes]
            genotype_dict[cs] = gt_counts

    # Convert to dataframe and calculate LD matrix
    genotype_df = pd.DataFrame(genotype_dict)
    ld_df = genotype_df.corr()
    print(ld_df)

    cs_to_drop = []

    # Loop through all credible set pairs
    for cs_i in credible_sets:
        for cs_j in credible_sets:
            if cs_i == cs_j:
                continue

            # Get LD value from LD matrix
            ld_value = abs(ld_df.loc[cs_i, cs_j])

            # If LD value > threshold, remove credible set with highest p-value
            if ld_value > ld_value_threshold:
                if top_variant_dict[cs_i]['p_value'] > top_variant_dict[cs_j]['p_value']:
                    cs_to_drop.append(cs_i)
                else:
                    cs_to_drop.append(cs_j)

    print(f'Credible sets to remove: {cs_to_drop}')

    # Remove credible sets that are in LD with another credible set
    df.loc[df['SusieRss_CS'].isin(cs_to_drop), 'SusieRss_CS'] = np.nan

    return df
    
def reorder_cs(df):
    '''Function that reorders credible sets from most significant to least significant'''
    
    no_cs_df = df[df['SusieRss_CS'].isna()].copy()
    cs_df = df[~df['SusieRss_CS'].isna()].copy()

    # No reordering needed if there are no credible sets
    if len(cs_df) == 0:
        return df

    cs_out = []
    min_p_out = []

    # For each credible set, get lowest p-value
    for cs in cs_df['SusieRss_CS'].unique():
        filtered = cs_df[cs_df['SusieRss_CS'] == cs]
        min_p = filtered['p_value'].min()
        cs_out.append(cs)
        min_p_out.append(min_p)

    # Create dataframe of credible set top variants and their p-values
    cs_order_df = pd.DataFrame({'cs': cs_out, 'min_p': min_p_out})
    cs_order_df.sort_values('min_p', inplace=True)

    cs_mapping_dict = {}
    cs_counter = 1

    # Create a mapping dictionary for reordering credible sets
    for _, row in cs_order_df.iterrows():
        cs_mapping_dict[row['cs']] = cs_counter
        cs_counter += 1

    df['SusieRss_CS'] = df['SusieRss_CS'].map(cs_mapping_dict)
    
    # Create empty list for columns
    columns = []
    cols_to_drop = []

    # Reorder columns
    for col in df.columns:
        if not col.startswith('lbf_cs'):
            columns.append(col)
            continue
        
        cs = int(col.split('_')[-1])

        if cs in cs_mapping_dict.keys():
            columns.append(f'lbf_cs_{cs_mapping_dict[cs]}')
            continue

        cols_to_drop.append(col)

    df.drop(columns=cols_to_drop, inplace=True)

    # Assign new reordered columns to dataframe
    df.columns = columns

    return df


def filter_columns(df):
    '''Function that filters relevant columns'''
    
    num_cs_requested = df.iloc[0]['L']

    cs_df = df[~df['SusieRss_CS'].isna()]
    credible_sets = cs_df['SusieRss_CS'].unique()

    # Drop LBF columns for which no credible set exists if asked for more than 1 credible set
    if num_cs_requested != 1:
        lbf_cols = [c for c in cs_df.columns if c.startswith('lbf_cs_')]
        keep_lbf_cols = [f'lbf_cs_{int(cs)}' for cs in credible_sets]
        drop_lbf_cols = [c for c in lbf_cols if c not in keep_lbf_cols]
        df.drop(columns=drop_lbf_cols, inplace=True)

    return df


# Load SuSiE output
df = pd.read_csv(susie_file, sep='\t')

# Don't process a locus if it has < 250 variants
if len(df) < 250:
    print(f'{locus} not included because only {len(df)} variants remain after filtering')
    exit(0)

# Remove credible sets that have no significant p-value
df = remove_insignificant_p_cs(df)

# Remove credible sets that have no LBF value > 2
df = remove_low_lbf_cs(df)

# Remove credible sets that are in LD with another credible set
df = remove_ld_cs(df)

# Order credible sets based on p-value
df = reorder_cs(df)

# Filter relevant columns
df = filter_columns(df)

# Sort by position
df.sort_values('position', inplace=True)

# Write to file
df.to_csv(f'{locus}_susie_filtered.txt.gz', sep='\t', index=False)

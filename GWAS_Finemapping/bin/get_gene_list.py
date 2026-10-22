#!/usr/bin/env python3 

import pandas as pd 
import sys

annotation_path = sys.argv[1]
all_genes_list_path = sys.argv[2]

# Load a list of all tested eQTL genes
with open(all_genes_list_path, 'r') as all_genes_list:
    all_tested_genes = [g.strip() for g in all_genes_list.readlines()]

# Load annotation dataframe
annotation_df = pd.read_csv(annotation_path, sep='\t')
annotation_df = annotation_df[annotation_df['Gene'].isin(all_tested_genes)]
print(annotation_df.head())

# Load the intervals from interval file
with open('intervals.txt', 'r') as interval_file:
    intervals = [i.strip() for i in interval_file.readlines()]

res_df = None

# Loop through intervals and extract genes from the annotation file to test
for interval in intervals:
    chromosome, start, end = interval.split('\t')
    chromosome = int(chromosome)
    start = int(start)
    end = int(end)

    filtered = annotation_df[
        (annotation_df['Chr'] == chromosome) &
        (annotation_df['ChrStart'] > start) &
        (annotation_df['ChrEnd'] < end)
        ]

    if len(filtered) > 0:
        if res_df is None:
            res_df = filtered
        else:
            res_df = pd.concat([res_df, filtered])

# If the resulting dataframe is not empty, write it to file
if res_df is not None:
    res_df = res_df[['Gene', 'Chr']]
    res_df.drop_duplicates('Gene', inplace=True)
    res_df.to_csv('gene_list.txt', sep='\t', index=False, header=None)
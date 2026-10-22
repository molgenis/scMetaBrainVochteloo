"""
Script to compare AlphaGenome predictions between SingleBrain and scMetaBrain top eQTLs.
Requires AlphaGenome API access.

Example usage:
python alphagenome_analysis.py --SINGLEBRAIN_PATH /path/to/SingleBrain_top_eqtls.tsv.gz --METABRAIN_PATH /path/to/scMetaBrain_top_eqtls.txt --METABRAIN_VCF_PATH /path/to/scMetaBrain_ref.vcf.gz --tissue "frontal cortex" --outdir "output_directory" --verbose

Author: Francisco Camoes Magalhaes
"""

#!/usr/bin/env python
from functions import score_variants, load_metabrain, plot_scatter, plot_histogram
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

import argparse

## Parse command line arguments
parser = argparse.ArgumentParser(description="Compare AlphaGenome scores between SingleBrain and scMetaBrain top eQTLs.")
parser.add_argument("--tissue", dest="tissue_name", default="frontal cortex", help="Tissue name used for AlphaGenome RNA-seq predictions.")
parser.add_argument("--outdir", dest="MAIN_DIR", default="singlebrain_test_67", help="Output directory for figures and intermediate files.")
parser.add_argument("--match-gene-mode", dest="match_gene_mode", choices=["none", "same_number", "same_genes"], default="none", help="How to match variants between datasets. 'none': use all variants, 'same_number': match dataset sizes, 'same_genes': restrict to shared genes.")
parser.add_argument("--singlebrain", required=True, dest="SINGLEBRAIN_PATH", help="Path to SingleBrain eQTL summary statistics.")
parser.add_argument("--metabrain", required=True, dest="METABRAIN_PATH", help="Path to scMetaBrain eQTL data.")
parser.add_argument("--metabrain-vcf", required=True, dest="METABRAIN_VCF_PATH", help="Path to the reference VCF used for scMetaBrain variants.")
parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")

args = parser.parse_args()

## Configurations
# Parsed arguments
tissue_name = args.tissue_name
MAIN_DIR = args.MAIN_DIR
match_gene_mode = args.match_gene_mode
verbose = args.verbose

SINGLEBRAIN_PATH = args.SINGLEBRAIN_PATH
METABRAIN_PATH = args.METABRAIN_PATH
METABRAIN_VCF_PATH = args.METABRAIN_VCF_PATH

## Logging setup
# Define a verbose print function
vprint = print if verbose else lambda *args, **kwargs: None

# Create output directories
figs_dir = os.path.join(MAIN_DIR, "figures")
data_dir = os.path.join(MAIN_DIR, "data")
os.makedirs(MAIN_DIR, exist_ok=True)
os.makedirs(figs_dir, exist_ok=True)
os.makedirs(data_dir, exist_ok=True)

# Define AlphaGenome scoring parameters
REQUESTED_SCORERS = ["AVI_SCORE", "RNA_SEQ", "AVI_SCORE_FEATURE_IMPORTANCE"]
DEFAULT_TISSUES = {tissue_name: tissue_name}
REQUIRED_COLUMNS = ["variant_id", "CHROM", "POS", "REF", "ALT"]

## Load SingleBrain top eQTL data
# SingleBrain configuration
SINGLEBRAIN_COI = ['feature', 'variant_id', 'chr', 'pos', 'ref', 'alt', 'fixed_beta', 'fixedfxqval']
SINGLEBRAIN_RENAME = {'chr': 'CHROM', 'pos': 'POS', 'ref': 'REF', 'alt': 'ALT', 'fixed_beta': 'BETA', 'fixedfxqval': 'FDR'}

# Load SingleBrain data
vprint(f"Loading SingleBrain data from\n{SINGLEBRAIN_PATH}")
df_singlebrain = pd.read_csv(SINGLEBRAIN_PATH, sep='\t')

## Filter and process SingleBrain data
vprint(f"Filtering dataframe...")
# Add a new column 'variant_id' that combines the gene and the variant_id
df_singlebrain['variant_id'] = df_singlebrain['feature'].str.split('.').str[0] + '_' + df_singlebrain['variant_id']
df_singlebrain = df_singlebrain[SINGLEBRAIN_COI]
# Subset to columns of interest and filter by FDR < 0.05
df_singlebrain = df_singlebrain.loc[df_singlebrain['fixedfxqval'] < 0.05]
# Drop rows with any missing values
df_singlebrain = df_singlebrain.dropna()
# Subset to only the top eQTL per gene, keeping the one with the lowest FDR
df_singlebrain = df_singlebrain.loc[df_singlebrain.groupby('feature')['fixedfxqval'].apply(lambda x: x.idxmin())]
# Rename columns to match the expected format
df_singlebrain = df_singlebrain.rename(columns=SINGLEBRAIN_RENAME)
# Subset to the final columns of interest
df_singlebrain = df_singlebrain[["variant_id", "CHROM", "POS", "REF", "ALT", "BETA"]]

# Save the eQTL data to a CSV file for reference
save_path= os.path.join(data_dir, 'SingleBrain_top_eqtls.csv')
vprint(f"Saving SingleBrain data to {save_path}")
df_singlebrain.to_csv(save_path, index=False)

vprint(f"Done.")

## Load scMetaBrain top eQTL data
# scMetaBrain configuration
METABRAIN_COI = ["ensembl_id", "snp_id", "qval", "beta", "beta_se", "snp_chromosome", "snp_position", "assessed_allele"]
METABRAIN_RENAME = {
    "snp_id": "variant_id",
    "qval": "FDR",
    "beta": "BETA",
    "beta_se": "SE",
    "snp_chromosome": "CHROM",
    "snp_position": "POS",
    "assessed_allele": "ALT",
    "ensembl_id": "gene_id",
}
COL_ORDER = ["variant_id", "CHROM", "POS", "REF", "ALT", "BETA", "SE", "FDR"]

# Load scMetaBrain data using the provided function
vprint(f"Loading scMetaBrain data from\n{METABRAIN_PATH}\nand\n{METABRAIN_VCF_PATH}")
df_metabrain = load_metabrain(METABRAIN_PATH, METABRAIN_VCF_PATH, cols_of_interest=METABRAIN_COI, rename_dict=METABRAIN_RENAME, col_order=COL_ORDER)
## Filter and process scMetaBrain data
vprint(f"Filtering dataframe...")
# Drop rows with any missing values in the specified columns
df_metabrain = df_metabrain.dropna(subset=['CHROM', 'POS', 'ALT', 'BETA', 'SE', 'FDR'])
# Extract gene_id from variant_id for grouping
df_metabrain["gene_id"] = df_metabrain["variant_id"].str.split('_').str[0]

# Subset to only the top eQTL per gene, keeping the one with the lowest FDR
df_metabrain = df_metabrain.loc[df_metabrain.groupby('gene_id')['FDR'].apply(lambda x: x.idxmin())]
df_metabrain = df_metabrain.drop("gene_id", axis=1)
# Save the eQTL data to a CSV file for reference
save_path= os.path.join(data_dir, 'scMetaBrain_top_eqtls.csv')
vprint(f"Saving scMetaBrain data to {save_path}")
df_metabrain.to_csv(save_path, index=False)


## Score the variants using AlphaGenome API
vprint(f"Scoring variants for tissues:\n{DEFAULT_TISSUES.keys()}")
df_metabrain = score_variants(df_metabrain, requested_scorers=REQUESTED_SCORERS, required_columns=REQUIRED_COLUMNS, tissues_of_interest_col_names=DEFAULT_TISSUES, batch_size=None, max_retries=500, sleep_time=5, return_skips=False, cached_loc=os.path.join(data_dir, "ag_rna_scores_cache.tsv.gz"))
df_singlebrain = score_variants(df_singlebrain, requested_scorers=REQUESTED_SCORERS, required_columns=REQUIRED_COLUMNS, tissues_of_interest_col_names=DEFAULT_TISSUES, batch_size=None, max_retries=500, sleep_time=5, return_skips=False, cached_loc=os.path.join(data_dir, "ag_rna_scores_cache.tsv.gz"))
df_metabrain = df_metabrain.dropna(subset=["AVI_SCORE", tissue_name, "BETA"])
df_singlebrain = df_singlebrain.dropna(subset=["AVI_SCORE", tissue_name, "BETA"])

## Subset the datasets based on the match_gene_mode
if match_gene_mode == "same_genes":
    # Extract gene_id from variant_id for both datasets
    df_metabrain["gene_id"] = df_metabrain["variant_id"].str.split('_').str[0]
    df_singlebrain["gene_id"] = df_singlebrain["variant_id"].str.split('_').str[0]

    common_genes = set(df_metabrain["gene_id"]).intersection(set(df_singlebrain["gene_id"]))
    df_metabrain = df_metabrain.loc[df_metabrain["gene_id"].isin(common_genes)].sort_values("gene_id")
    df_singlebrain = df_singlebrain.loc[df_singlebrain["gene_id"].isin(common_genes)].sort_values("gene_id")

elif match_gene_mode == "same_number":
    n_rows = len(df_singlebrain)
    df_metabrain = df_metabrain.nsmallest(n_rows, "FDR")

elif match_gene_mode == "none":
    pass
else:
    raise ValueError(f"Unknown match_gene_mode: {match_gene_mode}")



## Plot the results
# Plot in two subplots, side by side

vprint(f"Making plots...")
fig, axs = plt.subplots(3, 2, figsize=(18, 20))
# Plot scatter plots of BETA vs RNA-seq scores for both datasets
plot_scatter("BETA", tissue_name, df_metabrain, SOURCE="scMetaBrain", ax = axs[0,0], draw_regression_line=True)
plot_scatter("BETA", tissue_name, df_singlebrain, SOURCE="SingleBrain", ax = axs[0,1], draw_regression_line=True)

# Plot only top 5% of absolute RNA-seq scores for both datasets
# Subset
df_metabrain["absolute_RNA"] = df_metabrain[tissue_name].abs()
df_singlebrain["absolute_RNA"] = df_singlebrain[tissue_name].abs()
n_top = max(1, int(0.05 * len(df_metabrain)))
df_metabrain_subset = df_metabrain.nlargest(n_top, "absolute_RNA")
n_top = max(1, int(0.05 * len(df_singlebrain)))
df_singlebrain_subset = df_singlebrain.nlargest(n_top, "absolute_RNA")

# Plot
plot_scatter("BETA", tissue_name, df_metabrain_subset, title_suffix="Top 5% absolute RNA-seq", SOURCE="scMetaBrain", ax = axs[1,0], draw_regression_line=True)
plot_scatter("BETA", tissue_name, df_singlebrain_subset, title_suffix="Top 5% absolute RNA-seq", SOURCE="SingleBrain", ax = axs[1,1], draw_regression_line=True)


# Calculate Wilcoxon test between the two datasets for the absolute RNA-seq values
if match_gene_mode == "same_genes":
    df_metabrain = df_metabrain.sort_values("gene_id")
    df_singlebrain = df_singlebrain.sort_values("gene_id")
    stat, p = wilcoxon(df_metabrain["absolute_RNA"], df_singlebrain["absolute_RNA"])
    vprint(f"Wilcoxon test between scMetaBrain and SingleBrain for\nabsolute RNA-seq values: W={stat}, p={p}")
    suffix=f"Wilcoxon test between scMetaBrain and SingleBrain for\nabsolute RNA-seq values: W={stat}, p={p}"
else:
    suffix = ''

# Plot histograms of the RNA-seq scores for both datasets
plot_histogram(df_metabrain, tissue_name, title_suffix=suffix, ax=axs[2,0], SOURCE="scMetaBrain")
plot_histogram(df_singlebrain, tissue_name, title_suffix=suffix, ax=axs[2,1], SOURCE="SingleBrain")


# set y-axis of histograms to be the same for both subplots
max_y = max(axs[2,0].get_ylim()[1], axs[2,1].get_ylim()[1])
min_y = min(axs[2,0].get_ylim()[0], axs[2,1].get_ylim()[0])
axs[2,0].set_ylim(min_y, max_y)
axs[2,1].set_ylim(min_y, max_y)
plt.tight_layout()

# Save the figure
filename = f"direct_comparison_mode_{match_gene_mode}.pdf"
fileloc = os.path.join(figs_dir, filename)
vprint(f"Saving plots to {fileloc}")
plt.savefig(fileloc)
plt.close()
vprint("Saved plots for scMetaBrain and SingleBrain BETA vs RNA-seq scores and histograms.")
print("Done.")

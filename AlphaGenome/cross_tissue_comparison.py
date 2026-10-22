"""
Script to compare how scMetaBrain eQTLs of different tissues correlate with tissue-specific scores from AlphaGenome.
Requires AlphaGenome API access.

Example usage:
python cross_tissue_comparison.py --vcf path/to/scMetaBrain_ref.vcf.gz --cell_type_dir path/to/scMetaBrain_cell_type_results/ --corr_method spearman --save_intermediate --verbose

Author: Francisco Camoes Magalhaes
"""
#!/usr/bin/env python
import os
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from functions import (
    score_variants,
    load_metabrain,
    plot_summary_figure,
)

# argparse
parser = argparse.ArgumentParser(
    description="Score scMetaBrain variants and generate summary figures."
)

parser.add_argument(
    "--corr_method",
    type=str,
    default="spearman",
    choices=["spearman", "pearson"],
    help="Correlation method.",
)

parser.add_argument(
    "--o",
    type=str,
    default="cross_tissue_comparison",
    help="Directory for output files.",
)

parser.add_argument(
    "--save_intermediate",
    action="store_true",
    help="Save scored TSV files.",
)

parser.add_argument(
    "--verbose",
    action="store_true",
    help="Verbose output.",
)

parser.add_argument(
    "--vcf_path",
    type=str,
    required=True,
    help="Path to the VCF file containing the scMetaBrain variant metadata."
)

parser.add_argument(
    "--cell_type_dir",
    type=str,
    required=True,
    help="Path to the directory containing the scMetaBrain cell type results.",
)


args = parser.parse_args()

corr_method = args.corr_method
MAIN_DIR = args.o
save_intermediate = args.save_intermediate
vprint = print if args.verbose else lambda *a, **k: None
VCF_PATH = args.vcf_path
_BASE = args.cell_type_dir

os.makedirs(MAIN_DIR, exist_ok=True)


# Define Constants

RNA_TRACKS = {
    "frontal cortex": "Brain_Frontal_Cortex_RNA_Seq",
    "CD14-positive monocyte": "CD14_pos_monocyte_RNA_Seq",
    "fibroblast of lung": "Lung_fibroblast_RNA_Seq",
}

tissue_dict = {
    "Brain_Frontal_Cortex_RNA_Seq": "Frontal Cortex",
    "CD14_pos_monocyte_RNA_Seq": "Monocytes",
    "Lung_fibroblast_RNA_Seq": "Lung fibroblast",
}

tissues = list(tissue_dict.keys())

METABRAIN_COI = [
    "ensembl_id",
    "snp_id",
    "qval",
    "beta",
    "beta_se",
    "snp_chromosome",
    "snp_position",
    "assessed_allele",
]

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

COL_ORDER = [
    "variant_id",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "BETA",
    "SE",
    "FDR",
]

cell_types = ["AST", "EX", "IN", "MIC", "OLI", "OPC"]


# scoring

def load_and_score(
    name: str,
    loc: str,
    save_variants: bool = True,
) -> pd.DataFrame:

    vprint(f"\nLoading variants for {name} from {loc}...")

    df = load_metabrain(
        loc,
        VCF_PATH,
        cols_of_interest=METABRAIN_COI,
        rename_dict=METABRAIN_RENAME,
        col_order=COL_ORDER,
    )

    vprint("Scoring variants...")

    df = score_variants(
        df,
        tissues_of_interest_col_names=RNA_TRACKS,
    )

    df = df.dropna(subset=RNA_TRACKS.values())

    if save_variants:
        outfile = os.path.join(MAIN_DIR, f"{name}_scored.tsv")
        vprint(f"Saving scored variants to {outfile}")
        df.to_csv(outfile, sep="\t", index=False)

    vprint("Done.")
    return df



# load data

loc_dict = {
    f"scMetaBrain_{ct}":
    f"{_BASE}/{ct}/top_qtl_results_all_MultiTest.txt"
    for ct in cell_types
}

dfs = {}

for name, loc in loc_dict.items():
    region = name.replace("scMetaBrain_", "")

    vprint(f"Processing {region}")

    df = load_and_score(
        name=name,
        loc=loc,
        save_variants=save_intermediate,
    )

    df = df.dropna()
    df = df.drop_duplicates(subset=["variant_id"])

    dfs[region] = df



# Plotting

vprint("Subsetting dataframes to the same number of rows.")

min_size = min(len(df) for df in dfs.values())

for region, df in dfs.items():
    dfs[region] = df.nsmallest(min_size, "FDR")

vprint("Calculating correlations between BETA and tissue-specific scores.")

results = []

for threshold in [0.1, 0.05]:
    for source in dfs.keys():

        df = dfs[source]

        vprint(
            f"Correlations for {source} "
            f"at top {threshold * 100:.0f}%"
        )

        for tissue in tissues:

            df[f"abs_{tissue}"] = df[tissue].abs()

            top_n = int(len(df) * threshold)

            top_df = df.nlargest(top_n, f"abs_{tissue}")

            corr = top_df["BETA"].corr(
                top_df[tissue],
                method=corr_method,
            )

            vprint(
                f"    {tissue}: {corr} "
                f"(n={len(top_df)})"
            )

            results.append(
                {
                    "source": source,
                    "tissue": tissue,
                    "threshold": threshold,
                    "n": len(top_df),
                    f"{corr_method}_corr": corr,
                }
            )

plot_summary_figure(
    dfs,
    results,
    os.path.join(
        MAIN_DIR,
        f"summary_figure_{corr_method}.pdf",
    ),
    threshold=0.05,
    corr_method=corr_method,
    tissue_dict=tissue_dict,
)

print("Done generating summary figure.")
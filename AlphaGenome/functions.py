#!/usr/bin/env python
import gc
import os
import re
import time
from typing import Any

import anndata as ad
import grpc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score
import tabix
from tqdm import tqdm

from alphagenome.atlas import atlas
from alphagenome.data import genome

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

cell_type_colors = {
    "AST": "#ba611b",
    "END": "#cb79a6",
    "EX": "#5bb3e4",
    "IN": "#166c9c",
    "MIC": "#e59f24",
    "OLI": "#148a6a",
    "OPC": "#80d6b9",
    "PER": "#f7c4e3"
}


def _batch(
    batch: list[genome.Variant],
    batched_results: dict[str, list[ad.AnnData]],
    atlas_client: Any,
    requested_scorers: dict[str, str],
    skips: int = 0,
    max_retries: int = 50,
    sleep_time: int = 60,
):
    success = False
    retry_count = 0

    while not success and retry_count < max_retries:
        try:
            batch_result = atlas_client.query_variants(batch, requested_scorers=requested_scorers)
            for scorer in requested_scorers:
                adata = batch_result[scorer]
                adata.obs_names_make_unique()
                batched_results[scorer].append(adata)
            success = True
        except Exception as e:
            if isinstance(e, grpc.RpcError):
                retry_count += 1
                if e.code() in (
                    grpc.StatusCode.RESOURCE_EXHAUSTED,
                    grpc.StatusCode.UNAVAILABLE,
                ):
                    print(f"Quota exceeded. Sleeping {sleep_time} seconds before retry.")
                    time.sleep(sleep_time)
                else:
                    raise

            elif isinstance(e, ValueError):
                if len(batch) == 1:
                    skips += 1
                    print(f"Skipping variant {batch[0].name} due to error: {e}. Total skips: {skips}")
                else:
                    mid = len(batch) // 2
                    batched_results, skips = _batch(batch[:mid], batched_results, atlas_client, requested_scorers, skips, max_retries=max_retries, sleep_time=sleep_time)
                    batched_results, skips = _batch(batch[mid:], batched_results, atlas_client, requested_scorers, skips,max_retries=max_retries, sleep_time=sleep_time)
                success = True
    if not success:
        raise RuntimeError(f"Failed to process batch after {max_retries} retries.")

    return batched_results, skips



def score_variants(
    vcf: pd.DataFrame,
    cached_loc: str = "data/rna_scores_cache.tsv.gz",
    tissues_of_interest_col_names: dict[str, str] | None = None,
    required_columns: list[str] = [
        "variant_id",
        "CHROM",
        "POS",
        "REF",
        "ALT",
    ],
    batch_size: int | None = 1000,
    max_retries: int = 500,
    sleep_time: int = 60,
    requested_scorers: list[str] = [
        "AVI_SCORE",
        "RNA_SEQ",
        "AVI_SCORE_FEATURE_IMPORTANCE",
    ],
    return_skips: bool = False,
):


    if tissues_of_interest_col_names is None:
        tissues_of_interest_col_names = {"frontal cortex": "Brain_Frontal_Cortex_RNA-Seq"}
    requested_scorers = list(requested_scorers)


    for column in required_columns:
        if column not in vcf.columns:
            raise ValueError(f'VCF file is missing required column: {column}.')

    # remove non-SNVs
    vcf = vcf[(vcf['REF'].str.len() == 1) & (vcf['ALT'].str.len() == 1)].copy()
    snv_vcf = vcf.copy()
    input_cols = set(snv_vcf.columns)
 
    API_KEY = os.environ.get("ALPHAGENOME_API_KEY")
    if API_KEY is None:
        raise ValueError("ALPHAGENOME_API_KEY environment variable is not set. Please set it to your AlphaGenome API key.")
    atlas_client = atlas.create(API_KEY)

    desired_variants = vcf['variant_id'].tolist()
 
    if os.path.exists(cached_loc):
        cached_df = pd.read_csv(cached_loc, sep='\t')
 
        # check if all the requested RNA tracks are in the cached DataFrame. If not,
        # we need to reqeust them again from the API, so we cant use them.
        if all(col in cached_df.columns for col in tissues_of_interest_col_names.values()):
            # Remove the cached variants from the vcf DataFrame
            cached_variants = set(cached_df['variant_id'])
            vcf = vcf[~vcf['variant_id'].isin(cached_variants)]
        else:
            print("Cached DataFrame does not contain all requested RNA tracks. Re-querying API for all variants.")
            cached_df = None
    else:
        cached_df = None
 
    to_score = []
    for index, row in vcf.iterrows():
        variant = genome.Variant(
            chromosome=row['CHROM'],
            position=int(row['POS']),
            reference_bases=row['REF'],
            alternate_bases=row['ALT'],
            name=row['variant_id'],
        )
        to_score.append(variant)


    if len(to_score) == 0:
        print("No new variants to score. Exiting.")
        if cached_df is None:      # empty input
            out = snv_vcf.reset_index(drop=True)
        else:
            out = snv_vcf.merge(cached_df, on='variant_id', how='left').reset_index(drop=True)
        
        if return_skips:
            return out, None
        return out
            
            
 
    batched_results = {scorer: [] for scorer in requested_scorers}

    # if batch_size is None, set it to the length of to_score to process all variants a single batch
    if batch_size is None:
        batch_size = len(to_score)
 
    for i in tqdm(range(0, len(to_score), batch_size), desc="Querying API batches", smoothing=0):
        batch = to_score[i:i + batch_size]
        batched_results, skips = _batch(batch, batched_results, atlas_client, requested_scorers, max_retries=max_retries, sleep_time=sleep_time)
 
    final_results = {}
    for scorer in requested_scorers:
        final_results[scorer] = ad.concat(batched_results[scorer], merge="first", index_unique="-")
 
    results = final_results
 
    rna_adata = results["RNA_SEQ"]
    avi_adata = results["AVI_SCORE"]
    feat_adata = results["AVI_SCORE_FEATURE_IMPORTANCE"]
 
    vcf["extracted_gene_id"] = vcf["variant_id"].str.extract(r'(ENSG\d+)', expand=False)
    vcf["constructed_variant_str"] = vcf["CHROM"] + ":" + vcf["POS"].astype(str) + ":" + vcf["REF"] + ">" + vcf["ALT"]
 
    rna_adata.obs["clean_gene_id"] = rna_adata.obs["gene_id"].astype(str).str.extract(r'(ENSG\d+)', expand=False)
    rna_adata.obs["variant_str"] = rna_adata.obs["variant"].astype(str)
    avi_adata.obs["variant_str"] = avi_adata.obs["variant"].astype(str)
    feat_adata.obs["variant_str"] = feat_adata.obs["variant"].astype(str)
 
    with np.errstate(all='ignore'):
        avi_vals = avi_adata.X[:, 0]
        avi_adata.obs["AVI_SCORE"] = avi_vals
        phred_scores = -10 * np.log10(1 - avi_vals)
        phred_scores[avi_vals == 1.0] = np.inf
        avi_adata.obs["AVI_SCORE_PHRED"] = phred_scores
 
    avi_df = avi_adata.obs[["variant_str", "AVI_SCORE", "AVI_SCORE_PHRED"]].drop_duplicates()
    del avi_adata, avi_vals, phred_scores
 
    feat_df = pd.DataFrame(feat_adata.X, columns=feat_adata.var["name"].values)
    feat_df["variant_str"] = feat_adata.obs["variant_str"].values
    feat_df = feat_df.drop_duplicates(subset=["variant_str"])
    del feat_adata
    gc.collect()
 
    vcf = vcf.merge(avi_df, how="left", left_on=["constructed_variant_str"], right_on=["variant_str"]).drop(columns=["variant_str"])
    vcf = vcf.merge(feat_df, how="left", left_on=["constructed_variant_str"], right_on=["variant_str"]).drop(columns=["variant_str"])
    del avi_df, feat_df
    gc.collect()

    # --- Tissue-specific RNA scores: one column per tissue ---
    for tissue, col_name in tissues_of_interest_col_names.items():
        var_mask = (rna_adata.var["biosample_name"] == tissue) & (rna_adata.var["biosample_life_stage"] == 'adult')

        if var_mask.sum() == 0:
            raise ValueError(f"No RNA tracks found for tissue '{tissue}' after filtering.")
        
        subset_rna = rna_adata[:, var_mask].copy()
 
        with np.errstate(all='ignore'):
            subset_rna.obs[col_name] = np.nanmean(subset_rna.X, axis=1)
 
        rna_df = subset_rna.obs[["clean_gene_id", "variant_str", col_name]].drop_duplicates()
        del subset_rna

        vcf = vcf.merge(rna_df, how="left", left_on=["extracted_gene_id", "constructed_variant_str"], right_on=["clean_gene_id", "variant_str"])
        vcf = vcf.drop(columns=["clean_gene_id", "variant_str"])
        
        del rna_df
        gc.collect()
        

    # rna_adata is only released once every tissue has been processed
    del rna_adata
    gc.collect()
 
    # --- Merge with previously scored variants and update the cache (once) ---
    score_cols = [c for c in vcf.columns if c not in list(input_cols) + ['extracted_gene_id', 'constructed_variant_str']]
    cache_new = vcf[['variant_id'] + score_cols]

    if cached_df is not None:
        cache_out = (pd.concat([cache_new, cached_df], ignore_index=True)
                    .drop_duplicates(subset=['variant_id'], keep='first'))  # fresh wins
    else:
        cache_out = cache_new
    tmp = cached_loc + ".tmp"
    cache_out.to_csv(tmp, sep="\t", index=False, compression='gzip')
    os.replace(tmp, cached_loc)

    vcf = (snv_vcf[snv_vcf['variant_id'].isin(desired_variants)]
            .merge(cache_out, on='variant_id', how='left')
            .reset_index(drop=True))
    
    if return_skips:
        return vcf, skips
    return vcf


def _populate_ref_alleles(
    df: pd.DataFrame,
    vcf_file: str,
):
    tb = tabix.open(vcf_file)

    def get_ref(row):
        chrom = str(row["CHROM"]).replace("chr", "")
        pos = int(row["POS"])
        alt = str(row["ALT"])
        row["REF"] = None  # stays None if no match is found

        try:
            for rec in tb.query(chrom, pos - 1, pos):
                # rec[1] = POS, rec[3] = REF, rec[4] = ALT
                if int(rec[1]) != pos:
                    continue
                if rec[4] == alt:
                    row["REF"] = rec[3]
                    return row
                if rec[3] == alt:
                    # Assessed allele is the VCF REF: swap alleles and flip BETA
                    row["REF"] = rec[3]
                    row["ALT"] = rec[4]
                    row["BETA"] *= -1
                    return row
        except tabix.TabixError as e:
            print(e)

        return row

    return df.apply(get_ref, axis=1)


def load_metabrain(
    path: str,
    vcf_path: str,
    cols_of_interest: list[str],
    rename_dict: dict[str, str],
    col_order: list[str],
    fdr_threshold: float = 0.05,
):
    df = pd.read_csv(path, sep="\t")
    df = df[cols_of_interest].rename(columns=rename_dict)

    df["REF"] = "N"
    df["CHROM"] = "chr" + df["CHROM"].astype(str)

    df = df.loc[df["FDR"] < fdr_threshold].copy()

    df["variant_id"] = df["gene_id"] + "_" + df["variant_id"]
    df = df.drop(columns=["gene_id"])
    df = df[col_order]

    df = _populate_ref_alleles(df, vcf_path)

    return df


def plot_scatter(
    x_col: str,
    y_col: str,
    df: pd.DataFrame,
    title_suffix: str = "",
    draw_regression_line: bool = False,
    main_loc: str = "",
    ax: Axes | None = None,
    SOURCE: str = "",
    savefile_name: str | None = None,
):
    pearson_corr = df[[x_col, y_col]].corr(method='pearson').iloc[0, 1]
    spearman_corr = df[[x_col, y_col]].corr(method='spearman').iloc[0, 1]
    directional_concordance = np.sign(df[x_col]) == np.sign(df[y_col])
    directional_concordance_rate = directional_concordance.mean()
    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)
    plt.scatter(df[x_col], df[y_col], alpha=0.5, c="#8888FF", edgecolors="black", linewidth=0.3, s=50)
    plt.title(f'{SOURCE}\nScatter plot of {y_col} vs {x_col}\nPearson correlation: {pearson_corr:.2f}, Spearman correlation: {spearman_corr:.2f}, Directional concordance: {directional_concordance_rate:.2%}\nn={len(df)}, {title_suffix}')
    plt.xlabel(x_col)
    plt.ylabel(f"AlphaGenome {y_col} Differential RNA-seq Prediction")
    plt.grid(True)
    if draw_regression_line:
        m, b = np.polyfit(df[x_col], df[y_col], 1)
        reg_x = np.linspace(df[x_col].min(), df[x_col].max(), 4)
        reg_y = m * reg_x + b
        plt.plot(reg_x, reg_y, 'r--')

    if ax is None:
        if savefile_name is None:
            savefile_name = f"scatter_{SOURCE}_{y_col}_vs_{x_col}_{title_suffix.lower().replace(' ', '_')}.pdf"
        plt.savefig(os.path.join(main_loc, savefile_name))
        plt.close()



def plot_histogram(
    df: pd.DataFrame,
    col: str,
    title_suffix: str = "",
    ax: Axes | None = None,
    SOURCE: str = "",
    bin_width: float = 0.05,
):
    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)

    bins = np.arange(0, df[col].abs().max() + bin_width, bin_width)
    plt.hist(df[col].abs(), bins=bins, color="#2222FF", edgecolor="black", linewidth=0.5)
    plt.title(f'{SOURCE}\nHistogram of absolute {col}\nn={len(df)}, {title_suffix}')
    plt.xlabel(f'Absolute {col}')
    plt.ylabel('Count')
    plt.yscale('log')
    plt.grid(True)


def plot_summary_figure(dfs, results, outfile, threshold=0.05,
                        scatter_sources=("EX", "MIC"),
                        scatter_tissues=("Brain_Frontal_Cortex_RNA_Seq",
                                         "CD14_pos_monocyte_RNA_Seq"),
                        corr_method="spearman",
                        tissue_dict = {},
                        bar_tissues=None):

    bar_tissues = bar_tissues or list(tissue_dict.keys())

    fig = plt.figure(figsize=(16, 18))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 1], hspace=0.3, wspace=0.2)

    # --- 2x2 scatterplots (all share the same x and y limits) ---
    ax0 = None
    for r, source in enumerate(scatter_sources):
        for c, tissue in enumerate(scatter_tissues):
            if ax0 is None:
                ax = ax0 = fig.add_subplot(gs[r, c])
            else:
                ax = fig.add_subplot(gs[r, c], sharex=ax0, sharey=ax0)

            plot_correlation_by_top_percent(
                dfs[source], "BETA", tissue,
                xlabel=f"{source} BETA",
                ylabel=f"{tissue_dict.get(tissue, tissue)} effect",
                title=f"{source} Beta vs {tissue_dict.get(tissue, tissue)}",
                ax=ax,
                legend_fontsize=8,
                legend_title_fontsize=9,
            )

    # --- barplot spanning the bottom row ---
    ax_bar = fig.add_subplot(gs[2, :])
    sources = list(dfs.keys())
    x = np.arange(len(bar_tissues))
    width = 0.8 / len(sources)
    for i, source in enumerate(sources):
        corrs = [
            next(res[f"{corr_method}_corr"] for res in results
                 if res["source"] == source and res["tissue"] == t
                 and res["threshold"] == threshold)
            for t in bar_tissues
        ]
        ax_bar.bar(x + i * width, corrs, width, label=source, color=cell_type_colors[source])
    ax_bar.set_xticks(x + width * (len(sources) - 1) / 2)
    ax_bar.set_xticklabels([tissue_dict.get(t, t) for t in bar_tissues], rotation=30, ha="right")
    if corr_method == "spearman":
        ax_bar.set_ylabel("Spearman Correlation")
    elif corr_method == "pearson":
        ax_bar.set_ylabel("Pearson Correlation")
    ax_bar.set_title(f"Correlation between BETA and tissue effect (top {int(threshold * 100)}%)")
    ax_bar.axhline(0, color="black", lw=0.8)
    ax_bar.legend(title="Source", ncol=len(sources))

    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


def plot_correlation_by_top_percent(
    df,
    x_col,
    y_col,
    top_percents=(100, 10, 5),
    xlabel=None,
    ylabel=None,
    title=None,
    colors=("grey", "royalblue", "crimson"),
    point_alphas=(0.05, 0.15, 0.3),
    figsize=(10, 8),
    save_path=None,
    show=False,
    ax=None,
    legend_fontsize=12,
    legend_title_fontsize=14,
):
    top_percents = sorted(top_percents, reverse=True)
    if any(p <= 0 or p > 100 for p in top_percents):
        raise ValueError("top_percents values must be in the range (0, 100].")
    if len(colors) < len(top_percents) or len(point_alphas) < len(top_percents):
        raise ValueError("Provide at least one color and alpha per threshold.")

    data = df.dropna(subset=[x_col, y_col]).copy()
    data["_abs_y"] = data[y_col].abs()

    # Convert percentages into counts of top variants to consider
    top_ns = []
    cutoffs = []
    for p in top_percents:
        if p >= 100:
            top_ns.append(len(data))
            cutoffs.append(-np.inf)  # Include all data for 100%
        else:
            top_ns.append(int(len(data) * p / 100))
            cutoffs.append(data.nlargest(top_ns[-1], "_abs_y")["_abs_y"].min())

    # Nested subsets used for regression lines and correlations
    subsets = [data.nlargest(n, "_abs_y") for n in top_ns]

    # Exclusive bins used for coloring the points (so they don't overplot)
    bins = []
    for i, c in enumerate(cutoffs):
        upper = cutoffs[i + 1] if i + 1 < len(cutoffs) else np.inf
        bins.append(data[(data["_abs_y"] >= c) & (data["_abs_y"] < upper)])

    # Use the provided axis, or create our own figure
    created_fig = ax is None
    if created_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # Points
    for i, b in enumerate(bins):
        sns.scatterplot(
            data=b, x=x_col, y=y_col, alpha=point_alphas[i], s=15,
            color=colors[i], zorder=1 + i, edgecolor="none", ax=ax,
        )

    # Regression lines
    for i, s in enumerate(subsets):
        sns.regplot(
            data=s, x=x_col, y=y_col, scatter=False, color=colors[i],
            line_kws={"linewidth": 2, "zorder": 10 + i}, ax=ax,
        )

    # Correlations + legend
    handles, labels, rows = [], [], []
    for i, (p, s) in enumerate(zip(top_percents, subsets)):
        rho_s, p_s = stats.spearmanr(s[x_col], s[y_col])
        rho_p, p_p = stats.pearsonr(s[x_col], s[y_col])
        rows.append({
            "top_percent": p, "n": len(s),
            "spearman_r": rho_s, "spearman_p": p_s,
            "pearson_r": rho_p, "pearson_p": p_p,
        })

        if p >= 100:
            head = f"All predictions (n={len(s)})"
        else:
            head = f"Top {p:g}% by |{ylabel or y_col}| (|y| ≥ {cutoffs[i]:.3g}, n={len(s)})"

        handles.append(Line2D([0], [0], color=colors[i], lw=2))
        labels.append(
            f"{head}\n"
            f"Spearman: {rho_s:.3f} (P={p_s:.1e})\n"
            f"Pearson: {rho_p:.3f} (P={p_p:.1e})"
        )

    ax.legend(
        handles, labels, title="Correlations by Threshold", loc="upper left",
        fontsize=legend_fontsize,
        title_fontproperties={"weight": "bold", "size": legend_title_fontsize},
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel(xlabel if xlabel is not None else x_col)
    ax.set_ylabel(ylabel if ylabel is not None else y_col)
    if title:
        ax.set_title(title)

    # Only manage the figure if we created it
    if created_fig:
        fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    if show:
        plt.show()

    return fig, ax, pd.DataFrame(rows)
#!/usr/bin/env python
# Author: M. Vochteloo, A. Kooijmans

import warnings
warnings.simplefilter("ignore", UserWarning)

import numpy as np
import pandas as pd
import scanpy
import h5py
import argparse
import gzip
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--poolsheet", required=True, type=str, help="")
parser.add_argument("--pool", required=True, type=str, help="")
parser.add_argument("--cell_annotation", required=True, type=str, help="")
parser.add_argument("--rb_genes", required=True, type=str, help="")
parser.add_argument("--mt_genes", required=True, type=str, help="")
parser.add_argument("--ancestry", required=False, type=str, default="EUR", help="")
parser.add_argument("--ncount_rna", required=False, type=int, default=500, help="")
parser.add_argument("--nfeature_rna", required=False, type=int, default=0, help="")
parser.add_argument("--complexity", required=False, type=int, default=100, help="")
parser.add_argument("--percent_rb", required=False, type=int, default=100, help="")
parser.add_argument("--percent_mt", required=False, type=int, default=5, help="")
parser.add_argument("--malat1", required=False, type=int, default=0, help="")
parser.add_argument("--cap_barcodes", required=False, type=int, default=None, help="")
parser.add_argument("--cr_barcodes", action="store_true", default=False, help="")
parser.add_argument("--feature_name", required=False, type=str, default="HGNC", choices=["HGNC", "ENSG", "HGNC_ENSG"], help="")
parser.add_argument("--aggregate_fun", required=False, type=str, default="sum", choices=["sum", "mean"], help="")
parser.add_argument("--save_filtered_h5", action="store_true", default=False, help="")
parser.add_argument("--out", required=True, type=str, help="")

args = parser.parse_args()

os.makedirs(args.out, exist_ok=True)

def gzopen(file, mode="r"):
    if file.endswith(".gz"):
        return gzip.open(file, mode + 't')
    else:
        return open(file, mode)

def load_file(fpath, sep="\t", header=0, index_col=None, must_contain=None):
    data = []
    columns = None
    index = 0
    with gzopen(fpath, mode="r") as fh:
        for index, line in enumerate(fh):
            if index == 0 or (index + 1) % 1e5 == 0:
                print("\tParsed {:,} lines".format(index + 1), end='\r')

            if index != header and must_contain is not None and must_contain not in line:
                continue

            values = line.rstrip("\n").split(sep)
            if index == header:
                columns = values
                continue
            data.append(values)
    print("\tParsed {:,} lines".format(index + 1), end='\n')

    if len(data) == 0:
        print("\tError, failed to load {}: no data".format(os.path.basename(fpath)))
        exit()

    df = pd.DataFrame(data, columns=columns)
    if index_col is not None:
        if len(set(df.columns)) != df.shape[1]:
            print("Error, not all columns are unique.")
            exit()

        column = df.columns[index_col]
        df.index = df[column]
        df.drop(column, axis=1, inplace=True)

    print("\tLoaded dataframe {} with shape: {}".format(os.path.basename(fpath), df.shape))
    return df

def load_metadata(barcodes_fpath):
    print("  Loading Barcodes")
    barcodes = load_file(barcodes_fpath, header=None).rename(columns={0: "Barcode"})
    barcodes.reset_index(drop=False, inplace=True)
    barcodes["Barcode"] = barcodes["Barcode"].str.split("-", n=1, expand=True)[0] + "_" + str(args.pool)
    barcodes["Pool"] = args.pool

    print("  Loading cell annotations")
    cell_annot = load_file(args.cell_annotation, must_contain=args.pool)
    cr_barcodes = None
    if "CellBender" in cell_annot:
        cr_barcodes = [barcode.split("_")[0] + "-1" for barcode in cell_annot.loc[cell_annot["CellBender"] == "False", "Barcode"].values]
    elif "CellRanger" in cell_annot:
        cr_barcodes = [barcode.split("_")[0] + "-1" for barcode in cell_annot.loc[cell_annot["CellRanger"] == "True", "Barcode"].values]
    else:
        print("  Warning, could not add CellRanger column.")

    merge_incl = []
    if "Dataset" in cell_annot.columns:
        print("  Warning, found 'Dataset' column. Will include this as key in merge.")
        merge_incl = ["Dataset"]

    print("\nMerging metadata ...")
    metadata = barcodes.merge(cell_annot, on="Barcode", how="inner")

    # Reorder back into original order.
    metadata = metadata.sort_values(by="index", ascending=True).drop(["index"], axis=1)

    # Check if all barcodes are unique.
    if len(metadata["Barcode"].unique()) != metadata.shape[0]:
        print("Error, not all barcodes are unique in the metadata.")
        exit()

    print("\tLoaded metadata with shape: {}".format(metadata.shape))
    return metadata, cr_barcodes

def load_counts(counts_fpath):
    adata = scanpy.read_10x_h5(counts_fpath)
    count_matrix = adata.X
    n_barcodes, n_genes = count_matrix.shape
    genes = adata.var_names.to_numpy()
    barcodes = adata.obs_names.to_numpy()
    if np.size(genes) != n_genes:
        print("Error, matrix size does not match gene annotation.")
        exit()
    if np.size(barcodes) != n_barcodes:
        print("Error, matrix size does not match barcodes annotation.")
        exit()
    print("\tLoaded matrix with {:,} barcodes and {:,} genes.".format(n_barcodes, n_genes))

    # Check for unique barcodes.
    if len(np.unique(barcodes)) != np.size(barcodes):
        print("Error, not all barcodes are unique.")
        exit()

    # Parse the gene info
    gene_ids = adata.var['gene_ids']
    gene_ids = gene_ids.reset_index(drop=False)
    gene_ids = gene_ids.to_numpy()
    if not (genes == gene_ids[:, 0]).all():
        print("Error, 'count_data.var_names' are expected to have the same order as 'count_data.var['gene_ids'].")
        exit()

    # Set gene names.
    features = get_features(m=gene_ids, indices={"HGNC": 0, "ENSG": 1})

    # Remove duplicates.
    u, c = np.unique(features, return_counts=True)
    dup = u[c > 1]
    features_mask = np.in1d(features, dup)
    if np.sum(features_mask) > 0:
        print("\tRemoving {:,} duplicate genes.".format(np.sum(features_mask)))
        adata = adata[:, ~features_mask]
        adata.var_names = features[~features_mask]

    print("\tLoaded expression with shape: {}".format(count_matrix.shape))
    return adata

def get_features(m, indices):
    if args.feature_name == "HGNC":
        return m[:, indices["HGNC"]]
    elif args.feature_name == "ENSG":
        return m[:, indices["ENSG"]]
    elif args.feature_name == "HGNC_ENSG":
        return np.char.add(np.char.add(m[:, indices["HGNC"]], "_"), m[:, indices["ENSG"]])
    else:
        print("Unexpected feature name '{}'".format(args.feature_name))
        exit()

def load_features(fpath):
    gene_ids = load_file(fpath=fpath).to_numpy()
    features = get_features(m=gene_ids, indices={"HGNC": 1, "ENSG": 0})
    return features

def get_malat1_feature():
    features = get_features(m=np.array([["ENSG00000251562", "MALAT1"]]), indices={"HGNC": 1, "ENSG": 0})
    return features[0]

def calc_barcode_qc(adata, cr_barcodes):
    print("  Adding UMI and feature counts")
    info = pd.DataFrame(np.hstack((np.sum(adata.X, axis=1), np.sum(adata.X != 0, axis=1))),
                        columns=["nCount_RNA", "nFeature_RNA"], index=adata.obs_names)

    print("  Adding complexity")
    info["complexity"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "complexity"] = (info.loc[mask, "nFeature_RNA"] / info.loc[mask, "nCount_RNA"]) * 100
    del mask

    print("  Adding ribosomal %")
    rb_features = load_features(args.rb_genes)
    rb_mask = np.in1d(adata.var_names, rb_features)
    info["rb"] = np.sum(adata.X[:, rb_mask], axis=1)

    info["percent.rb"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "percent.rb"] = (info.loc[mask, "rb"] / info.loc[mask, "nCount_RNA"]) * 100
    del rb_features, rb_mask, mask

    print("  Adding mitochondiral %")
    mt_features = load_features(args.mt_genes)
    mt_mask = np.in1d(adata.var_names, mt_features)
    info["mt"] = np.sum(adata.X[:, mt_mask], axis=1)

    info["percent.mt"] = np.nan
    mask = info["nCount_RNA"] > 0
    info.loc[mask, "percent.mt"] = (info.loc[mask, "mt"] / info.loc[mask, "nCount_RNA"]) * 100
    del mt_features, mt_mask, mask

    print("  Adding MALAT1")
    info["MALAT1"] = 0
    malat1_feature = get_malat1_feature()
    malat1_mask = adata.var_names == malat1_feature
    info["MALAT1"] = adata.X[:, malat1_mask].toarray()
    del malat1_feature, malat1_mask

    print("  Adding nCount_RNA index")
    info["nCount_RNAIndex"] = info[["nCount_RNA"]].rank(method="first", ascending=False)

    print("  Adding CellRanger flag")
    if cr_barcodes is None:
        info["CellRanger"] = True
    else:
        info["CellRanger"] = False
        info.loc[info.index.isin(cr_barcodes), "CellRanger"] = True

    return info

def apply_barcode_qc(barcode_qc):
    print(barcode_qc)

    n_barcodes = barcode_qc.shape[0]
    print("  Input: {:,} barcodes".format(n_barcodes))

    ncount_rna_mask = barcode_qc["nCount_RNA"] >= args.ncount_rna
    print("\tnCount_RNA >= {} yields {:,} [{:.2f}%] barcodes".format(args.ncount_rna, np.sum(ncount_rna_mask), (100 / n_barcodes) * np.sum(ncount_rna_mask)))

    nfeature_rna_mask = barcode_qc["nFeature_RNA"] >= args.nfeature_rna
    print("\tnFeature_RNA >= {} yields {:,} [{:.2f}%] barcodes".format(args.nfeature_rna, np.sum(nfeature_rna_mask), (100 / n_barcodes) * np.sum(nfeature_rna_mask)))

    complexity_mask = barcode_qc["complexity"] <= args.complexity
    print("\tcomplexity <= {} yields {:,} [{:.2f}%] barcodes".format(args.complexity, np.sum(complexity_mask), (100 / n_barcodes) * np.sum(complexity_mask)))

    percent_rb_mask = barcode_qc["percent.rb"] <= args.percent_rb
    print("\tpercent_rb <= {} yields {:,} [{:.2f}%] barcodes".format(args.percent_rb, np.sum(percent_rb_mask), (100 / n_barcodes) * np.sum(percent_rb_mask)))

    percent_mt_mask = barcode_qc["percent.mt"] <= args.percent_mt
    print("\tpercent_mt <= {} yields {:,} [{:.2f}%] barcodes".format(args.percent_mt, np.sum(percent_mt_mask), (100 / n_barcodes) * np.sum(percent_mt_mask)))

    malat1_mask = barcode_qc["MALAT1"] >= args.malat1
    print("\tMALAT1 >= {} yields {:,} [{:.2f}%] barcodes".format(args.malat1, np.sum(malat1_mask), (100 / n_barcodes) * np.sum(malat1_mask)))

    cap_barcodes = np.inf
    if args.cap_barcodes is not None:
        cap_barcodes = args.cap_barcodes
    cap_barcodes_mask = barcode_qc["nCount_RNAIndex"] <= cap_barcodes
    print("\tmax cells <= {} yields {:,} [{:.2f}%] barcodes".format(cap_barcodes, np.sum(cap_barcodes_mask), (100 / n_barcodes) * np.sum(cap_barcodes_mask)))

    cr_barcodes_mask = pd.Series(True, index=barcode_qc.index, name="index")
    if args.cr_barcodes:
        cr_barcodes_mask = barcode_qc["CellRanger"]
    print("\tCellRanger barcodes yields {:,} [{:.2f}%] barcodes".format(np.sum(cr_barcodes_mask), (100 / n_barcodes) * np.sum(cr_barcodes_mask)))

    mask = ncount_rna_mask & nfeature_rna_mask & complexity_mask & percent_rb_mask & percent_mt_mask & malat1_mask & cap_barcodes_mask & cr_barcodes_mask

    print("  Pass-QC: {:,} [{:.2f}%] barcodes".format(np.sum(mask), (100 / n_barcodes) * np.sum(mask)))

    barcode_qc["tag"] = "Outlier"
    barcode_qc.loc[mask, "tag"] = "NotOutlier"

    return barcode_qc

def filter_counts(adata,metadata):
    ncells = adata.X.shape[0]

    # Mask cells based on qc filters
    mask = (metadata["tag"] == "NotOutlier")
    ncells_pass_qc = np.sum(mask) 
    print(f'{ncells_pass_qc} / {ncells} cells passed qc')

    filtered_adata = adata=adata[mask, :]
    filtered_barcodes = filtered_adata.obs_names.to_numpy()
    if args.save_filtered_h5:
        save_filtered_counts_h5(fpath=os.path.join(args.out, f"{args.pool}.filtered.counts.h5"), adata=filtered_adata)
        # Save barcodes
        fout = open(os.path.join(args.out, f"{args.pool}.filtered.barcodes.csv"),'wt')
        for x in filtered_barcodes:
            fout.write(f'{x}\n')
        fout.close()

    return filtered_adata,filtered_barcodes


def save_filtered_counts_h5(fpath, adata):
    """
    Write hdf5 file from Cell Ranger v4 or later versions.
    Source: https://github.com/scverse/anndata/issues/595
    """
    if os.path.exists(fpath):
        raise FileExistsError(f"There already is a file `{fpath}`.")

    int_max = lambda x: int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    str_max = lambda x: max([len(i) for i in x])

    # Save.
    w = h5py.File(fpath, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("feature_type", data=np.array(adata.var.feature_types, dtype=f'|S{str_max(adata.var.feature_types)}'))
    ftrs.create_dataset("genome", data=np.array(adata.var.genome, dtype=f'|S{str_max(adata.var.genome)}'))
    ftrs.create_dataset("id", data=np.array(adata.var.gene_ids, dtype=f'|S{str_max(adata.var.gene_ids)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))

#############################################
print("\nLoading poolsheet ...")
fpaths_df = load_file(args.poolsheet, must_contain=args.pool)
if fpaths_df.shape[0] > 1:
    fpaths_df = fpaths_df.loc[fpaths_df["Pool"] == args.pool, :].reset_index()
fpaths = fpaths_df.to_dict("index")
del fpaths_df

if len(fpaths) != 1:
    print("Error, pool is not unique.")
    exit()
if fpaths[0]["Pool"] != args.pool:
    print("Error, unexpected Pool value.")
    exit()

counts_fpath = fpaths[0]["Counts"]
barcodes_fpath = fpaths[0]["Barcodes"]
print("  Counts fpath: " + counts_fpath)
print("  Barcodes fpath: " + barcodes_fpath)

print("\nCreating metadata ...")
metadata, cr_barcodes = load_metadata(barcodes_fpath=barcodes_fpath)

print("\nLoading counts matrix ...")
adata = load_counts(counts_fpath=counts_fpath)

print("\nCalculating barcode QC stats ...")
barcode_qc = calc_barcode_qc(adata=adata, cr_barcodes=cr_barcodes)

print("\nApplying barcode QC ...")
barcode_qc = apply_barcode_qc(barcode_qc=barcode_qc)
if barcode_qc.shape[0] == 0:
    print("Error, 0 barcodes passed QC.")
    exit()

print("\nMerging metadata and barcode QC ...")
# Add the Barcode_Pool column.
barcode_qc.reset_index(drop=False, inplace=True)
barcode_qc = barcode_qc.rename(columns={"index": "barcode"})
barcode_qc["Barcode"] = barcode_qc["barcode"].str.split("-", n=1, expand=True)[0] + "_" + str(args.pool)

# Merge on Barcode and reorder by the count matrix (barcode column).
metadata = metadata.merge(barcode_qc, on="Barcode", how="inner").set_index("barcode").loc[adata.obs_names, :]

print("\tSaving metadata file")
metadata.to_csv(os.path.join(args.out, str(args.pool) + ".full.metadata.tsv.gz"), sep="\t", header=True, index=False, compression="gzip")

print("\tFiltering counts")
filtered_adata,filtered_barcodes = filter_counts(adata,metadata)

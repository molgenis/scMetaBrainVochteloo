#!/usr/bin/env python
# Author: M. Vochteloo
import argparse
import pandas as pd
import json
import os
import subprocess

"""
Syntax: 
./validate_datasheet.py -h

./validate_datasheet.py --datasheet 2024-10-23-scMetaBrain_datasheet.json
"""



parser = argparse.ArgumentParser(description="")
parser.add_argument("--datasheet", required=True, type=str, help="")
parser.add_argument("--skip_md5sum", action='store_true', help="")
args = parser.parse_args()

print("Options in effect:")
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
print("")


def get_md5sum(fpath):
    output = subprocess.check_output(["md5sum", fpath], stderr=subprocess.STDOUT).decode()
    md5sum, _ = output.split("  ")
    return md5sum


########################################################################

print("Loading datasheet")
df = pd.DataFrame(json.load(open(args.datasheet, "r")))
df.set_index("Dataset", inplace=True)
df.index.name = None
print(df)

columns = df.columns.tolist()
datasets = df.index.tolist()

fhout = None
if not args.skip_md5sum:
    fhout = open(args.datasheet.replace(".json", ".md5sums"), "w")

for dataset in datasets:
    print("Processing dataset '{}'".format(dataset))
    for column in columns:
        print("  column '{}'".format(column))
        fpath = df.loc[dataset, column]
        if os.path.exists(fpath):
            if args.skip_md5sum:
                continue
            md5sum = get_md5sum(fpath=fpath)
            fhout.write(f"{md5sum}  {dataset}_{column}\n")
        else:
            print(f"File '{fpath}' does not exist")
            continue

        if column == "PoolSheet":
            poolsheet_df = pd.read_csv(fpath, sep="\t", header=0, index_col=None)
            for _, row in poolsheet_df.iterrows():
                for poolsheet_column in ["Counts", "Barcodes"]:
                    poolsheet_fpath = row[poolsheet_column]
                    if os.path.exists(poolsheet_fpath):
                        if args.skip_md5sum:
                            continue
                        md5sum = get_md5sum(fpath=poolsheet_fpath)
                        fhout.write(f"{md5sum}  {dataset}_{column}_{poolsheet_column}_{row['Pool']}\n")
                    else:
                        print(f"File '{poolsheet_fpath}' does not exist")
                        continue

if fhout is not None:
    fhout.close()
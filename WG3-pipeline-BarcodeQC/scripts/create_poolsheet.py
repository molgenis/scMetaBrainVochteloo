import json
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--pools", nargs="+", required=True, type=str, help="")
parser.add_argument("--outdir", required=True, type=str, help="")
parser.add_argument("--dataset", required=True, type=str, help="")
parser.add_argument("--outfile", required=True, type=str, help="")

args = parser.parse_args()
print("Options in effect:")
for arg in vars(args):
	print("  --{} {}".format(arg, getattr(args, arg)))
print("")

print("Creating poolsheet for {len(args.pools)} pools")
print("Writing to:")
print(args.outfile)
fout = open(args.outfile,'wt')
fout.write('Pool\tCounts\tBarcodes\tBam\tDataset\n')
for pool in args.pools:
	counts = f'{args.outdir}/{pool}/{pool}.filtered.counts.h5'
	barcodes = f'{args.outdir}/{pool}/{pool}.filtered.barcodes.csv'
	bam = ""
	fout.write(f"{pool}\t{counts}\t{barcodes}\t{bam}\t{args.dataset}\n")
fout.close()

print("Done.")

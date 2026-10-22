#!/usr/bin/env Rscript
# Adapted from drneavin https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/blob/master/Demultiplexing/scripts/Singlet_QC_Figures.R
.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-p", "--poolsheet", required=TRUE, help="")
parser$add_argument("-s", "--seurat", required=TRUE, help="")
parser$add_argument("-o", "--out", required=TRUE, help="The output directory where results will be saved.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

## make sure the directory exists ###
dir.create(args$out, recursive = TRUE)

print("Options in effect:")
for (name in names(args)) {
	print(paste0("  --", name, " ", args[[name]]))
}
print("")

# suppressMessages(suppressWarnings(library(scCustomize)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggforce)))
suppressMessages(suppressWarnings(library(ggnewscale)))

seurat <- readRDS(args$seurat)


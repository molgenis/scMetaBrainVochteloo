#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Merge seurat object
# author: Martijn Vochteloo
# date: 2023-12-14
# description:

#   ____________________________________________________________________________
#   Import libraries                                                        ####

.libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(progressr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(scCustomize)))

#   ____________________________________________________________________________
#   Set up parameter variables                                              ####

option_list <-  list(
  make_option("--seurat",
              type = "character",
              default = NULL,
              help = crayon::green(""),
              metavar = "character"),
  make_option("--out",
              type = "character",
              default = "split_object",
              help = crayon::green("Output file name [default= %default]"),
              metavar = "character"),
  make_option("--path",
              type = "character",
              default = ".",
              help = crayon::green("Output path to store results [default= %default]"),
              metavar = "character")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print("Options in effect:")
for (name in names(opt)) {
	print(paste0("  --", name, " ", opt[[name]]))
}
print("")

#   ____________________________________________________________________________
#   Helper functions

seurat_objects <- apply(opt$seurat, loadRDS)
seurat <- merge(seurat_objects[[1]], y = seurat_objects[2:length(seurat_objects)])
saveRDS(seurat, file.path(opt$path, paste0(opt$out, ".RDS")))
#!/usr/bin/env Rscript

library(coloc)
library(susieR)

args <- commandArgs(trailingOnly=TRUE)
gene <- args[1]
relFinemapPath <- args[2]
cellType <- args[3]
gwas <- args[4]

relFinemap = data.table::fread(paste0(relFinemapPath),data.table=F)

relFinemapSel <- relFinemap

##Get LBF columns for the GWAS
matches <- grep("lbf_cs_.*_gwas", colnames(relFinemapSel))

# Quit if no credible sets are found for the GWAS
if (length(matches) == 0) {
    message("No columns matching 'lbf_cs_.*_gwas' found. Exiting.")
    quit(save = "no", status = 0)
}

# Grep GWAS LBF columns
lbfCols <- grep(
    pattern = "^lbf_cs_[0-9]+_gwas$",
    x = colnames(relFinemapSel),
    value = TRUE
)

# Extract credible set number
csNumbers <- as.integer(
    sub("^lbf_cs_([0-9]+)_gwas$", "\\1", lbfCols)
)

# Sort columns by credible set number
lbfCols <- lbfCols[order(csNumbers)]

# Create bfSetR
bfSetR <- t(relFinemapSel[, lbfCols, drop = FALSE])
bfSetR <- as.matrix(bfSetR)

colnames(bfSetR) <- relFinemapSel$snp_formatted_sorted

print(str(bfSetR))
print(which(!is.na(unique(relFinemapSel$SusieRss_CS_gwas))))

bfSetR = as.matrix(bfSetR[,,drop=F])
colnames(bfSetR) = relFinemapSel$snp_formatted_sorted
print(str(bfSetR))


# Get LBF columns for the eQTLs
lbfCols <- grep(
    pattern = "^lbf_cs_[0-9]+_eqtl$",
    x = colnames(relFinemapSel),
    value = TRUE
)

# Extract credible set number
csNumbers <- as.integer(
    sub("^lbf_cs_([0-9]+)_eqtl$", "\\1", lbfCols)
)

# Sort columns by credible set number
lbfCols <- lbfCols[order(csNumbers)]

# Create bfSetR
bfSetO <- t(relFinemapSel[, lbfCols, drop = FALSE])
bfSetO <- as.matrix(bfSetO)

colnames(bfSetO) <- relFinemapSel$snp_formatted_sorted

# Run coloc
colocOut = coloc.bf_bf(bfSetR,bfSetO)
colocOutSummary = as.data.frame(colocOut$summary)
print(colocOutSummary)


# Fix the names of the relevant sets LBFS
for(r in 1:nrow(bfSetR)){
    colocOutSummary$idx1[which(colocOutSummary$gwas==r)] = rownames(bfSetR)[r]
    print(str(colocOutSummary))
}

# Fix the names of the other sets LBFS
for(r in 1:nrow(bfSetO)){
    colocOutSummary$idx2[which(colocOutSummary$eqtl==r)] = rownames(bfSetO)[r]
}

outMatrix = list()
outMatrix[[paste(gene)]] = colocOutSummary
outMatrix = as.data.frame(data.table::rbindlist(outMatrix,fill = T))

colnames(outMatrix)[colnames(outMatrix) == "hit1"] <- "sorted_snp_gwas"
colnames(outMatrix)[colnames(outMatrix) == "hit2"] <- "sorted_snp_eqtl"



print(outMatrix)

outMatrix$gwas_lambda <- relFinemap$SusieRss_lambda_gwas[1]
outMatrix$eqtl_lambda <- relFinemap$SusieRss_lambda_eqtl[1]
outMatrix$cell_type <- cellType
write.table(outMatrix,file = paste0(gwas, '_', cellType, '_', gene,"_mv_coloc_out.txt") ,quote=F, sep="\t",row.names=F)



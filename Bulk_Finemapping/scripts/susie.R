#!/usr/bin/env Rscript

library(susieR)
set.seed(1)

# Parse command line arguments
args <- commandArgs(trailingOnly=TRUE)
locus <- args[1]
data_path <- args[2]
ld_data_path <- args[3]
eqtl_n <- as.integer(args[4])
num_cs <- as.integer(args[5])

filterNomSig = NA

# Load data
df <- read.csv(data_path, sep = "\t")
print(str(df))

# Load LD matrix
ld_df <- read.table(ld_data_path, header = TRUE, sep = "\t", row.names = 1)
ld_matrix <- as.matrix(ld_df)
print(str(ld_matrix))

# Ensure order of snps between eqtls and LD matrix is the same
df <- df[df$variant_id_formatted %in% rownames(ld_df), ]
df <- df[match(rownames(ld_df), df$variant_id_formatted), ]
df <- df[!is.na(df$variant_id_formatted), ]


estimated_res_var = T

# Calculate lambda
df["SusieRss_lambda"] = estimate_s_rss(z=df[['MetaPZ']], R = as.matrix(ld_matrix),n=eqtl_n)
print(str(df))

# Create kriging plot
condz = kriging_rss(df[['MetaPZ']], as.matrix(ld_matrix), n=eqtl_n)
options(bitmapType = "cairo")
png(paste0(locus, "_zscore_comparison.png"), width = 800, height = 600)
print(condz$plot)
dev.off()

print(str(df))

fitted_rss2 = tryCatch({
    fitted_rss2 = susie_rss(z=df[['MetaPZ']], R = as.matrix(ld_matrix), n = eqtl_n, L = num_cs, estimate_residual_variance = TRUE)
},error = function(e){
    fitted_rss2 = susie_rss(z=df[['MetaPZ']], R = as.matrix(ld_matrix), n = eqtl_n, L = num_cs, estimate_residual_variance = FALSE)
    estimated_res_var = F
    return(fitted_rss2)
})

if(!fitted_rss2$converged){
    fitted_rss2 = susie_rss(z=df[['MetaPZ']], R = as.matrix(ld_matrix), n = eqtl_n, L = num_cs, estimate_residual_variance = FALSE)
    estimated_res_var = F
}

print = T
direct = T
formatted = T

if(formatted){
    outMat = NULL
    }

if(direct){
    print(summary(fitted_rss2))
}

df['L'] = num_cs

if(formatted){
    df["SusieRss_pip"] = fitted_rss2$pip 
    df["SusieRss_CS"] = NA
    df["SusieRss_ResVar"] = estimated_res_var

    if(length(fitted_rss2$sets[[1]])!=0){
        for(l in 1:length(fitted_rss2$sets[[1]])){
            row_index = as.numeric(fitted_rss2$sets[[1]][[l]])  # Ensure index is numeric
            df[["SusieRss_CS"]][row_index] = gsub("L", "", names(fitted_rss2$sets[[1]])[l])
        }
    }

    lbfOut = t(fitted_rss2$lbf_variable)
    if(ncol(lbfOut)<num_cs){
        lbfOut = as.data.frame(lbfOut)
        for(j in 1:(num_cs - ncol(lbfOut))){
            lbfOut[paste("lbf_cs",j,sep="_")] = NA
        }
    }
    colnames(lbfOut) = paste("lbf_cs",1:num_cs,sep="_")
    df = cbind(df,lbfOut)
    outMat = rbind(outMat,df)
}

# Write to file
write.table(outMat, paste(locus, "susie.txt", sep="_"), quote=F, sep="\t", row.names=F)

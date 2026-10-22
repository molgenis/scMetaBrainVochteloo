#!/usr/bin/env Rscript

library(susieR)
library(RcppCNPy)
set.seed(1)

# Parse command line arguments
args <- commandArgs(trailingOnly=TRUE)
locus <- args[1]
data_path <- args[2]
ld_data_path <- args[3]
gwas_n <- as.integer(args[4])
num_cs <- as.integer(args[5])

filterNomSig = NA
estimated_res_var = T

# Load locus file
df <- read.csv(data_path, sep = "\t")
print(str(df))

# Load LD matrix from numpy file
ld_matrix <- npyLoad(ld_data_path)

# Calculate z-scores for lambda calculation
df$zscore <- df$beta / df$standard_error

# Calculate predicted z-scores
condz = kriging_rss(df$zscore, as.matrix(ld_matrix), n=gwas_n)

# Add predicted z-score column to dataframe
df$z_pred = condz$conditional_dist$condmean

# Create column with absolute differences between observed z-score and predicted z-score
df$z_diff <- abs(df$zscore - df$z_pred)

# Calculate z-scores over zscores differences
df$z_diff_zscore <- (df$z_diff - mean(df$z_diff, na.rm = TRUE)) / sd(df$z_diff, na.rm = TRUE)

# Get indices of outlier variants based on above calculated z-scores
outlier_idx <- which(df$z_diff_zscore > 3)

# Remove outlier variants from LD matrix
ld_matrix <- ld_matrix[-outlier_idx, -outlier_idx]

# Remove outlier variants from dataframe
df <- df[-outlier_idx, ]

# Calculate lambda values after removal of outliers
lambda = estimate_s_rss(z=df$zscore, R = as.matrix(ld_matrix),n=gwas_n)
df["SusieRss_lambda"] = lambda

# Run SuSiE
fitted_rss2 = tryCatch({
    fitted_rss2 = susie_rss(z=df$zscore, R = as.matrix(ld_matrix), n = gwas_n, L = num_cs, estimate_residual_variance = TRUE)
},error = function(e){
    fitted_rss2 = susie_rss(z=df$zscore, R = as.matrix(ld_matrix), n = gwas_n, L = num_cs, estimate_residual_variance = FALSE)
    estimated_res_var = F
    return(fitted_rss2)
})

if(!fitted_rss2$converged){
    fitted_rss2 = susie_rss(z=df$zscore, R = as.matrix(ld_matrix), n = gwas_n, L = num_cs, estimate_residual_variance = FALSE)
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


if(formatted){
    df["SusieRss_pip"] = fitted_rss2$pip 
    df["SusieRss_CS"] = NA
    df["SusieRss_ResVar"] = estimated_res_var
    df['L'] = num_cs

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


write.table(outMat, paste0(locus, "_susie.txt"), quote=F, sep="\t", row.names=F)


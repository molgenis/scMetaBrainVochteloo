# Based on the mashr eQTL analysis outline vignette:
# https://stephenslab.github.io/mashr/articles/eQTL_outline.html

library(ashr)
library(mashr)
library(corrplot)

# Load strong effects
Bhat_strong <- as.matrix(read.table("../1-data/strong/strong_betas.txt", header = TRUE, row.names = 1, fill = TRUE))
Shat_strong <- as.matrix(read.table("../1-data/strong/strong_standard_errors.txt", header = TRUE, row.names = 1, fill = TRUE))

# Load random effects
Bhat_random <- as.matrix(read.table("../1-data/random/random_betas.txt", header = TRUE, row.names = 1, fill = TRUE))
Shat_random <- as.matrix(read.table("../1-data/random/random_standard_errors.txt", header = TRUE, row.names = 1, fill = TRUE))

# Estimate correlation structure
data.temp = mash_set_data(Bhat_random,Shat_random)
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

# Set up data driven covariances
data.random = mash_set_data(Bhat_random,Shat_random,V=Vhat)
data.strong = mash_set_data(Bhat_strong,Shat_strong, V=Vhat)

# PCA
U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)

# Fit mash model
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

# Compute posterior summaries
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

str(m2)

# Save mashr object to file
saveRDS(m2, file = "m2.rds")

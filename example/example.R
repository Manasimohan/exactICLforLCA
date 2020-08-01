library(exactICLforLCA)
library(BayesLCA)

# Read data
data(Alzheimer)

## set CONSTANTS
alpha_var <- 1
beta_var <- 1
# G <- as.integer(readline(prompt="Enter the number of groups: "))
G <- 3

# LCA fit
fit <- blca.em(Alzheimer, G)
Z <- unMAP(MAP(Zscore(Alzheimer, fit)))

Y <- as.matrix(Alzheimer)

res_final <- ICLFit(Z, Y, G, alpha_var, beta_var)
ICL_val <- res_final$ICL_val_max
Z <- res_final$Z_max

print(ICL_val)
print(ncol(Z))

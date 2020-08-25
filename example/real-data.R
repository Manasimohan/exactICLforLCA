library(exactICLforLCA)
library(BayesLCA)

# Read data
data(Alzheimer)

## set CONSTANTS
alpha_var <- 1
beta_var <- 1
delta_var <- 1
# G <- as.integer(readline(prompt="Enter the number of groups: "))
G <- 3

# LCA fit
fit <- blca.em(Alzheimer, G)
Z <- unMAP(MAP(Zscore(Alzheimer, fit)))

Y <- as.matrix(Alzheimer)

res <- ICLFit(Z, Y, G, alpha_var, beta_var, delta_var)
ICLsummary()

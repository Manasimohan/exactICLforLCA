library(exactICLforLCA)
library(BayesLCA)

# Read data
data(Alzheimer)

## set CONSTANTS
alpha_var <- 1
beta_var <- 1
# G <- as.integer(readline(prompt="Enter the number of groups: "))
G <- 9

# LCA fit
fit <- blca.em(Alzheimer, G)
Z <- unMAP(MAP(Zscore(Alzheimer, fit)))

Y <- as.matrix(Alzheimer)

# print(ICLGOne(Z, Y, G, alpha_var, beta_var))

res <- ICLFit(Z, Y, G, alpha_var, beta_var)
ICLsummary()

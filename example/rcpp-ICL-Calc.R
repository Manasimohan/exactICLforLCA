library(BayesLCA)
library(Rcpp)
sourceCpp("Rcpp\\ICLCalc.cpp")

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

res_final <- ICLCalc(alpha_var, beta_var, G, Y, Z, delta_var)

print(res_final)

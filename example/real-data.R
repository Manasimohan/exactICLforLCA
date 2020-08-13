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
ICL_original <- res_final$ICL_original

#print(ICL_val)
#print(ncol(Z))
#print(ICL_original)


print(paste('Oringinal ICL value : ', format(ICL_original, digits = 5)))
print(paste("Incial number of clusters : ", G))
print(paste("ICL value post processing : ", format(ICL_val, digits = 5)))
print(paste("Number of clusters post processing : ", ncol(Z)))


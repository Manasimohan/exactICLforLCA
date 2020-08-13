library(BayesLCA)
library(exactICLforLCA)

# Generate data
type1 <- c(0.9, 0.3, 0.4, 0.2, 0.6)
type2 <- c(0.1, 0.7, 0.6, 0.8, 0.4)
X<- rlca(500, rbind(type1, type2), c(0.7, 0.3))

## set CONSTANTS
alpha_var <- 1
beta_var <- 1
# G <- as.integer(readline(prompt="Enter the number of groups: "))
G <- 3

# LCA fit
fit <- blca.em(X, G)
Z <- unMAP(MAP(Zscore(X, fit)))


Y <- as.matrix(X)

res_final <- ICLFit(Z, Y, G, alpha_var, beta_var)
ICL_val <- res_final$ICL_val_max
Z <- res_final$Z_max
ICL_original <- res_final$ICL_original

print(paste('Oringinal ICL value : ', format(ICL_original, digits = 5)))
print(paste("Incial number of clusters : ", G))
print(paste("ICL value post appling sweep and merge : ", format(ICL_val, digits = 5)))
print(paste("Number of clusters post processing : ", ncol(Z)))


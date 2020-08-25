library(BayesLCA)
library(exactICLforLCA)
library(flexclust)

rlca_fun <- function (n, itemprob = 0.5, classprob = 1, fit = NULL)
{
  if (is.null(fit)) {
    itemprob <- as.matrix(itemprob)
    G <- nrow(itemprob)
    M <- ncol(itemprob)
  }
  else {
    itemprob <- fit$itemprob
    classprob <- fit$classprob
    G <- nrow(itemprob)
    M <- ncol(itemprob)
  }
  x <- matrix(runif(n * M), nrow = n)
  z <- rep(0, n)
  classvec <- as.vector(rmultinom(1, n, prob = classprob))
  ind <- c(0, cumsum(classvec))
  for (g in 1:G){
    x[(ind[g] + 1):ind[g + 1], ] <- t(t(x[(ind[g] + 1):ind[g + 1], ]) < itemprob[g, ]) * 1
    z[(ind[g] + 1):ind[g + 1]] <- g
  }
  return(list(x=x, z=z))
}

# set CONSTANTS
alpha_var <- 1
beta_var <- 1
delta_var <- 1
G <- 3
n <- 1000

# Generate data
type1 <- c(0.9, 0.3, 0.7,0.1)
type2 <- c(0.1, 0.7, 0.3,0.9)
sim_data_list <- rlca_fun(n, rbind(type1,type2), c(0.6,0.4))

# splitting data into x and z
x_sim <- sim_data_list$x
z_sim <- sim_data_list$z

# splitting X based on groups
x_1 <- x_sim[z_sim == 1, ]
x_2 <- x_sim[z_sim == 2, ]

# finding means of X groups
#apply(x_1, 2, mean)
#apply(x_2, 2, mean)
#table(z_sim)/length(z_sim)

# LCA fit
fit_new <- blca.em(x_sim, G )
Z_new <- MAP(Zscore(x_sim, fit_new))

# splitting Z based on groups
#X_1 <- x_sim[Z_new == 1, ]
#X_2 <- x_sim[Z_new == 2, ]

# finding means of Z groups
#apply(X_1, 2, mean)
#apply(X_2, 2, mean)
#table(Z_new)/length(Z_new)

# ICL computation
Y <- as.matrix(x_sim)
Z <- unMAP(Z_new)

#Outputs

#Obtaing ICL for Z matrix
res_final <- ICLFit(Z, Y, G, alpha_var, beta_var, delta_var)
ICLsummary()

Z_new <- MAP(res_final$Z_max)

#Comparing the two matrix
table(Z_new, z_sim)
randIndex(Z_new, z_sim)

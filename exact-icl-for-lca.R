install.packages('BayesLCA')
library(BayesLCA)

data(Alzheimer)

Alzheimer

fit <- blca.em(Alzheimer, 2)

Z <- unMAP(MAP(Zscore(Alzheimer, fit)))

Y <- as.matrix(Alzheimer)

##Variables

delta <- c(1, 1)

delta_prime <- delta + apply(Z, 2, sum)

alpha <- 1

beta <- 1

G <- 2

r <- ncol(Y)

alpha_gj <- beta_gj <- matrix(0, nrow = G, ncol = r)

##Function Call

log_beta_vec <- function(delta){
  
  l_num <- sum(lgamma(delta))
  
  l_denom <- lgamma(sum(delta))
  
  l_num - l_denom
  
}


for(j in 1:r){
  
  alpha_gj[, j] <- alpha + apply(Y[, j] * Z, 2, sum)
  
  beta_gj[, j] <- beta + apply( (1 - Y[, j]) * Z, 2, sum)
  
}


b_num <- 0

for (g in 1:G)
{
  for(j in 1:r)
  {
    
    b_num = b_num + lbeta(alpha_gj[g,j], beta_gj[g,j])
  }
}

b_num

b_denom <- G * r * lbeta(alpha, beta)

ICL <- ( log_beta_vec(delta_prime) - log_beta_vec(delta) ) + (  b_num - b_denom )

ICL

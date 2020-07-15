install.packages('BayesLCA')
library(BayesLCA)

data(Alzheimer)

Y <- as.matrix(Alzheimer)

r <- ncol(Y)

ICL_calc <- function() {
  
  alpha <- as.integer(readline(prompt="Enter alpha value: "))
  
  beta <- as.integer(readline(prompt="Enter beta value: "))
  
  G <- as.integer(readline(prompt="Enter the number of groups: "))
  
  delta <- matrix(0, nrow=1 , ncol = G)
  
  readinteger <- function()
  { 
    for(i in 1:G)
    {
      delta[1,i]<- as.integer(readline(prompt="Enter integer for delta: "))
    }
    delta
  }
  delta<-readinteger ()
  return(list("alpha"=alpha, "beta"=beta, "G"=G, "delta"=delta))
}

inputs <- ICL_calc()

alpha <- inputs$alpha
beta <- inputs$beta
G <- inputs$G
delta <- inputs$delta

fit <- blca.em(Alzheimer, G)

Z <- unMAP(MAP(Zscore(Alzheimer, fit)))

delta_prime <- delta + apply(Z, 2, sum)

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

ICL <- ( log_beta_vec(delta_prime) - log_beta_vec(delta) ) + (  b_num - (G * r * lbeta(alpha, beta)) )

ICL


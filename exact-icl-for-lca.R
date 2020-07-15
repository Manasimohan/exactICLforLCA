library(BayesLCA)

# Read data
data(Alzheimer)

## set CONSTANTS 
alpha <- 1
beta <- 1
G <- as.integer(readline(prompt="Enter the number of groups: "))

# LCA fit
fit <- blca.em(Alzheimer, G)
Z <- unMAP(MAP(Zscore(Alzheimer, fit)))

Y <- as.matrix(Alzheimer)

icl_calc_func <- function(alpha, beta, G, Y, Z) {
  # set variables
  delta <- rep(1 , G)
  r <- ncol(Y)
  
  # initalise alpha and beta matrix
  alpha_gj <- beta_gj <- matrix(0, nrow = G, ncol = r)
  
  # initalise delta matrix
  delta_prime <- delta + apply(Z, 2, sum)
  
  # calc alpha and beta of groups based on Y and Z
  for(j in 1:r){
    alpha_gj[, j] <- alpha + apply(Y[, j] * Z, 2, sum)
    beta_gj[, j] <- beta + apply( (1 - Y[, j]) * Z, 2, sum)
  }
  
  # log beta fun
  log_beta_vec <- function(delta) {
    l_num <- sum(lgamma(delta))
    l_denom <- lgamma(sum(delta))
    l_num - l_denom
  }
  
  # first eqn
  first_var <- log_beta_vec(delta_prime) - log_beta_vec(delta)
  
  # second eqn numerator value
  b_num <- 0
  for (g in 1:G)
  {
    for(j in 1:r)
    {
      b_num = b_num + lbeta(alpha_gj[g,j], beta_gj[g,j])
    }
  }
  
  # second eqn denomenarator value  
  b_denom <- G * r * lbeta(alpha, beta)
  
  # second eqn
  sec_var <- b_num - b_denom 
  
  # ICL calc
  ICL <- first_var + sec_var
  
  ICL
}

update_Z_of_obs_i <- function(Z, G, i, t = 2) {
  if(G == 2) {
    for (j in 1:G) {
      if(Z[i,j] == 0) {
        Z[i,j] <- 1
      } else {
        Z[i,j] <- 0
      }
    }
  } else {
    for(j in 1:G) {
      if(Z[i,j] == 1) {
        Z[i,j] <- 0
        k <- c(j)
      }
    }
    if (t == k) {
      Z[i, k] <- 1 
    } else {
      Z[i, t] <- 1 
    }
  }
  Z
}

icl_calc_func(alpha, beta, G, Y, Z)

for(i in  1 : nrow(Z)) {
  
  if (G == 2) {
    ICL_val1 <- icl_calc_func(alpha, beta, G, Y, Z)
    Z <- update_Z_of_obs_i(Z, G, i)
    ICL_val2 <- icl_calc_func(alpha, beta, G, Y, Z)
    if(ICL_val2 - ICL_val1 > 0) {
      # DO nothing keep the changed Z
    } else {
      # revert it back to original
      Z <- update_Z_of_obs_i(Z, G, i)
    }
  } else {
    # original ICL value without any changes
    ICL_val <- icl_calc_func(alpha, beta, G, Y, Z)
    # to find which value in the ith observation has 1
    # and store that in g
    for(j in 1:G) {
      if(Z[i,j] == 1) {
        g <- c(j)
      }
    }
    h_vals <- c(1:G)
    h_vals <- setdiff(h_vals, g)
    
    ICL_max <- ICL_val
    ICL_h <- g
    for(h in h_vals) {
      # changing the cluster of ith obervation from group g to group h
      Z <- update_Z_of_obs_i(Z, G, i, h)
      # calculating ICL value of the new combination
      ICL_val_of_h <- icl_calc_func(alpha, beta, G, Y, Z)
      # reverting back to original combination
      Z <- update_Z_of_obs_i(Z, G, i, g)
      
      ICL_del <- ICL_val_of_h - ICL_max
      
      if(ICL_del > 0) {
        ICL_max <- ICL_val_of_h
        ICL_h <- h
      }
    }
    # changing to the combination with highest ICL value
    Z <- update_Z_of_obs_i(Z, G, i, ICL_h)
  }
}

ICL_val_new <- icl_calc_func(alpha, beta, G, Y, Z)
ICL_val_new
View(Z)

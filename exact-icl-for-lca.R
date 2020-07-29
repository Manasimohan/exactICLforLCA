library(BayesLCA)

# Read data
data(Alzheimer)

## set CONSTANTS 
alpha <- 1
beta <- 1
# G <- as.integer(readline(prompt="Enter the number of groups: "))
G <- 9

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
    for(g in 1:G) {
      temp_alpha <- 0
      temp_beta <- 0
      for(i in 1:nrow(Z)) {
        temp_alpha <- temp_alpha + Z[i, g] * Y[i, j]
        temp_beta <- temp_beta + Z[i, g] * (1- Y[i, j])
      }
      alpha_gj[g, j] <- alpha + temp_alpha
      beta_gj[g, j] <- beta + temp_beta
    }
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

check_group_reduction <- function(Z, G) {
  group_red <- FALSE
  group_index <- -1
  for ( i in 1 : G) {
    if (length(Z[, i][Z[, i] == TRUE]) == 0) {
      group_red <- TRUE
      group_index <- i
      break
    }
  }
  if(group_red) {
    Z <- Z[ , -c(group_index)]
  }
  Z
}

group_reduction <- function(samp_df) {
  for(i in samp_df) {
    
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
        # check if the groups have reduced and del Z[, col] of that group
        Z1 <- check_group_reduction(Z, G)
        
        # if the group has reduced reduce G val
        if(ncol(Z) != ncol(Z1)) {
          G <- G - 1
        }
        
        # calculating ICL value of the new combination
        ICL_val_of_h <- icl_calc_func(alpha, beta, G, Y, Z1)
        
        # reverting back to original combination
        if(ncol(Z) != ncol(Z1)) {
          G <- G + 1
        } else {
          Z <- update_Z_of_obs_i(Z, G, i, g)
        }
        
        ICL_del <- ICL_val_of_h - ICL_max
        
        if(ICL_del > 0) {
          ICL_max <- ICL_val_of_h
          ICL_h <- h
        }
      }
      group_reduced <- FALSE
      # if the group has reduced reduce G val
      if(ncol(Z) != ncol(Z1)) {
        group_reduced <- TRUE
      }
      # changing to the combination with highest ICL value
      Z <- update_Z_of_obs_i(Z, G, i, ICL_h)
      Z <- check_group_reduction(Z, G)
      
      if(group_reduced) {
        G <- G - 1
        i <- 0
      }
    }
    #print(i)
    #print(G)
  }
  return(list("G"=G, "Z"=Z))
}

ICL_fit <- function(Z, Y, alpha, beta) {
  ICL_old <- icl_calc_func(alpha, beta, G, Y, Z)
  print(ICL_old)
  
  ICL_val_new <- 0
  
  iter <- 0
  while(TRUE){
    samp_ind <- sample(1:nrow(Z))
    
    results<-group_reduction(samp_ind)
    
    G <- results$G
    Z <- results$Z
    
    if(ICL_val_new != 0) {
      ICL_old <- ICL_val_new
    }
    ICL_val_new <- icl_calc_func(alpha, beta, G, Y, Z)
    
    print(ICL_val_new)
    
    del <- ICL_val_new - ICL_old
    
    #del > 0
    
    iter <- iter + 1
    print(del)
    if(del <= 0){
      break
    }
    if(iter==10){
      break
    }
  }
  return(list("G"=G, "Z"=Z, "ICL_val_new"=ICL_val_new))
}

check_group_merge <- function(Z, Y, g1, g2, alpha, beta) {
  
  Z[ , g1] <- Z[ , g1] + Z[ , g2]
  Z <- Z[ , -c(g2)]
  
  ICL_val_merge <- icl_calc_func(alpha, beta, ncol(Z), Y, Z)
  
  return(list("G"=G, "Z"=Z, "ICL_val_merge"=ICL_val_merge))
}

ICL_group_merge <- function(ICL_val_old, Z, Y, G, alpha, beta) {
  ICL_max <- ICL_val_old
  Z_max <- 0
  g1 <- 0
  g2 <- 0
  
  for(i in 1:(G-1)) {
    for(j in (i+1):G) {
      results <- check_group_merge(Z, Y, i, j, alpha, beta)
      
      G <- results$G
      Z1 <- results$Z
      ICL_val_merge <- results$ICL_val_merge
      
      
      ICL_del <- ICL_val_merge - ICL_max
      
      print(ICL_del)
      
      if(ICL_del > 0) {
        ICL_max <- ICL_val_merge
        g1 <- i
        g2 <- j
        Z_max <- Z1
      }
    }
  }
  return(list("Z_max"=Z_max, "ICL_max"=ICL_max))
}


results<-ICL_fit(Z, Y, alpha, beta)

G <- results$G
Z <- results$Z
ICL_val_old <- results$ICL_val_new

results <- ICL_group_merge(ICL_val_old, Z, Y, G, alpha, beta)

ICL_max <- results$ICL_max
Z_max <- results$Z_max

print(ICL_max)
print(ncol(Z_max))


#Z#View(Z)
#length(Z[, 4][Z[, 4] == TRUE])
G
#ncol(Z)



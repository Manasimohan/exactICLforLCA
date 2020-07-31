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

ICLCalc <- function(alpha_var, beta_var, G, Y, Z) {
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
      alpha_gj[g, j] <- alpha_var + temp_alpha
      beta_gj[g, j] <- beta_var + temp_beta
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
  b_denom <- G * r * lbeta(alpha_var, beta_var)

  # second eqn
  sec_var <- b_num - b_denom

  # ICL calc
  ICL <- first_var + sec_var

  ICL
}

updateZ <- function(Z, G, i, t = 2) {
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

checkGroupReduction <- function(Z, G) {
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

groupReduction <- function(samp_df, Z, G) {
  for(i in samp_df) {

    if (G == 2) {
      ICL_val1 <- ICLCalc(alpha_var, beta_var, G, Y, Z)
      Z <- updateZ(Z, G, i)

      ICL_val2 <- ICLCalc(alpha_var, beta_var, G, Y, Z)
      if(ICL_val2 - ICL_val1 > 0) {
        # DO nothing keep the changed Z
      } else {
        # revert it back to original
        Z <- updateZ(Z, G, i)
      }
    } else {
      # original ICL value without any changes
      ICL_val <- ICLCalc(alpha_var, beta_var, G, Y, Z)
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
        Z <- updateZ(Z, G, i, h)
        # check if the groups have reduced and del Z[, col] of that group
        Z1 <- checkGroupReduction(Z, G)

        # if the group has reduced reduce G val
        if(ncol(Z) > ncol(Z1)) {
          G <- G - 1
        }

        # calculating ICL value of the new combination
        ICL_val_of_h <- ICLCalc(alpha_var, beta_var, G, Y, Z1)

        # reverting back to original combination
        if(ncol(Z) > ncol(Z1)) {
          G <- G + 1
        } else {
          Z <- updateZ(Z, G, i, g)
        }

        ICL_del <- ICL_val_of_h - ICL_max

        if(ICL_del > 0) {
          ICL_max <- ICL_val_of_h
          ICL_h <- h
        }
      }
      group_reduced <- FALSE
      # if the group has reduced reduce G val
      if(ncol(Z) > ncol(Z1)) {
        group_reduced <- TRUE
      }
      # changing to the combination with highest ICL value
      Z <- updateZ(Z, G, i, ICL_h)
      Z <- checkGroupReduction(Z, G)

      if(group_reduced) {
        G <- G - 1
      }
    }
  }
  return(list("G"=G, "Z"=Z))
}

ICLSweep <- function(Z, Y, G, alpha_var, beta_var) {
  ICL_old <- ICLCalc(alpha_var, beta_var, G, Y, Z)

  ICL_val_new <- 0
  Z_max <- Z
  G_max <- G

  iter <- 0
  while(TRUE){
    samp_ind <- sample(1:nrow(Z))

    results<-groupReduction(samp_ind, Z, G)

    G <- results$G
    Z <- results$Z

    if(ICL_val_new != 0) {
      ICL_old <- ICL_val_new
    }
    ICL_val_new <- ICLCalc(alpha_var, beta_var, G, Y, Z)

    del <- ICL_val_new - ICL_old

    iter <- iter + 1
    if(del <= 0){
      break
    }
    if(iter==10){
      break
    }
    Z_max <- Z
    G_max <- G
  }
  return(list("G_max"=G_max, "Z_max"=Z_max, "ICL_old"=ICL_old))
}

checkGroupMerge <- function(Z, Y, g1, g2, alpha_var, beta_var) {

  Z[ , g1] <- Z[ , g1] + Z[ , g2]
  Z <- Z[ , -c(g2)]

  ICL_val_merge <- ICLCalc(alpha_var, beta_var, ncol(Z), Y, Z)

  return(list("Z"=Z, "ICL_val_merge"=ICL_val_merge))
}

ICLGroupMerge <- function(ICL_val, Z, Y, G, alpha_var, beta_var) {
  ICL_max <- ICL_val
  Z_max <- 0
  g1 <- 0
  g2 <- 0

  if(G == 2) {
    Z_max <- Z
    return(list("Z_max"=Z_max, "ICL_max"=ICL_max))
  }

  for(i in 1:(G-1)) {
    for(j in (i+1):G) {
      results <- checkGroupMerge(Z, Y, i, j, alpha_var, beta_var)

      Z1 <- results$Z
      ICL_val_merge <- results$ICL_val_merge


      ICL_del <- ICL_val_merge - ICL_max


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


ICLFit <- function(Z, Y, G, alpha_var, beta_var) {
  ICL_val_max <- 0
  Z_max <- 0
  iter <- 0
  while(TRUE){
    results<-ICLSweep(Z, Y, G, alpha_var, beta_var)

    G <- results$G_max
    Z <- results$Z_max
    ICL_val <- results$ICL_old

    res <- ICLGroupMerge(ICL_val, Z, Y, G, alpha_var, beta_var)

    ICL_val <- res$ICL_max
    Z_merge <- res$Z_max

    iter <- iter + 1
    Z <- Z_merge

    if(G == ncol(Z_merge)){
      ICL_val_max <- ICL_val
      Z_max <- Z_merge
      break
    }
    if(iter==10){
      ICL_val_max <- ICL_val
      Z_max <- Z_merge
      break
    }
    G <- ncol(Z)
  }
  return(list("Z_max"=Z_max, "ICL_val_max"=ICL_val_max))
}

res_final <- ICLFit(Z, Y, G, alpha_var, beta_var)
ICL_val <- res_final$ICL_val_max
Z <- res_final$Z_max

print(ICL_val)
print(ncol(Z))

source('R\\ICL.Calc.R')

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

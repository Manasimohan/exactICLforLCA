source('R\\ICL.Calc.R')

#' ICL Sweep
#'
#' @param ICL value, Z, Y, G, alpha and beta
#' @return Z_max, ICL_max
#'
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

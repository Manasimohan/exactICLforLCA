source('R\\ICL.Sweep.R')
source('R\\ICL.GroupMerge.R')

#' ICL Fit
#'
#' @param Z, Y, G, alpha and beta
#' @return ICL
#' @export
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

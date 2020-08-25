source('R\\ICL.Calc.R')
source('R\\ICL.GroupReduction.R')

#' ICL Sweep
#'
#' @param Z, Y, G, alpha and beta
#' @return G_max, Z_max, ICL_old
#'
ICLSweep <- function(Z, Y, G, alpha_var, beta_var, delta_var) {
  ICL_old <- ICLCalc(alpha_var, beta_var, G, Y, Z, delta_var)
  ICL_original <- ICL_old

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
    ICL_val_new <- ICLCalc(alpha_var, beta_var, G, Y, Z, delta_var)

    del <- ICL_val_new - ICL_old

    iter <- iter + 1
    if(del <= 0){
      break
    }
    if(iter==50){
      break
    }
    Z_max <- Z
    G_max <- G
  }
  return(list("G_max"=G_max, "Z_max"=Z_max, "ICL_old"=ICL_old, "ICL_original"=ICL_original))
}

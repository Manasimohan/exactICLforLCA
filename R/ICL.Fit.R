source('R\\ICL.Sweep.R')
source('R\\ICL.GroupMerge.R')

pkg <- new.env()
pkg$results <- list()

#' ICL Fit
#'
#' @param Z, Y, G, alpha and beta
#' @return ICL
#' @export
ICLFit <- function(Z, Y, G, alpha_var, beta_var) {
  ICL_val_max <- 0
  Z_max <- 0
  iter <- 0

  pkg$G <- G

  icl_list <- c()
  g_list <- c()

  while(TRUE){
    results <- ICLSweep(Z, Y, G, alpha_var, beta_var)

    if(iter == 0) {
      ICL_original <- results$ICL_original

      icl_list <- c(icl_list, ICL_original)
      g_list <- c(g_list, G)
    }
    G <- results$G_max
    Z <- results$Z_max
    ICL_val <- results$ICL_old


    res <- ICLGroupMerge(ICL_val, Z, Y, G, alpha_var, beta_var)

    ICL_val <- res$ICL_max
    Z_merge <- res$Z_max

    iter <- iter + 1
    Z <- Z_merge

    if(G == ncol(Z_merge)){
      icl_list <- c(icl_list, ICL_val)
      g_list <- c(g_list, ncol(Z_merge))

      ICL_val_max <- ICL_val
      Z_max <- Z_merge
      break
    }
    if(iter==50){
      icl_list <- c(icl_list, ICL_val)
      g_list <- c(g_list, ncol(Z_merge))

      ICL_val_max <- ICL_val
      Z_max <- Z_merge
      break
    }

    icl_list <- c(icl_list, ICL_val)
    g_list <- c(g_list, G)

    G <- ncol(Z)
  }

  results <- list("Z_max"=Z_max, "ICL_val_max"=ICL_val_max, "ICL_original"=ICL_original, "icl_list"=icl_list, "g_list"=g_list)
  pkg$results <- results

  return(results)
}

ICLsummary <- function() {
  results <- pkg$results
  G <- pkg$G
  print(paste('Oringinal ICL value : ', format(results$ICL_original, digits = 5)))
  print(paste("Incial number of clusters : ", G))
  print(paste("ICL value post processing : ", format(results$ICL_val_max, digits = 5)))
  print(paste("Number of clusters post processing : ", ncol(results$Z)))
  print(results$icl_list)
  print(results$g_list)
}

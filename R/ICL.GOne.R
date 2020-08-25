source('R\\ICL.Calc.R')

#' ICL GOne
#'
#' @param Z, Y, G, alpha and beta
#' @return ICL
#' @export
ICLGOne <- function(Z, Y, G, alpha_var, beta_var, delta_var) {
  Zrows <- as.matrix(rep(1, nrow(Z)))
  ICL_one <- ICLCalc(alpha_var, beta_var, ncol(Zrows), Y, Zrows, delta_var)
  return(ICL_one)
}

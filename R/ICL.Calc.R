#' ICL Calculation
#'
#' @param Z, Y, G, alpha and beta
#' @return ICL

ICLCalc <- function(alpha_var, beta_var, G, Y, Z, delta_var) {
  # set variables
  delta <- rep(delta_var , G)
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

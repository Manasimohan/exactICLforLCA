#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

NumericMatrix repeat(int val, int times) {
  NumericMatrix result(1, times);

  for(int i = 0; i < times; i++) {
    result(0, i) = val;
  }
  return result;
}

NumericMatrix applysum(NumericMatrix mat) {
  NumericMatrix result(1, mat.ncol());

  for(int i = 0; i < mat.nrow(); i++) {
    for(int j = 0; j < mat.ncol(); j++) {
      result(0, j) = result(0, j) + mat(i, j);
    }
  }
  return result;
}

float applysum_vec(NumericVector vec) {
  float result = 0;

  for(int i = 0; i < vec.length(); i++) {
    result = result + vec[i];
  }
  return result;
}

NumericMatrix matsum(NumericMatrix mat1, NumericMatrix mat2) {
  NumericMatrix result(mat1.nrow(), mat1.ncol());

  for(int i = 0; i < mat1.nrow(); i++) {
    for(int j = 0; j < mat1.ncol(); j++) {
      result(i, j) = mat1(i, j) + mat2(i, j);
    }
  }

  return result;
}

float log_beta_vec(NumericMatrix delta) {
  // for l_num
  float l_num = applysum_vec(lgamma(delta));

  // for l_denom
  float l_denom = lgamma(applysum_vec(delta));

  return l_num - l_denom;
}


// [[Rcpp::export]]
float ICLCalc(int alpha_var, int beta_var, int G, NumericMatrix Y, NumericMatrix Z, int delta_var) {

  NumericMatrix delta = repeat(delta_var, G);

  // find ncol(Y)
  int r = Y.ncol();

  // initalise alpha and beta matrix
  NumericMatrix alpha_gj( G , r );
  NumericMatrix beta_gj( G , r );

  //NumericMatrix delta_temp = ;
  NumericMatrix delta_prime = matsum(applysum(Z) , delta);

  // calc alpha and beta of groups based on Y and Z
  for(int j = 0; j <= r; j++){
    for(int g = 0; g <= G; g++) {
      int temp_alpha = 0;
      int temp_beta = 0;
      for(int i = 0; i <= Z.nrow(); i++) {
        temp_alpha = temp_alpha + Z(i, g) * Y(i, j);
        temp_beta = temp_beta + Z(i, g) * (1- Y(i, j));
      }
      alpha_gj(g, j) = alpha_var + temp_alpha;
      beta_gj(g, j) = beta_var + temp_beta;
    }
  }

  // first eqn
  float first_var = log_beta_vec(delta_prime) - log_beta_vec(delta);

  float b_num = 0;

  for (int g = 0; g<G; g++)
  {
    for(int j = 0; j<r; j++)
    {
      b_num = b_num + R::lbeta(alpha_gj(g,j), beta_gj(g,j));
    }
  }

  // second eqn denomenarator value
  float b_denom = G * r * R::lbeta(alpha_var, beta_var);

  // second eqn
  float sec_var = b_num - b_denom;

  //ICL calc
  float ICL = first_var + sec_var;

  return ICL;
}

// Exercise Rcpp. Translate determinant.R to Rcpp codes.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
double Rcpp_compute_determinant(mat A) {
  if (A.n_cols != A.n_rows) {
    Rcpp::stop("Only square matrix has determinant.");
  }
  
  double res = 0.0;
  if (A.n_rows == 2) {
    return A(0,0) * A(1,1) - A(0,1) * A(1,0);
  } else {
    for (unsigned int i = 0; i < A.n_cols; i++) {
      double sign = (i % 2 == 0 ? 1.0: -1.0);
      mat B = A;
      B.shed_row(0);
      B.shed_col(i);
      res = res + sign * A(0,i) * Rcpp_compute_determinant(B);
    }
  }
  return res;
}

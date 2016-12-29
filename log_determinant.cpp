// Exercise. Optimise function "log_det_pd" with Rcpp.

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
double cpp_log_det_pd(mat A, int m, int n, vec cheb_coeff, vec s) {
  int d = A.n_rows;
  double G = 0;
  
  for (int i = 0; i < m; i++) {
    vec v = Rcpp ::RcppArmadillo::sample(s, d, true);
    vec u = cheb_coeff(0) * v;
    if (n > 1) {
      vec w0 = v;
      vec w1 = A * v;
      u = u + cheb_coeff(1) * A * v;
      for (int j = 1; j < n; j++) {
        vec w2 = 2 * A * w1 - w0;
        u = u + cheb_coeff(j+1) * w2;
        w0 = w1;
        w1 = w2;
      }
    }
    G = G + dot(v,u) / m;
  }
  return G;
}

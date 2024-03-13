#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;
//using namespace arma;

// [[Rcpp::export]]
arma::mat arma_inv (arma::mat M){
  arma::mat  c_inv = inv_sympd(M);
  return c_inv;
}

// [[Rcpp::export]]
double arma_log_det(arma::mat M){
  double  ld  = log_det_sympd(M);
  return ld;
}

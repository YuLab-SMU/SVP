#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//// This script was from the CelliD package
//' Obtain the pair distance of row between \code{Ar} and \code{Br} matrix
//' @param Ar matrix which number of column should be equal to column number of \code{Br}.
//' @param Br matrix which number of column should be equal to column number of \code{Ar}.
//' @return a distance matrix of row feature in \code{Ar} and \code{Br}.
//[[Rcpp::export]]
NumericMatrix fastPDist(NumericMatrix Ar, NumericMatrix Br) {
    int m = Ar.nrow(), 
        n = Br.nrow(),
        k = Ar.ncol();
    arma::mat A = arma::mat(Ar.begin(), m, k, false);
    arma::mat B = arma::mat(Br.begin(), n, k, false); 
    arma::colvec An =  sum(square(A), 1);
    arma::colvec Bn =  sum(square(B), 1);
    arma::mat C = -2 * (A * B.t());
    C.each_col() += An;
    C.each_row() += Bn.t();
    return wrap(sqrt(C)); 
}

//[[Rcpp::export]]
arma::mat fusiondist(
        arma::mat s,
        arma::mat p,
        double alpha = 0.2,
        double beta = 0.1){
  double alpha_s = 1.0 - alpha;
  double beta_s = 1.0 - beta;
  arma::mat z = beta_s * (alpha_s * s + alpha * p) + beta * s % p;
  return (z);
}

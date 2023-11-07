#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//// This script was from the CelliD package
//' Obtain the pair distance of row between Ar and Br matrix
//' @param Ar matrix which number of column should be equal to column number of Br.
//' @param Br matrix which number of column should be equal to column number of Ar.
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

#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
using namespace arma;

//' Obtain the pair distance of row between \code{A} and \code{B} matrix
//' @param A matrix which number of column should be equal to column number of \code{B}.
//' @param B matrix which number of column should be equal to column number of \code{A}.
//' @return a distance matrix of row feature in \code{A} and \code{B}.
//[[Rcpp::export]]
SEXP fastPDist(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::VectorXd An = A.array().square().rowwise().sum();
    Eigen::VectorXd Bn = B.array().square().rowwise().sum();
    Eigen::MatrixXd C = -2 * (A * B.transpose());
    C.colwise() += An;
    C.rowwise() += Bn.transpose();

    return wrap(C.array().sqrt());

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

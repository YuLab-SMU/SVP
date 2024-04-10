#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace RcppParallel;
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


struct colorder : public Worker{
  const arma::mat& x;
  arma::umat& result;

  colorder(const arma::mat& x, arma::umat& result):
      x(x), result(result){}

  void operator()(std::size_t begin, std::size_t end){
      for (uword i = begin; i < end; i++){
          result.col(i) = arma::sort_index(x.col(i), "ascend") + 1;
      }
  }
};

// [[Rcpp::export]]
NumericMatrix ParallelColOrder(const arma::mat& x, int top_n){
    arma::umat ordmat = arma::umat(x.n_rows, x.n_cols);
    colorder runcolorder(x, ordmat);
    parallelFor(0, x.n_cols, runcolorder);
    return wrap(ordmat.head_rows(top_n));
}

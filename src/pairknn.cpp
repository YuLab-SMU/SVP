#include <RcppArmadillo.h>
//#include <RcppEigen.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
//using namespace Eigen;
using namespace arma;

//// Obtain the pair distance of row between \code{A} and \code{B} matrix
//// param A matrix which number of column should be equal to column number of \code{B}.
//// param B matrix which number of column should be equal to column number of \code{A}.
//// return a distance matrix of row feature in \code{A} and \code{B}.
////[[Rcpp::export]]
//SEXP fastPDist(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
//    Eigen::VectorXd An = A.array().square().rowwise().sum();
//    Eigen::VectorXd Bn = B.array().square().rowwise().sum();
//    Eigen::MatrixXd C = -2 * (A * B.transpose());
//    C.colwise() += An;
//    C.rowwise() += Bn.transpose();
//
//    return wrap(C.array().sqrt());
//
//}


////[[Rcpp::export]]
//arma::mat fusiondist(
//        arma::mat s,
//        arma::mat p,
//        double alpha = 0.2,
//        double beta = 0.1){
//  double alpha_s = 1.0 - alpha;
//  double beta_s = 1.0 - beta;
//  arma::mat z = beta_s * (alpha_s * s + alpha * p) + beta * s % p;
//  return (z);
//}

arma::vec fastPDist2(arma::mat A, arma::vec B, arma::vec An, double Bn){
    arma::vec C = -2 * (A * B);
    C = sqrt(C + An + Bn);
    C.replace(arma::datum::nan, 6.629066e-15);
    return (C);
}

struct CalPairKnn : public Worker{
  const arma::mat& x;
  const arma::mat& y;
  const arma::vec& an;
  const arma::vec& bn;
  const uword& k;
  arma::umat& res1;
  arma::mat& res2;

  CalPairKnn(const arma::mat& x, const arma::mat& y, const arma::vec& an,
             const arma::vec& bn, const uword& k, arma::umat& res1,
             arma::mat& res2
    ): x(x), y(y), an(an), bn(bn), k(k), res1(res1), res2(res2){}

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
       arma::vec tmp = fastPDist2(y, x.col(i), an, bn(i));
       arma::uvec tind = arma::sort_index(tmp, "ascend");
       arma::vec dd = tmp(tind);
       res1.col(i) = tind.subvec(0, k);
       res2.col(i) = dd.subvec(0, k);
    }
  }
};

//[[Rcpp::export]]
List pairKnnCpp(arma::mat x, arma::mat y, arma::uword topn = 2){
    arma::vec an = sum(square(y), 1);
    arma::vec bn = sum(square(x), 1);
    x = x.t();
    arma::umat res1(topn + 1, x.n_cols);
    arma::mat res2(topn + 1, x.n_cols);
    CalPairKnn runpairknn(x, y, an, bn, topn, res1, res2);
    parallelFor(0, x.n_cols, runpairknn);
    return List::create(Named("index") = res1 + 1,
                        Named("distance") = res2);
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

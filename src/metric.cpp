#include "autocorutils.h"
#include "progress.h"
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

struct calf1 : public Worker{
  const arma::sp_mat& x;
  const arma::sp_mat& y;
  simple_progress& p;
  mat& result;

  calf1(const arma::sp_mat& x, const arma::sp_mat& y,
        simple_progress& p, mat& result)
  : x(x), y(y), p(p), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        for (uword j = 0; j < y.n_cols; j ++){
            result(i, j) = calculateF1(x.col(i).as_dense(), y.col(j).as_dense());
            p.increment();
        }
    }
  }
};

//[[Rcpp::export]]
arma::mat CalF1Parallel(
  arma::sp_mat x,
  arma::sp_mat y
  ){
  y = y.t();
  arma::mat res(x.n_cols, y.n_cols);
  simple_progress p (x.n_cols * y.n_cols);

  calf1 runcalf1(x, y, p, res);
  parallelFor(0, x.n_cols, runcalf1);
  res.replace(arma::datum::nan, 0);
  return(res);
}

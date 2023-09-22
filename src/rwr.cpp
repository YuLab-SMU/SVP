#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppParallel)]]
using namespace RcppParallel;
using namespace arma;
using namespace std;

struct calrwr : public Worker{
  const arma::sp_mat& x;
  const arma::sp_mat& v;
  const double restart;
  const double stop_delta;
  const int stop_step;

  mat& result;

  calrwr(const arma::sp_mat& x, const arma::sp_mat& v, const double restart,
         const double stop_delta, const int stop_step, mat& result)
  : x(x), v(v), restart(restart), stop_delta(stop_delta),
  stop_step(stop_step), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        int step = 0;
        double delta = 1;
        arma::sp_mat pt = v.col(i);
        while((delta > stop_delta) && (step <= stop_step)){
            arma::sp_mat px = ((1 - restart) * (x * pt)) + (restart * v.col(i));
            delta = arma::accu(abs(px - pt));
            pt = px;
            step = step + 1;
        }
        result.col(i) = conv_to<mat>::from(pt);
    }
  }
};

// [[Rcpp::export]]
NumericMatrix parallelCalRWR(
              arma::sp_mat x,
              arma::sp_mat v,
              double restart = 0.75,
              double stop_delta = 0.000001,
              int stop_step = 2){

    uword n = v.n_cols;
    mat result(x.n_rows, n);

    calrwr runrwr(x, v, restart, stop_delta, stop_step, result);

    parallelFor(0, n, runrwr);
    
    return wrap(result);

}



//#include <RcppArmadillo.h>
//// [[Rcpp::depends(RcppArmadillo)]]
//
//using namespace Rcpp;
//using namespace arma;
//
//// [[Rcpp::export]]
//NumericMatrix calRWRCPP(arma::sp_mat x,
//	      arma::sp_mat v,
//	      double restart = .7, 
//	      double stop_delta = 1e-6,
//	      int stop_step = 50
//	     ){
//    int step = 0;
//    double delta = 1;
//    arma::sp_mat pt = v;
//    while((delta > stop_delta) && (step <= stop_step)){
//	arma::sp_mat px = ((1 - restart) * x * pt) + (restart * v);
//	delta = arma::accu(abs(px - pt));
//	pt = px;
//      step = step + 1;
//    }
//    arma::mat pt2 = arma::conv_to<mat>::from(pt);
//    return wrap(pt2);
//}
//
//// [[Rcpp::export]]
//NumericMatrix calRWRCPP2(
//	      arma::sp_mat x,
//              arma::sp_mat v,
//              double restart = .7,
//              double stop_delta = 1e-6,
//              int stop_step = 50
//             ){
//    int n = v.n_cols;
//    mat res(x.n_rows, n);
//    for (int i = 0; i < n; i++){
//      int step = 0;
//      double delta = 1;
//      arma::sp_mat pt = v.col(i);
//      while((delta > stop_delta) && (step <= stop_step)){
//          arma::sp_mat px = ((1 - restart) * x * pt) + (restart * v.col(i));
//          delta = arma::accu(abs(px - pt));
//          pt = px;
//          step = step + 1;
//      }
//	res(i) = arma::conv_to<mat>::from(pt);
//    }
//    return wrap(res);
//}


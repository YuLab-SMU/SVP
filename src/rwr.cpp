#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
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
        arma::mat pt = conv_to<mat>::from(v.col(i));
        arma::uvec ind = find(pt > 0);
        while((delta > stop_delta) && (step <= stop_step) && (ind.n_elem > 1)){
            arma::mat px = ((1 - restart) * (x * pt)) + (restart * v.col(i));
            delta = arma::accu(abs(px - pt));
            pt = px;
            step = step + 1;
        }
        result.col(i) = pt;
    }
  }
};


//' Computer the affinity score of all nodes in a graph to a seeds 
//' using Random Walk with Restart
//' @param x a adjacency matrix of a graph.
//' @param v a matrix define sets of starting seeds, each column 
//' corresponds to one set of seeds that a walker starts.
//' @param restart the restart probability used for RWR, it must be 
//' between 0 and 1, default is .75.
//' @param stop_delta minimum threshold to stop RWR, default is 1e-10.
//' @param stop_step step number to stop RWR, default is 50.
// [[Rcpp::export]]
NumericMatrix parallelCalRWR(
              arma::sp_mat x,
              arma::sp_mat v,
              double restart = 0.75,
              double stop_delta = 1e-10,
              int stop_step = 50){

    uword n = v.n_cols;
    mat result(x.n_rows, n);

    calrwr runrwr(x, v, restart, stop_delta, stop_step, result);

    parallelFor(0, n, runrwr);
    
    return wrap(result);

}


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
//
//// [[Rcpp::export]]
//NumericMatrix calRWRCPP2(
//	      arma::sp_mat x,
//              arma::sp_mat v,
//              double restart = .7,
//              double stop_delta = 1e-10,
//              int stop_step = 50
//             ){
//    int n = v.n_cols;
//    arma::mat res(x.n_rows, n);
//    for (int i = 0; i < n; i++){
//      int step = 0;
//      double delta = 1;
//      arma::mat pt = arma::conv_to<mat>::from(v.col(i));
//      arma::uvec ind = find(pt > 0);
//      while((delta > stop_delta) && (step <= stop_step) && (ind.n_elem > 1)){
//          arma::mat px = ((1 - restart) * x * pt) + (restart * v.col(i));
//          delta = arma::accu(abs(px - pt));
//          pt = px;
//          step = step + 1;
//      }
//	res.col(i) = pt;
//    }
//    return wrap(res);
//}


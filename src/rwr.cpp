#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix calRWRCPP(arma::sp_mat x,
	      arma::sp_mat v,
	      double restart = .7, 
	      double delta = 1,
	      int step = 0,
	      double stop_delta = 1e-6,
	      int stop_step = 50
	     ){
    arma::sp_mat pt = v;
    while((delta > stop_delta) && (step <= stop_step)){
	arma::sp_mat px = ((1 - restart) * x * pt) + (restart * v);
	//arma::mat tmp = arma::conv_to<mat>::from(px - pt);
	delta = arma::accu(abs(px - pt));
	pt = px;
        step = step + 1;
    }
    arma::mat pt2 = arma::conv_to<mat>::from(pt);
    return wrap(pt2);
}

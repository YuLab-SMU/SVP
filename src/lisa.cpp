#include "autocorutils.h"
#include "progress.h"
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

struct RunLocalG: public Worker{
    const arma::mat& xm;
    const arma::mat& w;
    const arma::vec& wi;
    const arma::vec& Wi2;
    simple_progress& p;
    const int n;
    arma::mat& result1;
    arma::mat& result2;
    arma::mat& result3;
    arma::mat& result4;
    arma::mat& result5;

    RunLocalG(const arma::mat& xm, const arma::mat& w, const arma::vec& wi,
            const arma::vec& Wi2, simple_progress& p, const int n, arma::mat& result1, 
	    arma::mat& result2, arma::mat& result3, arma::mat& result4, arma::mat& result5
    ):
    xm(xm), w(w), wi(wi), Wi2(Wi2), p(p), n(n), result1(result1), result2(result2),
    result3(result3), result4(result4), result5(result5){ }

    void operator()(std::size_t begin, std::size_t end){
        for (uword i = begin; i < end; i++){
            arma::mat res(n, 6);
            res = CalLocalGCpp(xm.col(i), w, wi, Wi2, n);
            result1.col(i) = res.col(0);
            result2.col(i) = res.col(1);
            result3.col(i) = res.col(2);
            result4.col(i) = res.col(3);
            result5.col(i) = res.col(4);
        }
    }
};

// [[Rcpp::export]]
Rcpp::List CalLocalGParallel(arma::sp_mat& x, arma::mat& w){
    arma::mat xm = conv_to<arma::mat>::from(x.t());
    int n = xm.n_cols;
    int m = xm.n_rows;

    arma::vec wi = sum(w, 1);
    arma::vec Wi2 = sum(pow(w, 2.0), 1);

    simple_progress p(n);

    arma::mat result1(m, n);
    arma::mat result2(m, n);
    arma::mat result3(m, n);
    arma::mat result4(m, n);
    arma::mat result5(m, n);

    RunLocalG runlocalg(xm, w, wi, Wi2, p, m, result1, result2,
            result3, result4, result5);

    parallelFor(0, n, runlocalg);

    List res(n);
    for (int i = 0; i < n; i++){
        arma::mat tmp(m, 5);
        res[i] = tidylocalg(result1.col(i), result2.col(i), result3.col(i),
                result4.col(i), result5.col(i), tmp);
        p.increment();
    }

    return (res);
}

struct RunLocalMoran : public Worker{
    const arma::mat& xm;
    const arma::mat& w;
    const arma::vec& wi;
    const arma::vec& Wi2;
    simple_progress& p;
    const int n;
    arma::mat& result1;
    arma::mat& result2;
    arma::mat& result3;
    arma::mat& result4;
    arma::mat& result5;
    arma::mat& result6;

    RunLocalMoran(const arma::mat& xm, const arma::mat& w, const arma::vec& wi,
            const arma::vec& Wi2, simple_progress& p, const int n, arma::mat& result1, arma::mat& result2, 
            arma::mat& result3, arma::mat& result4, arma::mat& result5, arma::mat& result6
    ):
    xm(xm), w(w), wi(wi), Wi2(Wi2), p(p), n(n), result1(result1), result2(result2), 
    result3(result3), result4(result4), result5(result5), result6(result6){ }

    void operator()(std::size_t begin, std::size_t end){
        for (uword i = begin; i < end; i++){
            arma::mat res(n, 6);
            res = CalLocalMoranCpp(xm.col(i), w, wi, Wi2, n);
            result1.col(i) = res.col(0);
            result2.col(i) = res.col(1);
            result3.col(i) = res.col(2);
            result4.col(i) = res.col(3);
            result5.col(i) = res.col(4);
            result6.col(i) = res.col(5);
            p.increment();
        }
    }
};


// [[Rcpp::export]]
Rcpp::List CalLocalMoranParallel(arma::sp_mat& x, arma::mat& w){
    arma::mat xm = conv_to<arma::mat>::from(x.t());
    int n = xm.n_cols;
    int m = xm.n_rows;

    arma::vec wi = sum(w, 1);
    arma::vec Wi2 = sum(pow(w, 2.0), 1);

    simple_progress p(n);
    
    arma::mat result1(m, n);
    arma::mat result2(m, n);
    arma::mat result3(m, n);
    arma::mat result4(m, n);
    arma::mat result5(m, n);
    arma::mat result6(m, n);

    RunLocalMoran runlocalmoran(xm, w, wi, Wi2, p, m, result1, result2,
            result3, result4, result5, result6);

    parallelFor(0, n, runlocalmoran);

    List res(n);
    for (int i = 0; i < n; i++){
        arma::mat tmp(m, 6);
        res[i] = tidylocalmoran(result1.col(i), result2.col(i), result3.col(i),
                result4.col(i), result5.col(i), result6.col(i), tmp);
    }

    return (res);
}


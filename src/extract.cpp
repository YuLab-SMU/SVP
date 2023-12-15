#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector SortNv(NumericVector x, bool decreasing = true) {
    IntegerVector idx = seq_along(x) - 1;
    if (decreasing){
        std::sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] > x[j];});
    }else{
        std::sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] < x[j];});
    }
    CharacterVector nm = x.names();
    NumericVector y = x[idx];
    y.names() = nm[idx];
    return y;
}

//' Extract the score of gene in each gene sets
//' @param x the score sparse matrix of gene for each gene sets.
//' @param rnm the row names of x matrix.
//' @param cnm the col names of x matrix.
//' @param g a list of gene set.
//' @return a list contained the score of gene in each gene sets
//[[Rcpp::export]]
List ExtractFeatureScoreCpp(arma::sp_mat& x,
              CharacterVector& rnm,
              CharacterVector& cnm,
              Rcpp::List& g
              ){
    arma::mat xx = conv_to<arma::mat>::from(x);
    uword n = rnm.length();
    Rcpp::List res(n);
    for (uword i = 0; i < n; i++){
        String gnm = rnm[i];
        CharacterVector gene = g[gnm];
        LogicalVector ind = in(cnm, gene);
        CharacterVector nm = cnm[ind];
        NumericVector f = as<NumericVector>(wrap(xx.row(i)));
        NumericVector fn = f[ind];
        fn.names() = nm;
        res[i] = SortNv(fn);
    }
    return wrap(res);
}

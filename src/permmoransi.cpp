#include <RcppParallel.h>
#include <RcppArmadillo.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

// this is refer to the Rfast
static double calc_pvalue_rnd_r(arma::vec& ds_c0, arma::vec& ds_c1, const unsigned int r,
                const double upper, const double lower, const double test_stat, const unsigned int nrows) {
        const unsigned int rnd_r = std::round(std::sqrt(r));
        arma::mat lh = arma::mat(ds_c0.size(), rnd_r, arma::fill::zeros);
        arma::mat rh = arma::mat(ds_c0.size(), rnd_r, arma::fill::zeros);
        for (unsigned int i = 0; i < rnd_r; ++i) {
                lh.col(i) = arma::shuffle(ds_c0);
                rh.col(i) = arma::shuffle(ds_c1);
        }
        arma::mat sxy = (arma::trans(lh) * rh);
        unsigned int sum = 1;
        for (unsigned int i = 0; i < sxy.n_rows; ++i) {
                for (unsigned int j = 0; j < sxy.n_cols; ++j) {
                        const double curr = (sxy(i, j) - upper) / lower;
                        const double curr_test_stat = std::abs(std::log((1 + curr) / (1 - curr)));
                        if (curr_test_stat > test_stat) {
                                ++sum;
                        }
                }
        }
        return sum / (double) (std::pow(rnd_r, 2) + 1);
}

// this is also refer to the Rfast
arma::rowvec calc_perm_cor(arma::vec& x, arma::vec& y, const unsigned int r) {
        const unsigned int nrows = x.size();

    double ds_c0_sum = 0;
    double ds_c1_sum = 0;
    double ds_c0_p2_sum = 0;
    double ds_c1_p2_sum = 0;
    double ds_sum = 0;
        for (unsigned int i = 0; i < nrows; i++) {
                const double curr_c0 = x[i];
                const double curr_c1 = y[i];
                ds_c0_sum += curr_c0;
                ds_c1_sum += curr_c1;
                ds_c0_p2_sum += std::pow(curr_c0, 2.0);
                ds_c1_p2_sum += std::pow(curr_c1, 2.0);
                ds_sum += curr_c0 * curr_c1;
        }
    const double upper = (ds_c0_sum * ds_c1_sum) / nrows;
    const double lower = std::sqrt((ds_c0_p2_sum - std::pow(ds_c0_sum, 2.0) / nrows) *
                        (ds_c1_p2_sum - std::pow(ds_c1_sum, 2.0) / nrows));

    const double cor = (ds_sum - upper) / lower;
    const double test_stat = std::abs(std::log((1 + cor) / (1 - cor)));
    const double pvalue = calc_pvalue_rnd_r(x, y, r, upper, lower, test_stat, nrows);
    arma::rowvec res = {cor, pvalue};
    return (res);
}


arma::rowvec permmoransi(arma::vec x, arma::mat weight, int permutation=999){
    arma::vec y = weight * x;
    arma::rowvec res = calc_perm_cor(y, x, permutation);
    res[0] = res[0] * var(y) / var(x);
    return(res);
}


struct RunPermMoransi : public Worker{
  const arma::mat& x;
  const arma::mat& weight;
  const int permutation;
  arma::mat& result;

  RunPermMoransi(const arma::mat& x, const arma::mat& weight,
          const int permutation, mat& result):
      x(x), weight(weight), permutation(permutation), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        result.row(i) = permmoransi(x.col(i), weight, permutation);
    }
  }
};

// [[Rcpp::export]]
arma::mat CalMoransiPermParallel(arma::sp_mat& x, arma::mat& weight, int permutation = 999){
  arma::mat xm =  conv_to<arma::mat>::from(x);
  xm = arma::trans(xm);
  
  int n = xm.n_cols;

  arma::colvec rowsumw = sum(weight, 1);
  rowsumw.replace(0, 1);
  weight = weight.each_col() / rowsumw;
  
  arma::mat result(n, 2);

  RunPermMoransi runpermmoransi(xm, weight, permutation, result);

  parallelFor(0, n, runpermmoransi);  

  return(result);

}

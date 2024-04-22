#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <algorithm>
#include <xoshiro.h>
#include <convert_seed.h>
#include <R_randgen.h>
#include "buildrand.h"

using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat outersubtractdot(arma::rowvec x){
  arma::mat xm = repelem(x, x.n_elem, 1);
  arma::mat res = pow(xm - xm.t(), 2);
  return (res);
}

double cal_gearysc(
          arma::rowvec x,
          arma::mat weight,
          double s,
          int n
        ){

    double m = mean(x);
    arma::rowvec y = x - m;
    arma::mat ym = outersubtractdot(x);
    double cv = accu(weight % ym);
    double v = accu(pow(y, 2));
    double obs = ((n - 1) * cv)/(s * v);

    return(obs);
}


arma::rowvec cal_gearysc_p_perm(
        arma::rowvec x,
        arma::mat weight,
        dqrng::xoshiro256plus rng,
        int permutation,
        double s,
        int n
        ){
    double obs = cal_gearysc(x, weight, s, n);
    arma::vec xmr = arma::vec(permutation, arma::fill::zeros);
    for (int i = 0; i < permutation; i++){
        arma::rowvec z = x;
        std::shuffle(std::begin(z), std::end(z), rng);
        xmr(i) = cal_gearysc(z, weight, s, n);
    }
    
    double expv = mean(xmr);
    double sdv = stddev(xmr);

    double pv = R::pnorm(obs, expv, sdv, 1, 0);

    arma::rowvec res = {obs, expv, sdv, pv};
    return(res);
  
}

struct RunGearysc : public Worker{
  const arma::mat& x;
  const arma::mat& weight;
  const uint64_t seed;
  const int permutation;
  const double s;
  const int n;
  arma::mat& result;

  RunGearysc(const arma::mat& x, const arma::mat& weight, const uint64_t seed, 
          const int permutation, const double s, const int n, mat& result):
      x(x), weight(weight), seed(seed), permutation(permutation), s(s),
      n(n), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    dqrng::xoshiro256plus rng(seed);
    for (uword i = begin; i < end; i++){
        dqrng::xoshiro256plus lrng(rng);
        lrng.long_jump(i + 1);
        result.row(i) = cal_gearysc_p_perm(x.row(i), weight, lrng, permutation, s, n);
    }
  }
};


// [[Rcpp::export]]
arma::mat CalGearyscParallel(arma::sp_mat& x, arma::mat& weight, int permutation = 999){
  arma::mat xm =  conv_to<arma::mat>::from(x);
  int n = x.n_rows;
  int m = x.n_cols;
 
  weight.diag().fill(0);
  arma::colvec rowsumw = sum(weight, 1);
  rowsumw.replace(0, 1);
  weight = weight.each_col() / rowsumw;

  double s = accu(weight);

  Rcpp::IntegerVector seed(2, dqrng::R_random_int);
  uint64_t seed2 = dqrng::convert_seed<uint64_t>(seed);

  arma::mat result(n, 4);
  RunGearysc rungearysc(xm, weight, seed2, permutation, s, m, result);
  parallelFor(0, n, rungearysc);

  //#ifdef _OPENMP
  //#pragma omp parallel for schedule(static)
  //#endif
  //for (unsigned int i = 0; i < n; i++){
  //    result.row(i) = cal_gearysc_p_perm(xm.row(i), weight, rmat, s, m);
  //}

  return(result);

}

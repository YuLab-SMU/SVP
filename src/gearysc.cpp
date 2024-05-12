#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <algorithm>
#include <xoshiro.h>
#include <convert_seed.h>
#include <R_randgen.h>
#include "progress.h"
#include "autocorutils.h"
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::rowvec cal_gearysc_p_perm(
        double obs,
        arma::rowvec x,
        arma::mat weight,
        arma::rowvec rowsumw,
        dqrng::xoshiro256plus rng,
        int permutation,
        double s,
        int n,
        int lower_tail
        ){
    arma::vec xmr = arma::vec(permutation, arma::fill::zeros);
    for (int i = 0; i < permutation; i++){
        arma::rowvec z = x;
        std::shuffle(std::begin(z), std::end(z), rng);
        xmr(i) = cal_gearysc(z, weight, s, n);
    }
    
    double expv = mean(xmr);
    double sdv = stddev(xmr);

    double z = (expv - obs) / sdv;

    double pv = R::pnorm(obs, expv, sdv, lower_tail, 0);

    arma::rowvec res = {obs, expv, sdv, z, pv};
    return(res);
}

arma::rowvec cal_gearysc_p_noperm(
        double obs,
        arma::rowvec x,
        double ei,
        double s,
        int n,
        double S1,
        double S2,
        int lower_tail
    ){
    
    arma::rowvec y = x - mean(x);
    //y = scaleCpp(y);
    double v = accu(pow(y, 2.0));

    double s_sq = pow(s, 2.0);
    double n2 = pow(n, 2.0);

    double k = (n * accu(pow(y, 4.0)))/pow(v, 2.0);
    double sdi = (n - 1.0) * S1 * (n2 - 3.0 * n + 3.0 - (n - 1.0) *k);
    sdi = sdi - (0.25 * (n - 1.0) * S2 * (n2 + 3.0 * n - 6.0 - (k * (n2 - n + 2.0))));
    sdi = sdi + s_sq * (n2 - 3.0 - pow(n - 1, 2.0) * k);
    sdi = sqrt(sdi/(n * (n - 2.0) * (n - 3.0) * s_sq));

    double z = (ei - obs)/sdi;

    double pv = R::pnorm(obs, ei, sdi, lower_tail, 0);
    arma::rowvec res = {obs, ei, sdi, z, pv};
    return (res);
}


arma::rowvec gearysc(
        arma::rowvec x,
        arma::mat weight,
        arma::rowvec rowsumw,
        dqrng::xoshiro256plus rng,
        int permutation,
        double S1,
        double S2,
        double s,
        int n,
        double ei,
	int lower_tail = 1
        ){

  double obs = cal_gearysc(x, weight, s, n);

  arma::rowvec res(5);
  if (permutation < 100){
      res = cal_gearysc_p_noperm(obs, x, ei, s, n, S1, S2, lower_tail);
  }else{
      res = cal_gearysc_p_perm(obs, x, weight, rowsumw, rng, permutation, s, n, lower_tail);
  }
  return(res);
}

struct RunGearysc : public Worker{
  const arma::mat& x;
  const arma::mat& weight;
  const arma::rowvec& rowsumw;
  simple_progress& p;
  const uint64_t seed;
  const int permutation;
  const double S1;
  const double S2;
  const double s;
  const int n;
  const double ei;
  const int lower_tail;
  arma::mat& result;

  RunGearysc(const arma::mat& x, const arma::mat& weight, const arma::rowvec& rowsumw, 
	  simple_progress& p, const uint64_t seed, const int permutation, const double S1,
	  const double S2, const double s, const int n, const double ei, const int lower_tail, mat& result):
      x(x), weight(weight), rowsumw(rowsumw), p(p), seed(seed), permutation(permutation),
      S1(S1), S2(S2), s(s), n(n), ei(ei), lower_tail(lower_tail), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    dqrng::xoshiro256plus rng(seed);
    for (uword i = begin; i < end; i++){
        dqrng::xoshiro256plus lrng(rng);
        lrng.long_jump(i + 1);
        result.row(i) = gearysc(x.row(i), weight, rowsumw, lrng, permutation, S1, S2, s, n, ei, lower_tail);
        p.increment();
    }
  }
};


// [[Rcpp::export]]
arma::mat CalGearyscParallel(arma::sp_mat& x, arma::mat& weight, int permutation = 999, int lower_tail = 1){
  arma::mat xm =  conv_to<arma::mat>::from(x);
  int n = x.n_rows;
  int m = x.n_cols;
  double ei = 1.0;
 
  arma::rowvec colsumw = sum(weight, 0);
  arma::colvec rowsumw = sum(weight, 1);
  arma::rowvec rowsumw2 = conv_to<arma::rowvec>::from(rowsumw);

  double S1 =  0.5 * accu(pow(weight + weight.t(), 2.0));
  double S2 = accu(pow(rowsumw2 + colsumw, 2.0));  
  double s = accu(weight);

  Rcpp::IntegerVector seed(2, dqrng::R_random_int);
  uint64_t seed2 = dqrng::convert_seed<uint64_t>(seed);

  simple_progress p(n);
  arma::mat result(n, 5);
  RunGearysc rungearysc(xm, weight, rowsumw2, p, seed2, permutation, S1, S2, s, m, ei, lower_tail, result);
  parallelFor(0, n, rungearysc);

  return(result);
}



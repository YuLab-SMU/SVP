#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <xoshiro.h>
#include <convert_seed.h>
#include <R_randgen.h>
#include "progress.h"
#include "autocorutils.h"
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::vec cal_global_lee_test(
  arma::rowvec x,
  arma::rowvec y,
  arma::mat w,
  dqrng::xoshiro256plus rng,
  double S2,
  int n,
  int permutation
  ){
  arma::vec res(permutation) ;
  for (int i = 0; i < permutation; i++){
      arma::rowvec x1 = x;
      arma::rowvec y1 = y;
      std::shuffle(std::begin(x1), std::end(x1), rng);
      std::shuffle(std::begin(y1), std::end(y1), rng);
      res[i] = cal_global_lee(x1, y1, w, S2, n); 
  }
  return(res);
}

struct RunGlobalLee : public Worker{
  const arma::mat& x;
  const arma::mat& w;
  const arma::uvec& f1;
  const arma::uvec& f2;
  simple_progress& p;
  const uword& nf2;
  const double S2;
  const int n;
  const uint64_t seed;
  const int permutation;
  const bool cal_pvalue;
  arma::mat& result;
  
  RunGlobalLee(const arma::mat& x, const arma::mat& w, const arma::uvec& f1,
          const arma::uvec& f2, simple_progress& p, const uword& nf2, 
          const double S2, const int n, const uint64_t seed, const int permutation,
          const bool cal_pvalue, arma::mat& result):
    x(x), w(w), f1(f1), f2(f2), p(p), nf2(nf2), S2(S2), n(n), seed(seed), 
    permutation(permutation), cal_pvalue(cal_pvalue), result(result){} 

  void operator()(std::size_t begin, std::size_t end){
      dqrng::xoshiro256plus rng(seed);
      for (uword i = begin; i < end; i++){
          dqrng::xoshiro256plus lrng(rng);
          for (uword j = 0; j < nf2; j ++){
              double obs = cal_global_lee(x.row(f1[i]), x.row(f2[j]), w, S2, n);
              if (cal_pvalue){
                  arma::vec resp = cal_global_lee_test(x.row(f1[i]), x.row(f2[j]), w, lrng, S2, n, permutation);
                  double pv = cal_permutation_p(resp, obs, permutation + 1);
                  result(i, j) = pv;
              }else{
                  result(i, j) = obs;
              }
              p.increment();
          }
      }
  } 
};

//[[Rcpp::export]]
arma::mat CalGlobalLeeParallel(
  arma::sp_mat x, 
  arma::mat w, 
  arma::uvec f1,
  arma::uvec f2,
  int permutation = 200,
  int alternative = 3,
  bool cal_pvalue = false
  ){
  arma::mat xm = conv_to<arma::mat>::from(x);
  int m = xm.n_cols;

  double S2 = accu(pow(sum(w, 1), 2.0));

  Rcpp::IntegerVector seed0(2, dqrng::R_random_int);
  uint64_t seed1 = dqrng::convert_seed<uint64_t>(seed0);
  
  int n1 = f1.n_elem;
  uword n2 = f2.n_elem;

  simple_progress p (n1 * n2);
  
  arma::mat res(n1, n2);

  RunGlobalLee rungloballeeperm(xm, w, f1, f2, p, n2, S2, m, seed1, permutation, cal_pvalue, res);
       
  parallelFor(0, n1, rungloballeeperm);
  
  return(res);
}

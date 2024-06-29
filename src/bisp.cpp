#include "autocorutils.h"
#include "xoshiro.h"
#include "convert_seed.h"
#include "R_randgen.h"
#include "progress.h"
#include <RcppParallel.h>
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
  const int alternative;
  arma::mat& result;
  arma::mat& presult;
  
  RunGlobalLee(const arma::mat& x, const arma::mat& w, const arma::uvec& f1,
          const arma::uvec& f2, simple_progress& p, const uword& nf2, 
          const double S2, const int n, const uint64_t seed, const int permutation,
          const bool cal_pvalue, const int alternative, arma::mat& result, arma::mat& presult):
    x(x), w(w), f1(f1), f2(f2), p(p), nf2(nf2), S2(S2), n(n), seed(seed), 
    permutation(permutation), cal_pvalue(cal_pvalue), alternative(alternative), 
    result(result), presult(presult){} 

  void operator()(std::size_t begin, std::size_t end){
      dqrng::xoshiro256plus rng(seed);
      for (uword i = begin; i < end; i++){
          dqrng::xoshiro256plus lrng(rng);
          for (uword j = 0; j < nf2; j ++){
              double obs = cal_global_lee(x.row(f1[i]), x.row(f2[j]), w, S2, n);
              result(i, j) = obs;
              if (cal_pvalue){
                  arma::vec resp = cal_global_lee_test(x.row(f1[i]), x.row(f2[j]), w, lrng, S2, n, permutation);
                  double pv = cal_permutation_p(resp, obs, permutation + 1, alternative);
                  presult(i, j) = pv;
              }
              p.increment();
          }
      }
  } 
};

//[[Rcpp::export]]
Rcpp::List CalGlobalLeeParallel(
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
  arma::mat pres(n1, n2);

  RunGlobalLee rungloballeeperm(xm, w, f1, f2, p, n2, S2, m, seed1, permutation, cal_pvalue, alternative, res, pres);
       
  parallelFor(0, n1, rungloballeeperm);
  List result = List::create(Named("Lee")=res, Named("pvalue")=pres);
  
  return(result);
}

//[[Rcpp::export]]
arma::vec RunLocalLee(
    arma::vec x,
    arma::vec y,
    arma::mat weight,
    double n
    ){

    arma::vec dx = x - mean(x);
    arma::vec dy = y - mean(y);

    double dx2 = accu(pow(dx, 2.0));
    double dy2 = accu(pow(dy, 2.0));

    arma::vec ldx = lagCpp(weight, dx);
    arma::vec ldy = lagCpp(weight, dy);

    arma::vec l = (n * ldx % ldy)/(sqrt(dx2) * sqrt(dy2));
    return(l);
}

arma::mat cal_local_moran_bv_perm(
    arma::vec x,
    arma::vec y,
    arma::mat w,
    int n,
    int perm
){
    arma::mat res(n, perm);
    for (int i = 0; i< perm; i++){
        arma::vec y1 = shuffle(y);
        res.col(i) = cal_local_moran_bv(x, y1, w);
    }
    return(res);
}

//[[Rcpp::export]]
arma::mat RunLocalMoranBvPerm(
  arma::vec x,
  arma::vec y,
  arma::mat w,
  int n,
  int permutation = 200
  ){
  
  x = scaleCpp2(x);
  y = scaleCpp2(y);  
  arma::mat result(n, 4);
  
  result.col(0) = cal_local_moran_bv(x, y, w);

  arma::mat res = cal_local_moran_bv_perm(x, y, w, n, permutation);
  result.col(1) = mean(res, 1);
  result.col(2) = var(res, 0, 1);
  result.col(3) = (result.col(0) - result.col(1))/sqrt(result.col(2));

  return(result);
}


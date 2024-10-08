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


arma::rowvec cal_moransi_p_noperm(
        double obs, 
        double ei, 
        double s, 
        arma::vec y, 
        double v, 
        int n, 
        double S1, 
        double S2,
	int lower_tail
  ){
  double s_sq = pow(s, 2.0);

  double k = (accu(pow(y, 4.0))*n)/pow(v, 2.0);
  double sdi = sqrt((n * ((pow(n, 2.0) - 3 * n + 3) * S1 - n * S2 + 3 * s_sq) -
              k * (n * (n - 1.0) * S1 - 2 * n * S2 + 6 * s_sq))/
          ((n - 1.0) * (n - 2.0) * (n - 3.0) * s_sq) - 1/(pow(n - 1.0, 2.0)));

  double z = (obs - ei)/sdi;

  double pv = R::pnorm(obs, ei, sdi, lower_tail, 0);
  arma::rowvec res = {obs, ei, sdi, z, pv};
  return(res);
}

arma::rowvec cal_moransi_p_perm(
        arma::vec x,
        arma::sp_mat weight,
        arma::vec rowsumw,
        dqrng::xoshiro256plus rng,
        int permutation,
        double s,
        int n,
        bool scaled = false,
	int lower_tail = 1
        ){
    double obs = cal_moransi(x, weight, rowsumw, s, n, scaled);
    arma::vec xmr = arma::vec(permutation, arma::fill::zeros);
    for (int i = 0; i < permutation; i++){
        arma::vec z = x;
        std::shuffle(std::begin(z), std::end(z), rng);
        xmr(i) = cal_moransi(z, weight, rowsumw, s, n, scaled);
    }
    
    double expv = mean(xmr);
    double sdv = stddev(xmr);

    double z = (obs - expv) / sdv;

    double pv = R::pnorm(obs, expv, sdv, lower_tail, 0);

    arma::rowvec res = {obs, expv, sdv, z, pv};
    return(res);
  
}

arma::rowvec moransi(
        arma::vec x, 
        arma::sp_mat weight,
        arma::vec rowsumw,
        dqrng::xoshiro256plus rng,
        int permutation,
        double S1,
        double S2, 
        double s,
        int n,
        double ei,
        bool scaled = false,
	int lower_tail = 1
        ){

  double m = mean(x);
  arma::vec y = x - m;
  double cv = moranouterdot(y, weight); 
  double v = accu(pow(y, 2));
  double obs = (n/s) * (cv/v);
  
  if (scaled){
    double imax = (n/s) * (stddev(rowsumw % y)/sqrt(v/(n -1)));
    obs = obs/imax;
  }

  arma::rowvec res(5);
  if (permutation <= 10){
      res = cal_moransi_p_noperm(obs, ei, s, y, v, n, S1, S2, lower_tail);
  }else{
      res = cal_moransi_p_perm(x, weight, rowsumw, rng, permutation, s, n, scaled, lower_tail);
  }
  return(res);
}


struct RunMoransi : public Worker{
  const arma::sp_mat& x;
  const arma::sp_mat& weight;
  const arma::vec& rowsumw;
  simple_progress& p;
  const uint64_t seed;
  const int permutation;
  const double S1;
  const double S2;
  const double s;
  const int n;
  const double ei;
  const bool scaled;
  const int lower_tail;
  arma::mat& result;

  RunMoransi(const arma::sp_mat& x, const arma::sp_mat& weight, const arma::vec& rowsumw, simple_progress& p,
          const uint64_t seed, const int permutation, const double S1, const double S2, 
          const double s, const int n, const double ei, const bool scaled, const int lower_tail, mat& result):
      x(x), weight(weight), rowsumw(rowsumw), p(p), seed(seed), permutation(permutation), 
      S1(S1), S2(S2), s(s), n(n), ei(ei), scaled(scaled), lower_tail(lower_tail), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    dqrng::xoshiro256plus rng(seed);
    for (uword i = begin; i < end; i++){
        dqrng::xoshiro256plus lrng(rng);
        lrng.long_jump(i + 1);
        result.row(i) = moransi(x.col(i).as_dense(), weight, rowsumw, lrng, permutation, S1, S2, s, n, ei, scaled, lower_tail);
        p.increment();
    }
  }
};


// [[Rcpp::export]]
arma::mat CalMoransiParallel(arma::sp_mat& x, arma::sp_mat& wm, bool scaled = false, 
        int permutation = 999, int lower_tail=1){
  arma::sp_mat xm =  x.t();
  arma::sp_mat weight = wm.t();
  int n = x.n_rows;
  int m = x.n_cols;
  double ei = -(1.0 / (m - 1));
  
  arma::vec colsumw = rowsumsp(weight);
  arma::vec rowsumw = rowsumsp(wm);
  
  double S1 =  0.5 * accu(powsp(wm + wm.t()));
  double S2 = accu(pow(rowsumw + colsumw, 2.0));
  double s = accu(weight);

  Rcpp::IntegerVector seed(2, dqrng::R_random_int);
  uint64_t seed2 = dqrng::convert_seed<uint64_t>(seed);

  arma::mat result(n, 5);

  simple_progress p(n);

  RunMoransi runmoransi(xm, weight, rowsumw, p, seed2, permutation, S1, S2, s, m, ei, scaled, lower_tail, result);

  parallelFor(0, n, runmoransi);

  return(result);
}

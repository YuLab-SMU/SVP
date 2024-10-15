#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include "progress.h"
#include "autocorutils.h"
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::rowvec cal_getisord_p_noperm(
    arma::vec x,
    double obs, 
    double B0, 
    double B1,
    double B2,
    double B3,
    double B4,
    double n3,
    double ei,
    int lower_tail
  ){
  double m1 = accu(x);
  arma::vec tmp = x % x;
  double m2 = accu(tmp);
  tmp = tmp % x;
  double m3 = accu(tmp);
  tmp = tmp % x;
  double m4 = accu(tmp);

  double pm1 = pow(m1, 2.0);
  double eg2 = B0 * pow(m2, 2.0) + B1 * m4 + B2 * pm1 * m2 + B3 * m1 * m3 + B4 * pow(pm1, 2.0);
  double mm12 = pm1 - m2;
  eg2 = eg2/(pow(mm12, 2.0) * n3);

  double sdi = sqrt(eg2 - pow(ei, 2.0));

  double z = (obs - ei)/sdi;
  
  double pv = R::pnorm(obs, ei, sdi, lower_tail, 0);

  arma::rowvec res = {obs, ei, sdi, z, pv};
  return(res);
}

arma::rowvec getisord(
  arma::vec x, 
  arma::sp_mat weight,
  double B0,
  double B1,
  double B2,
  double B3,
  double B4,
  double n3,
  double ei,
  int lower_tail = 1
  ){

  double obs = cal_getisord(x, weight);

  arma::rowvec res = cal_getisord_p_noperm(x, obs, B0, B1, B2, B3, B4, n3, ei, lower_tail);

  return(res);
}


struct RunGetisOrd : public Worker{
  const arma::sp_mat& x;
  const arma::sp_mat& weight;
  simple_progress& p;
  const double B0;
  const double B1;
  const double B2;
  const double B3;
  const double B4;
  const double n3;
  const double ei;
  const int lower_tail;
  arma::mat& result;

  RunGetisOrd(const arma::sp_mat& x, const arma::sp_mat& weight, simple_progress& p, 
          const double B0, const double B1, const double B2, const double B3, 
          const double B4, const double n3, const double ei, const int lower_tail, 
          mat& result):
      x(x), weight(weight), p(p), B0(B0), B1(B1), B2(B2), B3(B3), B4(B4), n3(n3), ei(ei), 
      lower_tail(lower_tail), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        result.row(i) = getisord(x.col(i).as_dense(), weight, B0, B1, B2, B3, B4, n3, ei, lower_tail);
        p.increment();
    }
  }
};


// [[Rcpp::export]]
arma::mat CalGetisOrdParallel(arma::sp_mat& x, arma::sp_mat& wm, int lower_tail=1){
  arma::sp_mat xm = x.t();
  arma::sp_mat weight = wm.t();
  int m = x.n_rows;
  int n = x.n_cols;

  arma::vec colsumw = rowsumsp(weight);
  arma::vec rowsumw = rowsumsp(wm);
  
  double S1 =  0.5 * accu(powsp(wm + wm.t()));
  double S2 = accu(pow(rowsumw + colsumw, 2.0));
  double s = accu(weight);

  double ei = s/(n * (n - 1.0));

  double SS0 = pow(s, 2.0);
  double nn = pow(n, 2.0);
  double n3 = n * (n - 1.0) * (n - 2.0) * (n - 3.0);

  double B0 = (nn - 3.0 * n + 3.0) *S1 - n * S2 + 3.0 * SS0;
  double B1 = -((nn - n) * S1 - 2.0 * n * S2 + 6 * SS0);
  double B2 = -(2.0 * n * S1 - (n + 3.0) * S2 + 6.0 * SS0);
  double B3 = 4.0 * (n - 1.0) * S1 - 2.0 * (n + 1.0) * S2 + 8.0 * SS0;
  double B4 = S1 - S2 + SS0;

  arma::mat result(m, 5);

  simple_progress p(m);

  RunGetisOrd rungetisord(xm, weight, p, B0, B1, B2, B3, B4, n3, ei, lower_tail, result);

  parallelFor(0, m, rungetisord);

  return(result);
}

#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include "progress.h"
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat outermultidot(arma::rowvec x){
  arma::mat xm = repelem(x, x.n_elem, 1);
  xm.diag().fill(0.0);
  arma::mat res = xm % xm.t();
  return (res);
}

arma::rowvec cal_getisord_p_noperm(
    arma::rowvec x,
    double obs, 
    double ei, 
    double s, 
    int n, 
    double S1, 
    double S2,
    int lower_tail
  ){
  double SS0 = pow(s, 2.0);
  double nn = pow(n, 2.0);
  double n3 = n * (n - 1.0) * (n - 2.0) * (n - 3.0);
  
  double B0 = (nn - 3.0 * n + 3.0) *S1 - n * S2 + 3.0 * SS0;
  double B1 = -((nn - n) * S1 - 2.0 * n * S2 + 6 * SS0);
  double B2 = -(2.0 * n * S1 - (n + 3.0) * S2 + 6.0 * SS0);
  double B3 = 4.0 * (n - 1.0) * S1 - 2.0 * (n + 1.0) * S2 + 8.0 * SS0;
  double B4 = S1 - S2 + SS0;
  double m1 = accu(x);
  double m2 = accu(pow(x, 2.0));
  double m3 = accu(pow(x, 3.0));
  double m4 = accu(pow(x, 4.0));

  double eg2 = B0 * pow(m2, 2.0) + B1 * m4 + B2 * pow(m1, 2.0) * m2 + B3 * m1 * m3 + B4 * pow(m1, 4.0);
  double mm12 = pow(m1, 2.0) - m2;
  eg2 = eg2/(pow(mm12, 2.0) * n3);

  double sdi = sqrt(eg2 - pow(ei, 2.0));

  double z = (obs - ei)/sdi;
  
  double pv = R::pnorm(obs, ei, sdi, lower_tail, 0);

  arma::rowvec res = {obs, ei, sdi, z, pv};
  return(res);
}

double cal_getisord(
  arma::rowvec x,
  arma::mat weight
  ){
  arma::mat ym = outermultidot(x);
  double cv = accu(weight % ym);
  double v = accu(ym);
  double obs = cv/v;

  return(obs);
}

arma::rowvec getisord(
  arma::rowvec x, 
  arma::mat weight,
  double S1,
  double S2, 
  double s,
  int n,
  double ei,
  int lower_tail = 1
  ){

  double obs = cal_getisord(x, weight);

  arma::rowvec res(5);
  
  res = cal_getisord_p_noperm(x, obs, ei, s, n, S1, S2, lower_tail);

  return(res);
}


struct RunGetisOrd : public Worker{
  const arma::mat& x;
  const arma::mat& weight;
  simple_progress& p;
  const double S1;
  const double S2;
  const double s;
  const int n;
  const double ei;
  const int lower_tail;
  arma::mat& result;

  RunGetisOrd(const arma::mat& x, const arma::mat& weight, simple_progress& p, 
	  const double S1, const double S2, const double s, const int n, 
	  const double ei, const int lower_tail, mat& result):
      x(x), weight(weight), p(p), S1(S1), S2(S2), s(s), n(n), ei(ei), 
      lower_tail(lower_tail), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        result.row(i) = getisord(x.row(i), weight, S1, S2, s, n, ei, lower_tail);
        p.increment();
    }
  }
};


// [[Rcpp::export]]
arma::mat CalGetisOrdParallel(arma::sp_mat& x, arma::mat& weight, int lower_tail=1){
  arma::mat xm =  conv_to<arma::mat>::from(x);
  int n = x.n_rows;
  int m = x.n_cols;

  arma::rowvec colsumw = sum(weight, 0);
  arma::colvec rowsumw = sum(weight, 1);
  arma::rowvec rowsumw2 = conv_to<arma::rowvec>::from(rowsumw);
  
  double S1 =  0.5 * accu(pow(weight + weight.t(), 2.0));
  double S2 = accu(pow(rowsumw2 + colsumw, 2.0));
  double s = accu(weight);

  double ei = s/(m * (m - 1.0));

  arma::mat result(n, 5);

  simple_progress p(n);

  RunGetisOrd rungetisord(xm, weight, p, S1, S2, s, m, ei, lower_tail, result);

  parallelFor(0, n, rungetisord);

  return(result);
}

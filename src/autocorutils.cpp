#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat outerdot(arma::rowvec x){
  arma::mat xm = repelem(x, x.n_elem, 1);
  arma::mat res = xm % xm.t();
  return (res);
}

arma::mat outersubtractdot(arma::rowvec x){
  arma::mat xm = repelem(x, x.n_elem, 1);
  arma::mat res = pow(xm - xm.t(), 2.0);
  return (res);
}

arma::rowvec scaleCpp(arma::rowvec x){
  double f = sqrt(accu(pow(x, 2.0))/max(1.0, x.n_elem - 1.0));
  if (f > 0.0){
     x = x / f;
  }
  return(x);
}

arma::mat outermultidot(arma::rowvec x){
  arma::mat xm = repelem(x, x.n_elem, 1);
  xm.diag().fill(0.0);
  arma::mat res = xm % xm.t();
  return (res);
}

arma::vec lagCpp(arma::mat w, arma::vec x){
    arma::mat xm = repelem(x, 1, x.n_elem);
    xm.diag().fill(0.0);
    w = w % xm.t();
    arma::vec res = sum(w, 1);
    return(res);
}

double cal_moransi(
          arma::rowvec x,
          arma::mat weight,
          arma::rowvec rowsumw,
          double s,
          int n,
          bool scaled = false
        ){

    double m = mean(x);
    arma::rowvec y = x - m;
    arma::mat ym = outerdot(y);
    double cv = accu(weight % ym);
    double v = accu(pow(y, 2));
    double obs = (n/s) * (cv/v);

    if (scaled){
      double imax = (n/s) * (stddev(rowsumw % y)/sqrt(v/(n -1)));
      obs = obs/imax;
    }

    return(obs);
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

double cal_gearysc(
          arma::rowvec x,
          arma::mat weight,
          double s,
          int n
        ){

    double m = mean(x);
    arma::rowvec y = x - m;
    //y = scaleCpp(y);
    arma::mat ym = outersubtractdot(x);
    double cv = accu(weight % ym);
    double v = accu(pow(y, 2.0));
    double obs = ((n - 1.0) * cv)/(2.0 * s * v);

    return(obs);
}

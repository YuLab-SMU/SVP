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

arma::vec scaleCpp2(arma::vec x){
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

arma::vec lagCpp1(arma::mat w, arma::rowvec x){
    arma::mat xm = repelem(x, x.n_elem, 1);
    xm.diag().fill(0.0);
    w = w % xm;
    arma::vec res = sum(w, 1);
    return(res);
}

arma::rowvec lagCpp2(arma::mat w, arma::rowvec x){
    arma::vec res = lagCpp1(w, x);
    return(res.t());
}

//[[Rcpp::export]]
arma::vec cal_local_moran_bv(
    arma::vec x, 
    arma::vec y,
    arma::mat weight
    ){
    arma::vec ldy = lagCpp(weight, y);
    arma::vec res = x % ldy;
    return(res);
}

double cal_global_lee(
    arma::rowvec x, 
    arma::rowvec y, 
    arma::mat weight, 
    double S2,
    int n
    ){
    
    arma::rowvec dx = x - mean(x);
    arma::rowvec dy = y - mean(y);
    
    double dx2 = accu(pow(dx, 2.0));
    double dy2 = accu(pow(dy, 2.0));

    arma::rowvec ldx = lagCpp2(weight, dx);
    arma::rowvec ldy = lagCpp2(weight, dy);

    double L = (n/S2) * (sum(ldx % ldy))/(sqrt(dx2) * sqrt(dy2));
    
    return(L);
}

double cal_permutation_p(
        arma::vec x,
        double obs,
        int permutation,
        int alternative = 3
  ){
  double pv = 0.0;
  double n = 0.0;

  if (alternative == 1){
      n = 1.0 * abs(sum(obs >= x) - permutation/2.0 + 1.0)/permutation;
      pv = R::punif(n, 0.0, 0.5, 0, 0);
      return(pv);
  }else if (alternative == 3){
      n = 1.0 * sum(x > obs)/permutation;
      pv = R::punif(n, 0.0, 1.0, 1, 0);
  }else{
      n = 1.0 * sum(x > obs)/permutation;
      pv = R::punif(n, 0.0, 1.0, 0, 0);
  }
  return(pv);
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

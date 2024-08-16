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

arma::mat CalLocalGCpp(
    arma::vec x,
    arma::mat w,
    arma::vec wi,
    arma::vec S1i,
    int n
    ){

    double x_star = sum(x);
    arma::vec xx = pow(x, 2.0);
    arma::vec Gi = lagCpp(w, x);
    arma::vec xibar = (sum(x) - x) / (n - 1);
    arma::vec si2 = (sum(xx) - xx)/(n - 1) - pow(xibar, 2.0);
    arma::vec EG = wi % xibar;
    arma::vec ZG = Gi - EG;
    arma::vec VG = si2 % (((n - 1) * S1i - pow(wi, 2.0))/(n - 2));
    ZG = ZG/sqrt(VG);

    arma::vec scale = x_star - x;
    arma::mat res(n, 5);

    res.col(0) = Gi/scale;
    res.col(1) = EG/scale;
    res.col(2) = VG/(pow(scale, 2.0));
    res.col(3) = ZG;
    res.col(4) = x;

    return(res);
}

arma::mat tidylocalg(
            arma::vec res1,
            arma::vec res2,
            arma::vec res3,
            arma::vec res4,
            arma::vec res5,
            arma::mat res
        ){
    res.col(0) = res1;
    res.col(1) = res2;
    res.col(2) = res3;
    res.col(3) = res4;
    res.col(4) = res5;
    return(res);
}

arma::mat CalLocalMoranCpp(
    arma::vec x,
    arma::mat w,
    arma::vec wi,
    arma::vec Wi2,
    int n
    ){

    arma::mat res(n, 6);
    arma::vec z = x - mean(x);
    arma::vec lz = lagCpp(w, z);
    arma::vec zz2 = pow(z, 2.0);
    double m2 = sum(zz2)/n;
    arma::vec I = (z/m2) % lz;
    arma::vec EI = -(zz2 % wi) / ((n - 1.0) * m2);
    arma::vec wwi2 = pow(wi, 2.0);
    arma::vec VI = pow(z/m2, 2.0) * (n / (n - 2.0));
    VI = VI % (Wi2 - (wwi2/ (n- 1.0)));
    VI = VI % (m2 - (zz2/ (n - 1.0)));
    arma::vec ZI = (I - EI)/sqrt(VI);

    res.col(0) = I;
    res.col(1) = EI;
    res.col(2) = VI;
    res.col(3) = ZI;
    res.col(4) = z;
    res.col(5) = lz;

    return(res);
}

arma::mat tidylocalmoran(
        arma::vec res1,
        arma::vec res2,
        arma::vec res3,
        arma::vec res4,
        arma::vec res5,
        arma::vec res6,
        arma::mat res
        ){
    res.col(0) = res1;
    res.col(1) = res2;
    res.col(2) = res3;
    res.col(3) = res4;
    res.col(4) = res5;
    res.col(5) = res6;
    return(res);
}

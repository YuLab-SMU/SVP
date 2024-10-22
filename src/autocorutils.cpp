#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


double moranouterdot(arma::vec x, arma::sp_mat w){
  double res = 0.0;
  for (size_t i = 0; i < w.n_cols; i++){
    res += accu(x(i) * x % w.col(i));
  }
  return(res);
}

double gearyouterdot(arma::vec x, arma::sp_mat w){
  double res = 0.0;
  for (size_t i = 0; i < w.n_cols; i++){
    res += accu(pow((x(i) - x), 2.0) % w.col(i));
  }
  return(res);
}

arma::vec getisordouterdot(arma::vec x, arma::sp_mat w){
  arma::vec res(2);
  for (size_t i = 0; i < w.n_cols; i++){
    res(0) += accu(x(i) * x) - pow(x(i), 2.0);
    res(1) += accu(x(i) * x  % w.col(i)) - pow(x(i),2.0) * w(i,i);
  }
  return(res);
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

arma::vec lagCpp(arma::sp_mat w, arma::vec x){
    arma::vec res(x.n_elem);
    for (size_t i = 0; i < w.n_cols; i++){
        res(i) = accu(x % w.col(i));
    }
    return(res);
}

double lagCpp2(arma::sp_mat w, arma::vec x, arma::vec y){
    double res = 0.0;
    for (size_t i = 0; i < w.n_cols; i++){
        res += accu(x % w.col(i)) * accu(y % w.col(i));
    }
    return(res);
}

arma::vec lagCpp3(arma::sp_mat w, arma::vec x, arma::vec y){
    arma::vec res(x.n_elem);
    for (size_t i = 0; i < w.n_cols; i++){
        res(i) = accu(x % w.col(i)) * accu(y % w.col(i));
    }
    return(res);
}

//[[Rcpp::export]]
arma::vec cal_local_moran_bv(
    arma::vec x, 
    arma::vec y,
    arma::sp_mat weight
    ){
    arma::vec res = x % lagCpp(weight, y);
    return(res);
}

double cal_global_lee(
    arma::vec x, 
    arma::vec y, 
    arma::sp_mat weight, 
    double S2,
    int n
    ){
    
    double dx2 = accu(pow(x - mean(x), 2.0));
    double dy2 = accu(pow(y - mean(y), 2.0));

    double ldxy = lagCpp2(weight, x-mean(x), y-mean(y));

    double L = (n/S2) * ldxy/(sqrt(dx2) * sqrt(dy2));
    
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
          arma::vec x,
          arma::sp_mat weight,
          arma::vec rowsumw,
          double s,
          int n,
          bool scaled = false
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

    return(obs);
}

double cal_getisord(
  arma::vec x,
  arma::sp_mat weight
  ){
  arma::vec tmp = getisordouterdot(x, weight);
  double obs = tmp(1)/tmp(0);

  return(obs);
}

double cal_gearysc(
          arma::vec x,
          arma::sp_mat weight,
          double s,
          int n
        ){

    double m = mean(x);
    arma::vec y = x - m;
    //y = scaleCpp2(y);
    double cv = gearyouterdot(y, weight);
    double v = accu(pow(y, 2.0));
    double obs = ((n - 1.0) * cv)/(2.0 * s * v);

    return(obs);
}

arma::mat CalLocalGCpp(
    arma::vec x,
    arma::sp_mat w,
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
    arma::sp_mat w,
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

arma::vec rowsumsp(arma::sp_mat x, bool flag = false){
   if (flag){
     for (auto it = x.begin(); it != x.end(); ++it) {
       *it = std::pow(*it, 2.0);
     }
   }
   arma::vec res = x * ones(x.n_cols, 1);
   return(res);
}

arma::sp_mat powsp(arma::sp_mat x){
  for (auto it = x.begin(); it != x.end(); ++it) {
    *it = std::pow(*it, 2.0);
  }
  return (x);
}

double calculateF1(const arma::vec& predictions, const arma::vec& actuals){
    double TP = arma::sum((predictions == 1) % (actuals == 1));
    double FP = arma::sum((predictions == 1) % (actuals == 0));
    double FN = arma::sum((predictions == 0) % (actuals == 1));
    double precision = TP / (TP + FP);
    double recall = TP / (TP + FN);
    double F1 = 2 * (precision * recall) / (precision + recall);

    return (F1);
}

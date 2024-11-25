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

arma::vec cal_quant(
  arma::vec diagM,
  arma::vec diagMt,
  arma::vec MII
  ){
  double Fon0 = accu(diagM);
  double Foff0 = accu(MII) - Fon0;
  double Fon1 = accu(pow(diagM, 2.0));
  double Foff1 = accu(diagMt) - Fon1;
  double Foff2 = accu(pow(MII - diagM, 2.0));
  double Fall2 = accu(pow(MII, 2.0));
  arma::vec res = {Foff0, Fon0, Foff1, Fon1, Foff2, Fall2};
  return(res);
}

arma::vec cal_Qquant(
  arma::vec x,
  arma::vec y,
  int n
){
  arma::vec xbar = x - mean(x);
  arma::vec ybar = y - mean(y);
  double t = ((n - 1.0)/n) * arma::stddev(x) * arma::stddev(y);
  arma::vec diagM = xbar % ybar/t;
  arma::vec MII(n);
  arma::vec diagMt(n);
  for (int i = 0; i < n; i++){
    MII(i) = accu(xbar(i) * ybar/t);
    diagMt(i) = accu(pow(ybar(i) * xbar/t, 2.0));
  }
  arma::vec res = cal_quant(diagM, diagMt, MII);
  return(res);
}

arma::vec cal_EGamma(
  arma::vec P,
  arma::vec Q,
  int n
  ){
  double EGammaoff = P(0) * Q(0)/(n * (n - 1.0));
  double EGammaon = P(1) * Q(1)/n;
  arma::vec res = {EGammaoff, EGammaon};
  return(res);
}

arma::vec cal_VarGamma(
  arma::vec P,
  arma::vec Q,
  arma::vec EG,
  int n
  ){
  double varGammaoff = 2.0 * P(2) * Q(2)/(n * (n - 1.0));
  varGammaoff = varGammaoff + 4.0 * (P(4) - P(2)) * (Q(4) - Q(2))/(n * (n - 1.0) * (n - 2.0));
  varGammaoff = varGammaoff +
                (pow(P(0), 2.0) + 2.0 * P(2) - 4.0 * P(4)) *
                (pow(Q(0), 2.0) + 2.0 * Q(2) - 4.0 * Q(4))/(n * (n - 1.0) * (n - 2.0) * (n - 3.0));
  varGammaoff = varGammaoff - pow(EG(0), 2.0);

  double varGammaon = 1.0 * P(3) * Q(3)/n;
  varGammaon = varGammaon + (pow(P(1), 2.0) - P(3)) * (pow(Q(1), 2.0) - Q(3))/(n * (n - 1.0));
  varGammaon = varGammaon - pow(EG(1), 2.0);

  double varGammaonoff = (P(5) - P(3) - P(4)) * (Q(5) - Q(3) - Q(4))/(2 * n * (n - 1.0));
  varGammaonoff = varGammaonoff +
                  (P(1) * P(0) - (P(5) - P(3) - P(4))) *
                  (Q(1) * Q(0) - (Q(5) - Q(3) - Q(4)))/(n * (n - 1.0) * (n - 2.0));
  varGammaonoff = varGammaonoff - EG(0) * EG(1);

  arma::vec res = {varGammaon, varGammaoff, varGammaonoff};
  return(res);
}

double cal_pval_lee(
  double L, 
  arma::vec dEG, 
  arma::vec dVarG, 
  int alternative=1
  ){
  double EL = dEG(1) + dEG(0);
  double VL = sqrt(dVarG(0) + dVarG(1) + 2.0 * dVarG(2));
  double ZL = (L - EL)/VL;
  double pval = 0.0;
  if (alternative == 3){
    pval = R::pnorm(ZL, 0, 1, 0, 0);
  }else if (alternative == 1){
    pval = 2 * R::pnorm(std::abs(ZL), 0, 1, 0, 0);
  }else if (alternative == 2){
    pval = R::pnorm(ZL, 0, 1, 1, 0);
  }
  return(pval);
}

double cal_pval_lee_pipeline(
  double L,
  arma::vec x,
  arma::vec y,
  arma::vec P,
  int n,
  int alternative = 1
  ){
  arma::vec Q = cal_Qquant(x, y, n);
  arma::vec dEG = cal_EGamma(P, Q, n);
  arma::vec dVarG = cal_VarGamma(P, Q, dEG, n);
  double pval = cal_pval_lee(L, dEG, dVarG, alternative);
  return(pval);
}

#include <RcppArmadillo.h>
using namespace arma;

#ifndef autocorutils_H
#define autocorutils_H

double moranouterdot(arma::vec x, arma::sp_mat w);

arma::rowvec scaleCpp(arma::rowvec x);
arma::vec scaleCpp2(arma::vec x);

arma::vec lagCpp(arma::sp_mat w, arma::vec x);

arma::vec lagCpp3(arma::sp_mat w, arma::vec x, arma::vec y);

arma::vec rowsumsp(arma::sp_mat x, bool flag = false);

arma::sp_mat powsp(arma::sp_mat x);

double calculateF1(const arma::vec& predictions, const arma::vec& actuals);

double cal_moransi(arma::vec x, arma::sp_mat weight, 
                   arma::vec rowsumw, double s, 
                   int n, bool scaled = false); 

double cal_getisord(arma::vec x, arma::sp_mat weight);

double cal_permutation_p(
   arma::vec x,
   double obs,
   int permutation,
   int alternative = 3
);

double cal_gearysc(
   arma::vec x, 
   arma::sp_mat weight, 
   double s, 
   int n
);

double cal_global_lee(
  arma::vec x, 
  arma::vec y, 
  arma::sp_mat weight, 
  double S2, 
  int n
);

arma::vec cal_local_moran_bv(
    arma::vec x,
    arma::vec y,
    arma::sp_mat weight
);

arma::mat CalLocalGCpp(
    arma::vec x,
    arma::sp_mat w,
    arma::vec wi,
    arma::vec S1i,
    int n
);

arma::mat tidylocalg(
    arma::vec res1,
    arma::vec res2,
    arma::vec res3,
    arma::vec res4,
    arma::vec res5,
    arma::mat res
);

arma::mat CalLocalMoranCpp(
    arma::vec x,
    arma::sp_mat w,
    arma::vec wi,
    arma::vec Wi2,
    int n
);

arma::mat tidylocalmoran(
    arma::vec res1,
    arma::vec res2,
    arma::vec res3,
    arma::vec res4,
    arma::vec res5,
    arma::vec res6,
    arma::mat res
);
#endif

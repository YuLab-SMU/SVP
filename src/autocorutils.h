#include <RcppArmadillo.h>
using namespace arma;

#ifndef autocorutils_H
#define autocorutils_H

arma::mat outerdot(arma::rowvec x);

arma::mat outersubtractdot(arma::rowvec x);

arma::mat outermultidot(arma::rowvec x); 

arma::rowvec scaleCpp(arma::rowvec x);
arma::vec scaleCpp2(arma::vec x);

arma::vec lagCpp(arma::mat w, arma::vec x);

double cal_moransi(arma::rowvec x, arma::mat weight, 
                   arma::rowvec rowsumw, double s, 
                   int n, bool scaled = false); 

double cal_getisord(arma::rowvec x, arma::mat weight);

double cal_permutation_p(
   arma::vec x,
   double obs,
   int permutation,
   int alternative = 3
);

double cal_gearysc(
   arma::rowvec x, 
   arma::mat weight, 
   double s, 
   int n
);

double cal_global_lee(
  arma::rowvec x, 
  arma::rowvec y, 
  arma::mat weight, 
  double S2, 
  int n
);

arma::vec cal_local_moran_bv(
    arma::vec x,
    arma::vec y,
    arma::mat weight
);

#endif

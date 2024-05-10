#include <RcppArmadillo.h>
using namespace arma;

#ifndef autocorutils_H
#define autocorutils_H

arma::mat outerdot(arma::rowvec x);

arma::mat outersubtractdot(arma::rowvec x);

arma::mat outermultidot(arma::rowvec x); 

arma::mat scaleCpp(arma::rowvec x);

arma::vec lagCpp(arma::mat w, arma::vec x);

double cal_moransi(arma::rowvec x, arma::mat weight, 
                   arma::rowvec rowsumw, double s, 
                   int n, bool scaled = false); 

double cal_getisord(arma::rowvec x, arma::mat weight);

double cal_gearysc(arma::rowvec x, arma::mat weight, double s, int n);

#endif

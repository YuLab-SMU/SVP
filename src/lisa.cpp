#include <RcppArmadillo.h>
#include "autocorutils.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
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


// [[Rcpp::export]]
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

#include <RcppParallel.h>
#include <RcppArmadillo.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat outersubtractdot(arma::rowvec x){
  arma::mat xm = repelem(x, x.n_elem, 1);
  arma::mat res = pow(xm - xm.t(), 2);
  return (res);
}

double cal_gearysc(
          arma::rowvec x,
          arma::mat weight,
          double s,
          int n
        ){

    double m = mean(x);
    arma::rowvec y = x - m;
    arma::mat ym = outersubtractdot(x);
    double cv = accu(weight % ym);
    double v = accu(pow(y, 2));
    double obs = ((n - 1) * cv)/(s * v);

    return(obs);
}


arma::rowvec cal_gearysc_p_perm(
        arma::rowvec x,
        arma::mat weight,
        double s,
        int n,
        int permutation = 999
        ){
    double obs = cal_gearysc(x, weight, s, n);
    arma::vec xmr = arma::vec(permutation, arma::fill::zeros);
    for (int i = 0; i < permutation; i++){
        xmr(i) = cal_gearysc(arma::shuffle(x), weight, s, n);
    }
    
    double expv = mean(xmr);
    double sdv = stddev(xmr);

    double pv = R::pnorm(obs, expv, sdv, 1, 0);

    arma::rowvec res = {obs, expv, sdv, pv};
    return(res);
  
}

struct RunGearysc : public Worker{
  const arma::mat& x;
  const arma::mat& weight;
  const double s;
  const int n;
  const int permutation;
  arma::mat& result;

  RunGearysc(const arma::mat& x, const arma::mat& weight, const double s, const int n,
             const int permutation, mat& result):
      x(x), weight(weight), s(s), n(n), permutation(permutation), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        result.row(i) = cal_gearysc_p_perm(x.row(i), weight, s, n, permutation);
    }
  }
};


// [[Rcpp::export]]
arma::mat CalGearyscParallel(arma::sp_mat& x, arma::mat& weight, int permutation = 999){
  arma::mat xm =  conv_to<arma::mat>::from(x);
  int n = x.n_rows;
  int m = x.n_cols;
 
  weight.diag().fill(0);
  arma::colvec rowsumw = sum(weight, 1);
  rowsumw.replace(0, 1);
  weight = weight.each_col() / rowsumw;

  double s = accu(weight);

  arma::mat result(n, 4);
  RunGearysc rungearysc(xm, weight, s, m, permutation, result);

  parallelFor(0, n, rungearysc);

  return(result);

}

#include <RcppParallel.h>
#include <RcppArmadillo.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat outerdot(arma::rowvec x){
  arma::mat xm = repelem(x, x.n_elem, 1);
  arma::mat res = xm % xm.t();
  return (res);
}

arma::rowvec cal_moransi_p_noperm(
        double obs, 
        double ei, 
        double s, 
        arma::rowvec y, 
        double v, 
        int n, 
        double S1, 
        double S2
  ){
  double s_sq = pow(s, 2);

  double k = (accu(pow(y,4))/n)/pow(v/n, 2);
  double sdi = sqrt((n * ((pow(n, 2) - 3 * n + 3) * S1 - n * S2 + 3 * s_sq) -
              k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s_sq))/
          ((n - 1) * (n - 2) * (n - 3) * s_sq) - 1/(pow(n - 1, 2)));

  double pv = R::pnorm(obs, ei, sdi, 1, 0);
  arma::rowvec res = {obs, ei, sdi, pv};
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


arma::rowvec cal_moransi_p_perm(
        arma::rowvec x,
        arma::mat weight,
        arma::rowvec rowsumw,
        double s,
        int n,
        bool scaled = false,
        int permutation = 999
        ){
    double obs = cal_moransi(x, weight, rowsumw, s, n, scaled);
    arma::vec xmr = arma::vec(permutation, arma::fill::zeros);
    for (int i = 0; i < permutation; i++){
        xmr(i) = cal_moransi(arma::shuffle(x), weight, rowsumw, s, n, scaled);
    }
    
    double expv = mean(xmr);
    double sdv = stddev(xmr);

    double pv = R::pnorm(obs, expv, sdv, 1, 0);

    arma::rowvec res = {obs, expv, sdv, pv};
    return(res);
  
}

arma::rowvec moransi(
        arma::rowvec x, 
        arma::mat weight,
        arma::rowvec rowsumw,
        double S1,
        double S2, 
        double s,
        int n,
        double ei,
        bool scaled = false,
        int permutation = 999
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

  arma::rowvec res(4);
  if (permutation <= 10){
      res = cal_moransi_p_noperm(obs, ei, s, y, v, n, S1, S2);
  }else{
      res = cal_moransi_p_perm(x, weight, rowsumw, s, n, scaled, permutation);
  }
  return(res);
}


struct RunMoransi : public Worker{
  const arma::mat& x;
  const arma::mat& weight;
  const arma::rowvec& rowsumw;
  const double S1;
  const double S2;
  const double s;
  const int n;
  const double ei;
  const bool scaled;
  const int permutation;
  arma::mat& result;

  RunMoransi(const arma::mat& x, const arma::mat& weight, const arma::rowvec& rowsumw,
          const double S1, const double S2, const double s, const int n, const double ei,
          const bool scaled, const int permutation, mat& result):
      x(x), weight(weight), rowsumw(rowsumw), S1(S1), S2(S2), s(s), n(n), ei(ei), 
      scaled(scaled), permutation(permutation), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        result.row(i) = moransi(x.row(i), weight, rowsumw, S1, S2, s, n, ei, scaled, permutation);
    }
  }
};


// [[Rcpp::export]]
arma::mat CalMoransiParallel(arma::sp_mat& x, arma::mat& weight, bool scaled = false, int permutation = 999){
  arma::mat xm =  conv_to<arma::mat>::from(x);
  int n = x.n_rows;
  int m = x.n_cols;
  double ei = -(1.0 / (m - 1));
  
  arma::colvec rowsumw = sum(weight, 1);
  rowsumw.replace(0, 1);
  weight = weight.each_col() / rowsumw;

  arma::rowvec colsumw = sum(weight, 0);
  rowsumw = sum(weight, 1);
  arma::rowvec rowsumw2 = conv_to<arma::rowvec>::from(rowsumw);
  
  double S1 =  0.5 * accu(square(weight + weight.t()));
  double S2 = accu(pow(rowsumw2 + colsumw, 2));
  double s = accu(weight);

  arma::mat result(n, 4);
  RunMoransi runmoransi(xm, weight, rowsumw2, S1, S2, s, m, ei, scaled, permutation, result);

  parallelFor(0, n, runmoransi);

  return(result);

}

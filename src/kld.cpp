#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;


NumericVector Quantile(NumericVector x, NumericVector probs) {
  const size_t n=x.size(), np=probs.size();
  if (n==0) return x;
  if (np==0) return probs;
  NumericVector index = (n-1.)*probs, y=x.sort(), x_hi(np), qs(np);
  NumericVector lo = Rcpp::floor(index), hi = Rcpp::ceiling(index);

  for (size_t i=0; i<np; ++i) {
    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];
    if ((index[i]>lo[i]) && (x_hi[i] != qs[i])) {
      double h;
      h = index[i]-lo[i];
      qs[i] = (1.-h)*qs[i] + h*x_hi[i];
    }
  }
  return qs;
}

double BandwidthNrdCpp(NumericVector x){
    NumericVector p = {0.25, 0.75};
    ////the arma::quantile might have bug.
    ////vec r = arma::quantile(as<arma::vec>(x), p);
    NumericVector r = Quantile(x, {0.25, 0.75});
    double h = (r[1] - r[0])/1.34;
    double w = pow(x.length(), -0.2);
    double s = sqrt(var(x));
    double v = 4 * 1.06 * std::min(s, h) * w;
    return (v);
}

arma::uvec findIntervalCpp(arma::vec x, arma::vec breaks) {
  uvec out(x.size());

  vec::iterator it, pos;
  uvec::iterator out_it;

  for(it = x.begin(), out_it = out.begin(); it != x.end();
      ++it, ++out_it) {
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    *out_it = std::distance(breaks.begin(), pos);
  }
  return (out);
}

NumericVector extractDensity(NumericMatrix x, NumericVector gx, NumericVector gy, arma::mat z){
  NumericVector coordx = x( _ , 0);
  NumericVector coordy = x( _ , 1);
  arma::vec oldx = as<arma::vec>(coordx);
  arma::vec oldy = as<arma::vec>(coordy);
  arma::uvec newx = findIntervalCpp(oldx, as<arma::vec>(gx));
  arma::uvec newy = findIntervalCpp(oldy, as<arma::vec>(gy));
  arma::mat sz = z.submat(newx, newy);
  arma::vec tmp = sz.diag();
  NumericVector res = as<NumericVector>(wrap(tmp));
  return (res);
}

NumericVector Kde2dWeightedCpp(NumericMatrix x,
                   NumericVector w, 
		   NumericVector gx,
		   NumericVector gy,
		   NumericVector h
                   ){
    
    int nx = x.nrow();
    int n = gx.length();
    w = w/sum(w) * w.length();

    NumericMatrix ax = outer(gx, x( _ , 0), std::minus<double>());
    NumericMatrix ay = outer(gy, x( _ , 1), std::minus<double>());    

    ax = ax / h[0];
    ay = ay / h[1];

    NumericVector v = rep_each(w, n);
    NumericVector dax = Rcpp::dnorm(as<NumericVector>(ax));

    NumericVector day = Rcpp::dnorm(as<NumericVector>(ay)) * v;
    day.attr("dim") = Dimension(n, nx);
    NumericMatrix daym = as<NumericMatrix>(day);

    daym = transpose(daym);

    v = v * dax;
    v.attr("dim") = Dimension(n, nx);
    
    NumericMatrix u = as<NumericMatrix>(v);

    arma::mat z = (as<arma::mat>(u) * as<arma::mat>(daym))/(sum(v) * h[0] * h[1]);
    
    //NumericMatrix tmp = as<NumericMatrix>(wrap(z));
    //NumericVector res = as<NumericVector>(tmp);

    NumericVector res = extractDensity(x, gx, gy, z);
    return (res);
}

//[[Rcpp::export]]
NumericVector CalRandSpatialKld(
        NumericMatrix coords, 
        NumericVector w, 
	NumericVector gx,
	NumericVector gy,
	NumericVector h,
        NumericVector bg,
        int random_times = 200,
	unsigned int seed = 1024
        ){
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    NumericVector bootkld(random_times);

    set_seed_r(seed);

    for (int j = 0; j < random_times; j++){
        NumericVector s = sample(w, w.size());
        NumericVector z = Kde2dWeightedCpp(coords, s, gx, gy, h);
        bootkld[j] = log(sum(z * log(z / bg)));
    }
    return(bootkld);
}

NumericVector CalKldPvalue(NumericVector boot, double x){
    boot = boot[is_finite(boot)];
    double bmean = mean(boot);
    double bsd = sd(boot);
    double pval = R::pnorm(x, bmean, bsd, 0, 0);
    NumericVector res = {x, bmean, bsd, pval};
    
    return (res);
}

//' Compute Background 2D Kernel Density
//' @param coords coordinate matrix.
//' @param gx Vector grid points in x direction, see(seq(lims[1], lims[2], length.out=200)).
//' @param gy Vector grid points in y direction, see(seq(lims[3], lims[4], length.out=200)).
//' @param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
//' (see MASS::bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
// [[Rcpp::export]]
NumericVector CalBgSpatialKld(NumericMatrix coords, 
	NumericVector gx, 
	NumericVector gy, 
        NumericVector h){
    double prop = 1.0 / coords.nrow();
    NumericVector bgw = rep(prop, coords.nrow());
    NumericVector bgm = Kde2dWeightedCpp(coords, bgw, gx, gy, h);
    NumericVector bgkld = bgm + 1e-300;
    return (bgkld);
}

//' Compute the Kullback–Leibler Divergence using 2D Kernel Density Estimation 
//' With Weighted and Statistical Test With Permutation for single weight vector.
//' @param coords coordinate matrix.
//' @param d the weight vector (the expression of gene or score of pathway).
//' @param bgkld the kernel density of background (the result of CalBgSpatialKld).
//' @param gx Vector grid points in x direction, see(seq(lims[1], lims[2], length.out=100)).
//' @param gy Vector grid points in y direction, see(seq(lims[3], lims[4], length.out=100)).
//' @param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
//' (see bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
//' @param random_times the permutation numbers for each weight to test whether 
//' it is significantly, default is 200.
//' @param seed The random seed to use to evaluate, default 123.
// [[Rcpp::export]]
NumericVector CalSpatialKld(NumericMatrix coords, 
                            NumericVector d,
                            NumericVector bgkld, 
			    NumericVector gx,
			    NumericVector gy,
			    NumericVector h,
                            int random_times = 200,
			    unsigned int seed = 123
			    ){

    NumericVector z = Kde2dWeightedCpp(coords, d, gx, gy, h);
    double kld = log(sum(z * log(z / bgkld)));
    if (! std::isnormal(kld)){
        kld = 0.0;
    }

    NumericVector bootkld = CalRandSpatialKld(coords, d, gx, gy, h, bgkld, random_times, seed);
    
    NumericVector res = CalKldPvalue(bootkld, kld);
    return(res);

}

//' Compute the Kullback–Leibler Divergence using 2D Kernel Density Estimation 
//' With Weighted And Statistical Test With Permutation.
//' @param coords coordinate matrix.
//' @param d matrix (the expression of gene or score of pathway).
//' @param l The limits of the rectangle covered by the grid as c(xl, xu, yl, yu).
//' @param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
//' (see bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
//' @param n the Number of grid points in each direction, default is 100.
//' @param random_times the permutation numbers for each weight to test whether
//' it is significantly, default is 200.
//' @param seed The random seed to use to evaluate, default 123.
// [[Rcpp::export]]
NumericMatrix CalSpatialKldCpp(NumericMatrix coords, arma::sp_mat d, NumericVector l, Nullable<NumericVector> h, 
	int n = 100, int random_times = 200, unsigned int seed = 123){

    NumericMatrix w = as<NumericMatrix>(wrap(d));

    NumericVector gx = wrap(arma::linspace(l[0], l[1], n));
    NumericVector gy = wrap(arma::linspace(l[2], l[3], n));

    NumericVector H(coords.ncol());
    if (h.isNull()){
      for (int j=0; j < coords.ncol();j ++){
         H[j] = BandwidthNrdCpp(coords(_,j)) / 4;
      }
    }else{
      H = h;
    }
    
    NumericVector bgkld = CalBgSpatialKld(coords, gx, gy, H);

    int num = w.nrow();

    NumericMatrix result(num, 4);
    for (int i = 0; i < num; i++){
       result(i, _ ) = CalSpatialKld(coords, w( i, _ ), bgkld, gx, gy, H, random_times, seed);
    }

    return (result);
}


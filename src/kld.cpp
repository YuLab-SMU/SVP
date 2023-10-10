#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

double BandwidthNrdCpp(NumericVector x){
    vec p = {0.25, 0.75};
    vec r = arma::quantile(as<arma::vec>(x), p);
    //NumericVector r = Quantile(x, {0.25, 0.75});
    double h = (r[1] - r[0])/1.34;
    double w = pow(x.length(), -0.2);
    double s = sqrt(var(x));
    double v = 4 * 1.06 * std::min(s, h) * w;
    return (v);
}

NumericMatrix Kde2dWeightedCpp(NumericMatrix x,
                   NumericVector w, 
                   int n, 
                   Nullable<NumericVector> lims){
    
    int nx = x.nrow();
    NumericVector l;
    
    if (lims.isNotNull()){
        l = lims;
    }else{
        NumericVector t1 = range(x( _ , 0));
        NumericVector t2 = range(x( _ , 1));
        l = {t1[0], t1[1], t2[0], t2[1]};
    }

    double h1 = BandwidthNrdCpp(x( _ , 0)) / 4;
    double h2 = BandwidthNrdCpp(x( _ , 1)) / 4;
    
    NumericVector gx = wrap(arma::linspace(l[0], l[1], n));
    NumericVector gy = wrap(arma::linspace(l[2], l[3], n));

    NumericMatrix ax = outer(gx, x( _ , 0), std::minus<double>());
    NumericMatrix ay = outer(gy, x( _ , 1), std::minus<double>());    

    ax = ax / h1;
    ay = ay / h2;

    NumericVector v = rep_each(w, n);
    NumericVector dax = Rcpp::dnorm(as<NumericVector>(ax));

    NumericVector day = Rcpp::dnorm(as<NumericVector>(ay));
    day.attr("dim") = Dimension(n, nx);
    NumericMatrix daym = as<NumericMatrix>(day);

    daym = transpose(daym) / (sum(w) * h1 * h2);

    v = v * dax;
    v.attr("dim") = Dimension(n, nx);
    
    NumericMatrix u = as<NumericMatrix>(v);

    arma::mat z = as<arma::mat>(u) * as<arma::mat>(daym);

    return wrap(z);
}

double CalRandSpatialKld(
        NumericMatrix coords, 
        NumericVector w, 
        NumericVector bg, 
        int n = 10){
    NumericVector s = sample(w, w.size());
    NumericVector z = as<NumericVector>(Kde2dWeightedCpp(coords, s, n, R_NilValue));
    double res = log(sum(z * log(z / bg + 1.0)));
    return (res);
}

NumericVector CalKldPvalue(NumericVector boot, double x){
    double bmean = mean(boot);
    double bsd = sd(boot);
    double pval = R::pnorm(x, bmean, bsd, 0, 0);
    NumericVector res = {x, bmean, bsd, pval};
    
    return (res);
}

//' Compute Background 2D Kernel Density
//' @param coords coordinate matrix.
//' @param n the Number of grid points in each direction, default is 25.
// [[Rcpp::export]]
NumericVector CalBgSpatialKld(NumericMatrix coords, int n = 25){
    double prop = 1.0 / coords.nrow();
    NumericVector bgw = rep(prop, coords.nrow());
    NumericMatrix bgm = Kde2dWeightedCpp(coords, bgw, n, R_NilValue);
    NumericVector bgkld = as<NumericVector>(bgm) + 1e-200;
    return (bgkld);
}

//' Compute the Kullback–Leibler Divergence using 2D Kernel Density Estimation 
//' With Weighted and Statistical Test With Permutation for single weight vector.
//' @param coords coordinate matrix.
//' @param d the weight vector (the expression of gene or score of pathway).
//' @param bgkld the kernel density of background (the result of CalBgSpatialKld).
//' @param n the Number of grid points in each direction, default is 25.
//' @param random_times the permutation numbers for each weight to test whether 
//' it is significantly, default is 999.
// [[Rcpp::export]]
NumericVector CalSpatialKld(NumericMatrix coords, 
                            NumericVector d,
                            NumericVector bgkld, 
                            int n = 25, 
                            int random_times = 999){

    NumericVector k = as<NumericVector>(Kde2dWeightedCpp(coords, d, n, R_NilValue));
    double kld = log(sum(k * log(k / bgkld + 1.0)));

    NumericVector bootkld(random_times);

    for (int j = 0; j < random_times; j++){
        bootkld[j] = CalRandSpatialKld(coords, d, bgkld, n);
    }
    
    NumericVector res = CalKldPvalue(bootkld, kld);
    return(res);

}

//' Compute the Kullback–Leibler Divergence using 2D Kernel Density Estimation 
//' With Weighted And Statistical Test With Permutation.
//' @param coords coordinate matrix.
//' @param d matrix (the expression of gene or score of pathway).
//' @param n the Number of grid points in each direction, default is 25.
//' @param random_times the permutation numbers for each weight to test whether
//' it is significantly, default is 999.
// [[Rcpp::export]]
NumericMatrix CalSpatialKldCpp(NumericMatrix coords, NumericMatrix d, int n = 25, int random_times = 999){
    
    NumericVector bgkld = CalBgSpatialKld(coords, n);

    int num = d.nrow();

    NumericMatrix result(num, 4);
    for (int i = 0; i < num; i++){
       result(i, _ ) = CalSpatialKld(coords, d( i, _ ), bgkld, n, random_times);
    }

    return (result);
}


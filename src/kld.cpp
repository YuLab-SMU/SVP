#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <algorithm>
#include <xoshiro.h>
#include <convert_seed.h>
#include <R_randgen.h>
#include "buildrand.h"
#include "progress.h"
using namespace RcppParallel;
using namespace arma;
using namespace Rcpp;
using namespace std;

arma::vec Quantile(arma::vec x, arma::vec probs) {
  const size_t n=x.n_elem, np=probs.n_elem;
  if (n==0) return x;
  if (np==0) return probs;
  arma::vec index = (n-1.0)*probs, y=sort(x), x_hi(np), qs(np);
  arma::vec lo = arma::floor(index), hi = arma::ceil(index);

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

double BandwidthNrdCpp(arma::vec x){
    arma::vec p = {0.25, 0.75};
    arma::vec r = Quantile(x, p);
    double h = (r[1] - r[0])/1.34;
    double w = pow(x.n_elem, -0.2);
    double s = sqrt(var(x));
    double v = 4 * 1.06 * std::min(s, h) * w;
    return (v);
}

// Find Interval Numbers or Indices in C++
// param x numeric vector (orignial).
// param breaks numeric vector (new).
// return the vector of length \code{length(x)} with values in \code{0:N}
// this is like \code{findInterval()} of R base, but result of this is can be used 
// to the C++
//[[Rcpp::export]]
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

arma::vec Kde2dWeightedCpp(
                   arma::rowvec w,
                   arma::mat ax,
                   arma::mat ay,
                   arma::vec h,
                   arma::uvec indx,
                   arma::uvec indy
                   ){
    int n = ax.n_rows;

    w = w/sum(w) * w.n_elem;

    ax = ax / h[0];
    ay = ay / h[1];
    arma::mat v = repelem(w, n, 1);
    arma::mat u = arma::normpdf(ax) % v;
    arma::mat day = arma::normpdf(ay) % v;
    arma::mat daym = day.t();

    arma::mat z = (u * daym)/(accu(v) * h[0] * h[1]);

    arma::mat sz = z.submat(indx, indy);
    arma::vec res = sz.diag();
    return (res);
}

struct RunWkde : public Worker{
  const arma::mat& w;
  const arma::mat& ax;
  const arma::mat& ay;
  const arma::vec& H;
  const arma::uvec& indx;
  const arma::uvec& indy;
  simple_progress& p;

  arma::mat& result;

  RunWkde(const arma::mat& w, const arma::mat& ax, const arma::mat& ay, 
          const arma::vec& H, const arma::uvec& indx, const arma::uvec& indy, 
          simple_progress& p, mat& result): 
      w(w), ax(ax), ay(ay), H(H), indx(indx), indy(indy), p(p), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        result.col(i) = Kde2dWeightedCpp(w.row(i), ax, ay, H, indx, indy);
        p.increment();
    }
  }
};


// Obtaion the difference between the grid points and original points
// param grid the grid points in one direction
// param x the original points in one direction
// return a matrix of the difference between the grid points and original points
//[[Rcpp::export]]
arma::mat outergrid(arma::vec grid, arma::vec x){
    arma::mat gxm = repelem(grid, 1, x.n_elem);
    arma::mat xm = repelem(x, 1, grid.n_elem);

    arma::mat ax = gxm - xm.t();

    return(ax);
}

// Compute the Kullback–Leibler Divergence by permutating a weight vector.
// param w the weight vector (the expression of gene or score of pathway).
// param bg the kernel density of background (the result of CalBgSpatialKld).
// param axm matrix the difference between the original point and grid points in x direction.
// param aym matrix the difference between the original point and grid points in y direction.
// param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
// (see bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
// param indx the index of original point by mapping to the grid points in x direction.
// param indy the index of original point by mapping to the grid points in y direction.
// param random_times the permutation numbers for each weight to test whether
// it is significantly, default is 200.
// return a vector of Kullback–Leibler Divergence with permutation.
////[[Rcpp::export]]
arma::vec CalRandSpatialKld(
    arma::rowvec w,
    arma::vec bg,
    arma::mat axm,
    arma::mat aym,
    arma::vec h,
    arma::uvec indx,
    arma::uvec indy,
    dqrng::xoshiro256plus rng,
    int random_times
    ){
    
    arma::vec bootkld(random_times);

    for (int j = 0; j < random_times; j++){
        arma::rowvec s = w;
        std::shuffle(std::begin(s), std::end(s), rng);
        arma::vec z = Kde2dWeightedCpp(s, axm, aym, h, indx, indy);
        z = z + 1e-300;
        bootkld[j] = log(sum(z % log(z / bg)));
    }
    return(bootkld);
}

arma::rowvec CalKldPvalue(arma::vec boot, double x){
    double bmean = arma::mean(boot);
    double bsd = arma::stddev(boot);
    double pval = R::pnorm(x, bmean, bsd, 0, 0);
    arma::rowvec res = {x, bmean, bsd, pval};

    return (res);
}

// Compute Background 2D Kernel Density
// param coords coordinate matrix.
// param axm matrix the difference between the original point and grid points in x direction.
// param aym matrix the difference between the original point and grid points in y direction.
// param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
// param indx the index of original point by mapping to the grid points in x direction.
// param indy the index of original point by mapping to the grid points in y direction.
// (see MASS::bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
// return a vector of 2D weighted kernel density value of background without spatial variability. 
// [[Rcpp::export]]
arma::vec CalBgSpatialKld(
        arma::mat coords,
        arma::mat axm,
        arma::mat aym,
        arma::vec h,
        arma::uvec indx,
        arma::uvec indy
    ){
    double prop = 1.0 / coords.n_rows;
    arma::rowvec bgw = rep(prop, coords.n_rows);
    arma::vec bgm = Kde2dWeightedCpp(bgw, axm, aym, h, indx, indy);
    arma::vec bgkld = bgm + 1e-300;
    return (bgkld);
}

// Compute the Kullback–Leibler Divergence using 2D Kernel Density Estimation 
// With Weighted and Statistical Test With Permutation for single weight vector.
// param d the weight vector (the expression of gene or score of pathway).
// param bgkld the kernel density of background (the result of CalBgSpatialKld).
// // param axm matrix the difference between the original point and grid points in x direction. 
// // param aym matrix the difference between the original point and grid points in y direction.
// // param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
// // (see bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
// // param indx the index of original point by mapping to the grid points in x direction.
// // param indy the index of original point by mapping to the grid points in y direction.
// // param random_times the permutation numbers for each weight to test whether 
// // it is significantly, default is 200.
// // return a vector of input features about the statistical test value with 2D weighted kernel density 
// //  and Kullback–Leibler Divergence.
// // [[Rcpp::export]]
// arma::rowvec CalSpatialKld(
//                         arma::rowvec d,
//                         arma::vec bgkld,
//                         arma::mat axm,
//                         arma::mat aym,
//                         arma::vec h,
//                         arma::uvec indx,
//                         arma::uvec indy,
//                         dqrng::xoshiro256plus rng,
//                         int random_times,
//                       ){
// 
//     arma::vec z = Kde2dWeightedCpp(d, axm, aym, h, indx, indy);
//     z = z + 1e-300;
//     double kld = log(sum(z % log(z / bgkld)));
// 
//     arma::vec bootkld(rmat.n_rows);
// 
//     bootkld = CalRandSpatialKld(d, bgkld, axm, aym, h, indx, indy, rng, random_times);
// 
//     arma::rowvec res = CalKldPvalue(bootkld, kld);
// 
//     return(res);
// }

struct SpatialKldCalWorker : public Worker{
  const arma::mat& w;
  const arma::vec& bgkld;
  const arma::mat& axm;
  const arma::mat& aym;
  const arma::vec& H;
  const arma::uvec& indx;
  const arma::uvec& indy;
  simple_progress& p;
  const uint64_t seed;
  const int random_times;
  mat& result;

  SpatialKldCalWorker(const arma::mat& w, const arma::vec& bgkld,
      const arma::mat& axm, const arma::mat& aym, const arma::vec& H, const arma::uvec& indx,
         const arma::uvec& indy, simple_progress& p, const uint64_t seed, const int random_times, mat& result)
  : w(w), bgkld(bgkld), axm(axm), aym(aym), H(H), indx(indx), indy(indy), p(p), seed(seed),
    random_times(random_times), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    dqrng::xoshiro256plus rng(seed);  
    for (uword i = begin; i < end; i++){
        arma::vec z = Kde2dWeightedCpp(w.row(i), axm, aym, H, indx, indy);
        double kld = log(sum(z % log(z / bgkld)));
        dqrng::xoshiro256plus lrng(rng);
        lrng.long_jump(i + 1);
        arma::vec bootkld = CalRandSpatialKld(w.row(i), bgkld, axm, aym, H, indx, indy, lrng, random_times);

        result.row(i) = CalKldPvalue(bootkld, kld);
        p.increment();
    }
  }
};


// Two-Dimensional Weighted Kernel Density Estimation And Mapping the Result To Original Dimension
// param x The 2-D coordinate matrix
// param w The weighted sparse matrix, the number columns the same than the number rows than x.
// param l The limits of the rectangle covered by the grid as c(xl, xu, yl, yu)
// param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
// (see bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
// param adjust numeric value to adjust to bandwidth, default is 1.
// param n number of grid points in the two directions, default is 400.
// return a matrix of 2D Weighted Kernel Density Estimation
// [[Rcpp::export]]
arma::sp_mat CalWkdeParallel(arma::mat& x, arma::sp_mat& w, arma::vec& l, Nullable<NumericVector> h,
        double adjust = 1.0, int n = 400) {

  arma::mat wv = conv_to<arma::mat>::from(w);

  arma::mat result(x.n_rows, w.n_rows);

  arma::vec gx = arma::linspace(l[0], l[1], n);
  arma::vec gy = arma::linspace(l[2], l[3], n);

  arma::vec H(x.n_cols);
  if (h.isNull()){
    for (uword j=0; j < x.n_cols;j ++){
       H[j] = BandwidthNrdCpp(x.col(j)) / 4 * adjust;
    }
  }else{
    H = as<arma::vec>(h);
  }

  //mapping to original coords
  arma::uvec indx = findIntervalCpp(x.col(0), gx) - 1;
  arma::uvec indy = findIntervalCpp(x.col(1), gy) - 1;

  arma::mat ax = outergrid(gx, x.col(0));
  arma::mat ay = outergrid(gy, x.col(1));

  uword num = wv.n_rows;
  simple_progress p(num);
  RunWkde runWkde(wv, ax, ay, H, indx, indy, p, result);
  parallelFor(0, num, runWkde);

  arma::sp_mat res = conv_to<arma::sp_mat>::from(result.t());

  return (res);
}


// Compute the Kullback–Leibler Divergence using 2D Kernel Density Estimation 
// With Weighted And Statistical Test With Permutation.
// param coords coordinate matrix.
// param d matrix (the expression of gene or score of pathway).
// param l The limits of the rectangle covered by the grid as c(xl, xu, yl, yu).
// param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
// (see bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
// param n the Number of grid points in each direction, default is 100.
// param random_times the permutation numbers for each weight to test whether
// it is significantly, default is 100.
// [[Rcpp::export]]
arma::mat CalSpatialKldCpp(arma::mat coords, arma::sp_mat d, arma::vec l, Nullable<NumericVector> h,
        int n = 100, int random_times = 100){

    arma::mat w = conv_to<arma::mat>::from(d);

    arma::mat gx = arma::linspace(l[0], l[1], n);
    arma::mat gy = arma::linspace(l[2], l[3], n);

    arma::vec H(coords.n_cols);
    if (h.isNull()){
      for (uword j=0; j < coords.n_cols;j ++){
         H[j] = BandwidthNrdCpp(coords.col(j)) / 4;
      }
    }else{
      H = as<arma::vec>(h);
    }

    uword num = w.n_rows;

    arma::mat result(num, 4);

    //mapping to original coords
    arma::uvec indx = findIntervalCpp(coords.col(0), gx) - 1;
    arma::uvec indy = findIntervalCpp(coords.col(1), gy) - 1;

    arma::mat axm = outergrid(gx, coords.col(0));
    arma::mat aym = outergrid(gy, coords.col(1));

    arma::vec bgkld = CalBgSpatialKld(coords, axm, aym, H, indx, indy);

    //arma::umat rmat = generate_random_permuation(w.n_cols, random_times);

    simple_progress p(num);

    Rcpp::IntegerVector seed(2, dqrng::R_random_int);
    uint64_t seed2 = dqrng::convert_seed<uint64_t>(seed);

    SpatialKldCalWorker spatialKldCalWorker(w, bgkld, axm, aym, H, indx, indy, p, seed2, random_times, result);
    parallelFor(0, num, spatialKldCalWorker);


    //#ifdef _OPENMP
    //#pragma omp parallel for schedule(static)
    //#endif
    //for (uword i = 0; i < num; i++){
    //    result.row(i) = CalSpatialKld(w.row(i), bgkld, axm, aym, H, indx, indy, rmat);
    //}

    return (result);
}


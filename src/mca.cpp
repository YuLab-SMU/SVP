//// This was from the CelliD package,
//// modified to support sparse matrix as input,
//// and speed up using Eigen library.
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace arma;
using namespace Eigen;

// [[Rcpp::export]]
List MCAStep1(const arma::sp_mat& X) {
    arma::mat AM = conv_to<arma::mat>::from(X);
    arma::vec rmin = arma::min(AM,1);
    arma::vec rmax = arma::max(AM,1);
    arma::vec range = (rmax - rmin);
    AM.each_col() -= rmin;
    AM.each_col() /= range;
    AM = join_cols(AM, 1 - AM);
    //AM.clear();
    long total = arma::accu(AM);
    arma::rowvec colsum = arma::sum(AM,0);
    arma::colvec rowsum = arma::sum(AM,1);
    AM.each_row() /= sqrt(colsum);
    AM.each_col() /= sqrt(rowsum);
    arma::colvec Dc = 1/(sqrt(rowsum/total));
    return List::create(Named("Z") = wrap(AM),
                        Named("Dc") = wrap(Dc));
}

//[[Rcpp::export]]
List MCAStep2(const Eigen::MatrixXd Z, const Eigen::MatrixXd V, const Eigen::VectorXd Dc){
    Eigen::MatrixXd FeaturesCoordinates = Z * V;
    int Zcol = Z.cols();
    FeaturesCoordinates = FeaturesCoordinates.array().colwise() * Dc.array();
    return List::create(Named("cellsCoordinates") = wrap(std::sqrt(Zcol) * V),
            Named("featuresCoordinates") = wrap(FeaturesCoordinates.topRows(FeaturesCoordinates.rows()/2)));
}


// //[[Rcpp::export]]
// List RunSVD(arma::mat x, int nv, int nu){
//     arma::mat U;
//     arma::vec s;
//     arma::mat V;
// 
//     svd(U, s, V, x);
//     return List::create(Named("d") = s.head(nv), 
//                         Named("u") = U.head_cols(nu), 
//                         Named("v") = V.head_cols(nv));
// 
// }

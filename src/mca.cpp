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
    arma::colvec rmin = arma::min(AM,1);
    arma::colvec rmax = arma::max(AM,1);
    arma::colvec range = (rmax -rmin);
    AM.each_col() -= rmin;
    AM.each_col() /= range;
    arma::mat FM = join_cols(AM, 1 - AM);
    AM.clear();
    long total = arma::accu(FM);
    arma::rowvec colsum = arma::sum(FM,0);
    arma::colvec rowsum = arma::sum(FM,1);
    FM.each_row() /= sqrt(colsum);
    FM.each_col() /= sqrt(rowsum);
    arma::colvec Dc = 1/(sqrt(rowsum/total));
    return List::create(Named("Z") = wrap(FM),
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

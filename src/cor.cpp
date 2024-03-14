// This script refer to the https://systematicinvestor.github.io/Correlation-Rcpp
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace RcppEigen;
using namespace arma;
using namespace std;

arma::vec CalWeight(const arma::vec& x){
    double medx = median(x);
    double madx = median(abs(x - medx));
    if (madx == 0.0){
        madx = median(abs(x - mean(x)));
    }
    arma::vec ui = (x - medx)/(9*madx);
    arma::uvec I = (1 - abs(ui))>0;
    arma::vec w = pow((1 - pow(ui, 2)), 2) % I;
    //w.replace(datum::nan, 0.0001);
    return (w);
}

// [[Rcpp::export]]
SEXP MatMultCpp(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

struct cal_bicor : public Worker {
    const arma::mat& mat;
    arma::mat& rmat;
    cal_bicor(const arma::mat& mat,
            arma::mat& rmat):
        mat(mat), rmat(rmat){ }

    void operator()(size_t begin, size_t end){
        for (size_t c1 = begin; c1 < end; c1++){
            arma::vec wX = CalWeight(mat.col(c1));
            double medX = median(mat.col(c1));
            for (size_t c2 = 0; c2 < c1; c2++){
                arma::vec wY = CalWeight(mat.col(c2));
                double medY = median(mat.col(c2));
                arma::vec xW = (mat.col(c1) - medX) % wX;
                arma::vec yW = (mat.col(c2) - medY) % wY;
                double sXX =  sqrt(accu(pow(xW, 2)));
                double sYY = sqrt(accu(pow(yW, 2)));
                rmat(c1, c2) = accu((xW/sXX) % (yW/sYY));
                rmat(c2, c1) = rmat(c1, c2);
            }
        }
    }
};


struct cal_bicor2 : public Worker {
    const arma::mat& mat1;
    const arma::mat& mat2;
    const size_t nc2;
    arma::mat& rmat;
    cal_bicor2(const arma::mat& mat1, const arma::mat& mat2, const size_t nc2,
            arma::mat& rmat):
        mat1(mat1), mat2(mat2), nc2(nc2), rmat(rmat){ }

    void operator()(size_t begin, size_t end){
        for (size_t c1 = begin; c1 < end; c1++){
            arma::vec wX = CalWeight(mat1.col(c1));
            double medX = median(mat1.col(c1));
            for (size_t c2 = 0; c2 < nc2; c2++){
                arma::vec wY = CalWeight(mat2.col(c2));
                double medY = median(mat2.col(c2));
                arma::vec xW = (mat1.col(c1) - medX) % wX;
                arma::vec yW = (mat2.col(c2) - medY) % wY;
                double sXX =  sqrt(accu(pow(xW, 2)));
                double sYY = sqrt(accu(pow(yW, 2)));
                rmat(c1, c2) = accu((xW/sXX) % (yW/sYY));
            }
        }
    }
};

// pre-compute sum and stdev
struct cor_step1 :  public Worker {
    const RcppParallel::RMatrix<double> mat;
    const int rstart, rend, nperiod;

    RcppParallel::RVector<double> rsum, rstdev;

    cor_step1(const Rcpp::NumericMatrix& mat, const int rstart, const int rend,
            Rcpp::NumericVector rsum, Rcpp::NumericVector rstdev)
        : mat(mat), rstart(rstart), rend(rend), nperiod(rend - rstart),
        rsum(rsum), rstdev(rstdev) {  }

    void operator() (size_t begin, size_t end) {
        for (size_t c = begin; c < end; ++c) {
            double sum, sum2;
            sum = sum2 = 0;

            for (int r = rstart; r < rend; ++r) {
                double d = mat(r,c);
                sum += d;
                sum2 += pow(d,2);
            }
            rsum[c] = sum;
            rstdev[c] = sqrt(nperiod * sum2 - pow(sum,2));
        }
    }
};

// compute correlation
struct cor_step2 : public Worker {
    const RcppParallel::RMatrix<double> mat;
    const int rstart, rend, nperiod;
    const RcppParallel::RVector<double> sum, stdev;

    RcppParallel::RMatrix<double> rmat;

    cor_step2(const Rcpp::NumericMatrix& mat, const int rstart, const int rend,
            const Rcpp::NumericVector& sum, const Rcpp::NumericVector& stdev,
            Rcpp::NumericMatrix rmat)
        : mat(mat), rstart(rstart), rend(rend), nperiod(rend - rstart),
        sum(sum), stdev(stdev), rmat(rmat) {  }

    void operator()(size_t begin, size_t end) {
        for (size_t c1 = begin; c1 < end; c1++) {
            for (size_t c2 = 0; c2 < c1; c2++) {
                double sXY = 0;
                for (int r = rstart; r < rend; r++) {
                    sXY += mat(r, c1) * mat(r, c2);
                }

                double sX = sum[c1];
                double sY = sum[c2];
                double SDx = stdev[c1];
                double SDy = stdev[c2];
                rmat(c1, c2) = (nperiod * sXY - sX * sY) / (SDx * SDy);
                rmat(c2, c1) = rmat(c1, c2);
            }
        }
    }
};

NumericMatrix CalCor(const NumericMatrix& mat, const int rstart,
        const int rend) {
    int nc = mat.ncol();
    NumericVector rsum(nc), rstdev(nc);

    cor_step1 cor_step1(mat, rstart, rend, rsum, rstdev);
    parallelFor(0, nc, cor_step1);

    NumericMatrix rmat(nc, nc);

    cor_step2 cor_step2(mat, rstart, rend, rsum, rstdev, rmat);
    parallelFor(0, nc, cor_step2);

    return (rmat);
}

//[[Rcpp::export]]
arma::mat corCpp(arma::sp_mat x, arma::sp_mat y){
    arma::mat z = arma::cor(arma::conv_to<mat>::from(x), arma::conv_to<mat>::from(y));
    return(z);
}


// [[Rcpp::export]]
NumericMatrix CalParallelCor(
    arma::sp_mat& x 
  ){
    arma::mat m = trans(arma::conv_to<arma::mat>::from(x));
    NumericMatrix y = wrap(m);
    NumericMatrix res = CalCor(y, 0, y.nrow());
    res.fill_diag(1.0);
    return (res);
}

arma::mat CalBiCor(const arma::mat& mat){
    int nc = mat.n_cols;
    arma::mat res(nc, nc);
    cal_bicor cal_bicor(mat, res);
    parallelFor(0, nc, cal_bicor);
    return(res);
}

arma::mat CalBiCorTwoMatrix(const arma::mat& mat1, const arma::mat& mat2){
    size_t nc1 = mat1.n_cols;
    size_t nc2 = mat2.n_cols;
    arma::mat res(nc1, nc2);
    cal_bicor2 cal_bicor2(mat1, mat2, nc2, res);
    parallelFor(0, nc1, cal_bicor2);
    return(res);
}

//[[Rcpp::export]]
arma::mat CalParallelBiCor(arma::sp_mat& x){
    arma::mat m = trans(arma::conv_to<arma::mat>::from(x));
    arma::mat res = CalBiCor(m);
    res.diag().fill(1.0);
    return(res);
}


//[[Rcpp::export]]
arma::mat CalParallelBiCorTwoMatrix(arma::sp_mat& x, arma::sp_mat& y){
    arma::mat m1 = trans(arma::conv_to<arma::mat>::from(x));
    arma::mat m2 = trans(arma::conv_to<arma::mat>::from(y));
    arma::mat res = CalBiCorTwoMatrix(m1, m2);
    return(res);
}


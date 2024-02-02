// This script refer to the https://systematicinvestor.github.io/Correlation-Rcpp
#include <RcppParallel.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;
using namespace std;

arma::vec CalWeight(const arma::vec& x){
    double medx = median(x);
    double madx = median(abs(x - medx));
    arma::vec ui = (x - medx)/(9*madx);
    arma::uvec I = (1 - abs(ui))>0;
    arma::vec w = pow((1 - pow(ui, 2)), 2) % I;
    return (w);
}

struct cal_bicor : public Worker {
    const arma::mat& mat;
    const int rstart, rend;
    arma::mat& rmat;
    cal_bicor(const arma::mat& mat, const int rstart, const int rend,
            arma::mat& rmat):
        mat(mat), rstart(rstart), rend(rend), rmat(rmat){ }

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

arma::mat CalBiCor(const arma::mat& mat, const int rstart, const int rend){
    int nc = mat.n_cols;
    arma::mat res(nc, nc);
    cal_bicor cal_bicor(mat, rstart, rend, res);
    parallelFor(0, nc, cal_bicor);
    return(res);
}

//[[Rcpp::export]]
arma::mat CalParallelBiCor(arma::sp_mat& x){
    arma::mat m = trans(arma::conv_to<arma::mat>::from(x));
    arma::mat res = CalBiCor(m, 0, m.n_rows);
    res.diag().fill(1.0);
    return(res);
}

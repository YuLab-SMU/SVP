#include <RcppArmadillo.h>
using namespace arma;

arma::umat generate_random_permuation(int ncol, int permutation){
    arma::umat rmat(permutation, ncol);

    for (unsigned int i = 0; i < permutation; i ++){
        arma::urowvec tmp = arma::linspace<arma::urowvec>(0, ncol - 1, ncol);
        rmat.row(i) = arma::shuffle(tmp);
    }

    return(rmat);
}


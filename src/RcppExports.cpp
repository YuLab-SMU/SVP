// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ExtractFeatureScoreCpp
List ExtractFeatureScoreCpp(arma::sp_mat& x, CharacterVector& rnm, CharacterVector& cnm, Rcpp::List& g);
RcppExport SEXP _SVP_ExtractFeatureScoreCpp(SEXP xSEXP, SEXP rnmSEXP, SEXP cnmSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type rnm(rnmSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type cnm(cnmSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(ExtractFeatureScoreCpp(x, rnm, cnm, g));
    return rcpp_result_gen;
END_RCPP
}
// findIntervalCpp
arma::uvec findIntervalCpp(arma::vec x, arma::vec breaks);
RcppExport SEXP _SVP_findIntervalCpp(SEXP xSEXP, SEXP breaksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type breaks(breaksSEXP);
    rcpp_result_gen = Rcpp::wrap(findIntervalCpp(x, breaks));
    return rcpp_result_gen;
END_RCPP
}
// Kde2dWeightedCpp
arma::vec Kde2dWeightedCpp(arma::rowvec w, arma::mat ax, arma::mat ay, arma::vec h, arma::uvec indx, arma::uvec indy);
RcppExport SEXP _SVP_Kde2dWeightedCpp(SEXP wSEXP, SEXP axSEXP, SEXP aySEXP, SEXP hSEXP, SEXP indxSEXP, SEXP indySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ax(axSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ay(aySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indx(indxSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indy(indySEXP);
    rcpp_result_gen = Rcpp::wrap(Kde2dWeightedCpp(w, ax, ay, h, indx, indy));
    return rcpp_result_gen;
END_RCPP
}
// outergrid
arma::mat outergrid(arma::vec grid, arma::vec x);
RcppExport SEXP _SVP_outergrid(SEXP gridSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(outergrid(grid, x));
    return rcpp_result_gen;
END_RCPP
}
// CalRandSpatialKld
arma::vec CalRandSpatialKld(arma::rowvec w, arma::vec bg, arma::mat axm, arma::mat aym, arma::vec h, arma::uvec indx, arma::uvec indy, int random_times, double seed);
RcppExport SEXP _SVP_CalRandSpatialKld(SEXP wSEXP, SEXP bgSEXP, SEXP axmSEXP, SEXP aymSEXP, SEXP hSEXP, SEXP indxSEXP, SEXP indySEXP, SEXP random_timesSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bg(bgSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type axm(axmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aym(aymSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indx(indxSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indy(indySEXP);
    Rcpp::traits::input_parameter< int >::type random_times(random_timesSEXP);
    Rcpp::traits::input_parameter< double >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(CalRandSpatialKld(w, bg, axm, aym, h, indx, indy, random_times, seed));
    return rcpp_result_gen;
END_RCPP
}
// CalBgSpatialKld
arma::vec CalBgSpatialKld(arma::mat coords, arma::mat axm, arma::mat aym, arma::vec h, arma::uvec indx, arma::uvec indy);
RcppExport SEXP _SVP_CalBgSpatialKld(SEXP coordsSEXP, SEXP axmSEXP, SEXP aymSEXP, SEXP hSEXP, SEXP indxSEXP, SEXP indySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type axm(axmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aym(aymSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indx(indxSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indy(indySEXP);
    rcpp_result_gen = Rcpp::wrap(CalBgSpatialKld(coords, axm, aym, h, indx, indy));
    return rcpp_result_gen;
END_RCPP
}
// CalSpatialKld
arma::rowvec CalSpatialKld(arma::rowvec d, arma::vec bgkld, arma::mat axm, arma::mat aym, arma::vec h, arma::uvec indx, arma::uvec indy, int random_times, double seed);
RcppExport SEXP _SVP_CalSpatialKld(SEXP dSEXP, SEXP bgkldSEXP, SEXP axmSEXP, SEXP aymSEXP, SEXP hSEXP, SEXP indxSEXP, SEXP indySEXP, SEXP random_timesSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bgkld(bgkldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type axm(axmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aym(aymSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indx(indxSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indy(indySEXP);
    Rcpp::traits::input_parameter< int >::type random_times(random_timesSEXP);
    Rcpp::traits::input_parameter< double >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(CalSpatialKld(d, bgkld, axm, aym, h, indx, indy, random_times, seed));
    return rcpp_result_gen;
END_RCPP
}
// MCAStep1
List MCAStep1(arma::sp_mat& X);
RcppExport SEXP _SVP_MCAStep1(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MCAStep1(X));
    return rcpp_result_gen;
END_RCPP
}
// MCAStep2
List MCAStep2(NumericMatrix Z, NumericMatrix V, NumericVector Dc);
RcppExport SEXP _SVP_MCAStep2(SEXP ZSEXP, SEXP VSEXP, SEXP DcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dc(DcSEXP);
    rcpp_result_gen = Rcpp::wrap(MCAStep2(Z, V, Dc));
    return rcpp_result_gen;
END_RCPP
}
// parallelCalRWR
NumericMatrix parallelCalRWR(arma::sp_mat x, arma::sp_mat v, double restart, double stop_delta, int stop_step);
RcppExport SEXP _SVP_parallelCalRWR(SEXP xSEXP, SEXP vSEXP, SEXP restartSEXP, SEXP stop_deltaSEXP, SEXP stop_stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type restart(restartSEXP);
    Rcpp::traits::input_parameter< double >::type stop_delta(stop_deltaSEXP);
    Rcpp::traits::input_parameter< int >::type stop_step(stop_stepSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelCalRWR(x, v, restart, stop_delta, stop_step));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SVP_ExtractFeatureScoreCpp", (DL_FUNC) &_SVP_ExtractFeatureScoreCpp, 4},
    {"_SVP_findIntervalCpp", (DL_FUNC) &_SVP_findIntervalCpp, 2},
    {"_SVP_Kde2dWeightedCpp", (DL_FUNC) &_SVP_Kde2dWeightedCpp, 6},
    {"_SVP_outergrid", (DL_FUNC) &_SVP_outergrid, 2},
    {"_SVP_CalRandSpatialKld", (DL_FUNC) &_SVP_CalRandSpatialKld, 9},
    {"_SVP_CalBgSpatialKld", (DL_FUNC) &_SVP_CalBgSpatialKld, 6},
    {"_SVP_CalSpatialKld", (DL_FUNC) &_SVP_CalSpatialKld, 9},
    {"_SVP_MCAStep1", (DL_FUNC) &_SVP_MCAStep1, 1},
    {"_SVP_MCAStep2", (DL_FUNC) &_SVP_MCAStep2, 3},
    {"_SVP_parallelCalRWR", (DL_FUNC) &_SVP_parallelCalRWR, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_SVP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

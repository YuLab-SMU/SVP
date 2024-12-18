// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cal_local_moran_bv
arma::vec cal_local_moran_bv(arma::vec x, arma::vec y, arma::sp_mat weight);
RcppExport SEXP _SVP_cal_local_moran_bv(SEXP xSEXP, SEXP ySEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_local_moran_bv(x, y, weight));
    return rcpp_result_gen;
END_RCPP
}
// CalGlobalLeeParallel
Rcpp::List CalGlobalLeeParallel(arma::sp_mat& x, arma::sp_mat& wm, arma::urowvec f1, arma::urowvec f2, int permutation, int alternative, bool cal_pvalue);
RcppExport SEXP _SVP_CalGlobalLeeParallel(SEXP xSEXP, SEXP wmSEXP, SEXP f1SEXP, SEXP f2SEXP, SEXP permutationSEXP, SEXP alternativeSEXP, SEXP cal_pvalueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type wm(wmSEXP);
    Rcpp::traits::input_parameter< arma::urowvec >::type f1(f1SEXP);
    Rcpp::traits::input_parameter< arma::urowvec >::type f2(f2SEXP);
    Rcpp::traits::input_parameter< int >::type permutation(permutationSEXP);
    Rcpp::traits::input_parameter< int >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< bool >::type cal_pvalue(cal_pvalueSEXP);
    rcpp_result_gen = Rcpp::wrap(CalGlobalLeeParallel(x, wm, f1, f2, permutation, alternative, cal_pvalue));
    return rcpp_result_gen;
END_RCPP
}
// RunLocalLee
arma::vec RunLocalLee(arma::vec& x, arma::vec& y, arma::sp_mat& wm, double n);
RcppExport SEXP _SVP_RunLocalLee(SEXP xSEXP, SEXP ySEXP, SEXP wmSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type wm(wmSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(RunLocalLee(x, y, wm, n));
    return rcpp_result_gen;
END_RCPP
}
// RunLocalMoranBvPerm
arma::mat RunLocalMoranBvPerm(arma::vec& x, arma::vec& y, arma::sp_mat& wm, int n, int permutation);
RcppExport SEXP _SVP_RunLocalMoranBvPerm(SEXP xSEXP, SEXP ySEXP, SEXP wmSEXP, SEXP nSEXP, SEXP permutationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type wm(wmSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type permutation(permutationSEXP);
    rcpp_result_gen = Rcpp::wrap(RunLocalMoranBvPerm(x, y, wm, n, permutation));
    return rcpp_result_gen;
END_RCPP
}
// MatMultCpp
Eigen::MatrixXd MatMultCpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _SVP_MatMultCpp(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(MatMultCpp(A, B));
    return rcpp_result_gen;
END_RCPP
}
// SpMatElemMultiMat
arma::sp_mat SpMatElemMultiMat(const arma::sp_mat x, const arma::mat y);
RcppExport SEXP _SVP_SpMatElemMultiMat(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(SpMatElemMultiMat(x, y));
    return rcpp_result_gen;
END_RCPP
}
// SpMatElemMultiSpMat
arma::sp_mat SpMatElemMultiSpMat(const arma::sp_mat x, const arma::sp_mat y);
RcppExport SEXP _SVP_SpMatElemMultiSpMat(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(SpMatElemMultiSpMat(x, y));
    return rcpp_result_gen;
END_RCPP
}
// MatElemMultiMat
arma::mat MatElemMultiMat(const arma::mat x, const arma::mat y);
RcppExport SEXP _SVP_MatElemMultiMat(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(MatElemMultiMat(x, y));
    return rcpp_result_gen;
END_RCPP
}
// corCpp
arma::mat corCpp(arma::sp_mat x, arma::sp_mat y);
RcppExport SEXP _SVP_corCpp(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(corCpp(x, y));
    return rcpp_result_gen;
END_RCPP
}
// CalParallelCor
NumericMatrix CalParallelCor(arma::sp_mat& x);
RcppExport SEXP _SVP_CalParallelCor(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(CalParallelCor(x));
    return rcpp_result_gen;
END_RCPP
}
// CalParallelBiCor
arma::mat CalParallelBiCor(arma::sp_mat& x);
RcppExport SEXP _SVP_CalParallelBiCor(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(CalParallelBiCor(x));
    return rcpp_result_gen;
END_RCPP
}
// CalParallelBiCorTwoMatrix
arma::mat CalParallelBiCorTwoMatrix(arma::sp_mat& x, arma::sp_mat& y);
RcppExport SEXP _SVP_CalParallelBiCorTwoMatrix(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(CalParallelBiCorTwoMatrix(x, y));
    return rcpp_result_gen;
END_RCPP
}
// ExtractFeatureScoreCpp
List ExtractFeatureScoreCpp(NumericMatrix& x, CharacterVector& rnm, CharacterVector& cnm, Rcpp::List& g);
RcppExport SEXP _SVP_ExtractFeatureScoreCpp(SEXP xSEXP, SEXP rnmSEXP, SEXP cnmSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type rnm(rnmSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type cnm(cnmSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(ExtractFeatureScoreCpp(x, rnm, cnm, g));
    return rcpp_result_gen;
END_RCPP
}
// CalGearyscParallel
arma::mat CalGearyscParallel(arma::sp_mat& x, arma::sp_mat& wm, int permutation, int lower_tail);
RcppExport SEXP _SVP_CalGearyscParallel(SEXP xSEXP, SEXP wmSEXP, SEXP permutationSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type wm(wmSEXP);
    Rcpp::traits::input_parameter< int >::type permutation(permutationSEXP);
    Rcpp::traits::input_parameter< int >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(CalGearyscParallel(x, wm, permutation, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// CalGetisOrdParallel
arma::mat CalGetisOrdParallel(arma::sp_mat& x, arma::sp_mat& wm, int lower_tail);
RcppExport SEXP _SVP_CalGetisOrdParallel(SEXP xSEXP, SEXP wmSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type wm(wmSEXP);
    Rcpp::traits::input_parameter< int >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(CalGetisOrdParallel(x, wm, lower_tail));
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
// CalWkdeParallel
arma::sp_mat CalWkdeParallel(arma::mat& x, arma::sp_mat& w, arma::vec& l, Nullable<NumericVector> h, double adjust, int n);
RcppExport SEXP _SVP_CalWkdeParallel(SEXP xSEXP, SEXP wSEXP, SEXP lSEXP, SEXP hSEXP, SEXP adjustSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type l(lSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type adjust(adjustSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(CalWkdeParallel(x, w, l, h, adjust, n));
    return rcpp_result_gen;
END_RCPP
}
// CalSpatialKldCpp
arma::mat CalSpatialKldCpp(arma::mat coords, arma::sp_mat d, arma::vec l, Nullable<NumericVector> h, int n, int random_times);
RcppExport SEXP _SVP_CalSpatialKldCpp(SEXP coordsSEXP, SEXP dSEXP, SEXP lSEXP, SEXP hSEXP, SEXP nSEXP, SEXP random_timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type l(lSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type random_times(random_timesSEXP);
    rcpp_result_gen = Rcpp::wrap(CalSpatialKldCpp(coords, d, l, h, n, random_times));
    return rcpp_result_gen;
END_RCPP
}
// CalLocalGParallel
Rcpp::List CalLocalGParallel(arma::sp_mat& x, arma::sp_mat& wm);
RcppExport SEXP _SVP_CalLocalGParallel(SEXP xSEXP, SEXP wmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type wm(wmSEXP);
    rcpp_result_gen = Rcpp::wrap(CalLocalGParallel(x, wm));
    return rcpp_result_gen;
END_RCPP
}
// CalLocalMoranParallel
Rcpp::List CalLocalMoranParallel(arma::sp_mat& x, arma::sp_mat& wm);
RcppExport SEXP _SVP_CalLocalMoranParallel(SEXP xSEXP, SEXP wmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type wm(wmSEXP);
    rcpp_result_gen = Rcpp::wrap(CalLocalMoranParallel(x, wm));
    return rcpp_result_gen;
END_RCPP
}
// MCAStep1
List MCAStep1(const arma::sp_mat& X);
RcppExport SEXP _SVP_MCAStep1(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MCAStep1(X));
    return rcpp_result_gen;
END_RCPP
}
// MCAStep2
List MCAStep2(const Eigen::MatrixXd Z, const Eigen::MatrixXd V, const Eigen::VectorXd Dc);
RcppExport SEXP _SVP_MCAStep2(SEXP ZSEXP, SEXP VSEXP, SEXP DcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type V(VSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Dc(DcSEXP);
    rcpp_result_gen = Rcpp::wrap(MCAStep2(Z, V, Dc));
    return rcpp_result_gen;
END_RCPP
}
// CalF1Parallel
arma::mat CalF1Parallel(arma::sp_mat x, arma::sp_mat y);
RcppExport SEXP _SVP_CalF1Parallel(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(CalF1Parallel(x, y));
    return rcpp_result_gen;
END_RCPP
}
// CalMoransiParallel
arma::mat CalMoransiParallel(arma::sp_mat& x, arma::sp_mat& wm, bool scaled, int permutation, int lower_tail);
RcppExport SEXP _SVP_CalMoransiParallel(SEXP xSEXP, SEXP wmSEXP, SEXP scaledSEXP, SEXP permutationSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type wm(wmSEXP);
    Rcpp::traits::input_parameter< bool >::type scaled(scaledSEXP);
    Rcpp::traits::input_parameter< int >::type permutation(permutationSEXP);
    Rcpp::traits::input_parameter< int >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(CalMoransiParallel(x, wm, scaled, permutation, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
// fastPDist
SEXP fastPDist(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _SVP_fastPDist(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(fastPDist(A, B));
    return rcpp_result_gen;
END_RCPP
}
// pairKnnCpp
List pairKnnCpp(arma::mat x, arma::mat y, arma::uword topn);
RcppExport SEXP _SVP_pairKnnCpp(SEXP xSEXP, SEXP ySEXP, SEXP topnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uword >::type topn(topnSEXP);
    rcpp_result_gen = Rcpp::wrap(pairKnnCpp(x, y, topn));
    return rcpp_result_gen;
END_RCPP
}
// colKnnCpp
List colKnnCpp(const arma::sp_mat& x, arma::uword k, bool weight);
RcppExport SEXP _SVP_colKnnCpp(SEXP xSEXP, SEXP kSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(colKnnCpp(x, k, weight));
    return rcpp_result_gen;
END_RCPP
}
// parallelCalRWR
NumericMatrix parallelCalRWR(arma::sp_mat& x, arma::sp_mat& v, double restart, double stop_delta, int stop_step);
RcppExport SEXP _SVP_parallelCalRWR(SEXP xSEXP, SEXP vSEXP, SEXP restartSEXP, SEXP stop_deltaSEXP, SEXP stop_stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type restart(restartSEXP);
    Rcpp::traits::input_parameter< double >::type stop_delta(stop_deltaSEXP);
    Rcpp::traits::input_parameter< int >::type stop_step(stop_stepSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelCalRWR(x, v, restart, stop_delta, stop_step));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SVP_cal_local_moran_bv", (DL_FUNC) &_SVP_cal_local_moran_bv, 3},
    {"_SVP_CalGlobalLeeParallel", (DL_FUNC) &_SVP_CalGlobalLeeParallel, 7},
    {"_SVP_RunLocalLee", (DL_FUNC) &_SVP_RunLocalLee, 4},
    {"_SVP_RunLocalMoranBvPerm", (DL_FUNC) &_SVP_RunLocalMoranBvPerm, 5},
    {"_SVP_MatMultCpp", (DL_FUNC) &_SVP_MatMultCpp, 2},
    {"_SVP_SpMatElemMultiMat", (DL_FUNC) &_SVP_SpMatElemMultiMat, 2},
    {"_SVP_SpMatElemMultiSpMat", (DL_FUNC) &_SVP_SpMatElemMultiSpMat, 2},
    {"_SVP_MatElemMultiMat", (DL_FUNC) &_SVP_MatElemMultiMat, 2},
    {"_SVP_corCpp", (DL_FUNC) &_SVP_corCpp, 2},
    {"_SVP_CalParallelCor", (DL_FUNC) &_SVP_CalParallelCor, 1},
    {"_SVP_CalParallelBiCor", (DL_FUNC) &_SVP_CalParallelBiCor, 1},
    {"_SVP_CalParallelBiCorTwoMatrix", (DL_FUNC) &_SVP_CalParallelBiCorTwoMatrix, 2},
    {"_SVP_ExtractFeatureScoreCpp", (DL_FUNC) &_SVP_ExtractFeatureScoreCpp, 4},
    {"_SVP_CalGearyscParallel", (DL_FUNC) &_SVP_CalGearyscParallel, 4},
    {"_SVP_CalGetisOrdParallel", (DL_FUNC) &_SVP_CalGetisOrdParallel, 3},
    {"_SVP_findIntervalCpp", (DL_FUNC) &_SVP_findIntervalCpp, 2},
    {"_SVP_outergrid", (DL_FUNC) &_SVP_outergrid, 2},
    {"_SVP_CalBgSpatialKld", (DL_FUNC) &_SVP_CalBgSpatialKld, 6},
    {"_SVP_CalWkdeParallel", (DL_FUNC) &_SVP_CalWkdeParallel, 6},
    {"_SVP_CalSpatialKldCpp", (DL_FUNC) &_SVP_CalSpatialKldCpp, 6},
    {"_SVP_CalLocalGParallel", (DL_FUNC) &_SVP_CalLocalGParallel, 2},
    {"_SVP_CalLocalMoranParallel", (DL_FUNC) &_SVP_CalLocalMoranParallel, 2},
    {"_SVP_MCAStep1", (DL_FUNC) &_SVP_MCAStep1, 1},
    {"_SVP_MCAStep2", (DL_FUNC) &_SVP_MCAStep2, 3},
    {"_SVP_CalF1Parallel", (DL_FUNC) &_SVP_CalF1Parallel, 2},
    {"_SVP_CalMoransiParallel", (DL_FUNC) &_SVP_CalMoransiParallel, 5},
    {"_SVP_fastPDist", (DL_FUNC) &_SVP_fastPDist, 2},
    {"_SVP_pairKnnCpp", (DL_FUNC) &_SVP_pairKnnCpp, 3},
    {"_SVP_colKnnCpp", (DL_FUNC) &_SVP_colKnnCpp, 3},
    {"_SVP_parallelCalRWR", (DL_FUNC) &_SVP_parallelCalRWR, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_SVP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// achr
Rcpp::List achr(const Rcpp::List& model, const Rcpp::List& state, const arma::mat& warmupPnts, const int nPnts, const int stepsPerPnt);
RcppExport SEXP _gembox_achr(SEXP modelSEXP, SEXP stateSEXP, SEXP warmupPntsSEXP, SEXP nPntsSEXP, SEXP stepsPerPntSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type state(stateSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type warmupPnts(warmupPntsSEXP);
    Rcpp::traits::input_parameter< const int >::type nPnts(nPntsSEXP);
    Rcpp::traits::input_parameter< const int >::type stepsPerPnt(stepsPerPntSEXP);
    rcpp_result_gen = Rcpp::wrap(achr(model, state, warmupPnts, nPnts, stepsPerPnt));
    return rcpp_result_gen;
END_RCPP
}
// mepc
Rcpp::List mepc(const Rcpp::List& model, const double beta, const double damp, const unsigned int maxIter, const double dlb, const double dub, const double epsil, const bool ff, const unsigned int fIdx, double fMeans, double fVars);
RcppExport SEXP _gembox_mepc(SEXP modelSEXP, SEXP betaSEXP, SEXP dampSEXP, SEXP maxIterSEXP, SEXP dlbSEXP, SEXP dubSEXP, SEXP epsilSEXP, SEXP ffSEXP, SEXP fIdxSEXP, SEXP fMeansSEXP, SEXP fVarsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type damp(dampSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const double >::type dlb(dlbSEXP);
    Rcpp::traits::input_parameter< const double >::type dub(dubSEXP);
    Rcpp::traits::input_parameter< const double >::type epsil(epsilSEXP);
    Rcpp::traits::input_parameter< const bool >::type ff(ffSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type fIdx(fIdxSEXP);
    Rcpp::traits::input_parameter< double >::type fMeans(fMeansSEXP);
    Rcpp::traits::input_parameter< double >::type fVars(fVarsSEXP);
    rcpp_result_gen = Rcpp::wrap(mepc(model, beta, damp, maxIter, dlb, dub, epsil, ff, fIdx, fMeans, fVars));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gembox_achr", (DL_FUNC) &_gembox_achr, 5},
    {"_gembox_mepc", (DL_FUNC) &_gembox_mepc, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_gembox(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
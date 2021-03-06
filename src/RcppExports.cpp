// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mgrwc
List mgrwc(const arma::colvec y, const arma::mat X, const arma::colvec weight, const int maxit, const double tol, double alpha, double mu, const bool mufixed);
RcppExport SEXP _mommix_mgrwc(SEXP ySEXP, SEXP XSEXP, SEXP weightSEXP, SEXP maxitSEXP, SEXP tolSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP mufixedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const bool >::type mufixed(mufixedSEXP);
    rcpp_result_gen = Rcpp::wrap(mgrwc(y, X, weight, maxit, tol, alpha, mu, mufixed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mommix_mgrwc", (DL_FUNC) &_mommix_mgrwc, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_mommix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

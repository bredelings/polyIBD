// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// runMCMC_cpp
Rcpp::List runMCMC_cpp(Rcpp::List args, Rcpp::List args_functions);
RcppExport SEXP _polyIBD_runMCMC_cpp(SEXP argsSEXP, SEXP args_functionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args_functions(args_functionsSEXP);
    rcpp_result_gen = Rcpp::wrap(runMCMC_cpp(args, args_functions));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests();

static const R_CallMethodDef CallEntries[] = {
    {"_polyIBD_runMCMC_cpp", (DL_FUNC) &_polyIBD_runMCMC_cpp, 2},
    {"run_testthat_tests",   (DL_FUNC) &run_testthat_tests,   0},
    {NULL, NULL, 0}
};

RcppExport void R_init_polyIBD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

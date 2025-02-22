// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// learn
NumericVector learn(int rep_size, List reps, NumericVector locs, double loc, NumericVector t_x, NumericVector syl_counter, int num_dems, double a, double innov, double p_att);
RcppExport SEXP _TransmissionBias_learn(SEXP rep_sizeSEXP, SEXP repsSEXP, SEXP locsSEXP, SEXP locSEXP, SEXP t_xSEXP, SEXP syl_counterSEXP, SEXP num_demsSEXP, SEXP aSEXP, SEXP innovSEXP, SEXP p_attSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type rep_size(rep_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< double >::type loc(locSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t_x(t_xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type syl_counter(syl_counterSEXP);
    Rcpp::traits::input_parameter< int >::type num_dems(num_demsSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type innov(innovSEXP);
    Rcpp::traits::input_parameter< double >::type p_att(p_attSEXP);
    rcpp_result_gen = Rcpp::wrap(learn(rep_size, reps, locs, loc, t_x, syl_counter, num_dems, a, innov, p_att));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TransmissionBias_learn", (DL_FUNC) &_TransmissionBias_learn, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_TransmissionBias(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ExpMat_C
arma::mat ExpMat_C(double t, arma::mat A);
RcppExport SEXP _PhenoTraits_ExpMat_C(SEXP tSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(ExpMat_C(t, A));
    return rcpp_result_gen;
END_RCPP
}
// vExpMat_C
arma::mat vExpMat_C(NumericVector t, arma::mat Q);
RcppExport SEXP _PhenoTraits_vExpMat_C(SEXP tSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(vExpMat_C(t, Q));
    return rcpp_result_gen;
END_RCPP
}
// Exp_2Q_C
arma::vec Exp_2Q_C(double t1, double t2, double Tau, arma::mat Q1, arma::mat Q2);
RcppExport SEXP _PhenoTraits_Exp_2Q_C(SEXP t1SEXP, SEXP t2SEXP, SEXP TauSEXP, SEXP Q1SEXP, SEXP Q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< double >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q1(Q1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q2(Q2SEXP);
    rcpp_result_gen = Rcpp::wrap(Exp_2Q_C(t1, t2, Tau, Q1, Q2));
    return rcpp_result_gen;
END_RCPP
}
// vExp_2Q_t1_C
arma::mat vExp_2Q_t1_C(NumericVector t1, double t2, double Tau, arma::mat Q1, arma::mat Q2);
RcppExport SEXP _PhenoTraits_vExp_2Q_t1_C(SEXP t1SEXP, SEXP t2SEXP, SEXP TauSEXP, SEXP Q1SEXP, SEXP Q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< double >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q1(Q1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q2(Q2SEXP);
    rcpp_result_gen = Rcpp::wrap(vExp_2Q_t1_C(t1, t2, Tau, Q1, Q2));
    return rcpp_result_gen;
END_RCPP
}
// Exp_2Q_scen3_C
arma::vec Exp_2Q_scen3_C(double t1, double Tau, arma::mat Q1, arma::mat Q2);
RcppExport SEXP _PhenoTraits_Exp_2Q_scen3_C(SEXP t1SEXP, SEXP TauSEXP, SEXP Q1SEXP, SEXP Q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q1(Q1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q2(Q2SEXP);
    rcpp_result_gen = Rcpp::wrap(Exp_2Q_scen3_C(t1, Tau, Q1, Q2));
    return rcpp_result_gen;
END_RCPP
}
// vExp_2Q_scen3_C
arma::mat vExp_2Q_scen3_C(NumericVector t1, double Tau, arma::mat Q1, arma::mat Q2);
RcppExport SEXP _PhenoTraits_vExp_2Q_scen3_C(SEXP t1SEXP, SEXP TauSEXP, SEXP Q1SEXP, SEXP Q2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< double >::type Tau(TauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q1(Q1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q2(Q2SEXP);
    rcpp_result_gen = Rcpp::wrap(vExp_2Q_scen3_C(t1, Tau, Q1, Q2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PhenoTraits_ExpMat_C", (DL_FUNC) &_PhenoTraits_ExpMat_C, 2},
    {"_PhenoTraits_vExpMat_C", (DL_FUNC) &_PhenoTraits_vExpMat_C, 2},
    {"_PhenoTraits_Exp_2Q_C", (DL_FUNC) &_PhenoTraits_Exp_2Q_C, 5},
    {"_PhenoTraits_vExp_2Q_t1_C", (DL_FUNC) &_PhenoTraits_vExp_2Q_t1_C, 5},
    {"_PhenoTraits_Exp_2Q_scen3_C", (DL_FUNC) &_PhenoTraits_Exp_2Q_scen3_C, 4},
    {"_PhenoTraits_vExp_2Q_scen3_C", (DL_FUNC) &_PhenoTraits_vExp_2Q_scen3_C, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_PhenoTraits(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

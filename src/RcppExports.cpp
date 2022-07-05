// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CHcpp
Rcpp::List CHcpp(const Rcpp::NumericMatrix& animals, const Rcpp::NumericMatrix& traps, const Rcpp::NumericMatrix& Tsk, const Rcpp::NumericMatrix& gsb, const Rcpp::NumericVector& NT, const int detectfn, const int detectorcode, const int nontargetcode, const int btype, const int Markov, const Rcpp::IntegerVector& binomN);
RcppExport SEXP _ipsecr_CHcpp(SEXP animalsSEXP, SEXP trapsSEXP, SEXP TskSEXP, SEXP gsbSEXP, SEXP NTSEXP, SEXP detectfnSEXP, SEXP detectorcodeSEXP, SEXP nontargetcodeSEXP, SEXP btypeSEXP, SEXP MarkovSEXP, SEXP binomNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type animals(animalsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type traps(trapsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type gsb(gsbSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type NT(NTSEXP);
    Rcpp::traits::input_parameter< const int >::type detectfn(detectfnSEXP);
    Rcpp::traits::input_parameter< const int >::type detectorcode(detectorcodeSEXP);
    Rcpp::traits::input_parameter< const int >::type nontargetcode(nontargetcodeSEXP);
    Rcpp::traits::input_parameter< const int >::type btype(btypeSEXP);
    Rcpp::traits::input_parameter< const int >::type Markov(MarkovSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type binomN(binomNSEXP);
    rcpp_result_gen = Rcpp::wrap(CHcpp(animals, traps, Tsk, gsb, NT, detectfn, detectorcode, nontargetcode, btype, Markov, binomN));
    return rcpp_result_gen;
END_RCPP
}
// rpsvcpp
Rcpp::NumericVector rpsvcpp(const Rcpp::IntegerMatrix& sk, const Rcpp::NumericMatrix& traps);
RcppExport SEXP _ipsecr_rpsvcpp(SEXP skSEXP, SEXP trapsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type sk(skSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type traps(trapsSEXP);
    rcpp_result_gen = Rcpp::wrap(rpsvcpp(sk, traps));
    return rcpp_result_gen;
END_RCPP
}
// popcpp
Rcpp::NumericMatrix popcpp(const Rcpp::NumericMatrix& mask, Rcpp::NumericVector& prob, double& maskspacing, int& N);
RcppExport SEXP _ipsecr_popcpp(SEXP maskSEXP, SEXP probSEXP, SEXP maskspacingSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< double& >::type maskspacing(maskspacingSEXP);
    Rcpp::traits::input_parameter< int& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(popcpp(mask, prob, maskspacing, N));
    return rcpp_result_gen;
END_RCPP
}
// popevencpp
Rcpp::NumericMatrix popevencpp(const Rcpp::NumericMatrix& bounds, int& N);
RcppExport SEXP _ipsecr_popevencpp(SEXP boundsSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type bounds(boundsSEXP);
    Rcpp::traits::input_parameter< int& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(popevencpp(bounds, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ipsecr_CHcpp", (DL_FUNC) &_ipsecr_CHcpp, 11},
    {"_ipsecr_rpsvcpp", (DL_FUNC) &_ipsecr_rpsvcpp, 2},
    {"_ipsecr_popcpp", (DL_FUNC) &_ipsecr_popcpp, 4},
    {"_ipsecr_popevencpp", (DL_FUNC) &_ipsecr_popevencpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ipsecr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

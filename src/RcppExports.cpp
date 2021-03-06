// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// jmcs_main
Rcpp::List jmcs_main(SEXP tol, SEXP k, SEXP n1, SEXP p1, SEXP p2, SEXP p1a, SEXP maxiter, SEXP point, SEXP xs, SEXP ws, SEXP yfile, SEXP cfile, SEXP mfile, SEXP Betasigmafile, SEXP Sigcovfile, SEXP gamma1, SEXP gamma2, SEXP trace);
RcppExport SEXP _FastJM_jmcs_main(SEXP tolSEXP, SEXP kSEXP, SEXP n1SEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP p1aSEXP, SEXP maxiterSEXP, SEXP pointSEXP, SEXP xsSEXP, SEXP wsSEXP, SEXP yfileSEXP, SEXP cfileSEXP, SEXP mfileSEXP, SEXP BetasigmafileSEXP, SEXP SigcovfileSEXP, SEXP gamma1SEXP, SEXP gamma2SEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< SEXP >::type k(kSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p1a(p1aSEXP);
    Rcpp::traits::input_parameter< SEXP >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< SEXP >::type point(pointSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ws(wsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yfile(yfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cfile(cfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mfile(mfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Betasigmafile(BetasigmafileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Sigcovfile(SigcovfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gamma1(gamma1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type gamma2(gamma2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(jmcs_main(tol, k, n1, p1, p2, p1a, maxiter, point, xs, ws, yfile, cfile, mfile, Betasigmafile, Sigcovfile, gamma1, gamma2, trace));
    return rcpp_result_gen;
END_RCPP
}
// jmcsf_main
Rcpp::List jmcsf_main(SEXP tol, SEXP k, SEXP n1, SEXP p1, SEXP p2, SEXP p1a, SEXP maxiter, SEXP point, SEXP xs, SEXP ws, SEXP yfile, SEXP cfile, SEXP mfile, SEXP Betasigmafile, SEXP Sigcovfile, SEXP gammafile, SEXP trace);
RcppExport SEXP _FastJM_jmcsf_main(SEXP tolSEXP, SEXP kSEXP, SEXP n1SEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP p1aSEXP, SEXP maxiterSEXP, SEXP pointSEXP, SEXP xsSEXP, SEXP wsSEXP, SEXP yfileSEXP, SEXP cfileSEXP, SEXP mfileSEXP, SEXP BetasigmafileSEXP, SEXP SigcovfileSEXP, SEXP gammafileSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< SEXP >::type k(kSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type p1a(p1aSEXP);
    Rcpp::traits::input_parameter< SEXP >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< SEXP >::type point(pointSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ws(wsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yfile(yfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cfile(cfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mfile(mfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Betasigmafile(BetasigmafileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Sigcovfile(SigcovfileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type gammafile(gammafileSEXP);
    Rcpp::traits::input_parameter< SEXP >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(jmcsf_main(tol, k, n1, p1, p2, p1a, maxiter, point, xs, ws, yfile, cfile, mfile, Betasigmafile, Sigcovfile, gammafile, trace));
    return rcpp_result_gen;
END_RCPP
}
// SimData
Rcpp::List SimData(SEXP k_val, SEXP p1_val, SEXP p1a_val, SEXP p2_val, SEXP g_val, SEXP truebeta, SEXP truegamma, SEXP truevee1, SEXP truevee2, SEXP randeffect, SEXP yfn, SEXP cfn, SEXP mfn);
RcppExport SEXP _FastJM_SimData(SEXP k_valSEXP, SEXP p1_valSEXP, SEXP p1a_valSEXP, SEXP p2_valSEXP, SEXP g_valSEXP, SEXP truebetaSEXP, SEXP truegammaSEXP, SEXP truevee1SEXP, SEXP truevee2SEXP, SEXP randeffectSEXP, SEXP yfnSEXP, SEXP cfnSEXP, SEXP mfnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_val(k_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p1_val(p1_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p1a_val(p1a_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p2_val(p2_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type g_val(g_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type truebeta(truebetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type truegamma(truegammaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type truevee1(truevee1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type truevee2(truevee2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type randeffect(randeffectSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yfn(yfnSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cfn(cfnSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mfn(mfnSEXP);
    rcpp_result_gen = Rcpp::wrap(SimData(k_val, p1_val, p1a_val, p2_val, g_val, truebeta, truegamma, truevee1, truevee2, randeffect, yfn, cfn, mfn));
    return rcpp_result_gen;
END_RCPP
}
// SimDataSF
Rcpp::List SimDataSF(SEXP k_val, SEXP p1_val, SEXP p1a_val, SEXP p2_val, SEXP g_val, SEXP truebeta, SEXP truegamma, SEXP truevee, SEXP randeffect, SEXP yfn, SEXP cfn, SEXP mfn);
RcppExport SEXP _FastJM_SimDataSF(SEXP k_valSEXP, SEXP p1_valSEXP, SEXP p1a_valSEXP, SEXP p2_valSEXP, SEXP g_valSEXP, SEXP truebetaSEXP, SEXP truegammaSEXP, SEXP trueveeSEXP, SEXP randeffectSEXP, SEXP yfnSEXP, SEXP cfnSEXP, SEXP mfnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type k_val(k_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p1_val(p1_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p1a_val(p1a_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type p2_val(p2_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type g_val(g_valSEXP);
    Rcpp::traits::input_parameter< SEXP >::type truebeta(truebetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type truegamma(truegammaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type truevee(trueveeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type randeffect(randeffectSEXP);
    Rcpp::traits::input_parameter< SEXP >::type yfn(yfnSEXP);
    Rcpp::traits::input_parameter< SEXP >::type cfn(cfnSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mfn(mfnSEXP);
    rcpp_result_gen = Rcpp::wrap(SimDataSF(k_val, p1_val, p1a_val, p2_val, g_val, truebeta, truegamma, truevee, randeffect, yfn, cfn, mfn));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FastJM_jmcs_main", (DL_FUNC) &_FastJM_jmcs_main, 18},
    {"_FastJM_jmcsf_main", (DL_FUNC) &_FastJM_jmcsf_main, 17},
    {"_FastJM_SimData", (DL_FUNC) &_FastJM_SimData, 13},
    {"_FastJM_SimDataSF", (DL_FUNC) &_FastJM_SimDataSF, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_FastJM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

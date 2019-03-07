// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// nodeIniter
IntegerVector nodeIniter(IntegerVector s0, double p);
RcppExport SEXP _SimInfInference_nodeIniter(SEXP s0SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(nodeIniter(s0, p));
    return rcpp_result_gen;
END_RCPP
}
// set_seed
void set_seed(unsigned int seed);
RcppExport SEXP _SimInfInference_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    set_seed(seed);
    return R_NilValue;
END_RCPP
}
// csample_num
NumericVector csample_num(NumericVector x, int size, bool replace, NumericVector prob);
RcppExport SEXP _SimInfInference_csample_num(SEXP xSEXP, SEXP sizeSEXP, SEXP replaceSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type replace(replaceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(csample_num(x, size, replace, prob));
    return rcpp_result_gen;
END_RCPP
}
// testShuffle
IntegerVector testShuffle(IntegerVector array);
RcppExport SEXP _SimInfInference_testShuffle(SEXP arraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type array(arraySEXP);
    rcpp_result_gen = Rcpp::wrap(testShuffle(array));
    return rcpp_result_gen;
END_RCPP
}
// sample_pools_Rcpp
IntegerVector sample_pools_Rcpp(int S, int I);
RcppExport SEXP _SimInfInference_sample_pools_Rcpp(SEXP SSEXP, SEXP ISEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    rcpp_result_gen = Rcpp::wrap(sample_pools_Rcpp(S, I));
    return rcpp_result_gen;
END_RCPP
}
// predict_env_sample_SISe_Rcpp
bool predict_env_sample_SISe_Rcpp(float whpp);
RcppExport SEXP _SimInfInference_predict_env_sample_SISe_Rcpp(SEXP whppSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type whpp(whppSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_env_sample_SISe_Rcpp(whpp));
    return rcpp_result_gen;
END_RCPP
}
// sample_herd_Rcpp
bool sample_herd_Rcpp(int s, int i);
RcppExport SEXP _SimInfInference_sample_herd_Rcpp(SEXP sSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_herd_Rcpp(s, i));
    return rcpp_result_gen;
END_RCPP
}
// sample_single_node_Rcpp
LogicalVector sample_single_node_Rcpp(IntegerVector svec, IntegerVector ivec);
RcppExport SEXP _SimInfInference_sample_single_node_Rcpp(SEXP svecSEXP, SEXP ivecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type svec(svecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ivec(ivecSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_single_node_Rcpp(svec, ivec));
    return rcpp_result_gen;
END_RCPP
}
// sample_nodes_Rcpp
LogicalVector sample_nodes_Rcpp(IntegerMatrix resMat, int seed);
RcppExport SEXP _SimInfInference_sample_nodes_Rcpp(SEXP resMatSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type resMat(resMatSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_nodes_Rcpp(resMat, seed));
    return rcpp_result_gen;
END_RCPP
}
// findLast
int findLast(IntegerMatrix obs, int start, int node);
RcppExport SEXP _SimInfInference_findLast(SEXP obsSEXP, SEXP startSEXP, SEXP nodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type node(nodeSEXP);
    rcpp_result_gen = Rcpp::wrap(findLast(obs, start, node));
    return rcpp_result_gen;
END_RCPP
}
// timeMatch
IntegerVector timeMatch(IntegerMatrix sim, IntegerMatrix obs, IntegerVector nodes);
RcppExport SEXP _SimInfInference_timeMatch(SEXP simSEXP, SEXP obsSEXP, SEXP nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type sim(simSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nodes(nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(timeMatch(sim, obs, nodes));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SimInfInference_nodeIniter", (DL_FUNC) &_SimInfInference_nodeIniter, 2},
    {"_SimInfInference_set_seed", (DL_FUNC) &_SimInfInference_set_seed, 1},
    {"_SimInfInference_csample_num", (DL_FUNC) &_SimInfInference_csample_num, 4},
    {"_SimInfInference_testShuffle", (DL_FUNC) &_SimInfInference_testShuffle, 1},
    {"_SimInfInference_sample_pools_Rcpp", (DL_FUNC) &_SimInfInference_sample_pools_Rcpp, 2},
    {"_SimInfInference_predict_env_sample_SISe_Rcpp", (DL_FUNC) &_SimInfInference_predict_env_sample_SISe_Rcpp, 1},
    {"_SimInfInference_sample_herd_Rcpp", (DL_FUNC) &_SimInfInference_sample_herd_Rcpp, 2},
    {"_SimInfInference_sample_single_node_Rcpp", (DL_FUNC) &_SimInfInference_sample_single_node_Rcpp, 2},
    {"_SimInfInference_sample_nodes_Rcpp", (DL_FUNC) &_SimInfInference_sample_nodes_Rcpp, 2},
    {"_SimInfInference_findLast", (DL_FUNC) &_SimInfInference_findLast, 3},
    {"_SimInfInference_timeMatch", (DL_FUNC) &_SimInfInference_timeMatch, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SimInfInference(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

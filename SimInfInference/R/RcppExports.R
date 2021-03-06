# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Init the nodes in u0.
#' p proportion should be infected.
#' @export
nodeIniter <- function(s0, p) {
    .Call(`_SimInfInference_nodeIniter`, s0, p)
}

#' a vector shuffle
NULL

#' set the seed
#' @export
set_seed <- function(seed) {
    invisible(.Call(`_SimInfInference_set_seed`, seed))
}

#' draw random sample
#' @export
csample_num <- function(x, size, replace, prob = as.numeric( c())) {
    .Call(`_SimInfInference_csample_num`, x, size, replace, prob)
}

#' a tester for the shuffle
#' @export
testShuffle <- function(array) {
    .Call(`_SimInfInference_testShuffle`, array)
}

#' Sample the pools
#' @export
sample_pools_Rcpp <- function(S, I) {
    .Call(`_SimInfInference_sample_pools_Rcpp`, S, I)
}

#' Predict the outcome of the env. sample
#' @export
predict_env_sample_SISe_Rcpp <- function(whpp) {
    .Call(`_SimInfInference_predict_env_sample_SISe_Rcpp`, whpp)
}

#' Sample the herd
#' @export
sample_herd_Rcpp <- function(s, i) {
    .Call(`_SimInfInference_sample_herd_Rcpp`, s, i)
}

#' Sample a node
#' @export
sample_single_node_Rcpp <- function(svec, ivec) {
    .Call(`_SimInfInference_sample_single_node_Rcpp`, svec, ivec)
}

#' Sample multiple nodes
#' @export
sample_nodes_Rcpp <- function(resMat, seed) {
    .Call(`_SimInfInference_sample_nodes_Rcpp`, resMat, seed)
}

#' Helper function for timeMatch.
#' Find the last entry in the vector of the same entry as
#' in the "starting" position.
#' @export
findLast <- function(obs, start, node) {
    .Call(`_SimInfInference_findLast`, obs, start, node)
}

#' match element that we should extract.
#' @export
timeMatch <- function(sim, obs, nodes) {
    .Call(`_SimInfInference_timeMatch`, sim, obs, nodes)
}


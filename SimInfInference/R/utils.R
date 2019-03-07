## Robin Eriksson 2018

##' idea from R-blogger
##' By R on roryverse, September 29, 2018
##' @export
get_seed <- function() {
    out <- sample.int(.Machine$integer.max, 1)
}

##' A wrapped version of the node sampler including a seed
##'
##' @export
sample_nodes_Rcpp_seed <- function(resMat) {
    sample_nodes_Rcpp(resMat, get_seed())
}

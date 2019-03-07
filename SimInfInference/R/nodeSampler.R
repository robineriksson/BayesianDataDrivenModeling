##library(SimInf)
##library(zoo) # for date handeling.
##library(pvclust)
##library(parallel)


##' Sample pools. Group individuals three and three in pools.
##' Randomly determine the status of each pool based on the number of
##' positive samples in each pool and the sensitivity.
##'
##' @param S_n Number of susceptible individuals
##' @param I_n Mumber of infected individuals
##' @return data.frame with the number of positive pools (pos) and the
##'     total number of pools (n).
##' @export
sample_pools <- function(S_n, I_n) {
    stopifnot(identical(length(S_n), length(I_n)))
    do.call("rbind", apply(cbind(S_n, I_n), 1, function(x) {
        if (sum(x) < 1)
            return(data.frame(pos = 0, n = 0))

        ## Create pools 3 and 3 from the vector of individuals. Create
        ## a matrix nrow x 3 where each row is one pool. Ignore that
        ## the number of individuals are not always divisble by three.
        samples <- sample(c(rep(0, x[1]), rep(1, x[2])), sum(x))
        suppressWarnings(m <- matrix(samples, ncol = 3))

        ## Calculate the number of positive samples in each pool.
        pools <- rowSums(m)

        ## Randomly determine the status of each pool based on the
        ## number of positive samples in each pool and the
        ## sensitivity.
        ##
        ## Mark E Arnold, Johanne Ellis-Iversen, Alasdair J C Cook,
        ## Robert H Davies, Ian M McLaren, Anthony C S Kay, and Geoff
        ## C Pritchard. Investigation into the effectiveness of pooled
        ## fecal samples for detection of verocytotoxin-producing
        ## escherichia coli o157 in cattle. Journal of Veterinary
        ## Diagnostic Investigation, 20(1):21-27, Jan 2008.
        ##
        ## a <- exp(-3.69 + pos * 8.20 / 3)
        ## sensitivity <- a / (1 + a)
        sensitivity <- c(0, 0.2775461, 0.8552848, 0.9891212)
        status <- runif(length(pools)) < sensitivity[pools + 1]

        data.frame(pos = sum(status), n = length(status))
    }))
}

##' Computes the pool prevalence
##'
##' @param S_n Number of susceptible individuals
##' @param I_n Mumber of infected individuals
##' @return Pool prevalence
##' @export
pool_prevalence <- function(S_n, I_n){
    if ((S_n + I_n) < 1)
      return(0)
    ## Create pools 3 and 3 from the vector of individuals. Create a
    ## matrix nrow x 3 where each row is one pool. Ignore that the
    ## number of individuals are not always divisble by three.
    suppressWarnings(m <- matrix(sample(c(rep(0, S_n), rep(1, I_n)),
                                        S_n + I_n),
                                 ncol = 3))
    ## Randomly determine the status of each pool based on the number
    ## of positive samples in each pool and the sensitivity.
    pools <- apply(m, 1, function(x) {
        sensitivity <- c(0, 0.28, 0.86, 0.99)
        return(runif(1) < sensitivity[sum(x) + 1])
    })
    100 * sum(pools) / length(pools)
}

##' Predict outcome from a overshoe and pooled pat sample in a
##' simulated herd
##'
##' There must be at least one pool positive in the herd (whpp > 0)
##' Parameters from: Widgren, S.; Eriksson, E.; Aspan, A.; Emanuelson,
##' U.; Alenius, S. & Lindberg, A. Environmental sampling for
##' evaluating verotoxigenic Escherichia coli O157: H7 status in dairy
##' cattle herds. Journal of Veterinary Diagnostic Investigation
##' @param wgpp Vector with the within-group pool prevalence
##' @param whpp Vector with the within-herd pool prevalence
##' @return logical vector
##' @export
predict_env_sample <- function(wgpp, whpp) {
    a <- exp(-2.14 + 0.06 * (wgpp - whpp) + 0.11 * whpp)
    (runif(length(wgpp)) < (a / (1 + a))) & (whpp > 0)
}

##' like the above, but trying to emulate only using 1 group per herd.
##' @param whpp withing herd prevalence
##' @export
predict_env_sample_single <- function(whpp){
    a <- exp(-2.14 + 0.11 * whpp)
    runif(1) <  a / (1 + a)
}

##' Given the within herd prevalence determine if the node is pos or neg.
##'
##' @param whpp withing herd prevalence
##' @export
predict_env_sample_SISe <- function(whpp){
    ## constants come from data in the above mentioned article,
    ## when applying a logistic regression (supervised learning)
    ## to te herd only data.
    a <- exp(-382.32431 + 91.49417 * whpp)
    prob <- a / (1+a)
    if(is.na(prob))
        prob <- 1
    runif(1) < prob

}



##' Sample a herd (node) to to determine if pos or false.
##' @param node The index of the herd to sample
##' @param week The index of the week to sample
##' @param S scalar with Susceptible
##' @param I scalar with Infected
##' @export
##' @return TRUE or FALSE
sample_herd <- function(S, I){
    ## Sample pools
    pools <- sample_pools(S, I)


    ## Determine within-group pool prevalence
    whpp <- 100 * pools$pos / ifelse(pools$n > 0, pools$n, 1)


    ## ## Determine within-herd pool prevalence
    ##whpp <- 100 * pools$pos / pools$n
    ##whpp <- wgpp ## only 1 group.

    ## Predict outcome from sampling a simulated herd. The status is
    ## positive if the overshoe and pooled pat sample is positive in
    ## either the calves or the young stock age categories.
    ##predict_env_sample(wgpp, whpp)
    ##predict_env_sample_single(whpp)
    predict_env_sample_SISe(whpp)

}

##' Simulate environmental sampling in the nodes (in parallel!)
##' @param result The result from simulating the model, list of dataframes for each nodes to sample.
##' @param cl0 the number of cores (threads) to use.
##' @return sampled results. Positive or negative at the time vectors.
##' @export
sample_nodes_par<- function(result, cl0=NULL){
    if(is.null(cl0)){
        cNum <- parallel::detectCores()
        cl <- parallel::makeCluster(getOption("cl.cores", cNum))
        parallel::clusterExport(cl=cl,
                      varlist=c("sample_herd","predict_env_sample_SISe","pool_prevalence","sample_pools"))
    }else
        cl <- cl0

    sampling <- parallel::parLapply(cl, X = result, function(X){
        ##sampling <- lapply(X = result, function(X){
        ## sample each node at each timestamp
        outcome <- apply(X,1,function(data){
            ## Determine the number of susceptible and infected
            S <- data["S"]
            I <- data["I"]

            sample_herd(S,I)
        })
        times <- X$time ##as.Date(res[,"date"])
        node <- X$node[1]
        data.frame(sample = as.numeric(outcome), ## {1, 0} instead of {true, false}
                   time = times,
                   node = rep(node, length(times)))
    })

    if(is.null(cl0))
        parallel::stopCluster(cl)


    post_sampling <- data.table::rbindlist(lapply(sampling, function(x){x}))

    return(post_sampling)
}

##' Simulate environmental sampling in the nodes (in parallel!)
##' @param result The result from simulating the model, list of dataframes for each nodes to sample.
##' @param cl0 the number of cores (threads) to use.
##' @return sampled results. Positive or negative at the time vectors.
##' @export
sample_nodes_par_Rcpp<- function(result, cl0=NULL){
    if(is.null(cl0)){
        cNum <- parallel::detectCores()
        cl <- parallel::makeCluster(getOption("cl.cores", cNum))
        parallel::clusterExport(cl=cl,
                      varlist=c("sample_herd","predict_env_sample_SISe","pool_prevalence","sample_pools"))
    }else
        cl <- cl0

    sampling <- parallel::parLapply(cl, X = result, function(X){
        ##sampling <- lapply(X = result, function(X){
        ## sample each node at each timestamp
        outcome <- apply(X,1,function(data){
            ## Determine the number of susceptible and infected
            S <- data["S"]
            I <- data["I"]

            sample_herd_Rcpp(S,I)
        })
        times <- X$time ##as.Date(res[,"date"])
        node <- X$node[1]
        data.frame(sample = as.numeric(outcome), ## {1, 0} instead of {true, false}
                   time = times,
                   node = rep(node, length(times)))
    })

    if(is.null(cl0))
        parallel::stopCluster(cl)


    post_sampling <- data.table::rbindlist(lapply(sampling, function(x){x}))

    return(post_sampling)
}

##' Simulate environmental sampling in the nodes (in parallel!)
##' @param result The result from simulating the model, list of dataframes for each nodes to sample.
##' @param cl0 the number of cores (threads) to use.
##' @return sampled results. Positive or negative at the time vectors.
##' @export
sample_nodes<- function(result){
    sampling <- lapply(result, function(X){
        ##sampling <- lapply(X = result, function(X){
        ## sample each node at each timestamp
        outcome <- apply(X,1,function(data){
            ## Determine the number of susceptible and infected
            S <- data["S"]
            I <- data["I"]

            sample_herd(S,I)
        })
        times <- X$time ##as.Date(res[,"date"])
        node <- X$node[1]
        data.frame(sample = as.numeric(outcome), ## {1, 0} instead of {true, false}
                   time = times,
                   node = rep(node, length(times)))
        ## list(sample = as.numeric(outcome), ## {1, 0} instead of {true, false}
        ##      time = times,
        ##      nodeId = node)
    })

    ## Combine the sampled results into a larger data.frame, like the input
    ## rbindlist exists, and might be faster.
    ## post_sampling <- do.call(rbind,lapply(sampling, function(x){
    ##     data.frame(sample = x$sample,
    ##                date = x$date,
    ##                node = rep(x$nodeId,length(x$sample)))
    ## }))

    ## post_sampling <- data.table::rbindlist(
    ##                                  lapply(sampling, function(x){
    ##                                      data.frame(sample = x$sample,
    ##                                                 time = x$time,
    ##                                                 node = rep(x$nodeId,length(x$sample)))
    ##                                  })
    ##                              )

    post_sampling <- data.table::rbindlist(lapply(sampling, function(x){x}))
    return(post_sampling)
}

##' Simulate environmental sampling in the nodes (in parallel!)
##' @param result The result from simulating the model, list of dataframes for each nodes to sample.
##' @param cl0 the number of cores (threads) to use.
##' @return sampled results. Positive or negative at the time vectors.
##' @export
sample_nodes_semi_Rcpp <- function(result){
    sampling <- lapply(result, function(X){
        ##sampling <- lapply(X = result, function(X){
        ## sample each node at each timestamp
        outcome <- apply(X,1,function(data){
            ## Determine the number of susceptible and infected
            S <- data["S"]
            I <- data["I"]

            sample_herd_Rcpp(S,I)
        })
        times <- X$time ##as.Date(res[,"date"])
        node <- X$node[1]
        data.frame(sample = as.numeric(outcome), ## {1, 0} instead of {true, false}
                   time = times,
                   node = rep(node, length(times)))
        ## list(sample = as.numeric(outcome), ## {1, 0} instead of {true, false}
        ##      time = times,
        ##      nodeId = node)
    })

    ## Combine the sampled results into a larger data.frame, like the input
    ## rbindlist exists, and might be faster.
    ## post_sampling <- do.call(rbind,lapply(sampling, function(x){
    ##     data.frame(sample = x$sample,
    ##                date = x$date,
    ##                node = rep(x$nodeId,length(x$sample)))
    ## }))

    ## post_sampling <- data.table::rbindlist(
    ##                                  lapply(sampling, function(x){
    ##                                      data.frame(sample = x$sample,
    ##                                                 time = x$time,
    ##                                                 node = rep(x$nodeId,length(x$sample)))
    ##                                  })
    ##                              )

    post_sampling <- data.table::rbindlist(lapply(sampling, function(x){x}))
    return(post_sampling)
}

##Rcpp::sourceCpp("~/Gits/BPD/SimInfInference/src/nodeSampler.cpp")
cppTest <- function(){
    hello(2)

    m <- matrix(rnorm(9),3,3)
    print(m)

    row_max(m)
}

RcppBench <- function(N = 1e3, S = 40, I = 10){
    t1 <- Sys.time()
    a <- replicate(n = N, sample_herd(S=S,I=I))
    t1 <- Sys.time() - t1
    print(t1)


    t2 <- Sys.time()
    b <- replicate(n = N, sample_herd_Rcpp(s=S,i=I))
    t2 <- Sys.time() - t2
    print(t2)

    p <- sapply(list(a,b), mean)
    return(p)
}

##' Sample one node given a vector of S and I values.
##' @export
sample_single_node <- function(Svec, Ivec) {
    df <- data.frame(S=Svec, I = Ivec);
    out <- apply(X=df, 1, function(X) {
        sample_herd_Rcpp(X[1], X[2])
        })
    return(out)
}

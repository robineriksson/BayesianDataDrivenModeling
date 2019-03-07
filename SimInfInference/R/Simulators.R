##' This file holds all Simulator functions.
##' Robinit Eriksson @2017
## library(Matrix)
## library(data.table)
## source("~/Gits/BPD/R/methods/nodeSampler.R")

##' SimInf included data (model)
##'
##' @param theta named vector of model parameters
##' @param extraArgs (runSeed, threads, nObs, tObs, tspan, includeTrue)
##'  runSeed the seed used in the actual simulation
##'  threads number of threads used, if null all available are used
##'  nObs the nodes that are observed durring measurements. (needs to be a subset of the avaiable nodes)
##'  tObs at what times that the nodes are observed. (needs to be a subset of tspan)
##'  tspan time (integer) vector used in simulator. When the state has to be recorded.
##'  includeTrue true if the underlying true state should be returned aswell.
##' @return result matrix in a list. Columns are nodes, and rows are timestamps.
## theta = c(upsilon_1 = 0.01, upsilon_2 = 0.01, upsilon_3 = 0.01, beta_t1 = 0.095, beta_t2 = 0.012, beta_t3 = 0.1, beta_t4 = 0.15)
##runSeed = NULL, threads = NULL, nObs = NULL, tObs = NULL, tspan = seq(1,4*365, by= 7), includeTrue = FALSE){
##' @export
SimInfSimulator_sandbox <- function(theta, extraArgs){
      ## Make sure the package is loaded
    ##stopifnot(require(SimInf, quietly = TRUE))
    stopifnot(c("runSeed", "threads", "nObs", "tObs", "tspan", "solver","prevLevel", "prevHerds", "events",
                "u0", "binary", "phiLevel", "nSim") %in% names(extraArgs))
    ## extract all the needed variables from the extraArgs argument
    runSeed <- extraArgs$runSeed
    threads <- extraArgs$threads
    nObs <- extraArgs$nObs
    tObs <- extraArgs$tObs
    tspan <- extraArgs$tspan
    includeTrue <- extraArgs$includeTrue
    solver <- extraArgs$solver
    prevLevel <- extraArgs$prevLevel
    prevHerds <- extraArgs$prevHerds
    events <- extraArgs$events
    u0 <- extraArgs$u0
    binary <- extraArgs$binary
    phiLevel <- extraArgs$phiLevel
    nSim <- extraArgs$nSim


    data("nodes", package = "SimInf")

    ## distance matrix
    d <- SimInf::distance_matrix(x = nodes$x, y = nodes$y, cutoff = 2500)

    ## make sure theta is correct format.
    theta <- setTheta(theta)

    ## create and initialize model
    model <- create_model_SISe(theta = theta, u0 = u0, distance = d,
                             tspan = tspan, events = events)


    ## extra initiation ...
    model <- init_model_widgren(model = model, theta = theta, tspan = tspan, p = prevLevel,
                                phi = phiLevel)
    ##model@U_sparse <- U_sparse

    ## run model nSim times and store in a list
    ## traj <- list()
    ## for (n in seq_len(nSim)) {
    ##     res <- SimInf::run(model = model, threads = threads, solver = solver)
    ##     ## prevalence(res, cases = I ~ S*I)
    ##     fulltraj <- SimInf::trajectory(model = res, compartments = I~I+S, node = as.numeric(nObs))
    ##     traj[[n]] <- fulltraj[fulltraj$time %in% tObs,]
    ## }
    ##traj <- traj[order(traj$node),]

    ## out <- lapply(traj, function(x){
    ##     if (binary){
    ##         ## create a list. 1 element for each node.
    ##         traj.list <- list()
    ##         for(i in seq_len(length(nObs)))
    ##             traj.list[[i]] = x[x$node == nObs[i],]
    ##         out <- sample_nodes_par(traj.list, cl0 = cl)
    ##         out <- as.data.frame(out)
    ##     } else
    ##         out <- x
    ## })

    ## run model nSim times and store in a list
    tcols <- which(tspan %in% tObs)
    nodeNames <- rep(as.numeric(nObs), each = length(tcols))

    traj <- list()

    for (n in seq_len(nSim)) {
        res <- SimInf::run(model = model, threads = threads, solver = solver)
        fulltraj <- SimInf::trajectory(model = res, compartments = c("S","I"), node = as.numeric(nObs), as.is = TRUE)
        reducedTraj <- fulltraj[,tcols]
        if(binary){
            ## do node sampler!
            sample <- sample_nodes_Rcpp_seed(as.matrix(reducedTraj))
            ##sample <- sample_nodes_Rcpp(as.matrix(reducedTraj), seed = n)
        } else {
            si <- seq(1,2*length(nObs),2)
            ii <- seq(2,2*length(nObs),2)
            sample <- data.frame(S = c(t(reducedTraj[si,])), I = c(t(reducedTraj[ii,])))
        }

        traj.df <- data.frame(node = nodeNames,
                              time = rep(tObs, length(nObs)),
                              sample)

        ## ## This needs to be done only ones!
        ## if(flag){
        ##     ##obsTimeElem <- sapply(as.list(seq_len(length(nObs))), function(X) {
        ##     ##obsTimeElem <- sapply(as.list(nObs), function(X) {
        ##     jOld <- 1
        ##     jNew <- 0
        ##     for(i in 1:length(nObs)){
        ##         sim <- traj.df[traj.df$node == nObs[i],]
        ##         obs <- obsDates[obsDates$node == nObs[i],]
        ##         found <- which(sim$time %in% obs$numTime)
        ##         jNew <- jOld + length(found)
        ##         ##if(i == 1) jOld <- 1;
        ##         obsTimeElem[jOld:(jNew-1)] <- length(tObs)*(i-1) + found
        ##         jOld <- jNew
        ##     }

        ##     flag <- FALSE

        ##traj.df$time <- zoo::as.Date(traj.df$time, origin = "2005-01-01")
        traj[[n]] <- traj.df
        ##traj[[n]] <- reducedTraj ##fulltraj[fulltraj$time %in% tObs,]
    }

    out <- traj



    return(out)
}


##' SimInf included data (model)
##'
##' @param theta named vector of model parameters
##' @param extraArgs (runSeed, threads, nObs, tObs, tspan, includeTrue)
##'  runSeed the seed used in the actual simulation
##'  threads number of threads used, if null all available are used
##'  nObs the nodes that are observed durring measurements. (needs to be a subset of the avaiable nodes)
##'  tObs at what times that the nodes are observed. (needs to be a subset of tspan)
##'  tspan time (integer) vector used in simulator. When the state has to be recorded.
##'  includeTrue true if the underlying true state should be returned aswell.
##' @return result matrix in a list. Columns are nodes, and rows are timestamps.
## theta = c(upsilon_1 = 0.01, upsilon_2 = 0.01, upsilon_3 = 0.01, beta_t1 = 0.095, beta_t2 = 0.012, beta_t3 = 0.1, beta_t4 = 0.15)
##runSeed = NULL, threads = NULL, nObs = NULL, tObs = NULL, tspan = seq(1,4*365, by= 7), includeTrue = FALSE){
##' @export
SimInfSimulator_real <- function(theta, extraArgs){
      ## Make sure the package is loaded
    ##stopifnot(require(SimInf, quietly = TRUE))
    stopifnot(c("runSeed", "threads", "nObs", "tObs", "tspan","solver","prevLevel", "prevHerds", "events",
                "u0", "binary", "model", "phiLevel", "nSim", "obsDates", "useSMHI") %in% names(extraArgs))
    ## extract all the needed variables from the extraArgs argument
    runSeed <- extraArgs$runSeed
    threads <- extraArgs$threads
    nObs <- extraArgs$nObs
    tObs <- extraArgs$tObs
    tspan <- extraArgs$tspan
    solver <- extraArgs$solver
    prevLevel <- extraArgs$prevLevel
    prevHerds <- extraArgs$prevHerds
    events <- extraArgs$events
    u0 <- extraArgs$u0
    binary <- extraArgs$binary
    model <- extraArgs$model
    phiLevel <- extraArgs$phiLevel
    nSim <- extraArgs$nSim
    obsDates <- extraArgs$obsDates
    useSMHI <- extraArgs$useSMHI




    ## if the prior knowledge have been applied to the upsilon parameters. Adjust theta.
    theta <- setTheta(theta)

    ## create and initialize model

    ## load "model"
    if(is.null(model)) {
        if(useSMHI)
            load("~/Gits/BPD/R/DATA/secret/SISe_SMHI.rda")
        else
            load("~/Gits/BPD/R/DATA/secret/SISe.rda")
    }

    ##model <- init_model_SISe_real(model = model, theta = theta, tspan = tspan, prevLevel = prevLevel, prevHerds = prevHerds)
    model <- init_model_widgren(model = model, theta = theta, tspan = tspan, p = prevLevel, phi = phiLevel)

    ## run model nSim times and store in a list
    tcols <- which(tspan %in% tObs)
    nodeNames <- rep(as.numeric(nObs), each = length(tcols))

    traj <- list()
    obsTimeElem <- numeric(dim(obsDates)[1])
    flag <- TRUE
    for (n in seq_len(nSim)) {
        res <- SimInf::run(model = model, threads = threads, solver = solver)
        fulltraj <- SimInf::trajectory(model = res, compartments = c("S","I"), node = as.numeric(nObs), as.is = TRUE)
        reducedTraj <- fulltraj[,tcols]
        if(binary){
            ## do node sampler!
            sample <- sample_nodes_Rcpp_seed(as.matrix(reducedTraj))
            ##sample <- sample_nodes_Rcpp(as.matrix(reducedTraj))##, seed = n)
        } else {
            si <- seq(1,2*length(nObs),2)
            ii <- seq(2,2*length(nObs),2)
            sample <- data.frame(S = c(t(reducedTraj[si,])), I = c(t(reducedTraj[ii,])))
        }

        traj.df <- data.frame(node = nodeNames,
                              time = rep(tObs, length(nObs)),
                              sample)

        ## This needs to be done only ones!
        if(flag){
            ##obsTimeElem <- sapply(as.list(seq_len(length(nObs))), function(X) {
            ##obsTimeElem <- sapply(as.list(nObs), function(X) {
            jOld <- 1
            jNew <- 0
            for(i in 1:length(nObs)){
                sim <- traj.df[traj.df$node == nObs[i],]
                obs <- obsDates[obsDates$node == nObs[i],]
                found <- which(sim$time %in% obs$numTime)
                jNew <- jOld + length(found)
                ##if(i == 1) jOld <- 1;
                obsTimeElem[jOld:(jNew-1)] <- length(tObs)*(i-1) + found
                jOld <- jNew
            }

            flag <- FALSE
        }
        ##traj.df$time <- zoo::as.Date(traj.df$time, origin = "2005-01-01")
        traj[[n]] <- traj.df[obsTimeElem,]
        ##traj[[n]] <- reducedTraj ##fulltraj[fulltraj$time %in% tObs,]
    }

    out <- traj

    ## out <- lapply(traj, function(x){
    ##     if (binary){
    ##         ## create a list. 1 element for each node.
    ##         traj.list <- list()
    ##         for(i in seq_len(length(nObs))){
    ##             sim <- x[x$node == nObs[i],]
    ##             obs <- obsDates[obsDates$node == nObs[i],]

    ##             traj.list[[i]] = sim[sim$time %in% obs$numTime,]
    ##         }
    ##         out <- sample_nodes_par(traj.list, cl0 = cl)
    ##         out <- as.data.frame(out)
    ##     } else
    ##         out <- x
    ## })

    ## out <- parallel::parLapply(cl, X = traj, function(X){
    ##     if (binary){
    ##         ## create a list. 1 element for each node.
    ##         traj.list <- list()
    ##         for(i in seq_len(length(nObs))){
    ##             sim <- X[X$node == nObs[i],]
    ##             obs <- obsDates[obsDates$node == nObs[i],]

    ##             traj.list[[i]] = sim[sim$time %in% obs$numTime,]
    ##         }
    ##         out <- sample_nodes(traj.list)
    ##         out <- as.data.frame(out)
    ##     } else
    ##         out <- X
    ## })


    return(out)
}


##' Re-format theta to model format
##'
##' @param theta named parameter vector
##' @return named parameter vector of correct model format.
##' @export
setTheta <- function(theta){
       ## if the prior knowledge have been applied to the upsilon parameters. Adjust theta.
    if(length(theta) == 6){
        theta.out <- c(upsilon = theta["upsilon"], beta_t1 = theta["beta_t1"],
                       beta_t2 = theta["beta_t2"] , beta_t3 = theta["beta_t1"] ,
                       beta_t4 = theta["beta_t3"], gamma = theta["gamma"],
                       prev = theta["prev"])
        names(theta.out) <- c("upsilon",
                              "beta_t1","beta_t2","beta_t3","beta_t4",
                              "gamma", "prev")
    } else if(length(theta) == 5){
        theta.out <- c(upsilon = theta["upsilon"], beta_t1 = theta["beta_t1"],
                       beta_t2 = theta["beta_t2"] , beta_t3 = theta["beta_t3"] ,
                       beta_t4 = theta["beta_t1"], gamma = theta["gamma"])
        names(theta.out) <- c("upsilon",
                              "beta_t1","beta_t2","beta_t3","beta_t4",
                              "gamma")
    }
    else if(length(theta) == 3){
        theta.out <- c(upsilon = theta["upsilon"], beta_t1 = theta["beta"],
                       beta_t2 = theta["beta"] , beta_t3 = theta["beta"] ,
                       beta_t4 = theta["beta"], gamma = theta["gamma"])
        names(theta.out) <- c("upsilon",
                              "beta_t1","beta_t2","beta_t3","beta_t4",
                              "gamma")
    } else if(length(theta) == 4){
        theta.out <- c(upsilon = theta["upsilon"], beta_t1 = theta["beta_t1"],
                       beta_t2 = theta["beta_t2"] , beta_t3 = theta["beta_t2"] ,
                       beta_t4 = theta["beta_t1"], gamma = theta["gamma"])
        names(theta.out) <- c("upsilon",
                              "beta_t1","beta_t2","beta_t3","beta_t4",
                              "gamma")
    } else
        theta.out <- theta
    return(theta.out)
}

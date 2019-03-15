library(SimInfInference)

##' Performes SLAM on the SISe dataset 1600 network
##'
##' @param nStop number of accepted parameter proposals
##' @param nSim number of simulations for each evaluation
##' @param debug if we want debug messages from the estimator
##' @param solver the stochastic solver to use in SimInf
##' @param obsspan how far between each observation
##' @param binary if the data is binary filtered
##' @param S SLAM parameter 1
##' @param e SLAM parameter 2
##' @param normalize if the SS should be normalized to observation
##' @param pertubate if to pertubate the initial guess or not.
##' @param thetaTrue The truth used for the generation of the obs.
##' @param threads the number of threads to use in SimInf
##' @param bs if bootstrapped is to be used
##' @param obsSeed the seed used for the observation
##' @param seed the seed used for the param. estimation
##' @param logParam observe the log-param space?
SLAMInference <- function(nStop = 100, nSim = 25,
                          debug = FALSE, solver = "ssm", obsspan = 60,
                          binary = FALSE,
                          Sw = NULL, S = 1e-3, e = NULL,
                          normalize = TRUE, pertubate = TRUE,
                          thetaTrue = c(upsilon = 0.0075,
                                        beta_t1 = 0.050,
                                        beta_t2 = 0.085,
                                        gamma = 0.1),
                          threads = NULL,
                          bs = TRUE, seed = 0,
                          logParam = FALSE){

    if(!is.null(Sw))
        S <- Sw
    if(is.null(e))
        e <- S/100
    ## S = 2.5e-4 for 3 beta model? e = 0.5
    ## S = 1e-3 for 2 beta model
    set.seed(seed) ## set up simulator
    tspan <- seq(1,4*365, 1)
    tObs <- seq(200,4*365, obsspan)
    data("nodes", package = "SimInf")
    nObs <- sample(1:length(rownames(nodes)),100)    ## The Simulator.
    Simulator <- SimInfSimulator_sandbox

    if(binary){
        ##for parallell sampling
        column <- "sample"
        logical <- TRUE
    } else {
        column <- "I"
        logical <- FALSE
    }

    cl <- NULL

    ## load data that is included the SimInf package
    events <- SimInf::events_SISe()
    u0 <- SimInf::u0_SISe()

    ##phiLevel = 0.5
    phiLevel <- "local"
    if("prev" %in% names(thetaTrue))
        prevLevel <- NULL
    else
        prevLevel <- 0.1

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs,
                               nObs = nObs, runSeed = NULL, threads = threads,
                               includeTrue = FALSE, solver = solver, prevLevel = prevLevel,
                               prevHerds = prevLevel, phiLevel = phiLevel, u0 = u0,
                               events = events, binary = binary, cl = cl, nSim = nSim)

    ## The Summary statistics
    SummaryStatistics <- aggWilkinson
    extraArgsSummaryStatistics <- list(column = column, fun = mean,
                                       qtr = FALSE, bs = bs, B = 200*nSim,
                                       useW = FALSE, logical = logical)

    ## The Proposal When only estimating upsilon.
    Proposal <- NULL
    extraArgsProposal <- NULL

    ## The Estimator
    Estimator <- SLAM


    ## True parameter
    ## upsilon <- 0.0075
    ## beta_t1 <- 0.05
    ## beta_t2 <- 0.085 ## whic is = 1.69*beta_t1
    ## gamma <- 0.1



    ##thetaTrue <- c(upsilon = upsilon, beta_t1= beta_t1, beta_t2 = beta_t2, gamma = gamma)


    ## pertubate starting guess
    if(pertubate)
        pert <- runif(length(thetaTrue), 0.95, 1.05)
    else
        pert <- rep(1, length(thetaTrue))

    if(logParam)
        theta0 <- log(pert * thetaTrue)
    else
        theta0 <- pert*thetaTrue

    ## if(is.null(S))
    ##     S = 2.4^2/length(thetaTrue) ## optimal from authors.


    extraArgsEstimator <- list(parameters = names(thetaTrue),
                               nStop = nStop, nSim = nSim, debug = debug, theta0 = theta0,
                               normalize = normalize, e = e, C0 = 1e-2, S = S, n0 = 1, accVec = NULL,
                               reinit = FALSE, xbar = NULL, sigma = NULL, logParam = logParam)


    ## Create object!
    infe <- Inference$new(Simulator = Simulator, SummaryStatistics = SummaryStatistics,
                          Proposal = Proposal, Estimator = Estimator,
                          extraArgsSimulator = extraArgsSimulator,
                          extraArgsSummaryStatistics = extraArgsSummaryStatistics,
                          extraArgsEstimator = extraArgsEstimator,
                          extraArgsProposal = extraArgsProposal)

    ## make observation
    infe$changeExtraArgsSimulator(nSim = 1)
    obs <- infe$runSimulator(thetaTrue)
    infe$setObservation(obs)
    infe$changeExtraArgsSimulator(nSim = nSim)

    ## explore the posterior!
    infe$runEstimation()



    return(infe)
}

##' Continue the already created inference class
##'
##' @param infe the inference class
##' @param nStop number of added parameters
##' @param theta a proposed theta0
##' @return inference class
contInference <- function(infe, nStop = 100){
    accVec <- as.matrix(infe$getPosterior())
    d <- dim(accVec)
    if(infe$getExtraArgsEstimator()$logParam)
        accVec[,1:(d[2]-1)] <- log(accVec[,1:(d[2]-1)])

    infe$changeExtraArgsEstimator(nStop = nStop, accVec = accVec, reinit = TRUE)



    infe$runEstimation()

    return(infe)
}

##' pertubationSL - explore the 1d parameter SL space.
##' We want to see if the minimas are well defined.
##'
##' @param thetalength the length of the theta vector
##' @param nSim the number of simulations per SL eval.
##' @param pert how much do we perturbe from the truth?
##' @param solver which stochastic solver
##' @param binary apply filter or not
##' @param obsspan how far between each observation
##' @param bs bootstrap or not
##' @param multiSS how many times to evaluate the point.
pertubationSL <- function(thetalength = 21,
                          nSim = 20, pert = c(0.8,1.2),
                          solver = "ssm", binary = FALSE,
                          obsspan = 60,
                          bs = TRUE,
                          multiSS = 20){

    set.seed(0)
    ## set up simulator
    tspan <- seq(1,4*365, 1)
    tObs <- seq(100,4*365, obsspan)
    ##tObs <- seq(1,4*365, 42)
    data("nodes", package = "SimInf")

    nObs <- sample(1:length(rownames(nodes)),100)
    ##nObs <- 1:1600


    ## The Simulator.
    Simulator <- SimInfSimulator_sandbox

    if(binary){
        ##for parallell sampling
        column <- "sample"
        logical <- TRUE
    } else {
        column <- "I"
        logical <- FALSE
    }


    ## load data that is included the SimInf package
    events <- SimInf::events_SISe()
    u0 <- SimInf::u0_SISe()

    ##phiLevel = 0.5
    phiLevel <- "local"
    prevLevel <- 0.1

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs,
                               nObs = nObs, runSeed = NULL, threads = NULL,
                               includeTrue = FALSE, solver = solver, prevLevel = prevLevel,
                               prevHerds = prevLevel, phiLevel = phiLevel, u0 = u0,
                               events = events, binary = binary, cl = NULL, nSim = nSim)

    ## The Summary statistics
    ##SummaryStatistics <- SimInfStatistics_sandbox_full
    SummaryStatistics <- aggWilkinson
    extraArgsSummaryStatistics <- list(column = column, fun = mean,
                                       qtr = FALSE, bs = bs, B = 500,
                                       useW = FALSE, logical = logical)

    ## SummaryStatistics <- aggregateTS
    ## extraArgsSummaryStatistics <- list(column = column)
    ##extraArgsSummaryStatistics <- list(nodes = nodes, numClusters = numClusters, SI = SI, span = SSspan,
    ##betalag = betalag, betapower = betapower)

    ## The Proposal
    ## When only estimating upsilon
    Proposal <- NULL
    extraArgsProposal <- NULL


    ## The Estimator
    Estimator <- slGrid

    ## True parameter
    thetaTrue = c(upsilon = 0.0075,
                  beta_t1 = 0.050,
                  beta_t2 = 0.085,
                  gamma = 0.1)


    extraArgsEstimator <- list(theta0 = thetaTrue,
                               pert = pert,
                               thetalength = thetalength,
                               allDim = FALSE,
                               logParam = FALSE,
                               nSim = nSim,
                               normalize = TRUE,
                               multiSS = multiSS)

    ## Create object!
    infe <- Inference$new(Simulator = Simulator, SummaryStatistics = SummaryStatistics,
                          Proposal = Proposal, Estimator = Estimator,
                          extraArgsSimulator = extraArgsSimulator,
                          extraArgsSummaryStatistics = extraArgsSummaryStatistics,
                          extraArgsEstimator = extraArgsEstimator,
                          extraArgsProposal = extraArgsProposal)


    ## make observation
    infe$changeExtraArgsSimulator(nSim = 1)
    obs <- infe$runSimulator(thetaTrue)
    infe$setObservation(obs)
    infe$changeExtraArgsSimulator(nSim = nSim)

    infe$setObservation(obs)

    infe$runEstimation()

    return(infe)
}

restartBugTest <- function(N=1000, n = 25, divisions = 4, S = 1e-4, solver = "ssm", debug = FALSE){

    sAll <- SLAMInference(N, n, S = S, solver = solver, debug = debug)
    dAll <- sAll$getPosterior()

    sShort <- SLAMInference(N/divisions, n, S = S, solver = solver, debug = debug)
    for(i in seq_len(divisions-1)){
        sShort <- contInference(sShort,N/divisions)
    }
    dShort <- sShort$getPosterior()


    d <- data.frame(all = dAll[,5], divided = dShort[,5])
    rate <- apply(d,2,function(x){length(unique(x))})/ N
    print(rate)
    return(list(d = d, sAll = sAll, sShort = sShort))
}

foraNight <- function(N, S){
    sfull <- SLAMInference(nStop = N, S = S, pertubate = TRUE)
    sbin <- SLAMInference(nStop = N, S = S, column = "sample", binary = TRUE, pertubate = TRUE)

    return(list(sfull = sfull, sbin = sbin))
}

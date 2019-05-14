library(SimInfInference)

##' Performes SLAM on the SISe dataset with real observations
##'
##' @param nStop number of accepted parameter proposals
##' @param nSim number of simulations for each evaluation
##' @param debug if we want debug messages from the estimator
##' @param solver the stochastic solver to use in SimInf
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
##' @param useW use weighted Summary Statistics
##' @param B the number of bootstrap replicates
##' @param dataDir the directory holding the epi-modeldata
SLAMInference <- function(nStop = 1e3, nSim = 10,
                          debug = FALSE, solver = "ssm",
                          Sw = NULL,
                          S = 1e-4, e = NULL,
                          normalize = TRUE, pertubate = TRUE,
                          thetaTrue = c(upsilon = 0.01563841,
                                        beta_t1 = 0.14309175,
                                        beta_t2 = 0.13366590,
                                        beta_t3 = 0.15541914,
                                        gamma = 0.09720055,
                                        prev = 0.02537256),
                          threads = NULL,
                          bs = TRUE, seed = 0,
                          logParam = FALSE, useW = TRUE, B = 100,
                          dataDir){
    if(!is.null(Sw))
        S <- Sw
    if(is.null(e))
        e <- S/100
    ## S = 1e-3 -> 20% accrate for inc. prev.
    ## e = 0.5
    ## S = 5e-5 -> 3.56%
    ## S = 5e-4 -> 1.4% or 2.2% acc.rate
    set.seed(seed) ## set up simulator


    paths <- NULL
    paths["obs"] <- paste(dataDir, "obs.rda", sep="")
    paths["model"] <- paste(dataDir, "SISe.rda", sep="")
    paths["nObs"] <- paste(dataDir, "nObs.RData", sep="")


    load(paths["obs"]) ## loads: obs
    obs$sample <- as.numeric(obs$status)


    ## We only to run the simulation for the maximum time that we have observations for.
    ## But we should start at the same time as the movements as they will affect
    ## the state of the distribution of induvidials at nodes.
    load(paths["model"]) ## loads: model


    tspan0 <- seq(head(model@events@time,1), tail(model@events@time,1), 1)
    realDates <- zoo::as.Date(tspan0, origin = "2005-01-01")#"2007-07-01")
    tinObs <- sort(unique(obs$time))
    tObs <- as.numeric(tinObs) - (as.numeric(tinObs[1]) - tspan0[which(realDates == tinObs[1])])
    obs$numTime <- as.numeric(obs$time) - as.numeric(tinObs[1]) +
        tspan0[which(realDates == tinObs[1])]

    tspan <- tspan0[-which(realDates > tail(tinObs,1))]


    ## load names of the observed nodes -> nObs
    load(paths["nObs"]) ## load: nObs


    ## The Simulator.
    Simulator <- SimInfSimulator_real

    events <- NULL
    u0 <- NULL

    ## loads variable: model
    ##phiLevel <- 0.05
    phiLevel <- "local"
    if("prev" %in% names(thetaTrue))
        prevLevel = NULL
    else
        prevLevel = 0.02

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, solver = solver,
                               prevLevel = prevLevel, prevHerds = prevLevel, phiLevel = phiLevel,
                               u0 = u0, events = events, binary = TRUE, model = model,
                               nSim = nSim, obsDates = obs, dataDir = dataDir)


    ## The Summary statistics
    SummaryStatistics <- aggWilkinson
    if(is.null(B))
        B <- nSim*20
    fun = sum
    extraArgsSummaryStatistics <- list(column = "sample", fun = fun, qtr = TRUE,  bs = bs,
                                       B = B, useW = useW, logical = FALSE)

    Proposal <- NULL
    extraArgsProposal <- NULL

    ## The Estimator
    Estimator <- SLAM

    ## pertubate starting guess
    if(pertubate)
        pert <- runif(length(thetaTrue), 0.95, 1.05)
    else
        pert <- rep(1, length(thetaTrue))

    if(logParam)
        theta0 <- log(pert * thetaTrue)
    else
        theta0 <- pert*thetaTrue

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


    ## load in observation
    infe$setObservation(list(obs))

    ## explore the posterior!
    infe$runEstimation()

    ##parallel::stopCluster(cl)


    return(infe)
}

##' Continue the already created inference class
##'
##' @param infe the inference class
##' @param nStop number of added parameters
##' @param theta a proposed theta0
##' @return inference class
contInference <- function(infe, nStop = 100){

    ## infe <- reinitCluster(infe)

    accVec <- as.matrix(infe$getPosterior())
    d <- dim(accVec)
    if(infe$getExtraArgsEstimator()$logParam)
        accVec[,1:(d[2]-1)] <- log(accVec[,1:(d[2]-1)])

    infe$changeExtraArgsEstimator(nStop = nStop, accVec = accVec, reinit = TRUE)

    infe$runEstimation()

    return(infe)
}

##' 1D marginal pertubation exploration
##' @param nSim the number of simulations used for each SL
##' @param pert how big the pertubation should be
##' @param thetalength how many points should be explored
##' @param solver what stochastic numerical solver to use
##' @param bs if bootstrap should be used
##' @param multiSS how many times SL should be evaluated at each point
##' @param thetaTrue the centroid to start from
##' @param useW if weighted samples should be used
##' @param allDim if we should do 1D or 2D marginal exploration
##' @param B how many bootstrap samples to use.
##' @param dataDir the directory holding the observation/model data
pertubationSL <- function(nSim = 20, pert = c(0.8,1.2),
                          thetalength = 11,
                          solver = "ssm",
                          bs = TRUE,
                          multiSS = 1,
                          thetaTrue = c(upsilon = 0.01563841,
                                        beta_t1 = 0.14309175,
                                        beta_t2 = 0.13366590,
                                        beta_t3 = 0.15541914,
                                        gamma = 0.09720055,
                                        prev = 0.02537256),
                          useW = TRUE,
                          allDim = FALSE,
                          B = NULL,
                          dataDir
                          ){

    paths <- NULL
    paths["obs"] <- paste(dataDir, "obs.rda", sep="")
    paths["model"] <- paste(dataDir, "SISe.rda", sep="")
    paths["nObs"] <- paste(dataDir, "nObs.RData", sep="")

    set.seed(0)
    ## set up simulator
    load(paths["obs"]) ## loads: obs
    obs$sample <- as.numeric(obs$status)

    load(paths["model"]) ## loads: model

    tspan0 <- seq(head(model@events@time,1), tail(model@events@time,1), 1)
    realDates <- zoo::as.Date(tspan0, origin = "2005-01-01")#"2007-07-01")
    tinObs <- sort(unique(obs$time))
    tObs <- as.numeric(tinObs) - (as.numeric(tinObs[1]) - tspan0[which(realDates == tinObs[1])])
    obs$numTime <- as.numeric(obs$time) - as.numeric(tinObs[1]) +
        tspan0[which(realDates == tinObs[1])]

    tspan <- tspan0[-which(realDates > tail(tinObs,1))]

    ## load names of the observed nodes -> nObs
    load(paths["nObs"]) ## load: nObs


    ## The Simulator.
    Simulator <- SimInfSimulator_real

    ## load data that is included the SimInf package
    events <- NULL
    u0 <- NULL

    phiLevel <- "local"
    if("prev" %in% names(thetaTrue))
        prevLevel = NULL
    else
        prevLevel = 0.02


    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, solver = solver,
                               prevLevel = prevLevel, prevHerds = prevLevel, phiLevel = phiLevel,
                               u0 = u0, events = events, binary = TRUE, model = model,
                               nSim = nSim, obsDates = obs)

    ## The Summary statistics
    SummaryStatistics <- aggWilkinson


    if(is.null(B))
        B <- nSim*20

    fun <- sum
    extraArgsSummaryStatistics <- list(column = "sample", fun = fun, qtr = TRUE,  bs = bs,
                                       B = B, useW = useW, logical = FALSE)

    ## The Proposal
    ## When only estimating upsilon
    Proposal <- NULL
    extraArgsProposal <- NULL


    ## The Estimator
    Estimator <- slGrid




    extraArgsEstimator <- list(theta0 = thetaTrue,
                               pert = pert,
                               thetalength = thetalength,
                               allDim = allDim,
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
    infe$setObservation(list(obs))

    ## run estimation
    infe$runEstimation()

    return(infe)
}

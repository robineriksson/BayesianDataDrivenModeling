library(SimInfInference)

##' Performes SLAM on the SISe dataset
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
##' @param B the number of boostrap replicates
##' @param dataDir path to the directory holding data used in simulation.
SLAMInference <- function(nStop = 1e3, nSim = 20,
                          debug = FALSE, solver = "ssm",
                          binary = FALSE,
                          Sw = 1e-3,
                          S = 1e-3, e = NULL,
                          normalize = TRUE, pertubate = TRUE,
    			  thetaTrue = c(upsilon = 0.01563841,
                                        beta_t1 = 0.14309175,
                                        beta_t2 = 0.13366590,
                                        beta_t3 = 0.15541914,
                                        gamma = 0.09720055,
                                        prev = 0.02537256),
                          threads = NULL,
                          bs = TRUE, obsSeed = 0, seed = 0,
                          logParam = FALSE, useW = TRUE, B = 100,
                          dataDir){
    if(!is.null(Sw))
        S <- Sw
    if(is.null(e))
        e <- S/100

    ## filepaths
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


    ## load d (distance between the observed nodes)
    load(paths["nObs"]) ## load: nObs

    ## The Simulator.
    Simulator <- SimInfSimulator_real

    if(binary){
        column <- "sample"
    } else {
        column <- "I"
    }
    logical <- FALSE


    cl <- NULL





    ## load data that is included the SimInf package
    events <- NULL
    u0 <- NULL


    ## initial phi level
    phiLevel <- "local"

    ## initial prevalence level
    if("prev" %in% names(thetaTrue))
        prevLevel = NULL
    else
        prevLevel = 0.02

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, solver = solver,
                               prevLevel = prevLevel, prevHerds = prevLevel, phiLevel = phiLevel,
                               u0 = u0, events = events, binary = binary, model = model,
                               nSim = nSim, obsDates = obs, dataDir = dataDir)


    ## The Summary statistics

    ##SummaryStatistics <- aggregateTS
    if(is.null(B))
        B <- nSim*20
    fun = sum
    extraArgsSummaryStatistics <- list(column = column, fun = fun, qtr = TRUE, bs = bs, B = B,
                                       useW = useW, logical = logical,
                                       fftcoeff = c(2,7))

    ## Proposal function, adaptive give NULL
    Proposal <- NULL
    extraArgsProposal <- NULL

    ## The Estimator
    Estimator <- SLAM

    set.seed(seed) ## set seed for parameter search.
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

    ## make observation
    set.seed(obsSeed) ## set seed for observation
    infe$changeExtraArgsSimulator(nSim = 1)
    obs <- infe$runSimulator(thetaTrue)
    infe$setObservation(obs)
    infe$changeExtraArgsSimulator(nSim = nSim)

    ## explore the posterior!
    set.seed(seed) ## reset seed for param estimation
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

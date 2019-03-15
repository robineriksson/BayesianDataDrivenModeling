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
##' @param useSMHI use SMHI based seasons.
SLAMInference <- function(nStop = 1e3, nSim = 20,
                          debug = FALSE, solver = "ssm",
                          Sw = NULL,
                          S = 1e-3, e = NULL,
                          normalize = TRUE, pertubate = TRUE,
                          thetaTrue = c(upsilon = 1.75e-2,
                                        beta_t1 = 1.57e-1,
                                        beta_t2 = 1.44e-1,
                                        beta_t3 = 1.50e-1,
                                        gamma = 0.1,
                                        prev = 0.02),
                          threads = NULL,
                          bs = TRUE, seed = 0,
                          logParam = FALSE, useW = TRUE, useSMHI = TRUE){
    if(!is.null(Sw))
        S <- Sw
    if(is.null(e))
        e <- S/100
    ## S = 1e-3 -> 20% accrate for inc. prev.
    ## e = 0.5
    ## S = 5e-5 -> 3.56%
    ## S = 5e-4 -> 1.4% or 2.2% acc.rate
    set.seed(seed) ## set up simulator

    load("~/Gits/BPD/R/DATA/secret/obsCleanDates.RData") ## loads: observation.dates

    ## We only to run the simulation for the maximum time that we have observations for.
    ## But we should start at the same time as the movements as they will affect
    ## the state of the distribution of induvidials at nodes.
    if(useSMHI)
        load("~/Gits/BPD/R/DATA/secret/SISe_smhi.rda") ## loads: model
    else
        load("~/Gits/BPD/R/DATA/secret/SISe.rda") ## loads: model

    tspan0 <- seq(head(model@events@time,1), tail(model@events@time,1), 1)
    realDates <- zoo::as.Date(tspan0, origin = "2005-01-01")#"2007-07-01")
    tinObs <- sort(unique(observation.dates$time))
    tObs <- as.numeric(tinObs) - (as.numeric(tinObs[1]) - tspan0[which(realDates == tinObs[1])])
    observation.dates$numTime <- as.numeric(observation.dates$time) - as.numeric(tinObs[1]) +
        tspan0[which(realDates == tinObs[1])]

    tspan <- tspan0[-which(realDates > tail(tinObs,1))]

    ## load names of the observed nodes -> nObs
    load("~/Gits/BPD/R/DATA/secret/nObs.RData") ## load: nObs


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
                               nSim = nSim, obsDates = observation.dates, useSMHI = useSMHI)


    ## The Summary statistics
    SummaryStatistics <- aggWilkinson

    extraArgsSummaryStatistics <- list(column = "sample", fun = mean, qtr = TRUE,  bs = bs,
                                       B = nSim*20, useW = useW, logical = FALSE)

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
    infe$setObservation(list(observation.dates))

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



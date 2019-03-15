library(SimInfInference)

##' Bootstrapped ABC (SLAM) estimate
##'
##' Generating Mboot chains,  then compute the mean bias.
##' @param nStop the length of the chain
##' @param Mboot the number of bootstraps
generateBootstrap <- function(nStop = 14000, Mboot = 10) {
	## generate the bootstrap samples
	boots = list()
	for (s in 1:Mboot) {
		infe <- bootSLreal(nStop, seed = s) 
		boots[[length(boots)+1]] <- infe
	}
	
	b = estBias(Mboot, 1000, 120)
	return(b)	
}


##' BootSLreal - a script that computes a SLAM chain
##' with the goal to estimate the bias of the "true" chain.
##'
##' @param nStop number of accepted parameter proposals
##' @param nSim number of simulations for each evaluation
##' @param debug if we want debug messages from the estimator
##' @param solver the stochastic solver to use in SimInf
##' @param binary if the data is binary filtered
##' @param S SLAM parameter 1
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
bootSLreal <- function(nStop = 100, nSim = 10,
                       debug = FALSE, solver = "ssm",
                       binary = TRUE,
                       normalize = TRUE,
                       thetaTrue = c(upsilon = 0.01748061,
                                     beta_t1 = 0.15210166,
                                     beta_t2 = 0.13849330,
                                     beta_t3 = 0.14714038,
                                     gamma = 0.10302909,
                                     prev = 0.01928063),
                       threads = NULL,
                       bs = TRUE, seed = 0,
                       logParam = FALSE, useW = TRUE, useSMHI = TRUE,
                       EstName = "SLAM",
                       S = 1e-3, trace = 10
                       ) {
    ## estimator specific parameters
    if(EstName == "SLAM") {
        e <- S/100
        EstimatorExtra= list(S = S, e = e, nStop = nStop,
                             C0 = 1e-2, n0 = 1, accVec = NULL,
                             reinit = FALSE, xbar = NULL, sigma = NULL)
    } else if(EstName == "maxSL") {
        EstimatorExtra = list(trace = trace, maxit = nStop)
    }


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

    ## old
    ##tObs <- seq(182+obsspan, 3287, obsspan) ## what nodes we observe
    ##


    ## load d (distance between the observed nodes)
    load("~/Gits/BPD/R/DATA/secret/d.RData") ## loads: distance matrix d
    nObs <- colnames(d)
    ## all nodes
    ##nodes <- 1:37221


    Simulator <- SimInfSimulator_real

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
    events <- NULL
    u0 <- NULL

    ##phiLevel = 0.5
    phiLevel <- "local"
    if("prev" %in% names(thetaTrue))
        prevLevel <- NULL
    else
        prevLevel <- 0.02

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, solver = solver,
                               prevLevel = prevLevel, prevHerds = prevLevel, phiLevel = phiLevel,
                               u0 = u0, events = events, binary = binary, model = model,
                               nSim = nSim, obsDates = observation.dates, useSMHI =  useSMHI)

    ## The Summary statistics
    SummaryStatistics <- aggWilkinson

    extraArgsSummaryStatistics <- list(column = column, fun = mean, qtr = TRUE, bs = bs, B = nSim*20,
                                       useW = useW, logical = logical)

    ## The Proposal When only estimating upsilon.
    Proposal <- NULL
    extraArgsProposal <- NULL

    ## The Estimator
    Estimator <- eval(get(EstName))

    ## if log of parameters
    if(logParam)
        theta0 <- log(thetaTrue)
    else
        theta0 <- thetaTrue



    extraArgsEstimator <- list(parameters = names(thetaTrue),
                               debug = debug, theta0 = theta0, nSim = nSim,
                               normalize = normalize, logParam = logParam)

    for (n in names(EstimatorExtra))
        extraArgsEstimator[n] = EstimatorExtra[n]



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

estBias <- function(B = 5, burnin = 200, thinning = 20) {
    thetaTrue <- c(upsilon = 0.01748061,
                   beta_t1 = 0.15210166,
                   beta_t2 = 0.13849330,
                   beta_t3 = 0.14714038,
                   gamma = 0.10302909,
                   prev = 0.01928063)

    dirname = "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/boot/sbin/"
    filename <- paste(dirname, "boot_seed_1.RData", sep = "")
    load(filename)
    s <- as.matrix(sboot_bin$getPosterior())
    sdim <- dim(s)
    s <- s[seq(burnin, sdim[1], thinning),]

    thetaVec <- array(0, dim = c(dim(s), B))
    thetaVec[,,1] <- s

    for (seed in 2:B) {
        filename <- paste(dirname, "boot_seed_", seed, ".RData", sep = "")
        load(filename)
        s <- as.matrix(sboot_bin$getPosterior())
        s <- s[seq(burnin, sdim[1], thinning),]
        thetaVec[,,seed] <- s
    }

    colnames(thetaVec) <- names(sboot_bin$getPosterior())

    meanEst <- t(apply(thetaVec[,1:length(thetaTrue),], c(2,3), mean))
    biasEst <- apply(meanEst, 2, mean) - thetaTrue


    thetaAll <- data.frame()
    for(b in 1:B)
        thetaAll <- rbind(thetaAll, thetaVec[,,b])



    return(list(bias = biasEst, mmean = meanEst, theta = thetaVec, theta0 = thetaTrue, thetaAll = thetaAll))
}

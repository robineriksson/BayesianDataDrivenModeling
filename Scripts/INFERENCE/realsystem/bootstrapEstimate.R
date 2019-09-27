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
		sboot_bin <- bootSLreal(nStop, seed = s)
		boots[[length(boots)+1]] <- sboot_bin
	}

	## estimate the bias using the list of runs.
	b = estBias(boots, Mboot, 1000, 120)
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
##' @param B the number of bootstrap replicates
##' @param dataDir the directory holding the epi-modeldata
bootSLreal <- function(nStop = 100, nSim = 10,
                       debug = FALSE, solver = "ssm",
                       binary = TRUE,
                       normalize = TRUE,
                       thetaTrue = c(upsilon = 0.01518903,
                                      beta_t1 = 0.14057693,
                                      beta_t2 = 0.12646374,
                                      beta_t3 = 0.15532180,
                                      gamma = 0.09650432,
                                      prev = 0.02471147),
                       threads = NULL,
                       bs = TRUE, seed = 0,
                       logParam = FALSE, useW = TRUE,
                       EstName = "SLAM",
                       S = 1e-3, trace = 10, B = 100, dataDir
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


    ## filepaths
    paths <- NULL
    paths["obs"] <- paste(dataDir, "obs.rda", sep="")
    paths["model"] <- paste(dataDir, "SISe.rda", sep="")
    paths["nObs"] <- paste(dataDir, "nObs.RData", sep="")



    load(paths["obs"]) ## loads: observation.dates
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

    tspan <- tspan0[-which(realDates > tail(tinObs,1))]

    ## load d (distance between the observed nodes)
    load(paths["nObs"]) ## load: nObs



    Simulator <- SimInfSimulator_real

    if(binary){
        ##for parallell sampling
        column <- "sample"

    } else {
        column <- "I"

    }
    logical <- FALSE

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
                               nSim = nSim, obsDates = obs, dataDir = dataDir)

    ## The Summary statistics
    SummaryStatistics <- aggWilkinson


    ## The Summary statistics
    SummaryStatistics <- aggWilkinson
    if(is.null(B))
        B <- nSim*20
    fun = sum
    extraArgsSummaryStatistics <- list(column = "sample", fun = fun, qtr = TRUE,  bs = bs,
                                       B = B, useW = useW, logical = FALSE, fftcoeff = c(2,7))

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

##' Estimate the bias from a list of inference runs.
##' @param infe.list a list of inference runs
##' @param B  the nubmer of inference runs to inclue
##' @param burnin what burnin to use for the Markov Chain
##' @param thinning what thinning to use for the Markov Chain
estBias <- function(infe.list, B = 10, burnin = 1000, thinning = 100) {
    thetaTrue <- c(upsilon = 0.01518903,
                   beta_t1 = 0.14057693,
                   beta_t2 = 0.12646374,
                   beta_t3 = 0.15532180,
                   gamma = 0.09650432,
                   prev = 0.02471147)

    sboot_bin <- infe.list[[1]]
    s <- as.matrix(sboot_bin$getPosterior())
    sdim <- dim(s)
    s <- s[seq(burnin, sdim[1], thinning),]

    thetaVec <- array(0, dim = c(dim(s), B))
    thetaVec[,,1] <- s

	## error check.
	if(B > length(infe.list)
		B <- length(infe.list)

	## loop for clarity.
    for (seed in 2:B) {
		sboot_bin <- infe.list[[seed]]
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

library(SimInfInference)


##' Performes SL MCMC on the SISe3 dataset
##' @param nStop number of accepted parameter proposals
##' @param nSim number of simulations for each evaluation
##' @param plotres true if result is plotted
SLMCMCInference <- function(nStop = 100, nSim = 25, plotres = FALSE, debug = FALSE,
                            solver = "aem", obsspan = 60, cl = NULL, column = "I",
                            binary = FALSE, Scaling = NULL,
                            ScalingVec =  c(0.03,0.005,0.03,0.03,0.2,0.2,0.01,0.01),
                            normalize = TRUE){
    set.seed(0) ## set up simulator
    tspan <- seq(1,4*365, 1)
    tObs <- seq(100,4*365, obsspan)
    data("nodes", package = "SimInf")
    nObs <- sample(1:length(rownames(nodes)),100)

    ## The Simulator.
    Simulator <- SimInfSimulator_sandbox

    ## for parallell sampling cNum <- parallel::detectCores() cl <-
    ## parallel::makeCluster(getOption("cl.cores", cNum))
    ## parallel::clusterExport(cl=cl,
    ## varlist=c("sample_herd","predict_env_sample_SISe","pool_prevalence","sample_pools"))


    ## load data that is included the SimInf package
    events <- SimInf::events_SISe()
    u0 <- SimInf::u0_SISe()

    phiLevel = 0.5
    prevLevel = 0.1
    extraArgsSimulator <- list(tspan = tspan, tObs = tObs,
                               nObs = nObs, runSeed = NULL, threads = NULL,
                               includeTrue = FALSE, solver = solver, prevLevel = prevLevel,
                               prevHerds = prevLevel, phiLevel = phiLevel, u0 = u0,
                               events = events, binary = binary, cl = cl)

    ## The Summary statistics
    SummaryStatistics <- SimInfStatistics_sandbox_full
    SummaryStatistics <- aggregateTS
    extraArgsSummaryStatistics <- list(column = column)

    ## The Proposal When only estimating upsilon.  S = 0.01 ->
    ## acceptance: ~0.31% S = 0.1, rho = 0.75, rhoBeta = 0.75 -> 5-11%
    ## S = 0.1, -> acceptance: ~ 0-28% S = 0.01 -> acceptance: ~ 40% S
    ## = 0.001 -> acceptance: ~ 22-93%
    Proposal <- Proposal_kernel_SISe_random_lognormal
    extraArgsProposal <- list(S = Scaling, Svec = ScalingVec, rho = 1, rhoBeta = 1)

    ## The Estimator
    Estimator <- SLMCMC


    ## True parameter
    ## upsilon <- 0.009
    ## beta_t1 <- 0.075
    ## beta_t2 <- 1.075*beta_t1
    ## gamma <- 0.1
    upsilon <- 0.0075
    beta_t1 <- 0.05
    beta_t2 <- 0.085 ## whic is = 1.69*beta_t1
    gamma <- 0.1


    thetaTrue <- c(upsilon = upsilon, beta_t1= beta_t1, beta_t2 = beta_t2, gamma = gamma)


    extraArgsEstimator <- list(parameters = c("upsilon", "beta_t1", "beta_t2", "gamma"),
                               nStop = nStop, nSim = nSim, debug = debug, theta0 = thetaTrue,
                               normalize = normalize)


    ## Create object!
    infe <- Inference$new(Simulator = Simulator, SummaryStatistics = SummaryStatistics,
                          Proposal = Proposal, Estimator = Estimator,
                          extraArgsSimulator = extraArgsSimulator,
                          extraArgsSummaryStatistics = extraArgsSummaryStatistics,
                          extraArgsEstimator = extraArgsEstimator,
                          extraArgsProposal = extraArgsProposal)

    ## make observation
    obs <- infe$runSimulator(thetaTrue)
    infe$setObservation(obs)

    ## explore the posterior!
    infe$runEstimation()

    return(infe)
}

weekend <- function(N1 = 3e3, N2 = 3e3, nSim = 25, Sw = 1e-3, pertubate = TRUE, bs = TRUE){
    if(N1 > 0){
        sfull <- SLAMInference(nStop = N1, nSim = nSim, binary = FALSE, Sw = Sw, pertubate = pertubate, bs = bs)
        save(sfull, file = "~/Gits/BPD/R/INFERENCE/1600system/output/slam/sfull.RData")
    }

    if(N2 > 0){
        sbin <- SLAMInference(nStop = N2, nSim = nSim, binary = TRUE, Sw = Sw, pertubate = pertubate, bs = bs)
        save(sbin, file = "~/Gits/BPD/R/INFERENCE/1600system/output/slam/sbin.RData")
    }
}

night <- function(N1 = 1e3, N2 = 1e3){
    if(N1 > 0){
        filename <- "~/Gits/BPD/R/INFERENCE/1600system/output/slam/sfull.RData"
        load(filename)
        sfull <- contInference(sfull, N1)
        save(sfull, file = filename)
    }

    if(N2 > 0){
        filename <- "~/Gits/BPD/R/INFERENCE/1600system/output/slam/sbin.RData"
        load(filename)
        sbin <- contInference(sbin, N2)
        save(sbin, file = filename)
    }
}

##' Performes SL MCMC on the SISe3 dataset
##' @param nStop number of accepted parameter proposals
##' @param nSim number of simulations for each evaluation
##' @param plotres true if result is plotted
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
                                       qtr = FALSE, bs = bs, B = 500,
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
    ## if(is.null(theta)){
    ##     posterior <- infe$getPosterior()
    ##     row <- dim(posterior)[1]
    ##     theta <- as.numeric(posterior[row, c("upsilon", "beta_t1", "beta_t2", "gamma")])
    ##     names(theta) <- c("upsilon", "beta_t1", "beta_t2", "gamma")
    ##     if(toLog)
    ##         theta <- log(theta)
    ## }

    accVec <- as.matrix(infe$getPosterior())
    d <- dim(accVec)
    if(infe$getExtraArgsEstimator()$logParam)
        accVec[,1:(d[2]-1)] <- log(accVec[,1:(d[2]-1)])

    infe$changeExtraArgsEstimator(nStop = nStop, accVec = accVec, reinit = TRUE)

    ## if(infe$getExtraArgsSimulator()$binary){
    ##     ##for parallell sampling
    ##     cNum <- parallel::detectCores()
    ##     cl <- parallel::makeCluster(getOption("cl.cores", cNum))
    ##     parallel::clusterExport(cl=cl,
    ##                             varlist=c("sample_herd",
    ##                                       "predict_env_sample_SISe",
    ##                                       "pool_prevalence",
    ##                                       "sample_pools"))
    ##     infe$changeExtraArgsSimulator(cl = cl)
    ## }

    infe$runEstimation()

    ## if(infe$getExtraArgsSimulator()$binary)
    ##     parallel::stopCluster(cl)

    return(infe)
}

pertubationSL <- function(thetalength = 21,
                          nSim = 20, pert = c(0.8,1.2),
                          thetaLength = 21,
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
    ##thetaTrue <- c(upsilon = 0.0050000, beta_t1 = 0.2302585, beta_t2 = 0.2302585, beta_t3 = 0.2302585, beta_t4 = 0.2302585)
    ##thetaTrue <- c(upsilon = 0.01, beta_t1 = 0.095, beta_t2 = 0.012, beta_t3 = 0.1, beta_t4 = 0.15)




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

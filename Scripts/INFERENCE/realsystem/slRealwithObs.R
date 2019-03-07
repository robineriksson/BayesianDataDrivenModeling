library(SimInfInference)

nightNodeRun <- function(N = 1e3, seed = NULL) {
    if(is.null(seed))
        seed <- sample(1:100,1)

    ##    upsilon    beta_t1    beta_t2    beta_t3      gamma       prev
    ##    0.01748061 0.15210166 0.13849330 0.14714038 0.10302909 0.01928063
    thetaTrue <- c(upsilon = 0.01748061,
                   beta_t1 = 0.15210166,
                   beta_t2 = 0.13849330,
                   beta_t3 = 0.14714038,
                   gamma = 0.10302909,
                   prev = 0.01928063)
    nSim <- 20

    Sw <- 1e-3
    bs <- TRUE


    dirname <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/sobs/"
    filename <- paste(dirname, "obs_seed_",seed,".RData", sep = "")
    print(filename)

    sobs <- SLAMInference(nStop = N, nSim = nSim, Sw = Sw, bs = bs, seed = seed)
    save(sobs,file=filename)



}

contNightNodeRun <- function(N = 5e2, seed) {
    dirname <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/sobs/"
    filename <- paste(dirname, "obs_seed_",seed,".RData", sep = "")
    print(filename)
    load(filename)
    sobs <- contInference(sobs, N)
    save(sobs,file=filename)
}

weekend <- function(N = 1e3, nSim = 20, normalize = TRUE, Sw = 1e-3, bs = TRUE, useW = TRUE ){
    if(N > 0) {
        sobs <- SLAMInference(nStop = N, nSim = nSim, normalize = normalize, Sw = Sw, bs = bs)
        filename <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sobs.RData"
        save(sobs, file = filename)
    }
}


night <- function(N = 1e3){
    if(N > 0){
        filename <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sobs.RData"
        load(filename)
        sobs <- contInference(sobs, N)
        save(sobs, file = filename)
    }
}


franken <- function(odd = FALSE, burnIn = 750) {
    dirname = "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/sobs/"
    if(odd)
        seeds <- c(2:6,8:9)
        #seeds <- c(3:11)
    else
        seeds <- c(2:10)

    filename0 <- paste(dirname, "obs_seed_1.RData", sep = "")
    load(filename0)
    dAll <- sobs$getPosterior()

    ## remove burnin
    dims <- dim(dAll)
    burnin <- burnIn:dims[1]

    dAll <- dAll[burnin,]
    for (s in seeds) {
        filename <- paste(dirname, "obs_seed_", s, ".RData", sep = "")
        load(filename)
        d <- sobs$getPosterior()[burnin,]
        dAll <- rbind(dAll, d)
    }

    return(dAll)
}

changePosterior <- function(infe, post){
    dims <- dim(post)
    post <- as.matrix(post)

    ## alter xbar
    xbar <- apply(post[,1:(dims[2]-1)], 2, mean)

    ## alter sigma
    cov.est <- function(X){
        d <- dim(X)
        if(is.null(d)){
            k <- length(X)
            xbar <- 1/(k+1) * sum(X)
            cov.est <- 1/k * (sum( X^2) - (k+1)*xbar^2)
        } else {
            k <- d[1]
            xbar <- 1/(k+1) * apply(X,2,sum)
            s <- t(X) %*% X - ((k+1)*xbar %*% t(xbar))
            cov.est <- s/k
        }
        return(cov.est)
    }

    C <- cov.est(post[,1:(dims[2]-1)])
    S <- infe$getExtraArgsEstimator()$S
    e <- infe$getExtraArgsEstimator()$e

    sigma <- S * C + S*e*diag(dims[2]-1)

    ## set the new parameters
    infe$changeExtraArgsEstimator(accVec = post, xbar = xbar, sigma = sigma, nStop = 1, reinit = TRUE)
    infe$runEstimation()

    return(infe)
}


##' Performes SL MCMC on the SISe3 dataset
##' @param nStop number of accepted parameter proposals
##' @param nSim number of simulations for each evaluation
##' @param plotres true if result is plotted
SLAMInference <- function(nStop = 100, nSim = 20,
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
    ##' TODO: rewrite to instead use:
    ## load("~/Gits/BPD/R/DATA/observation_aggregates.RData") # loads: data, aggregated observation.dates
    ## load("~/Gits/BPD/R/DATA/tinObs.RData") # load: tinObs, unique times at which we observe nodes.
    ## load("~/Gits/BPD/R/DATA/nodeTime.RData") # load: nodeTime, node and sample time.




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

    ## load d (distance between the observed nodes)
    ## load("~/Gits/BPD/R/DATA/secret/d.RData")
    ## nObs <- colnames(d)

    ## load names of the observed nodes -> nObs
    load("~/Gits/BPD/R/DATA/secret/nObs.RData") ## load: nObs


    ## The Simulator.
    Simulator <- SimInfSimulator_real

    ##for parallell sampling
    ## cNum <- parallel::detectCores()
    ## cl <- parallel::makeCluster(getOption("cl.cores", cNum))
    ## parallel::clusterExport(cl=cl,
    ##                         varlist=c("sample_herd",
    ##                                   "predict_env_sample_SISe",
    ##                                   "pool_prevalence",
    ##                                   "sample_pools"))
    ## cl <- NULL

    ## load data that is included the SimInf package
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

    return(infe)
}

reinitCluster <- function(infe){
    ##for parallell sampling
    cNum <- parallel::detectCores()
    cl <- parallel::makeCluster(getOption("cl.cores", cNum))
    parallel::clusterExport(cl=cl,
                            varlist=c("sample_herd",
                                      "predict_env_sample_SISe",
                                      "pool_prevalence",
                                      "sample_pools"))

    infe$changeExtraArgsSimulator(cl = cl)
    infe$changeExtraArgsSummaryStatistics(cl = cl)

    return(infe)
}

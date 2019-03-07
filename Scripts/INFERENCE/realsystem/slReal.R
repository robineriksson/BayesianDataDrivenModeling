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
    tspan <- seq(182, 3287, 1)
    tObs <- seq(182+obsspan, 3287, obsspan) ## what nodes we observe
    ## load d (distance between the observed nodes)
    load("~/Gits/BPD/R/DATA/secret/d.RData")
    nObs <- colnames(d)
    ## all nodes
    nodes <- 1:37221

    ## The Simulator.
    Simulator <- SimInfSimulator_real


    ## load data that is included the SimInf package
    events <- NULL
    u0 <- NULL

    ## ## for parallell sampling
    ## cNum <- parallel::detectCores()
    ## cl <- parallel::makeCluster(getOption("cl.cores", cNum))
    ## parallel::clusterExport(cl=cl,
    ##                         varlist=c("sample_herd","predict_env_sample_SISe","pool_prevalence","sample_pools"))

    ## loads variable: model
    load("~/Gits/BPD/R/DATA/secret/SISe.rda")
    phiLevel <- 0.05

    if("prev" %in% names(thetaTrue))
        prevLevel = NULL
    else
        prevLevel = 0.02

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, solver = solver,
                               prevLevel = prevLevel, prevHerds = prevLevel, phiLevel = phiLevel,
                               u0 = u0, events = events, binary = binary, cl = cl, model = model)


    ## The Summary statistics
    ##SummaryStatistics <- SimInfStatistics_sandbox_full
    SummaryStatistics <- aggQtrWilkinson
    extraArgsSummaryStatistics <- list(column = column, fun = mean)

    ## The Proposal
    Proposal <- Proposal_kernel_SISe_random_lognormal
    extraArgsProposal <- list(S = Scaling, Svec = ScalingVec, rho = 1, rhoBeta = 1)

    ## The Estimator
    Estimator <- SLMCMC

    ## make observation
    ## True parameter
    upsilon <- 0.0135
    beta_t1 <- 0.15
    beta_t2 <- 0.46*beta_t1 ##0.069
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

nightNodeRun <- function(N = 1e3, seed = NULL, binary = TRUE) {
    if(is.null(seed))
        seed <- sample(1:100,1)

    ##    upsilon    beta_t1    beta_t2    beta_t3      gamma       prev
    ##    0.01748061 0.15210166 0.13849330 0.14714038 0.10302909 0.01928063
    ## thetaTrue <- c(upsilon = 0.01748061,
    ##                beta_t1 = 0.15210166,
    ##                beta_t2 = 0.13849330,
    ##                beta_t3 = 0.14714038,
    ##                gamma = 0.10302909,
    ##                prev = 0.01928063)

     thetaTrue <- c(upsilon = 0.01754161,
                    beta_t1 = 0.15860932,
                    beta_t2 = 0.14621449,
                    beta_t3 = 0.15045281,
                    gamma = 0.10044851,
    		    prev = 0.02027267)


    nSim <- 20

    Sw <- 1e-3
    bs <- TRUE

    if(binary)
        dirname <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/sbin/"
    else
        dirname <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/sfull/"
    filename <- paste(dirname, "seed_",seed,".RData", sep = "")
    print(filename)

    s <- SLAMInference(nStop = N, nSim = nSim, binary = binary, Sw = Sw, bs = bs, seed = seed)
    if(binary) {
        sbin <- s
        save(sbin,file=filename)
    } else {
        sfull <- s
        save(sfull,file=filename)
    }

}

contNightNodeRun <- function(N = 5e2, seed) {
    dirname <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/sbin/"
    filename <- paste(dirname, "seed_",seed,".RData", sep = "")
    print(filename)
    load(filename)
    sbin <- contInference(sbin, N)
    save(sbin,file=filename)
}



weekend <- function(N1 = 3e3, N2= 3e3, nSim = 20, pertubate = TRUE, Sw = 1e-3, bs = TRUE, useW = TRUE, seed = 0){

    if(N1 > 0) {
        sfull <- SLAMInference(nStop = N1, nSim = nSim, binary = FALSE, pertubate = pertubate, Sw = Sw, bs = bs,
                               useW = useW, seed = seed)
        filename = "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sfull.RData"
        save(sfull, file = filename)
    }

    if(N2 > 0) {
        sbin <- SLAMInference(nStop = N2, nSim = nSim, binary = TRUE, pertubate = pertubate, Sw = Sw, bs = bs,
                              useW = useW, seed = seed)
        filename = "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sbin.RData"
        save(sbin, file = filename)
    }
}

night <- function(N1 = 1e3, N2 = 1e3){
    if(N1 > 0){
        filename = "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sfull.RData"
        load(filename)
        sfull <- contInference(sfull, N1)
        save(sfull, file = filename)
    }

    if(N2 > 0){
        filename = "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sbin.RData"
        load(filename)
        sbin <- contInference(sbin, N2)
        save(sbin, file = filename)
    }
}

multirun <- function(N) {
    s1 <- SLAMInference(N, nSim = 20, bs = TRUE)
    s2 <- SLAMInference(N, nSim = 20, binary = TRUE, bs = TRUE)
    return(list(sfull = s1, sbin = s2))
}


franken <- function(binary = TRUE, burnIn = 850){
    dirname = "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/"
    if(binary) {
        dirname <- paste(dirname,"sbin/",sep="")
    }
    else {
        dirname <- paste(dirname,"sfull/",sep="")
    }

    seeds <-  2:10

    filename0 <- paste(dirname, "seed_1.RData", sep = "")
    load(filename0)

    if(binary)
        infe <- sbin
    else
        infe <- sfull
    dAll <- infe$getPosterior()

    ## remove burnin
    dims <- dim(dAll)
    burnin <- burnIn:dims[1]

    dAll <- dAll[burnin,]
    for (s in seeds) {
        filename <- paste(dirname, "seed_", s, ".RData", sep = "")
        load(filename)

        if(binary)
            infe <- sbin
        else
            infe <- sfull

        d <- infe$getPosterior()[burnin,]
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
##' @param debug if we want debug messages from the estimator
##' @param column what the name of the data column is
##' @param binary if the data is binary filtered
##' @param S SLAM parameter 1
##' @param e SLAM parameter 2
##' @param normalize normalize the summary statistics before comparison
##' @param pertubate if to pertubate the initial guess or not.
SLAMInference <- function(nStop = 100, nSim = 20,
                          debug = FALSE, solver = "ssm",
                          binary = FALSE,
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
                          bs = FALSE, obsSeed = 0, seed = 0,
                          logParam = FALSE, useW = TRUE, useSMHI = TRUE){
    if(!is.null(Sw))
        S <- Sw
    if(is.null(e))
        e <- S/100




    ## old S = 1e-4 -->0.85% acceptance (using the divided approach)
    ## S = 1e-3 --> 60%
    ## S = 1e-2 --> 30-40%


    load("~/Gits/BPD/R/DATA/secret/obsCleanDates.RData") ## loads: observation.dates
    ##tinObs <- sort(unique(observation.dates$time))
    ##tObs <- as.numeric(tinObs) - 1xb3514
    ##tspan <- seq(182, 3287, 1)

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

    ## The Simulator.
    Simulator <- SimInfSimulator_real

    if(binary){
        column <- "sample"
    } else {
        column <- "I"
    }
    logical <- FALSE

    ##     if(is.null(threads))
    ##         cNum <- parallel::detectCores()
    ##     else
    ##         cNum <- threads
    ##     cl <- parallel::makeCluster(getOption("cl.cores", cNum))
    ##     parallel::clusterExport(cl=cl,
    ##                             varlist=c("sample_herd",
    ##                                       "predict_env_sample_SISe",
    ##                                       "pool_prevalence",
    ##                                       "sample_pools"))
    ## } else {
    ##     if(is.null(column))
    ##         column <- "I"
    cl <- NULL





    ## load data that is included the SimInf package
    events <- NULL
    u0 <- NULL

    ## loads variable: model
    load("~/Gits/BPD/R/DATA/secret/SISe.rda")
    ##phiLevel <- 0.05
    phiLevel <- "local"
    if("prev" %in% names(thetaTrue))
        prevLevel = NULL
    else
        prevLevel = 0.02

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, solver = solver,
                               prevLevel = prevLevel, prevHerds = prevLevel, phiLevel = phiLevel,
                               u0 = u0, events = events, binary = binary, model = model,
                               nSim = nSim, obsDates = observation.dates, useSMHI =  useSMHI)


    ## The Summary statistics

    ##SummaryStatistics <- aggregateTS
    SummaryStatistics <- aggWilkinson
    extraArgsSummaryStatistics <- list(column = column, fun = mean, qtr = TRUE, bs = bs, B = nSim*20,
                                       useW = useW, logical = logical)

    ## The Proposal When only estimating upsilon.  S = 0.01 ->
    ## acceptance: ~0.31% S = 0.1, rho = 0.75, rhoBeta = 0.75 -> 5-11%
    ## S = 0.1, -> acceptance: ~ 0-28% S = 0.01 -> acceptance: ~ 40% S
    ## = 0.001 -> acceptance: ~ 22-93%
    Proposal <- NULL
    extraArgsProposal <- NULL

    ## The Estimator
    Estimator <- SLAM


     ## make observation
    ## True parameter
    ## upsilon <- 0.0135
    ## beta_t1 <- 0.16
    ## beta_t2 <- 0.055
    ## gamma <- 0.1

    ##thetaTrue <- c(upsilon = upsilon, beta_t1= beta_t1, beta_t2 = beta_t2, gamma = gamma)

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


    ## r <- infe$runSimulator(thetaTrue)
    ## z <- infe$runStatistics(r)
    ## browser()
    ## browser()
    ## plotME(r, "I")


    ## explore the posterior!
    set.seed(seed) ## reset seed for param estimation
    infe$runEstimation()


    ##parallel::stopCluster(cl)

    return(infe)
}

##' plots the last accepted parameter pair
##' against a simulation using theta0
##' @param infe the inference object.
plotME<- function(infe) {
    nSim <- infe$getExtraArgsSimulator()$nSim

    thetaTrue <- c(upsilon = 1.75e-2,
                   beta_t1 = 1.57e-1,
                   beta_t2 = 1.44e-1,
                   beta_t3 = 1.50e-1,
                   gamma = 0.1,
                   prev = 0.02)

    infe$changeExtraArgsSimulator(nSim = 1)
    obs <- infe$runSimulator(thetaTrue)
    obs.p <- agg(obs, infe$getExtraArgsSummaryStatistics())

    infe$changeExtraArgsSimulator(nSim = nSim)
    thetaPost <- as.numeric(tail(infe$getPosterior(),1)[1:length(thetaTrue)])
    names(thetaPost) <- names(thetaTrue)
    post <- infe$runSimulator(thetaPost)
    post.p <- agg(post, infe$getExtraArgsSummaryStatistics())
    post.p <- matrix(post.p, ncol = nSim)


    ## plot!
    t = 1:length(obs.p)
    obsdata <- data.frame(t = t, observed = obs.p)
    simdata <- data.frame(t = t, post.p)
    simdata.mf <- reshape::melt(simdata, id.var="t")

    dribbon = data.frame(t = t, ymin = apply(post.p,1,min), ymax = apply(post.p,1,max))

    library(dplyr, warn.conflicts = FALSE)
    pplot <- ggplot2::ggplot()
    pplot <- pplot + ggplot2::geom_ribbon(data = dribbon,
                                          ggplot2::aes(x = t, ymin = ymin, ymax = ymax),
                                          alpha = 0.3)

    mapping <- ggplot2::aes(x = t, y = value, color = variable)

    pplot <- pplot + ggplot2::geom_line(ggplot2::aes(x=t, y=observed),
                                        color = "red", data = obsdata)
    pplot <- pplot + ggplot2::geom_point(ggplot2::aes(x=t, y=value, color = variable),
                                         colour = "blue", alpha = 0.3, data = simdata.mf)

    plot(pplot)

    return(list(obs = obs.p, post = post.p, pplot = pplot))
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

    ## if(infe$getExtraArgsSimulator()$binary)
    ##     parallel::stopCluster(cl)

    return(infe)
}

## reinitCluster <- function(infe){
##     ##for parallell sampling
##     cNum <- parallel::detectCores()
##     cl <- parallel::makeCluster(getOption("cl.cores", cNum))
##     parallel::clusterExport(cl=cl,
##                             varlist=c("sample_herd",
##                                       "predict_env_sample_SISe",
##                                       "pool_prevalence",
##                                       "sample_pools"))

##     infe$changeExtraArgsSimulator(cl = cl)
##     infe$changeExtraArgsSummaryStatistics(cl = cl)

##     return(infe)
## }

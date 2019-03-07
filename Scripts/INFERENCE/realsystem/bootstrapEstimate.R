library(SimInfInference)

nightNodeRun <- function(nStop = 100, seed = NULL, method = "SLAM") {
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

    nSim <- 10
    binary = TRUE
    bs <- TRUE

    dirname <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/boot/sbin/"

    filename <- paste(dirname, "boot_seed_",seed,".RData", sep = "")
    print(filename)
    s <- bootSLreal(nSim = nSim, nStop = nStop, binary = binary, bs = bs, EstName = method,
                    seed = seed)

    sboot_bin <- s
    save(sboot_bin,file=filename)
}

contNightNodeRun <- function(N = 5e2, seed) {
    dirname <- "~/Gits/BPD/R/INFERENCE/realsystem/output/slam/batchRun/boot/sbin/"
    filename <- paste(dirname, "boot_seed_",seed,".RData", sep = "")
    print(filename)
    load(filename)
    sboot_bin <- contInference(sboot_bin, N)
    save(sboot_bin,file=filename)
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

multiEstimate <- function(B = 5, nSim = 10, nStop = 25, method = "SLAM") {
    thetaTrue = c(upsilon = 0.0075,
                  beta_t1 = 0.050,
                  beta_t2 = 0.085,
                  gamma = 0.1)

    if(method == "SLAM")
        thetaVec <- array(0, dim = c(nStop, length(thetaTrue)+1, B))
    else {
        ## method = "maxSL"
        thetaVec <- matrix(0, nrow = B, ncol = length(thetaTrue)+1)
    }

    for (seed in 1:B) {
        s <- bootSL1600(nSim = nSim, nStop = nStop, thetaTrue = thetaTrue, seed = seed, EstName = method)

        if(method == "SLAM") {
            thetaVec[,,seed] <- as.matrix(s$getPosterior())
        } else if(method == "maxSL") {
            thetaVec[seed,] <- as.numeric(s$getPosterior())
        }

    }

    colnames(thetaVec) <- names(s$getPosterior())


    if(method == "SLAM") {
        meanEst <- t(apply(thetaVec[,1:length(thetaTrue),], c(2,3), mean))
        biasEst <- apply(meanEst, 2, mean) - thetaTrue
    } else if(method == "maxSL") {
        meanEst <- apply(thetaVec[,1:length(thetaTrue)], 2, mean)
        biasEst <- meanEst - thetaTrue
    }

    return(list(bias = biasEst, mmean = meanEst, theta = thetaVec, theta0 = thetaTrue))
}

plotBootPost <- function(fulldata) {
    data <- fulldata$theta
    dataAll <- fulldata$thetaAll
    truth <- fulldata$theta0

    dims <- dim(data)
    rows <- dims[1]
    cols <- dims[2]-1
    B <- dims[3]

    dataMean <- apply(data, 2, mean)

    par(mfrow = c(3,2))
    for(col in 1:cols) {

        xlim <- c( min(data[,col,]), max(data[,col,]))

        d <- density(dataAll[,col])
        ylim <- c(0,max(d$y)*8)

        plot(d, lwd = 2, xlim = xlim, ylim = ylim,
             main = colnames(data)[col])
        if(!is.null(truth))
            abline(v = truth[col], lwd = 3)

        abline(v = dataMean[col], lwd = 3, lty = 2, col = "red")

        for(b in 1:B) {
            lines(density(data[,col,b]),col = b+1, lty = b)
        }
    }
}

plotBootPost2 <- function(data) {
    mean.norm <- sweep(data$mmean,2,data$theta0,"/")
    dims <- dim(data$theta)
    sl <- data$theta[,dims[2],]

    par(mfrow = c(2,1))
    boxplot(mean.norm, main = "posterior mean / true <-- over 10 bootstrap estimates")

    plot(c(0,dim(sl)[1]), c(min(apply(sl,2,min)), max(apply(sl,2,max))) ,type = "n",
         xlab = "index", ylab = "(log) SL", main = "SL trajectories for all bootstrap estimates")
    for(s in 1:dim(sl)[2])
        lines(sl[,s], col = s)


}


bootSL1600 <- function(nSim = 20,
                       debug = FALSE, solver = "ssm", obsspan = 60,
                       binary = TRUE,
                       normalize = TRUE,
                       thetaTrue = c(upsilon = 0.0075,
                                     beta_t1 = 0.050,
                                     beta_t2 = 0.085,
                                     gamma = 0.1),
                   threads = NULL,
                   bs = TRUE, seed = 0, logParam = FALSE,
                   EstName = "maxSL",
                   S = 1e-3, nStop = 100, trace = 10
                   ) {

    if(EstName == "SLAM") {
        e <- S/100
        EstimatorExtra= list(S = S, e = e, nStop = nStop,
                              C0 = 1e-2, n0 = 1, accVec = NULL,
                              reinit = FALSE, xbar = NULL, sigma = NULL)
    } else if(EstName == "maxSL") {
        EstimatorExtra = list(trace = trace, maxit = nStop)
    }

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
    extraArgsSummaryStatistics <- list(column = column, fun = mean, qtr = FALSE, bs = bs, B = nSim*20,
                                       useW = FALSE, logical = logical)

    ## The Proposal When only estimating upsilon.
    Proposal <- NULL
    extraArgsProposal <- NULL

    ## The Estimator
    Estimator <- eval(get(EstName))



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

bootSLreal <- function(nSim = 20,
                       debug = FALSE, solver = "ssm", obsspan = 60,
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
                       EstName = "maxSL",
                       S = 1e-3, nStop = 100, trace = 10
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

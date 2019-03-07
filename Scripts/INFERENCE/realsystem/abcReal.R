library(SimInfInference)

ABCInference <- function(nStop = 100, epsilon = 0.1, method = "rejection", debug = FALSE, plotres = FALSE, solver = "aem", column = "I", obsspan = 60, cl = NULL, binary = FALSE){
    set.seed(0)
    ## set up simulator
    tspan <- seq(182, 3287, 1)
    tObs <- seq(182+obsspan, 3287, obsspan)

    ## what nodes we observe
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
    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, solver = solver,
                               prevLevel = NULL, prevHerds = NULL, u0 = u0,
                               events = events, binary = binary, cl = cl, model = model)
    ## The Summary statistics
    ##SummaryStatistics <- SimInfStatistics_sandbox_full
    SummaryStatistics <- aggregateTS
    extraArgsSummaryStatistics <- list(column = column)

    ## The Proposal
    ## When only estimating upsilon.
    Proposal <- dummy_Proposal_SISe()
    ##Proposal <- Proposal_SISe3sp_fix()

    ## The Estimator
    Estimator <- ABC_pkg
    extraArgsEstimator <- list(parameters = c("upsilon", "beta_t1", "beta_t2", "gamma"),
                               nStop = nStop, debug = debug, epsilon = epsilon, method = method)



    ## Create object!
    infe <- Inference$new(Simulator = Simulator, SummaryStatistics = SummaryStatistics,
                          Proposal = Proposal, Estimator = Estimator,
                          extraArgsSimulator = extraArgsSimulator,
                          extraArgsSummaryStatistics = extraArgsSummaryStatistics,
                          extraArgsEstimator = extraArgsEstimator)

    ## make observation
    ## True parameter
    upsilon <- 0.0120
    beta_t1 <- 0.1
    beta_t2 <- 1.075*beta_t1
    gamma <- 0.1

    thetaTrue <- c(upsilon = upsilon, beta_t1= beta_t1, beta_t2 = beta_t2, gamma = gamma)

    sim0 <- infe$runSimulator(thetaTrue)

    infe$setObservation(sim0)

    ## run inference!
    infe$runEstimation()

    ## plot the posterior!
    if(plotres)
        infe$plot()

    return(infe)
}


##' Test the area around the true parameter space.
##' Idea is that we will find an "extrema" that
##' will be possible for the inference algorithms to fnid.
pertubationABC <- function(thetalength = 10, thetaVec = NULL, epsilon = 1, method = "rejection",
                           pert = c(0.4,0.4,0.4,0.4), solver = "aem", numClusters = 1, SI = "SInomeanS",
                           obsspan = 60, SSspan = 90, betalag = c(2,3,4), betapower = c(1,1,1)){
    set.seed(0)
    ## set up simulator
    tspan <- seq(182, 3287, 1)
    tObs <- seq(182+obsspan, 3287, obsspan)

    ## what nodes we observe
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

    ## for parallell sampling
    cNum <- parallel::detectCores()
    cl <- parallel::makeCluster(getOption("cl.cores", cNum))
    parallel::clusterExport(cl=cl,
                            varlist=c("sample_herd","predict_env_sample_SISe","pool_prevalence","sample_pools"))

    ## loads variable: model
    load("~/Gits/BPD/R/DATA/secret/SISe.rda")
    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, solver = solver,
                               prevLevel = NULL, prevHerds = NULL, u0 = u0,
                               events = events, binary = FALSE, cl = cl, model = model)
    ## The Summary statistics
    SummaryStatistics <- aggregateTS ##SimInfStatistics_real_full

    ##extraArgsSummaryStatistics <- list(nodes = nObs, numClusters = numClusters, SI = SI, span = SSspan,
    ##                                   betalag = betalag, betapower = betapower, distance = d)

    extraArgsSummaryStatistics <- list(column = "I")
    ## The Proposal
    ## When only estimating upsilon
    Proposal <- NULL

    ## The Estimator
    Estimator <- NULL
    extraArgsEstimator <- NULL

    ## Create object!
    infe <- Inference$new(Simulator = Simulator, SummaryStatistics = SummaryStatistics,
                          Proposal = Proposal, Estimator = Estimator,
                          extraArgsSimulator = extraArgsSimulator,
                          extraArgsSummaryStatistics = extraArgsSummaryStatistics,
                          extraArgsEstimator = extraArgsEstimator)

    ## make observation
    thetaTrue <- c(upsilon = 0.0270000, beta_t1 = 0.2302585, beta_t2 = 0.2302585, beta_t3 = 0.2302585, beta_t4 = 0.2302585)
    ##thetaTrue <- c(upsilon = 0.0108, beta_t1 = 0.095, beta_t2 = 0.012, beta_t3 = 0.1, beta_t4 = 0.15)

    s1 <- system.time(
        sim0 <- infe$runSimulator(thetaTrue)
    )

    s2 <- system.time(
        ss0 <- infe$runStatistics(sim0)
    )
    ##s2 <- NULL

    print("simulator time")
    print(s1)
    print("summary time")
    print(s2)


    ## create grid of all possible combinations.
    if(is.null(thetaVec)){
        ## do only for 1 beta. We assume that they behave similarly
        upsilon <- seq(from = (1-pert[1])*thetaTrue[["upsilon"]], to = (1+pert[2])*thetaTrue[["upsilon"]],
                       length.out = thetalength)
        beta <- seq(from =  (1-pert[3])*thetaTrue[["beta_t1"]], to =  (1+pert[4])*thetaTrue[["beta_t1"]],
                    length.out = thetalength)

        ## create a parameter grid to traverse.
        thetaVec <- expand.grid(upsilon, beta)
        colnames(thetaVec) <- c("upsilon", "beta")
    }

    len <- dim(thetaVec)[1]
    print(sprintf("Will do %d iterations", len))
    print(sprintf("Estimated time to completion: %f hours (or %f minutes)", (s1+s2)[3]*len/60/60, (s1+s2)[3]*len/60))


    ss <- matrix(numeric(length(ss0)*len), ncol = length(ss0))
    ## create simulations for all combinations.
    for(i in 1:len){
        theta <- c(upsilon = thetaVec$upsilon[i], beta_t1 = thetaVec$beta[i], beta_t2 = 0.012, beta_t3 = 0.1, beta_t4 = 0.15)
        sim <- infe$runSimulator(theta)
        ss[i,] <- infe$runStatistics(sim)
    }

    ## stop cluster.
    parallel::stopCluster(cl)

    colnames(ss) <-  names(ss0)

    abcObj <- abc::abc(target = ss0, param = thetaVec, sumstat = ss, method = "rejection", tol = epsilon)
    dist <- abcObj$dist

    dat <- data.frame(thetaVec, dist)
    ## require(lattice)
    ## contourplot(dist ~ upsilon*beta, data = dat, region = TRUE)
    return(list(dat = dat, ss0 = ss0))
}




## plot3dStuff <-  function(data, invert = FALSE, ex = TRUE){
##     require(akima); require(rgl); require(plot3D);
##     x = data[,1]
##     y = data[,2]
##      if(ex)
##         data[,3] <- exp(data[,3])
##     if(invert){
##         z = 1/data[,3]
##         if(Inf %in% z)
##             z = 1/(1+data[,3])
##         zlab <- "1/distance"
##     } else {
##         z = data[,3]
##         zlab <- "distance"
##     }
##     s=interp(x,y,z)
##     persp3D(s$x,s$y,s$z,
##             xlab = expression(upsilon), ylab = expression(beta[1]), zlab = zlab)



## }

## plot2dStuff <- function(data){
##     require(lattice)
##     contourplot(dist ~ upsilon*beta, data = data, region = TRUE)
## }

## plot2dotherStuff <-  function(data, invert = FALSE, ex = TRUE){
##     require(akima);require(plot3D);
##     x = data[,1]
##     y = data[,2]
##      if(ex)
##         data[,3] <- exp(data[,3])
##     if(invert){
##         z = 1/data[,3]
##         if(Inf %in% z)
##             z = 1/(1+data[,3])
##     } else
##         z = data[,3]

##     s=interp(x,y,z)
##     image2D(x = s$x, y = s$y, z = s$z,
##             contoure = list(levels = 1, col = "black", lwd = 2),
##             xlab = expression(upsilon), ylab = expression(beta[1]), main = "1/distance")
##     rect(0.0095, 0.090, 0.0105, 0.10)
## }

## plot2and3 <-  function(data, invert = TRUE, ex = FALSE){
##     require(akima); require(rgl); require(plot3D);

##     x = data[,1]
##     y = data[,2]

##     if(ex)
##         data[,3] <- exp(data[,3])

##     if(invert){
##         z = 1/data[,3]
##         if(Inf %in% z)
##             z = 1/(1+data[,3])
##         zlab <- "1/distance"
##     } else {
##         z = data[,3]
##         zlab <- "distance"
##     }

##     s=interp(x,y,z)
##     ##X11()
##     par(mfrow = c(1,2), mar = c(10,1,10,1))
##     image2D(x = s$x, y = s$y, z = s$z,
##             contoure = list(levels = 1, col = "black", lwd = 2),
##             xlab = expression(upsilon), ylab = expression(beta[1]), main = "1/distance")
##     rect(0.0095, 0.090, 0.0105, 0.10)
##     persp3D(s$x,s$y,s$z,
##             xlab = expression(upsilon), ylab = expression(beta[1]), zlab = zlab)

## }

##' when applying -inf the vector
##' seems to become of factor form.
##' This function cleans this away.
cleanData <- function(data){
    data$dist <- as.numeric(levels(data$dist)[data$dist])
    data
}

##library(SimInf)
library(SimInfInference)

weekend <- function(N1 = 3e3, N2 = 3e3, epsilon = NULL){
    if(N1 > 0){
        abcfull <- ABCInference(nStop = N1, column = "I", binary = FALSE, epsilon = epsilon)
        save(abcfull, file = "~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcfull.RData")
    }

    if(N2 > 0){
        abcbin <- ABCInference(nStop = N2, column = "sample", binary = TRUE, epsilon = epsilon)
        save(abcbin, file = "~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcbin.RData")
    }
}

night <- function(N1 = 1e3, N2 = 1e3){
    if(N1 > 0){
        filename <- "~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcfull.RData"
        load(filename)
        abcfull <- contInference(abcfull, N1)
        save(abcfull, file = filename)
    }

    if(N2 > 0){
        filename <- "~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcbin.RData"
        load(filename)
        abcbin <- contInference(abcbin, N2)
        save(abcbin, file = filename)
    }
}

ABCInference <- function(nStop = 100, epsilon = 0.1,
                         method = "rejection", debug = FALSE,
                         solver = "ssm", column = "I", obsspan = 60,
                         binary = FALSE, distance = euclid,
                         thetaTrue = c(upsilon <- 0.009
                                       beta_t1 <- 0.075
                                       beta_t2 <- 0.085
                                       gamma <- 0.1)
                         ){
    set.seed(0)
    ## set up simulator
    tspan <- seq(1,4*365, 1)
    tObs <- seq(100,4*365, obsspan)
    ##tObs <- seq(1,4*365, 42)
    data("nodes", package = "SimInf")
    ##nObs <- sample(as.numeric(rownames(nodes)),100)
    nObs <- sample(1:length(rownames(nodes)),100)

    ## The Simulator.
    Simulator <- SimInfSimulator_sandbox

    if(binary){
        ##for parallell sampling
        cNum <- parallel::detectCores()
        cl <- parallel::makeCluster(getOption("cl.cores", cNum))
        parallel::clusterExport(cl=cl,
                                varlist=c("sample_herd",
                                          "predict_env_sample_SISe",
                                          "pool_prevalence",
                                          "sample_pools"))
    } else
        cl <- NULL

    ## load data that is included the SimInf package
    events <- SimInf::events_SISe()
    u0 <- SimInf::u0_SISe()

    phiLevel = 0.5
    prevLevel = NULL
    extraArgsSimulator <- list(tspan = tspan, tObs = tObs,
                               nObs = nObs, runSeed = NULL, threads = NULL,
                               includeTrue = FALSE, solver = solver, prevLevel = prevLevel,
                               prevHerds = prevLevel, phiLevel = phiLevel, u0 = u0,
                               events = events, binary = binary, cl = cl, nSim = 1)
    ## The Summary statistics
    ##SummaryStatistics <- SimInfStatistics_sandbox_full
    SummaryStatistics <- aggregateTS
    extraArgsSummaryStatistics <- list(column = column)

    ## The Proposal
    ## When only estimating upsilon.
    Proposal <- Proposal_SISe_allUniformNarrowBeta3Prev
    extraArgsProposal <- list()
    ##Proposal <- Proposal_SISe3sp_fix()

    ## The Estimator
    Estimator <- ABC##_pkg

    ## make observation
    ## True parameter
    ## upsilon <- 0.009
    ## beta_t1 <- 0.075
    ## beta_t2 <- 1.075*beta_t1
    ## gamma <- 0.1
    ## upsilon <- 0.0075
    ## beta_t1 <- 0.05
    ## beta_t2 <- 0.085 ## whic is = 1.69*beta_t1
    ## gamma <- 0.1
    ##thetaTrue <- c(upsilon = upsilon, beta_t1 = beta_t1, beta_t2 = beta_t2, gamma = gamma)

    extraArgsEstimator <- list(parameters = names(thetaTrue),
                               nStop = nStop, debug = debug, epsilon = epsilon, method = method,
                               distance = distance, theta0 = thetaTrue, accVec = NULL)


    ## Create object!
    infe <- Inference$new(Simulator = Simulator, SummaryStatistics = SummaryStatistics,
                          Proposal = Proposal, Estimator = Estimator,
                          extraArgsSimulator = extraArgsSimulator,
                          extraArgsSummaryStatistics = extraArgsSummaryStatistics,
                          extraArgsEstimator = extraArgsEstimator)

    sim0 <- infe$runSimulator(thetaTrue)

    infe$setObservation(sim0)

    ## run inference!
    infe$runEstimation()

    return(infe)
}

modCanberra <- function(x,y){2 * sqrt( sum( ((x-y)/(x+y))^2) )}

euclid <- function(x,y){sqrt( sum(((x-y)/x)^2) )}


##' Continue the already created inference class
##'
##' @param infe the inference class
##' @param nStop number of added parameters
##' @param theta a proposed theta0
##' @return inference class
contInference <- function(infe, nStop = 100){
    accVec <- as.matrix(infe$getPosterior())

    infe$changeExtraArgsEstimator(nStop = nStop, accVec = accVec)

    if(infe$getExtraArgsSimulator()$binary){
        ##for parallell sampling
        cNum <- parallel::detectCores()
        cl <- parallel::makeCluster(getOption("cl.cores", cNum))
        parallel::clusterExport(cl=cl,
                                varlist=c("sample_herd",
                                          "predict_env_sample_SISe",
                                          "pool_prevalence",
                                          "sample_pools"))
        infe$changeExtraArgsSimulator(cl = cl)
    }

    infe$runEstimation()

    if(infe$getExtraArgsSimulator()$binary)
        parallel::stopCluster(cl)

    return(infe)
}

##' Test the area around the true parameter space.
##' Idea is that we will find an "extrema" that
##' will be possible for the inference algorithms to fnid.
pertubationABC <- function(thetalength = 10, thetaVec = NULL, epsilon = 1, method = "rejection",
                           pert = c(0.4,0.4,0.4,0.4), solver = "aem", numClusters = 1,
                           obsspan = 60, column = "I"){
    set.seed(0)
    ## set up simulator
    tspan <- seq(1,4*365, 1)
    tObs <- seq(200,4*365, obsspan)
    ##tObs <- seq(1,4*365, 42)
    data("nodes", package = "SimInf")
    ##nObs <- sample(as.numeric(rownames(nodes)),100)
    nObs <- sample(1:length(rownames(nodes)),100)
    ##nObs <- 1:1600


    ## The Simulator.
    Simulator <- SimInfSimulator_sandbox

    ## for parallell sampling
    cNum <- parallel::detectCores()
    cl <- parallel::makeCluster(getOption("cl.cores", cNum))
    parallel::clusterExport(cl=cl,
                  varlist=c("sample_herd","predict_env_sample_SISe","pool_prevalence","sample_pools"))


    ## load data that is included the SimInf package
    events <- SimInf::events_SISe()
    u0 <- SimInf::u0_SISe()

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                               runSeed = NULL, threads = NULL, includeTrue = FALSE,
                               solver = solver, prevLevel = NULL, prevHerds = NULL,
                               u0 = u0, events = events, binary = FALSE, cl = cl)
    ## The Summary statistics
    ##SummaryStatistics <- SimInfStatistics_sandbox_full
    SummaryStatistics <- aggregateTS
    extraArgsSummaryStatistics <- list(column = column)
    ##extraArgsSummaryStatistics <- list(nodes = nodes, numClusters = numClusters, SI = SI, span = SSspan,
    ##betalag = betalag, betapower = betapower)

    ## The Proposal
    ## When only estimating upsilon
    Proposal <- NULL

    ## The Estimator
    Estimator <- abcGrid

    ## True parameter
    upsilon <- 0.009
    beta_t1 <- 0.075
    beta_t2 <- 1.075*beta_t1
    gamma <- 0.1

    thetaTrue <- c(upsilon = upsilon, beta_t1= beta_t1, beta_t2 = beta_t2, gamma = gamma)

    extraArgsEstimator <- list(theta0 = thetaTrue, method = "rejection", epsilon = epsilon, pert = pert, thetalength = thetalength)

    ## Create object!
    infe <- Inference$new(Simulator = Simulator, SummaryStatistics = SummaryStatistics,
                          Proposal = Proposal, Estimator = Estimator,
                          extraArgsSimulator = extraArgsSimulator,
                          extraArgsSummaryStatistics = extraArgsSummaryStatistics,
                          extraArgsEstimator = extraArgsEstimator)

    ## make observation
    #thetaTrue <- c(upsilon = 0.0050000, beta_t1 = 0.2302585, beta_t2 = 0.2302585, beta_t3 = 0.2302585, beta_t4 = 0.2302585)

    ##thetaTrue <- c(upsilon = 0.01, beta_t1 = 0.095, beta_t2 = 0.012, beta_t3 = 0.1, beta_t4 = 0.15)


    sim0 <- infe$runSimulator(thetaTrue)

    infe$setObservation(sim0)

    infe$runEstimation()

    return(infe)


}

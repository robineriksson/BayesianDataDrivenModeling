library(SimInfInference)


##' Perform ABC rejection sampling on the 1600 network.
##'
##' @param nStop the number of samples
##' @param epsilon the proportion of accepted proposals
##' @param solver the stochastic solver
##' @param obsspan how far between are the observations
##' @param binary is the data subjected to the filter or not?
##' @param distance the ABC distance function
##' @param thetaTrue the true parameter used for observation generation.
##' @return inference object.
ABCInference <- function(nStop = 100, epsilon = 0.1,
                         method = "rejection", debug = FALSE,
                         solver = "ssm", obsspan = 60,
                         binary = FALSE, distance = euclid,
                         thetaTrue = c(upsilon = 0.005,
                                       beta_t1 = 0.025,
                                       beta_t2 = 0.058,
                                       gamma = 0.1)
                         ){
    set.seed(0)

    ## set up simulator
    tspan <- seq(1,4*365, 1)
    tObs <- seq(100,4*365, obsspan)
    data("nodes", package = "SimInf")
    nObs <- sample(1:length(rownames(nodes)),100)

    ## The Simulator.
    Simulator <- SimInfSimulator_sandbox

    if(binary){
        column <- "sample"
    } else {
        column <- "I"
    }

    logical <- FALSE
    cl <- NULL


    ## load data that is included the SimInf package
    events <- SimInf::events_SISe()
    u0 <- SimInf::u0_SISe()

    phiLevel <- "local"
    if("prev" %in% names(thetaTrue))
        prevLevel <- NULL
    else
        prevLevel <- 0.1

    extraArgsSimulator <- list(tspan = tspan, tObs = tObs,
                               nObs = nObs, runSeed = NULL, threads = NULL,
                               includeTrue = FALSE, solver = solver, prevLevel = prevLevel,
                               prevHerds = prevLevel, phiLevel = phiLevel, u0 = u0,
                               events = events, binary = binary, cl = cl, nSim = 1)


    ## The Summary statistics
    SummaryStatistics <- aggWilkinson
    extraArgsSummaryStatistics <- list(column = column, fun = mean,
                                        qtr = FALSE, bs = FALSE, B = NULL,
                                        useW = FALSE, logical = logical)

    ## The Proposal
    ## When only estimating upsilon.
    Proposal <- Proposal_SISe_allUniformNarrow
    extraArgsProposal <- list()


    ## The Estimator
    Estimator <- ABC##_pkg

    ## make observation
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

## a modified Canberra distance function
modCanberra <- function(x,y){2 * sqrt( sum( ((x-y)/(x+y))^2) )}

## a normalized euclidian distance function
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

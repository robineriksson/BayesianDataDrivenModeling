##' This file holds the (R6) inference class for
##' parameters estimation and posterior exploration.
##'



##library(R6)

##' An inference class that will hold the functionality of all the functions.
##'
##' @export
Inference <- R6::R6Class("Inference",
                         lock_objects = FALSE,
                         public = list(
                             initialize = function(Simulator, SummaryStatistics, Estimator, Proposal,
                                                   extraArgsSimulator = NULL,
                                                   extraArgsSummaryStatistics = NULL,
                                                   extraArgsEstimator = NULL,
                                                   extraArgsProposal = NULL){

                                 private$Simulator <- Simulator
                                 private$SummaryStatistics <- SummaryStatistics
                                 private$Estimator <- Estimator
                                 private$Proposal <- Proposal

                                 private$setNames(Simulator, SummaryStatistics, Estimator, Proposal)

                                 private$extraArgsSimulator <- extraArgsSimulator
                                 private$extraArgsSummaryStatistics <- extraArgsSummaryStatistics
                                 private$extraArgsEstimator <- extraArgsEstimator
                                 private$extraArgsProposal <- extraArgsProposal
                             },
                             ## Check if summary statistics are normally distributed for N simulations.
                             normalTest = function(theta, N = NULL){
                                 ## Get the nSim variable from the estimator argument
                                 nSim <- private$extraArgsEstimator$nSim

                                 ## Check for normality!
                                 ss.list <- list()
                                 ## Simulate data
                                 if(!is.null(N))
                                     self$changeExtraArgsSimulator(nSim = N)
                                 else
                                     N <- private$extraArgsSimulator$nSim

                                 data <- self$runSimulator(theta=theta)
                                 ## Extract Summary Statistics from data!
                                 ss.c <- private$SummaryStatistics(data = data, extraArgs = private$extraArgsSummaryStatistics)
                                 if(!is.null(N))
                                     self$changeExtraArgsSimulator(nSim = nSim)


                                 ss.m <- t(matrix(ss.c, ncol = N))

                                 shap <- apply(ss.m,2, function(x){shapiro.test(x)$p.value} > 0.05)
                                 print(shap)
                                 ##jb <- apply(ss.m,2, function(x){jb.norm.test(x)$p.value} > 0.05)
                                 ## QQ plots
                                 len <- dim(ss.m)[2]
                                 numCols = ceiling(len/2)
                                 par(mfrow = c(2,numCols))
                                 for(i in 1:len){
                                     qqnorm(ss.m[,i], main = names(ss.m)[i])
                                     qqline(ss.m[,i])
                                 }
                                 return(shap.test = shap)
                             },
                             ## use the estimator with the supplied functions and arguments.
                             runEstimation = function(){
                                 self$tic <- Sys.time()
                                 est.out <- private$Estimator(observation = private$observation,
                                                              Simulator = private$Simulator,
                                                              SummaryStatistics = private$SummaryStatistics,
                                                              Proposal = private$Proposal,
                                                              extraArgsEstimator = private$extraArgsEstimator,
                                                              extraArgsSimulator = private$extraArgsSimulator,
                                                              extraArgsSummaryStatistics = private$extraArgsSummaryStatistics,
                                                              extraArgsProposal = private$extraArgsProposal
                                                              )
                                 self$toc <- Sys.time()

                                 if(private$first)
                                     private$first <- FALSE
                                 ## if(!private$first){
                                 ##     if("posterior" %in% names(est.out)){
                                 ##         ##est.out$posterior<- rbind(private$Posterior, est.out$posterior)
                                 ##         est.out$posterior<- est.out$posterior
                                 ##     } else {
                                 ##         private$Posterior <- rbind(private$Posterior, est.out)
                                 ##     }
                                 ##     if("extra" %in% names(est.out)) {
                                 ##         private$EstimatorOutput <- rbind(private$EstimatorOutput, est.out$extra)
                                 ##     }

                                 ## } else {
                                 ##     private$first <- FALSE
                                 ## }

                                 private$Posterior <- est.out$posterior

                                 if("slamOut" %in% names(est.out))
                                     self$changeExtraArgsEstimator(xbar = est.out$slamOut$xbar,
                                                                   sigma = est.out$slamOut$sigma)

                                 if("extra" %in% names(est.out))
                                     private$EstimatorOutput <- est.out$extra

                             },
                             ## plot the posterior.
                             ## Histograms for all the parameters in the posterior.
                             hist = function(){
                                 size <- dim(private$Posterior)
                                 len <- size[2]-1 ## don't mind the last column.
                                 numRows = ceiling(len/2)
                                 par(mfrow = c(numRows,2))
                                 for(i in 1:len)
                                     hist(private$Posterior[,i], main = "", xlab = names(private$Posterior)[i])
                             },
                             ## Plot the full data frame. Scatterplots of all the
                             plot = function(n = 1, N = NULL){
                                 if(private$first){
                                     print("No posterior estimated yet")
                                     print("Run: runEstimation() first")
                                     return()
                                 }
                                 Posterior <- private$Posterior
                                 dims <- dim(Posterior)
                                 len <- dims[2] - 1

                                 if(is.null(N))
                                     N <- dims[1]

                                 data <- Posterior[n:N,]

                                 ## for abc, ordering can be of importance.
                                 data <- data[order(data[,dims[2]]),]
                                 data <- data[1:N,1:len]


                                 g <- GGally::ggpairs(data = data[n:N,],
                                                      columnLabels = c("upsilon", "beta[1]", "beta[2]", "gamma"),
                                                      labeller = "label_parsed",
                                                      upper = list(continuous = "density"),
                                                      lower = list(continuous = GGally::wrap("points", alpha = 0.6)))
                                 g

                                 ## ## want to show the histograms on the diagonal.
                                 ## panel.hist <- function(x, ...){
                                 ##     usr <- par("usr"); on.exit(par(usr))
                                 ##     par(usr = c(usr[1:2], 0, 1.5) )
                                 ##     h <- hist(x, plot = FALSE)
                                 ##     breaks <- h$breaks; nB <- length(breaks)
                                 ##     y <- h$counts; y <- y/max(y)
                                 ##     rect(breaks[-nB], 0, breaks[-1], y, ...)
                                 ## }
                                 ## pairs(x=data, diag.panel = panel.hist)

                             },
                             ## Summary of the attached functions.
                             summary = function(){
                                 ## Simulator
                                 cat("Simulator:\n")
                                 print(private$Simulator.name)
                                 cat("Simulator arguments:\n")
                                 print(names(private$extraArgsSimulator))
                                 cat("\n")
                                 ## Summary
                                 cat("SummaryStatistics:\n")
                                 print(private$SummaryStatistics.name)
                                 cat("SummaryStatistics arguments:\n")
                                 print(names(private$extraArgsSummaryStatistics))
                                 cat("\n")
                                 ## Estimator
                                 cat("Estimator:\n")
                                 print(private$SummaryStatistics.name)
                                 cat("Estimator arguments:\n")
                                 print(names(private$extraArgsEstimator))
                                 cat("\n")
                                 ## Proposal
                                 cat("Proposal:\n")
                                 print(private$Proposal.name)
                                 ## Tic Toc
                                 cat("Last estimation took the following time:\n")
                                 print(self$toc - self$tic)
                             },
                             ## Execute a simulation for a specific theta
                             runSimulator = function(theta){
                                 return(private$Simulator(theta = theta, extraArgs = private$extraArgsSimulator))
                             },
                             ## From a given simulation, extract summary statistics
                             runStatistics = function(data){
                                 return(private$SummaryStatistics(data = data, extraArgs = private$extraArgsSummaryStatistics))
                             },
                             ## Set the observation
                             setObservation = function(observation){
                                 private$observation <- observation
                             },
                             ## Get the posterior
                             getPosterior = function(){
                                 return(private$Posterior)
                             },
                             ## Get Simulator extraarguments
                             getExtraArgsSimulator = function(...){
                                 return(private$extraArgsSimulator)
                             },
                             ## Get Summary Statistics extraarguments
                             getExtraArgsSummaryStatistics = function(...){
                                 return(private$extraArgsSummaryStatistics)
                             },
                             ## Get Estimator extraarguments
                             getExtraArgsEstimator = function(...){
                                 return(private$extraArgsEstimator)
                             },
                             ## Get Proposal extraarguments
                             getExtraArgsProposal = function(...){
                                 return(private$extraArgsProposal)
                             },
                             ## Get Extra output from Estimator
                             getEstimatorOutput = function(...){
                                 return(private$EstimatorOutput)
                             },
                             ## Get the acceptance rate from Estimator
                             getAccRate = function(){
                                 if(private$first)
                                     return(0)
                                 else {
                                     d <- self$getPosterior()
                                     dims <- dim(d)
                                     rows <- dims[1]
                                     cols <- dims[2]
                                     len <- length(unique(d[,cols]))

                                     accRate <- len / rows

                                     return(accRate)
                                 }
                             },
                             ##--------------------------
                             ## Change private attributes!
                             ##--------------------------
                             ## Change the estimator function.
                             changeEstimator = function(Estimator){
                                 private$Estimator <- Estimator
                             },
                             ## Change the Simulator extraargument.
                             changeExtraArgsSimulator = function(...){
                                 items <- list(...)
                                 for(nam in names(items))
                                     private$extraArgsSimulator[[nam]] <- items[[nam]]
                             },
                             ## Change the Summary Statistics extraargument.
                             changeExtraArgsSummaryStatistics = function(...){
                                 items <- list(...)
                                 for(nam in names(items))
                                     private$extraArgsSummaryStatistics[[nam]] <- items[[nam]]
                             },
                             ## Change the Estimator extraargument.
                             changeExtraArgsEstimator = function(...){
                                 items <- list(...)
                                 for(nam in names(items))
                                     private$extraArgsEstimator[[nam]] <- items[[nam]]
                             },
                             ##--------------------------
                             ## Timing variables
                             ##--------------------------
                             ## Tic
                             tic = NULL,
                             ## Toc
                             toc = NULL
                         ),
                         private = list(
                             ## set the function names
                             setNames = function(Simulator, SummaryStatistics, Estimator, Proposal){
                                 private$Simulator.name <- as.character(match.call()[1])
                                 private$SummaryStatistics.name <- as.character(match.call()[2])
                                 private$Estimator.name <- as.character(match.call()[3])
                                 private$Proposal.name <- as.character(match.call()[4])
                             },
                             ## The simulator
                             Simulator = NULL,
                             Simulator.name = NULL,
                             ## Function to summarize simulated data into a vector of summary statistics
                             SummaryStatistics = NULL,
                             SummaryStatistics.name = NULL,
                             ## The estimator to use, i.e, ABC or SLMCMC.
                             Estimator = NULL,
                             Estimator.name = NULL,
                             ## Paramter proposal function
                             Proposal = NULL,
                             Proposal.name = NULL,
                             ## obeservation data
                             observation = NULL,
                             ## Posterior distribution
                             Posterior = NULL,
                             ## Arguments for the estimator
                             extraArgsEstimator = NULL,
                             ## Extra output that can be given from the estimator.
                             EstimatorOutput = NULL,
                             ## Arguments for the simulator
                             extraArgsSimulator = NULL,
                             ## Argument for the summary statistics transformator.
                             extraArgsSummaryStatistics = NULL,
                             ## Arguments for the proposal function
                             extraArgsProposal = NULL,
                             ## flag to check if the estimator has been runned before.
                             first = TRUE
                         )
                         )

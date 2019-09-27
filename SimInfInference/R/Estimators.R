##' Robin Eriksson 2018
##' Selection of parameter posterior estimation methods.

##' --------------------------------------------------------------
##' Algortihm:
##' While  n < Nstop
##' --- t <- propose parameters
##' --- d <- Simulate data(t)
##' --- s <- transform into summary statistics(d)
##' --- if( || s_org - s || < epsilon )
##' --- --- store t in theta
##' --- --- n <- n + 1
##' --- else
##' --- --- n <- n
##' return theta
##' --------------------------------------------------------------
##' ABC estimator. Will perfome
##' Approximate Bayesian Computations in order
##' to approximate the posterior distribution
##' of the model parameters.
##' @param observation Observed data, that parameters are to be fit to.
##' @param Simulator (function) Model simulator. Theta dependent.
##' @param SummaryStatistics (function) transformer of the output data from the simulator. Output is a vector.
##' @param Proposal proposal functions for the parameters.
##' @param extraArgsABC extra argruments for ABC: epsilon, nStop, debug
##' @param extraArgsSimulator extra arguments for Simulator
##' @param extraArgsSummaryStatistics extra arguments for summary statistics transformer.
##' @param extraArgsProposal extra arguments for the proposal function
##' @return data frame with samples from the posterior distribution and their "ABC" distance.
##' @export
ABC <- function(observation, Simulator, SummaryStatistics, Proposal,
                extraArgsEstimator, extraArgsSimulator, extraArgsSummaryStatistics,
                extraArgsProposal){
    stopifnot(c("epsilon", "nStop", "debug", "distance", "theta0", "accVec") %in% names(extraArgsEstimator))

    epsilon <- extraArgsEstimator$epsilon
    nStop <- extraArgsEstimator$nStop
    debug <- extraArgsEstimator$debug
    distance <- extraArgsEstimator$distance
    theta0 <- extraArgsEstimator$theta0
    accVec <- extraArgsEstimator$accVec

    ## create summary statistics from observations.
    ss_obs <- SummaryStatistics(data = observation, extraArgs = extraArgsSummaryStatistics)

    ## output. The set of parameters drawn from the approximate posterior.
    ## dimensions are govern by the number of parameters.
    if(is.null(accVec)){
        theta <- matrix(0, nrow = nStop, ncol = length(theta0))
        delta <- numeric(nStop)

        nStart <- 1
    }
    else{
        dims <- dim(accVec)
        theta <- matrix(0, nrow = nStop + dims[1], ncol = dims[2]-1)
        theta[1:dims[1],] <- accVec[,1:(dims[2]-1)]
        delta <- numeric(nStop + dims[1])
        delta[1:dims[1]] <- accVec[,dims[2]]

        nStart <- dims[1] + 1
        nStop <- dims[1] + nStop
    }

    colnames(theta) <- names(theta0)


    ## progressbar
    pb <- pbapply::timerProgressBar(width = 50)

    for(n in nStart:nStop){
        ## Propose parameters
        theta_prop <- Proposal(theta = theta[n,], extraArgs = extraArgsProposal)

        if(debug) {
                cat("\n***************\ndebug info\n***************\n")
                cat(theta_prop, "\n")
        }

        dist <- computeDelta(Simulator, SummaryStatistics,
                             extraArgsSimulator, extraArgsSummaryStatistics,
                             distance, theta_prop, ss_obs)

        theta[n,] <- theta_prop
        delta[n] <- dist

        if(debug) {
            cat("\n")
            cat("dist = ", dist, "\n")
        }


        ## }
        ## progressbar
        pbapply::setTimerProgressBar(pb, n/nStop)
    }
    close(pb)

    if(!is.null(epsilon))
        accRows <- which(delta < epsilon)
    else
        accRows <- 1:nStop
    theta.m <- theta[accRows,]
    delta.m <- delta[accRows]


    ## organize output
    ##theta.m <- matrix(unlist(theta), nrow = length(theta_prop))

    ## organize output
    theta.d <- data.frame(theta, delta)
    names(theta.d) <- c(names(theta0), "delta")
    return(list(posterior = theta.d))
}

##' Simulate data and return the synthetic log likelihood
##' @param Simulator the simulator
##' @param SummaryStatistics the summary statistics computer
##' @param extraArgsSimulator the arguments for the simulator
##' @param extraArgsSummaryStatistics the arguments for the statistics
##' @param theta the model parameter
##' @param ss_obs the observed summary statistics
##' @return the synthetic likelihood of the proposed parameter theta
##' @export
computeDelta <- function(Simulator, SummaryStatistics,
                         extraArgsSimulator,extraArgsSummaryStatistics,
                         distance, theta, ss_obs) {

    ## Simulate data using the proposed parameters
    data <- try(
        Simulator(theta = theta, extraArgs = extraArgsSimulator)
    )

    ## Calculate the summary statistics.
    if(is(data, "try-error")){
        ss.m <- t(matrix(NA, nrow = 1, ncol = length(ss_obs)))
    } else {
        ss.m <- t(matrix(SummaryStatistics(data = data,
                                           extraArgs = extraArgsSummaryStatistics),
                         nrow = 1))
    }

    ## Check distance between SS
    delta <- try(
        distance(ss_obs, ss.m)
    )

    if(is(dist, "try-error") || !is.numeric(delta) || !is.finite(delta))
        delta <- Inf

    return(delta)
}



##' --------------------------------------------------------------
##' SLAM (Synthetic likelihood Adaptive Metropolis  estimator.
##' Computes the synthetic likelihood (SL) (Wood 2010) and the performes a Metropolis-Hasting algortihm
##' with SL instead of the true likelihood
##' @param observation Observed data, that parameters are to be fit to.
##' @param Simulator (function) Model simulator. Theta dependent.
##' @param SummaryStatistics (function) transformer of the output data from the simulator. Output is a vector.
##' @param Proposal proposal functions for the parameters.
##' @param extraArgsEstimator extra argruments for the estimator: nStop, nSim, debug
##' @param extraArgsSimulator extra arguments for Simulator
##' @param extraArgsSummaryStatistics extra arguments for summary statistics transformer.
##' @param extraArgsProposal extra arguments for the proposal function
##' @return samples from the posterior distribution.
##' @export
SLAM <- function(observation, Simulator, SummaryStatistics, Proposal,
                 extraArgsEstimator, extraArgsSimulator,
                 extraArgsSummaryStatistics, extraArgsProposal){
    stopifnot(c("nStop", "nSim", "debug", "theta0", "parameters",
                "normalize", "C0", "e", "S", "reinit", "n0", "logParam") %in% names(extraArgsEstimator))
    reinit <- extraArgsEstimator$reinit


    nStop <- extraArgsEstimator$nStop
    nSim <- extraArgsEstimator$nSim
    debug <- extraArgsEstimator$debug
    theta0 <- extraArgsEstimator$theta0
    parameters <- extraArgsEstimator$parameters
    normalize <- extraArgsEstimator$normalize
    C0 <- extraArgsEstimator$C0
    e <- extraArgsEstimator$e
    S <- extraArgsEstimator$S

    n0 <- extraArgsEstimator$n0
    logParam <- extraArgsEstimator$logParam

    ## create summary statistics from observations.
    ## Do this ones here, in order to save potentially "stly" computations.
    ss_obs <- SummaryStatistics(data = observation, extraArgs = extraArgsSummaryStatistics)


    if(reinit) {
        reVals <- SLAMreInit(observation, Simulator, SummaryStatistics, Proposal,
                              extraArgsEstimator, extraArgsSimulator,
                              extraArgsSummaryStatistics, extraArgsProposal)
        d <- reVals$d
        theta <- reVals$theta
        syntl <- reVals$syntl
        numAcc <- length(unique(syntl[1:d]))
        xbar <- extraArgsEstimator$xbar
        sigma <- extraArgsEstimator$sigma
        sl_old <- syntl[d]
    } else {
        ## output. The set of parameters drawn from the approximate posterior.
        ## dimensions are govern by the number of parameters.

        ## ------------------------------------ ##
        ## --- Set-up for new and re-init ----- ##
        ## ------------------------------------ ##

        theta <- matrix(numeric(length(theta0)*nStop), ncol = length(theta0))
        syntl <- numeric(nStop)
        theta[1,] <- theta0


        colnames(theta) <- names(theta0)


        ## init synthetic likelihood
        if(debug){
            print("Preparing starting s-likelihood")
            if(logParam)
                print(exp(theta[1,]))
            else
                print(theta[1,])
        }


        ## ---------------------------- ##
        ## Compute synthetic likelihood ##
        ## ---------------------------- ##
        sl_old <- computeSL(Simulator, SummaryStatistics,
                            extraArgsSimulator,
                            extraArgsSummaryStatistics,
                            extraArgsEstimator,
                            theta0, ss_obs)



        syntl[1] <- sl_old
        d <- 1

        ## report acceptance rate
        numAcc <- 1
    }

    if(debug)
            print(sprintf("sl = %f", sl_old))





    if(d > 1)
        lenStop <- d+nStop-1
    else
        lenStop <- nStop-1

    if(d > 1) {
        if(!(d %% 100))
            init <- TRUE
        else
            init <- FALSE
    } else
        init <- TRUE


    ## progressbar
    pb <- pbapply::timerProgressBar(width = 50)

    ## ----------------------------- ##
    ## --- Main SLAM algorithm ----- ##
    ## ----------------------------- ##
    for (n in d:lenStop){
        ## ---------------------------- ##
        ## --- Propose parameters ----- ##
        ## ---------------------------- ##

        ## Estimate covariance
        if(n <= n0){
            sigma <- diag(rep(C0,length(theta[n,])))
        } else if(init) {
            C <- cov.est(theta[1:n,])
            sigma <- S * C + S*e*diag(length(theta[n,]))
            xbar <- apply(theta[1:n,],2,mean)
            init <- FALSE
        } else {
            C.rec <- cov.rec(theta[n,], xbar = xbar, C = sigma, k = n, s = S, e = e)
            sigma <- C.rec$C
            xbar <- C.rec$xbar

            ## to make sure that the recursive covariance is in line with
            ## the correct value, re-init the covariance every 100th value.
            if(!(n %% 100))
                init <- TRUE
        }


        ## stepping
        eps <- MASS::mvrnorm(1, mu = numeric(length(theta[n,])), Sigma = sigma)
        theta_prop <- theta[n,] + eps

        ## ---------------------------- ##
        ## Compute synthetic likelihood ##
        ## ---------------------------- ##
        sl_new <- computeSL(Simulator, SummaryStatistics,
                            extraArgsSimulator,
                            extraArgsSummaryStatistics,
                            extraArgsEstimator,
                            theta_prop, ss_obs)


        ## ---------------------------- ##
        ## ---acceptance probability--- ##
        ## ---------------------------- ##
        alpha <- min(1, exp(sl_new - sl_old))

        u <- runif(n = 1)

        if(debug){
            cat("\n***************\ndebug info\n***************\n")
            if(logParam)
                cat(exp(theta_prop), "\n")
            else
                cat(theta_prop, "\n")
            cat("\n")
            cat("sl = ", sl_new, "\n")
            cat("alpha = ", alpha, "\n")
            cat("u = ", u, "\n")
        }

        if(u <= alpha){
            if(debug){
                cat("hit\n")
            }
            numAcc <- numAcc + 1
            theta[n+1,] <- theta_prop
            syntl[n+1] <- sl_new
            sl_old <- sl_new

        } else {
            theta[n+1,] <- theta[n,]
            syntl[n+1] <- syntl[n]
        }

        ## update progressbar
        if(d > 1)
            frac <- (n-(d-1))/(lenStop-(d-1))
        else
            frac <- n/lenStop
        pbapply::setTimerProgressBar(pb, frac)

    }
    close(pb)

    cat("\nAcceptance rate: ", numAcc/(lenStop+1), "\n")

    ## organize output
    if(logParam)
        thetaOut <- exp(theta)
    else
        thetaOut <- theta

    ## organize output
    theta.d <- data.frame(theta, syntl)
    names(theta.d) <- c(parameters, "Synthetic likelihood")
    return(list(posterior = theta.d, slamOut = list(xbar = xbar, sigma = sigma)))
}

##' helper function for SLAM when reinit the Markov chain.
SLAMreInit <- function(observation, Simulator, SummaryStatistics, Proposal,
                       extraArgsEstimator, extraArgsSimulator,
                       extraArgsSummaryStatistics, extraArgsProposal){
    stopifnot(c("nStop", "nSim", "debug", "theta0", "parameters",
                "normalize", "C0", "e", "S", "accVec", "n0", "xbar",
                "sigma", "logParam") %in% names(extraArgsEstimator))

    nStop <- extraArgsEstimator$nStop
    nSim <- extraArgsEstimator$nSim
    debug <- extraArgsEstimator$debug
    theta0 <- extraArgsEstimator$theta0
    parameters <- extraArgsEstimator$parameters
    normalize <- extraArgsEstimator$normalize
    C0 <- extraArgsEstimator$C0
    e <- extraArgsEstimator$e
    S <- extraArgsEstimator$S
    accVec <- extraArgsEstimator$accVec
    n0 <- extraArgsEstimator$n0
    logParam <-extraArgsEstimator$logParam

    ## extract values from the previous states in the chain.
    xbar <- extraArgsEstimator$xbar
    sigma <- extraArgsEstimator$sigma

    ## create summary statistics from observations.
    ## Do this ones here, in order to save potentially "costly" computations.
    ss_obs <- SummaryStatistics(data = observation, extraArgs = extraArgsSummaryStatistics)

    ## output. The set of parameters drawn from the approximate posterior.
    ## dimensions are govern by the number of parameters.

    ## ------------------------------------ ##
    ## --- Set-up for re-init ----- ##
    ## ------------------------------------ ##

    d <- dim(accVec)[1]
    theta <- matrix(numeric(length(theta0)*(d+nStop)), ncol = length(theta0))
    theta[1:d,1:length(theta0)] <- accVec[,1:length(theta0)] ## last element is SL
    syntl <- numeric(d+nStop)
    syntl[1:d] <- accVec[,length(theta0)+1]


    colnames(theta) <- names(theta0)

    sl_old <- syntl[d]

    lenStop <- d+nStop-1

    return(list(d = d, theta = theta, syntl = syntl))
}

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

mean.rec <- function(xbar, x, t){
    xbarnew <- 1/t*((t-1) * xbar + x)
}

cov.rec <- function(X,xbarOld,C, k, s, e){

    xbar <- mean.rec(xbarOld, X, k)

    C <-  (k-1)/k * C + s/k * (k * xbarOld %*% t(xbarOld) -
                               (k+1) * xbar %*% t(xbar) +
                               X %*% t(X) +
                               e*diag(length(X)))

    return(list(C = C, xbar = xbar))
}

empCov.wood <- function(X){
    sigma <- tryCatch(synlik::robCov(t(X))$COV,
                      error = function(e){
                          ## compute mu_hat (mean of SS)
                          muHat <- apply(X,2,mean)

                          ## compute S (s - mu_hat)
                          S <- apply(X,1,function(X,y){X-y},y=muHat)

                          ## compute sigma_hat
                          sigmaHat <- S %*% t(S) / (dim(S)[1]-1)
                          return(sigmaHat)
                      })
    return(sigma)
}

##' Simulate data and return the synthetic log likelihood
##' @param Simulator the simulator
##' @param SummaryStatistics the summary statistics computer
##' @param extraArgsSimulator the arguments for the simulator
##' @param extraArgsSummaryStatistics the arguments for the statistics
##' @param extraArgsEstimator the arguments for the estimator
##' @param theta the model parameter
##' @param nSim the number of simulations
##' @param ss_obs the observed summary statistics
##' @param normalize if the summary statistics are to be normalized or not.
##' @return the synthetic likelihood of the proposed parameter theta
##' @export
computeSL <- function(Simulator, SummaryStatistics,
                      extraArgsSimulator, extraArgsSummaryStatistics,
                      extraArgsEstimator,
                      theta, ss_obs){
    if("multiSS" %in% names(extraArgsEstimator)) {
        multiSS <- extraArgsEstimator$multiSS
        nSS <- multiSS
    } else
        nSS <- 1

    normalize <- extraArgsEstimator$normalize
    nSim <- extraArgsEstimator$nSim
    logParam <- extraArgsEstimator$logParam
    ## this might be removed
    bs <- extraArgsSummaryStatistics$bs

    if(logParam)
        thetaProp <- exp(theta)
    else
        thetaProp <- theta

    if(bs)
        ncol <- extraArgsSummaryStatistics$B
    else
        ncol <- nSim

    ## Calculate the summary statistics.
    sl <- numeric(nSS)
    for(n in seq_len(nSS)) {

        ## Simulate data
        data <- try(
            Simulator(theta = thetaProp, extraArgs = extraArgsSimulator)
        )
        if(is(data, "try-error"))
            return(-Inf)


        ##ss.m <- t(matrix(0, ncol = ncol, nrow = length(ss_obs)))
        ss <- SummaryStatistics(data = data,
                                extraArgs = extraArgsSummaryStatistics)




        if(class(ss) == "matrix"){
            if(dim(ss)[1] == ncol)
                ss.m <- ss
            else
                ss.m <- t(ss)
        } else {
            ss.m <- t(matrix(ss, ncol = ncol, nrow = length(ss_obs)))
        }



        ## if we normalize or not.
        if(normalize){
            ss.m <- sweep(ss.m, 2, ss_obs, "/")
            y <- rep(1, length(ss_obs))
        }

        ## Compute synthetic likelihood
        sl[n] <- try(
            SyntheticLogLikelihood_wood(y = y, X = ss.m)
        )

        if(is(sl[n], "try-error") || !is.numeric(sl[n]) || !is.finite(sl[n]))
            sl <- -Inf

        ## make sure that sl_new stays numeric and not factor (a bug)
        if(class(sl[n]) == "factor")
            sl[n] <- as.numeric(levels(sl[n])[sl[n]])

        ## cat("\n")
        ## print(round(as.numeric(apply(ss,1,mean)),4))
        ## cat("\n")
        ## print(sl[n])
    }

    return(sl)
}

##' Compute the synthetic likelihood
##' @param ss_obs the observed synthetic likelihood
##' @param ss simulated summary statistics as a matrix, each row is a simulaton.
##' @export
SyntheticLogLikelihood <- function(ss_obs, ss){
    ## compute mu_hat (mean of SS)
    muHat <- apply(ss,2,mean)

    ## compute S (s - mu_hat)
    S <- apply(ss,1,function(x,y){x-y},y=muHat)

    ## compute sigma_hat

    sigmaHat <- S %*% t(S) / (dim(S)[1]-1)
    invSigma <- MASS::ginv(sigmaHat)
    logDet <- 0.5*log(det(sigmaHat))
    ## Robust estimation of sigma hat
    ## rc <- robCov(S)
    ## E <- rc$E
    ## invSigma <- t(E) %*% E
    ## logDet <- rc$half.ldet.V
    ##sigmaHat <- rc$COV

    ## l for each element
    vec <- array(ss_obs-muHat)
    ##lsOld <- -0.5*t(vec) %*% ginv(sigmaHat) %*% vec - 0.5*log(det(sigmaHat))
    ls <- -0.5*t(vec) %*% invSigma %*% vec - logDet

    ## if(is.infinite(ls) || is.nan(ls))
    ##     return(-1e10)
    ## browser()
    ##L <- sum(ls)
    ls

}

##' SL grid exploration
##'
##' Compute the SL on a "given" grid.
##' @param observation Observed data, that parameters are to be fit to.
##' @param Simulator (function) Model simulator. Theta dependent.
##' @param SummaryStatistics (function) transformer of the output data from the simulator. Output is a vector.
##' @param extraArgsEstimator extra argruments for the estimator: nStop, nSim, debug
##' @param extraArgsSimulator extra arguments for Simulator
##' @param extraArgsSummaryStatistics extra arguments for summary statistics transformer.
##' @param proposal null
##' @return grid parameters with distance
##' @export
slGrid <- function(observation, Simulator, SummaryStatistics, Proposal,
                   extraArgsEstimator, extraArgsSimulator,
                   extraArgsSummaryStatistics, extraArgsProposal){
    stopifnot(c("theta0", "pert", "thetalength", "normalize", "nSim", "logParam", "allDim", "multiSS")
              %in% names(extraArgsEstimator))


    theta0 <- extraArgsEstimator$theta0
    pert <- extraArgsEstimator$pert
    thetalength <- extraArgsEstimator$thetalength
    allDim <- extraArgsEstimator$allDim
    multiSS <- extraArgsEstimator$multiSS

    ## Observation!
    ss_obs <- SummaryStatistics(data = observation, extraArgs = extraArgsSummaryStatistics)

    ## create grid of all possible combinations.
    thetaVec <- genGrid(theta0 = theta0, pert = pert,
                        thetalength = thetalength, allDim = allDim)



    len <- dim(thetaVec)[1]
    ## Computed synthetic likelihood vector .
    sl_vec <- matrix(NA, nrow = len, ncol = multiSS)

    theta_prop <- theta0
    pb <- pbapply::timerProgressBar(width = 50)
    for(i in 1:len){
        for(j in seq_len(length(theta0)))
            theta_prop[j] <- thetaVec[i,j]
        ## theta_prop <- c(upsilon = thetaVec$upsilon[i],
        ##                 beta_t1 = thetaVec$beta_t1[i],
        ##                 beta_t2 = thetaVec$beta_t2[i],
        ##                 gamma = thetaVec$gamma[i])

        sl_new <- computeSL(Simulator, SummaryStatistics,
                            extraArgsSimulator,
                            extraArgsSummaryStatistics,
                            extraArgsEstimator,
                            theta_prop, ss_obs)

        sl_vec[i,] <- sl_new
        pbapply::setTimerProgressBar(pb, i/len)
    }
    close(pb)

    theta.d <- data.frame(thetaVec, dist = sl_vec)
    return(list(posterior = theta.d))
}

##' generate parameter grid to explore
##' @param theta0 the center point
##' @param thetalength how fine the grid should be
##' @param pert by how much should the center be pertubated?
##' @param allDim 1d or all d?
##' @export
genGrid <- function(theta0, thetalength, pert, allDim){

    ## create grid of all possible combinations.
    theta.list <- list()
    for (i in seq_len(length(theta0))) {
        theta.i <- seq(from = pert[1]*theta0[i],
                       to = pert[2]*theta0[i],
                       length.out = thetalength)
        theta.list[[length(theta.list)+1]] <- theta.i
    }

    ## create a parameter grid to traverse.
    if(allDim) {

        thetaVec <- expand.grid(theta.list)
        colnames(thetaVec) <- names(theta0)
    } else {
        thetaVec <- as.data.frame(matrix(0, nrow = thetalength*length(theta0),
                                         ncol = length(theta0)))
        names(thetaVec) <- names(theta0)

        for (i in seq_len(length(theta0))) {
            prebox <- rep(theta0[i], thetalength*(i-1))
            postbox <- rep(theta0[i], thetalength*(length(theta0)-i))

            thetaVec[,i] <- c(prebox, theta.list[[i]], postbox)
        }
    }
    return(thetaVec)
}

##' --------------------------------------------------------------
##' adaptiveMIS (Synthetic likelihood Adaptive Metropolis independent sampler)
##' Computes the synthetic likelihood (SL) (Wood 2010) and the performes a (adaptive ) MIS
##' with SL instead of the true likelihood
##' @param observation Observed data, that parameters are to be fit to.
##' @param Simulator (function) Model simulator. Theta dependent.
##' @param SummaryStatistics (function) transformer of the output data from the simulator. Output is a vector.
##' @param Proposal proposal functions for the parameters.
##' @param extraArgsEstimator extra argruments for the estimator: nStop, nSim, debug
##' @param extraArgsSimulator extra arguments for Simulator
##' @param extraArgsSummaryStatistics extra arguments for summary statistics transformer.
##' @param extraArgsProposal extra arguments for the proposal function
##' @return samples from the posterior distribution.
##' @export
adaptiveMIS <- function(observation, Simulator, SummaryStatistics, Proposal,
                 extraArgsEstimator, extraArgsSimulator,
                 extraArgsSummaryStatistics, extraArgsProposal){
    stopifnot(c("nStop", "nSim", "debug", "theta0", "parameters",
                "normalize", "C0", "e", "S", "reinit", "n0", "logParam") %in% names(extraArgsEstimator))
    reinit <- extraArgsEstimator$reinit


    nStop <- extraArgsEstimator$nStop
    nSim <- extraArgsEstimator$nSim
    debug <- extraArgsEstimator$debug
    theta0 <- extraArgsEstimator$theta0
    parameters <- extraArgsEstimator$parameters
    normalize <- extraArgsEstimator$normalize
    C0 <- extraArgsEstimator$C0
    e <- extraArgsEstimator$e
    S <- extraArgsEstimator$S

    n0 <- extraArgsEstimator$n0
    logParam <- extraArgsEstimator$logParam

    ## create summary statistics from observations.
    ## Do this ones here, in order to save potentially "stly" computations.
    ss_obs <- SummaryStatistics(data = observation, extraArgs = extraArgsSummaryStatistics)



    if(reinit) {
        reVals <- SLAMreInit(observation, Simulator, SummaryStatistics, Proposal,
                              extraArgsEstimator, extraArgsSimulator,
                              extraArgsSummaryStatistics, extraArgsProposal)
        d <- reVals$d
        theta <- reVals$theta
        syntl <- reVals$syntl
        numAcc <- length(unique(syntl[1:d]))
        xbar <- extraArgsEstimator$xbar
        sigma <- extraArgsEstimator$sigma
        sl_old <- syntl[d]
    } else {
        ## output. The set of parameters drawn from the approximate posterior.
        ## dimensions are govern by the number of parameters.

        ## ------------------------------------ ##
        ## --- Set-up for new and re-init ----- ##
        ## ------------------------------------ ##

        theta <- matrix(numeric(length(theta0)*nStop), ncol = length(theta0))
        syntl <- numeric(nStop)
        theta[1,] <- theta0


        colnames(theta) <- names(theta0)


        ## init synthetic likelihood
        if(debug){
            print("Preparing starting s-likelihood")
            if(logParam)
                print(exp(theta[1,]))
            else
                print(theta[1,])
        }


        ## ---------------------------- ##
        ## Compute synthetic likelihood ##
        ## ---------------------------- ##
        sl_old <- computeSL(Simulator, SummaryStatistics,
                            extraArgsSimulator,
                            extraArgsSummaryStatistics,
                            extraArgsEstimator,
                            theta0, ss_obs)



        syntl[1] <- sl_old
        d <- 1

        sigma <- NULL
        xbar <- NULL
        ## report acceptance rate
        numAcc <- 1
    }

    ## the weight of the previous proposal
    if(is.null(xbar) | is.null(sigma))
        w_old <- exp( sl_old - emdbook::dmvnorm(1, mu = theta0, Sigma = diag(0.1,dimlen,dimlen), log = TRUE))
    else
        w_old <- exp( sl_old - emdbook::dmvnorm(theta0, mu = xbar, Sigma = sigma, log = TRUE))

    if(debug)
        print(w_old)

    if(d > 1)
        lenStop <- d+nStop-1
    else
        lenStop <- nStop-1




    ## progressbar
    pb <- pbapply::timerProgressBar(width = 50)

    ## ----------------------------- ##
    ## --- Main MIS algorithm ----- ##
    ## ----------------------------- ##
    dimlen <- length(theta[,1])
    for (n in d:lenStop){
        ## ---------------------------- ##
        ## --- Propose parameters ----- ##
        ## ---------------------------- ##

        ## do MIS from some prior.
        if(is.null(xbar) | is.null(sigma)) {
            theta_prop <- MASS::mvrnorm(1, mu = theta[1,], Sigma = diag(0.1,dimlen,dimlen))
        } else {
            ##  do adaptive MIS
            theta_prop <- MASS::mvrnorm(1, mu = xbar, Sigma = sigma)
        }

        ## ---------------------------- ##
        ## Compute synthetic likelihood ##
        ## ---------------------------- ##

        sl_new <- computeSL(Simulator, SummaryStatistics,
                            extraArgsSimulator,
                            extraArgsSummaryStatistics,
                            extraArgsEstimator,
                            theta_prop, ss_obs)

        ## compute importance weight for the sample
        # w <- likelihood / prior
        w_new <- exp( sl_new - emdbook::dmvnorm(theta_prop, mu = xbar, Sigma = sigma, log = TRUE))

        ## ---------------------------- ##
        ## ---acceptance probability--- ##
        ## ---------------------------- ##
        alpha <- min(1, w_new / w_old)

        u <- runif(n = 1)

        if(debug){
            cat("\n***************\ndebug info\n***************\n")
            if(logParam)
                cat(exp(theta_prop), "\n")
            else
                cat(theta_prop, "\n")
            cat("\n")
            cat("sl = ", sl_new, "\n")
            cat("w = ", w_new, "\n")
            cat("alpha = ", alpha, "\n")
            cat("u = ", u, "\n")
        }

        if(u <= alpha){
            if(debug){
                cat("hit\n")
            }
            numAcc <- numAcc + 1
            theta[n+1,] <- theta_prop
            syntl[n+1] <- sl_new
            sl_old <- sl_new
            w_old <- w_new
        } else {
            theta[n+1,] <- theta[n,]
            syntl[n+1] <- syntl[n]
        }

        ## update progressbar
        if(d > 1)
            frac <- (n-(d-1))/(lenStop-(d-1))
        else
            frac <- n/lenStop
        pbapply::setTimerProgressBar(pb, frac)

    }
    close(pb)

    cat("\nAcceptance rate: ", numAcc/(lenStop+1), "\n")

    ## organize output
    if(logParam)
        thetaOut <- exp(theta)
    else
        thetaOut <- theta

    ## organize output
    theta.d <- data.frame(theta, syntl)
    names(theta.d) <- c(parameters, "Synthetic likelihood")
    return(list(posterior = theta.d, slamOut = list(xbar = xbar, sigma = sigma)))



}

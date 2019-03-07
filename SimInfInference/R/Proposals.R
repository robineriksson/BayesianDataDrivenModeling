
##' Proposal functions for all parameters in the model
##' @return names list with names lists for each parameter.
##' @export
Proposal_SISe3sp <- function(){
    ## we say that 2*up1 = 2*up2 = up3
    upsilon <- function(){runif(n = 1, min = 1e-4, max = 1e-1)}
    beta_t1 <- function(){runif(n = 1, min = 1e-4, max = 5e-1)}
    beta_t2<-  function(){runif(n = 1, min = 1e-4, max = 5e-1)}
    beta_t3<-  function(){runif(n = 1, min = 1e-4, max = 5e-1)}
    beta_t4<-  function(){runif(n = 1, min = 1e-4, max = 5e-1)}

    return(list(upsilon = upsilon, beta_t1 = beta_t1, beta_t2 = beta_t2, beta_t3 = beta_t3, beta_t4 = beta_t4))
}

##' Proposal functions for all parameters in the SISe_sp model
##'
##' Propose a new parameter set. Independent of the lastly proposed
##' and in all dimensions at once.
##' @param theta the old parameter vector
##' @param extraArgs not in use.
##' @return proposed parameter vector
##' @export
Proposal_SISe_allUniformNarrow <- function(theta, extraArgs){
    all <- function(theta){
        theta["upsilon"] <- runif(n = 1, min = 0, max = 0.015)
        theta["beta_t1"] <- runif(n = 1, min = 0, max = 0.1)
        theta["beta_t2"] <- runif(n = 1, min = 0, max = 0.17)
        theta["gamma"] <- runif(n = 1, min = 0, max = 0.2)
        return(list(theta = theta, id = 1))
    }

    probVec <- rep(1,1)
    probVec <- probVec / sum(probVec)

    ## just a check to make sure!
    stopifnot( sum(probVec) == 1)

    funVec <- c(all)
    fun <- sample(x = funVec, size = 1, prob = probVec)

    fun.out <- fun[[1]](theta = theta)
    return(fun.out$theta)
}


##' Proposal functions for all parameters in the SISe_sp model
##'
##' Propose a new parameter set. Independent of the lastly proposed
##' and in all dimensions at once.
##' @param theta the old parameter vector
##' @param extraArgs not in use.
##' @return proposed parameter vector
##' @export
Proposal_SISe_allUniformNarrowBeta3 <- function(theta, extraArgs){
    all <- function(theta){
        theta["upsilon"] <- runif(n = 1, min = 0, max = 0.0252)
        theta["beta_t1"] <- runif(n = 1, min = 0, max = 0.21)
        theta["beta_t2"] <- runif(n = 1, min = 0, max = 0.1926)
        theta["beta_t3"] <- runif(n = 1, min = 0, max = 0.2)
        theta["gamma"] <- runif(n = 1, min = 0, max = 0.2)
        return(list(theta = theta, id = 1))
    }

    probVec <- rep(1,1)
    probVec <- probVec / sum(probVec)

    ## just a check to make sure!
    stopifnot( sum(probVec) == 1)

    funVec <- c(all)
    fun <- sample(x = funVec, size = 1, prob = probVec)

    fun.out <- fun[[1]](theta = theta)
    return(fun.out$theta)
}

##' Proposal functions for all parameters in the SISe_sp model
##'
##' Propose a new parameter set. Independent of the lastly proposed
##' and in all dimensions at once.
##' @param theta the old parameter vector
##' @param extraArgs not in use.
##' @return proposed parameter vector
##' @export
Proposal_SISe_allUniformNarrowBeta3Prev <- function(theta, extraArgs){
    all <- function(theta){
        theta["upsilon"] <- runif(n = 1, min = 0, max = 0.0252)
        theta["beta_t1"] <- runif(n = 1, min = 0, max = 0.21)
        theta["beta_t2"] <- runif(n = 1, min = 0, max = 0.1926)
        theta["beta_t3"] <- runif(n = 1, min = 0, max = 0.2)
        theta["gamma"] <- runif(n = 1, min = 0, max = 0.2)
        theta["prev"] <- runif(n = 1, min = 0, max = 0.2)
        return(list(theta = theta, id = 1))
    }

    probVec <- rep(1,1)
    probVec <- probVec / sum(probVec)

    ## just a check to make sure!
    stopifnot( sum(probVec) == 1)

    funVec <- c(all)
    fun <- sample(x = funVec, size = 1, prob = probVec)

    fun.out <- fun[[1]](theta = theta)
    return(fun.out$theta)
}

##' Proposal functions for all parameters in the SISe_sp model
##'
##' Propose a new parameter set. Independent of the lastly proposed
##' and in all dimensions at once.
##' @param theta the old parameter vector
##' @param extraArgs not in use.
##' @return proposed parameter vector
##' @export
Proposal_SISe_allUniformWide <- function(theta, extraArgs){

     all <- function(theta){
        theta["upsilon"] <- runif(n = 1, min = 0, max = 0.25)
        theta["beta_t1"] <- runif(n = 1, min = 0, max = 0.25)
        theta["beta_t2"] <- runif(n = 1, min = 0, max = 0.25)
        theta["gamma"] <- runif(n = 1, min = 0, max = 0.25)
        return(list(theta = theta, id = 1))
    }

    probVec <- rep(1,1)
    probVec <- probVec / sum(probVec)

    ## just a check to make sure!
    stopifnot( sum(probVec) == 1)

    funVec <- c(all)
    fun <- sample(x = funVec, size = 1, prob = probVec)

    fun.out <- fun[[1]](theta = theta)
    return(fun.out$theta)
}

##' Proposal functions for all parameters in the SISe_sp model
##' @return names list with names lists for each parameter.
##' @export
Proposal_SISe_real<- function(){
    upsilon <- function(){runif(n = 1, min = 0.007, max = 0.017)}
    beta_t1 <- function(){runif(n = 1, min = 0.05, max = 0.15)}
    beta_t2 <- function(){runif(n = 1, min = 0.05, max = 0.15)}
    gamma   <- function(){runif(n = 1, min = 0.05, max = 0.15)}

    return(list(upsilon = upsilon, beta_t1 = beta_t1, beta_t2 = beta_t2, gamma = gamma))
}

##' Proposal functions for all parameters in the model
##' @return names list with names lists for each parameter.
##' @export
Proposal_SISe3sp_fix<- function(){
    ## we say that 2*up1 = 2*up2 = up3
    upsilon <- function(){runif(n = 1, min = 1e-4, max = 1e-1)}
    beta_t1 <- function(){0.095}
    beta_t2<-  function(){0.012}
    beta_t3<-  function(){0.1}
    beta_t4<-  function(){0.15}
    return(list(upsilon = upsilon, beta_t1 = beta_t1, beta_t2 = beta_t2,
                beta_t3 = beta_t3, beta_t4 = beta_t4))

}

##' Proposal functions for all parameters in the model
##' Kernel function. Takes a value and alters it with a log-normal random walk.
##' @return names list with names lists for each parameter.
##' @export
Proposal_SISe3sp_Kernel <- function(){
    S <- 0.3
    upsilon <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["upsilon"] <- theta["upsilon"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    beta_t1 <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t1"] <- theta["beta_t1"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    beta_t2 <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t2"] <- theta["beta_t2"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    beta_t3 <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t3"] <- theta["beta_t3"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    beta_t4 <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t4"] <- theta["beta_t4"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    all <- function(theta){
        fac <- exp(rnorm(5)*S)
        qratio <- prod(fac)
        theta <- theta*fac
        return(list(theta = theta, qratio = qratio))
    }
    return(list(upsilon = upsilon, beta_t1 = beta_t1, beta_t2 = beta_t2, beta_t3 = beta_t3, beta_t4 = beta_t4))
}

##' Proposal functions for all parameters in the model
##' Kernel function. Takes a value and alters it with a log-normal random walk.
##' @return names list with names lists for each parameter.
##' @export
Proposal_SISe_Kernel <- function(S = 0.3){
    upsilon <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["upsilon"] <- theta["upsilon"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    beta_t1 <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t1"] <- theta["beta_t1"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    beta_t2 <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t2"] <- theta["beta_t2"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    gamma <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["gamma"] <- theta["gamma"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    all <- function(theta){
        fac <- exp(rnorm(5)*S)
        qratio <- prod(fac)
        theta <- theta*fac
        return(list(theta = theta, qratio = qratio))
    }
    return(list(upsilon = upsilon, beta_t1 = beta_t1, beta_t2 = beta_t2, gamma = gamma, all = all))
}

##' Proposal functions for all parameters in the model
##' Kernel function. Takes a value and alters it with a log-normal random walk.
##' @return names list with names lists for each parameter.
##' @export
Proposal_SISe_Kernel_extended <- function(S = 0.3){
    upsilon <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["upsilon"] <- theta["upsilon"]*qratio
        return(list(theta = theta, qratio = qratio, id = 1))
    }
    beta_t1 <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t1"] <- theta["beta_t1"]*qratio
        return(list(theta = theta, qratio = qratio, id = 2))
    }
    beta_t2 <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t2"] <- theta["beta_t2"]*qratio
        return(list(theta = theta, qratio = qratio, id = 3))
    }
    gamma <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["gamma"] <- theta["gamma"]*qratio
        return(list(theta = theta, qratio = qratio, id = 4))
    }
    upsilonBeta <- function(theta){
        fac <- exp(rnorm(3)*S)
        qratio <- prod(fac)
        theta[c("upsilon","beta_t1", "beta_t2")] <- theta[c("upsilon","beta_t1", "beta_t2")]*fac
        return(list(theta = theta, qratio = qratio, id = 5))
    }
    upsilonGamma <- function(theta){
        fac <- exp(rnorm(2)*S)
        qratio <- prod(fac)
        theta[c("upsilon","gamma")] <- theta[c("upsilon","gamma")]*fac
        return(list(theta = theta, qratio = qratio, id = 5))
    }
    betaGamma <- function(theta){
        fac <- exp(rnorm(3)*S)
        qratio <- prod(fac)
        theta["gamma"] <- theta["gamma"]*fac[1]
        theta[c("beta_t1", "beta_t2")] <- theta[c("beta_t1", "beta_t2")]/fac[2:3]
        return(list(theta = theta, qratio = qratio, id = 6))
    }
    all <- function(theta){
        fac <- exp(rnorm(4)*S)
        qratio <- prod(fac)
        theta <- theta*fac
        return(list(theta = theta, qratio = qratio, id = 7))
    }
    return(list(upsilon = upsilon, beta_t1 = beta_t1, beta_t2 = beta_t2, gamma = gamma,
                upsilonBeta = upsilonBeta, upsilonGamma = upsilonGamma, betaGamma = betaGamma,
                all = all))
}

##' Proposal functions for all parameters in the model
##' Kernel function. Takes a value and alters it with a normal random walk.
##' @return names list with names lists for each parameter.
##' @export
Proposal_kernel_SISe_random_lognormal <- function(theta, extraArgs){
    stopifnot(c("S","Svec","rho","rhoBeta") %in% names(extraArgs))
    S <- extraArgs$S
    Svec <- extraArgs$Svec
    rho <- extraArgs$rho
    rhoBeta <- extraArgs$rhoBeta


    upsilon <- function(S, rho, rhoBeta, theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["upsilon"] <- theta["upsilon"]*qratio
        return(list(theta = theta, qratio = qratio, id = 1))
    }
    beta_t1 <- function(S, rho, rhoBeta, theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t1"] <- rnorm(1, mean = theta["beta_t1"], sd = S)
        return(list(theta = theta, qratio = qratio, id = 2))
    }
    beta_t2 <- function(S, rho, rhoBeta, theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["beta_t2"] <- theta["beta_t2"]*qratio
        return(list(theta = theta, qratio = qratio, id = 3))
    }
    gamma <- function(S, rho, rhoBeta, theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["gamma"] <- theta["gamma"]*qratio
        return(list(theta = theta, qratio = qratio, id = 4))
    }
    upsilonBeta <- function(S, rho, rhoBeta, theta){
        step <- MASS::mvrnorm(n = 1,
                              mu = rep(0,3),
                              Sigma = matrix(c(1, -rho, -rho,
                                               -rho, 1, rhoBeta,
                                               -rho, rhoBeta, 1),ncol=3))
        expStep <- exp(step*S)
        qratio <- prod(expStep)

        theta[c("upsilon","beta_t1", "beta_t2")] <- theta[c("upsilon","beta_t1", "beta_t2")]*qratio
        return(list(theta = theta, qratio = qratio, id = 5))
    }
    upsilonGamma <- function(S, rho, rhoBeta, theta){
        step <- MASS::mvrnorm(n = 1,
                              mu = rep(0,2),
                              Sigma = matrix(c(1, rho,
                                               rho, 1), ncol=2))
        expStep <- exp(step*S)
        qratio <- prod(expStep)

        theta[c("upsilon","gamma")] <- theta[c("upsilon","gamma")]*expStep
        return(list(theta = theta, qratio = qratio, id = 6))
    }
    betaGamma <- function(S, rho, rhoBeta, theta){
        step <- MASS::mvrnorm(n = 1,
                              mu = rep(0,3),
                              Sigma = matrix(c(1, rho, rho,
                                               rho, 1, rhoBeta,
                                               rho, rhoBeta, 1),ncol=3))
        expStep <- exp(step*S)
        qratio <- prod(expStep)

        theta[c("beta_t1", "beta_t2", "gamma")] <- theta[c("beta_t1", "beta_t2", "gamma")]*expStep

        return(list(theta = theta, qratio = qratio, id = 7))
    }
    all <- function(S, rho, rhoBeta, theta){
        expStep <- exp(rnorm(n = 4)*S)
        qratio <- prod(expStep)
        theta[c("upsilon","beta_t1", "beta_t2", "gamma")] <-
            theta[c("upsilon","beta_t1", "beta_t2", "gamma")]*expStep
        return(list(theta = theta, qratio = qratio, id = 8))
    }

    ## probVec <- c(0.05,0.05,0.05,0.05, ## 0.2
    ##              0.2,0.2,0.2, ## 0.6
    ##              0.2) ## 0.2
    probVec <- rep(1,8)
    probVec <- probVec / sum(probVec)

    ## just a check to make sure!
    stopifnot( sum(probVec) == 1)

    funVec <- c(upsilon, beta_t1, beta_t2, gamma, upsilonBeta, upsilonGamma, betaGamma, all)
    funElement <- sample(x = seq_len(length(funVec)), size = 1, prob = probVec)
    func <- funVec[[funElement]]

    if(is.null(S))
        S <- Svec[funElement]
    fun.out <- func(S = S, rho = rho, rhoBeta = rhoBeta, theta = theta)

    return(list(theta = fun.out$theta, qratio = fun.out$qratio, id = fun.out$id))
}


##' Proposal functions for all parameters in the model
##' Kernel function. Takes a value and alters it with a normal random walk.
##' @return names list with names lists for each parameter.
##' @export
Proposal_SISe_Kernel_normalRandom <- function(theta, extraArgs){
    stopifnot(c("S") %in% names(extraArgs))
    S <- extraArgs$S

    upsilon <- function(S, theta){
        qratio <- 1
        theta["upsilon"] <- rnorm(1, mean = theta["upsilon"], sd = S)
        return(list(theta = theta, qratio = qratio, id = 1))
    }
    beta_t1 <- function(S, theta){
        qratio <- 1
        theta["beta_t1"] <- rnorm(1, mean = theta["beta_t1"], sd = S)
        return(list(theta = theta, qratio = qratio, id = 2))
    }
    beta_t2 <- function(S, theta){
        qratio <- 1
        theta["beta_t2"] <- rnorm(1, mean = theta["beta_t2"], sd = S)
        return(list(theta = theta, qratio = qratio, id = 3))
    }
    gamma <- function(S, theta){
        qratio <- 1
        theta["gamma"] <- rnorm(1, mean = theta["gamma"], sd = S)
        return(list(theta = theta, qratio = qratio, id = 4))
    }
    upsilonBeta <- function(S, theta){
        qratio <- 1
        rho <- 1/3
        step <- MASS::mvrnorm(n = 1,
                        mu = theta[c("upsilon","beta_t1", "beta_t2")],
                        Sigma = matrix(c(S^2, -rho*S^2, -rho*S^2,
                                         -rho*S^2, S^2, -rho*S^2,
                                         -rho*S^2, -rho*S^2,S^2),ncol=3))
        theta[c("upsilon","beta_t1", "beta_t2")] <- step
        return(list(theta = theta, qratio = qratio, id = 5))
    }
    upsilonGamma <- function(S, theta){
        qratio <- 1
        rho <- 1/2
        step <- MASS::mvrnorm(n = 1,
                        mu = theta[c("upsilon","gamma")],
                        Sigma = matrix(c(S^2,rho*S^2,rho*S^2,S^2),ncol=2))
        theta[c("upsilon","gamma")] <- step
        return(list(theta = theta, qratio = qratio, id = 5))
    }
    betaGamma <- function(S, theta){
        qratio <- 1
        rho <- 1/3
        step <- MASS::mvrnorm(n = 1,
                        mu = theta[c("beta_t1", "beta_t2","gamma")],
                        Sigma = matrix(c(S^2, -rho*S^2, -rho*S^2,
                                          -rho*S^2, S^2, -rho*S^2,
                                          -rho*S^2, -rho*S^2,S^2),ncol=3))
        theta[c("beta_t1", "beta_t2", "gamma")] <- step

        return(list(theta = theta, qratio = qratio, id = 6))
    }
    all <- function(S, theta){
        step <- rnorm(n = 4, mean = theta, sd = S)
        qratio <- 1
        theta[c("upsilon","beta_t1", "beta_t2", "gamma")] <- step
        return(list(theta = theta, qratio = qratio, id = 7))
    }

    probVec <- c(0.075,0.075,0.075,0.075, ## 0.3
                 0.20,0.20,0.20, ## 0.6
                 0.1) ## 0.1
    probVec <- probVec / sum(probVec)

    ## just a check to make sure!
    stopifnot( sum(probVec) == 1)

    funVec <- c(upsilon, beta_t1, beta_t2, gamma, upsilonBeta, upsilonGamma, betaGamma, all)
    fun <- sample(x = funVec, size = 1, prob = probVec)

    fun.out <- fun[[1]](S = S, theta = theta)

    return(list(theta = fun.out$theta, qratio = fun.out$qratio, id = fun.out$id))
}

##' Proposal functions for all parameters in the model
##' Kernel function. Takes a value and alters it with a log-normal random walk.
##' @return names list with names lists for each parameter.
##' @export
Proposal_SISe3sp_Kernel_fix <- function(){
    S <- 0.05
    upsilon <- function(theta){
        qratio <- exp(rnorm(n=1)*S)
        theta["upsilon"] <- theta["upsilon"]*qratio
        return(list(theta = theta, qratio = qratio))
    }
    return(list(upsilon = upsilon))
}

## Robin Eriksson 2018

##' Create the model before running a trajectory
##' @param theta (named) model parameters {up1,up2,up3,beta1,beta2,beta3,beta4}
##' @param u0 initial infection data
##' @param dist distance matrix connecting nodes
##' @param tspan time (integer) vector used in simulator
##' @return an initialised model
##' @export
create_model_SISe <- function(theta, u0, distance, tspan, events)
{

    numberNodes = dim(u0)[1]
    mod <- SimInf::SISe_sp(u0 = u0, tspan = tspan, phi = rep(0,numberNodes),
                   events = events,
                   upsilon = theta["upsilon"],
                   gamma = theta["gamma"],
                   beta_t1 = theta["beta_t1"], beta_t2 = theta["beta_t2"],
                   beta_t3 = theta["beta_t3"], beta_t4 = theta["beta_t4"],
                   end_t1 = 91, end_t2 = 182, end_t3 = 273, end_t4 = 365,
                   distance = distance, alpha = 1,
                   coupling = 0)

    mod
}

##' Implementation by Stefan Widgren with alterations
##' Init the SISe Model
##'
##' Every 1/p-th animal is moved to infected state.
##' @export
init_model_widgren <- function(model, theta, tspan, p = NULL, phi = NULL)
{
    if("prev" %in% names(theta))
        p <- theta["prev"]
    else if (is.null(p))
        p <- 0.04
    ##betaVec <- c(theta["beta_t1"], theta["beta_t2"], theta["beta_t3"],theta["beta_t4"])

    ##phi <- as.numeric(1/meanbetaVec) - theta["gamma"] / theta["upsilon"])
    ##phi <- 0
    if(is.null(phi))
        phi <- max(c(0,as.numeric(1/theta["beta_t1"] - theta["gamma"] / theta["upsilon"])))

    ##phi <- 0

    ##print(phi)
    ##phi <- 0.1
    ## Set tspan
    ##stopifnot(tspan[1] >= model@tspan[1] && tail(tspan,1) <= tail(model@tspan,1))
    model@tspan <- tspan

    ## Move animals to the infected compartment according to the
    ## prevalence by moving every 1/p animal from S to I
    model@u0["S", ] <- model@u0["S", ] + model@u0["I", ] -> Sigma
    model@u0["I", ] <- 0L
    i <- which(model@u0["S", ] > 0)

    ## pinv <- abs(1/p) ## bug somewhere. This maybe fixes it.
    ## u <- seq(from = sample(seq_len(pinv), 1),
    ##          to = sum(model@u0["S", i]),
    ##          by = pinv)
    ## w <- cumsum(model@u0["S", i])
    ## index <- 1
    ## k <- numeric(length(u))
    ## for (j in seq_len(length(u))) {
    ##     repeat {
    ##         if (u[j] <= w[index])
    ##             break
    ##         index <- index + 1
    ##     }
    ##     k[j] <- i[index]
    ## }
    ## k <- table(k)
    ## model@u0["I", as.integer(names(k))] <- as.integer(k)
    ircpp <- nodeIniter(as.numeric(model@u0["S",]), p)
    model@u0["I", ] <- ircpp

    model@u0["S", ] <- model@u0["S", ] - model@u0["I", ]

    model@gdata["upsilon"] <- theta["upsilon"]
    model@gdata["gamma"] <- theta["gamma"]
    model@gdata["beta_t1"] <- theta["beta_t1"]
    model@gdata["beta_t2"] <- theta["beta_t2"]
    model@gdata["beta_t3"] <- theta["beta_t3"]
    model@gdata["beta_t4"] <- theta["beta_t4"]
    model@gdata["coupling"] <- 0 ## no diffusion between nodes.
    ##model@v0[1, as.integer(names(k))] <- phi
    if(phi == "local") {
        locPhi <- ircpp / (theta["beta_t1"] * Sigma)
        locPhi[which(is.nan(locPhi))] <- 0
        model@v0["phi", ] <- locPhi
    } else
        model@v0["phi", ] <- phi

    model

}

library(SimInfInference)



##' Get the "top" posterior from inference class.
##' @param infe inference class
##' @param abc is the inference class an ABC?
##' @param tag the method that was used.
getTop <- function(infe, rows, abc = TRUE, tag = NA){
    d <- infe$getPosterior()
    cols <- dim(d)[2]
    if(abc){
        o <- order(d$delta)
        top <- d[o[seq_len(length(rows))],1:(cols-1)]
    } else {
        top <- d[rows, 1:(cols-1)]
    }
    top$method <- rep(tag, dim(top)[1])

    return(top)
}

getFun <- function(name, post, fun = mean) {
    cols <- dim(post)[2]

    out <- apply(post[post$method == name, seq_len(cols-1)], 2, function(x){fun(as.numeric(x))})

    return(out)
}



main <- function(rows = seq(1,14000,100), bootEstimate = TRUE) {
    post <- data.frame()

    ## load data
    load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sfull.RData")
    sfull <- getTop(sfull, rows, FALSE, "SLAM")
    post <- rbind(post, sfull)

    load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sbin.RData")
    sbin <- getTop(sbin, rows, FALSE, "SLAM / filtered")
    post <- rbind(post, sbin)

    load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sobs.RData")
    sobs <- getTop(sobs, rows, FALSE, "SLAM / observations")
    post <- rbind(post, sobs)




    cols <- dim(post)[2]
    methodnames <- as.list(unique(post$method))
    paramnames <- names(post[,1:(cols-1)])


    thetaTrue <- c(upsilon = 0.017154161,
                   beta_t1 = 0.15860932,
                   beta_t2 = 0.14621449,
                   beta_t3 = 0.15045281,
                   gamma = 0.10044851,
                   prev = 0.02027267)

    post <- rbind(post, c(as.numeric(thetaTrue), "theTrue"))


    ## compute Bias
    means <- c()
    for(name in methodnames) {
        means <-  c(means,getFun(name, post, mean))
    }

    means <- matrix(means, nrow = length(paramnames))


    bias <- sweep(means, 1, thetaTrue, '-')
    rownames(bias) <- paramnames
    colnames(bias) <- as.character(methodnames)


    ## input bias value
    if(bootEstimate) {
        load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/bias.RData")
        bias[,which(colnames(bias) == "SLAM / observations")] = b$bias
    } else
        bias[,which(colnames(bias) == "SLAM / observations")] = bias[,which(colnames(bias) == "SLAM / filtered")]


    ## compute variance
    variance <- c()
    for(name in methodnames) {
        variance <-  c(variance, getFun(name, post, var))
    }
    variance.mat <- matrix(variance, nrow = length(paramnames))

    rownames(variance.mat) <- paramnames
    colnames(variance.mat) <- as.character(methodnames)

    CoV <- sqrt(variance.mat) / means

    MSE <- variance.mat + bias^2


    NMRSE <- sqrt(MSE)/means

    return(list(bias = bias,
                variance = variance.mat,
                MSE = MSE,
                NMRSE = NMRSE,
                CoV = CoV)
           )
}

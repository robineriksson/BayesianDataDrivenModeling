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



main <- function(rows = seq(100,10000,100)) {
    post <- data.frame()

    ## load data
    load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/publish/beta2/sfull.RData")
    sfull <- getTop(sfull, rows, FALSE, "SLAM")
    post <- rbind(post, sfull)

    load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/publish/beta2/sbin.RData")
    sbin <- getTop(sbin, rows, FALSE, "SLAM / filtered")
    post <- rbind(post, sbin)

    load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/toPublish/abcfull.RData")
    abcfull <- getTop(abcfull, rows, TRUE, "ABC")
    post <- rbind(post, abcfull)

    load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/toPublish/abcbin.RData")
    abcbin <- getTop(abcbin, rows, TRUE, "ABC / filtered")
    post <- rbind(post, abcbin)


    cols <- dim(post)[2]
    methodnames <- as.list(unique(post$method))
    paramnames <- names(post[,1:(cols-1)])

    thetaTrue = c(upsilon = 0.0075,
                  beta_t1 = 0.05,
                  beta_t2 = 0.085,
                  gamma = 0.1)


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

    ## compute variance
    variance <- c()
    for(name in methodnames) {
        variance <-  c(variance, getFun(name, post, var))
    }
     variance.mat <- matrix(variance, nrow = length(paramnames))

    rownames(variance.mat) <- paramnames
    colnames(variance.mat) <- as.character(methodnames)

    variance.dimless <- (sqrt(variance.mat)/means)^2

    MSE <- variance.mat + bias^2

    NMRSE <- sqrt(MSE)/means

    return(list(bias = bias,
                variance = variance.mat,
                variance.dimless = variance.dimless,
                MSE = MSE,
                NMRSE = NMRSE)
           )
}

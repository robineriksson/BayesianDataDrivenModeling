## Robin Eriksson 2019
library(SimInfInference)

library(ggplot2)

##' Find the quantiles of the data x for the Bayesian CI
##' @param x the data to find the CI in
##' @param alpha parameter for deciding the CI region.
findQuantile <- function(x,alpha = 0.05) {
    xs <- sort(x)
    cpxx <- cumsum(xs) / sum(xs)
    upper = 1-alpha/2
    lower = 1-upper
    low <- xs[which(cpxx >= lower)[1]]   # lower boundary
    up <- xs[which(cpxx >= upper)[1]-1] # upper boundary

    return(c("low"=low,"up"=up))
}

##' Exploit the findQuantile function to extract the quantiles
##' for an entire matrix.
quantileHelper <- function(sample) {
    ## on the data, extract
    ## (1) mean
    ## (2) 95% Credible Interval (upper and lower)

    out <- data.frame("mean" = apply(sample, 2, mean))
    quant  <- apply(sample, 2, findQuantile)
    out$low  <- quant["low",]
    out$up  <- quant["up",]

    return(out)
}

##' Helper function to plot the posterior fit
plotit <- function(output, obs.y, textsize = 10) {

    x <- 1:length(names(obs.y))

    ## prepare sample data
    output.mat <- do.call("rbind",output)

    ## extract necessary information
    data <- quantileHelper(output.mat)

    data$x <- seq_len(length(names(obs.y)))
    data$obs <- obs.y

    g <- ggplot(data)


    ## ribbon
    g <- g + geom_ribbon(mapping = aes(x=x, ymin=low, ymax=up, color = "95% CI"), alpha=0.2)

    ## mean
    g <- g + geom_line(mapping = aes(x=x, y=mean, color = "Mean"), lwd = 2)

    ## observations
    g <- g + geom_point(mapping = aes(x=x, y=obs, color = "Observation"), pch = 17, cex = 4)#, col="cornflowerblue")


    ## change the x-axis test
    g <- g + scale_x_discrete(limits = names(obs.y))

    ## labels
    g <- g + xlab("Time") + ylab("Positive samples") + labs(color="Source")

    ## add publication theme
    g <- g + theme_Publication(textsize) +
        theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
              legend.position = "top",
              text = ggplot2::element_text((family="LM Roman 10"))) +
        scale_color_manual(values=c("grey80", "black", "cornflowerblue"))

    return(g)
}


##' Evaluate how well the posterior fits the data
##'
##' @param N number of posterior sampels
##' @param filename  name and filepath of the saved posterior
##' @param datadir path to the directory holding data used in simulation.
evalPosterior <- function(N = 250, filename = "/posterior/real/mis_obs.rda",
                          dataDir = "../DATA") {
    paths <- NULL
    ## load observation
    paths["obs"] <- paste(dataDir, "/VTEC/obs.rda", sep="")
    load(paths["obs"]) ## loads: obs
    obs$sample <- as.numeric(obs$status)

    ## aggregate
    aggregate <- function(obs) {
        obs.q <- qtrStore(obs, "time", "sample", logical=FALSE)
        obs.m <- matrix(obs.q$sample, ncol = 126)
        obs.y <- apply(obs.m, 1, function(x){sum(x, na.rm=TRUE)})
        names(obs.y) <- unique(obs.q$time)
        return(obs.y)
    }
    obs.y <- aggregate(obs)

    ## load posterior
    paths["post"] <- paste(dataDir, filename, sep="")
    load(paths["post"])
    mobs$changeExtraArgsSimulator(nSim=1)
    post <- mobs$getPosterior()

    ## sample from the parameter space
    rows <- seq(1,13000,10)
    tosample <- sample(rows, N, replace=TRUE)
    post <- post[tosample, -7]

    ## simulate 1 realisation per parameter sample
    output <- list()
    for(i in seq_len(N)) {
        theta <- as.numeric(post[i,])
        names(theta) <- colnames(post)
        sim <- mobs$runSimulator(theta)
        sim.y <- aggregate(sim[[1]])
        output[[i]] <- sim.y
    }


    p <- plotit(output, obs.y)

    return(list("sample" = output, "obs" = obs.y, "plot" = p))
}

## Robin Eriksson 2019
library(SimInfInference)

##' Evaluate how well the posterior fits the data
##'
##' @param N number of posterior sampels
##' @param filename  name and filepath of the saved posterior
##' @param datadir path to the directory holding data used in simulation.
evalPosterior <- function(N = 250, filename, datadir) {
    paths <- NULL
    ## load observation
    paths["obs"] <- paste(dataDir, "obs.rda", sep="")
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
    paths["post"] <- filename
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


##' Helper function to plot the posterior fit
plotit <- function(output, obs.y, textsize = 10) {
    x <- 1:length(names(obs.y))

    ## prepare sample data
    output.mat <- do.call("rbind",output)

    ## extract necessary information
    data <- quantileHelper(output.mat)

    data$x <- seq_len(length(names(obs.y)))
    data$obs <- obs.y

    data.m <- reshape2::melt(data, id.var = "x")

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

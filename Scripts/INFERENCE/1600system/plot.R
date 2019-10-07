library(SimInfInference)
library(ggplot2)
library(gridExtra)
library(car)
library(cowplot)
library(GGally)
library(ggpubr)

##' Helper function, diagonal plot
m1density <- function(data, mapping, thetaTrue, lims = NULL, ...) {
    g <- ggplot2::ggplot(data = data, mapping = mapping)

    g <- g + geom_density(adjust = 2, lwd = 0.25)

    if(!is.null(lims)){
        x <- lims[,colnames(lims) == rlang::get_expr(mapping[[1]])]
        g <- g + xlim(x)
    }


    if(!is.null(thetaTrue))
        g <- g + geom_vline(xintercept = thetaTrue[as.character(rlang::get_expr(mapping[[1]]))],
                            linetype = 2, lwd = 0.25)


    cols <- getColors(length(unique(data$method)),"1600")

    return(g)
}

##' Helper function, lower triangle plots
myEllips <- function(data, mapping, thetaTrue, lims = NULL, noellips = TRUE, ...){
    g <- ggplot(data = data, mapping = mapping) +
        scale_shape(solid=T)


    alpha <- 0.5
    size <- 1

    g <- g + geom_point(mapping = aes(shape = method), alpha = alpha, size = size)


    if(!noellips)
        g <- g +  stat_ellipse(level = 0.95, lwd = 0.25)



    if(!is.null(lims)){
        x <- lims[,colnames(lims) == rlang::get_expr(mapping[[1]])]
        y <- lims[,colnames(lims) == rlang::get_expr(mapping[[2]])]
        g <- g + xlim(x) + ylim(y)
    }

    if(!is.null(thetaTrue)) {
        g <- g + geom_vline(xintercept = thetaTrue[as.character(rlang::get_expr(mapping[[1]]))],
                            linetype = 2,lwd = 0.25)
        g <- g + geom_hline(yintercept = thetaTrue[as.character(rlang::get_expr(mapping[[2]]))],
                            linetype = 2,lwd = 0.25)
    }



    cols <- getColors(length(unique(data$method)),"1600")


    return(g)
}

##' Helper function
##'
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

##' Helper function
##'
##' Extracts the posteriors, calls the correct functions. returns plot object
##' @param names the names of the inference objects to load
##' @param rows thinning of the posterior
##' @param thetaTrue a "true" theta to plot
##' @param selectedParams what parameters should be plot?
##' @param includeEllips Gaussian ellipses around the scatterplots.
plotMultiDens <- function(names, rows = seq(100,10000,100), thetaTrue = NULL, selectedParams = NULL,
                          noellips = TRUE){



    post <- data.frame()
    for(name in names){
        ## load data
        dirname <- "../../DATA/posterior/1600/"

        if(name == "sfull") {
            ##load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/publish/beta2/sfull.RData")
            load(paste(dirname,name,".RData",sep=""))

            load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/sfull.RData")
            sfull <- getTop(sfull, rows, FALSE, "SLAM / unfiltered")
            post <- rbind(post, sfull)
        } else if(name == "sbin") {
            ##load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/publish/beta2/sbin.RData")
            load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/sbin.RData")
            sbin <- getTop(sbin, rows, FALSE, "SLAM / filtered")
            post <- rbind(post, sbin)
        } else if(name == "abcfull") {
            load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcfull.RData")
            abcfull <- getTop(abcfull, rows, TRUE, "ABC / unfiltered")
            post <- rbind(post, abcfull)
        } else if(name == "abcbin") {
            load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcbin.RData")
            abcbin <- getTop(abcbin, rows, TRUE, "ABC / filtered")
            post <- rbind(post, abcbin)
        } else {
            print("name not accounted for in code")
            stopifnot(FALSE)
        }
    }





    numParam <- dim(post)[2]-1

    if(is.null(selectedParams))
        columns <- colnames(post)[1:numParam]
    else
        columns <- selectedParams
    ## selected parameters + SL column
    post <- post[,c(which(colnames(post) %in% columns), numParam + 1)]
    rename <- c()
    for(name in columns) {
        if(name == "beta_t1")
            rename <- c(rename, "beta[1]")
        else if(name == "beta_t2")
            rename <- c(rename, "beta[2]")
        else
            rename <- c(rename, name)
    }
    columnLabels <- rename

    numParam <- length(columns)

    limits <- rbind(as.numeric(apply(post,2,min)[1:numParam])*0.99,
                    as.numeric(apply(post,2,max)[1:numParam])*1.01)

    colnames(limits) <- columns

    post$method  <- factor(post$method, levels = c("ABC / unfiltered",
                                                   "ABC / filtered",
                                                   "SLAM / unfiltered",
                                                   "SLAM / filtered"))



    p <- GGally::ggpairs(post,
                         mapping = aes(color = method, linetype = method),
                         columns = columns,
                         columnLabels = columnLabels,
                         labeller = "label_parsed",
                         legend = c(1,1),
                         upper = list(continuous = "blank"),
                         diag = list(continuous = wrap(m1density, thetaTrue = thetaTrue, lims = limits)),
                         lower = list(continuous = wrap(myEllips, thetaTrue = thetaTrue, lims = limits,
                                                        noellips = noellips))
                         ) +
        theme_Publication() +
        theme(legend.position="top")


    return(p)
}

##' Main function
##'
##' Creates the ABC posterior plot.
main <- function(){
    ## which posteriors to use
    names <- c("abcfull", "abcbin", "sfull", "sbin")

    ## how to thin the plot
    rows <- seq(15000,30000,100)

    ## Only plot some parameter dimensions?
    selectedParams <- NULL
    ## show ellips or not on scatter plot
    noellips <- TRUE

    ## true value.
    thetaTrue <- c(upsilon = 0.005,
                  beta_t1 = 0.025,
                  beta_t2 = 0.058,
                  gamma = 0.1)

    ## create the plot object
    p <- plotMultiDens(names, rows, thetaTrue, selectedParams, noellips=noellips)

    p <- p + theme(text = element_text(family="LM Roman 10"))
    return(p)
}

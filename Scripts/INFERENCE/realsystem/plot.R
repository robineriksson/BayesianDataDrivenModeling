##' Robin Eriksson 2019

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


    cols <- getColors(length(unique(data$method)),"real")#[c(1,3,2)]
    g <- g + scale_color_manual(values=cols)

    if(!is.null(thetaTrue))
        g <- g + geom_vline(xintercept = thetaTrue[as.character(rlang::get_expr(mapping[[1]]))],
                            linetype = 2, lwd = 0.25)

    return(g)
}

##' Helper function, lower triangle plots
myEllips <- function(data, mapping, thetaTrue, lims = NULL, includeEllips = FALSE, ...){
    g <- ggplot(data = data, mapping = mapping)

    alpha  <- 0.7
    size <- 0.5

    g <- g + geom_point(mapping = aes(shape = method), alpha = alpha, size = size)

    if(includeEllips)
        g <- g +  stat_ellipse(level = 0.95, lwd = 0.25)

    g <- g + scale_shape(solid=F)

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

    cols <- getColors(length(unique(data$method)),"real")
    g <- g + scale_color_manual(values=cols)

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
##' @param infe if we don't give names, we supply the objects directly.
##' @param includeEllips Gaussian ellipses around the scatterplots.
plotMultiDens <- function(names, rows = seq(100,10000,100), thetaTrue = NULL,
                          selectedParams = NULL, infe = NULL,
                          includeEllips=TRUE){


    post <- data.frame()
    if(is.null(infe)) {
        dirname <- "../../DATA/posterior/real/"
        for(name in names){
            ## load data
            if(name == "mis_full") {
                load(paste(dirname,name,".rda",sep=""))
                mfull <- getTop(mfull, rows, FALSE, "SLAM / filtered")
                post <- rbind(post, mfull)
            } else if(name == "mis_bin") {
                load(paste(dirname,name,".rda",sep=""))
                mbin <- getTop(mbin, rows, FALSE, "SLAM / unfiltered")
                post <- rbind(post, mbin)
            } else if(name == "mis_obs") {
                load(paste(dirname,name,".rda",sep=""))
                mobs <- getTop(mobs, rows, FALSE, "SLAM / observations")
                post <- rbind(post, mobs)
            } else {
                print("name not accounted for in code")
                stopifnot(FALSE)
            }
        }
    } else {
        if(is.list(infe)) {
            for (i in 1:length(infe)) {
                name <- names(infe)[i]
                post <- rbind(post, getTop(infe[[i]], rows, FALSE, name))
            }
        } else
            post <- getTop(infe, rows, FALSE, "SLAM")
    }

    ## remove potentials NaN
    nans <- which(is.na(post))
    if(!(length(nans) == 0))
        post <- post[-nans,]

    numParam <- dim(post)[2]-1

    ## which columns and columnlabels ...
    if(is.null(selectedParams)) {
        if(numParam == 4) { ## 2 beta
            columns = c("upsilon", "beta_t1", "beta_t2", "gamma")
            columnLabels = c("upsilon", "beta[1]", "beta[2]", "gamma")
        } else if(numParam == 5) { ## 3 beta
            columns = c("upsilon", "beta_t1", "beta_t2", "beta_t3", "gamma")
            columnLabels = c("upsilon", "beta[1]", "beta[2]", "beta[3]", "gamma")
        } else if(numParam == 6) { ## 3 beta & prev
            columns = c("upsilon", "beta_t1", "beta_t2", "beta_t3", "gamma", "prev")
            columnLabels = c("upsilon", "beta[1]", "beta[2]", "beta[4]", "gamma", "p[0]")
        } else {stopifnot(FALSE)}
    } else {
        columns <- selectedParams
        ## selected parameters + SL column
        post <- post[,c(which(colnames(post) %in% columns), numParam + 1)]
        rename <- c()
        old  <- TRUE
        if(old) {
            for(name in columns) {
                if(name == "beta_t1")
                    rename <- c(rename, "beta[1]")
                else if(name == "beta_t2")
                    rename <- c(rename, "beta[2]")
                ## this one is special
                else if(name == "beta_t3")
                    rename <- c(rename, "beta[4]")
                else if(name == "prev")
                    rename <- c(rename, "p[0]")
                else
                    rename <- c(rename, name)
            }
        } else {
            for(name in columns) {
                if(name == "beta_t1") {
                    rename <- c(rename, TeX("$\\beta_1$"))
                } else if(name == "beta_t2") {
                    rename <- c(rename, TeX("$\\beta_2$"))
                    ## this one is special
                } else if(name == "beta_t3") {
                    rename <- c(rename, TeX("$\\beta_4$"))
                } else if(name == "prev") {
                    rename <- c(rename, TeX("$p_0$"))
                } else if(name == "upsilon") {
                    rename <- c(rename, TeX("$\\upsilon$"))
                } else {
                    rename <- c(rename, name)
                }
            }
        }
        columnLabels <- rename

        numParam <- length(columns)


    }


    lim.bot <- 0.95
    lim.top <- 1.05
    limits <- rbind(as.numeric(apply(post,2,min)[1:numParam])*lim.bot,
                    as.numeric(apply(post,2,max)[1:numParam])*lim.top)
    colnames(limits) <- columns

    post$method  <- factor(post$method, levels = c("SLAM / unfiltered",
                                                   "SLAM / filtered",
                                                   "SLAM / observations"))

    p <- GGally::ggpairs(post,
                         mapping = aes(color = method, linetype = method),
                         columns = columns,
                         columnLabels = columnLabels,
                         labeller = "label_parsed",
                         legend = c(1,1),
                         upper = list(continuous = "blank"),
                         diag = list(continuous = wrap(m1density, thetaTrue = thetaTrue, lims = limits)),
                         lower = list(continuous = wrap(myEllips, thetaTrue = thetaTrue,
                                                        lims = limits, includeEllips = includeEllips))
                         ) +
        theme_Publication() +
        theme(legend.position="top")


    return(p)
}

##' Main function.
##'
##' Creates the multi posterior plot.
main <- function() {
    ## which posteriors to use
    names <- c("mis_full", "mis_bin", "mis_obs")

    ## what thinning to use.
    rows <- seq(1,13000,10)

    ## what parameters to plot, NULL if all
    selectedParams <- c("upsilon", "beta_t1", "gamma", "prev")

    ## Can be modified so that the user supplies inference objects and not loads them.
    infe <- NULL

    ## If to use the Gaussian ellipses or not.
    includeEllips <- TRUE

    ## "true" parameters
    thetaTrue <- c(upsilon = 0.01518903,
                   beta_t1 = 0.14057693,
                   beta_t2 = 0.12646374,
                   beta_t3 = 0.15532180,
                   gamma = 0.09650432,
                   prev = 0.02471147)

    ## create the plot.
    p <- plotMultiDens(names, rows, thetaTrue =thetaTrue, selectedParams = selectedParams, infe = infe,
                       includeEllips = includeEllips)



    p <- p + theme(text = element_text(family="LM Roman 10"))

    return(p)
}

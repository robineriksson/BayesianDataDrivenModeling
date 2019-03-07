##' Robin Eriksson 2018
##' Plot the distributions in a comparative manner.
##'
##'

library(SimInfInference)
library(ggplot2)
library(gridExtra)
library(car)
library(cowplot)
library(GGally)


m1density <- function(data, mapping, thetaTrue, lims = NULL, ...) {
    g <- ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_density(adjust = 2, lwd = 0.25)

    if(!is.null(lims)){
        x <- lims[,colnames(lims) == mapping[[1]]]
        g <- g + xlim(x)
    }


    if(!is.null(thetaTrue))
        g <- g + geom_vline(xintercept = thetaTrue[as.character(mapping[[1]])], linetype = 2, lwd = 0.25)
    return(g)
}

## not used.
myPoints <- function(data, mapping, ...){
    g <- ggplot(data = data, mapping = mapping) +
        ##ggplot2::geom_point(alpha = 0.25) +
        stat_density_2d(geom = "contour", bins = 5)

    return(g)
}

myEllips <- function(data, mapping, thetaTrue, lims = NULL,...){
    g <- ggplot(data = data, mapping = mapping) +
        geom_point(aes(shape = method), alpha = 0.7, size = 0.5) +
        scale_shape(solid=F)
        ## stat_density2d(geom="density_2d",
        ##                aes(color = method, alpha=10*..level..),
        ##                #size=2,
        ##                contour=TRUE,
        ##                lwd = 0.25)
    #g <- g +  stat_ellipse(level = 0.95, lwd = 0.25)



    if(!is.null(lims)){
        x <- lims[,colnames(lims) == mapping[[1]]]
        y <- lims[,colnames(lims) == mapping[[2]]]
        g <- g + xlim(x) + ylim(y)
    }

    if(!is.null(thetaTrue)) {
        ## g <- g + geom_point(aes(x = thetaTrue[as.character(mapping[[1]])[1]],
        ##                         y = thetaTrue[as.character(mapping[[2]])[2]]),
        ##                     colour = "black")
        g <- g + geom_vline(xintercept = thetaTrue[as.character(mapping[[1]])],  linetype = 2,lwd = 0.25)
        g <- g + geom_hline(yintercept = thetaTrue[as.character(mapping[[2]])],  linetype = 2,lwd = 0.25)
      }
    return(g)
}

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


##' Get the "top" posterior from inference class.
##' @param infe inference class
##' @param abc is the inference class an ABC?
##' @param tag the method that was used.
getTopSL <- function(infe, rows, abc = TRUE, tag = NA){
    d <- infe$getPosterior()
    cols <- dim(d)[2]
    if(abc){
        o <- order(d$delta)
        top <- d[o[seq_len(length(rows))], cols]
    } else {
        top <- d[rows, cols]
    }

    top <- data.frame(t = 1:length(rows), sl = top, method = rep(tag, length(top)))

    return(top)
}



plotMultiDens <- function(names, rows = seq(100,10000,100), thetaTrue = NULL, selectedParams = NULL){


    post <- data.frame()
    for(name in names){
        ## load data
       if(name == "sobs") {
            load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sobs.RData")
            sobs <- getTop(sobs, rows, FALSE, "SLAM / observations")
            post <- rbind(post, sobs)
       } else if(name == "sbin") {
           load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sbin.RData")
           sbin <- getTop(sbin, rows, FALSE, "SLAM / filtered")
           post <- rbind(post, sbin)
       } else if(name == "sfull") {
           load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sfull.RData")
           sfull <- getTop(sfull, rows, FALSE, "SLAM")
           post <- rbind(post, sfull)
       } else if(name == "sboot") {
           load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sboot.RData")
           sboot <- as.data.frame(sboot[seq(1,dim(sboot)[1],1),1:6])
           sboot$method <- rep("boot - SLAM", dim(sboot)[1])
           post <- rbind(post, sboot)
       } else {
           print("name not accounted for in code")
           stopifnot(FALSE)
       }
    }

    ##post <- rbind(sl1, sl2, sl3)

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
        columnLabels <- rename

        numParam <- length(columns)


    }



    limits <- rbind(as.numeric(apply(post,2,min)[1:numParam])*0.85,
                    as.numeric(apply(post,2,max)[1:numParam])*1.15
                    )
    colnames(limits) <- columns
    #limits <- NULL


    p <- GGally::ggpairs(post,
                         mapping = aes(color = method, linetype = method),
                         columns = columns,
                         columnLabels = columnLabels,
                         labeller = "label_parsed",
                         legend = c(1,1),
                         upper = list(continuous = "blank"),
                         ##upper = list(continuous = "density"),
                         diag = list(continuous = wrap(m1density, thetaTrue = thetaTrue, lims = limits)),
                         lower = list(continuous = wrap(myEllips, thetaTrue = thetaTrue, lims = limits))
                         ) +
        theme_Publication() +
        theme(legend.position="top")


    return(p)
}


## plotDensitites <- function(){
##     ## load and re-store the data.
##     load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sfull.RData")
##     sl1 <- sfull
##     load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sbin.RData")
##     sl2 <- sbin
##     load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sobs.RData")
##     sl3 <- sobs



##     getExtrema <- function(listOfInfe, func = min, extra = 1){
##         data <- do.call("rbind",
##                         lapply(listOfInfe, function(x){x$getPosterior()[,1:4]})
##                         )
##         extrema <- apply(data, 2, func)
##         extrema <- extrema * extra
##         return(extrema)
##     }



##     posteriorPlot <- function(infe, n = 1, N = NULL, by = NULL, ordered = FALSE, limits = NULL){
##         data <- infe$getPosterior()
##         dims <- dim(data)

##         if(is.null(N))
##             N <- dims[1]

##         if(ordered){
##             data <- data[order(data[,dims[2]]),]
##         }

##         if(!is.null(by))
##             data <- data[seq(n,N,by),]
##         else
##             data <- data[n:N,]

##         my2dDensity <- function(data, mapping, lims = NULL,  ...){
##             g <- ggplot2::ggplot(data = data, mapping = mapping) +
##                 ggplot2::stat_density_2d(geom = "contour", bins = 5,
##                                          ggplot2::aes(color = ..level..))
##             if(!is.null(lims)){
##                 x <- lims[,colnames(lims) == mapping[[1]]]
##                 y <- lims[,colnames(lims) == mapping[[2]]]

##                 g <- g + xlim(x) + ylim(y)
##             }
##             return(g)
##         }

##         my1dDensity <- function(data, mapping, lims = NULL, ...) {
##             g <- ggplot2::ggplot(data = data, mapping = mapping) +
##                 ggplot2::geom_density()

##             if(!is.null(lims)){
##                 x <- lims[,colnames(lims) == mapping[[1]]]
##                 g <- g + xlim(x)
##             }
##             return(g)

##         }

##         myPoints <- function(data, mapping, lims = NULL, ...){
##             g <- ggplot2::ggplot(data = data, mapping = mapping) +
##                 ggplot2::geom_point()

##             if(!is.null(lims)){
##                 x <- lims[,colnames(lims) == mapping[[1]]]
##                 y <- lims[,colnames(lims) == mapping[[2]]]

##                 g <- g + xlim(x) + ylim(y)
##             }
##             return(g)
##         }



##         g <- GGally::ggpairs(data = data[,seq(1,dims[2]-1)],
##                              columnLabels = c("upsilon", "beta[1]", "beta[2]", "gamma"),
##                              labeller = "label_parsed",
##                              upper = list(continuous = GGally::wrap(my2dDensity, lims = limits)),
##                              diag = list(continuous = GGally::wrap(my1dDensity, lims = limits)),
##                              lower = list(continuous = GGally::wrap(myPoints, lims = limits))
##                              ) +
##             theme_Publication()


##     }

##     mins <- getExtrema(list(sl1, sl2, sl3), min, 0.9)
##     maxs <- getExtrema(list(sl1, sl2, sl3), max, 1.1)
##     limits <- rbind(mins, maxs)

##     p1 <- posteriorPlot(sl1, 250, 3000, 10, FALSE, limits = limits)
##     p2 <- posteriorPlot(sl2, 250, 3000, 10, FALSE, limits = limits)
##     p3 <- posteriorPlot(sl3, 250, 3000, 10, FALSE, limits = limits)

##     return(list(p1 = p1, p2 = p2, p3 = p3))
## }


plotEllips <- function(N = 250, level = 0.95){
    ## load and restoredata
    load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sfull.RData")
    sl1 <- sfull
    load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sbin.RData")
    sl2 <- sbin

    ## extract the posteriors
    ## full state sl
    d.sl1 <- sl1$getPosterior()[seq(500,3000,1),c("upsilon", "gamma")]
    n.sl1 <- sample(1:dim(d.sl1)[1], N, replace = TRUE)
    post.sl1 <- d.sl1[n.sl1,]
    post.sl1$type <- rep("sl-Full", N)

    ## binary state sl
    d.sl2 <- sl2$getPosterior()[seq(500,3000,1),c("upsilon", "gamma")]
    n.sl2 <- sample(1:dim(d.sl2)[1], N, replace = TRUE)
    post.sl2 <- d.sl2[n.sl2,]
    post.sl2$type <- rep("sl-binary", N)

    ## true value
    true.df <- data.frame(upsilon = rep(0.0075,N), gamma = rep(0.1,N), type = rep("actual",N))

    data <- data.frame(upsilon = c(post.sl1$upsilon, post.sl2$upsilon),
                       gamma = c(post.sl1$gamma, post.sl2$gamma),
                       type = c(post.sl1$type, post.sl2$type)
                       )
    data$type <- factor(data$type)

    ## plot!
    p <- ggplot()
    p <- p + SimInfInference::theme_Publication() +
        ##theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
        labs(x = expression(upsilon), y = expression(gamma)) #+
    ## ggtitle(sprintf("1600 system. %.0f posterior samples, %.0f%% confidence interval", N, 100*level)) +
    ## theme(plot.title = element_text(size = 18))

    p <- p + geom_point(data = data,
                        aes(x = upsilon, y = gamma, colour = type), size = 1) +
        stat_ellipse(data= data, aes(x = upsilon, y = gamma, colour = type), size = 1.5, level = level) +
        geom_point(aes(x = 0.0135, y = 0.1), shape = 13, size = 6)

    return(p)
}

plotSL <- function(names = c("sfull", "sbin", "sobs"), rows = seq(1,14000,100),
                   save = FALSE, saveInline = FALSE) {



    post <- data.frame()
    for(name in names){
        ## load data
        if(name == "sobs") {
            load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sobs.RData")
            sobs <- getTopSL(sobs, rows, FALSE, "sobs")
            post <- rbind(post, sobs)
        } else if(name == "sbin") {
            load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sbin.RData")
            sbin <- getTopSL(sbin, rows, FALSE, "sbin")
            post <- rbind(post, sbin)
        } else if(name == "sfull") {
            load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sfull.RData")
            sfull <- getTopSL(sfull, rows, FALSE, "sfull")
            post <- rbind(post, sfull)
        } else {
            print("name not accounted for in code")
            stopifnot(FALSE)
        }
    }


    post$method <- as.character(post$method)
    post[which(post$method == "sfull"),"method"] <-"SLAM"
    post[which(post$method == "sbin"),"method"] <- "SLAM / filtered"
    post[which(post$method == "sobs"),"method"] <- "SLAM / observations"



    g <- ggplot2::ggplot(post, ggplot2::aes(x = t, y = sl, color = method))
    g <- g + ggplot2::geom_line(ggplot2::aes(linetype = method), lwd = 0.25)

    colors <- c("#f8766dff","#00BA38", "#619CFF") #FFC34B
    linetypes <- c(1,3,2)


    ## for(i in seq_len(length(colors))){
    ##     subpost <- post[which(post$method == names[i]),]
    ##     g <- g + ggplot2::geom_line(data = subpost,
    ##                                 ggplot2::aes(x = t, y = sl),
    ##                                 color = colors[i], linetype = linetypes[i], lwd = 0.25)
    ## }

    g <- g + theme_Publication(6) + ggplot2::theme(legend.position="top",#c(1,0.5),
                                                   legend.direction = "horizontal",
                                        #legend.background = element_rect(fill="lightblue",
                                        #                                 size=1, linetype="solid"),
                                        #legend.key.width = unit(0.5,"cm"),
                                                   axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                                                   axis.title = ggplot2::element_text(face = "plain",size = ggplot2::rel(1)),
                                                   plot.margin=grid::unit(c(1,1,1,1),"mm"),
                                                   legend.margin=margin(l = -0.5, unit='cm'))
    ## plot.margin=grid::unit(c(0,0,0,0),"mm"),
    ## axis.line = element_line(size = 0.25),
    ## axis.ticks = element_line(colour = "black", size = 0.25),

    ##                                               )
    g <- g + ggplot2::labs(x = "Sample number", y = "Syntethic log likelihood")


    #g <- g + ggplot2::labs(x = "", y = "")


    if(save){
        dirname <- "~/Gits/BPD/PLOTS/toPublish/real/"
        if(saveInline) {
            g <- g + ggplot2::labs(x = "Sample", y = "SL")
            g <- g + ggplot2::theme(legend.position = "none",
                                    plot.margin=grid::unit(c(0,0,0,0),"mm"),
                                    axis.line = element_line(size = 0.25),
                                    axis.ticks = element_line(colour = "black", size = 0.25),
                                    axis.title.y = element_text(margin = margin(t = 0,
                                                                                r = 0,
                                                                                b = 0,
                                                                                l = 1,
                                                                                unit = "pt")),
                                    axis.title.x = element_text(margin = margin(t = -4,
                                                                                r = 0,
                                                                                b = 1,
                                                                                l = 0,
                                                                                unit = "pt"))
                                    )



            ggsave(paste(dirname, "realMultiSeries.pdf", sep = ""), g,  width = 3.8, height = 1.8, units = "cm")
        } else
            ggsave(paste(dirname, "realMultiSeries.pdf", sep = ""), g,  width = 8.7, height = 4.3, units = "cm")
    }



    return(g)


}

computeSLmoments <- function(names = c("sfull", "sbin", "sobs"), rows = seq(1,14000,100), save = FALSE) {
    post <- data.frame()
    for(name in names){
        ## load data
        if(name == "sobs") {
            load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sobs.RData")
            sobs <- getTopSL(sobs, rows, FALSE, "sobs")
            post <- rbind(post, sobs)
        } else if(name == "sbin") {
            load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sbin.RData")
            sbin <- getTopSL(sbin, rows, FALSE, "sbin")
            post <- rbind(post, sbin)
        } else if(name == "sfull") {
            load("~/Gits/BPD/R/INFERENCE/realsystem/output/slam/sfull.RData")
            sfull <- getTopSL(sfull, rows, FALSE, "sfull")
            post <- rbind(post, sfull)
        } else {
            print("name not accounted for in code")
            stopifnot(FALSE)
        }
    }


    post$method <- as.character(post$method)
    sl.full <- post[which(post$method == "sfull"),"sl"]
    sl.bin <- post[which(post$method == "sbin"),"sl"]
    sl.obs <- post[which(post$method == "sobs"),"sl"]

    sl.df <- data.frame(full = sl.full, bin = sl.bin, obs = sl.obs)
    means <- apply(sl.df, 2, mean)
    sd <- apply(sl.df, 2, sd)

    return(list(means = means, sd = sd))
}

## main <- function(){

##     ##ellips <- plotEllips(N = 5e2, level = 0.95)
##     densities <- plotDensities()
##     dirname <- "~/Gits/BPD/PLOTS/toPublish/real/"

##     ##ggsave(paste(dirname, "ellipses.pdf", sep = ""), ellips)
##     ##ggsave(paste(dirname, "abcPosterior1600FullState.pdf", sep = ""), densities$p1)
##     ##ggsave(paste(dirname, "abcPosterior1600BinaryState.pdf", sep = ""), densities$p2)
##     ggsave(paste(dirname, "slRealFullStatePosterior.pdf", sep = ""), densities$p1)
##     ggsave(paste(dirname, "slRealBinaryStatePosterior.pdf", sep = ""), densities$p2)
##     ggsave(paste(dirname, "slRealObsPosterior.pdf", sep = ""), densities$p3)


## }

main <- function(names = c("sfull", "sbin", "sobs"), rows = seq(1,14000,100),
                 save = FALSE, widgrenParam = FALSE,
                 selectedParams = c("upsilon", "beta_t1", "gamma", "prev")){
    if(widgrenParam)
        thetaTrue = c(upsilon = 1.75e-2,
                      beta_t1 = 1.57e-1,
                      beta_t2 = 1.44e-1,
                      beta_t3 = 1.50e-1,
                      gamma = 0.1,
                      prev = 0.02)
    else
        thetaTrue <- c(upsilon = 0.017154161,
                       beta_t1 = 0.15860932,
                       beta_t2 = 0.14621449,
                       beta_t3 = 0.15045281,
                       gamma = 0.10044851,
                       prev = 0.02027267)
    p <- plotMultiDens(names, rows, thetaTrue = thetaTrue, selectedParams = selectedParams)
    if(save){
        p <- p + theme_Publication(6) +
            theme(plot.margin=grid::unit(c(0,0,0,0),"mm"),
                  legend.margin=margin(t = 0.1,b=-0.25, unit='cm'),
                  axis.line = element_line(size = 0.25),
                  axis.ticks = element_line(colour = "black", size = 0.25),
                  panel.grid.major = ggplot2::element_line(colour="#f0f0f0", size = 0.25),
                  legend.position = "top")

        dirname <- "~/Gits/BPD/PLOTS/toPublish/real/"
        ggsave(paste(dirname, "realMultiPosterior.pdf", sep = ""), p,  width = 8.7, height = 8.7, units = "cm")
    }
    return(p)
}

## mainMulti <- function(names = c("sfull", "sbin", "sobs"), rows = seq(1,10000,100),
##                       save = FALSE) {
##     library(egg)
##     post <- main3(names,rows)
##     traj <- plotSL(names,rows)

##     plot.with.inset <- post + annotation_custom(ggplotGrob(traj), xmin = 5, xmax = 7, ymin = 30, ymax = 44)

##     plot.with.inset


## }

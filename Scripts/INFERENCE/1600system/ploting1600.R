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
library(ggpubr)

## plotMultiDens_old<- function(){

##     ## load data
##     load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcfull.RData")
##     abc1 <- getTop(abcfull, TRUE, "abcfull")
##     load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcbin.RData")
##     abc2 <- getTop(abcbin, TRUE, "abcbin")
##     load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/sfull.RData")
##     sl1 <- getTop(sfull, FALSE, "sfull")
##     load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/sbin.RData")
##     sl2 <- getTop(sbin, FALSE, "sbin")

##     post <- rbind(abc1, abc2, sl1, sl2)
##     ## post <- rbind(post, sl1)
##     ## post <- rbind(post, sl2)



##     p <- GGally::ggpairs(post,
##                  mapping = aes(color = method, linetype = method),
##                  columns = c("upsilon", "beta_t1", "beta_t2", "gamma"),
##                  columnLabels = c("upsilon", "beta[1]", "beta[2]", "gamma"),
##                  labeller = "label_parsed",
##                  legend = c(1,1),
##                  upper = list(continuous = "blank"),
##                  ##upper = list(continuous = "density"),
##                  diag = list(continuous = m1density),
##                  lower = list(continuous = myEllips)
##                  ) +
##         theme_Publication() +
##         theme(legend.position="top")


##     return(p)
## }



m1density <- function(data, mapping, thetaTrue, lims = NULL, ...) {
    g <- ggplot2::ggplot(data = data, mapping = mapping) +
        ggplot2::geom_density(adjust = 2, lwd = 0.25)

    if(!is.null(lims)){
        x <- lims[,colnames(lims) == mapping[[1]]]
        g <- g + xlim(x)
    }


    if(!is.null(thetaTrue))
        g <- g + geom_vline(xintercept = thetaTrue[as.character(mapping[[1]])], linetype = 2, lwd = 0.25)


    g <- g + scale_color_manual(values=c("#FFC34B", "#619CFF", "#f8766dff", "#00BA38"))

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
        scale_shape(solid=F)


    dSLAM <- data[which(data$method == "SLAM"),]
    dSLAM.filt <- data[which(data$method == "SLAM / filtered"),]
    dABC <- data[-grep("SLAM",data$method),]




    g <- g + geom_point(data=dABC, mapping = aes(shape = method), alpha = 0.2, size = 0.5)
    g <- g + geom_point(data=dSLAM.filt, mapping = aes(shape = method), alpha = 0.6, size = 0.5)
    g <- g + geom_point(data=dSLAM, mapping = aes(shape = method), alpha = 0.6, size = 0.5)

        ## stat_density2d(geom="density_2d",
        ##                aes(color = method, alpha=10*..level..),
        ##                #size=2,
        ##                contour=TRUE,
        ##                lwd = 0.25)
    g <- g +  stat_ellipse(level = 0.95, lwd = 0.25)



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

    g <- g + scale_color_manual(values=c("#FFC34B", "#619CFF", "#f8766dff", "#00BA38"))
    g <- g + scale_shape_manual(values = c(2:5))#,3,4))


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

plotMultiDens <- function(names, rows = seq(100,10000,100), thetaTrue = NULL, selectedParams = NULL){



    post <- data.frame()
    for(name in names){
        ## load data

        if(name == "sfull") {
            load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/publish/beta2/sfull.RData")
            sfull <- getTop(sfull, rows, FALSE, "SLAM")
            post <- rbind(post, sfull)
        } else if(name == "sbin") {
            load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/publish/beta2/sbin.RData")
            sbin <- getTop(sbin, rows, FALSE, "SLAM / filtered")
            post <- rbind(post, sbin)
        } else if(name == "abcfull") {
            load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/toPublish/abcfull.RData")
            abcfull <- getTop(abcfull, rows, TRUE, "ABC")
            post <- rbind(post, abcfull)
        } else if(name == "abcbin") {
            load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/toPublish/abcbin.RData")
            abcbin <- getTop(abcbin, rows, TRUE, "ABC / filtered")
            post <- rbind(post, abcbin)
        } else {
            print("name not accounted for in code")
            stopifnot(FALSE)
        }
    }

    ##post <- rbind(sl1, sl2, sl3)




    numParam <- dim(post)[2]-1
    ## ## which columns and columnlabels ...
    ## if(numParam == 4) { ## 2 beta
    ##     columns = c("upsilon", "beta_t1", "beta_t2", "gamma")
    ##     columnLabels = c("upsilon", "beta[1]", "beta[2]", "gamma")
    ## } else if(numParam == 5) { ## 3 beta
    ##     columns = c("upsilon", "beta_t1", "beta_t2", "beta_t3", "gamma")
    ##     columnLabels = c("upsilon", "beta[1]", "beta[2]", "beta[3]", "gamma")
    ## } else if(numParam == 6) { ## 3 beta & prev
    ##     columns = c("upsilon", "beta_t1", "beta_t2", "beta_t3", "gamma", "prev")
    ##     columnLabels = c("upsilon", "beta[1]", "beta[2]", "beta[3]", "gamma", "p[0]")
    ## } else {stopifnot(FALSE)}

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

    limits <- rbind(as.numeric(apply(post,2,min)[1:numParam])*0.85,
                    as.numeric(apply(post,2,max)[1:numParam])*1.15)

    colnames(limits) <- columns

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




plotDensitites <- function(){
    ## load and re-store the data.
    ## load("~/Gits/BPD/R/INFERENCE/old/1600Full/output/abInference.RData")
    ## abc1 <- a
    load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcfull.RData")
    abc1 <- abcfull
    load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcbin.RData")
    abc2 <- abcbin
    load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/publish/beta2/sfull.RData")
    sl1 <- sfull
    load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/publish/beta2/sbin.RData")
    sl2 <- sbin







    posteriorPlot <- function(infe, n = 1, N = NULL, by = NULL, ordered = FALSE, limits = NULL){
        data <- infe$getPosterior()
        dims <- dim(data)

        if(is.null(N))
            N <- dims[1]

        if(ordered){
            data <- data[order(data[,dims[2]]),]
        }

        if(!is.null(by))
            data <- data[seq(n,N,by),]
        else
            data <- data[n:N,]

        my2dDensity <- function(data, mapping, lims = NULL,  ...){
            g <- ggplot2::ggplot(data = data, mapping = mapping) +
                ggplot2::stat_density_2d(geom = "contour", bins = 5,
                                         ggplot2::aes(color = ..level..))
            if(!is.null(lims)){
                x <- lims[,colnames(lims) == mapping[[1]]]
                y <- lims[,colnames(lims) == mapping[[2]]]

                g <- g + xlim(x) + ylim(y)
            }
            return(g)
        }

        my1dDensity <- function(data, mapping, lims = NULL, ...) {
            g <- ggplot2::ggplot(data = data, mapping = mapping) +
                ggplot2::geom_density()

            if(!is.null(lims)){
                x <- lims[,colnames(lims) == mapping[[1]]]
                g <- g + xlim(x)
            }
            return(g)

        }

        myPoints <- function(data, mapping, lims = NULL, ...){
            g <- ggplot2::ggplot(data = data, mapping = mapping) +
                ggplot2::geom_point()

            if(!is.null(lims)){
                x <- lims[,colnames(lims) == mapping[[1]]]
                y <- lims[,colnames(lims) == mapping[[2]]]

                g <- g + xlim(x) + ylim(y)
            }
            return(g)
        }



        g <- GGally::ggpairs(data = data[,seq(1,dims[2]-1)],
                             columnLabels = c("upsilon", "beta[1]", "beta[2]", "gamma"),
                             labeller = "label_parsed",
                             upper = list(continuous = GGally::wrap(my2dDensity, lims = limits)),
                             diag = list(continuous = GGally::wrap(my1dDensity, lims = limits)),
                             lower = list(continuous = GGally::wrap(myPoints, lims = limits))
                             ) +
            theme_Publication()
        g
    }

    mins <- getExtrema(list(abc1, abc2, sl1, sl2), min, 0.9)
    maxs <- getExtrema(list(abc1, abc2, sl1, sl2), max, 1.1)
    limits <- rbind(mins, maxs)

    p1 <- posteriorPlot(abc1, 1, 500, 1, TRUE, limits = limits)
    p2 <- posteriorPlot(abc2, 1, 500, 1, TRUE, limits = limits)
    p3 <- posteriorPlot(sl1, 500, 10000, 10, FALSE, limits = limits)
    p4 <- posteriorPlot(sl2, 500, 10000, 10, FALSE, limits = limits)

    ##return(list(p1, p2, p3, p4))
    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

getExtrema <- function(listOfInfe, func = min, extra = 1){
    data <- do.call("rbind",
                    lapply(listOfInfe, function(x){x$getPosterior()[,1:4]})
                    )
    extrema <- apply(data, 2, func)
    extrema <- extrema * extra
    return(extrema)
}

plotEllips <- function(N = 250, level = 0.95){
    ## load and restoredata
    load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcfull.RData")
    abc1 <- abcfull
    load("~/Gits/BPD/R/INFERENCE/1600system//output/ABC/abcbin.RData")
    abc2 <- abcbin
    load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/sfull.RData")
    sl1 <- sfull
    load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/sbin.RData")
    sl2 <- sbin

    ## extract the posteriors
    ## full state abc
    d.abc1 <- abc1$getPosterior()
    d.abc1 <- d.abc1[order(d.abc1[,5])[1:1000],c("upsilon", "gamma")]
    n.abc1 <- sample(1:dim(d.abc1)[1], N, replace = TRUE)
    post.abc1 <- d.abc1[n.abc1,]
    post.abc1$type <- rep("abc-full", N)

    ## binary state abc
    d.abc2 <- abc2$getPosterior()
    d.abc2 <- d.abc2[order(d.abc2[,5])[1:1000],c("upsilon", "gamma")]
    n.abc2 <- sample(1:dim(d.abc2)[1], N, replace = TRUE)
    post.abc2 <- d.abc2[n.abc2,]
    post.abc2$type <- rep("abc-binary", N)


    ## full state sl
    d.sl1 <- sl1$getPosterior()[seq(500,3000,1),c("upsilon", "gamma")]
    n.sl1 <- sample(1:dim(d.sl1)[1], N, replace = TRUE)
    post.sl1 <- d.sl1[n.sl1,]
    post.sl1$type <- rep("sl-full", N)

    ## binary state sl
    d.sl2 <- sl2$getPosterior()[seq(500,3000,1),c("upsilon", "gamma")]
    n.sl2 <- sample(1:dim(d.sl2)[1], N, replace = TRUE)
    post.sl2 <- d.sl2[n.sl2,]
    post.sl2$type <- rep("sl-binary", N)

    ## true value
    true.df <- data.frame(upsilon = rep(0.0075,N), gamma = rep(0.1,N), type = rep("actual",N))

    data <- data.frame(upsilon = c(post.abc1$upsilon, post.abc2$upsilon,
                                   post.sl1$upsilon, post.sl2$upsilon),
                       gamma = c(post.abc1$gamma, post.abc2$gamma,
                                 post.sl1$gamma, post.sl2$gamma),
                       type = c(post.abc1$type, post.abc2$type,
                                post.sl1$type, post.sl2$type)
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
        geom_point(aes(x = 0.0075, y = 0.1), shape = 13, size = 6)

    return(p)
}

abstractPlot <- function(){
    ## retrieve abc (bin) posterior
    load("~/Gits/BPD/R/INFERENCE/1600system/output/ABC/abcbin.RData")
    abc <- abcbin
    d.abc <- abc$getPosterior()
    post.abc <- d.abc[order(d.abc[,5])[1:250],1:4]
    dens.abc <- apply(post.abc, 2, density)

    ## retrieve sl (bin) posterior
    load("~/Gits/BPD/R/INFERENCE/1600system/output/slam/sbin.RData")
    sl <- sbin
    d.sl <- sl$getPosterior()
    post.sl <- d.sl[seq(500,3000,10),1:4]
    dens.sl <- apply(post.sl, 2, density)

    ## get limits
    maxx.sl <- sapply(dens.sl, function(val){max(val$x)})
    maxx.abc <- sapply(dens.abc, function(val){max(val$x)})
    maxx <- apply(do.call("rbind", list(maxx.sl, maxx.abc)), 2, max)

    maxy.sl <- sapply(dens.sl, function(val){max(val$y)})
    maxy.abc <- sapply(dens.abc, function(val){max(val$y)})
    maxy <- apply(do.call("rbind", list(maxy.sl, maxy.abc)), 2, max)

    minx.sl <- sapply(dens.sl, function(val){min(val$x)})
    minx.abc <- sapply(dens.abc, function(val){min(val$x)})
    minx <- apply(do.call("rbind", list(minx.sl, minx.abc)), 2, min)

    miny.sl <- sapply(dens.sl, function(val){min(val$y)})
    miny.abc <- sapply(dens.abc, function(val){min(val$y)})
    miny <- apply(do.call("rbind", list(miny.sl, miny.abc)), 2, min)



    ## plot
    par(oma = c(3, 2, 1, 0))
    par(mfrow = c(2,2), mar = c(3.9, 1.1, 2.1, 1.1))


    ## upsilon
    variable <- "upsilon"
    plot(x = c(minx[variable], maxx[variable]), y = c(miny[variable], maxy[variable]), type = "n",
         xlab = expression(upsilon), ylab = "Density")
    polygon(dens.sl$upsilon, col=rgb(0, 0, 1, 0.5))
    polygon(dens.abc$upsilon, col=rgb(1, 0, 0, 0.5))
    abline(v=0.0075, col = rgb(1,0,1,1), lwd = 3, lty = 2)


    ## beta_t1
    variable <- "beta_t1"
    plot(x = c(minx[variable], maxx[variable]), y = c(miny[variable], maxy[variable]), type = "n",
         xlab = expression(beta[t[1]]), ylab = "Density")
    polygon(dens.sl$beta_t1, col=rgb(0, 0, 1, 0.5))
    polygon(dens.abc$beta_t1, col=rgb(1, 0, 0, 0.5))
    abline(v=0.05, col = rgb(1,0,1,1), lwd = 3, lty = 2)


    ## beta_t2
    variable <- "beta_t2"
    plot(x = c(minx[variable], maxx[variable]), y = c(miny[variable], maxy[variable]), type = "n",
         xlab = expression(beta[t[2]]), ylab = "Density")
    polygon(dens.sl$beta_t2, col=rgb(0, 0, 1, 0.5))
    polygon(dens.abc$beta_t2, col=rgb(1, 0, 0, 0.5))
    abline(v=0.085, col = rgb(1,0,1,1), lwd = 3, lty = 2)


    ## gamma
    variable <- "gamma"
    plot(x = c(minx[variable], maxx[variable]), y = c(miny[variable], maxy[variable]), type = "n",
         xlab = expression(gamma), ylab = "Density")
    polygon(dens.sl$gamma, col=rgb(0, 0, 1, 0.5))
    polygon(dens.abc$gamma, col=rgb(1, 0, 0, 0.5))
    abline(v=0.1, col = rgb(1,0,1,1), lwd = 3, lty = 2)


    par(fig = c(0, 1, 0, 1), oma = c(0, 3.5, 3, 0), mar = c(0, 0, 0, 0), new = TRUE)

    legend("center",legend = c("ABC", "SLMCMC", "True"), col = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5), rgb(1,0,1,1)),
           horiz = TRUE, lty = c(1,1,2), lwd = c(5,5,3), xpd = TRUE, bty = "n")


}

gridPlot <- function(data = NULL, logy = FALSE, save = FALSE) {
    if(is.null(data)) {
        ## load data into post
        load("~/Gits/BPD/R/INFERENCE/1600system/output/slgrid/sgrid10.RData")
        post <- sgrid$getPosterior()
    }
    else
        post <- data$getPosterior()

    ## if(expIt) {
    ##     post$dist <- exp(post$dist)
    ##     ylab = expression(-1 %*% "SL")
    ## } else {
    ##     ylab = expression(-1 %*% "log SL")
    ## }
    ## extract each dim analysis
    dims <- dim(post)
    thetaLen <- 4
    distLen <- dims[2]-thetaLen
    postLen <- dims[1]
    dataSec <- seq(1,dims[1]+1, by = postLen/thetaLen)




    #upsilon <- post[dataSec[1]:(dataSec[2]-1), c(1, (thetaLen+1):distLen)]

    ##* UPSILON
    upsilon <- data.frame(value = rep(NA, postLen/thetaLen))
    upsilon$value <- post[dataSec[1]:(dataSec[2]-1),"upsilon"]
    ## normalize
    upsilonAll <- -1*post[dataSec[1]:(dataSec[2]-1), (thetaLen+1):distLen]
    upsilonAll.norm <- (upsilonAll - min(upsilonAll))/ (max(upsilonAll) - min(upsilonAll))
    ##
    upsilon$mean <- apply(upsilonAll.norm,1, mean)
    upsilon$sd <- apply(upsilonAll.norm, 1, sd)
    ## within the area!
    upsilon$category <- rep("outside", postLen/thetaLen)
    category <- which( (upsilon$mean-2*upsilon$sd) <= (upsilon$mean[11] + 2*upsilon$sd[11]))
    upsilon$category[category] <- rep("inside", length(category))

    ##* BETA_T1
    beta_t1 <- data.frame(value = rep(NA, postLen/thetaLen))
    beta_t1$value <- post[dataSec[2]:(dataSec[3]-1),"beta_t1"]
    ## normalize
    beta_t1All <- -1*post[dataSec[2]:(dataSec[3]-1), (thetaLen+1):distLen]
    beta_t1All.norm <- (beta_t1All - min(beta_t1All))/ (max(beta_t1All) - min(beta_t1All))
    ##
    beta_t1$mean <- apply(beta_t1All.norm,1, mean)
    beta_t1$sd <- apply(beta_t1All.norm, 1, sd)
    ## within the area!
    beta_t1$category <- rep("outside", postLen/thetaLen)
    category <- which( (beta_t1$mean-2*beta_t1$sd) <= (beta_t1$mean[11] + 2*beta_t1$sd[11]))
    beta_t1$category[category] <- rep("inside", length(category))

    ##* BETA_T2
    beta_t2 <- data.frame(value = rep(NA, postLen/thetaLen))
    beta_t2$value <- post[dataSec[3]:(dataSec[4]-1), "beta_t2"]
    ## normalize
    beta_t2All <- -1*post[dataSec[3]:(dataSec[4]-1), (thetaLen+1):distLen]
    beta_t2All.norm <- (beta_t2All - min(beta_t2All))/ (max(beta_t2All) - min(beta_t2All))
    ##
    beta_t2$mean <- apply(beta_t2All.norm,1, mean)
    beta_t2$sd <- apply(beta_t2All.norm, 1, sd)
    ## within the area!
    beta_t2$category <- rep("outside", postLen/thetaLen)
    category <- which( (beta_t2$mean-2*beta_t2$sd) <= (beta_t2$mean[11] + 2*beta_t2$sd[11]))
    beta_t2$category[category] <- rep("inside", length(category))

    ##* GAMMA
    gamma <- data.frame(value = rep(NA, postLen/thetaLen))
    gamma$value <- post[dataSec[4]:(dataSec[5]-1),"gamma"]
    ## normalize
    gammaAll <- -1*post[dataSec[4]:(dataSec[5]-1), (thetaLen+1):distLen]
    gammaAll.norm <- (gammaAll - min(gammaAll))/ (max(gammaAll) - min(gammaAll))
    ##
    gamma$mean <- apply(gammaAll.norm,1, mean)
    gamma$sd <- apply(gammaAll.norm, 1, sd)
    ## within the area!
    gamma$category <- rep("outside", postLen/thetaLen)
    category <- which( (gamma$mean-2*gamma$sd) <= (gamma$mean[11] + 2*gamma$sd[11]))
    gamma$category[category] <- rep("inside", length(category))


     if(save) {
         dirname <- "~/Gits/BPD/PLOTS/toPublish/1600/"
         if(logy)
             filename <- paste(dirname, "gridSemiLog.pdf", sep="")
         else
             filename <- paste(dirname, "grid.pdf", sep="")
         ##pdf(filename)
     }


    ylab = expression(-1 %*% "log SL")

    if(save)
        base_size <- 6
    else
        base_size <- 14

    p1 <- ggplot(upsilon, aes(x=value, y=mean)) +
        geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), fill="grey80")+
        geom_line() +
        geom_vline(xintercept = mean(upsilon$value), linetype = 2) +
        geom_point(aes(color = category)) +
        ##geom_pointrange(aes(ymin=mean-2*sd, ymax=mean+2*sd)) +
        ##geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), alpha = 0.3) +
        theme_Publication(base_size) +
        labs(x = expression(upsilon), y = ylab)


    p2 <- ggplot(beta_t1, aes(x=value, y=mean)) +
        geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), fill="grey80")+
        geom_line() +
        geom_vline(xintercept = mean(beta_t1$value), linetype = 2) +
        geom_point(aes(color = category)) +
        #geom_pointrange(aes(ymin=mean-2*sd, ymax=mean+2*sd)) +
        ##geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), alpha = 0.3) +
        theme_Publication(base_size) +
        labs(x = expression(beta[1]), y = ylab)


    p3 <- ggplot(beta_t2, aes(x=value, y=mean)) +
        geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), fill="grey80")+
        geom_line() +
        geom_point(aes(color = category)) +
        geom_vline(xintercept = mean(beta_t2$value), linetype = 2) +
        #geom_pointrange(aes(ymin=mean-2*sd, ymax=mean+2*sd)) +
        theme_Publication(base_size) +
        labs(x = expression(beta[2]), y = ylab)



    p4 <- ggplot(gamma, aes(x=value, y=mean)) +
        geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), fill="grey80")+
        geom_line() +
        geom_vline(xintercept = mean(gamma$value), linetype = 2) +
        geom_point(aes(color = category)) +
        ## geom_pointrange(aes(ymin=mean-2*sd, ymax=mean+2*sd)) +
        theme_Publication(base_size) +
        labs(x = expression(gamma), y = ylab)

    if(save) {
        ## margin fix
        p1 <- p1 + theme(plot.margin=grid::unit(c(1,1,1,1),"mm"),
                         legend.margin=margin(l = -0.5, unit='cm'),
                         legend.title =element_blank())

        p2 <- p2 + theme(plot.margin=grid::unit(c(1,1,1,1),"mm"),
                         legend.margin=margin(l = -0.5, unit='cm'))

        p3 <- p3 + theme(plot.margin=grid::unit(c(1,1,1,1),"mm"),
                         legend.margin=margin(l = -0.5, unit='cm'))

        p4 <- p4 + theme(plot.margin=grid::unit(c(1,1,1,1),"mm"),
                         legend.margin=margin(l = -0.5, unit='cm'))

    }



    #p <- grid.arrange(p1,p2,p3,p4, ncol = 2)

    p <- ggarrange(plotlist = list(p1, p2, p3, p4),
                   common.legend = TRUE,
                   legend = "top")

    if(save)
        ggsave(filename = filename, p, width = 13.05, height = 8.7, units = "cm")

    return(p)


}

main <- function(names = c("abcfull", "abcbin", "sfull", "sbin"), rows = seq(100,10000,100),
                  save = FALSE, selectedParams = NULL){ ## c("upsilon","beta_t2","gamma")){
    thetaTrue = c(upsilon = 0.0075,
                  beta_t1 = 0.05,
                  beta_t2 = 0.085,
                  gamma = 0.1)

    p <- plotMultiDens(names, rows, thetaTrue, selectedParams)
    if(save) {
        p <- p + theme_Publication(6) +
            theme(plot.margin=grid::unit(c(0,0,0,0),"mm"),
                  legend.margin=margin(t = 0.1,b=-0.25, unit='cm'),
                  axis.line = element_line(size = 0.25),
                  axis.ticks = element_line(colour = "black", size = 0.25),
                  panel.grid.major = ggplot2::element_line(colour="#f0f0f0", size = 0.25),
                      legend.position = "top")
        dirname <- "~/Gits/BPD/PLOTS/toPublish/1600/"
        ggsave(paste(dirname, "1600MultiPosterior.pdf", sep = ""), p, width = 8.7, height = 8.7, units = "cm")
    }



    return(p)
}

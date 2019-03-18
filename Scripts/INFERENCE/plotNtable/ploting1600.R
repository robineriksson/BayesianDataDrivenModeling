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

gridPlot <- function(data = NULL, logy = FALSE, save = FALSE) {
    if(is.null(data)) {
        ## load data into post
        load("~/Gits/BPD/R/INFERENCE/1600system/output/slgrid/sgrid10.RData")
        post <- sgrid$getPosterior()
    }
    else
        post <- data$getPosterior()

    dims <- dim(post)
    thetaLen <- 4
    distLen <- dims[2]-thetaLen
    postLen <- dims[1]
    dataSec <- seq(1,dims[1]+1, by = postLen/thetaLen)


    ##* UPSILON
    upsilon <- data.frame(value = rep(NA, postLen/thetaLen))
    upsilon$value <- post[dataSec[1]:(dataSec[2]-1),"upsilon"]
    ## normalize
    upsilonAll <- -1*post[dataSec[1]:(dataSec[2]-1), (thetaLen+1):dims[2]]
    upsilonAll.norm <- (upsilonAll - min(upsilonAll))/ (max(upsilonAll) - min(upsilonAll))
    ##
    upsilon$mean <- apply(upsilonAll.norm,1, mean)
    upsilon$sd <- apply(upsilonAll.norm, 1, sd)
    ## within the area!
    upsilon$category <- rep("outside", postLen/thetaLen)
    smallest <- which(upsilon$mean == min(upsilon$mean))
    category <- which( (upsilon$mean-2*upsilon$sd) <= (upsilon$mean[smallest] + 2*upsilon$sd[smallest]))
    upsilon$category[category] <- rep("inside", length(category))

    ##* BETA_T1
    beta_t1 <- data.frame(value = rep(NA, postLen/thetaLen))
    beta_t1$value <- post[dataSec[2]:(dataSec[3]-1),"beta_t1"]
    ## normalize
    beta_t1All <- -1*post[dataSec[2]:(dataSec[3]-1), (thetaLen+1):dims[2]]
    beta_t1All.norm <- (beta_t1All - min(beta_t1All))/ (max(beta_t1All) - min(beta_t1All))
    ##
    beta_t1$mean <- apply(beta_t1All.norm,1, mean)
    beta_t1$sd <- apply(beta_t1All.norm, 1, sd)
    ## within the area!
    beta_t1$category <- rep("outside", postLen/thetaLen)
    smallest <- which(beta_t1$mean == min(beta_t1$mean))
    category <- which( (beta_t1$mean-2*beta_t1$sd) <= (beta_t1$mean[smallest] + 2*beta_t1$sd[smallest]))
    beta_t1$category[category] <- rep("inside", length(category))

    ##* BETA_T2
    beta_t2 <- data.frame(value = rep(NA, postLen/thetaLen))
    beta_t2$value <- post[dataSec[3]:(dataSec[4]-1), "beta_t2"]
    ## normalize
    beta_t2All <- -1*post[dataSec[3]:(dataSec[4]-1), (thetaLen+1):dims[2]]
    beta_t2All.norm <- (beta_t2All - min(beta_t2All))/ (max(beta_t2All) - min(beta_t2All))
    ##
    beta_t2$mean <- apply(beta_t2All.norm,1, mean)
    beta_t2$sd <- apply(beta_t2All.norm, 1, sd)
    ## within the area!
    beta_t2$category <- rep("outside", postLen/thetaLen)
    smallest <- which(beta_t2$mean == min(beta_t2$mean))
    category <- which( (beta_t2$mean-2*beta_t2$sd) <= (beta_t2$mean[smallest] + 2*beta_t2$sd[smallest]))
    beta_t2$category[category] <- rep("inside", length(category))

    ##* GAMMA
    gamma <- data.frame(value = rep(NA, postLen/thetaLen))
    gamma$value <- post[dataSec[4]:(dataSec[5]-1),"gamma"]
    ## normalize
    gammaAll <- -1*post[dataSec[4]:(dataSec[5]-1), (thetaLen+1):dims[2]]
    gammaAll.norm <- (gammaAll - min(gammaAll))/ (max(gammaAll) - min(gammaAll))
    ##
    gamma$mean <- apply(gammaAll.norm,1, mean)
    gamma$sd <- apply(gammaAll.norm, 1, sd)
    ## within the area!
    gamma$category <- rep("outside", postLen/thetaLen)
    smallest <- which(gamma$mean == min(gamma$mean))
    category <- which( (gamma$mean-2*gamma$sd) <= (gamma$mean[smallest] + 2*gamma$sd[smallest]))
    gamma$category[category] <- rep("inside", length(category))


     if(save) {
         dirname <- "~/Gits/BPD/PLOTS/toPublish/1600/"
         if(logy)
             filename <- paste(dirname, "gridSemiLog.pdf", sep="")
         else
             filename <- paste(dirname, "grid.pdf", sep="")
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
        theme_Publication(base_size) +
        labs(x = expression(upsilon), y = ylab)


    p2 <- ggplot(beta_t1, aes(x=value, y=mean)) +
        geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), fill="grey80")+
        geom_line() +
        geom_vline(xintercept = mean(beta_t1$value), linetype = 2) +
        geom_point(aes(color = category)) +
        theme_Publication(base_size) +
        labs(x = expression(beta[1]), y = ylab)


    p3 <- ggplot(beta_t2, aes(x=value, y=mean)) +
        geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), fill="grey80")+
        geom_line() +
        geom_point(aes(color = category)) +
        geom_vline(xintercept = mean(beta_t2$value), linetype = 2) +
        theme_Publication(base_size) +
        labs(x = expression(beta[2]), y = ylab)



    p4 <- ggplot(gamma, aes(x=value, y=mean)) +
        geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), fill="grey80")+
        geom_line() +
        geom_vline(xintercept = mean(gamma$value), linetype = 2) +
        geom_point(aes(color = category)) +
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

library(SimInf)
library(SimInfInference)
library(ggplot2)

##' Run the simulation study
##' return table with the population prevalence
##' includes Credible interval.
##'
##' @param N number of posterior samples
##' @param nodes the number of nodes looked at.
##' @param CI how many SD to include in the CI
##' @param extend how many days to relax the mathc with.
##' @param PosteriorFilePath the filepath to the computed posterior.
##' @param dataDir the path to the directory holding the data
main <- function(N = 100, nodes = 250, CI = 2, extend = 1,
                 PosteriorFilePath, dataDir) {
    r <- runSimulation(nodes = nodes, N = N, extend = extend,
                       filename = PosteriorFilePath, dataDir = dataDir)

    qq <- apply(prev, 2, function(x){findQuantile(x,0.05)})
    prev <- r$prev
    means <- apply(prev, 2, mean)
    sds <- apply(prev, 2, sd)

    df <- data.frame(mean = means, sd = sds, qq.low = qq["low",], qq.up = qq["up",])

    credible <- apply(df,1,function(x,CI){x[1] + CI*c(-x[2], x[2])}, CI = CI)
    df$low <- credible[1,]
    df$high <- credible[2,]

    rownames(df) <- colnames(prev)
    return(list(df = df, r = r))
}

##' Find the Credible Interval (CI)
##' @param x the data to look at
##' @param alpha Classic CI parameter.
findQuantile <- function(x,alpha = 0.05) {
    xs <- sort(x)
    cpxx <- cumsum(xs) / sum(xs)
    upper = 1-alpha/2
    lower = 1-upper
    low <- xs[which(cpxx >= lower)[1]]   # lower boundary
    up <- xs[which(cpxx >= upper)[1]-1] # upper boundary

    return(c("low"=low,"up"=up))
}


loadPosterior <- function(filename, N, by = 100) {
    ## load posterior
    load(filename)
    posterior.All <- mobs$getPosterior()
    dims <- dim(posterior.All)
    rows <- seq(1,dims[1],by)
    cols <- seq_len(dims[2]-1)
    posterior.thinned <- posterior.All[rows,cols]

    ## if we do 100 draws and do not use replace,
    ## then we're using all the accepted parameter pairs
    theta.post <- posterior.thinned[sample(seq_len(length(rows)), N, replace = T),]
    return(theta.post)
}


##' Compute the population prevalence for the number of nodes in input.
##' Looks at two different intervals, mentioned in the paper.
##'
##' @param nodes the nodes to observe
##' @param N number of posterior draws
##' @param extend if we look outside of the exact transfer day date.
##' @param filename path to filename
##' @param dataDir path to data directory
runSimulation <- function(nodes = 250, N = 100, extend = NULL,
                          filename, dataDir {
    set.seed(0)

    ## load posterior
    theta.post <- loadPosterior(filename, N)

    ## define tspan vector
    sise_path <- paste(dataDir, "SISe.rda", sep = "")
    load(sise_path) ## load model
    tspan0 <- seq(head(model@events@time,1), tail(model@events@time,1), 1)

    ## 0. Define date A and B, find all nodes that have exit events on both dates.
    datetime <- as.Date(tspan0, origin = "2004-12-31")

    ## For what timespan should we track the temporal mean?
    ## dateOI = dates of interest
    dateOI <- c("2008-10-01","2009-01-01","2009-04-01","2009-07-01")
    dateOI <- c(dateOI, "2011-10-01","2012-01-01","2012-04-01","2012-07-01")
    dateOI.elem <- c()
    for(d in dateOI)
        dateOI.elem <- c(dateOI.elem, which(datetime == d))
    dateOI.num <- tspan0[dateOI.elem]

    dateB <- "2012-11-01"
    dateB.elem <- which(datetime == dateB)
    dateB.num <- tspan0[dateB.elem]

    tspan.record <- tspan0[dateOI.elem]

    tspan <- tspan0[1:dateOI.elem[length(dateOI.elem)]]

    ## find the noteset
    events <- model@events
    events.df <- data.frame(event = events@event,
                            time = events@time,
                            node = events@node,
                            n = events@n)

    if(is.null(extend))
        eventsOI <- events.df[which(events.df$time %in% dateOI.num),]
    else {
        extension <- c()
        for(e in 1:extend)
            extension <- c(extension, dateOI.num - e, dateOI.num + e)
        eventsOI <- events.df[which(events.df$time %in% extension),]
    }

    eventsOI.exit <- eventsOI[which(eventsOI$event == 0),]


    print("Run Simulations!")
    prev.out <- matrix(NA, nrow = N, ncol = 2)
    colnames(prev.out) <- c("08-09", "11-12")

    det.mean <- numeric(N)

    pb <- pbapply::timerProgressBar(width = 50)
    for(n in seq_len(N)) {
        ## extract theta(n)
        ## + class massaging
        t <- as.numeric(theta.post[n,])
        names(t) <- names(theta.post)
        theta <- setTheta(t)

        ## reload model
        load(sise_path)

        ## initialize model
        model <- init_model_widgren(model = model, theta = theta,
                                    tspan = tspan,
                                    phi = "local")

        res <- SimInf::run(model, threads = NULL, solver = "ssm")


        ##* Sample the nodes ****

        ## what nodes to sample?
        nodes.selected <- sampleNodes(eventsOI.exit, nodes, extend)


        ## within-node-prevalence for all the selected nodes for all time
        prev.all <- SimInf::prevalence(res, I ~ S + I, "wnp", node = unlist(c(nodes.selected)))

        ## wnp for selected times
        prev.sel <- prev.all[which(prev.all$time %in%  tspan.record),]
        ## remove NaN, replace with 0
        prev.sel$prevalence[which(is.na(prev.sel$prevalence))] = 0


        ## add column which indicates what yearly average to include in
        prev.sel$year <- rep(NA, dim(prev.sel)[1])
        for(i in 1:dim(nodes.selected)[2]) {
            prev.sel$year[which(prev.sel$node %in% nodes.selected[,i])] = i
        }

        ## add column that states how many samples we should draw from each node
        prev.sel$n <- rep(NA, dim(prev.sel)[1])
        for(i in 1:dim(nodes.selected)[2]) {
            tab <- table(nodes.selected[,i])

            for(j in 1:length(tab)) {
                prev.sel$n[which(prev.sel$node == names(tab[j]))] = tab[[j]]
            }
        }

        ## get the size of the nodes
        prev.sel <- getNodeSize(res, prev.sel, tspan.record, nodes.selected)

        ## swab test the nodes!
        prev.sel$swab <- swabNodes(prev.sel)

        ## evaluate the swabs
        prev.out[n,] <- evalSwab(prev.sel)

        ## how many nodes where detected?
        det.mean[n] <- numberOfInfNodes(prev.sel, tspan.record)
        pbapply::setTimerProgressBar(pb, n/N)
    }
    close(pb)

    return(list(prev = prev.out, tspan.record = tspan.record, det.mean = det.mean))
}

##' Get the size of the population in the node
getNodeSize <- function(res, prev.sel, tspan.record, nodes.selected) {
    ## ## size of the population at the node and time
    u <- SimInf::trajectory(res, c("S","I"),as.is=T, node = unlist(c(nodes.selected)))

    ## selected dates
    u.sel.dates <- u[,as.character(tspan.record)]
    ## from S and I to Sigma = S+I
    i1 <- seq(1,dim(u)[1],2)
    u.sel <- u.sel.dates[i1,] + u.sel.dates[i1+1,]
    rownames(u.sel) <- unique(unlist(c(nodes.selected)))

    prev.sel$size <- rep(NA, dim(prev.sel)[1])
    nodenames <- unique(unlist(c(nodes.selected)))
    for(m in 1:length(nodenames)) {
        node <- nodenames[m]
        prev.sel$size[which(prev.sel$node == node)] = u.sel[m,]
    }

    return(prev.sel)
}

##' Find and randomly select node that are active
##' in the given timeframe.
sampleNodes <- function(eventList, N, extend) {
    ## for each yearset, select N nodes!
    if(is.null(extend))
        extend <- 0

    udates <- unique(eventList$time)
    nodeset1 <- eventList$node[which(eventList$time >= (udates[1] - extend) &
                                    eventList$time <= (udates[length(udates)/2] + extend) )]

    nodeset2 <- eventList$node[which(eventList$time >= (udates[length(udates)/2 + 1] - extend) &
                                     eventList$time <= (udates[length(udates)] + extend))]


    nodes <- data.frame(n1 = sample(nodeset1, N, replace = T),
                        n2 = sample(nodeset2, N, replace = T))

    return(nodes)
}

##' Artifically swab (test for infection) the nodes
swabNodes <- function(nodeData) {
    swab <- apply(nodeData, 1, function(x){
        sum(rbinom(n = x["n"], size = 1, prob = x["prevalence"]))
    })
    return(swab)
}

##' Determine the number of "Inf" in the node.
numberOfInfNodes <- function(nodeData, time) {
    years = unique(nodeData$year)
    detected <- numeric(length(years)*4)
    for(u in years) {
        atyear <- nodeData[nodeData$year == u,]
        for(q in 1:4) { ## quarter of a year
            tab <- table(atyear[atyear$time == time[q + 4*(u-1)], "swab"])
            if(length(tab) > 1) {
                ## if multiple  findings are made at the same node,
                ## this should still count as 1, as they're from the same node.
                detected[q + 4*(u-1)] <- sum(tab[-1])/sum(tab)
            } else
                detected[q + 4*(u-1)] <- 0
        }
    }
    det.mean <- mean(detected)
    return(det.mean)
}

##' Evaluate the swab test. Does is test positive or negative?
evalSwab <- function(nodeData) {
    years <- unique(nodeData$year)
    evaled <- c()
    for(y in years) {
        swabTab <- table(nodeData$swab[nodeData$year == y])

        if(length(swabTab) > 1) {
            neg <- swabTab[[1]]
            pos <- 0
            for(i in 2:length(swabTab))
                pos <- pos + ((i-1)*swabTab[[i]])
            evaled <- c(evaled, pos/(neg + pos))
        } else ## all are negative
            evaled <- c(evaled, 0)
    }

    return(evaled)
}

##' Set the size of the node.
setsize <- function(result, nodes, tspan.record, prev) {
    u <- SimInf::trajectory(result, c("S","I"),as.is=T, node = nodes)
    u.sel <- u[, as.character(tspan.record)]
    ## from S and I to Sigma = S+I
    i1 <- seq(1,dim(u.sel)[1],2)
    u.size <- u.sel[i1,] + u.sel[i1+1,]
    nodenames <- unique(nodes)
    rownames(u.size) <- nodenames

    ## paste to the prevalence dataframe
    prev$size.median <- apply(u.size, 2, median)
    prev$size.mean <- apply(u.size, 2, mean)

    return(prev)
}

##' Node prevalence for interval 2008-2009
##'
##' @param nodes the nodes to observe
##' @param N number of posterior draws
##' @param extend if we look outside of the exact transfer day date.
##' @param filename path to filename
##' @param dataDir path to data directory
trueNOPsim <- function(nodes = 250, N = 100, extend = NULL, filename, dataDir) {
    set.seed(0)

    ## load posterior
    theta.post <- loadPosterior(filename, N)

    ## define tspan vector
    sise_path <- paste(dataDir, "SISe.rda", sep = "")
    load(sise_path) ## load model
    tspan0 <- seq(head(model@events@time,1), tail(model@events@time,1), 1)

    ## 0. Define date A and B, find all nodes that have exit events on both dates.
    datetime <- as.Date(tspan0, origin = "2004-12-31")

    ## For what timespan should we track the temporal mean?
    ## dateOI = dates of interest
    dateOI <- c("2008-10-01","2009-01-01","2009-04-01","2009-07-01")
    ##dateOI <- c(dateOI, "2011-10-01","2012-01-01","2012-04-01","2012-07-01")
    dateOI.elem <- c()
    for(d in dateOI)
        dateOI.elem <- c(dateOI.elem, which(datetime == d))
    dateOI.num <- tspan0[dateOI.elem]

    tspan.record <- tspan0[dateOI.elem]

    tspan <- tspan0[1:dateOI.elem[length(dateOI.elem)]]

    ## find the noteset
    events <- model@events
    events.df <- data.frame(event = events@event,
                            time = events@time,
                            node = events@node,
                            n = events@n)

    if(is.null(extend))
        eventsOI <- events.df[which(events.df$time %in% dateOI.num),]
    else {
        extension <- c()
        for(e in 1:extend)
            extension <- c(extension, dateOI.num - e, dateOI.num + e)
        eventsOI <- events.df[which(events.df$time %in% extension),]
    }

    eventsOI.exit <- eventsOI[which(eventsOI$event == 0),]

    ## what nodes to sample?
    nodes.events <- sampleNodes(eventsOI.exit, nodes, extend)
    nodes.all <- NULL
    nObs_path = paste(dataDir, "nObs.rda", sep ="")
    load(nObs_path)
    nodes.n126 <- as.numeric(nObs)


    print("Run Simulations!")
    prev.out <- matrix(NA, nrow = N, ncol = 2)
    colnames(prev.out) <- c("08-09", "11-12")

    prev <- data.frame(time = NA, prevalence = NA, size.median = NA, size.mean = NA, selection = NA, N = NA)


    pb <- pbapply::timerProgressBar(width = 50)
    for(n in seq_len(N)) {
        ## extract theta(n)
        ## + class massaging
        t <- as.numeric(theta.post[n,])
        names(t) <- names(theta.post)
        theta <- setTheta(t)

        ## reload model
        load(sise_path)

        ## initialize model
        model <- init_model_widgren(model = model, theta = theta,
                                    tspan = tspan,
                                    phi = "local")

        res <- SimInf::run(model, threads = NULL, solver = "ssm")

        ## within-node-prevalence for all the selected nodes for all time
        prev.events <- SimInf::prevalence(res, I ~ S + I, "nop", node = unlist(c(nodes.events)))
        prev.events <- prev.events[which(prev.events$time %in% tspan.record), ]

        prev.all <- SimInf::prevalence(res, I ~ S + I, "nop", node = nodes.all)
        prev.all <- prev.all[which(prev.all$time %in% tspan.record), ]

        prev.n126 <- SimInf::prevalence(res, I ~ S + I, "nop", node = nodes.n126)
        prev.n126 <- prev.n126[which(prev.n126$time %in% tspan.record), ]

        ## get the average node population for the nodeset observed!

        ## ## size of the population at the node and time



        ## prev.sel$size <- rep(NA, dim(prev.sel)[1])
        ## nodenames <- unique(unlist(c(nodes.selected)))
        ## for(m in 1:length(nodenames)) {
        ##     node <- nodenames[m]
        ##     prev.sel$size[which(prev.sel$node == node)] = u.sel[m,]
        ## }


        ## get the size of the nodes
        prev.events <- setsize(res, unlist(c(nodes.events)), tspan.record, prev.events)
        prev.all <- setsize(res, nodes.all, tspan.record, prev.all)
        prev.n126 <- setsize(res, nodes.n126, tspan.record, prev.n126)

        prev.rbind <- rbind(prev.events, prev.all,prev.n126)
        prev.rbind$selection <- rep(c("event", "all", "n126"), each = length(tspan.record))
        prev.rbind$N <- rep(n, dim(prev.rbind)[1])

        prev <- rbind(prev, prev.rbind)
        pbapply::setTimerProgressBar(pb, n/N)
    }
    close(pb)

    prev <- prev[-1,]
    return(prev)
}

plotit <- function(res, size = 14) {

    df <- res$r$prev
    df.r <- data.frame(prevalence = unlist(c(df)))
    df.r$period <- rep(colnames(df), each = dim(df)[1])

    p <- ggplot(df.r, aes(prevalence, fill = period)) +
        geom_density(alpha = 0.7) +
        scale_x_continuous(labels = scales::percent,
                           breaks = seq(0.015, 0.05, 0.0075)) +
        theme_Publication(size)

    return(p)
}

##' Code that runs and check the node prevalence for the 126
##' herds and checks the node prevalence for simulation and observation.
##' @param runs the number of iterations to use for the "posterior"
##' @param filename name of the posterior file
##' @param dataDir the directory holding the observation/model data.
n126nodePrevalence <- function(runs = 1, filename, dataDir) {
    set.seed(0)

    ## draw samples from the posterior
    theta.post <- loadPosterior(filename, runs)

  	## what nodes
    nObs_path <- paste(dataDir, "nObs.rda", sep = "")
    load(nObs_path) ## load: nObs

    ## same for all runs
    sise_path <- paste(dataDir, "SISe.rda", sep = "")
    load(sise_path) ## load model


	obs_path <- paste(dataDir, "obs.rda", sep = "")
    load(obs_path) ## load: obs
    obs$sample <- as.numeric(obs$status)




    ## ** FOR OUTPUT
    tab.obs <- numPos(table(obs$sample))

    obs.qtr <- qtrStore(as.data.frame(obs), "time", "sample")

    tab.obs.qtr <- numPos(table(obs.qtr$sample))
    wm.obs <- meanNodePrev(obs.qtr)
    ## **


    ##p1 <- plotAggregate2(obs, fun = fun)
    ##t <- 1:length(p1)

    tspan0 <- seq(head(model@events@time,1), tail(model@events@time,1), 1)
    realDates <- zoo::as.Date(tspan0, origin = "2005-01-01")#"2007-07-01")
    tinObs <- sort(unique(obs$time))
    tObs <- as.numeric(tinObs) - (as.numeric(tinObs[1]) - tspan0[which(realDates == tinObs[1])])
    obs$numTime <- as.numeric(obs$time) - as.numeric(tinObs[1]) +
        tspan0[which(realDates == tinObs[1])]

    tspan <- tspan0[-which(realDates > tail(tinObs,1))]

    tab <- data.frame(matrix(0, ncol = 3, nrow = runs))
    names(tab) <- c("neg", "pos", "prop")

    tab.qtr <- tab

    wm <- numeric(runs)


    ## p <- data.frame(matrix(nrow=length(p1),ncol=runs))
    ## colnames(p) <- seq_len(runs)
    ## do simulation for all samples.
    ##ssResObs <- list()
    pb <- pbapply::timerProgressBar(width = 50)
    for(k in 1:runs) {
        tt <- as.numeric(theta.post[k,])
        names(tt) <- names(theta.post)
        theta <- setTheta(tt)

        extraArgs <- list(tspan = tspan, tObs = tObs, nObs = nObs,
                                   runSeed = NULL, threads = NULL, solver = "ssm",
                                   prevLevel = NULL, prevHerds = NULL, phiLevel = "local",
                                   u0 = NULL, events = NULL, binary = TRUE, model = model,
                                   nSim = 1, obsDates = obs)

        res <- SimInfSimulator_real(theta, extraArgs)[[1]]

        ## ** FOR OUTPUT
        tab[k, ] <- numPos(table(res$sample))

        res.qtr <- qtrStore(res, "time", "sample")
        tab.qtr[k, ] <- numPos(table(res.qtr$sample))

        wm[k] <- meanNodePrev(res.qtr)

        pbapply::setTimerProgressBar(pb, k/runs)
    }
    close(pb)

    ## summarizing the data
    mean.all <- mean(tab$prop)
    sd.all <- sd(tab$prop)

    qt.all <- findQuantile(tab$prop,0.05)

    mean.qtr <- mean(tab.qtr$prop)
    sd.qtr <- sd(tab.qtr$prop)
    qt.qtr <- findQuantile(tab.qtr$prop,0.05)

    df <- data.frame(mean = c(mean.all, mean.qtr),
                     sd = c(sd.all, sd.qtr),
                     qt.low = c(qt.all["low"], qt.qtr["low"]),
                     qt.up = c(qt.all["up"], qt.qtr["up"]))

    rownames(df) <- c("all", "qtr")

    qt.w <- findQuantile(wm)
    df.w <- data.frame(obs = wm.obs, mean = mean(wm), sd = sd(wm), qtLow = qt.w["low"], qtUp = qt.w["up"])

    ##


    return(list(obs = tab.obs,  obs.qtr = tab.obs.qtr, tab = tab, tab.qtr = tab.qtr,
                df = df, df.w = df.w))
}

##' number of positive nodes
numPos <- function(tab) {
    neg <- tab[[1]]
    pos <- sum(tab[-1])

    return(data.frame(neg = neg, pos = pos, prop = pos/(pos + neg)))
}

##' Find the average node prevalence.
meanNodePrev <- function(res) {
    ## transform into a matrix
    m <- matrix(c(res$sample), nrow = length(unique(res$time)))

    ## compute the node prevalence at each quarter
    nprev <- apply(m, 1, function(x){
        tab <- table(x)
        pos <- sum(tab[-1])
        neg <- tab[[1]]
        nprev <- pos / (pos + neg)
    })

    ## compte the weights for each quarter measure
    res.n <- matrix(res$number, nrow = length(unique(res$time)))

    w <- computeWeights(res.n)

    ## compute the (weighted) temporal mean
    wm <- weighted.mean(nprev, w)
    return(wm)
}

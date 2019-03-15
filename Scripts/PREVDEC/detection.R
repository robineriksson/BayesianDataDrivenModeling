library(SimInfInference)
library(SimInf)
library(ggplot2)

## % *** intro here, look back to the 1st paragraph of the abstract
## % (``public security'')

## % *** something about how bad EHEC/VTEC is to the public (see siminf1)

## % sample N = 100 parameters
## % simulate independent simulations for, e.g., 10 years, taking
## % measurements only
## % at the final 5 years or so
## % RECORD the nodes: phi and/or the highest P = I/(S+I)

## % AFTER the simulation, RANK the nodes and find the 10--100 or so most
## % likely to be detected (see below)

## % RE-RUN with this information and add ~2--3 other sets of ``good
## % nodes'' {e.g., in-degree, #individuals, random}

## % discuss: formulate as a hypothesis (e.g., "set 1 is best at
## % detecting"), and a null hypothesis ("there is no difference")

## % *** Fig here: probability of detection vs. time (so no statistical
## % test)

## % Note: probability of detection at time t = P_detect(t) =
## % 1-prod_{node i in set} (1-P_i(t)), P_i(t) = sigmoid(phi) or
## % sigmoid(P), where sigmoid(.) goes from 0 to 1 and is controlled by a
## % cut-off c and a sharpness parameter eps. The change of value goes
## % from 0 to 1 around the value c, and it does so in an interval of
## % length about eps.


##' Run the detection experient found in the paper
##' @param N the number of posterior samples
##' @param M the number of top nodes to pick
##' @param rankP detect on prevalence or phi
##' @param toSave save the output plot
##' @param PosteriorFilePath the directory that holds the posterior
##' @param dataDir path to data directory
main <- function(N = 100, M = 10, rankP = FALSE, PosteriorFilePath, dataDir) {
    set.seed(0)

    ## perform node network check here so that we get the same random
    ## nodes each time.
    networkNodes <- networkEval(M, dataDir)


    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## % sample N = 100 parameters
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print("Draw posterior")
    theta.post <- drawPosterior(N, PosteriorFilePath, dataDir)


    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## simulate independent simulations for, e.g., 10 years, taking
    ## measurements only at the final 5 years or so RECORD the nodes:
    ## phi and/or the highest P = I/(S+I)
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print("Init run")
    run <- runSimulation(init = TRUE, N = N, theta.post = theta.post,
                       nodeSelection = NULL, dataDir = dataDir)



    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## AFTER the simulation, RANK the nodes and find the 10--100 or so
    ## most likely to be detected add ~2--3 other sets of ``good
    ## nodes'' {e.g., in-degree, #individuals, random}
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print("Rank nodes")
    topNodes <- rankNodes(run, M, rankP = rankP)
    nodeSelection <- data.frame("Simulation" = topNodes, networkNodes)



    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## RE-RUN with this information. (Now we want the mean + std
    ## result)
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print("Re run")
    reRun <- runSimulation(init = FALSE, N = N, theta.post = theta.post,
                           nodeSelection = nodeSelection, dataDir = dataDir)


    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Fig idea: probability of detection vs. time
    ##
    ## where probability of detection at time t:
    ## P_detect(t) = 1-prod_{node i in set} (1-P_i(t)),
    ## P_i(t) = sigmoid(phi) or sigmoid(P)
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p <- detectionProb(reRun, nodeSelection, rankP = FALSE, na.rm = FALSE)

    return(list(run = run, nodeSelection = nodeSelection, reRun = reRun, theta = theta.post, plot = p))
}

##' Draw parameters from posterior
##' @param N number of draws
drawPosterior <- function(N, filename) {
    ## load posterior
    load(filename)
    posterior.All <- sobs$getPosterior()
    dims <- dim(posterior.All)
    rows <- seq(1,dims[1],100)
    cols <- seq_len(dims[2]-1)

    posterior.thinned <- posterior.All[rows,cols]

    ## if we do 100 draws and do not use replace, then we're using all the accepted parameter pairs
    theta.post <- posterior.thinned[sample(seq_len(length(rows)), N, replace = T),]

    return(theta.post)
}

##' Run initial simulation for ranking of nodes
##' @param init is it the initial run?
##' @param N number of simulation/post draws
##' @param theta.post sampled posterior parameters
##' @param dataDir the path to the directory that holds the data
runSimulation <- function(init = TRUE, N, theta.post, nodeSelection = NULL, dataDir) {
    ## define tspan vector
    sise_smhi_path <- paste(dataDir, "SISe_smhi.rda", sep = "")
    load(sise_smhi_path) ## load model
    tspan <- seq(head(model@events@time,1), tail(model@events@time,1), 1)

    tmeasurment <- seq(round(median(tspan)),tail(tspan,1),15)



    cols <- tmeasurment
    ncol <- length(cols)
    if(init)
        rows <- 1:dim(model@ldata)[2]
    else
        rows <- as.numeric(unlist(nodeSelection))

    ## SimInf::trajectory sorts (increasing) the nodes, so we need to
    ## keep track on which order we want to keep.
    rows.order <- order(rows)

    nrow <- length(rows)

    v0 <- array(numeric(N*nrow*ncol), c(nrow, ncol, N))
    prev0 <- array(numeric(N*nrow*ncol), c(nrow, ncol, N))

    for(n in seq_len(N)) {
        ## extract theta(n)
        ## + class massaging
        t <- as.numeric(theta.post[n,])
        names(t) <- names(theta.post)
        theta <- setTheta(t)

        ## reload model
        load(sise_smhi_path)

        ## initialize model
        model <- init_model_widgren(model = model, theta = theta,
                                    tspan = tspan,
                                    phi = "local")

        res <- SimInf::run(model, threads = NULL, solver = "ssm")

        ## ## extract final state, to be used as
        ## u <- SimInf::trajectory(res, c("S","I"), as.is = TRUE)
        ## size <- dim(u)
        ## u0.m <- matrix(u[,size[2]], nrow = 2)
        ## rownames(u0.m) <- c("S","I")



        v <- SimInf::trajectory(res, c("phi"), node = rows, as.is = TRUE)
        v.m <- matrix(v[,which(as.numeric(colnames(v)) %in% cols)], ncol = ncol)


        v0[rows.order,,n] <- v.m


        p <- SimInf::trajectory(res, c("S","I"), node = rows, as.is = TRUE)
        p.m <- p[,which(as.numeric(colnames(p)) %in% cols)]

        S <- p.m[which(rownames(p.m) == "S"),]
        I <- p.m[which(rownames(p.m) == "I"),]
        P <- I / (S+I)

        prev0[rows.order,,n] <- as.matrix(P)

    }

    return(list(prev0 = prev0, v0 = v0, time = tmeasurment))

}

##' Rank the nodes from the simulation
##' @param run output from main function
##' @param M number of selected nodes
##' @param rankP do ranking on P or phi
rankNodes <- function(run, M = 10, rankP = FALSE, na.rm = FALSE) {
    if(rankP)
        value <- run$prev0
    else {
        value <- run$v0
        value <- (value - min(value)) / (max(value) - min(value))
    }

    ## mean and sd over repeats
    meanVals <- apply(X = value, c(1,2), function(x){mean(x, na.rm = na.rm)})
    sdVals <- apply(value, c(1,2), function(x){sd(x, na.rm = na.rm)})

    ## 99.6%(?) of all values included
    maxVals <- meanVals + 3*sdVals

    ## temporal mean over each node
    meanInRows <- apply(maxVals, 1, function(x){mean(x, na.rm = na.rm)})

    ## find which node that has the largest temporal mean
    maxOrder <- order(meanInRows, decreasing = TRUE)

    rowsMax <- maxOrder[1:M]
    ## rows100 <- maxOrder[101:(100+M)]
    ## rows1000 <- maxOrder[1001:(1000+M)]
    return("Simulation" = rowsMax)
    ##return(data.frame(rowsMax, rows100 = rows100, rows1000 = rows1000))
}

##' Evaluate the network trafic for possible candidates
##' @param M number of selected nodes
##' @param dataDir the path to the drectory that holds the data
networkEval <- function(M = 10, dataDir) {
    sise_smhi_path <- paste(dataDir, "SISe_smhi.rda", sep = "")
    load(sise_smhi_path) ## loads: model
    events <- model@events



    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## indegree
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## here, only count external transfers in indegree
    ## which events are out?
    send <- which(events@event == 3)
    #3 sent to node dest[send]
    send.node = events@dest[send]
    send.n = events@n[send]
    send.num <- numeric(sum(send.n)) ## put all in row for easy table!
    last <- 1
    for(i in 1:length(send.n)) {
        num <- send.n[i]
        send.num[last:(last+num)] <- send.node[i]
        last <- last+num
    }
    send.table <- table(send.num)
    indegree <- head(order(send.table, decreasing = TRUE), M)

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## outdegree
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## here, only count external transfers in outdegree
    ## sent from node[send]
    send.node = events@node[send]
    send.n = events@n[send]
    send.num <- numeric(sum(send.n)) ## put all in row for easy table!
    last <- 1
    for(i in 1:length(send.n)) {
        num <- send.n[i]
        send.num[last:(last+num)] <- send.node[i]
        last <- last+num
    }
    send.table <- table(send.num)
    outdegree <- head(order(send.table, decreasing = TRUE), M)

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## largest number of animals (u0)
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u0 <- model@u0["I",] ## the animals are all in I compartment
    maxAnim <- head(order(u0, decreasing = TRUE), M)


    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## random selection
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## build vector of all unique nodes
    nodes <- seq_len(dim(model@ldata)[2])
    randomNodes <- sample(nodes, M)


    ## collect in dataframe
    output <- data.frame("Indegree" = indegree,
                         "Outdegree" = outdegree,
                         "Largest" = maxAnim,
                         "Random" = randomNodes)
    return(output)
}


plotReRun <- function(run, nodeSelection, na.rm = FALSE, rankP = FALSE) {
    if(rankP)
        value <- run$prev0
    else {
        value <- run$v0
        value <- (value - min(value)) / (max(value) - min(value))
    }
    M1 <- dim(nodeSelection)[1]
    M2 <- dim(nodeSelection)[2]


    meanVec <- c()
    sdVec <- c()
    nameVec <- c()

    for(i in 0:(M2-1)) {
        m <- i*10 + (1:10)

        ## mean over nodes
        value.mean <- t(apply(value[m,,], c(2,3), function(x){mean(x, na.rm = na.rm)}))

        ## mean and sd over replicates
        pm <- apply(value.mean, 2, function(x){mean(x, na.rm = na.rm)})
        psd <- apply(value.mean, 2, function(x){sd(x, na.rm = na.rm)})

        meanVec <- c(meanVec,pm)
        sdVec <- c(sdVec,psd)
        nameVec <- c(nameVec, rep(names(nodeSelection)[i+1], length(pm)))

    }

    df <- data.frame(time = rep(1:length(pm), M2),
                     ymean = meanVec,
                     ymin = meanVec - sdVec,
                     ymax = meanVec + sdVec,
                     name = nameVec)



    yminPos <- numeric(length(df$ymin))
    for(i in 1:length(yminPos))
        yminPos[i] <- max(0, df$ymin[i])

    df$yminPos <- yminPos




    gp <- ggplot2::ggplot(data = df)
    ## mean lines
    gp <- gp + ggplot2::geom_line(ggplot2::aes(x = time, y = ymean, color = name), lwd = 0.5)

    ## ribbon lines
    gp <- gp + ggplot2::geom_ribbon(ggplot2::aes(x = time, ymin = yminPos,
                                                 ymax = ymax, fill = name),
                                    alpha = 0.1)#, fill = colors[n])

    gp <- gp + theme_Publication()

    gp

    return(list(df = df, gp = gp))


}

##' Sigmoid(.) goes from 0 to 1 and is controlled by a
##' cut-off c and a sharpness parameter k. The change of value goes
##' from 0 to 1 around the value p0, and it does so in an interval of
##' length about eps.
##' @param p value to be evaluated by the function
##' @param p0 the x-value of the sigmoid's midpoint
##' @param k the logistic growth rate
sigmoid <- function(p, p0 = 0.3, k = 20) {
    F = 1 / ( 1 + exp(-k*(p-p0)))
    return(F)
}


plotSig <- function(p0 = 0.3, k = 20) {
    x <- seq(0,1,length.out = 1e6)
    s <- sigmoid(x, p0 = p0, k = k)

    ## 95% interval
    prange <- c(0.025, 0.975)
    inside <- which(s > prange[1] & s < prange[2])
    edge <- c(inside[1], tail(inside,1))

    print(sprintf("95%% conf. int.: [%f, %f]", x[edge[1]], x[edge[2]]))

    ## plot line
    plot(x,s, type = "l")

    ## plot area
    lines(x[inside],s[inside], col = "red")
    abline(h = prange, v = x[edge])

    ## plot s = 0.5
    abline(h = 0.5, v = p0, col = "green")

    ## plot derivative (CDF -> PDF)
    l <- 2:length(x)
    d <- (s[l]-s[l-1]) / (x[l] - x[l-1])
    d <- d/max(d)
    lines(x[l],d, col = "blue")

}

##' probability of detection at time t = P_detect(t) =
##' 1-prod_{node i in set} (1-P_i(t)), P_i(t) = sigmoid(phi) or
##' sigmoid(P),
##' @param p matrix for node set
detectionProb <- function(run, nodeSelection, rankP = FALSE, na.rm = FALSE, cutLevel = 1, numSD = 1, alt1 = TRUE) {
    if(rankP)
        value <- run$prev0
    else {
        value <- run$v0
        ## normalize value
        value <- value/max(value)
    }

    dims <- dim(value)

    sets <- dim(nodeSelection)[2]
    setSize <- dim(nodeSelection)[1]


    if(alt1) {
        prob <- array(0, c(sets,dims[2],dims[3]))

        ## for each nodeset we run apply the sigmoid
        for(i in 0:(sets-1)) {
            ii <- 1:10 + i*10
            sig <- sigmoid(value[ii,,])
            prob[i+1,,] <- 1 - apply(sig, c(2,3), function(x){prod(1 - x, na.rm = na.rm)})

        }

        ## Compute the probability of detecting each nodeset
        allProbs <- as.data.frame(matrix(NA, nrow = sets*dims[2], ncol = 4))
        names(allProbs) <- c("mean", "max", "min", "name")

        allProbs$mean <- c(t(apply(prob, c(1,2), mean)))
        probSD <- c(t(apply(prob, c(1,2), sd)))
        allProbs$max <- allProbs$mean + numSD*probSD
        allProbs$min <- allProbs$mean - numSD*probSD

        allProbs$name <- rep(names(nodeSelection), each = dims[2])

    } else {
        ## First we take the mean over all posterior samples
        postMean <- apply(value, c(1,2), mean)

        ## then add the sd for X % CI
        postSD <- apply(value, c(1,2), sd)
        postMean.up <- postMean + numSD*postSD
        postMean.down <- postMean - numSD*postSD

        ## then apply the detection filter
        postMean.sig <- sigmoid(postMean)
        postMean.up.sig <- sigmoid(postMean.up)
        postMean.down.sig <- sigmoid(postMean.down)

        ## Compute the probability of detecting each nodeset
        allProbs <- as.data.frame(matrix(NA, nrow = sets*dims[2], ncol = 4))
        names(allProbs) <- c("mean", "max", "min", "name")

        for(i in 0:(sets-1)) {
            ii <- 1:setSize + i*setSize
            sig.mean <- postMean.sig[ii,]
            sig.up = postMean.up.sig[ii,]
            sig.down = postMean.down.sig[ii,]


            allProbs$mean[1:dims[2] + i*dims[2]] <- 1 - apply(sig.mean, 2, function(x){prod(1 - x)})
            allProbs$max[1:dims[2] + i*dims[2]] <- 1 - apply(sig.up, 2, function(x){prod(1 - x)})
            allProbs$min[1:dims[2] + i*dims[2]] <- 1 - apply(sig.down, 2, function(x){prod(1 - x)})
            allProbs$name[1:dims[2] + i*dims[2]] <- rep(names(nodeSelection)[i+1], dims[2])
        }



    }


    ymin <- allProbs$min
    ymax <- allProbs$max
    ymean <- allProbs$mean


    yminPos <- numeric(length(ymin))
    ymaxPos <- numeric(length(ymax))
    for(i in 1:length(yminPos)) {
        yminPos[i] <- max(0, allProbs$min[i])
        ymaxPos[i] <- min(cutLevel, allProbs$max[i])
    }


    time = as.Date(run$time, origin = "2004-12-31")
    allProbs$time <- rep(time, sets)

    allProbs$yminPos <- yminPos
    allProbs$ymaxPos <- ymaxPos


    gp <- ggplot2::ggplot(data = allProbs)
    ## mean lines
    gp <- gp + ggplot2::geom_line(ggplot2::aes(x = time, y = ymean, color = name), lwd = 0.5)

    ## ribbon lines
    gp <- gp + ggplot2::geom_ribbon(ggplot2::aes(x = time, ymin = yminPos,
                                                 ymax = ymaxPos, fill = name),
                                    alpha = 0.1)#, fill = colors[n])
    ## percentages on y axis
    gp <- gp + ggplot2::scale_y_continuous(labels = scales::percent)

    ## theme-ing
    gp <- gp + theme_Publication(6) +
        labs(x = "Year", y = "Probability of detection") +
        labs(color="Node selection", fill = "Node selection") +
        theme(plot.margin=grid::unit(c(1,1,1,1),"mm"),
              legend.margin=margin(t = -0.1, unit='cm'),
              legend.position="top",
              #axis.title.y = element_text(margin = margin(t = 0, r = -5, b = 0, l = 0)),
              axis.text.x = element_text(angle = 0, hjust = 0.5),
              axis.title = element_text(face = "plain",size = rel(1))) +
        guides(color=guide_legend(nrow=1,byrow=TRUE), fill=guide_legend(nrow=1,byrow=TRUE))

    return(gp)

}

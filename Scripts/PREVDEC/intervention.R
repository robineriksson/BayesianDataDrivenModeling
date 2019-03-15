## Robin Eriksson
library(SimInf)
library(SimInfInference)

library(ggplot2)



##' main function, that employs all subfunctions
##'
##' @param N the number of posterior samples
##' @param savePlot save the output plot
##' @param PosteriorFilePath filepath to the posterior.
##' @param dataDir path to data directory
##' @return list with dataframe and plot handle
main <- function(N, PosteriorFilePath, dataDir) {
    res <- runSimulator(N, PosteriorFilePath, dataDir)
    ## create more informative names for each intervention
    oldnames <- names(res)
    newnames <- numeric(length(oldnames))
    for(n in seq_len(length(oldnames))) {
        if(oldnames[n] == "baseline")
            newnames[n] <- "No intervention"
        else if(oldnames[n] == "transrate")
            newnames[n] <- "10% reduced indirect transmission rate"
        else if(oldnames[n] == "decayrate")
            newnames[n] <- "10% increased bacterial decay rate"
        else if(oldnames[n] == "noInfectTransp")
            newnames[n] <- "No transmission by transport"
        else
            newnames[n] <- oldnames[n]
    }
    names(res) <- newnames

    ## create data frame of each N simulations.
    prev <- lapply(res,prevSimulator)

    ## extract values post 2008
    prev.2008 <- lapply(prev, post2008)

    ## thin out the data
    prev.2008t <- lapply(prev.2008, postThin)

    ## create a ribbon, mean & sd
    ribbon <- do.call("rbind", lapply(prev.2008t, postRibbon))
    ## drop rownames
    rownames(ribbon) <- c()

    ## add "factor name"
    resnames <- names(res)
    internames <- numeric(length(resnames))
    for(i in seq_len(length(resnames))) {
        internames[i] <- paste(i,") ", resnames[i], sep = "")
    }
    ribbon$name <- rep(internames, each = dim(ribbon)[1]/length(resnames))
    p <- plotSimulator(ribbon)

    return(list(fulldata = prev, ribbon = ribbon, p = p))
}

##' Perform the statistical test.
##' Are there any difference in the prevention methods?
##'
##' @param prev the full data produced by the main function
##' @param ribbon the ribbon produced by the main function.
statTest <- function(prev, logit = TRUE, B = 1000) {
    ## in order not to recompute every time
    ncols <- length(prev)

    ## prepare the data
    prev.2008 <- lapply(prev, post2008)
    prev.2008t <- lapply(prev.2008, postThin)
    data <- do.call("rbind",prev.2008t)

    ## get max time
    mtime <- max(data$time)

    ## extract the max time data
    maxdata <- t(data[which(data$time == mtime),])
    maxdata <- matrix(as.numeric(maxdata[-1,]), ncol = ncols)

    ## logit?
    if(logit) {
        maxdata <- log(maxdata)
    }

    ## Do we have normality?
    swVal <- c()

    ## plot histograms with fitted normal curve
    par(mfrow = c(ceiling(ncols/2),2))
    for(i in seq_len(ncols)) {
        if(B > 0) {
            ## non parametric bootstrap
            xx <- sample(maxdata[,i], size = B, replace = T)
        } else
            xx <-  maxdata[,i]

        swVal <- c(swVal,shapiro.test(xx)$p.value)

        hist(xx, freq = F)
        ## think they are in the same order.
        curve(dnorm(x, mean = mean(xx), sd = sd(xx)), col = 1+i, add =T)
    }

    ## perform normality test
    #swVal <- apply(maxdata, 2, function(x){shapiro.test(x)$p.value})

    alpha <- 0.05
    swTest <- swVal > alpha
    return(swTest)
}


##' run N simulations with N drawn posterior parameter samples.
##'
##' @param N number of simulations
##' @param filename the directory that holds the posterior
##' @param dataDir path to data directory
##' @return result
runSimulator <- function(N, filename, dataDir) {
    set.seed(0)

    ## load posterior
    load(filename)
    posterior.All <- sobs$getPosterior()
    dims <- dim(posterior.All)
    rows <- seq(1,dims[1],100)
    cols <- seq_len(dims[2]-1)

    posterior.thinned <- posterior.All[rows,cols]

    ## if we do 100 draws and do not use replace, then we're using all the accepted parameter pairs
    theta.post <- posterior.thinned[sample(seq_len(length(rows)), N, replace = T),]


    ## define tspan vector
    sise_smhi_path <- paste(dataDir, "SISe_smhi.rda", sep = "")
    load(sise_smhi_path) ## load model
    tspan0 <- seq(head(model@events@time,1), tail(model@events@time,1), 1)

    ## 0. Pre run until intervetion day
    ## save these models
    ## then continue them with the desired prevention method.

    datetime <- as.Date(tspan0, origin = "2004-12-31")

    ## u0 is based on the day before we apply
    u0Date <- "2008-12-31"
    u0Date.elem <- which(datetime == u0Date)
    u0Date.num <- tspan0[u0Date.elem]


    tspan1 <- tspan0[1:u0Date.elem]
    tspan2 <- tspan0[(u0Date.elem+1):length(tspan0)]


    print("pre intervention")
    u0 <- list()
    v0 <- list()
    prev0 <- list()
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
                                    tspan = tspan1,
                                    phi = "local")

        res <- SimInf::run(model, threads = NULL, solver = "ssm")

        ## extract final state, to be used as
        u <- SimInf::trajectory(res, c("S","I"), as.is = TRUE)
        size <- dim(u)
        u0.m <- matrix(u[,size[2]], nrow = 2)
        rownames(u0.m) <- c("S","I")

        u0[[n]] <- u0.m

        v <- SimInf::trajectory(res, c("phi"), as.is = TRUE)
        v0.m <- matrix(v[,size[2]], nrow = 1)
        rownames(v0.m) <- c("phi")

        v0[[n]] <- v0.m

        prev0[[n]] <- SimInf::prevalence(res, I ~., "pop")
    }


    ## the events need also to be transfered to the fast-forwarded models.
    events <- model@events
    e0 <- which(events@time == (u0Date.num + 1))[1]
    elen <- length(events@event)


    events@event <- events@event[e0:elen]
    events@time <- events@time[e0:elen]
    events@node <- events@node[e0:elen]
    events@dest <- events@dest[e0:elen]
    events@n <- events@n[e0:elen]
    events@proportion <- events@proportion[e0:elen]
    events@select <- events@select[e0:elen]
    events@shift <- events@shift[e0:elen]



    print("Baseline")
    ## 1.
    ## Run baseline case
    baseline <- list()
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
                                    tspan = tspan2,
                                    phi = "local")


        model@u0 <- u0[[n]]
        model@v0 <- v0[[n]]
        model@events <- events
        res <- try(SimInf::run(model, threads = NULL, solver = "ssm"))


        baseline[[n]] <- rbind(prev0[[n]], SimInf::prevalence(res, I ~., "pop"))
    }






    ## ***************************************************
    ## Evaluating prevention
    ## ***************************************************

    ## 2. Reduce upsilon by 10%, (indirect trans. rate)
    print("Reduce upsilon")
    transrate <- list()
    for(n in seq_len(N)) {
        ## extract theta(n)
        ## + class massaging
        t <- as.numeric(theta.post[n,])
        names(t) <- names(theta.post)

        t["upsilon"] <- t["upsilon"]*0.9
        theta <- setTheta(t)

         ## reload model
        load(sise_smhi_path)

        ## initialize model
        model <- init_model_widgren(model = model, theta = theta,
                                    tspan = tspan2,
                                    phi = "local")


        model@u0 <- u0[[n]]
        model@v0 <- v0[[n]]
        model@events <- events
        res <- try(SimInf::run(model, threads = NULL, solver = "ssm"))

        transrate[[n]] <- rbind(prev0[[n]], SimInf::prevalence(res, I ~., "pop"))
    }

    ## 3. Increase beta by 10%, (decay of env. inf. pressure per day)
    print("Increase beta")
    decayrate <- list()
    for(n in seq_len(N)) {
        ## extract theta(n)
        ## + class massaging
        t <- as.numeric(theta.post[n,])
        names(t) <- names(theta.post)

        t["beta_t1"] <- t["beta_t1"]*1.1
        t["beta_t2"] <- t["beta_t2"]*1.1
        t["beta_t3"] <- t["beta_t3"]*1.1
        theta <- setTheta(t)

         ## reload model
        load(sise_smhi_path)

        ## initialize model
        model <- init_model_widgren(model = model, theta = theta,
                                    tspan = tspan2,
                                    phi = "local")


        model@u0 <- u0[[n]]
        model@v0 <- v0[[n]]
        model@events <- events
        res <- try(SimInf::run(model, threads = NULL, solver = "ssm"))


        decayrate[[n]] <- rbind(prev0[[n]], SimInf::prevalence(res, I ~., "pop"))
    }

    ## 4. Infected induviduals are "cured" on transport
    events <- model@events
    nmat <- matrix(c(0,-1),nrow=2)
    mode(nmat) <- "integer"
    rownames(nmat) <- c("S","I")
    colnames(nmat) <- c("1")
    events@N <- nmat

    events@shift <- as.integer(events@shift + 1)



    print("No infection transports")
    noInfectTransp <- list()
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
                                    tspan = tspan2,
                                    phi = "local")


        model@u0 <- u0[[n]]
        model@v0 <- v0[[n]]
        model@events <- events
        res <- try(SimInf::run(model, threads = NULL, solver = "ssm"))

        noInfectTransp[[n]] <- rbind(prev0[[n]], SimInf::prevalence(res, I ~., "pop"))
    }



    print("done")
    return(list(baseline = baseline, noInfectTransp = noInfectTransp,
                decayrate = decayrate, transrate = transrate))
}


##' Extract the prevalence from the simulator result
##'
##' @param results a list of model results
##' @return a dataframe with the prevalence and time
prevSimulator <- function(results) {
    ## extract the aggregated prevalence vector from each result
    ## pmat <- lapply(results, function(X) {
    ##     SimInf::prevalence(model = X, formula = I ~., type = type)
    ## })
    df <- results[[1]]
    oldNames <- c("time", "p1")
    for(i in 2:length(results)) {
        df[,i+1] <- results[[i]]$prevalence
        names(df) <- c(oldNames, paste("p",i,sep=""))
        oldNames <- names(df)
    }
    df$time <- as.Date(df$time, origin = "2004-12-31")

    return(df)
}

##' Cut the timeseries to only show post 2008.
##'
##' @param df dataframe with prevalence from the result
post2008 <- function(df) {
    post2008 <- df[which(df$time > "2007-12-31"), ]
    return(post2008)
}

##' Thin out the dataseries
##'
##' @param df dataframe with prevalence from th result
postThin <- function(df, n = 25) {
    postThin <- df[seq(1,dim(df)[1], n),]
    return(postThin)
}


postRibbon <- function(df) {
    dims <- dim(df)
    dribbon = data.frame(time = df$time,
                         ymean = apply(df[,2:dims[2]],1,mean),
                         ysd = apply(df[,2:dims[2]],1,sd))

    return(dribbon)

}
##' Plot the timeseries from the dataframe
##'
##' @param df lift the dataframes with prevalence from the result
##' @return a plot handle
plotSimulator <- function(df, publish = FALSE, numSD = 2) {

    newnames <- character(4)

    newnames[1] <- "No intervention"
    newnames[2] <- "10% reduced indirect transmission rate"
    newnames[3] <- "10% incr bacterial decay rate"
    newnames[4] <- "No transmission by transport"
    ## browser()
    ## names(df) <- newnames

    gp <- ggplot2::ggplot(data = df)
    ## mean lines

    ## ribbon lines
    gp <- gp + ggplot2::geom_ribbon(ggplot2::aes(x = time, ymin = ymean - numSD*ysd,
                                                 ymax = ymean + numSD*ysd, fill = name),
                                    alpha = 0.1)#, fill = colors[n])

    ## get end for legend
    d_end <- df[which(df$t == tail(sort(df$t),1)),]
    sidename <- as.list(1:dim(d_end)[1])
    d_end$name <- sapply(sidename, function(x){paste(x,")",sep="")})

    ##d_end[d_end$name == "2)", "ymean"] = d_end[d_end$name == "2)","ymean"]*0.1
    gp <- gp + ggplot2::geom_text(data = d_end,
                                  ggplot2::aes(x = time,
                                               y = ymean,
                                               label=name),
                                  hjust = -1, vjust = 0,
                                  fontface = "bold")

    gp <- gp + ggplot2::geom_vline(xintercept = 14230) ## 2009-01-01 in numeric 14245
    gp <- gp + ggplot2::geom_vline(xintercept = 14610, lty = 2, color = "#00B283") ## 2009-01-01 in numeric 14245
    gp <- gp + ggplot2::geom_vline(xintercept = 15340, lty = 3, color = "#FF3B00") ## 2009-01-01 in numeric 14245

    d_times <- data.frame(time = as.Date(c(1445,1825,2555), origin = "2004-12-31"),
                          level = max(df$ymean+df$ysd),
                          name = c("Year 0", "Year 1", "Year 3"))




    ##
    ##gp <- gp + ggthemes::theme_economist_white()
    if(publish) {
        gp <- gp + theme_Publication(6)
        gp <- gp + ggplot2::geom_line(ggplot2::aes(x = time, y = ymean, color = name), lwd = 0.5)
        gp <- gp + ggplot2::geom_text(data = d_end,
                                      ggplot2::aes(x = time,
                                                   y = ymean,
                                                   label=name),
                                      hjust = -1, vjust = 0,
                                      fontface = "bold", size = 2)
        gp <- gp + ggplot2::geom_text(data = d_times,
                                      ggplot2::aes(x = time,
                                                   y = level,
                                                   label=name),
                                      hjust = -0.1,
                                      vjust = 1,
                                      fontface = "bold", size = 2
                                      )
        gp <- gp + theme(legend.position="top",#c(1,0.5),
                         legend.direction = "horizontal",
                         #legend.background = element_rect(fill="lightblue",
                         #                                 size=1, linetype="solid"),
                         #legend.key.width = unit(0.5,"cm"),
                         axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                         axis.title = ggplot2::element_text(face = "plain",size = ggplot2::rel(1)),
                         plot.margin=grid::unit(c(1,1,1,1),"mm"),
                         legend.margin=margin(l = -0.5, unit='cm'))
        #gp <- gp + legend.key(label=abbreviate)

    } else {
        gp <- gp + theme_Publication()
        gp <- gp + ggplot2::geom_line(ggplot2::aes(x = time, y = ymean, color = name), size = 1)
        gp <- gp + ggplot2::geom_text(data = d_end,
                                      ggplot2::aes(x = time,
                                                   y = ymean,
                                                   label=name),
                                      hjust = -1, vjust = 0,
                                      fontface = "bold")
        gp <- gp + ggplot2::geom_text(data = d_times,
                                      ggplot2::aes(x = time,
                                                   y = level,
                                                   label=name),
                                      hjust = -0.1,
                                      vjust = 1,
                                      fontface = "bold"

                                      )

    }

    gp <- gp + ggplot2::scale_y_continuous(labels = scales::percent)
    gp <- gp + ggplot2::labs(x = "Year", y = "Population prevalence",
                             colour = "Intervention:", fill ="Intervention:")
    gp <- gp + ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE),
                               fill = ggplot2::guide_legend(nrow=2,byrow=TRUE))

    #gp <- gp + ggplot2::theme(legend.position = "top")


    return(gp)
}




##' Perform the statistical test.
##' Are there any difference in the prevention methods?
##'
##' @param prev the full data produced by the main function
intCompare <- function(mainoutput) {
    prev <- mainoutput$fulldata
    ## in order not to recompute every time
    ncols <- length(prev)

    ## prepare the data
    prev.2008 <- lapply(prev, post2008)
    ##prev.2008t <- lapply(prev.2008, postThin)
    ##data <- do.call("rbind",prev.2008t)
    data <- do.call("rbind",prev.2008)

    ## get two comparative year slots
    t1 <- "2010-01-01"
    t2 <- "2012-01-01"


    t1data <- t(data[which(data$time == t1),])
    t1matrix <- matrix(as.numeric(t1data[-1,]), ncol = ncols)
    t1matrix <- sweep(t1matrix,1,t1matrix[,1],"/")

    mean1 <- apply(t1matrix,2,mean)
    sd1 <- apply(t1matrix,2,mean)
    tf1data <- data.frame(t1matrix, as.character(rep("Year 1", dim(t1matrix)[1])))
    names(tf1data) <- c(names(prev),"time")


    t2data <- t(data[which(data$time == t2),])
    t2matrix <- matrix(as.numeric(t2data[-1,]), ncol = ncols)
    t2matrix <- sweep(t2matrix,1,t2matrix[,1],"/")
    tf2data <- data.frame(t2matrix, as.character(rep("Year 3", dim(t2matrix)[1])))

    names(tf2data) <- c(names(prev),"time")




    tfdata <- data.frame(rbind(tf1data,tf2data))
    names(tfdata) <- c("(1)","(2)","(3)","(4)","time")
    ##names(tfdata) <- c(unique(mainoutput$ribbon$name),"time")

    means1 <- apply(tfdata[tfdata$time == "Year 1", 2:4], 2, mean)
    sd1 <- apply(tfdata[tfdata$time == "Year 1", 2:4], 2, sd)

    means2 <- apply(tfdata[tfdata$time == "Year 3", 2:4], 2, mean)
    sd2 <- apply(tfdata[tfdata$time == "Year 3", 2:4], 2, sd)

    sumData <- data.frame(means1, sd1, means2, sd2)
    rownames(sumData) <- names(tfdata)[2:4]

    tfdata <- tfdata[,-1]
    tfdata.melt <- reshape::melt(tfdata, id.var = c("time"))


    ## boxplot(tfdata[tfdata$time == 2, 1:ncols],col= scales::alpha("red", 0.5), las = 2)
    ## boxplot(tfdata[tfdata$time == 1, 1:ncols],col= scales::alpha("blue", 0.5), add = T, las = 2)

    p <- ggplot2::ggplot(tfdata.melt, ggplot2::aes(x = variable, y = value, color = time)) +
        ggplot2::geom_boxplot(lwd = 0.1, outlier.size=NA, fatten = 0.5) + #coord_flip() +
        theme_Publication(6) +
        ggplot2::labs(x = "Intervention method", y = "Reduction factor", color = "time") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                       legend.position = c(0.85,0.8),
                       legend.direction = "vertical",
                       ##legend.key.size= grid::unit(0.1, "cm"),
                       legend.background = element_rect(colour='grey'),
                       panel.grid.major.x = element_blank()) +
        ##scale_y_continuous(position = "right") +
        scale_color_manual(values = c("#00B283","#FF3B00")) +
    guides(color=ggplot2::guide_legend(title.vjust = -1, title = NULL, label.vjust = 2, key.vjust = 2))

    return(list(p = p, sumData = sumData))

}

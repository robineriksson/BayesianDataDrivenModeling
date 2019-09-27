##' Siminf summary statistics
##' input should be
##' ## simulation data
##' ## nodes to sample
##' ## time stamps to sample at (measurements)
##' output should be
##' ### vector where each element is a summary statistic
##'
##' Different "main" functions here should act at different
##' vectors returned, i.e, different summary statistic combinations.
##' Could be usefull if we use different models, as SS choice is
##' often heuristic.
##'
##' Robin Eriksson @ 2017

##' aggregate the time series from the simulated data.
##' and extract the Wilkinson statistics
##' @param data simulated data frame
##' @param extraArgs column in which data i found
##' @export
aggWilkinson <- function(data, extraArgs){
    stopifnot(c("bs", "useW") %in% names(extraArgs))
    bs <- extraArgs$bs
    useW <- extraArgs$useW

    aggSS.all <- agg(data, extraArgs)
    aggSS <- aggSS.all$data
    name <- aggSS.all$name
    w <- aggSS.all$w

    sumSS.mat <- matrix(aggSS, ncol = length(data))
    if(bs && length(data) != 1)
        sumSS.mat <- bootstrap(sumSS.mat, extraArgs)
    if(!useW)
        w <- NULL

    rownames(sumSS.mat) <- name

    wilkSS <- apply(sumSS.mat, 2, wilkinsonStatistics, w = w)
    return(wilkSS)
}

##' Naive bootstrap draw new samples!
##' @param mat the data matrix
##' @param B the number of bootstrap samples
##' @export
bootstrap <- function(mat,extraArgs){
    stopifnot(c("B") %in% names(extraArgs))
    B <- extraArgs$B

    dims <- dim(mat)
    rows <- dims[1]
    cols <- dims[2]

    boot <- matrix(0, nrow = rows, ncol = B)

    for (row in seq_len(rows)) {
        ## draw samples from the row: row
        draw <- sample(x = mat[row,], size = B, replace = TRUE)
        ## these are the bootstrap draws for the element.
        boot[row,] <- draw
    }


    return(boot)

}

##' basic plotting for the bootstrap replications.
aplot <- function(mat1, mat2){
    t <- 1:dim(mat1)[1]
    maxmax <- max(max(apply(mat1,2,max)),max(apply(mat2,2,max)))

    y0 <- seq(0,maxmax,length.out=length(t))

    plot(t,y0, type = "n")

    ## the bootstrapped
    for(l in 1:dim(mat2)[2]) {
        lines(t, mat2[,l], col = scales::alpha(rgb(0,0,1), 0.005))
    }
    ## mean of the bootstrap
    lines(t, apply(mat2,1,mean), col = scales::alpha(rgb(0,0,1), 0.9))


    ## the original
    for(l in 1:dim(mat1)[2]) {
        points(t, mat1[,l], pch = 16, col = scales::alpha(rgb(1,0,0), 0.1))
    }
    ## mean of the original
    lines(t, apply(mat1,1,mean), col = scales::alpha(rgb(1,0,0), 0.9))



}


##' aggregate the time series from the simulated data, also over quarters!
##' @param data simulated data frame
##' @param extraArgs column in which data i found
##' @export
agg <- function(data, extraArgs){
    stopifnot(c("column", "fun", "qtr", "logical") %in% names(extraArgs))
    column <- extraArgs$column
    fun <- extraArgs$fun
    qtr <- extraArgs$qtr
    logical <- extraArgs$logical

    agged <- lapply(data, function(data){
        ## prepare timeseries vector
        agged <- aggregateHelper(column = column, data = data, fun = fun, qtr = qtr, logical = logical)

        agged.unlist <- unlist(agged)
        data.frame(data = agged.unlist[grep("data", names(agged.unlist))],
                   name = agged.unlist[grep("name", names(agged.unlist))],
                   w = agged.unlist[grep("w", names(agged.unlist))])
    })

    ## same w for all
    w <- agged[[1]]$w
    name <- agged[[1]]$name
    data <- c(sapply(agged, function(x){x$data}))
    return(list(data = data, name = name, w = w))
}


##' A helper function
##' @export
aggregateHelper <- function(column, data, fun, qtr, logical){
    nodes <- unique(data$node)
    if(qtr) {
        ## aggregate into qtrs first
        data.qtrs <- qtrStore(data, "time", column, logical)
        data.mat <- matrix(data.qtrs[,column], ncol = length(nodes))
        numb.mat <- matrix(data.qtrs[,"number"], ncol = length(nodes))

        w <- computeWeights(numb.mat)

    } else {
        data.mat <- matrix(data[order(data$node),column], ncol = length(nodes))
        w <- rep(NA, dim(data.mat)[1])
    }

    data.fun <- apply(data.mat, 1, function(x){fun(x, na.rm = TRUE)})

    if(qtr)
        name <- unique(data.qtrs$time)
    else
        name <- unique(data$time)
    return(list(data = data.fun, name = name, w = w))
}


##' Store data in quarters.
##'
##' @param data as list
##' @return a data.frame from the data in the list.
##' @export
qtrStore <- function(indata, timeCol, dataCol, logical = FALSE){
    qtr <- zoo::as.yearqtr(c("2009 Q4",
                             "2010 Q1", "2010 Q2", "2010 Q3", "2010 Q4",
                             "2011 Q1", "2011 Q2", "2011 Q3", "2011 Q4",
                             "2012 Q1", "2012 Q2", "2012 Q3", "2012 Q4")
                           )

    ##nodes <- as.list(unique(indata$node))
    nodes <- unique(indata$node)
    df.list <- list()
    for (node in nodes)
        df.list[[length(df.list)+1]] <- qtrHelp(node, indata, timeCol, dataCol, logical, qtr)

    df <- data.table::rbindlist(df.list)
    return(df)
}

qtrHelp <-  function(node, indata, timeCol, dataCol, logical = FALSE, qtr){
    dat <- indata[indata$node == node,]
    dataqtr <- zoo::as.yearqtr(as.Date(dat[,timeCol], origin = "2005-01-01"))

    outside <- which(!(qtr %in% dataqtr))
    inside <- which(qtr %in% dataqtr)

    status <- numeric(length(qtr))
    numMeas <- numeric(length(qtr))


    for(i in inside){
        k <- which(dataqtr == qtr[i])
        ## if(dataCol == "I")
        ##     status[i] <- mean(dat[k,dataCol])
        ## else
        if(logical)
            status[i] <- as.logical(sum(dat[k,dataCol]))
        else
            status[i] <- sum(dat[k,dataCol])
        numMeas[i] <- length(dat[k,dataCol])
    }

    status[outside] <- NA

    df <- data.frame(node = rep(node, length(qtr)),
                     time = qtr,
                     sample = as.numeric(status),
                     number = numMeas)
    names(df) <- c("node", timeCol, dataCol, "number")
    return(df)
}

##' Compute autoregression constants, mean, and log of variance
##' @param vec a named list with lists of vectors to extract summary statistics from
##' @param w weigts
##' @export
wilkinsonStatistics <- function(vec, w = NULL){
    if (is.list(vec)) {
        ac23 <- sapply(vec, function(x){acf(x, lag.max = 2, plot = FALSE)$acf[2:3]})
        SS <- c(ac23, sapply(vec,mean), sapply(vec,function(x){log(var(x) + 1)}))
    } else {
        ## this one is not from Wilkinson
        if(is.null(w)) {
            SS <- c()
            ## mean statistics
            ms <- halfyearMean(vec)
            SS <- c(SS,ms)


            ## Good for Upsilon
            ## fft statistics
            w <- rep(1,length(vec))/length(vec)
            if("fftcoeff" %in% names(extra)) {
                fftcoef <- extra$fftcoef
                fs <- fftAbs(vec,w,fftcoef)
                SS <- c(SS,fs)
            }

            ##ac23 <- acf(vec, lag.max = 2, plot = FALSE)$acf[c(2,3)]
            ##SS <-  c(ac23, mean(vec), log(var(vec) + 1))
        } else {
            SS <- c()

            ## quarter statistics
            qs <- qtrStatistics(vec,w)
            SS <- c(SS,qs)

            ## Good for Upsilon
            ## fft statistics
            w <- rep(1,length(vec))/length(vec)
            if("fftcoeff" %in% names(extra)) {
                fftcoef <- extra$fftcoef
                fs <- fftAbs(vec,w,fftcoef)
                SS <- c(SS,fs)
            }
        }

    }
    return(SS)
}

##' Compute the two period means.
##' @param vec summarized data in vector (as sum or mean from matrix)
##' @export
halfyearMean <- function(vec){
    ## we have two periods, t1 and t2.
    times <- as.numeric(names(vec))
    ## make into yearqtr
    times.qt <- zoo::as.yearqtr(as.Date(times, origin = "2005-01-01"))

    qvec <- c("Q1","Q2","Q3","Q4")
    qmean <- numeric(length(qvec))
    for (qi in seq_len(length(qvec))) {
        q <- qvec[qi]
        ## which are in quarter q?
        in.q <- grep(q,times.qt)
        vec.q <- vec[in.q]
        ## compute sum
        qmean[qi] <- mean(vec.q)
    }

    ## Q1 and Q2
    mean.t1 <- (qmean[1] + qmean[2])/2
    ## Q3 and Q4
    mean.t2 <- (qmean[3] + qmean[4])/2

    return(c(m1 = mean.t1, m2 = mean.t2))
}

##' Compute fft statistics (weighted)
##' @param vec summarized data in vector (as sum or mean from matrix)
##' @param w the weights
##' @param numCoef the number of coefficients to save.
##' @export
fftAbs <- function(vec,w,numCoef=c(4,5)) {
    ## Weight the data
    y <- vec*w

    ## transform the data
    y.fft <- fft(y)

    z <- abs(y.fft)[numCoef]

    ## find the numCoef largest coefficients (amplitudes)
    ## ampl.uni <- unique(abs(y.fft))
    ## ampl.ord <- order(ampl.uni,decreasing=TRUE)
    ## ampl.larg <- ampl.uni[ampl.ord[1:numCoef]]

    ## name the output
    name <- as.character(numCoef)
    for (i in 1:length(numCoef))
        name[i] <- paste("fft_coef", numCoef[i], sep="")

    ##names(ampl.larg) <- name
    names(z) <- name
    ##return(ampl.larg)
    return(z)
}

##' Compute mean and variance (weighted)
##' @param vec summarized data in vector (as sum or mean from matrix)
##' @param w the weights
##' @export
qtrStatistics <- function(vec, w) {
    qt <- zoo::as.yearqtr(as.numeric(names(vec)))

    ## one for each quarter
    valsMean <- numeric(4)
    ##valsVar <- numeric(4)
    qvec <- c("Q1","Q2","Q3","Q4")
    for (qi in seq_len(length(qvec))) {
        q <- qvec[qi]
        ## which are in quarter q?
        vec.q <- vec[grep(q,qt)]
        w.q <- w[grep(q,qt)]
        ## compute the weigted mean
        valsMean[qi] <- Weighted.Desc.Stat::w.mean(vec.q,w.q)

    }

    names(valsMean) = sapply(qvec,function(x){paste("mean",x,sep="")})
    return(valsMean)
}


##' Compute mean and variance (weighted)
##' @param vec summarized data in vector (as sum or mean from matrix)
##' @param w the weights
##' @export
acStatistics <- function(vec,w) {
                                        # weight the vector.
    vec.w <- vec*w

    # for looping
    len <- length(vec)
    val1 <- 0
    val2 <- 0


    for (i in seq(1,len-3,3)) {
        ## rolling windows
        ac23 <- acf(vec.w[i:(i+3)], lag.max = 2, plot = FALSE)$acf[c(2,3)]
        val1 <- val1 + ac23[1]
        val2 <- val2 + ac23[2]
    }

    ## compute the mean
    ac2 <- val1/floor(len/3)
    ac3 <- val2/floor(len/3)

    return(c("ac2.r" = ac2, "ac3.r" = ac3))
}

##' Compute mean and variance (weighted)
##' @param vec summarized data in vector (as sum or mean from matrix)
##' @param wmean the weighted mean
##' @export
countStatistics <- function(vec, wmean) {
    aboveMean <- length(which(vec > wmean))
    belowMean <- length(which(vec < wmean))
    return(c("above" = aboveMean, "below" = belowMean))
}


##' Compute the weights used in the weighted
##' @param data Observation data in matrix form
##' @export
computeWeights <- function(data){
    w <- apply(data, 1, sum)
    ## normalize the weights
    w <- w/sum(w)
    return(w)
}

##' Compute mean and variance (weighted)
##' @param vec summarized data in vector (as sum or mean from matrix)
##' @param w the weights
##' @export
weightedStat <- function(vec, w){
    wm <- weighted.mean(vec, w)
    wvar <- weighted.var(vec, w)

    ##return(c(wm, wvar))
    ## mean
    wm <- Weighted.Desc.Stat::w.mean(vec,w)
    ## var
    wvar <- Weighted.Desc.Stat::w.var(vec,w)
    ## kurtosis
    wk <- Weighted.Desc.Stat::w.kurtosis(vec,w)
    ## skewness
    ws <- Weighted.Desc.Stat::w.skewness(vec,w)

    return(c("mean" = wm, "logvar" = log(wvar+1), "logkurt" = log(wk+2), "skew" = ws))
}

##' Compute mean and variance (weighted)
##' @param vec summarized data in vector (as sum or mean from matrix)
##' @param w the weights
##' @export
lagStat <- function(vec,w) {
    vec.w <- vec*w
    vec.0 <- vec.w[1:(length(vec)-1)]
    vec.1 <- vec.w[2:length(vec)]

    lagdiff <- vec.0 - vec.1

    lagmean <- mean(lagdiff)
    ##lagvar <- var(lagdiff)
    return(c("lagmean" = lagmean))
}


##' Weighted Variance
##' Credit to Dr. Gavin Simpson +44 (0)20 7679 0522
##' @export
weighted.var <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2,
                                       na.rm = na.rm)
}

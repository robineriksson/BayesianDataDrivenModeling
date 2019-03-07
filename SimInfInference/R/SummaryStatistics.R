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
    w <- aggSS.all$w

    sumSS.mat <- matrix(aggSS, ncol = length(data))
    if(bs && length(data) != 1)
        sumSS.mat <- bootstrap(sumSS.mat, extraArgs)
    if(!useW)
        w <- NULL
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


    for(l in 1:dim(mat2)[2]) {
        lines(t, mat2[,l], col = "blue")
    }

     for(l in 1:dim(mat1)[2]) {
        lines(t, mat1[,l], col = "red")
    }
    ## df1 <- data.frame(t = t, mat1)
    ## df1.m <- reshape::melt(df1, id.var = "t")

    ## df2 <- data.frame(t = t, mat2)
    ## df2.m <- reshape::melt(df2, id.var = "t")

    ## p <- ggplot2::ggplot(df1.m, ggplot2::aes(x = t, y = value, col = "red")) +
    ##     ggplot2::geom_line() +
    ##     ggplot2::geom_line(data = df2.m, mapping = ggplot2::aes(x = t, y = value, colour = "blue"))

    ## p



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
        agged <- aggregateHelper(column = column, data = data, fun = fun, qtr = qtr)
                   ##   function(column, data){
                   ##       nodes <- unique(data$node)
                   ##       if(qtr) {
                   ##           ## aggregate into qtrs first
                   ##           data.qtrs <- qtrStore(data, "time", column, logical)
                   ##           data.mat <- matrix(data.qtrs[,column], ncol = length(nodes))
                   ##       } else {
                   ##           data.mat <- matrix(data[order(data$node),column], ncol = length(nodes))
                   ##       }

                   ##       data.fun <- apply(data.mat, 1, function(x){fun(x, na.rm = TRUE)})
                   ##   }
        ## , data = data)

        agged.unlist <- unlist(agged)
        data.frame(data = agged.unlist[grep("data", names(agged.unlist))],
                   w = agged.unlist[grep("w", names(agged.unlist))])
    })

    ## same w for all
    w <- agged[[1]]$w
    data <- c(sapply(agged, function(x){x$data}))
    return(list(data = data, w = w))
}


##' A helper function
##' @export
aggregateHelper <- function(column, data, fun, qtr){
    nodes <- unique(data$node)
    if(qtr) {
        ## aggregate into qtrs first
        data.qtrs <- qtrStore(data, "time", column)
        data.mat <- matrix(data.qtrs[,column], ncol = length(nodes))

        w <- computeWeights(data.mat, column)

    } else {

        data.mat <- matrix(data[order(data$node),column], ncol = length(nodes))
        w <- rep(NULL, dim(data.mat)[2])
    }

    data.fun <- apply(data.mat, 1, function(x){fun(x, na.rm = TRUE)})

    return(list(data = data.fun, w = w))
}


##' Store data in quarters.
##'
##' @param data as list
##' @return a data.frame from the data in the list.
##' @export
qtrStore <- function(indata, timeCol, dataCol, logical = FALSE){
    ## potential outliers are removed!
    ## qtr <- zoo::as.yearqtr(c("2009 Q4",
    ##                          ##"2010 Q1",
    ##                          "2010 Q2",
    ##                          ##"2010 Q3",
    ##                          ##"2010 Q4",
    ##                          "2011 Q1", "2011 Q2", "2011 Q3", "2011 Q4",
    ##                          ##"2012 Q1",
    ##                          ##"2012 Q2",
    ##                          "2012 Q3", "2012 Q4")
    ##                        )
    qtr <- zoo::as.yearqtr(c("2009 Q4",
                             "2010 Q1", "2010 Q2", "2010 Q3", "2010 Q4",
                             "2011 Q1", "2011 Q2", "2011 Q3", "2011 Q4",
                             "2012 Q1", "2012 Q2", "2012 Q3", "2012 Q4")
                           )

    nodes <- as.list(unique(indata$node))
    df.list <- lapply(nodes, function(node){
        dat <- indata[indata$node == node,]
        dataqtr <- zoo::as.yearqtr(as.Date(dat[,timeCol], origin = "2005-01-01"))

        outside <- which(!(qtr %in% dataqtr))
        inside <- which(qtr %in% dataqtr)

        status <- numeric(length(qtr))

        for(i in inside){
            k <- which(dataqtr == qtr[i])
            ## if(dataCol == "I")
            ##     status[i] <- mean(dat[k,dataCol])
            ## else
            if(logical)
                status[i] <- as.logical(sum(dat[k,dataCol]))
            else
                status[i] <- sum(dat[k,dataCol])
        }

        status[outside] <- NA

        df <- data.frame(node = rep(node, length(qtr)),
                         time = qtr,
                         sample = as.numeric(status))
        names(df) <- c("node", timeCol, dataCol)
        return(df)
    })

    df <- data.table::rbindlist(df.list)
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
        ac23 <- acf(vec, lag.max = 2, plot = FALSE)$acf[2:3]
        ## this one is not from Wilkinson
        ##lcoef <- as.numeric(lm(vec~seq_len(length(vec)))$coef)
        if(is.null(w))
            SS <-  c(ac23, mean(vec), log(var(vec) + 1))
        else {
            ws <- weightedStat(vec,w)
            SS <- c(ac23,ws)
        }
    }
    return(SS)
}


##' Compute the weights used in the weighted
##' @param data Observation data in matrix form
##' @export
computeWeights <- function(data, column){
    ## compute weightes
    if(column == "I") {
        ## If we track #I we would have to write down the weights earlier. This
        ## alternative method gives us a hinge on how the should be.
        w <- apply(data, 1, function(x){length(x) - length(which(is.na(x)))})
    } else {
        w <- apply(data, 1, function(x){
            tab <- table(x)
            sums <- 0
            if("0" %in% names(tab))
                sums <- tab[["0"]]
            sums <- sums + sum(x, na.rm = TRUE)
            return(sums)
        })
    }

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

    return(c(wm, log(wvar + 1)))
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

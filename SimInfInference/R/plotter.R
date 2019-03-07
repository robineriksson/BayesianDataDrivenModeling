
##' Plot the posterior for all the parameters.
##'
##' @param infe the inference calss
##' @param N the number of values to plot.
##' @param by "thin" out the chain by every by element
##' @return a ggpairs plot object
##' @export
posteriorPlot <- function(infe, n = 1, N = NULL, by = NULL, ordered = FALSE){
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

    myDensity <- function(data, mapping, ...) {
        g <- ggplot2::ggplot(data = data, mapping = mapping) +
            ggplot2::stat_density_2d(geom = "contour", bins = 5,
                                     ggplot2::aes(color = ..level..))
    }

    if(dim(data)[2] == 5)
        columnLabels = c("upsilon", "beta[1]", "beta[2]", "gamma")
    else if(dim(data)[2] == 6)
        columnLabels = c("upsilon", "beta[1]", "beta[2]", "beta[3]", "gamma")
    else
        columnLabels = names(data)[1:(dims[2]-1)]

    g <- GGally::ggpairs(data = data[,seq(1,dims[2]-1)],
                         columnLabels = columnLabels,
                         labeller = "label_parsed",
                         upper = list(continuous = myDensity),
                         lower = list(continuous = GGally::wrap("points", alpha = 0.6))
                         ) +
        theme_Publication()
    g
}


##' theme:ing for ggplots
##' Credit goes to:
##' Koundinya Desiraju (04/07/2015)
##' @export
theme_Publication <- function(base_size=14) {
      (ggthemes::theme_foundation(base_size=base_size)
       + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                         size = ggplot2::rel(1.2), hjust = 0.5),
                        text = ggplot2::element_text(),
                        panel.background = ggplot2::element_rect(colour = NA),
                        plot.background = ggplot2::element_rect(colour = NA),
                        panel.border = ggplot2::element_rect(colour = NA),
                        axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                        axis.title.y = ggplot2::element_text(angle= 90,vjust =2),
                        axis.title.x = ggplot2::element_text(vjust = -0.2),
                        axis.text = ggplot2::element_text(),
                        axis.text.x = ggplot2::element_text(angle = 270, hjust = 1),
                        axis.line = ggplot2::element_line(colour="black"),
                        axis.ticks = ggplot2::element_line(),
                        panel.grid.major = ggplot2::element_line(colour="#f0f0f0"),
                        panel.grid.minor = ggplot2::element_blank(),
                        legend.key = ggplot2::element_rect(colour = NA),
                        legend.position = "bottom",
                        legend.direction = "horizontal",
                        legend.key.size= grid::unit(0.2, "cm"),
                        legend.spacing = grid::unit(0, "cm"),
                        legend.title = ggplot2::element_text(face="italic"),
                        plot.margin=grid::unit(c(10,5,5,5),"mm"),
                        strip.background= ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                        strip.text = ggplot2::element_text(face="bold")
                        ))

}

##' Plot mean aggregated data
##'
##' Plot aggregated the data that is used when computing SS
##' @param data distance matrix data
##' @param column which column to "aggregate"
##' @return the aggregates (mean) time series.
##' @export
plotAggregate_mean<- function(data, column = "sample"){
    nodes <- unique(data$node)
    data.mat <- matrix(data[order(data$node),column], ncol = length(nodes))
    data.mean <- apply(data.mat, 1, function(x){mean(x, na.rm = TRUE)})
    plot(data.mean, type = "b")
    return(data.mean)
}

##' Plot mean aggregated data
##'
##' Plot aggregated the data that is used when computing SS
##' @param data distance matrix data
##' @param column which column to "aggregate"
##' @return the aggregates (mean) time series.
##' @export
plotAggregate_sum<- function(data, column = "sample"){
    nodes <- unique(data$node)
    data.mat <- matrix(data[order(data$node),column], ncol = length(nodes))
    data.mean <- apply(data.mat, 1, function(x){sum(x, na.rm = TRUE)})
    plot(data.mean, type = "b")
    return(data.mean)
}

##' Plot mean aggregated data
##'
##' Plot aggregated the data that is used when computing SS
##' @param data distance matrix data
##' @param column which column to "aggregate"
##' @return the aggregates (mean) time series.
##' @export
plotAggregate_sum_qtr<- function(listOfData, column = "sample"){
    y <- sapply(listOfData, function(data){
        data.qtr <- qtrStore(data, "time", column)
        nodes <- unique(data.qtr$node)
        data.qtr.mat <- matrix(data.qtr[order(data.qtr$node),column], ncol = length(nodes))
        data.qtr.sum <- apply(data.qtr.mat, 1, function(x){sum(x, na.rm = TRUE)})
    })

    plot(y, type = "b")
    return(y)
}

##' plot posterior in 3d
##'
##' plot
##' @param data distance matrix data
##' @param invert invert the distance/likelihood
##' @param ex exponentiate the distance/likelihood
##' @export
plot3dStuff <-  function(data, invert = FALSE, ex = TRUE){
    ##require(akima); require(rgl); require(plot3D);
    x = data[,1]
    y = data[,2]
     if(ex)
        data[,3] <- exp(data[,3])
    if(invert){
        z = 1/data[,3]
        if(Inf %in% z)
            z = 1/(1+data[,3])
        zlab <- "1/distance"
    } else {
        z = data[,3]
        zlab <- "distance"
    }
    s=akima::interp(x,y,z)
    plot3D::persp3D(s$x,s$y,s$z,
            xlab = expression(upsilon), ylab = expression(beta[1]), zlab = zlab)



}

##' plot 2d lattice
##'
##' plot
##' @param data distance matrix data
##' @export
plot2dStuff <- function(data){
    lattice::contourplot(dist ~ upsilon*beta, data = data, region = TRUE)
}

##' plot posterior in 2d, with interpolated values.
##'
##' plot
##' @param data distance matrix data
##' @param invert invert the distance/likelihood
##' @param ex exponentiate the distance/likelihood
##' @export
plot2dotherStuff <-  function(data, type = "exp", theta0 = NULL, rec = c(0.0095, 0.090, 0.0105, 0.10)){
       ##require(akima);require(plot3D);
    cols <- length(data)
    if(cols == 3){
        x = data[,1]
        y = data[,2]
        if(ex)
            data[,3] <- exp(data[,3])
        if(invert){
            z = 1/data[,3]
            if(Inf %in% z)
                z = 1/(1+data[,3])
        } else
            z = data[,3]

        s=akima::interp(x,y,z)
        plot3D::image2D(x = s$x, y = s$y, z = s$z,
                        contoure = list(levels = 1, col = "black", lwd = 2),
                        xlab = expression(upsilon), ylab = expression(beta[1]), main = "1/distance")
        rect(0.0095, 0.090, 0.0105, 0.10)
    } else if(cols > 3){
        cb <- combn(cols-1,2)
        cb.len <- dim(cb)[2]
        print(cb.len)
        par(mfrow = c(ceiling(cb.len/2),2), mar = c(4,7,4,7))

        for (i in seq_len(cb.len)){
            static <- !(unique(as.numeric(cb)) %in% cb[,i])


            rows.mat <- round(data[,static], digits = 7) == round(theta0[static], digits = 7)
            rows <- which(apply(rows.mat, 1, sum) > 1)

            x = data[rows,cb[1,i]]
            y = data[rows,cb[2,i]]
            if(type == "exp")
                z <- exp(data[rows, cols])
            else if(type == "invert"){
                z = 1/data[rows,cols]
                if(Inf %in% z)
                    z = 1/(1+data[rows,cols])
            } else
                z = data[rows,cols] ## last column is always the "distance"

            s=akima::interp(x,y,z)
            plot3D::image2D(x = s$x, y = s$y, z = s$z,
                            contoure = list(levels = 1, col = "black", lwd = 2),
                            xlab = names(theta0[cb[1,i]]), ylab = names(theta0[cb[2,i]]), main = "1/distance")
            ##rect(rec[1],rec[2],rec[3],rec[4])



        }
    }
}

##' plot posterior in 2d, with interpolated values.
##'
##' plot
##' @param data distance matrix data
##' @param invert invert the distance/likelihood
##' @param ex exponentiate the distance/likelihood
##' @export
plot2Class <-  function(inf, type = "exp", theta0 = NULL, rec = c(0.0095, 0.090, 0.0105, 0.10)){
    data <- inf$getPosterior()
    theta0 <- inf$getExtraArgsEstimator()$theta0

    ##require(akima);require(plot3D);
    cols <- length(data)
    if(cols == 3){
        x = data[,1]
        y = data[,2]
        if(ex)
            data[,3] <- exp(data[,3])
        if(invert){
            z = 1/data[,3]
            if(Inf %in% z)
                z = 1/(1+data[,3])
        } else
            z = data[,3]

        s=akima::interp(x,y,z)
        plot3D::image2D(x = s$x, y = s$y, z = s$z,
                        contoure = list(levels = 1, col = "black", lwd = 2),
                        xlab = expression(upsilon), ylab = expression(beta[1]), main = "1/distance")
        rect(0.0095, 0.090, 0.0105, 0.10)
    } else if(cols > 3){
        cb <- combn(cols-1,2)
        cb.len <- dim(cb)[2]
        print(cb.len)
        par(mfrow = c(ceiling(cb.len/2),2), mar = c(4,7,4,7))

        for (i in seq_len(cb.len)){
            static <- !(unique(as.numeric(cb)) %in% cb[,i])


            rows.mat <- round(data[,static], digits = 7) == round(theta0[static], digits = 7)
            rows <- which(apply(rows.mat, 1, sum) > 1)

            x = data[rows,cb[1,i]]
            y = data[rows,cb[2,i]]
            if(type == "exp")
                z <- exp(data[rows, cols])
            else if(type == "invert"){
                z = 1/data[rows,cols]
                if(Inf %in% z)
                    z = 1/(1+data[rows,cols])
            } else
                z = data[rows,cols] ## last column is always the "distance"

            s=akima::interp(x,y,z)
            plot3D::image2D(x = s$x, y = s$y, z = s$z,
                            contoure = list(levels = 1, col = "black", lwd = 2),
                            xlab = names(theta0[cb[1,i]]), ylab = names(theta0[cb[2,i]]), main = "1/distance")
            ##rect(rec[1],rec[2],rec[3],rec[4])



        }
    }
}

##' plot posterior in 2d and 3d side by side
##'
##' plot
##' @param data distance matrix data
##' @param invert invert the distance/likelihood
##' @param ex exponentiate the distance/likelihood
##' @export
plot2and3 <-  function(data, invert = TRUE, ex = FALSE){
    ##require(akima); require(rgl); require(plot3D);

    x = data[,1]
    y = data[,2]

    if(ex)
        data[,3] <- exp(data[,3])

    if(invert){
        z = 1/data[,3]
        if(Inf %in% z)
            z = 1/(1+data[,3])
        zlab <- "1/distance"
    } else {
        z = data[,3]
        zlab <- "distance"
    }

    s=akima::interp(x,y,z)
    ##X11()
    par(mfrow = c(1,2), mar = c(10,1,10,1))
    plot3D::image2D(x = s$x, y = s$y, z = s$z,
            contoure = list(levels = 1, col = "black", lwd = 2),
            xlab = expression(upsilon), ylab = expression(beta[1]), main = "1/distance")
    rect(0.0095, 0.090, 0.0105, 0.10)
    plot3D::persp3D(s$x,s$y,s$z,
            xlab = expression(upsilon), ylab = expression(beta[1]), zlab = zlab)

}

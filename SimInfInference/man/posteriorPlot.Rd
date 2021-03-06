% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotter.R
\name{posteriorPlot}
\alias{posteriorPlot}
\title{Plot the posterior for all the parameters.}
\usage{
posteriorPlot(infe, n = 1, N = NULL, by = NULL, ordered = FALSE)
}
\arguments{
\item{infe}{the inference calss}

\item{N}{the number of values to plot.}

\item{by}{"thin" out the chain by every by element}
}
\value{
a ggpairs plot object
}
\description{
Plot the posterior for all the parameters.
}

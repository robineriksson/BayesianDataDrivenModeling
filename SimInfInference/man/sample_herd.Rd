% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nodeSampler.R
\name{sample_herd}
\alias{sample_herd}
\title{Sample a herd (node) to to determine if pos or false.}
\usage{
sample_herd(S, I)
}
\arguments{
\item{S}{scalar with Susceptible}

\item{I}{scalar with Infected}

\item{node}{The index of the herd to sample}

\item{week}{The index of the week to sample}
}
\value{
TRUE or FALSE
}
\description{
Sample a herd (node) to to determine if pos or false.
}

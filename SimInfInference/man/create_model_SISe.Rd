% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initmodel.R
\name{create_model_SISe}
\alias{create_model_SISe}
\title{Create the model before running a trajectory}
\usage{
create_model_SISe(theta, u0, distance, tspan, events)
}
\arguments{
\item{theta}{(named) model parameters {up1,up2,up3,beta1,beta2,beta3,beta4}}

\item{u0}{initial infection data}

\item{tspan}{time (integer) vector used in simulator}

\item{dist}{distance matrix connecting nodes}
}
\value{
an initialised model
}
\description{
Create the model before running a trajectory
}

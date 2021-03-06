% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/full_models.R
\name{tsc_objectivea}
\alias{tsc_objectivea}
\title{Objective function for full of tscopula plus margin model}
\usage{
tsc_objectivea(theta, modelspec, modeltype, dens, cdf, y)
}
\arguments{
\item{theta}{vector of parameter values}

\item{modelspec}{list containing model specification}

\item{modeltype}{character string giving type of model}

\item{dens}{marginal density function}

\item{cdf}{marginal cdf}

\item{y}{vector of data values}
}
\value{
Value of objective function at parameters.
}
\description{
Objective function for full of tscopula plus margin model
}
\keyword{internal}

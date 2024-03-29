% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sarmacopula.R
\name{sarmacopula_objective}
\alias{sarmacopula_objective}
\title{Objective function for SARMA copula process}
\usage{
sarmacopula_objective(theta, modelspec, u)
}
\arguments{
\item{theta}{vector of parameters of DARMA process}

\item{modelspec}{vector containing model order (p,q)(P,Q)D}

\item{u}{vector of data}
}
\value{
Value of objective function at parameters.
}
\description{
Objective function for SARMA copula process
}
\keyword{internal}

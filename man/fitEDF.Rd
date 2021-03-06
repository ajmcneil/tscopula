% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/full_models.R
\name{fitEDF}
\alias{fitEDF}
\title{Fit tscm using empirical distribution function}
\usage{
fitEDF(x, y, tsoptions, control)
}
\arguments{
\item{x}{an object of class \linkS4class{tscm}.}

\item{y}{a vector or time series of data.}

\item{tsoptions}{a list of parameters passed to fitting.
This differs according to the class of x.}

\item{control}{list of control parameters to be passed to the
\code{\link[stats]{optim}} function.}
}
\value{
An object of class \linkS4class{tscmfit}.
}
\description{
Fit tscm using empirical distribution function
}
\keyword{internal}

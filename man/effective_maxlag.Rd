% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dvinecopula3.R
\name{effective_maxlag}
\alias{effective_maxlag}
\title{Find effective maximum lag}
\usage{
effective_maxlag(tau, tautol, maxlag)
}
\arguments{
\item{tau}{values of Kendall partial autocorrelation function.}

\item{tautol}{value at which kpacf should be truncated.}

\item{maxlag}{manual override of maximum lag.}
}
\value{
a value for the effective maximum lag floored at one.
}
\description{
Find effective maximum lag
}
\keyword{internal}

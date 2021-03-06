% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vtransforms.R
\name{vgradient}
\alias{vgradient}
\title{Calculate gradient of v-transform}
\usage{
vgradient(x, u)
}
\arguments{
\item{x}{an object of class \linkS4class{Vtransform}.}

\item{u}{a vector or time series with values in [0, 1].}
}
\value{
A vector or time series of values of gradient.
}
\description{
Calculate gradient of v-transform
}
\examples{
vgradient(Vsymmetric(), c(0, 0.25, 0.5, 0.75, 1))
}

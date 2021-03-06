% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dvinecopula.R
\name{dvinecopula}
\alias{dvinecopula}
\title{Constructor function for dvinecopula process}
\usage{
dvinecopula(family = "indep", pars = list(NULL), rotation = 0)
}
\arguments{
\item{family}{a vector of family names}

\item{pars}{a list containing the parameters of each lag}

\item{rotation}{a vector of rotations}
}
\value{
An object of class \linkS4class{dvinecopula}.
}
\description{
Constructor function for dvinecopula process
}
\examples{
dvinecopula(family = c("joe", "gauss", "t"), pars = list(3, .5, c(1, 2)), rotation = c(180, 0, 0))
}

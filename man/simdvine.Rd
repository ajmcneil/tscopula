% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dvinecopula.R
\name{simdvine}
\alias{simdvine}
\title{D-vine simulation helper function}
\usage{
simdvine(pc_list, n, innov, start)
}
\arguments{
\item{pc_list}{a list of pair copulas.}

\item{n}{number of data to be simulated.}

\item{innov}{vector of innovations of length n.}

\item{start}{vector of start values with length equal to order of process.}
}
\value{
a vector of length n.
}
\description{
D-vine simulation helper function
}
\keyword{internal}

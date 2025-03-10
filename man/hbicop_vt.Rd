% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dvinecopulavt.R
\name{hbicop_vt}
\alias{hbicop_vt}
\title{h-function for linear inverse-v-transformed copula}
\usage{
hbicop_vt(u, cond_var, family, vt1, vt2, inverse = FALSE)
}
\arguments{
\item{u}{two-column matrix of data at which h-function is evaluated}

\item{cond_var}{identity of conditioning variable (1/2)}

\item{family}{name of copula family}

\item{vt1}{first v-transform}

\item{vt2}{second v-transform}

\item{inverse}{logical variable specifying whether inverse is taken}
}
\value{
vector of values of h-function
}
\description{
h-function for linear inverse-v-transformed copula
}
\keyword{internal}

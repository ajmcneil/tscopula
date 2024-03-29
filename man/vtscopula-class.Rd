% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting_vtscopula.R
\docType{class}
\name{vtscopula-class}
\alias{vtscopula-class}
\alias{show,vtscopula-method}
\alias{coef,vtscopula-method}
\alias{predict,vtscopula-method}
\alias{sim,vtscopula-method}
\alias{kendall,vtscopula-method}
\title{Time series copula processes with v-transforms}
\usage{
\S4method{show}{vtscopula}(object)

\S4method{coef}{vtscopula}(object)

\S4method{predict}{vtscopula}(object, data, x, type = "df")

\S4method{sim}{vtscopula}(object, n = 1000)

\S4method{kendall}{vtscopula}(object, lagmax = 20)
}
\arguments{
\item{object}{an object of the class.}

\item{data}{vector of past data values.}

\item{x}{vector of arguments of prediction function.}

\item{type}{type of prediction function ("df" for density, "qf" for quantile function
or "dens" for density).}

\item{n}{length of realization.}

\item{lagmax}{maximum value of lag.}
}
\description{
Class of objects for v-transformed time series copula processes.
}
\section{Methods (by generic)}{
\itemize{
\item \code{show(vtscopula)}: Show method for vtscopula objects

\item \code{coef(vtscopula)}: Coef method for vtscopula class

\item \code{predict(vtscopula)}: Prediction method for vtscopula class

\item \code{sim(vtscopula)}: Simulation method for vtscopula class

\item \code{kendall(vtscopula)}: Calculate Kendall's tau values for vtscopula model

}}
\section{Slots}{

\describe{
\item{\code{Vcopula}}{object of class \linkS4class{tscopulaU}.}

\item{\code{Vtransform}}{object of class \linkS4class{Vtransform}.}

\item{\code{Wcopula}}{object of class \linkS4class{tscopula}.}
}}

\examples{
copobject <- armacopula(pars = list(ar = 0.6, ma = 0.2))
sim(vtscopula(copobject, Vtransform = V2p()))
mod <- vtscopula(armacopula(list(ar = 0.95, ma = -0.85)))
kendall(mod)
}

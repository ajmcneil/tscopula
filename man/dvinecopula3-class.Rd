% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dvinecopula3.R
\docType{class}
\name{dvinecopula3-class}
\alias{dvinecopula3-class}
\alias{coef,dvinecopula3-method}
\alias{kendall,dvinecopula3-method}
\alias{show,dvinecopula3-method}
\alias{sim,dvinecopula3-method}
\alias{predict,dvinecopula3-method}
\title{D-vine copula processes of type 3}
\usage{
\S4method{coef}{dvinecopula3}(object)

\S4method{kendall}{dvinecopula3}(object, lagmax = 20)

\S4method{show}{dvinecopula3}(object)

\S4method{sim}{dvinecopula3}(object, n = 1000)

\S4method{predict}{dvinecopula3}(object, data, x, type = "df")
}
\arguments{
\item{object}{an object of the class.}

\item{lagmax}{maximum value of lag to be considered.}

\item{n}{length of realization.}

\item{data}{vector of past data values.}

\item{x}{vector of arguments of prediction function.}

\item{type}{type of prediction function ("df" for density, "qf" for quantile function
or "dens" for density).}
}
\description{
Class of objects for d-vine copula processes. See \link{dvinecopula3} for more details.
}
\section{Methods (by generic)}{
\itemize{
\item \code{coef(dvinecopula3)}: Coef Method for dvinecopula3 class

\item \code{kendall(dvinecopula3)}: Calculate Kendall's tau values for pair copulas in type 3 d-vine copula

\item \code{show(dvinecopula3)}: Show method for dvinecopula3 class

\item \code{sim(dvinecopula3)}: Simulation method for dvinecopula3 class

\item \code{predict(dvinecopula3)}: Prediction method for dvinecopula2 class

}}
\section{Slots}{

\describe{
\item{\code{name}}{name of the d-vine copula process.}

\item{\code{modelspec}}{list containing the family, rotation, and name of KPACF}

\item{\code{pars}}{list comprising of the parameters.}
}}


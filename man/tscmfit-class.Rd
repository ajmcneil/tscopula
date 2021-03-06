% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/full_models.R
\docType{class}
\name{tscmfit-class}
\alias{tscmfit-class}
\alias{logLik,tscmfit-method}
\alias{resid,tscmfit-method}
\title{Fitted tscm model}
\usage{
\S4method{logLik}{tscmfit}(object)

\S4method{resid}{tscmfit}(object, trace = FALSE)
}
\arguments{
\item{object}{an object of the class.}

\item{trace}{extract trace instead of residuals.}
}
\description{
Class of objects for fitted \linkS4class{tscm} models.
}
\section{Methods (by generic)}{
\itemize{
\item \code{logLik}: method for tscmfit class

\item \code{resid}: Residual method for tscmfit class
}}

\section{Slots}{

\describe{
\item{\code{tscopula}}{an object of class \linkS4class{tscopula}.}

\item{\code{margin}}{an object of class \linkS4class{margin}.}

\item{\code{data}}{a vector or time series of data to which process has been fitted.}

\item{\code{fit}}{a list containing details of the fit.}
}}


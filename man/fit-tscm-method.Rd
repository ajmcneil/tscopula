% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/full_models.R
\name{fit,tscm-method}
\alias{fit,tscm-method}
\title{Fit Method for tscm Class}
\usage{
\S4method{fit}{tscm}(
  x,
  y,
  tsoptions = list(),
  control = list(warn.1d.NelderMead = FALSE, trace = FALSE, maxit = 5000),
  method = "IFM"
)
}
\arguments{
\item{x}{an object of class \linkS4class{tscm}.}

\item{y}{a vector or time series of data.}

\item{tsoptions}{a list of parameters passed to fitting.}

\item{control}{list of control parameters to be passed to the
\code{\link[stats]{optim}} function.}

\item{method}{character string specifying method.}
}
\value{
an object of class \linkS4class{tscmfit}.
}
\description{
Fit Method for tscm Class
}
\examples{
mod <- tscm(dvinecopula(family = "gauss", pars = 0.5), margin("weibull2"))
y <- sim(mod)
fit(mod, y)
}
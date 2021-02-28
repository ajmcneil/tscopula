#' Time Series Copula Processes
#'
#' Class of objects for time series copula processes.
#'
#'
#' @export
#'
#' @import stats
#' @import methods
#' @import graphics
#' @import utils
#'
setClass("tscopula", contains = c("VIRTUAL"))

#' New Generic for Simulating Time Series Models
#'
#' Methods are available for objects of class \linkS4class{tscopula},
#' \linkS4class{margin} and \linkS4class{tscm}.
#'
#' @param x an object of the model class.
#' @param ... further arguments to be passed on.
#'
#' @return A realization from the time series model.
#' @export
#'
#'
setGeneric("sim", function(x, ...) {
  standardGeneric("sim")
})

#' Strict White Noise Copula Process
#'
#' @export
#'
setClass("swncopula", contains = "tscopula")

#' Constructor Function for Strict White Noise Copula
#'
#' @return The strict white noise copula process.
#' @export
#'
#' @examples
#' swncopula()
swncopula <- function() {
  new("swncopula")
}

#' Simulation Method for Strict White Noise Copula
#'
#' @param x an object of class \linkS4class{swncopula}.
#' @param n numeric value for length of simulated realisation.
#'
#' @return A realisation of strict white noise of length n.
#' @export
#'
#' @examples
#' sim(swncopula())
setMethod("sim", "swncopula", function(x, n = 1000) {
  runif(n)
})

#' Coef Method for Strict White Noise Copula
#'
#' @param object an object of class \linkS4class{swncopula}.
#'
#' @return Null object
#' @export
#'
setMethod("coef", "swncopula", function(object) {
  NULL
})

#' Show Method for Strict White Noise Copula
#'
#' @param object an object of class \linkS4class{swncopula}.
#'
#' @return A summary of a \linkS4class{swncopula}.
#' @export
#'
setMethod("show", "swncopula", function(object) {
  cat("SWN \n")
})

#' New Generic for Kendall correlations
#'
#' Methods are available for objects of class \linkS4class{armacopula}
#'  and \linkS4class{dvinecopula}.
#'
#' @param x an object of the model class.
#' @param ... further arguments to be passed on.
#'
#' @return A realization from the time series model.
#' @export
#'
#'
setGeneric("kendall", function(x, ...) {
  standardGeneric("kendall")
})

#' New Generic for Generalized Lagging of Data
#'
#' Methods are available for objects of class \linkS4class{tscopulafit}
#'  and \linkS4class{tscmfit}.
#'
#' @param x an object of the model class.
#' @param ... further arguments to be passed on.
#'
#' @return A realization from the time series model.
#' @export
#'
#'
setGeneric("glag", function(x, ...) {
  standardGeneric("glag")
})



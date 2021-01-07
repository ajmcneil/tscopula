#' Time Series Copula Processes
#'
#' Class of objects for time series copula processes.
#'
#'
#' @return
#' @export
#'
setClass("tscopula", contains = c("VIRTUAL"))

#' New Generic for Simulating Time Series Models
#'
#' Methods are available for objects of class \linkS4class{tscopula},
#' \linkS4class{margin} and \linkS4class{tsc}.
#'
#' @param x an object of the model class.
#' @param ...
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
#' @return
#' @export
#'
setClass("swncopula", contains = "tscopula")

#' Constructor Function for Strict White Noise Copula
#'
#' @return the strict white noise copula process
#' @export
#'
#' @examples
#' swncopula()
swncopula <- function() {
  new("swncopula")
}

#' Simulation Method for Strict White Noise Copula
#'
#' @param swncopula an object of class \linkS4class{swncopula}.
#' @param n numeric value for length of simulated realisation.
#'
#' @return a realisation of strict white noise of length n
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
#' @return
#' @export
#'
setMethod("coef", "swncopula", function(object) {
  NULL
})

#' Show Method for Strict White Noise Copula
#'
#' @param object an object of class \linkS4class{swncopula}.
#'
#' @return
#' @export
#'
setMethod("show", "swncopula", function(object) {
  cat("SWN \n")
})







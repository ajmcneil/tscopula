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

#' New Generic for Estimating Time Series Models
#'
#' Methods are available for objects of class \linkS4class{tscopula},
#' \linkS4class{vtscopula}, \linkS4class{margin}
#' \linkS4class{tsc} and \linkS4class{gvtcopar}.
#'
#' @param x an object of the model class.
#' @param y a vector or time series of data.
#' @param ...
#'
#' @return An object of the fitted model class.
#' @export
#'
#'
setGeneric("fit", function(x, y, ...) {
  standardGeneric("fit")
})

#' New Generic for evaluating the joint density of a Time Series Models
#'
#' Methods are available for objects of class \linkS4class{tscopula},
#' \linkS4class{vtscopula}, \linkS4class{margin}
#' \linkS4class{tsc} and \linkS4class{gvtcopar}.
#'
#' @param x an object of the model class.
#' @param y a vector or time series of data.
#' @param ...
#'
#' @return joint density
#' @export
#'
#'
setGeneric("joint", function(x, y, ...) {
  standardGeneric("joint")
})

#' New Generic for carrying out a likelihood ratio test between two Time Series Models
#'
#' Methods are available for objects of class \linkS4class{tscopulafit},
#' \linkS4class{tscmfit}
#'
#' @param x,y objects of the model class.
#' @param ...
#'
#' @return a likelihood ratio test result.
#' @export
#'
#'
setGeneric("LRT", function(x, y, ...) {
  standardGeneric("LRT")
})

#' SARMA copula processes
#'
#' Class of objects for seasonal ARMA copula processes.
#'
#' @slot name name of seasonal ARMA copula process.
#' @slot modelspec vector containing number of AR, MA, SAR and SMA parameters as well
#' as the order D of seasonal differencing.
#' @slot pars list consisting of vector of AR parameters named `ar`
#' and vector of MA parameters named `ma`, SAR parameters named `sar`
#' and vector of SMA parameters named `sma`.
#'
#' @export
#'
setClass("sarmacopula", contains = "tscopula", slots = list(
  name = "character",
  modelspec = "numeric",
  pars = "list"
))

#' Constructor function for SARMA copula process
#'
#' @param pars list consisting of vector of AR parameters named `ar`
#' and vector of MA parameters named `ma`, SAR parameters named `sar`
#' and vector of SMA parameters named `sma`.
#' @param period period of seasonal model.
#'
#' @return An object of class \linkS4class{sarmacopula}.
#' @export
#'
#' @examples
#' sarmacopula(list(ar = 0.5, ma = 0.4, sar = 0.2, sma = 0.6), period = 4)
sarmacopula <- function(pars = list(ar = 0, ma = 0, sar = 0, sma = 0), period = 4) {
  p <- length(pars$ar)
  q <- length(pars$ma)
  P <- length(pars$sar)
  Q <- length(pars$sma)
  if (p > 0){
    if (non_stat(pars$ar))
      stop("Non-stationary AR part")
    names(pars$ar) <- paste("ar", 1:p, sep = "")
  }
  if (q > 0){
    if (non_invert(pars$ma))
      stop("Non-stationary MA part")
    names(pars$ma) <- paste("ma", 1:q, sep = "")
  }
  if (P > 0){
    if (non_stat(pars$sar))
      stop("Non-stationary SAR part")
    names(pars$sar) <- paste("sar", 1:P, sep = "")
  }
  if (Q > 0){
    if (non_invert(pars$sma))
      stop("Non-stationary SMA part")
    names(pars$sma) <- paste("sma", 1:Q, sep = "")
  }
  if ((p == 0) & (q == 0) & (P == 0) & (Q == 0)) {
    stop("Specify named ar, ma, sar or sma parameters")
  }
  new("sarmacopula",
      name = paste("SARMA(", p, ",", q, ")(", P, ",", Q,")", period, sep = ""),
      modelspec = c(p = p, q = q, P = P, Q = Q, period = period),
      pars = pars
  )
}

#' @describeIn sarmacopula Coef method for SARMA copula class
#'
#' @param object an object of the class.
#'
#' @export
setMethod("coef", "sarmacopula", function(object) {
  c(object@pars$ar, object@pars$ma, object@pars$sar, object@pars$sma)
})


#' Turn vector of SARMA parameters into list
#'
#' @param theta vector of SARMA model parameters
#' @param order order of model
#'
#' @return a list containing SARMA parameters in components ar, ma, sar and sma
#' @export
#'
sarmavec2list <- function(theta, order){
  p <- order[1]
  q <- order[2]
  P <- order[3]
  Q <- order[4]
  output <- list()
  if (p > 0)
    output$ar <- as.numeric(theta[1:p])
  if (q > 0)
    output$ma <- as.numeric(theta[(p + 1):(p + q)])
  if (P > 0)
    output$sar <- as.numeric(theta[(p + q + 1):(p + q + P)])
  if (Q > 0)
    output$sma <- as.numeric(theta[(p + q + P + 1):(p + q + P + Q)])
  output
}

#' @describeIn sarmacopula Show method for SARMA copula process
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", c(object = "sarmacopula"), function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("parameters: \n")
  print(coef(object))
})

#' Transform a sarmacopula object into an armacopula object
#'
#' @param object an object of class \linkS4class{sarmacopula}.
#'
#' @return An object of class \linkS4class{armacopula}.
#' @export
#'
#' @examples
#' sarma2arma(sarmacopula(list(ar = 0.5, ma = 0.4, sar = 0.2, sma = 0.6), period = 4))
sarma2arma <- function(object){
  if (!(is(object, "sarmacopula")))
    stop("Not sarmacopula object")
  period <- object@modelspec["period"]
  ar <- as.numeric(object@pars$ar) # Null to numeric(0) if necessary
  ma <- as.numeric(object@pars$ma)
  if (object@modelspec[3] > 0)
    ar <- expand_ar(ar, object@pars$sar, period)
  if (object@modelspec[4] > 0)
    ma <- expand_ma(ma, object@pars$sma, period)
  output <- list()
  if (length(ar) > 0)
    output$ar <- ar
  if (length(ma) > 0)
    output$ma <- ma
  armacopula(pars = output)
}

#' Expand AR coefficients to include SAR coefficients of SARMA model
#'
#' @param ar vector of AR coefficients
#' @param sar vector of SAR coefficients
#' @param period period of SARMA model
#'
#' @return vector of AR coefficients in equivalent ARMA model
#' @keywords internal
#'
expand_ar <- function(ar, sar, period){
  output <- rep(0, length(ar) + period*length(sar))
  sar <- as.vector(sapply(sar, function(x, d){c(rep(0, d-1), x)}, d = period))
  newvals <- -coefficients(polynom::polynomial(c(1, -sar)) *
                             polynom::polynomial(c(1, -ar)))[-1]
  if (length(newvals) > 0)
    output[1:length(newvals)] <- newvals
  output
}

#' Expand MA coefficients to include SMA coefficients of SARMA model
#'
#' @param ma vector of MA coefficients
#' @param sma vector of SMA coefficients
#' @param period period of SARMA model
#'
#' @return vector of MA coefficients in equivalent ARMA model
#' @keywords internal
#'
expand_ma <- function(ma, sma, period){
  output <- rep(0, length(ma) + period*length(sma))
  sma <- as.vector(sapply(sma, function(x, d){c(rep(0, d-1), x)}, d = period))
  newvals <- coefficients(polynom::polynomial(c(1, sma)) *
                            polynom::polynomial(c(1, ma)))[-1]
  if (length(newvals) > 0)
    output[1:length(newvals)] <- newvals
  output
}

#' @describeIn sarmacopula Simulation method for sarmacopula class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
#' @examples
#' sim(sarma2arma(sarmacopula(list(ar = 0.5, ma = 0.4, sar = 0.2, sma = 0.6), period = 4)))
setMethod("sim", c(object = "sarmacopula"), function(object, n = 1000) {
sim(sarma2arma(object), n = n)
})

#' Objective function for SARMA copula process
#'
#' @param theta vector of parameters of DARMA process
#' @param modelspec vector containing model order (p,q)(P,Q)D
#' @param u vector of data
#'
#' @return Value of objective function at parameters.
#' @keywords internal
#'
sarmacopula_objective <- function(theta, modelspec, u) {
  period <- modelspec[5]
  pars <- sarmavec2list(theta, modelspec)
  ar <- as.numeric(pars$ar) # Null to numeric(0) if necessary
  ma <- as.numeric(pars$ma)
  if (modelspec[3] > 0)
    ar <- expand_ar(ar, pars$sar, period)
  if (modelspec[4] > 0)
    ma <- expand_ma(ma, pars$sma, period)
  armacopula_objective(c(ar, ma), c(length(ar), length(ma)), u)
}

#' Residual function for sarmacopula object
#'
#' @param object a fitted sarmacopula object.
#' @param data the data to which copula is fitted.
#' @param trace extract trace instead of residuals.
#'
#' @return vector of model residuals
#' @keywords internal
#'
resid_sarmacopula <- function(object, data = NA, trace = FALSE){
  resid_armacopula(sarma2arma(object), data, trace)
}

#' @describeIn sarmacopula Calculate Kendall's tau values for sarmacopula model
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag.
#'
#' @export
#'
#' @examples
#' mod <- sarmacopula(list(ar = 0.5, ma = 0.4, sar = 0.2, sma = 0.6), period = 4)
#' kendall(mod)
setMethod("kendall", c(object = "sarmacopula"), function(object, lagmax = 20){
  kendall(sarma2arma(object), lagmax)
}
)

#' Generalized lagging for fitted sarmacopula objects
#'
#' @param copula a sarmacopula object.
#' @param data the data to which copula is fitted.
#' @param lagmax the maximum lag value.
#' @param glagplot logical value indicating generalized lag plot.
#'
#' @return If \code{glagplot} is \code{TRUE} a list of generalized lagged datasets
#' of maximum length 9 is returned to facilitate a generalized lagplot.
#' If \code{glagplot} is \code{FALSE} a vector of length \code{lagmax} containing
#' the Kendall rank correlations for the generalized lagged datasets is returned.
#'
#' @keywords internal
glag_for_sarmacopula <- function(copula, data, lagmax, glagplot = FALSE) {
  glag_for_armacopula(sarma2arma(copula), data, lagmax, glagplot)
}

#' @describeIn sarmacopula Prediction method for sarmacopula class
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#'
#' @export
#'
setMethod("predict", c(object = "sarmacopula"), function(object, data, x, type = "df") {
  predict(sarma2arma(object), data, x, type)
})

#' Transform a sarmacopula into a dvinecopula2 object
#'
#' @param object an object of class \linkS4class{sarmacopula}.
#'
#' @return An object of class \linkS4class{dvinecopula2}.
#' @export
#'
#' @examples
#' sarma2dvine(sarmacopula(list(ar = 0.5, ma = 0.4, sar = 0.2, sma = 0.6), period = 4))
sarma2dvine <- function(object){
  if (!(is(object, "sarmacopula")))
    stop("Not sarmacopula object")
  period <- object@modelspec["period"]
  dvinecopula2(family = "gauss", pars = object@pars,
                          kpacf = paste("kpacf_sarma", period, sep = ""))
}


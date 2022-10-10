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
  if ("ar" %in% names(pars)) {
    arpars <- pars$ar
    if (is.null(arpars))
      stop("No NULL values for parameters; omit from list instead")
    if (non_stat(arpars))
      stop("Non-stationary AR part")
    p <- length(arpars)
    names(pars$ar) <- paste("ar", 1:p, sep = "")
  }
  else {
    p <- 0
  }
  if ("ma" %in% names(pars)) {
    mapars <- pars$ma
    if (is.null(mapars))
      stop("No NULL values for parameters; omit from list instead")
    if (non_invert(mapars))
      stop("Non-invertible MA part")
    q <- length(mapars)
    names(pars$ma) <- paste("ma", 1:q, sep = "")
  }
  else {
    q <- 0
  }
  if ("sar" %in% names(pars)) {
    sarpars <- pars$sar
    if (is.null(sarpars))
      stop("No NULL values for parameters; omit from list instead")
    if (non_stat(sarpars))
      stop("Non-stationary SAR part")
    P <- length(sarpars)
    names(pars$sar) <- paste("sar", 1:P, sep = "")
  }
  else {
    P <- 0
  }
  if ("sma" %in% names(pars)) {
    smapars <- pars$sma
    if (is.null(smapars))
      stop("No NULL values for parameters; omit from list instead")
    if (non_invert(smapars))
      stop("Non-invertible SMA part")
    Q <- length(smapars)
    names(pars$sma) <- paste("sma", 1:Q, sep = "")
  }
  else {
    Q <- 0
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
  p <- object@modelspec[1]
  q <- object@modelspec[2]
  P <- object@modelspec[3]
  Q <- object@modelspec[4]
  if (p > 0) {
    arpars <- object@pars$ar
  } else {
    arpars <- NULL
  }
  if (q > 0) {
    mapars <- object@pars$ma
  } else {
    mapars <- NULL
  }
  if (P > 0) {
    sarpars <- object@pars$sar
  } else {
    sarpars <- NULL
  }
  if (Q > 0) {
    smapars <- object@pars$sma
  } else {
    smapars <- NULL
  }
  c(arpars, mapars, sarpars, smapars)
})

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
  D <- object@modelspec["period"]
  ar <- numeric()
  ma <- numeric()
  sar <- numeric()
  sma <- numeric()
  if ("ar" %in% names(object@pars))
    ar <- object@pars$ar
  if ("ma" %in% names(object@pars))
    ma <- object@pars$ma
  if ("sar" %in% names(object@pars)){
    sar <- object@pars$sar
    sar <- as.vector(sapply(sar, function(x, d){c(rep(0, d-1), x)}, d=D))
    ar <- -coefficients(polynom::polynomial(c(1, -sar)) * polynom::polynomial(c(1, -ar)))[-1]
  }
  if ("sma" %in% names(object@pars)){
    sma <- object@pars$sma
    sma <- as.vector(sapply(sma, function(x, d){c(rep(0, d-1), x)}, d=D))
    ma <- coefficients(polynom::polynomial(c(1, sma)) * polynom::polynomial(c(1, ma)))[-1]
  }
  output <- list()
  if (length(ar) > 0)
    output$ar <- ar
  if (length(ma >0))
    output$ma <- ma
  armacopula(pars = output)
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
  xdata <- qnorm(u)
  p <- modelspec[1]
  ar <- 0
  q <- modelspec[2]
  ma <- 0
  if (p > 0) {
    ar <- theta[1:p]
  }
  if (q > 0) {
    ma <- theta[(p + 1):(p + q)]
  }
  P <- modelspec[3]
  sar <- 0
  Q <- modelspec[4]
  sma <- 0
  if (P > 0) {
    sar <- theta[(p + q + 1):(p + q + P)]
  }
  if (Q > 0) {
    sma <- theta[(p + q + P + 1):(p + q + P + Q)]
  }
  D <- modelspec[5]
  if (non_stat(ar) | non_invert(ma) | non_stat(sar) | non_invert(sma)) {
    output <- NA
  } else {
    if (P > 0){
      sar <- as.vector(sapply(sar, function(x, d){c(rep(0, d-1), x)}, d=D))
      ar <- -coefficients(polynom::polynomial(c(1, -sar)) * polynom::polynomial(c(1, -ar)))[-1]
      if (length(ar) == 0)
        ar <- 0
    }
    if (Q > 0){
      sma <- as.vector(sapply(sma, function(x, d){c(rep(0, d-1), x)}, d=D))
      ma <- coefficients(polynom::polynomial(c(1, sma)) * polynom::polynomial(c(1, ma)))[-1]
      if (length(ma) == 0)
        ma <- 0
    }
    sp <- starmaStateSpace(ar, ma, c(p + D*P, q + D*Q))
    ans <- FKF::fkf(
      a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt, Zt = sp$Zt,
      HHt = sp$HHt, GGt = sp$GGt, yt = rbind(xdata)
    )
    output <- -ans$logLik + sum(log(dnorm(xdata)))
  }
  return(output)
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


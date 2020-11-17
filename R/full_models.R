#' Full models
#'
#' Class of objects for composite time series models consisting
#' of stationary copula processes and marginal distributions.
#'
#' @slot tscopula an object of class \linkS4class{tscopula}.
#' @slot margin an object of class \linkS4class{margin}.
#'
#' @return
#' @export
#'
setClass("tscm",
  slots = list(
    tscopula = "tscopula",
    margin = "margin"
  )
)

#' Show Method for tscm Class
#'
#' @param object an object of class \linkS4class{tscm}.
#'
#' @return
#' @export
#'
setMethod("show", "tscm", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("_______ \n")
  cat("MARGIN: \n")
  show(object@margin)
  cat("_______\n")
  cat("COPULA: \n")
  show(object@tscopula)
  if (is(object, "tscmfit")) {
    ests <- object@fit$par
    if ("skeleton" %in% names(attributes(ests))) {
      attributes(ests)$skeleton <- NULL
    }
    if (is.element("hessian", names(object@fit))) {
      ses <- safe_ses(object@fit$hessian)
      ests <- rbind(ests, ses)
      dimnames(ests)[[1]] <- c("par", "se")
    }
    cat("_________________________\n")
    cat("summary of all estimates:\n")
    print(ests)
    cat("convergence status:", object@fit$convergence, ", log-likelihood:",
      -object@fit$value, "\n",
      sep = " "
    )
  }
})

#' Coefficient Method for tscm Class
#'
#' @param object an object of class \linkS4class{tscm}.
#'
#' @return vector of coefficients of model.
#' @export
#'
setMethod("coef", "tscm", function(object) {
  c(coef(object@tscopula), coef(object@margin))
})

#' Constructor Function for time series
#'
#' @param tscopula an object of class \linkS4class{tscopula}.
#' @param margin an object of class \linkS4class{margin}.
#'
#' @return An object of class \linkS4class{tscm}.
#' @export
#'
#' @examples
#' tscm(dvinecopula(), margin("weibull2"))
tscm <- function(tscopula, margin = new("margin", name = "unif")) {
  if (is(tscopula, "tscopulafit")) {
    tscopula <- tscopula@tscopula
  }
  if (is(margin, "marginfit")) {
    margin <- margin@margin
  }
  new("tscm",
    tscopula = tscopula,
    margin = margin
  )
}

#' Simulation Method for tscm class
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param n length of realization.
#'
#' @return A realization of the time series of length n
#' @export
#'
#' @examples
#' mod <- tscm(dvinecopula(family = "gauss", pars = 0.5), margin("weibull2"))
#' sim(mod)
setMethod(
  "sim",
  c(x = "tscm"),
  function(x, n = 1000) {
    Utilde <- sim(x@tscopula, n)
    qmarg(x@margin, Utilde)
  }
)

#' Fitted tscm Model
#'
#' Class of objects for fitted \linkS4class{tscm} models.
#'
#' @slot tscopula an object of class \linkS4class{tscopula}.
#' @slot margin an object of class \linkS4class{margin}.
#' @slot data a vector or time series of data to which process has been fitted.
#' @slot fit a list containing details of the fit.
#'
#' @return
#' @export
#'
setClass(
  "tscmfit",
  contains = "tscm",
  slots = list(
    tscopula = "tscopula",
    margin = "margin",
    data = "ANY",
    fit = "list"
  )
)

#' Fit Method for tscm Class
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param y a vector or time series of data.
#' @param tsoptions a list of parameters passed to fitting.
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#' @param method character string specifying method.
#'
#'
#' @return an object of class \linkS4class{tscmfit}.
#' @export
#'
#' @examples
#' mod <- tscm(dvinecopula(family = "gauss", pars = 0.5), margin("weibull2"))
#' y <- sim(mod)
#' fit(mod, y)
setMethod(
  "fit", c(x = "tscm", y = "ANY"),
  function(x,
           y,
           tsoptions = list(),
           control = list(
             warn.1d.NelderMead = FALSE,
             trace = FALSE,
             maxit = 5000
           ),
           method = "IFM") {
    defaults <- list(hessian = FALSE, method = "Nelder-Mead", fulcrum = NA, avoidzero = TRUE)
    tsoptions <- setoptions(tsoptions, defaults)
    if (is(x, "tscmfit")) {
      tscopula <- x@tscopula
      margin <- x@margin
      if (is(tscopula, "tscopulafit")) {
        tscopula <- tscopula@tscopula
      }
      if (is(margin, "marginfit")) {
        margin <- margin@margin
      }
      x <- new("tscm", tscopula = tscopula, margin = margin)
    }
    if ((method == "full") & is(x@tscopula, "tscopulaU")) {
      method <- paste(method, "A", sep = "")
    }
    if ((method == "full") & is(x@tscopula, "vtscopula")) {
      method <- paste(method, "B", sep = "")
    }
    switch(method,
      empirical = fitEDF(x, y, tsoptions, control),
      IFM = fitSTEPS(x, y, tsoptions, control),
      fullA = fitFULLa(x, y, tsoptions, control),
      fullB = fitFULLb(x, y, tsoptions, control),
      stop("Not a known method")
    )
  }
)

#' Fit tscm Using Empirical Distribution Function
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param y a vector or time series of data.
#' @param tsoptions a list of parameters passed to fitting.
#' This differs according to the class of x.
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return an object of class \linkS4class{tscmfit}.
#' @keywords internal
#'
fitEDF <- function(x, y, tsoptions, control) {
  U <- strank(as.numeric(y))
  copfit <- fit(x@tscopula, U, tsoptions, control)
  new("tscmfit",
    tscopula = copfit@tscopula,
    margin = new("margin", name = "edf"),
    data = y,
    fit = copfit@fit
  )
}

#' Fit tscm in Two Steps
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param y a vector or time series of data.
#' @param tsoptions a list of parameters passed to fitting.
#' This differs according to the class of x.
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return an object of class \linkS4class{tscmfit}.
#' @keywords internal
#'
fitSTEPS <- function(x, y, tsoptions, control) {
  margfit <- fit(x@margin, y, tsoptions = tsoptions, control = control)
  U <- pmarg(margfit, y)
  copfit <- fit(x@tscopula, U, tsoptions = tsoptions, control = control)
  combinedfit <- list()
  names(margfit@fit$par) <- paste("margin.", names(margfit@fit$par), sep = "")
  combinedfit$par <- c(copfit@fit$par, margfit@fit$par)
  if ("hessian" %in% names(tsoptions)) {
    if (tsoptions$hessian) {
      combinedfit$hessian <- Matrix::bdiag(margfit@fit$hessian, copfit@fit$hessian)
    }
  }
  combinedfit$convergence <- sum(copfit@fit$convergence, margfit@fit$convergence)
  combinedfit$value <- sum(copfit@fit$value, margfit@fit$value)
  new("tscmfit",
    tscopula = copfit@tscopula,
    margin = margfit@margin,
    data = y,
    fit = combinedfit
  )
}

#' Fit tscm Jointly
#'
#' @param x an object of class \linkS4class{tscm}.
#' @param y a vector or time series of data.
#' @param tsoptions list of variables
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return an object of class \linkS4class{tscmfit}.
#' @keywords internal
#'
fitFULLa <- function(x, y, tsoptions, control) {
  dens <- eval(parse(text = paste("d", x@margin@name, sep = "")))
  cdf <- eval(parse(text = paste("p", x@margin@name, sep = "")))
  parlist <- x@tscopula@pars
  parlist$margin <- x@margin@pars[!x@margin@fixed]
  theta <- tsunlist(parlist, tsoptions$fulcrum)
  fit <- optim(
    par = theta,
    fn = tsc_objectivea,
    modelspec = x@tscopula@modelspec,
    modeltype = is(x@tscopula)[[1]],
    dens = dens,
    cdf = cdf,
    margfixed = x@margin@pars[x@margin@fixed],
    y = as.numeric(y),
    method = tsoptions$method,
    hessian = tsoptions$hessian,
    control = control
  )
  newpars <- tsrelist(fit$par, fulcrum = tsoptions$fulcrum)
  x@margin@pars[!x@margin@fixed] <- newpars$margin
  x@tscopula@pars <- newpars[names(newpars) != "margin"]
  new("tscmfit", tscopula = x@tscopula, margin = x@margin, data = y, fit = fit)
}



#' Objective Function for Full Fit
#'
#' @param theta
#' @param order
#' @param dens
#' @param cdf
#' @param margfixed
#' @param Vtrans
#' @param y
#' @param fulcrum
#' @param avoidzero
#'
#' @return
#' @keywords internal
#'
tsc_objectivea <-
  function(theta, modelspec, modeltype, dens, cdf, margfixed, y) {
    margpars <- theta[substring(names(theta), 1, 6) == "margin"]
    nonmargpars <- theta[substring(names(theta), 1, 6) != "margin"]
    names(margpars) <- substring(names(margpars), 8)
    margpars <- c(margpars, margfixed)
    dx <- do.call(dens, append(as.list(margpars), list(x = y, log = TRUE)))
    termA <- -sum(dx)
    if (is.na(termA)) {
      return(NA)
    }
    U <- do.call(cdf, append(margpars, list(q = y)))
    objective <- eval(parse(text = paste(modeltype, "_objective", sep = "")))
    termBC <-
      objective(nonmargpars, modelspec, U)
    return(termA + termBC)
  }

#' Convert tscopula Object to tscm Object
#'
#' @param from a \linkS4class{tscopula} object.
#' @param to a \linkS4class{tscm} object.
#'
#' @return a \linkS4class{tscm} object.
#' @export
#'
setMethod(
  "coerce", c(from = "tscopula", to = "tscm"),
  function(from, to = "tsc", strict = TRUE) {
    new("tscm", tscopula = from, margin = new("margin", name = "unif"))
  }
)

#' Convert tscopulafit Object to be tscmfit Object
#'
#' @param from a \linkS4class{tscopulafit} object.
#' @param to a \linkS4class{tscmfit} object.
#'
#' @return a \linkS4class{tscmfit} object.
#' @export
#'
setMethod(
  "coerce", c(from = "tscopulafit", to = "tscmfit"),
  function(from, to = "tscmfit", strict = TRUE) {
    new("tscmfit",
      tscopula = from,
      margin = new("margin", name = "unif"),
      data = from@data,
      fit = from@fit
    )
  }
)

#' logLik Method for tscmfit Class
#'
#' @param object an object of class \linkS4class{tscmfit}.
#'
#' @return an object of class logLik
#' @export
#'
setMethod("logLik", "tscmfit", function(object) {
  ll <- -object@fit$value[length(object@fit$value)]
  attr(ll, "nobs") <- length(object@data)
  attr(ll, "df") <- length(object@fit$par)
  class(ll) <- "logLik"
  ll
})




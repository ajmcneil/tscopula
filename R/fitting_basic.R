setClassUnion("tscopulaU", c("armacopula", "dvinecopula", "dvinecopula2"))

#' Calculate Standardized Ranks of Data
#'
#' @param x a vector or time series of data.
#'
#' @return A vector or time series of standardized ranks in the interval (0,1)
#' @export
#'
#' @examples
#' strank(rnorm(100))
strank <- function(x) {
  U <- rank(x, na.last = "keep", ties.method = "random") / (length(x) + 1)
  if (inherits(x, "zoo")) {
    attributes(U) <- attributes(x)
  }
  U
}

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

#' Fitted Time Series Copula Processes
#'
#' Class of objects for time series copula processes.
#'
#' @slot name name of time series copula process.
#'
#' @return
#' @export
#'
setClass("tscopulafit",
  contains = c("tscopula"),
  slots = list(
    tscopula = "tscopula",
    data = "ANY",
    fit = "list"
  )
)

#' Simulation Method for tscopulafit Class
#'
#' @param x an object of class \linkS4class{tscopulafit}.
#' @param n length of realization.
#'
#' @return A realization of a time series copula process.
#' @export
#'
#' @examples
#' ar1 <- armacopula(list(ar = 0.7))
#' data <- sim(ar1, 1000)
#' ar1fit <- fit(ar1, data)
#' sim(ar1fit)
setMethod("sim", c(x = "tscopulafit"), function(x, n = 1000) {
  sim(x@tscopula, n)
})

#' Coefficients for tscopulafit Class
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#'
#' @return Coefficients of fitted model.
#' @export
#'
setMethod("coef", c(object = "tscopulafit"), function(object) {
  coef(object@tscopula)
})

#' Show method for tscopulafit objects
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#'
#' @return parameters of tscopulafit model
#' @export
#'
setMethod("show", "tscopulafit", function(object) {
  if (is(object@tscopula, "vtscopula")) {
    show(object@tscopula)
  } else if (is(object@tscopula, "tscopulaU")) {
    cat("object class: ", is(object@tscopula)[[1]], "\n", sep = "")
    cat("name: ", object@tscopula@name, "\n", sep = "")
    if (is(object@tscopula, "dvinecopula")) {
      fams <- sapply(object@tscopula@modelspec, FUN = function(v) {
        tmp <- v$family
        if (v$rotation != 0) {
          tmp <- paste(tmp, v$rotation, sep = "")
        }
        if (is(object, "dvinecopulaNE")) {
          if (!(v$exchangeable)) {
            tmp <- paste(tmp, "_NE", sep = "")
          }
        }
        tmp
      })
      if (length(unique(fams)) == 1) {
        fams <- fams[1]
      }
      cat("copula family: ", fams, "\n")
    }
    if (is(object@tscopula, "dvinecopula2")) {
      modobject <- object@tscopula
      famname <- modobject@modelspec$family
      if (modobject@modelspec$rotation !=0)
        famname <- paste(famname, "with rotation", modobject@modelspec$rotation)
      cat("copula family: ", famname, "\n", sep = "")
      kpacf  <- modobject@modelspec$kpacf
      if (modobject@modelspec$maxlag != Inf)
        kpacf <- paste(kpacf, "with max lag", modobject@modelspec$maxlag)
      cat("KPACF: ", kpacf,"\n", sep = "")
    }
  }
  else {
    stop("Unknown tscopula type")
  }
  ests <- object@fit$par
  attributes(ests)$skeleton <- NULL
  if (is.element("hessian", names(object@fit))) {
    ses <- safe_ses(object@fit$hessian)
    ests <- rbind(ests, ses)
    dimnames(ests)[[1]] <- c("par", "se")
  }
  cat("_____________________\n")
  cat("Summary of estimates:\n")
  print(ests)
  cat("convergence status: ", object@fit$convergence, ", log-likelihood: ",
    -object@fit$value, "\n",
    sep = ""
  )
})

#' Calculate Standard Errors Safely
#'
#' @param hess a Hessian matrix from a model fit.
#'
#' @return a vector of standard errors.
#' @export
#'
safe_ses <- function(hess) {
  hessinverse <- tryCatch(solve(hess), error = function(e) {
    warning("Hessian can't be inverted")
    return(diag(rep(NA, length(diag(hess)))))
  })
  sqrt(abs(diag(hessinverse)))
}

#' Fit Method for tscopulafit Class
#'
#' @param x an object of class \linkS4class{tscopulafit}.
#' @param y vector or time series of data to which the copula process is to be fitted.
#' @param tsoptions list of options
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{tscopulafit}.
#' @export
#'
#' @examples
#' ar1 <- armacopula(list(ar = 0.7))
#' data <- sim(ar1, 1000)
#' ar1fit <- fit(fit(ar1, data), sim(ar1, 1000))
setMethod(
  "fit", c(x = "tscopulafit", y = "ANY"),
  function(x, y,
           tsoptions = list(),
           control = list(warn.1d.NelderMead = FALSE)) {
    x <- x@tscopula
    fit(x, y, tsoptions = tsoptions, control = control)
  }
)

#' Fit Method for tscopulaU Class
#'
#' @param x an object of class \linkS4class{tscopulaU}.
#' @param y vector or time series of data to which the copula process is to be fitted.
#' @param tsoptions list of options
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{tscopulafit}.
#' @export
#'
#' @examples
#' data <- sim(armacopula(list(ar = 0.5, ma = 0.4)), n = 1000)
#' fit(armacopula(list(ar = 0.5, ma = 0.4)), data)
setMethod(
  "fit", c(x = "tscopulaU", y = "ANY"),
  function(x, y,
           tsoptions = list(),
           control = list(warn.1d.NelderMead = FALSE)) {
    defaults <- list(hessian = FALSE, method = "Nelder-Mead")
    tsoptions <- setoptions(tsoptions, defaults)
    objective <- eval(parse(text = paste(is(x)[[1]], "_objective", sep = "")))
    fit <- optim(
      par = tsunlist(x@pars),
      fn = objective,
      modelspec = x@modelspec,
      u = as.numeric(y),
      method = tsoptions$method,
      hessian = tsoptions$hessian,
      control = control
    )
    x@pars <- tsrelist(fit$par)
    new("tscopulafit",
      tscopula = x,
      data = y,
      fit = fit
    )
  }
)

#' Unlist Function for tscopula Objects
#'
#' @param pars list of parameters
#' @param fulcrum numerical value for fixed fulcrum or NA
#'
#' @return relistable vector of parameters
#' @keywords internal
#'
tsunlist <- function(pars, fulcrum = NA) {
  if ("vt" %in% names(pars)) {
    if (is.null(pars$vt)) {
      pars <- pars[-length(pars)]
    }
  }
  if (!is.na(fulcrum)) {
    if (length(pars$vt) == 1) pars <- pars[!(names(pars) == "vt")]
    if (length(pars$vt) > 1) pars$vt <- pars$vt[-1]
  }
  unlist(as.relistable(pars))
}

#' Relist Function for tscopula Objects
#'
#' @param ests a relistable vector of parameter estimates
#' @param fulcrum numerical value for fixed fulcrum or NA
#'
#' @return list of parameters for \linkS4class{tscopulafit} object.
#' @keywords internal
#'
tsrelist <- function(ests, fulcrum = NA) {
  newpars <- relist(ests)
  attributes(newpars)$class <- NULL
  if (!is.na(fulcrum)) {
    newpars$vt <- c(delta = fulcrum, newpars$vt)
  }
  newpars
}

#' Set Optional Choices for tscopula Fitting
#'
#' @param tsoptions a list of options chosen by user.
#' @param defaults a list of defaults specified in function.
#'
#' @return the definitive list of options to be used in function.
#' @keywords internal
#'
setoptions <- function(tsoptions, defaults) {
  defnames <- names(defaults)
  defaults[(optnames <- names(tsoptions))] <- tsoptions
  defaults
}

#' logLik Method for tscopulafit Class
#'
#' @param object an object of class \linkS4class{tscopulafit}.
#'
#' @return an object of class logLik
#' @export
#'
setMethod("logLik", "tscopulafit", function(object) {
  logLik(as(object, "tscmfit"))
})

#' Plot Method for tscopulafit Class
#'
#' @param x an object of class \linkS4class{tscopulafit}.
#' @param y missing.
#' @param plotoption number of plot required.
#' @param plottype type of plot required.
#' @param bw logical variable specifying whether black-white options should be chosen.
#' @param klimit maximum lag value for dvinecopula2 cplots
#' @param vtransform logical variable specifying whether v-transform is plotted,
#'
#' @return
#' @export
#'
#' @examples
#' data <- sim(armacopula(list(ar = 0.5, ma = 0.4)), n = 1000)
#' fit <- fit(armacopula(list(ar = 0.5, ma = 0.4)), data)
#' plot(fit, plotoption = 1) # try plotoption 1 through 5
setMethod("plot", c(x = "tscopulafit", y = "missing"),
          function(x, plotoption = 1, plottype = "copula", bw = FALSE, klimit = 30) {
            if (plottype == "vtransform")
              plot(x@tscopula@Vtransform)
            else
              switch(is(x@tscopula)[1],
                     armacopula = plot_armacopula(x@tscopula, x@data, plotoption, bw),
                     dvinecopula = plot_dvinecopula(x@tscopula, x@data, plotoption, bw),
                     dvinecopula2 = plot_dvinecopula2(x@tscopula, x@data, plotoption, bw, klimit),
                     vtscopula = {
                       Vdata <- vtrans(x@tscopula@Vtransform, x@data, correction = TRUE)
                       plot(new("tscopulafit",
                                tscopula = x@tscopula@Vcopula,
                                data = Vdata,
                                fit = x@fit
                       ),
                       plotoption = plotoption,
                       bw = bw
                       )
                     })
})


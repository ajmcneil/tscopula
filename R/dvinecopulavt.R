#' D-vine copula processes with v-transforms
#'
#' Class of objects for d-vine copula processes. See \link{dvinecopulavt} for more details.
#'
#' @slot name name of the d-vine copula process.
#' @slot modelspec list containing the family, rotation, and name of KPACF
#' @slot pars list comprising of the parameters.
#'
#' @export
#'
setClass("dvinecopulavt", contains = "tscopula", slots = list(
  name = "character",
  modelspec = "list",
  pars = "list"
))

#' Constructor function for dvinecopulavt process
#'
#' This function sets up a stationary d-vine process of finite or infinite order based on a single
#' inverse-v-transformed copula family from a subset of those that can be implemented using
#' \code{\link[rvinecopulib]{bicop_dist}} in the \code{rvinecopulib} package.
#'
#' The permitted choices of base copula family are currently Joe, Gumbel, Frank or Clayton survival. If
#' Clayton is chosen, the \code{rotation} argument must be set to 180, while if Joe or Gumbel are chosen, the
#' \code{rotation} argument must be zero (which is the default); any other options will return an error
#'
#' The copulas are parameterized using the Kendall partial autocorrelation function (kpacf)
#' of the base copula sequence specified
#' by the \code{kpacf} argument. The default choice is the kpacf of a standard ARMA process which is
#' implemented in the function \code{\link{kpacf_arma}}. The parameters
#' of the kpacf should be set as a list using the \code{pars} argument; the required parameters should usually
#' be clear from the documentation of the chosen kpacf function and must be correctly named.
#'
#' For data showing stochastic volatility, we expect positive serial dependencies in the base copula sequence.
#' For this reason, we do not consider models where the kpacf takes negative values.
#'
#' The \code{maxlag} parameter specifies the maximum lag of the process; a finite number gives a finite-order
#' stationary d-vine process, but the parameter may also be set to \code{Inf} for an infinite-order process.
#'
#' @param family family name
#' @param rotation a scalar specifying the rotation (default is 0)
#' @param kpacf a character string giving the name of the Kendall pacf
#' @param pars a list containing the parameters of the model
#' @param maxlag a scalar specifying the maximum lag
#'
#' @return An object of class \linkS4class{dvinecopulavt}.
#' @export
#'
#' @examples
#' dvinecopulavt(family = "joe", kpacf = "kpacf_arma",
#' pars = list(ar = 0.95, ma = -0.85), maxlag = 30)
dvinecopulavt <- function(family = "joe",
                         rotation = 0,
                         kpacf = "kpacf_arma",
                         pars = list(ar = 0.1, ma = 0),
                         maxlag = Inf) {
  if (!(is(family, "character")))
    stop("copula family must be specified by name")
  fam <- tolower(family)
  permitted1 <- ((fam %in% c("joe", "gumbel", "frank")) & (rotation == 0))
  permitted2 <- ((fam %in% c("clayton","frank")) & (rotation == 180))
  if (!(permitted1 | permitted2))
    stop ("Illegal family and rotation choice: can have Joe or
          Gumbel with no rotation or Clayton with 180 degree rotation")
  if (is.null(names(pars)))
    stop("parameters should be named (p1 and p2 for exp/power)")
  if (!("delta1" %in% names(pars)))
    pars <- c(pars, delta1 = 0.5)
  if (!("delta2" %in% names(pars)))
    pars <- c(pars, delta2 = 0.5)
  modelspec <- list(family = fam,
                    rotation = rotation,
                    kpacf = kpacf,
                    maxlag = maxlag,
                    npar = length(unlist(pars)),
                    negtau = "none")
  new("dvinecopulavt",
      name = paste("d-vine-vt"),
      modelspec = modelspec,
      pars = pars
  )
}

#' @describeIn dvinecopulavt Coef Method for dvinecopulavt class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", c(object = "dvinecopulavt"), function(object) {
  return(unlist(object@pars))
})

#' @describeIn dvinecopulavt Show method for dvinecopulavt class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", c(object = "dvinecopulavt"), function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  famname <- object@modelspec$family
  if (object@modelspec$rotation !=0)
    famname <- paste(famname, "with rotation", object@modelspec$rotation)
  cat("copula family: ", famname, "\n", sep = "")
  cat("linear v-transforms with fulcrums at", object@pars$delta1, "and", object@pars$delta2, "\n")
  kpacf  <- object@modelspec$kpacf
  if (object@modelspec$maxlag != Inf)
    kpacf <- paste(kpacf, "with max lag", object@modelspec$maxlag)
  cat("KPACF: ", kpacf,"\n", sep = "")
  cat("parameters: \n")
  print(coef(object))
})

#' Objective function for dvinecopulavt process
#'
#' @param theta parameters of kpacf
#' @param modelspec list specifying model
#' @param u data
#'
#' @return Value of objective function at parameters.
#'
#' @keywords internal
#'
dvinecopulavt_objective <- function(theta, modelspec, u) {
  n <- length(u)
  kpacf <- eval(parse(text = modelspec$kpacf))
  tauvals <- kpacf((n-1), theta)
  if (is.na(sum(tauvals)))
    return(NA)
  k <- 1
  largetau <- (abs(tauvals) > .Machine$double.eps)
  if (sum(largetau) > 0)
    k <- max((1:(n-1))[largetau])
  k <- min(k, modelspec$maxlag)
  if (min(tauvals[1:k]) <0)
    return(NA)
  pc_list <- vector("list", k)
  for (i in 1:k) {
    fam <- tolower(modelspec$family)
    rot <- modelspec$rotation
    coppars <- ktau_to_par(
      family = fam,
      tau = tauvals[i]
    )
    if (is.na(coppars))
      return(NA)
    pc_list[[i]] <- tryCatch(rvinecopulib::bicop_dist(
      family = fam,
      rotation = rot,
      parameters = coppars
    ),
    error = function(e) {
      return(NA)
    }
    )
    if (is.na(pc_list[[i]][[1]])) {
      return(NA)
    }
  }
  v <- cbind(u[1:(n - 1)], u[2:n])
  delta <- as.numeric(c(theta["delta1"], theta["delta2"]))
  LL <- 0
  for (j in 1:k) {
    v_vt <- cbind(vtrans(Vlinear(delta[1]), v[,1]), vtrans(Vlinear(delta[2]), v[,2]))
    LL <- LL + sum(log(rvinecopulib::dbicop(u = v_vt, family = pc_list[[j]])))
    if (j == k) {
      return(-LL)
    }
    n <- dim(v)[1]
    v <- cbind(
      hbicop(v[(1:(n - 1)), ], cond_var = 2, family = pc_list[[j]], delta),
      hbicop(v[(2:n), ], cond_var = 1, family = pc_list[[j]], delta)
    )
  }
}

#' h-function for linear inverse-v-transformed copula
#'
#' @param u two-column matrix of data at which h-function is evaluated
#' @param cond_var identity of conditioning variable (1/2)
#' @param family name of copula family
#' @param delta vector of two fulcrum values
#'
#' @return vector of values of h-function
#'
#' @keywords internal
#'
hbicop <- function(u, cond_var, family, delta) {
  cond1 <- u[,2] > delta[2]
  cond2 <- u[,1] > delta[2]
  mult1 <- -(delta[2]^(1-cond1))*(delta[2] - 1)^cond1
  mult2 <- -(delta[1]^(1-cond2))*(delta[1] - 1)^cond2
  v <- cbind(vtrans(Vlinear(delta[1]), u[,1]), vtrans(Vlinear(delta[2]), u[,2]))
  switch(cond_var,
           mult1 * rvinecopulib::hbicop(u = v, cond_var = 1, family = family) + delta[2],
           mult2 * rvinecopulib::hbicop(u = v, cond_var = 2, family = family) + delta[1]
    )
}

#' @describeIn dvinecopulavt Calculate Kendall's tau values for core pair copulas
#' in d-vine copula model with v-transforms
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag.
#'
#' @export
#'
#' @examples
#' copmod <- dvinecopulavt(family = "joe", kpacf = "kpacf_arma",
#' pars = list(ar = 0.95, ma = -0.85), maxlag = 30)
#' kendall(copmod)
setMethod("kendall", c(object = "dvinecopulavt"), function(object, lagmax = 20) {
  kpacf <- eval(parse(text = object@modelspec$kpacf))
  kpacf(lagmax, object@pars)
}
)

#' Residual function for dvinecopulavt object
#'
#' @param object a fitted dvinecopulavt object.
#' @param data the data to which copula is fitted.
#' @param trace extract trace instead of residuals.
#'
#' @return vector of model residuals
#' @keywords internal
#'
resid_dvinecopulavt <- function(object, data = NA, trace = FALSE){
  browser()
  n <- length(data)
  pc_list <- mklist_dvine2(object, n-1, truncate = TRUE, tol = 1/3)
  k <- length(pc_list)
  if (trace)
    target <- rep(0.5, n)
  else
    target <- data
  res <- rep(NA, n)
  res[1] <- target[1]
  if (k >1){
    for (i in 2:k){
      pcs <- lapply(1:(i-1), function(j) {
        replicate(i - j, pc_list[[j]], simplify = FALSE)
      })
      vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure(i:1))
      vals <- c(data[1:(i-1)], target[i])
      res[i] <- rvinecopulib::rosenblatt(t(vals), vc_short)[i]
    }
  }
  pcs <- lapply(1:k, function(j) {
    replicate(k + 1 - j, pc_list[[j]], simplify = FALSE)
  })
  vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure((k + 1):1))
  for (i in ((k+1):n)){
    vals <- c(data[(i-k):(i-1)], target[i])
    res[i] <- rvinecopulib::rosenblatt(t(vals), vc_short)[k+1]
  }
  qnorm(res)
}

#' Rosenblatt Transforms for copulas with v-transforms
#'
#' @param latest argument of Rosenblatt functiom
#' @param previous conditioning vales
#' @param cv identity of conditioning variable
#' @param pcs list of pair copulas
#' @param delta vector fulcrum parameters
#'
#' @return
#' @export
#'
RT <- function(latest, previous, cv, pcs, delta) {
  if ((latest == 0) | (latest == 1)) {
    return(latest)
  }
  k <- length(previous)
  x <- c(previous, latest)
  if (k == 1) {
    return(hbicop(cbind(x[1], x[2]), cv, pcs[[1]], delta))
  } else {
    return(hbicop(cbind(
      RT(x[k], x[1:(k - 1)], cv = 2, pcs, delta),
      RT(x[k + 1], x[2:k], cv = 1, pcs, delta)
    ), cv, pcs[[k]], delta))
  }
}

#' Rosenblatt forward function with v-transforms
#'
#' @param x argument of Rosenblatt functiom
#' @param u conditioning values
#' @param pcs list of pair copulas
#' @param delta vector fulcrum parameters
#'
#' @return
#' @export
#'
Rforward <- function(x, u, pcs, delta) {
  k <- length(u)
  if (k == 1) {
    return(as.numeric(hbicop(cbind(u, x), 1, pcs[[1]], delta)))
  } else {
    return(hbicop(cbind(
      Rbackward(u[k], u[(k-1):1], pcs, delta),
      Rforward(x, u[1:(k-1)], pcs, delta)
    ), 1, pcs[[k]], delta))
  }
}

#' Rosenblatt backward function with v-transforms
#'
#' @param x argument of Rosenblatt functiom
#' @param u conditioning vales
#' @param pcs list of pair copulas
#' @param delta vector fulcrum parameters
#'
#' @return
#' @export
#'
Rbackward <- function(x, u, pcs, delta) {
  k <- length(u)
  if (k == 1) {
    return(as.numeric(hbicop(cbind(x, u), 2, pcs[[1]], delta)))
  } else {
    return(hbicop(cbind(
      Rbackward(x, u[1:(k-1)], pcs, delta),
      Rforward(u[k], u[(k-1):1], pcs, delta)
    ), 2, pcs[[k]], delta))
  }
}

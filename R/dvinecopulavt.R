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
#' The permitted choices of base copula family are currently Joe, Gumbel, Frank, ast or Clayton survival. If
#' Clayton is chosen, the \code{rotation} argument must be set to 180, while if Joe, Gumbel or ast are chosen, the
#' \code{rotation} argument must be zero (which is the default); any other options will return an error
#'
#' The copulas are parameterized using the Kendall partial autocorrelation function (kpacf)
#' of the base copula sequence specified
#' by the \code{kpacf} argument. The default choice is the kpacf of a standard ARMA process which is
#' implemented in the function \code{\link{kpacf_arma}}. The parameters
#' of the kpacf should be set as a list using the \code{pars} argument; the required parameters should usually
#' be clear from the documentation of the chosen kpacf function and must be correctly named.
#'
#' The arguments \code{vt1} and \code{vt2} are used to enter two parametric v-transforms which may be created, for example,
#' by \code{\link{Vlinear}} or \code{\link{V2p}}. However, the latter is very slow and the
#' variable \code{V2override} has to be set to \code{TRUE} if you want to include 2-parameter
#' v-transforms. While fitting is possible, residual analysis and simulation are almost always
#' prohibitively slow.
#'
#' For data showing stochastic volatility, we expect positive serial dependencies in the base copula sequence.
#' For this reason, we do not consider models where the kpacf takes negative values.
#'
#' In practice, the sequence of base copulas will be truncated at the last copula for which the kpacf exceeds \code{tautol}.
#' The \code{maxlag} parameter is typically used to force the truncation to take place at a lower lag (to increase speed).
#' This can also be achieved by increasing the value of \code{tautol}.
#'
#' @param family family name
#' @param rotation a scalar specifying the rotation (default is 0)
#' @param kpacf a character string giving the name of the Kendall pacf
#' @param pars a list containing the parameters of the model
#' @param vt1 first v-transform
#' @param vt2 second v-transform
#' @param tautol scalar value at which kpacf is truncated
#' @param maxlag a scalar which can be used to force a given value for maximum lag
#' @param V2override logical variable stating whether 2-parameter v-transform
#' should be permitted
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
                         vt1 = Vlinear(0.5),
                         vt2 = Vlinear(0.5),
                         tautol = 1e-04,
                         maxlag = Inf,
                         V2override = FALSE) {
  if (!(is(family, "character")))
    stop("copula family must be specified by name")
  fam <- tolower(family)
  permitted1 <- ((fam %in% c("joe", "gumbel", "frank", "ast")) & (rotation == 0))
  permitted2 <- ((fam %in% c("clayton","frank")) & (rotation == 180))
  if (!(permitted1 | permitted2))
    stop ("Illegal family and rotation choice: can have Joe, ast and
          Gumbel with no rotation or Clayton with 180 degree rotation")
  if (((vt1@name == "V2p") | (vt2@name == "V2p")) & (!V2override))
    stop("Use of 2-parameter v-transforms very slow. Set V2override to TRUE if you
         really want to use them.")
  vtpermitted <- (vt1@name %in% c("Vlinear", "V2p")) & (vt2@name %in% c("Vlinear", "V2p"))
  if (!(vtpermitted))
    stop("Only Vlinear and V2p v-transforms are currently permitted")
  vt1pars <- coef(vt1)
  names(vt1pars) <- paste(names(vt1pars),"1",sep="")
  if (!("delta1" %in% names(pars)))
    pars <- c(pars, vt1pars["delta1"])
  if ((vt1@name == "V2p") & (!("kappa1" %in% names(pars))))
    pars <- c(pars, vt1pars["kappa1"])
  vt2pars <- coef(vt2)
  names(vt2pars) <- paste(names(vt2pars),"2",sep="")
  if (!("delta2" %in% names(pars)))
    pars <- c(pars, vt2pars["delta2"])
  if ((vt2@name == "V2p") & (!("kappa2" %in% names(pars))))
    pars <- c(pars, vt2pars["kappa2"])
  modelspec <- list(family = fam,
                    rotation = rotation,
                    kpacf = kpacf,
                    tautol = tautol,
                    maxlag = maxlag,
                    npar = length(unlist(pars)),
                    vt1 = vt1,
                    vt2 = vt2,
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
  cat("v-transform 1: ", object@modelspec$vt1@name, "\n")
  cat("v-transform 2: ", object@modelspec$vt2@name, "\n")
  cat("KPACF: ", object@modelspec$kpacf,"\n", sep = "")
  cat("parameters: \n")
  print(coef(object))
})

#' Helper function for building lists of objects to describe copula sequence
#'
#' @param family name of copula family
#' @param rotation rotation for copula family
#' @param parameters parameters of copula
#'
#' @return a list containing specification of copula family
#'
#' @keywords internal
#'
tscopula_bicop_dist <- function(family, rotation, parameters){
  if (family == "ast")
    return(list(family = family, rotation = rotation, parameters = parameters))
  else
    return(tryCatch(rvinecopulib::bicop_dist(
      family = family,
      rotation = rotation,
      parameters = parameters
    ),
    error = function(e) {
      return(NA)
    }
    ))
}

#' Compute density of bivariate copula
#'
#' @param u matrix of data with two columns
#' @param family list describing copula family
#'
#' @return a vector containing values of density
#'
#' @keywords internal
#'
tscopula_dbicop <- function(u, family){
  if (family$family == "ast")
    return(copula::dCopula((1+u)/2, copula::tCopula(0, df = family$parameters)))
  else
    return(rvinecopulib::dbicop(u = u, family = family))
}

#' Compute h function of bivariate copula
#'
#' @param u matrix of data with two columns
#' @param cond_var conditioning variable (1 or 2)
#' @param family list describing copula family
#' @param inverse logical variable specifying whether inverse should be returned
#'
#' @return a vector containing values of density
#'
#' @keywords internal
#'
tscopula_hbicop <- function(u, cond_var, family, inverse = FALSE){
  if (family$family == "ast"){
    theta <- family$parameters
    mult <- switch(cond_var,
                    qt((1+u[,1])/2, df = theta),
                    qt((1+u[,2])/2, df = theta))
    ratio <- sqrt((theta + 1)/(theta + mult^2))
    if (inverse){
      mult2 <- switch(cond_var,
                      qt((1+u[,2])/2, df = theta + 1),
                      qt((1+u[,1])/2, df = theta + 1))
      return(2*pt(mult2 / ratio , df = theta) -1)
    }
    else{
      mult2 <- switch(cond_var,
                      qt((1+u[,2])/2, df = theta),
                      qt((1+u[,1])/2, df = theta))
      return(2*pt(ratio * mult2, df = theta + 1) -1)
    }
  }
  else
    return(rvinecopulib::hbicop(u = u, cond_var = cond_var, family = family, inverse = inverse))
}


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
  k <- effective_maxlag(tauvals, modelspec$tautol, modelspec$maxlag)
  if (min(tauvals[1:k]) < 0)
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
    pc_list[[i]] <- tscopula_bicop_dist(fam, rot, coppars)
    if (is.na(pc_list[[i]][[1]])) {
      return(NA)
    }
  }
  v <- cbind(u[1:(n - 1)], u[2:n])
  vt1 <- getvt(modelspec$vt1, theta, 1)
  vt2 <- getvt(modelspec$vt2, theta, 2)
  if ((vt1@pars["delta"] <= 0) | (vt1@pars["delta"] >= 1) |
      (vt2@pars["delta"] <= 0) | (vt2@pars["delta"] >= 1))
    return(NA)
  if (length(vt1@pars) > 1)
    if (vt1@pars["kappa"] <= 0)
      return(NA)
  if (length(vt2@pars) > 1)
    if (vt2@pars["kappa"] <= 0)
      return(NA)
  LL <- 0
  for (j in 1:k) {
    v_vt <- cbind(vtrans(vt1, v[,1]), vtrans(vt2, v[,2]))
    LL <- LL + sum(log(tscopula_dbicop(u = v_vt, family = pc_list[[j]])))
    if (j == k) {
      return(-LL)
    }
    n <- dim(v)[1]
    v <- cbind(
      hbicop_vt(v[(1:(n - 1)), ], cond_var = 2, family = pc_list[[j]],
             vt1 = vt1, vt2 = vt2),
      hbicop_vt(v[(2:n), ], cond_var = 1, family = pc_list[[j]],
             vt1 = vt1, vt2 = vt2)
    )
  }
}

#' Helper function for reparameterizing v-transforms
#'
#' @param vt a v-transform
#' @param theta vector of parameters containing parameters of v-transform
#' @param number number of v-transform (1 or 2)
#'
#' @return a v-transform
#'
#' @keywords internal
#'
getvt <- function(vt, theta, number){
  vtnms <- names(vt@pars)
  vtpars <- theta[paste(vtnms, number, sep = "")]
  names(vtpars) <- vtnms
  vt@pars <- vtpars
  vt
}

#' h-function for linear inverse-v-transformed copula
#'
#' @param u two-column matrix of data at which h-function is evaluated
#' @param cond_var identity of conditioning variable (1/2)
#' @param family name of copula family
#' @param vt1 first v-transform
#' @param vt2 second v-transform
#' @param inverse logical variable specifying whether inverse is taken
#'
#' @return vector of values of h-function
#'
#' @keywords internal
#'
hbicop_vt <- function(u, cond_var, family, vt1, vt2, inverse = FALSE) {
  if ((cond_var == 1) & (vt2@name == "Vlinear"))
  {
    delta2 <- vt2@pars["delta"]
    cond1 <- u[,2] > delta2
    mult1 <- (-1)*(delta2^(1-cond1))*((delta2 - 1)^cond1)
    v <- cbind(vtrans(vt1, u[,1]), vtrans(vt2, u[,2]))
    hval <- tscopula_hbicop(u = v, cond_var = 1,
                              family = family, inverse = inverse)
    return(mult1 * hval + delta2)
  }
  else if ((cond_var == 2) & (vt1@name == "Vlinear"))
  {
    delta1 <- vt1@pars["delta"]
    cond2 <- u[,1] > delta1
    mult2 <- (-1)*(delta1^(1-cond2))*((delta1 - 1)^cond2)
    v <- cbind(vtrans(vt1, u[,1]), vtrans(vt2, u[,2]))
    hval <- tscopula_hbicop(u = v, cond_var = 2,
                            family = family, inverse = inverse)
    return(mult2 * hval + delta1)
  }
  else if (cond_var == 1)
  {
    u[,1] <- pmax(u[,1], 1e-06) # lower bound on conditioning variable
    integrand <- function(x, y, family, vt1, vt2){
      tscopula_dbicop(u = cbind(vtrans(vt1, y), vtrans(vt2, x)), family = family)
    }
    tmp <- sapply(seq_along(u[,1]), function(i) integrate(integrand, lower = 0, upper = u[i,2], y = u[i,1],
                                                          family = family, vt1 = vt1, vt2 = vt2, stop.on.error=FALSE)$value)
    return(pmax(pmin(tmp,1),0))
  }
  else if (cond_var == 2)
  {
    u[,2] <- pmax(u[,2], 1e-06) # lower bound on conditioning variable
    integrand <- function(x, y, family, vt1, vt2){
      tscopula_dbicop(u = cbind(vtrans(vt1, x), vtrans(vt2, y)), family = family)
    }
    tmp <- sapply(seq_along(u[,2]), function(i) integrate(integrand, lower = 0, upper = u[i,1], y = u[i,2],
                                                          family = family, vt1 = vt1, vt2 = vt2, stop.on.error=FALSE)$value)
    return(pmax(pmin(tmp,1),0))
  }
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
  tau <- kpacf(lagmax, object@pars)
  k <- effective_maxlag(tau, object@modelspec$tautol, object@modelspec$maxlag)
  tau[which((1:lagmax) > k)] <- 0
  tau
}
)

#' Make list of pair copulas for dvinecopulavt object
#'
#' @param x an object of class dvinecopulavt
#' @param maxlag maximum possible lag to consider
#'
#' @return a list of pair copulas
#' @keywords internal
#'
mklist_dvinevt <- function(x, maxlag){
  kpacf <- eval(parse(text = x@modelspec$kpacf))
  tauvals <- kpacf(maxlag, x@pars)
  k <- effective_maxlag(tauvals, x@modelspec$tautol, x@modelspec$maxlag)
  pc_list <- vector("list", k)
  for (i in 1:k) {
    fam <- tolower(x@modelspec$family)
    rot <- x@modelspec$rotation
    coppars <- ktau_to_par(
      family = fam,
      tau = tauvals[i]
    )
    if (fam == "t")
      coppars <- c(coppars, x@pars$df)
    pc_list[[i]] <- tscopula_bicop_dist(
      family = fam,
      rotation = rot,
      parameters = coppars)
  }
  pc_list
}

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
  n <- length(data)
  pc_list <- mklist_dvinevt(object, n-1)
  k <- length(pc_list)
  if (trace)
    target <- rep(0.5, n)
  else
    target <- data
  res <- rep(NA, n)
  vt1 <- object@modelspec$vt1
  vt2 <- object@modelspec$vt2
  res[1] <- target[1]
  for (i in 2:k)
    res[i] <- Rforward(data[i], matrix(data[(i-1):1], ncol = i-1, nrow =1), pc_list, vt1, vt2)
  pastdata <- matrix(NA, ncol = k, nrow = n - k)
  for (i in ((k+1):n))
    pastdata[i-k,] <- data[(i-1):(i-k)]
  res[(k+1):n] <- Rforward(data[(k+1):n], pastdata, pc_list, vt1, vt2)
  qnorm(res)
}

#' Rosenblatt forward function with v-transforms
#'
#' @param x vector argument of Rosenblatt functiom
#' @param u matrix of conditioning values. Number of rows
#' must be either 1 or same length as x. Number of columns
#' should not be much more than 15 (due to repeated recursive
#' calling)
#' @param pcs list of pair copulas
#' @param vt1 first v-transform
#' @param vt2 second v-transform
#'
#' @return vector of same length as x
#'
#' @keywords internal
#'
Rforward <- function(x, u, pcs, vt1, vt2) {
  if (!(is.matrix(u)))
    stop("Must have matrix of conditioning variables")
  k <- dim(u)[2]
  if (k == 1) {
    return(as.numeric(hbicop_vt(cbind(as.vector(u), x), cond_var = 1,
                             family = pcs[[1]], vt1 = vt1, vt2 = vt2)))
  } else {
    return(hbicop_vt(cbind(
      Rbackward(u[, k], matrix(u[, ((k-1):1)], ncol = k-1), pcs, vt1, vt2),
      Rforward(x, matrix(u[, (1:(k-1))], ncol = k-1), pcs, vt1, vt2)
    ), cond_var = 1, family = pcs[[k]], vt1 = vt1, vt2 = vt2))
  }
}

#' Rosenblatt backward function with v-transforms
#'
#' @param x vector argument of Rosenblatt functiom
#' @param u matrix of conditioning values. Number of rows
#' must be either 1 or same length as x. Number of columns
#' should not be much more than 15 (due to repeated recursive
#' calling)
#' @param pcs list of pair copulas
#' @param vt1 first v-transform
#' @param vt2 second v-transform
#'
#' @return vector of same length as x
#'
Rbackward <- function(x, u, pcs, vt1, vt2) {
  if (!(is.matrix(u)))
    stop("Must have matrix of conditioning variables")
  k <- dim(u)[2]
  if (k == 1) {
    return(as.numeric(hbicop_vt(cbind(x, as.vector(u)), cond_var = 2,
                             pcs[[1]], vt1 = vt1, vt2 = vt2)))
  } else {
    return(hbicop_vt(cbind(
      Rbackward(x, matrix(u[, (1:(k-1))], ncol = k-1), pcs, vt1, vt2),
      Rforward(u[, k], matrix(u[, ((k-1):1)], ncol = k-1), pcs, vt1, vt2)
    ), cond_var = 2, family = pcs[[k]], vt1 = vt1, vt2 = vt2))
  }
}

#' Inverse Rosenblatt forward function with v-transforms
#'
#' @param x vector argument of Rosenblatt functiom
#' @param u matrix of conditioning values. Number of rows
#' must be either 1 or same length as x. Number of columns
#' should not be much more than 15 (due to repeated recursive
#' calling)
#' @param pcs list of pair copulas
#' @param vt1 first v-transform
#' @param vt2 second v-transform
#'
#' @return vector of same length as x
#'
RforwardI <- function(x, u, pcs, vt1, vt2) {
  if (!(is.matrix(u)))
    stop("Must have matrix of conditioning variables")
  k <- dim(u)[2]
  if (k == 1) {
    hval <- hbicop_vt(cbind(as.vector(u), x),
                      cond_var = 1, family = pcs[[1]], vt1 = vt1,
                      vt2 = vt2, inverse = TRUE)
    return(as.numeric(hval))
  } else {
    arg1 <- hbicop_vt(cbind(Rbackward(u[,k], matrix(u[,((k-1):1)], ncol = k-1), pcs, vt1, vt2), x),
                   cond_var = 1, family = pcs[[k]], vt1 = vt1, vt2 = vt2, inverse = TRUE)
    return(RforwardI(arg1, matrix(u[,(1:(k-1))], ncol = k-1), pcs, vt1, vt2))
  }
}

#' @describeIn dvinecopulavt Simulation method for dvinecopulavt class
#'
#' @param object an object of the class.
#' @param n length of realization.
#' @param forcetrunc logical parameter: TRUE truncates the copula sequence at lag 10
#' to accelerate simulation if copula sequence is longer; FALSE turns this feature off.
#'
#' @export
#'
setMethod("sim", c(object = "dvinecopulavt"), function(object, n = 1000, forcetrunc = TRUE) {
  pc_list <- mklist_dvinevt(object, n-1)
  k <- length(pc_list)
  if ((k > 10) & forcetrunc){
    k <- 10
    warning("Copula sequence has been truncated to length 10 for speed.
            Can be overridden if you are a patient person.")
  }
  vt1 <- object@modelspec$vt1
  vt2 <- object@modelspec$vt2
  sim <- numeric(n)
  sim[1] <- runif(1)
  for (t in 2:n) {
    pastvalues <- sim[(t-1):max(1, (t - k))]
    sim[t] <- RforwardI(runif(1), matrix(pastvalues, ncol = length(pastvalues)),
                        pc_list, vt1, vt2)
  }
  sim
})

#' Generalized lagging for fitted dvinecopulavt objects
#'
#' @param copula a dvinecopulavt object
#' @param data the data to which copula is fitted
#' @param lagmax the maximum lag value.
#' @param glagplot logical value indicating generalized lag plot.
#'
#' @return If \code{glagplot} is \code{TRUE} a list of generalized lagged datasets
#' of maximum length 9 is returned to facilitate a generalized lagplot.
#' If \code{glagplot} is \code{FALSE} a vector of length \code{lagmax} containing
#' the Kendall rank correlations for the generalized lagged datasets is returned.
#' @keywords internal
glag_for_dvinecopulavt <- function(copula, data, lagmax, glagplot = FALSE) {
  if (glagplot)
    lagmax <- min(lagmax, 9)
  pc_list <- mklist_dvinevt(copula, lagmax)
  k <- length(pc_list)
  n <- length(data)
  vt1 <- copula@modelspec$vt1
  vt2 <- copula@modelspec$vt2
  data <- cbind(as.numeric(data[1:(n - 1)]), as.numeric(data[2:n]))
  if (glagplot){
    output <- vector(mode = "list", length = k)
    output[[1]] <- data
  }
  else{
    output <- rep(NA, k)
    vdata <- cbind(vtrans(vt1, data[,1]),
                   vtrans(vt2, data[,2]))
    output[1] <- cor(vdata, method = "kendall")[1, 2]
  }
  if (k >1){
    for (i in 1:(k - 1)) {
      n <- dim(data)[1]
      data <-
        cbind(hbicop_vt(data[(1:(n - 1)), ], cond_var = 2, family = pc_list[[i]],
                     vt1 = vt1, vt2 = vt2),
              hbicop_vt(data[(2:n), ], cond_var = 1, family = pc_list[[i]],
                     vt1 = vt1, vt2 = vt2))
      if (glagplot)
        output[[i+1]] <- data
      else{
        vdata <- cbind(vtrans(vt1, data[,1]),
                       vtrans(vt2, data[,2]))
        output[i+1] <- cor(vdata, method = "kendall")[1, 2]
      }
    }
  }
  output
}

#' @describeIn dvinecopulavt Prediction method for dvinecopulavt class
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#'
#' @export
#'
setMethod("predict", c(object = "dvinecopulavt"), function(object, data, x, type = "df") {
  pc_list <- mklist_dvinevt(object, length(data))
  k <- length(pc_list)
  n <- length(data)
  vt1 <- object@modelspec$vt1
  vt2 <- object@modelspec$vt2
  lastkdata <- data[(n-k+1):n]
  switch(type,
         "df" = Rforward(x, matrix(rev(lastkdata), ncol = k), pc_list, vt1, vt2),
         "qf" = RforwardI(x, matrix(rev(lastkdata), ncol = k), pc_list, vt1, vt2),
         "dens" = {
            output <- rep(NA, length(x))
            for (i in 1:length(x)){
              data <- c(as.numeric(lastkdata), x[i])
              n <- length(data)
              v <- cbind(data[1:(n - 1)], data[2:n])
              output[i] <- 0
              for (j in 1:k) {
                v_vt <- cbind(vtrans(vt1, v[,1]), vtrans(vt2, v[,2]))
                tmp <- log(tscopula_dbicop(v_vt[n-1,], family = pc_list[[j]]))
                output[i] <- output[i] + sum(tmp)
                  if (j < k) {
                    n <- dim(v)[1]
                    arg1 <- matrix(v[(1:(n - 1)), ], ncol =2)
                    arg2 <- matrix(v[(2:n), ], ncol =2)
                    v <- cbind(
                 hbicop_vt(arg1, cond_var = 2, family = pc_list[[j]], vt1 = vt1, vt2 = vt2),
                 hbicop_vt(arg2, cond_var = 1, family = pc_list[[j]], vt1 = vt1, vt2 = vt2)
                  )
                }
              }
            }
           exp(output)
           })
})




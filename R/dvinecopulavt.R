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
#' @param vt1 first v-transform
#' @param vt2 second v-transform
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
                         vt1 = Vlinear(0.5),
                         vt2 = Vlinear(0.5),
                         maxlag = 15) {
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
  vt1pars <- coef(vt1)
  names(vt1pars) <- paste(names(vt1pars),"1",sep="")
  vt2pars <- coef(vt2)
  names(vt2pars) <- paste(names(vt2pars),"2",sep="")
  pars <- c(pars, vt1pars, vt2pars)
  modelspec <- list(family = fam,
                    rotation = rotation,
                    kpacf = kpacf,
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
  vt1 <- getvt(modelspec$vt1, theta, 1)
  vt2 <- getvt(modelspec$vt2, theta, 2)
  LL <- 0
  for (j in 1:k) {
    v_vt <- cbind(vtrans(vt1, v[,1]), vtrans(vt2, v[,2]))
    LL <- LL + sum(log(rvinecopulib::dbicop(u = v_vt, family = pc_list[[j]])))
    if (j == k) {
      return(-LL)
    }
    n <- dim(v)[1]
    v <- cbind(
      hbicop(v[(1:(n - 1)), ], cond_var = 2, family = pc_list[[j]], vt1, vt2),
      hbicop(v[(2:n), ], cond_var = 1, family = pc_list[[j]], vt1, vt2)
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
hbicop <- function(u, cond_var, family, vt1, vt2, inverse = FALSE) {
  if (vt1@name != "Vlinear")
    stop("first vt non-linear")
  if (vt2@name != "Vlinear")
    stop("second vt non-linear")
  delta1 <- vt1@pars["delta"]
  delta2 <- vt2@pars["delta"]
  cond1 <- u[,2] > delta1
  cond2 <- u[,1] > delta2
  mult1 <- -(delta2^(1-cond1))*(delta2 - 1)^cond1
  mult2 <- -(delta1^(1-cond2))*(delta1 - 1)^cond2
  v <- cbind(vtrans(vt1, u[,1]), vtrans(vt2, u[,2]))
  switch(cond_var,
           mult1 * rvinecopulib::hbicop(u = v, cond_var = 1,
                                        family = family, inverse = inverse) + delta2,
           mult2 * rvinecopulib::hbicop(u = v, cond_var = 2,
                                        family = family, inverse = inverse) + delta1
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
  n <- length(data)
  pc_list <- mklist_dvine2(object, n-1, truncate = TRUE, tol = 1/3)
  k <- length(pc_list)
  if (trace)
    target <- rep(0.5, n)
  else
    target <- data
  res <- rep(NA, n)
  delta <- as.numeric(unlist(object@pars)[c("delta1", "delta2")])
  res[1] <- target[1]
  for (i in 2:k)
    res[i] <- Rforward(data[i], matrix(data[(i-1):1], ncol = i-1, nrow =1), pc_list, delta)
  pastdata <- matrix(NA, ncol = k, nrow = n - k)
  for (i in ((k+1):n))
    pastdata[i-k,] <- data[(i-1):(i-k)]
  res[(k+1):n] <- Rforward(data[(k+1):n], pastdata, pc_list, delta)
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
#' @param delta vector of fulcrum parameters
#'
#' @return vector of same length as x
#'
#' @keywords internal
#'
Rforward <- function(x, u, pcs, delta) {
  if (!(is.matrix(u)))
    stop("Must have matrix of conditioning variables")
  k <- dim(u)[2]
  if (k == 1) {
    return(as.numeric(hbicop(cbind(as.vector(u), x), 1, pcs[[1]], delta)))
  } else {
    return(hbicop(cbind(
      Rbackward(u[, k], matrix(u[, ((k-1):1)], ncol = k-1), pcs, delta),
      Rforward(x, matrix(u[, (1:(k-1))], ncol = k-1), pcs, delta)
    ), 1, pcs[[k]], delta))
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
#' @param delta vector of fulcrum parameters
#'
#' @return vector of same length as x
#'
Rbackward <- function(x, u, pcs, delta) {
  if (!(is.matrix(u)))
    stop("Must have matrix of conditioning variables")
  k <- dim(u)[2]
  if (k == 1) {
    return(as.numeric(hbicop(cbind(x, as.vector(u)), 2, pcs[[1]], delta)))
  } else {
    return(hbicop(cbind(
      Rbackward(x, matrix(u[, (1:(k-1))], ncol = k-1), pcs, delta),
      Rforward(u[, k], matrix(u[, ((k-1):1)], ncol = k-1), pcs, delta)
    ), 2, pcs[[k]], delta))
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
#' @param delta vector of fulcrum parameters
#'
#' @return vector of same length as x
#'
RforwardI <- function(x, u, pcs, delta) {
  if (!(is.matrix(u)))
    stop("Must have matrix of conditioning variables")
  k <- dim(u)[2]
  if (k == 1) {
    return(as.numeric(hbicop(cbind(as.vector(u), x),
                             1, pcs[[1]], delta, inverse = TRUE)))
  } else {
    arg1 <- hbicop(cbind(Rbackward(u[,k], matrix(u[,((k-1):1)], ncol = k-1), pcs, delta),x),
                   1, pcs[[k]], delta, inverse = TRUE)
    return(RforwardI(arg1, matrix(u[,(1:(k-1))], ncol = k-1), pcs, delta))
  }
}

#' @describeIn dvinecopulavt Simulation method for dvinecopulavt class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
setMethod("sim", c(object = "dvinecopulavt"), function(object, n = 1000, forcetrunc = TRUE) {
  pc_list <- mklist_dvine2(object, n-1, truncate = TRUE, tol = 1/3)
  k <- length(pc_list)
  if ((k > 10) & forcetrunc){
    k <- 10
    warning("Copula sequence has been truncated to length 10 for speed.
            Can be overridden if you are a patient person.")
  }
  delta <- as.numeric(unlist(object@pars)[c("delta1", "delta2")])
  sim <- numeric(n)
  sim[1] <- runif(1)
  for (t in 2:n) {
    pastvalues <- sim[(t-1):max(1, (t - k))]
    sim[t] <- RforwardI(runif(1), matrix(pastvalues, ncol = length(pastvalues)),
                        pc_list, delta)
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
  pc_list <- mklist_dvine2(copula, lagmax, truncate = FALSE)
  k <- length(pc_list)
  n <- length(data)
  delta <- as.numeric(unlist(copula@pars)[c("delta1", "delta2")])
  data <- cbind(as.numeric(data[1:(n - 1)]), as.numeric(data[2:n]))
  if (glagplot){
    output <- vector(mode = "list", length = k)
    output[[1]] <- data
  }
  else{
    output <- rep(NA, k)
    vdata <- cbind(vtrans(Vlinear(delta[1]), data[,1]),
                   vtrans(Vlinear(delta[2]), data[,2]))
    output[1] <- cor(vdata, method = "kendall")[1, 2]
  }
  if (k >1){
    for (i in 1:(k - 1)) {
      n <- dim(data)[1]
      data <-
        cbind(hbicop(data[(1:(n - 1)), ], pc_list[[i]], cond_var = 2, delta),
              hbicop(data[(2:n), ], pc_list[[i]], cond_var = 1, delta))
      if (glagplot)
        output[[i+1]] <- data
      else{
        vdata <- cbind(vtrans(Vlinear(delta[1]), data[,1]),
                       vtrans(Vlinear(delta[2]), data[,2]))
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
  pc_list <- mklist_dvine2(object, length(data)-1, truncate = TRUE, tol = 1/3)
  k <- length(pc_list)
  n <- length(data)
  delta <- as.numeric(coef(object)[c("delta1", "delta2")])
  lastkdata <- data[(n-k+1):n]
  switch(type,
         "df" = Rforward(x, matrix(rev(lastkdata), ncol = k), pc_list, delta),
         "qf" = RforwardI(x, matrix(rev(lastkdata), ncol = k), pc_list, delta),
         "dens" = {
            output <- rep(NA, length(x))
            for (i in 1:length(x)){
              data <- c(as.numeric(lastkdata), x[i])
              n <- length(data)
              v <- cbind(data[1:(n - 1)], data[2:n])
              output[i] <- 0
              for (j in 1:k) {
                v_vt <- cbind(vtrans(Vlinear(delta[1]), v[,1]), vtrans(Vlinear(delta[2]), v[,2]))
                tmp <- log(rvinecopulib::dbicop(v_vt[n-1,], family = pc_list[[j]]))
                output[i] <- output[i] + sum(tmp)
                  if (j < k) {
                    n <- dim(v)[1]
                    arg1 <- matrix(v[(1:(n - 1)), ], ncol =2)
                    arg2 <- matrix(v[(2:n), ], ncol =2)
                    v <- cbind(
                 hbicop(arg1, cond_var = 2, family = pc_list[[j]], delta),
                 hbicop(arg2, cond_var = 1, family = pc_list[[j]], delta)
                  )
                }
              }
            }
           exp(output)
           })
})




#' Type 2 D-vine Copula Processes
#'
#' Class of objects for d-vine copula processes.
#'
#' @slot name name of the d-vine copula process.
#' @slot modelspec list containing the family, rotation, and name of KPACF
#' @slot pars list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("dvinecopula2", contains = "tscopula", slots = list(
  name = "character",
  modelspec = "list",
  pars = "list"
))

#' Constructor Function for dvinecopula2 process
#'
#' @param family family name
#' @param rotation rotation
#' @param kapcf
#' @param pars a list containing the parameters of each lag
#' @param rotation a vector of rotations
#'
#' @return An object of class \linkS4class{dvinecopula}.
#' @export
#'
#' @examples
dvinecopula2 <- function(family = "gauss",
                         rotation = 0,
                         kpacf = "kpacf_arma",
                         pars = list(ar = 0.1, ma = 0.1),
                         maxlag = Inf) {
  if (class(family) != "character")
    stop("copula family must be specified by name")
 modelspec <- list(family = tolower(family),
                   rotation = rotation,
                   kpacf = kpacf,
                   maxlag = maxlag)
  new("dvinecopula2",
      name = paste("type2-d-vine"),
      modelspec = modelspec,
      pars = pars
  )
}

#' Title
#'
#' @param k
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
kpacf_arma <- function(k, theta){
  if (is.list(theta))
    theta <- tsunlist(theta)
  ar <- numeric()
  ma <- numeric()
  nm <- substring(names(theta), 1, 2)
  if ("ar" %in% nm)
    ar <- theta[nm == "ar"]
  if ("ma" %in% nm)
    ma <- theta[nm == "ma"]
  if ((non_stat(ar)) | (non_invert(ma)))
    return(rep(NA, length(k)))
  rho <- ARMAacf(ar = ar, ma = ma, lag.max = max(k), pacf = TRUE)[k]
  (2/pi)*asin(rho)
}

#' Title
#'
#' @param k
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
kpacf_exp <- function(k, theta){
  if (is.list(theta))
    theta <- tsunlist(theta)
  arg <- theta[1] + theta[2] * k
  arg <- pmin(arg, 0)
  exp(arg)
}

#' Coef Method for dvinecopula2 Class
#'
#' @param object an object of class \linkS4class{dvinecopula2}.
#'
#' @return parameters of tscopula model
#' @export
#'
setMethod("coef", c(object = "dvinecopula2"), function(object) {
  if (length(object@pars) == 1) {
    return(object@pars[[1]])
  } else {
    return(unlist(object@pars))
  }
})

#' Show Method for dvinecopula2 Class
#'
#' @param object an object of class \linkS4class{dvinecopula2}.
#'
#' @return
#' @export
#'
setMethod("show", c(object = "dvinecopula2"), function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  famname <- object@modelspec$family
  if (object@modelspec$rotation !=0)
    famname <- paste(famname, "with rotation", object@modelspec$rotation)
  cat("copula family: ", famname, "\n", sep = "")
  kpacf  <- object@modelspec$kpacf
  if (object@modelspec$maxlag != Inf)
    kpacf <- paste(kpacf, "with max lag", object@modelspec$maxlag)
  cat("KPACF: ", kpacf,"\n", sep = "")
  cat("parameters: \n")
  print(coef(object))
})


#' Title
#'
#' @param theta
#' @param modelspec
#' @param u
#'
#' @return
#' @export
#'
#' @examples
dvinecopula2_objective <- function(theta, modelspec, u) {
  n <- length(u)
  kpacf <- eval(parse(text = modelspec$kpacf))
  tauvals <- kpacf(1:(n-1), theta)
  if (is.na(sum(tauvals)))
    return(NA)
  k <- 1
  largetau <- (abs(tauvals) > .Machine$double.eps)
  if (sum(largetau) > 0)
    k <- max((1:(n-1))[largetau])
  k <- min(k, modelspec$maxlag)
  pc_list <- vector("list", k)
  for (i in 1:k) {
    pc_list[[i]] <- tryCatch(rvinecopulib::bicop_dist(
      family = tolower(modelspec$family),
      rotation = modelspec$rotation,
      parameters = ktau_to_par(
        family = modelspec$family,
        tau = tauvals[i]
      )
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
  LL <- 0
  for (j in 1:k) {
    LL <- LL + sum(log(rvinecopulib::dbicop(u = v, family = pc_list[[j]])))
    if (j == k) {
      return(-LL)
    }
    n <- dim(v)[1]
    v <- cbind(
      rvinecopulib::hbicop(v[(1:(n - 1)), ], cond_var = 2, family = pc_list[[j]]),
      rvinecopulib::hbicop(v[(2:n), ], cond_var = 1, family = pc_list[[j]])
    )
  }
}

#' Title
#'
#' @param family
#' @param tau
#'
#' @return
#' @export
#'
#' @examples
ktau_to_par <- function(family, tau){
  if (family %in% c("joe", "gumbel", "clayton"))
    if (tau < 0)
      stop("Negative tau not allowed for this family")
  rvinecopulib::ktau_to_par(family, tau)
}

#' Simulation Method for dvinecopula2 Class
#'
#' @param x an object of class \linkS4class{dvinecopula2}.
#' @param n length of realization.
#'
#' @return A realization of a time series copula process.
#' @export
#'
#' @examples
setMethod("sim", c(x = "dvinecopula2"), function(x, n = 1000) {
  kpacf <- eval(parse(text = x@modelspec$kpacf))
  tauvals <- kpacf(1:(n-1), x@pars)
  kmax <- max((1:(n-1))[abs(tauvals) > .Machine$double.eps])
  k <- min(max(1, kmax), x@modelspec$maxlag)
  pc_list <- vector("list", k)
  for (i in 1:k) {
    pc_list[[i]] <- rvinecopulib::bicop_dist(
      family = tolower(x@modelspec$family),
      rotation = x@modelspec$rotation,
      parameters = rvinecopulib::ktau_to_par(
        family = x@modelspec$family,
        tau = tauvals[i]))
  }
  # bring pair-copulas in appropriate form for `rvinecopulib::vinecop_dist()`
  pcs <- lapply(seq_along(pc_list), function(i) {
    replicate(k - i + 1, pc_list[[i]], simplify = FALSE)
  })
  # set up k + 1 dimensional vine copula model
  vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure((k + 1):1))
  # initialize first steps of simulation
  sim <- numeric(n)
  u <- rvinecopulib::rvinecop(1, vc_short)
  sim[1:(k + 1)] <- u
  # conditional simulation for all future time points
  for (t in (k + 2):n) {
    u_new <- c(sim[(t - k):(t - 1)], 0.5) # (0.5 is a dummy)
    pit <- rvinecopulib::rosenblatt(t(u_new), vc_short)
    sim[t] <- rvinecopulib::inverse_rosenblatt(t(c(pit[-(k + 1)], runif(1))), vc_short)[k + 1]
  }
  sim
})

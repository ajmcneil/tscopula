#' Non-Exchangeable d-vine copula processes
#'
#' Class of objects for d-vine copula processes.
#'
#' @slot name name of the d-vine copula process.
#' @slot modelspec list containing the family, number of parameters, rotations and exchangeability indicator.
#' @slot pars list comprising of the parameters.
#'
#' @return
#' @export
#'
setClass("dvinecopulaNE", contains = "dvinecopula", slots = list(
  name = "character",
  modelspec = "list",
  pars = "list"
))

#' Constructor Function for dvinecopulaNE process
#'
#' @param family a vector of family names
#' @param pars a list containing the parameters of each lag
#' @param rotation a vector of rotations
#' @param exchangeable a vector of exchangeability indicators
#'
#' @return An object of class \linkS4class{dvinecopulaNE}.
#'
#' @export
#'
#' @examples
#' dvinecopulaNE(family = c("joe", "gauss", "t"), pars = list(3, .5, c(1, 2)), rotation = c(180, 0, 0), exchangeable = c(T, F, T))
dvinecopulaNE <- function(family = "indep", pars = list(NULL), rotation = 0, exchangeable = TRUE) {
  k <- length(pars)
  if ((length(family) > 1) & (k != length(family))) {
    stop("specify single family or vector with same length as parameter list")
  }
  if ((length(rotation) > 1) & (k != length(rotation))) {
    stop("specify single rotation value or vector with same length as parameter list")
  }
  if ((length(exchangeable) > 1) & (k != length(exchangeable))) {
    stop("specify single value for exchangeable or vector with same length as parameter list")
  }
  if (length(family) == 1) {
    family <- rep(family, k)
  }
  if (length(rotation) == 1) {
    rotation <- rep(rotation, k)
  }
  if (length(exchangeable) == 1) {
    exchangeable <- rep(exchangeable, k)
  }
  modelspec <- list()
  for (i in seq_along(family)) {
    modelspec[[i]] <- list(
      family = family[i], npars = length(pars[[i]]),
      rotation = rotation[i],
      exchangeable = exchangeable[i]
    )
  }
  pars <- lapply(pars, function(v) {
    if (length(v) > 0) {
      names(v) <- paste("p", (1:length(v)), sep = "")
    }
    v
  })
  names(pars) <- paste("cop", 1:k, sep = "")
  new("dvinecopulaNE",
    name = paste("d-vine(", k, ")", sep = ""),
    modelspec = modelspec,
    pars = pars
  )
}

#' Objective Function for non-exchangeable dvinecopula process
#'
#' @param theta parameters of copulas
#' @param modelspec list of families of copulas
#' @param udata data
#' @return
#' @keywords internal
#'
dvinecopulaNE_objective <- function(theta, modelspec, u) {
  k <- length(modelspec)
  npars <- sapply(modelspec, function(v) {
    v$npars
  })
  pars <- split(theta, rep(1:k, npars))
  n <- length(u)
  pc_list <- vector("list", k)
  extrapars <- vector("list", k)
  for (i in seq_along(modelspec)) {
    parsi <- pars[[i]]
    extrapars[[i]] <- c(NA, NA)
    if (!modelspec[[i]]$exchangeable) {
      thetatilde <- parsi[length(parsi)]
      if (thetatilde >= 0) {
        extrapars[[i]] <- c(1 - exp(-thetatilde), 0)
      }
      if (thetatilde < 0) {
        extrapars[[i]] <- c(0, 1 - exp(thetatilde))
      }
      parsi <- parsi[-length(parsi)]
    }
    pc_list[[i]] <- tryCatch(rvinecopulib::bicop_dist(
      family = tolower(modelspec[[i]]$family),
      rotation = modelspec[[i]]$rotation,
      parameters = parsi
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
    LL <- LL + sum(log(dbicop(u = v, family = pc_list[[j]], extrapars[[j]])))
    if (j == k) {
      return(-LL)
    }
    n <- dim(v)[1]
    v <- cbind(
      hbicop(v[(1:(n - 1)), ], cond_var = 2, family = pc_list[[j]], extrapars[[j]]),
      hbicop(v[(2:n), ], cond_var = 1, family = pc_list[[j]], extrapars[[j]])
    )
  }
}

#' Density for non-exchangeable pair copulas
#'
#' @param u
#' @param family
#' @param theta
#'
#' @return
#' @export
#'
dbicop <- function(u, family, theta = c(NA, NA)) {
  if (is.na(sum(theta))) {
    rvinecopulib::dbicop(u = u, family = family)
  } else {
    v1 <- u[, 1]^(1 - theta[1])
    v2 <- u[, 2]^(1 - theta[2])
    v <- cbind(v1, v2)
    rvinecopulib::dbicop(u = v, family = family) * (1 - theta[1]) * (1 - theta[2]) +
      rvinecopulib::hbicop(u = v, cond_var = 1, family = family) * theta[2] * (1 - theta[1]) / v2 +
      rvinecopulib::hbicop(u = v, cond_var = 2, family = family) * theta[1] * (1 - theta[2]) / v1 +
      rvinecopulib::pbicop(u = v, family = family) * theta[1] * theta[2] / (v1 * v2)
  }
}

#' H-function for non-exchangeable pair copulas
#'
#' @param u
#' @param cond_var
#' @param family
#' @param theta
#'
#' @return
#' @export
#'
hbicop <- function(u, cond_var, family, theta = c(NA, NA)) {
  if (is.na(sum(theta))) {
    rvinecopulib::hbicop(u = u, cond_var = cond_var, family = family)
  } else {
    v1 <- u[, 1]^(1 - theta[1])
    v2 <- u[, 2]^(1 - theta[2])
    v <- cbind(v1, v2)
    switch(cond_var,
           rvinecopulib::hbicop(u = v, cond_var = 1, family = family) * (1 - theta[1]) * u[, 2]^theta[2] +
             rvinecopulib::pbicop(u = v, family = family) * theta[1] * u[, 2] / (v1 * v2),
           rvinecopulib::hbicop(u = v, cond_var = 2, family = family) * (1 - theta[2]) * u[, 1]^theta[1] +
             rvinecopulib::pbicop(u = v, family = family) * theta[2] * u[, 1] / (v1 * v2)
    )
  }
}

#' Rosenblatt Transforms for non-exchangeable copulas
#'
#' @param latest
#' @param previous
#' @param cv
#' @param pcs
#' @param extrapars
#'
#' @return
#' @export
#'
RT <- function(latest, previous, cv, pcs, extrapars) {
  if ((latest == 0) | (latest == 1)) {
    return(latest)
  }
  k <- length(previous)
  x <- c(previous, latest)
  if (k == 1) {
    return(hbicop(cbind(x[1], x[2]), cv, pcs[[1]], extrapars[[1]]))
  } else {
    return(hbicop(cbind(
      RT(x[k], x[1:(k - 1)], cv = 2, pcs, extrapars),
      RT(x[k + 1], x[2:k], cv = 1, pcs, extrapars)
    ), cv, pcs[[k]], extrapars[[k]]))
  }
}

RT_vectorized <- Vectorize(RT, c("latest"))

#' Inverse Rosenblatt transform  for non-exchangeable pair copulas
#'
#' @param U
#' @param previous
#' @param pcs
#' @param extrapars
#'
#' @return
#' @export
#'
RTinverse <- function(U, previous, pcs, extrapars, tol = .Machine$double.eps^0.75) {
  root <- uniroot(function(x, U, previous, pcs, extrapars, tol) {
    U - RT_vectorized(x, previous, 1, pcs, extrapars)
  }, interval = c(0, 1), U = U, previous = previous, pcs = pcs, extrapars = extrapars, tol = tol)$root
  if (root <= 0) {
    root <- tol
  }
  if (root >= 1) {
    root <- 1 - tol
  }
  root
}

#' Simulation Method for dvinecopulaNE Class
#'
#' @param x an object of class \linkS4class{dvinecopulaNE}.
#' @param n length of realization.
#'
#' @return A realization of a time series copula process.
#' @export
#'
#' @examples
#' sim(dvinecopulaNE("gauss", list(c(0.5, 2)), exchangeable = F))
setMethod("sim", c(x = "dvinecopulaNE"), function(x, n = 1000) {
  k <- length(x@modelspec)
  pc_list <- vector("list", k)
  extrapars <- vector("list", k)
  for (i in seq_along(x@modelspec)) {
    extrapars[[i]] <- c(0, 0)
    parsi <- x@pars[[i]]
    if (!x@modelspec[[i]]$exchangeable) {
      thetatilde <- parsi[length(parsi)]
      if (thetatilde >= 0) {
        extrapars[[i]] <- c(1 - exp(-thetatilde), 0)
      }
      if (thetatilde < 0) {
        extrapars[[i]] <- c(0, 1 - exp(thetatilde))
      }
      parsi <- parsi[-length(parsi)]
    }
    pc_list[[i]] <- rvinecopulib::bicop_dist(
      family = tolower(x@modelspec[[i]]$family),
      rotation = x@modelspec[[i]]$rotation,
      parameters = parsi
    )
  }
  sim <- numeric(n)
  sim[1] <- runif(1)
  for (t in 2:n) {
    sim[t] <- RTinverse(runif(1), sim[max(1, (t - k)):(t - 1)], pc_list, extrapars)
  }
  sim
})




#' D-vine Copula Processes
#'
#' Class of objects for d-vine copula processes.
#'
#' @slot name name of the d-vine copula process.
#' @slot modelspec list containing the family, number of parameters and rotations
#' @slot pars list comprising of the parameters.
#'
#' @export
#'
setClass("dvinecopula", contains = "tscopula", slots = list(
  name = "character",
  modelspec = "list",
  pars = "list"
))

#' Constructor Function for dvinecopula process
#'
#' @param family a vector of family names
#' @param pars a list containing the parameters of each lag
#' @param rotation a vector of rotations
#'
#' @return An object of class \linkS4class{dvinecopula}.
#' @export
#'
#' @examples
#' dvinecopula(family = c("joe", "gauss", "t"), pars = list(3, .5, c(1, 2)), rotation = c(180, 0, 0))
dvinecopula <- function(family = "indep", pars = list(NULL), rotation = 0) {
  if (class(family) != "character")
    stop("families must be specified by names")
  else
    family <- tolower(family)
  k <- length(pars)
  if ((length(family) > 1) & (k != length(family))) {
    stop("specify single family or vector with same length as parameter list")
  }
  if ((length(rotation) > 1) & (k != length(rotation))) {
    stop("specify single rotation value or vector with same length as parameter list")
  }
  if (length(family) == 1) {
    family <- rep(family, k)
  }
  if (length(rotation) == 1) {
    rotation <- rep(rotation, k)
  }
  modelspec <- list()
  for (i in seq_along(family)) {
    tmp <- rvinecopulib::bicop_dist(family[i], rotation[i], pars[[i]])
    modelspec[[i]] <- list(
      family = family[i], npars = length(pars[[i]]),
      rotation = rotation[i]
    )
  }
  pars <- lapply(pars, function(v) {
    if (length(v) > 0) {
      names(v) <- paste("p", (1:length(v)), sep = "")
    }
    v
  })
  names(pars) <- paste("cop", 1:k, sep = "")
  new("dvinecopula",
      name = paste("d-vine(", k, ")", sep = ""),
      modelspec = modelspec,
      pars = pars
  )
}

#' Coef Method for dvinecopula Class
#'
#' @param object an object of class \linkS4class{dvinecopula}.
#'
#' @return parameters of tscopula model
#' @export
#'
setMethod("coef", c(object = "dvinecopula"), function(object) {
  if (length(object@pars) == 1) {
    return(object@pars[[1]])
  } else {
    return(unlist(object@pars))
  }
})

#' Show Method for dvinecopula Class
#'
#' @param object an object of class \linkS4class{dvinecopula}.
#'
#' @return A summary of an object of class \linkS4class{dvinecopula}.
#' @export
#'
setMethod("show", "dvinecopula", function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  fams <- sapply(object@modelspec, FUN = function(v) {
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
  cat("parameters: \n")
  print(coef(object))
})

#' Simulation Method for dvinecopula Class
#'
#' @param x an object of class \linkS4class{dvinecopula}.
#' @param n length of realization.
#'
#' @return A realization of a time series copula process.
#' @export
#'
#' @examples
#' sim(dvinecopula("gauss", 0.5))
setMethod("sim", c(x = "dvinecopula"), function(x, n = 1000) {
  pc_list <- mklist_dvine(x)
  simdvine(pc_list, n)
})

#' Make list of pair copulas for dvinecopula object
#'
#' @param x an object of class dvinecopula
#'
#' @return a list of pair copulas
#' @keywords internal
#'
mklist_dvine <- function(x){
  k <- length(x@modelspec)
  pc_list <- vector("list", k)
  for (i in seq_along(x@modelspec)) {
    pc_list[[i]] <- rvinecopulib::bicop_dist(
      family = tolower(x@modelspec[[i]]$family),
      rotation = x@modelspec[[i]]$rotation,
      parameters = x@pars[[i]][1:x@modelspec[[i]]$npars]
    )
  }
  pc_list
}

#' D-vine simulation helper function
#'
#' @param pc_list a list of pair copulas.
#' @param n number of data to be simulated.
#'
#' @return a vector of length n.
#' @export
#'
simdvine <- function(pc_list, n){
  # code template provided by Thomas Nagler
  k <- length(pc_list)
  pcs <- lapply(seq_along(pc_list), function(i) {
    replicate(k - i + 1, pc_list[[i]], simplify = FALSE)
  })
  vc_short <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure((k + 1):1))
  Z <- runif(n)
  U <- numeric(n)
  U[1:(k + 1)] <- rvinecopulib::inverse_rosenblatt(t(Z[1:(k+1)]), vc_short)
  for (t in (k + 2):n) {
    lastvals <- c(U[(t - k):(t - 1)], 0.5)
    rt <- rvinecopulib::rosenblatt(t(lastvals), vc_short)
    rt[k+1] <- Z[t]
    irt <- rvinecopulib::inverse_rosenblatt(t(rt), vc_short)
    U[t] <- irt[k + 1]
  }
  U
}


#' Objective Function for dvinecopula process
#'
#' @param theta parameters of copulas
#' @param modelspec list of families of copulas
#' @param udata data
#'
#' @return Value of objective function at parameters.
#'
#' @keywords internal
#'
dvinecopula_objective <- function(theta, modelspec, u) {
  k <- length(modelspec)
  npars <- sapply(modelspec, function(v) {
    v$npars
  })
  pars <- split(theta, rep(1:k, npars))
  n <- length(u)
  pc_list <- vector("list", k)
  for (i in seq_along(modelspec)) {
    pc_list[[i]] <- tryCatch(rvinecopulib::bicop_dist(
      family = tolower(modelspec[[i]]$family),
      rotation = modelspec[[i]]$rotation,
      parameters = pars[[i]]
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

#' Generalized lagging for fitted dvinecopula objects
#'
#' @param copula a dvinecopula object
#' @param data the data to which copula is fitted
#' @param lagmax the maximum lag value.
#' @param glagplot logical value indicating generalized lag plot.
#'
#' @keywords internal
glag_for_dvinecopula <- function(copula, data, lagmax, glagplot = FALSE) {
  pc_list <- mklist_dvine(copula)
  k <- min(length(pc_list), lagmax)
  if (glagplot)
    k <- min(k, 9)
  n <- length(data)
  data <- cbind(as.numeric(data[1:(n - 1)]), as.numeric(data[2:n]))
  if (glagplot){
    output <- vector(mode = "list", length = k)
    output[[1]] <- data
  }
  else{
  output <- rep(NA, k)
  output[1] <- cor(data, method = "kendall")[1, 2]
  }
  if (k > 1){
  for (i in 1:(k - 1)) {
    n <- dim(data)[1]
    data <-
      cbind(rvinecopulib::hbicop(data[(1:(n - 1)), ], pc_list[[i]], cond_var = 2),
            rvinecopulib::hbicop(data[(2:n), ], pc_list[[i]], cond_var = 1))
    if (glagplot)
      output[[i+1]] <- data
    else
      output[i+1] <- cor(data, method = "kendall")[1, 2]
  }
  }
  output
}


#' Calculate Kendall's tau values for pair copulas in d-vine copula
#'
#' @param x a \linkS4class{dvinecopula} object
#' @param lagmax maximum value of lag
#'
#' @return vector consisting of Kendall's tau values for each pair copula
#' @export
#'
#' @examples
#' mixmod <- dvinecopula(family = c("gumbel", "gauss"), pars = list(1.5, -0.6))
#' kendall(mixmod)
setMethod("kendall", c(x = "dvinecopula"), function(x, lagmax = 20) {
  tau <- vector("numeric",length(x@modelspec))
  if (length(tau) > lagmax)
    tau <- tau[1:lagmax]
  for (i in seq_along(tau)){
    model <- rvinecopulib::bicop_dist(
      family = tolower(x@modelspec[[i]]$family),
      rotation = x@modelspec[[i]]$rotation,
      parameters = x@pars[[i]][1:x@modelspec[[i]]$npars]
    )
    tau[i] <- rvinecopulib::par_to_ktau(model)
  }
  names(tau) <- sapply(x@modelspec, function(v){v$family})
  tau
}
)

#' Residual function for dvinecopula object
#'
#' @param object a fitted dvinecopula object.
#' @param data the data to which copula is fitted.
#' @param trace extract trace instead of residuals.
#'
#' @return vector of model residuals
#' @keywords internal
#'
resid_dvinecopula <- function(object, data = NA, trace = FALSE){
  pc_list <- mklist_dvine(object)
  k <- length(pc_list)
  n <- length(data)
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

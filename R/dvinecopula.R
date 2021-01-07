#' D-vine Copula Processes
#'
#' Class of objects for d-vine copula processes.
#'
#' @slot name name of the d-vine copula process.
#' @slot modelspec list containing the family, number of parameters and rotations
#' @slot pars list comprising of the parameters.
#'
#' @return
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
#' @return
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
  # code from Thomas Nagler
  k <- length(x@modelspec)
  pc_list <- vector("list", k)
  for (i in seq_along(x@modelspec)) {
    pc_list[[i]] <- rvinecopulib::bicop_dist(
      family = tolower(x@modelspec[[i]]$family),
      rotation = x@modelspec[[i]]$rotation,
      parameters = x@pars[[i]][1:x@modelspec[[i]]$npars]
    )
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

#' Objective Function for dvinecopula process
#'
#' @param theta parameters of copulas
#' @param modelspec list of families of copulas
#' @param udata data
#' @return
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

#' Plot Function for dvinecopula Objects
#'
#' @param copula a fitted dvinecopula object
#' @param data the data to which copula is fitted
#' @param plotoption number giving plot choice
#' @param bw logical for black-white plot
#'
#' @return
#' @export
plot_dvinecopula <- function(copula, data, plotoption, bw){
  data0 <- data
  n <- length(data)
  k <- length(copula@modelspec)
  kplotmax <- min(k, 9)
  datavecs <- vector(mode = "list", length = kplotmax)
  tau_empirical <- rep(NA, k)
  data <- cbind(data[1:(n - 1)], data[2:n])
  datavecs[[1]] <- data
  tau_empirical[1] <- cor(data, method = "kendall")[1, 2]
  for (i in 1:(k - 1)) {
    n <- dim(data)[1]
    model <- rvinecopulib::bicop_dist(
      family = tolower(copula@modelspec[[i]]$family),
      rotation = copula@modelspec[[i]]$rotation,
      parameters = copula@pars[[i]][1:copula@modelspec[[i]]$npars]
    )
    data <-
      cbind(rvinecopulib::hbicop(data[(1:(n - 1)), ], model, cond_var = 2),
            rvinecopulib::hbicop(data[(2:n), ], model, cond_var = 1))
    tau_empirical[i+1] <- cor(data, method = "kendall")[1, 2]
    if (i < kplotmax)
      datavecs[[i+1]] <- data
  }
  tau_theoretical <- get_tau(copula)
  colchoice <- ifelse(bw, "gray50", "red")
  switch(plotoption,
         {
           plot(1:k, tau_empirical, xlab ="k", ylab = "tau",
                ylim = range(tau_empirical, tau_theoretical), type = "h")
           lines(1:k, tau_theoretical, col = colchoice)
           abline(h = 0)},
         {
           pacf0 <- pacf(qnorm(data0), plot = FALSE)
           plot(2/pi * asin(pacf0$acf[,,1]), type = "h", ylab = "Gaussian KPACF",
                ylim = range(tau_theoretical, tau_empirical, pacf0$acf))
           lines(1:k, tau_theoretical, col = colchoice)
           lines(1:k, tau_empirical)
           abline(h = 0)
         },
         {
           lc <- ifelse(kplotmax > 4, 3, 2)
           lr <- ceiling(kplotmax / lc)
           default_par <- par(mfrow = c(lr, lc), mar = c(2.1, 2.1, 1.5, 0.5), oma = rep(2, 4),
               pty = "s", cex = 0.5)
           for (i in 1:kplotmax)
             plot(datavecs[[i]], main = paste("Lag ", i, sep = ""), asp = 1,
                  xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
           par(default_par)
         },
         stop("Not a plot option")
  )
}

#' Calculate Kendall's tau values for pair copulas in d-vine copula
#'
#' @param vinemodel a \linkS4class{dvinecopula} object
#'
#' @return vector consisting of Kendall's tau values for each pair copula
#' @export
#'
#' @examples
#' mixmod <- dvinecopula(family = c("gumbel", "gauss", "joe", "clayton"), pars = list(1.5, -0.6, 1.6, 2.1))
#' get_tau(mixmod)
get_tau <- function(vinemodel){
  if (is(vinemodel, "tscopulafit"))
    vinemodel <- vinemodel@tscopula
  if (!is(vinemodel, "dvinecopula"))
    stop("This function is for d-vine copulas")
  tau <- vector("numeric",length(vinemodel@modelspec))
  for (i in seq_along(tau)){
    model <- rvinecopulib::bicop_dist(
      family = tolower(vinemodel@modelspec[[i]]$family),
      rotation = vinemodel@modelspec[[i]]$rotation,
      parameters = vinemodel@pars[[i]][1:vinemodel@modelspec[[i]]$npars]
    )
    tau[i] <- rvinecopulib::par_to_ktau(model)
  }
  names(tau) <- sapply(vinemodel@modelspec, function(v){v$family})
  tau
}

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

#' Constructor Function for dvinecopula2 Process
#'
#' @param family family name
#' @param rotation rotation
#' @param kpacf character string giving name of Kendal pacf
#' @param pars a list containing the parameters of each lag
#' @param maxlag scalar specifying maximum lag
#'
#' @return An object of class \linkS4class{dvinecopula}.
#' @export
#'
#' @examples
#' dvinecopula2(family = "joe", kpacf = "kpacf_arma",
#' pars = list(ar = 0.95, ma = -0.85), maxlag = 30)
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
                   maxlag = maxlag,
                   npar = length(unlist(pars)))
  new("dvinecopula2",
      name = paste("type2-d-vine"),
      modelspec = modelspec,
      pars = pars
  )
}

#' KPACF of ARMA Process
#'
#' @param k vector of lags
#' @param theta list with components ar and ma specifying the ARMA parameters
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
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
    return(rep(NA, k))
  pacf <- ARMAacf(ar = ar, ma = ma, lag.max = k, pacf = TRUE)
  (2/pi)*asin(pacf)
}

#' Compute Partial Autocorrelations from Autocorrelations
#'
#' @param rho vector of autocorrelation values (excluding 1)
#'
#' @return vector of partial autocorrelation values
#' @export
#'
#' @examples
#' rho <- ARMAacf(ar = -0.9, ma = 0.8, lag.max = 50)[-1]
#' alpha <- acf2pacf(rho)
acf2pacf <- function(rho){
  drop(.Call(stats:::C_pacf1, c(1,rho), lag.max = length(rho)))
}

#' Compute Autocorrelations from Partial Autocorrelations
#'
#' @param alpha vector of partial autocorrelation values
#'
#' @return vector of autocorrelation values
#' @export
#'
#' @examples
#' alpha <- ARMAacf(ar = -0.9, ma = 0.8, lag.max = 50, pacf = TRUE)
#' rho <- pacf2acf(alpha)
pacf2acf <- function(alpha){
  n <- length(alpha)
  rho <- rep(alpha[1], n)
  if (n > 1){
    rho[2] <- rho[1]^2 + alpha[2]*(1-rho[1]^2)
    if (n >2){
      for (k in 3:length(alpha)){
        M <- diag(rep(1,(k-1)))
        for (i in 1:(k-2))
          for (j in (i+1):(k-1)){
            M[i,j] <- rho[j-i]
            M[j,i] <- M[i,j]
          }
        Mi <- solve(M)
        v <- rho[1:(k-1)]
        term1 <- t(v) %*% Mi %*% rev(v)
        D1 <- t(v) %*% Mi %*% v
        D2 <- t(rev(v)) %*% Mi %*% rev(v)
        rho[k] <- term1 + alpha[k]*sqrt((1-D1)*(1-D2))
      }
    }
  }
  rho
}

#' KPACF of ARFIMA Process
#'
#' @param k vector of lags
#' @param theta list with components ar, ma and H specifying the ARFIMA parameters
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_arfima1 <- function(k, theta){
  if (is.list(theta))
    theta <- tsunlist(theta)
  phi0 <- numeric()
  theta0 <- numeric()
  H <- numeric()
  nm <- substring(names(theta), 1, 5)
  if ("phi" %in% nm)
    phi0 <- theta[nm == "phi"]
  if ("theta" %in% nm)
    theta0 <- theta[nm == "theta"]
  if ("H" %in% nm)
    H <- plogis(theta[nm == "H"])
  acvf <- arfima::tacvfARFIMA(phi = phi0, theta = theta0, H = H, maxlag = k)
  if (is.null(acvf))
    return(rep(NA, k))
  acf <- acvf[-1]/acvf[1]
  pacf <- acf2pacf(acf)
  (2/pi)*asin(pacf)
}


#' KPACF of Fractional Brownian Noise
#'
#' @param k vector of lags
#' @param theta parameter of process
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_fbn <- function(k, theta){
  if (is.list(theta))
    theta <- tsunlist(theta)
  theta <- plogis(theta)
  acf <- (((1:k) + 1)^{2 * theta[1]} + abs((1:k) - 1)^{2 * theta[1]} - 2 * (1:k)^{2 * theta[1]})/2
  pacf <- acf2pacf(acf)
  return(suppressWarnings((2/pi)*asin(pacf)))
}

#' KPACF of Exponential Type
#'
#' @param k vector of lags
#' @param theta parameters of exponential decay function
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_exp <- function(k, theta){
  if (is.list(theta))
    theta <- tsunlist(theta)
  arg <- pmin(theta[1] + theta[2] * (1:k), 0)
  exp(arg)
}

#' KPACF of Power Type
#'
#' @param k vector of lags
#' @param theta parameters of power decay function
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_pow <- function(k, theta){
  if (is.list(theta))
    theta <- tsunlist(theta)
  arg <- pmin(theta[1] + theta[2] * log(1:k), 0)
  exp(arg)
}

#' KPACF of Transformed Exponential Type
#'
#' @param k vector of lags
#' @param theta parameters of exponential decay function
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_exp2 <- function(k, theta){
  if (is.list(theta))
    theta <- tsunlist(theta)
  arg <- pmin(theta[1] + theta[2] * (1:k), 0)
  acf <- exp(arg)
  acf2pacf(acf)
}

#' KPACF of Transformed Power Type
#'
#' @param k vector of lags
#' @param theta parameters of power decay function
#'
#' @return vector of Kendall partial autocorrelations for each lag k
#' @export
#'
kpacf_pow2 <- function(k, theta){
  if (is.list(theta))
    theta <- tsunlist(theta)
  arg <- pmin(theta[1] + theta[2] * log(1:k), 0)
  acf <- exp(arg)
  acf2pacf(acf)
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

#' Objective Function for dvinecopula2 process
#'
#' @param theta parameters of kpacf
#' @param modelspec list specifying model
#' @param u data
#' @return
#' @keywords internal
dvinecopula2_objective <- function(theta, modelspec, u) {
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

#' Transform Kendall's tau Values to Copula Parameters
#'
#' @param family name of copula family
#' @param tau value of Kendall's tau
#'
#' @return
#' @keywords internal
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
setMethod("sim", c(x = "dvinecopula2"), function(x, n = 1000) {
  kpacf <- eval(parse(text = x@modelspec$kpacf))
  tauvals <- kpacf((n-1), x@pars)
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

#' Plot Function for dvinecopula2 Objects
#'
#' @param copula a fitted dvinecopula2 object
#' @param data the data to which copula is fitted
#' @param plotoption number giving plot choice
#' @param bw logical for black-white plot
#' @param klimit maximum lag value for plots
#'
#' @return
#' @export
plot_dvinecopula2 <- function(copula, data, plotoption, bw, klimit){
  data0 <- data
  n <- length(data)
  kpacf <- eval(parse(text = copula@modelspec$kpacf))
  tauvals <- kpacf((n-1), copula@pars)
  kmax <- max((1:(n-1))[abs(tauvals) > .Machine$double.eps])
  k <- min(max(1, kmax), copula@modelspec$maxlag, klimit)
  kplotmax <- min(k, 9)
  datavecs <- vector(mode = "list", length = kplotmax)
  tau_empirical <- rep(NA, k)
  data <- cbind(as.numeric(data[1:(n - 1)]), as.numeric(data[2:n]))
  datavecs[[1]] <- data
  tau_empirical[1] <- cor(data, method = "kendall")[1, 2]
  for (i in 1:(k-1)) {
    n <- dim(data)[1]
    model <- rvinecopulib::bicop_dist(
      family = tolower(copula@modelspec$family),
      rotation = copula@modelspec$rotation,
      parameters = rvinecopulib::ktau_to_par(
        family = copula@modelspec$family,
        tau = tauvals[i]))
    data <-
      cbind(rvinecopulib::hbicop(data[(1:(n - 1)), ], model, cond_var = 2),
            rvinecopulib::hbicop(data[(2:n), ], model, cond_var = 1))
    tau_empirical[i+1] <- cor(data, method = "kendall")[1, 2]
    if (i < kplotmax)
      datavecs[[i+1]] <- data
  }
  tau_theoretical <- tauvals[1:k]
  colchoice <- ifelse(bw, "gray50", "red")
  switch(plotoption,
         {
           plot(1:k, tau_empirical, xlab ="k", ylab = "tau",
                ylim = range(tau_empirical, tau_theoretical), type = "h")
           lines(1:k, tau_theoretical, col = colchoice)
           abline(h = 0)},
         {
           pacf0 <- pacf(qnorm(data0), lag.max = k, plot = FALSE)
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

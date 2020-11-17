#' Generalized lag plot for dvinecopula object and data
#'
#' @param data a vector or time series of data
#' @param vinemodel an object of class \linkS4class{dvinecopula} or a fitted version thereof
#' @param kmax maximum number of lags in plot
#'
#' @return an invisble vector of Kendall's tau values
#' @export
#'
#' @examples
#' mixmod <- dvinecopula(family = c("gumbel", "gauss", "joe", "clayton"), pars = list(1.5, -0.6, 1.6, 2.1))
#' data <- sim(mixmod, n = 500)
#' glag.plot(data, mixmod)
glag.plot <- function(data, vinemodel, kmax = 9) {
  if (is(vinemodel, "tscopulafit"))
    vinemodel <- vinemodel@tscopula
  if (!is(vinemodel, "dvinecopula"))
    stop("This function is for d-vine copulas")
  k <- min(length(vinemodel@modelspec), kmax)
  kend <- rep(NA, k)
  if (k >= 9)
    lc <- 3
  else
    lc <- 2
  lr <- ceiling(k / lc)
  par(
    mfrow = c(lr, lc),
    mar = c(2.1, 2.1, 1.5, 0.5),
    oma = rep(2, 4),
    pty = "s",
    cex = 0.5
  )
  n <- length(data)
  data <- cbind(data[1:(n - 1)], data[2:n])
  kend[1] <- cor(data, method = "kendall")[1, 2]
  plot(
    data,
    main = paste("Lag 1 : tau = ", round(kend[1], 2), sep = ""),
    asp = 1,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "",
    ylab = ""
  )
  for (i in 1:(k - 1)) {
    n <- dim(data)[1]
    model <- rvinecopulib::bicop_dist(
      family = tolower(vinemodel@modelspec[[i]]$family),
      rotation = vinemodel@modelspec[[i]]$rotation,
      parameters = vinemodel@pars[[i]][1:vinemodel@modelspec[[i]]$npars]
    )
    data <-
      cbind(hbicop(data[(1:(n - 1)), ], model, cond_var = 2),
            hbicop(data[(2:n), ], model, cond_var = 1))
    kend[i + 1] <- cor(data, method = "kendall")[1, 2]
    plot(
      data,
      main = paste("Lag ", i + 1 , " : tau = ", round(kend[i + 1], 2), sep =
                     ""),
      asp = 1,
      xlim = c(0, 1),
      ylim = c(0, 1),
      xlab = "",
      ylab = ""
    )
  }
  par(mfrow = c(1, 1))
  invisible(kend)
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

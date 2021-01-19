#' Marginal Model for Time Series
#'
#' Class of objects for marginal models for stationary time series. The
#' object is given a name and there must exist functions pname, qname,
#' dname and rname. For example, the object could be named norm and
#' make use of \code{\link[stats]{pnorm}}, \code{\link[stats]{qnorm}},
#' \code{\link[stats]{dnorm}} and \code{\link[stats]{rnorm}}.
#' As well as the parameters of the distribution, dname must have the
#' logical argument log specifying whether log density should be computed.
#'
#' @slot name name of the marginal model class.
#' @slot pars a numeric vector containing the named parameters of the distribution
#' which are passed as arguments to pname, qname, dname and rname.
#'
#' @export
#'
#' @examples
#' new("margin", name = "norm", pars = c(mu = 0, sigma = 1))
setClass("margin", slots = list(
  name = "character",
  pars = "numeric"
))

#' Coef Method for margin Class
#'
#' @param object an object of class \linkS4class{margin}.
#'
#' @return Vector of coefficients of marginal distribution.
#' @export
#'
#'
setMethod("coef", "margin", function(object) {
  if (is(object, "marginfit")) {
    object <- object@margin
  }
  object@pars
})

#' Constructor Function for margin
#'
#' @param name character string giving name of distribution
#' @param pars parameters of the distribution
#'
#' @return An object of class \linkS4class{margin}.
#' @export
#'
#' @examples
#' margin("sst")
margin <- function(name, pars = NULL) {
  pfunc <- eval(parse(text = paste("p", name, sep = "")))
  defaults <- unlist(formals(pfunc)[-1])
  if (name == "norm") {
    defaults <- defaults[1:2]
  }
  if (!is.null(pars)) {
    parnames <- names(pars)
    parset <- parnames[parnames %in% names(defaults)]
    if (length(parset) == 0) {
      stop("Unrecognized parameter names")
    }
    defaults[parset] <- pars[parset]
  }
  new("margin", name = name, pars = defaults)
}

#' Compute CDF of Marginal Model
#'
#' Compute the cumulative distribution function of the marginal model.
#'
#' @param x an object of class \linkS4class{margin}.
#' @param q vector of values at which CDF should be computed.
#'
#' @return A vector of values for the CDF.
#' @export
#'
#' @examples
#' margmod <- margin("norm", pars = c(mean = 0, sd = 1))
#' pmarg(margmod, c(-2, 0, 2))
pmarg <- function(x, q) {
  if (is(x, "marginfit")) {
    x <- x@margin
  }
  func <- eval(parse(text = paste("p", x@name, sep = "")))
  do.call(func, append(x@pars, list(q = q)))
}

#' Compute Quantiles of Marginal Model
#'
#' Compute the quantile function of the marginal model.
#'
#' @param x an object of class \linkS4class{margin}.
#' @param p vector of probabilities for which quantiles should be computed.
#'
#' @return A vector of values for the quantile function.
#' @export
#'
#' @examples
#' margmod <- margin("norm", pars = c(mean = 0, sd = 1))
#' qmarg(margmod, c(0.05, 0.5, 0.95))
qmarg <- function(x, p) {
  if (is(x, "marginfit")) {
    x <- x@margin
  }
  func <- eval(parse(text = paste("q", x@name, sep = "")))
  do.call(func, append(x@pars, list(p = p)))
}

#' Compute Density of Marginal Model
#'
#' Compute the density function of the marginal model.
#'
#' @param x an object of class \linkS4class{margin}.
#' @param y vector of values for which density should be computed.
#' @param log logical variable specifying whether log density should be
#' returned.
#'
#' @return A vector of values for the density.
#' @export
#'
#' @examples
#' margmod <- margin("norm", pars = c(mean = 0, sd = 1))
#' dmarg(margmod, c(-2, 0, 2), log = TRUE)
dmarg <- function(x, y, log = FALSE) {
  if (is(x, "marginfit")) {
    x <- x@margin
  }
  func <- eval(parse(text = paste("d", x@name, sep = "")))
  do.call(func, append(x@pars, list(x = y, log = log)))
}

#' Density of Laplace Distribution
#'
#' @param x argument of density
#' @param mu location parameter
#' @param scale scale parameter
#' @param log flag for log density
#'
#' @return A vector of values for the density.
#' @keywords internal
#'
dlaplace <- function(x, mu = 0, scale = 1, log = FALSE){
  dsdoubleweibull(x, mu = mu, shape = 1, scale = scale, gamma = 1, log = log)
}

#' Distribution Function of Laplace Distribution
#'
#' @param q argument of distribution function
#' @param mu location parameter
#' @param scale scale parameter
#'
#' @return A vector of values for the distribution function.
#' @keywords internal
#'
plaplace <- function(q, mu = 0, scale = 1){
  psdoubleweibull(q, mu = mu, shape = 1, scale = scale, gamma = 1)
}

#' Quantile Function of Laplace Distribution
#'
#' @param p probability
#' @param mu location parameter
#' @param scale scale parameter
#'
#' @return A vector of values for the quantile function.
#' @keywords internal
#'
qlaplace <- function(p, mu = 0, scale = 1){
  qsdoubleweibull(p, mu = mu, shape = 1, scale = scale, gamma = 1)
}

#' Random Number Generation for Laplace Distribution
#'
#' @param n size
#' @param mu location parameter
#' @param scale scale parameter
#'
#' @return A vector of simulated values of length n.
#' @keywords internal
#'
rlaplace <- function(n, mu = 0, scale = 1){
  qlaplace(runif(n), mu, scale)
}

#' Density of Skew Laplace Distribution
#'
#' @param x argument of density
#' @param mu location parameter
#' @param scale scale parameter
#' @param skewness parameter
#' @param log flag for log density
#'
#' @return A vector of values for the density.
#' @keywords internal
#'
dslaplace <- function(x, mu = 0, scale = 1, gamma = 1, log = FALSE){
  dsdoubleweibull(x, mu = mu, shape = 1, scale = scale, gamma = gamma, log = log)
}

#' Distribution Function of Skew Laplace Distribution
#'
#' @param q argument of distribution function
#' @param mu location parameter
#' @param scale scale parameter
#' @param gamma skewness parameter
#'
#' @return A vector of values for the CDF.
#' @keywords internal
#'
pslaplace <- function(q, mu = 0, scale = 1, gamma = 1){
  psdoubleweibull(q, mu = mu, shape = 1, scale = scale, gamma = gamma)
}

#' Quantile Function of Skew Laplace Distribution
#'
#' @param p probability
#' @param mu location parameter
#' @param scale scale parameter
#' @param gamma skewness parameter
#'
#' @return A vector of values for the quantile function.
#' @keywords internal
#'
qslaplace <- function(p, mu = 0, scale = 1, gamma = 1){
  qsdoubleweibull(p, mu = mu, shape = 1, scale = scale, gamma = gamma)
}

#' Random Number Generation for Skew Laplace Distribution
#'
#' @param n size
#' @param mu location parameter
#' @param scale scale parameter
#' @param gamma skewness paramater
#'
#' @return A vector of simulated values of length n.
#' @keywords internal
#'
rslaplace <- function(n, mu = 0, scale = 1, gamma = 1){
  qslaplace(runif(n), mu, scale, gamma)
}

#' Density of Double Weibull Distribution
#'
#' @param x argument of density
#' @param mu location parameter
#' @param shape shape parameter
#' @param scale scale parameter
#' @param log flag for log density
#'
#' @return A vector of values for the density.
#' @keywords internal
#'
ddoubleweibull <- function(x, mu = 0, shape = 1, scale = 1, log = FALSE){
  dsdoubleweibull(x, mu = mu, shape = shape, scale = scale, gamma = 1, log = log)
}

#' Distribution Function of Double Weibull Distribution
#'
#' @param q argument of distribution function
#' @param mu location parameter
#' @param shape shape paramater
#' @param scale scale parameter
#'
#' @return A vector of values for the CDF.
#' @keywords internal
#'
pdoubleweibull <- function(q, mu = 0, shape = 1, scale = 1){
  psdoubleweibull(q, mu = mu, shape = shape, scale = scale, gamma = 1)
}

#' Quantile Function of Double Weibull Distribution
#'
#' @param p probability
#' @param mu location parameter
#' @param shape shape parameter
#' @param scale scale parameter
#'
#' @return A vector of values for the quantile function.
#' @keywords internal
#'
qdoubleweibull <- function(p, mu = 0, shape = 1, scale = 1){
  qsdoubleweibull(p, mu = mu, shape = shape, scale = scale, gamma = 1)
}

#' Random Number Generation for Double Weibull Distribution
#'
#' @param n size
#' @param mu location parameter
#' @param shape shape parameter
#' @param scale scale parameter
#'
#' @return A vector of simulated values of length n.
#' @keywords internal
#'
rdoubleweibull <- function(n, mu = 0, shape = 1, scale = 1){
  qdoubleweibull(runif(n), mu = mu, shape = shape, scale = scale)
}

#' Density of Skew Double Weibull Distribution
#'
#' @param x argument of density
#' @param mu location parameter
#' @param shape shape parameter
#' @param scale scale parameter
#' @param gamma skewness parameter
#' @param log flag for log density
#'
#' @return A vector of values for the density.
#' @keywords internal
#'
dsdoubleweibull <- function(x, mu = 0, shape = 1, scale = 1, gamma = 1, log = FALSE)
{
  if ((scale <= 0) | (shape <= 0)) {
    return(NA)
  }
  arg <- rep(NA, length(x))
  y <- (x - mu)/scale
  arg[y < 0] <- gamma * abs(y[y < 0])
  arg[y >= 0] <- y[y >=0]/gamma
  tmp <- (shape - 1)*log(arg) - arg^shape + log(shape) -log(scale) - log(gamma + 1/gamma)
  if (!log)
    tmp <- exp(tmp)
  tmp
}

#' Distribution Function of Skew Double Weibull Distribution
#'
#' @param q argument of distribution function
#' @param mu location parameter
#' @param shape shape paramater
#' @param scale scale parameter
#' @param gamma skewness parameter
#'
#' @return A vector of values for the CDF.
#' @keywords internal
#'
psdoubleweibull <- function(q, mu = 0, shape = 1, scale = 1, gamma = 1)
{
  arg <- rep(NA, length(q))
  y <- (q - mu)/scale
  arg[y < 0] <- gamma * abs(y[y < 0])
  arg[y >= 0] <- y[y >=0]/gamma
  cumy <- exp(-arg^shape)/(1+gamma^2)
  cumy[y > 0] <- (1-(gamma^2)*cumy[y > 0])
  cumy
}

#' Quantile Function of Skew Double Weibull Distribution
#'
#' @param p probability
#' @param mu location parameter
#' @param shape shape parameter
#' @param scale scale parameter
#' @param gamma skewness parameter
#'
#' @return A vector of values for the quantile function.
#' @keywords internal
#'
qsdoubleweibull <- function(p, mu = 0, shape = 1, scale = 1, gamma = 1)
{
  tmp <- rep(NA, length(p))
  wt <- (1+gamma^2)
  lower <- (p <= 1/wt)
  tmp[lower] <- mu - scale * ((-log(wt*p[lower]))^(1/shape))/gamma
  tmp[!lower] <- mu + scale * ((-log(wt*(1-p[!lower])/(wt-1)))^(1/shape))*gamma
  tmp
}

#' Random Number Generation for Skew Double Weibull Distribution
#'
#' @param n size
#' @param mu location parameter
#' @param shape shape parameter
#' @param scale scale parameter
#' @param gamma skewness parameter
#'
#' @return A vector of values of length n.
#' @keywords internal
#'
rsdoubleweibull <- function(n, mu = 0, shape = 1, scale = 1, gamma = 1)
{
  qsdoubleweibull(runif(n), mu, shape, scale, gamma)
}

#' Distribution Function of Student t
#'
#' @param q quantile
#' @param df degrees of freedom
#' @param mu location parameter
#' @param sigma scale parameter
#'
#' @return A vector of values for the CDF.
#' @keywords internal
#'
pst <- function(q, df = 10, mu = 0, sigma = 1) {
  pt((q - mu) / sigma, df)
}

#' Quantile Function of Student t
#'
#' @param p probability
#' @param df degrees of freedom
#' @param mu location parameter
#' @param sigma scale parameter
#'
#' @return A vector of values for the quantile function.
#' @keywords internal
#'
qst <- function(p, df, mu, sigma) {
  qt(p, df) * sigma + mu
}

#' Density of Student t
#'
#' @param x variable
#' @param df degrees of freedom
#' @param mu location parameter
#' @param sigma scale parameter
#' @param log
#'
#' @return A vector of values for the density.
#' @keywords internal
#'
dst <- function(x, df, mu, sigma, log = FALSE) {
  if ((sigma < 0) | (df < 0)) {
    return(NA)
  }
  dens <- dt((x - mu) / sigma, df = df, log = log)
  if (log) {
    return(dens - log(sigma))
  } else {
    return(dens / sigma)
  }
}

#' Random Number Generation for Student t
#'
#' @param n size
#' @param df degrees of freedom
#' @param mu location parameter
#' @param sigma scale parameter
#'
#' @return A vector of values of length n.
#' @keywords internal
#'
rst <- function(n, df, mu, sigma) {
  rt(n, df) * sigma + mu
}

#' Distribution Function of Skewed Student t
#'
#' @param q quantile
#' @param df degrees of freedom
#' @param mu location parameter
#' @param sigma scale parameter
#' @param gamma skewness parameter
#'
#' @return A vector of values for the CDF.
#' @keywords internal
#'
psst <- function(q, df = 10, gamma = 1, mu = 0, sigma = 1) {
  result <- rep(NA, length(q))
  x <- (q - mu) / sigma
  result[x < 0] <- 2 / (gamma^2 + 1) * pt(gamma * x[x < 0], df)
  result[x >= 0] <- 1 / (gamma^2 + 1) + 2 / (1 + (1 / gamma^2)) * (pt(x[x >= 0] / gamma, df) -
    1 / 2)
  result
}

#' Quantile Function of Skewed Student t
#'
#' @param p probability
#' @param df degrees of freedom
#' @param mu location parameter
#' @param sigma scale parameter
#' @param gamma skewness parameter
#'
#' @return A vector of values for the quantile function.
#' @keywords internal
#'
qsst <- function(p, df, gamma, mu, sigma) {
  result <- rep(NA, length(p))
  probzero <- 1 / (gamma^2 + 1)
  result[p < probzero] <- 1 / gamma * qt(((gamma^2 + 1) * p[p < probzero]) / 2, df)
  result[p >= probzero] <- gamma * qt((1 + 1 / gamma^2) / 2 * (p[p >= probzero] - probzero) +
    1 / 2, df)
  result * sigma + mu
}

#' Density of Skewed Student t
#'
#' @param x variable
#' @param df degrees of freedom
#' @param mu location parameter
#' @param sigma scale parameter
#' @param gamma skewness parameter
#' @param log
#'
#' @return A vector of values for the density.
#' @keywords internal
#'
dsst <- function(x, df, gamma, mu, sigma, log = FALSE) {
  if ((sigma < 0) | (df < 0)) {
    return(NA)
  }
  result <- rep(NA, length(x))
  x <- (x - mu) / sigma
  result[x < 0] <- dt(gamma * x[x < 0], df, log = log)
  result[x >= 0] <- dt(x[x >= 0] / gamma, df, log = log)
  if (log) {
    return(result + log(2 / (gamma + 1 / gamma)) - log(sigma))
  } else {
    return(result * (2 / (gamma + 1 / gamma)) / sigma)
  }
}

#' Random Number Generation for Skewed Student t
#'
#' @param n size
#' @param df degrees of freedom
#' @param mu location parameter
#' @param sigma scale parameter
#' @param gamma skewness parameter
#'
#' @return A vector of simulated values of length n.
#' @keywords internal
#'
rsst <- function(n, df, gamma, mu, sigma) {
  p <- runif(n)
  result <- rep(NA, n)
  probzero <- 1 / (gamma^2 + 1)
  result[p < probzero] <- 1 / gamma * qt(((gamma^2 + 1) * p[p < probzero]) / 2, df)
  result[p >= probzero] <- gamma * qt((1 + 1 / gamma^2) / 2 * (p[p >= probzero] - probzero) +
    1 / 2, df)
  result * sigma + mu
}

#' Simulation Method for margin Class
#'
#' @param x an object of class \linkS4class{margin}.
#' @param n length of realization.
#'
#' @return A vector of simulated values of length n.
#' @export
#'
#' @examples
#' margmod <- margin("norm", pars = c(mean = 0, sd = 1))
#' sim(margmod, n = 500)
setMethod("sim", c(x = "margin"), function(x, n = 1000) {
  func <- eval(parse(text = paste("r", x@name, sep = "")))
  do.call(func, append(x@pars, list(n = n)))
})

#' Fitted Marginal Model for Time Series
#'
#' @slot margin an object of class \linkS4class{margin}.
#' @slot data numeric vector or time series of data.
#' @slot fit a list containing details of the maximum likelihood fit.
#'
#' @export
#'
setClass("marginfit",
  contains = "margin",
  slots = list(
    margin = "margin",
    data = "ANY",
    fit = "list"
  )
)

#' Fit Method for margin Class
#'
#' @param x an object of class \linkS4class{margin}.
#' @param y a vector or time series of data.
#' @param tsoptions list of optional arguments:
#' hessian is logical variable specifying whether Hessian matrix should be returned;
#' start is vector od named starting values
#' @param control list of control parameters to be passed to the
#' \code{\link[stats]{optim}} function.
#'
#' @return An object of class \linkS4class{marginfit}.
#' @export
#'
#' @examples
#' margmod <- margin("norm", pars = c(mean = 0, sd = 1))
#' data <- sim(margmod, n = 500)
#' fit(margmod, data)
setMethod("fit", c(x = "margin", y = "ANY"), function(x, y,
                                                      tsoptions = list(),
                                                      control = list(maxit = 1000)) {
  defaults <- list(hessian = FALSE, method = "Nelder-Mead")
  tsoptions <- setoptions(tsoptions, defaults)
  dens <- eval(parse(text = paste("d", x@name, sep = "")))
  objective <- function(theta, dens, y) {
    dx <- do.call(dens, append(theta, list(x = y, log = TRUE)))
    -sum(dx)
  }
  fit <- optim(x@pars,
    fn = objective,
    dens = dens,
    y = y,
    method = tsoptions$method,
    hessian = tsoptions$hessian,
    control = control
  )
  x@pars <- fit$par
  new("marginfit", margin = x, data = y, fit = fit)
})

#' Show Method for margin Class
#'
#' @param object an object of class \linkS4class{margin}.
#'
#' @return A summary of an object of class \linkS4class{margin}.
#' @export
#'
#'
setMethod("show", "margin", function(object) {
  if (!(is(object, "marginfit"))){
    cat("name: ", object@name, "\n", sep = "")
    cat("parameters: \n")
    print(coef(object))
  }
  if (is(object, "marginfit")) {
    cat("name: ", object@margin@name, "\n", sep = "")
    ests <- object@fit$par
    if (is.element("hessian", names(object@fit))) {
      ses <- safe_ses(object@fit$hessian)
      ests <- rbind(ests, ses)
      dimnames(ests)[[1]] <- c("par", "se")
    }
    cat("estimates:\n")
    print(ests)
    cat("convergence status: ", object@fit$convergence, ", log-likelihood: ", -object@fit$value,
      "\n",
      sep = ""
    )
  }
})

#' logLik Method for marginfit Class
#'
#' @param object an object of class \linkS4class{marginfit}.
#'
#' @return An object of class logLik.
#' @export
#'
setMethod("logLik", "marginfit", function(object) {
  ll <- -object@fit$value
  attr(ll, "nobs") <- length(object@data)
  attr(ll, "df") <- length(object@fit$par)
  class(ll) <- "logLik"
  ll
})

#' Plot Method for marginfit Class
#'
#' @param x an object of class \linkS4class{marginfit}.
#' @param plotoption number of plot required.
#' @param bw logical variable specifying whether black-white options should be chosen.
#'
#' @export
#'
setMethod("plot", c(x = "marginfit", y = "missing"),
          function(x, plotoption = 1L, bw = FALSE) {
            colchoice <- ifelse(bw, "gray50", "red")
            switch(plotoption,
             {
              qus <- qmarg(x@margin, ppoints(x@data))
              plot(qus,
                   sort(as.numeric(x@data)),
                   xlab = "Theoretical",
                   ylab = "data"
              )
              abline(0, 1, col = colchoice)
            },
            zoo::plot.zoo(x@data, xlab = "", ylab = "data", type = "h"),
            stop("Not a plot option")
            )
            })



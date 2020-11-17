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
#' @slot fixed logical vector specifying any parameters that are held fixed
#'
#' @return
#' @export
#'
#' @examples
#' new("margin", name = "norm", pars = list(mu = 0, sigma = 1))
setClass("margin", slots = list(
  name = "character",
  pars = "numeric",
  fixed = "logical"
))

#' Coef Method for margin Class
#'
#' @param object an object of class \linkS4class{margin}.
#'
#' @return
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
#' @param fixed fixed parameters, for estimation purposes
#'
#' @return an object of class \linkS4class{margin}
#' @export
#'
#' @examples
#' margin("sst")
#' margin("weibull2")
margin <- function(name, pars = NULL, fixed = FALSE) {
  pfunc <- eval(parse(text = paste("p", name, sep = "")))
  defaults <- unlist(formals(pfunc)[-1])
  if (name == "norm") {
    defaults <- defaults[1:2]
  }
  if (sum(fixed) == 0) {
    fixed <- rep(FALSE, length(defaults))
  }
  if (!is.null(pars)) {
    parnames <- names(pars)
    parset <- parnames[parnames %in% names(defaults)]
    if (length(parset) == 0) {
      stop("Unrecognized parameter names")
    }
    defaults[parset] <- pars[parset]
  }
  if (length(fixed) != length(defaults)) {
    stop("Incorrect fixed specification")
  }
  new("margin", name = name, pars = defaults, fixed = fixed)
}

#' Compute CDF of Marginal Model
#'
#' Compute the cumulative distribution function of the marginal model.
#'
#' @param x an object of class \linkS4class{margin}.
#' @param q vector of values at which CDF should be computed.
#'
#' @return
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
#' @return
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
#' @return
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

#' Empirical Distribution Function
#'
#' Compute a version of the empirical distribution function that is
#' standardized by (n+1) to lie strictly in (0,1).
#'
#' @param q vector of values at which to evaluate the empirical distribution
#' function.
#' @param data vector of data values specifying the empirical distribution
#' function
#'
#' @return
#' @export
#'
#' @examples
#' data <- rnorm(100)
#' pedf(c(-2, -1, 0, 1, 2), data)
pedf <- function(q, data) {
  n <- length(data)
  (ecdf(data)(q) * n + 0.5) / (n + 1)
}

#' Quantiles of Empirical Distribution Function
#'
#' @param p vector of probabilities at which to evaluate the quantiles
#' of the empirical distribution function.
#' @param data vector of data values specifying the empirical distribution
#' function
#'
#' @return
#' @export
#'
#' @examples
#' data <- rnorm(100)
#' qedf(c(0.25, 0.5, 0.75), data)
qedf <- function(p, data) {
  quantile(data, p, type = 6)
}

#' Density of Laplace Distribution
#'
#' @param x argument of density
#' @param mu location parameter
#' @param scale scale parameter
#' @param log flag for log density
#'
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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
#' @return
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

#' Distribution Function of NIG
#'
#' @param q quantile
#' @param alphabar
#' @param gamma
#' @param mu
#' @param sigma
#'
#' @return
#' @keywords internal
#'
pNIG <- function(q, alphabar = 1, gamma = 0, mu = 0, sigma = 1) {
  obj <- ghyp::NIG(alpha.bar = alphabar, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::pghyp(q, object = obj)
}

#' Quantile Function of NIG
#'
#' @param p probability
#' @param alphabar
#' @param gamma
#' @param mu
#' @param sigma
#'
#' @return
#' @keywords internal
#'
qNIG <- function(p, alphabar, gamma, mu, sigma) {
  obj <- ghyp::NIG(alpha.bar = alphabar, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::qghyp(p, object = obj)
}

#' Density of NIG
#'
#' @param x variable
#' @param alphabar
#' @param gamma
#' @param mu
#' @param sigma
#' @param log
#'
#' @return
#' @keywords internal
#'
dNIG <- function(x, alphabar, gamma, mu, sigma, log = FALSE) {
  if ((alphabar < 0) | (sigma < 0)) {
    return(NA)
  }
  obj <- ghyp::NIG(alpha.bar = alphabar, mu = mu, sigma = sigma, gamma = gamma)
  out <- ghyp::dghyp(x, object = obj, logvalue = TRUE)
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}

#' Random Number Generation for NIG
#'
#' @param n size
#' @param alphabar
#' @param gamma
#' @param mu
#' @param sigma
#'
#' @return
#' @keywords internal
#'
rNIG <- function(n, alphabar, gamma, mu, sigma) {
  obj <- ghyp::NIG(alpha.bar = alphabar, mu = mu, sigma = sigma, gamma = gamma)
  ghyp::rghyp(n, object = obj)
}

#' Distribution Function of weibull2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pweibull2 <- function(q, sh1 = 1, sc1 = 1, sh2 = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * pweibull(q[q > 0], shape = sh1, scale = sc1)
  result[q <= 0] <- (1 - p0) * pweibull(-q[q <= 0], shape = sh2, scale = sc2, lower.tail = F)
  return(result)
}

#' Quantile Function of weibull2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qweibull2 <- function(p, sh1 = 1, sc1 = 1, sh2 = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- qweibull((p[p > (1 - p0)] - (1 - p0)) / p0, shape = sh1, scale = sc1)
  result[p <= (1 - p0)] <- -qweibull(1 - p[p <= (1 - p0)] / (1 - p0), shape = sh2, scale = sc2)
  return(result)
}

#' Density of weibull2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dweibull2 <- function(x, sh1 = 1, sc1 = 1, sh2 = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * dweibull(x[x > 0], shape = sh1, scale = sc1)
  result[x <= 0] <- (1 - p0) * dweibull(-x[x <= 0], shape = sh2, scale = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for weibull2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rweibull2 <- function(n, sh1 = 1, sc1 = 1, sh2 = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- rweibull(sum(mix), shape = sh1, scale = sc1)
  result[mix == 0] <- -rweibull(sum(1 - mix), shape = sh2, scale = sc2)
  return(result)
}

#' Distribution Function of burr2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pburr2 <- function(q, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * actuar::pburr(q[q > 0], shape1 = sh1, shape2 = sh1b, scale = sc1)
  result[q <= 0] <- (1 - p0) * actuar::pburr(-q[q <= 0], shape1 = sh2, shape2 = sh2b, scale = sc2, lower.tail = F)
  return(result)
}

#' Quantile Function of burr2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qburr2 <- function(p, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- actuar::qburr((p[p > (1 - p0)] - (1 - p0)) / p0, shape1 = sh1, shape2 = sh1b, scale = sc1)
  result[p <= (1 - p0)] <- -actuar::qburr(1 - p[p <= (1 - p0)] / (1 - p0), shape1 = sh2, shape2 = sh2b, scale = sc2)
  return(result)
}

#' Density of burr2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dburr2 <- function(x, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * actuar::dburr(x[x > 0], shape1 = sh1, shape2 = sh1b, scale = sc1)
  result[x <= 0] <- (1 - p0) * actuar::dburr(-x[x <= 0], shape1 = sh2, shape2 = sh2b, scale = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for burr2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rburr2 <- function(n, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- actuar::rburr(sum(mix), shape1 = sh1, shape2 = sh1b, scale = sc1)
  result[mix == 0] <- -actuar::rburr(sum(1 - mix), shape1 = sh2, shape2 = sh2b, scale = sc2)
  return(result)
}

#' Distribution Function of lgamma2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shapelog parameter of right branch
#' @param r1 ratelog parameter of right branch
#' @param sh2 shapelog parameter of left branch
#' @param r2 ratelog parameter of left branch
#'
#' @return
#' @keywords internal
#'
plgamma2 <- function(q, sh1 = 1, r1 = 2, sh2 = 1, r2 = 2, p0 = 0) {
  sh1 <- abs(sh1)
  r1 <- abs(r1)
  sh2 <- abs(sh2)
  r2 <- abs(r2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * actuar::plgamma(q[q > 0] + 1, shapelog = sh1, ratelog = r1)
  result[q <= 0] <- (1 - p0) * actuar::plgamma(-q[q <= 0] + 1, shapelog = sh2, ratelog = r2, lower.tail = F)
  return(result)
}

#' Quantile Function of lgamma2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shapelog parameter of right branch
#' @param r1 ratelog parameter of right branch
#' @param sh2 shapelog parameter of left branch
#' @param r2 ratelog parameter of left branch
#'
#' @return
#' @keywords internal
#'
qlgamma2 <- function(p, sh1 = 1, r1 = 2, sh2 = 1, r2 = 2, p0 = 0) {
  sh1 <- abs(sh1)
  r1 <- abs(r1)
  sh2 <- abs(sh2)
  r2 <- abs(r2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- actuar::qlgamma((p[p > (1 - p0)] - (1 - p0)) / p0, shapelog = sh1, ratelog = r1) - 1
  result[p <= (1 - p0)] <- -actuar::qlgamma(1 - p[p <= (1 - p0)] / (1 - p0), shapelog = sh2, ratelog = r2) + 1
  return(result)
}

#' Density of lgamma2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shapelog parameter of right branch
#' @param r1 ratelog parameter of right branch
#' @param sh2 shapelog parameter of left branch
#' @param r2 ratelog parameter of left branch
#'
#' @return
#' @keywords internal
#'
dlgamma2 <- function(x, sh1 = 1, r1 = 2, sh2 = 1, r2 = 2, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  r1 <- abs(r1)
  sh2 <- abs(sh2)
  r2 <- abs(r2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * actuar::dlgamma(x[x > 0] + 1, shapelog = sh1, ratelog = r1)
  result[x <= 0] <- (1 - p0) * actuar::dlgamma(-x[x <= 0] + 1, shapelog = sh2, ratelog = r2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for lgamma2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shapelog parameter of right branch
#' @param r1 ratelog parameter of right branch
#' @param sh2 shapelog parameter of left branch
#' @param r2 ratelog parameter of left branch
#'
#' @return
#' @keywords internal
#'
rlgamma2 <- function(n, sh1 = 1, r1 = 2, sh2 = 1, r2 = 2, p0 = 0) {
  sh1 <- abs(sh1)
  r1 <- abs(r1)
  sh2 <- abs(sh2)
  r2 <- abs(r2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- actuar::rlgamma(sum(mix), shapelog = sh1, ratelog = r1)
  result[mix == 0] <- -actuar::rlgamma(sum(1 - mix), shapelog = sh2, ratelog = r2)
  return(result)
}

#' Distribution Function of gweibull2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pgweibull2 <- function(q, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * rmutil::pgweibull(q[q > 0], s = sh1, m = sh1b, f = sc1)
  result[q <= 0] <- (1 - p0) * (1 - rmutil::pgweibull(-q[q <= 0], s = sh2, m = sh2b, f = sc2))
  return(result)
}

#' Quantile Function of gweibull2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qgweibull2 <- function(p, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- rmutil::qgweibull((p[p > (1 - p0)] - (1 - p0)) / p0, s = sh1, m = sh1b, f = sc1)
  result[p <= (1 - p0)] <- -rmutil::qgweibull(1 - p[p <= (1 - p0)] / (1 - p0), s = sh2, m = sh2b, f = sc2)
  return(result)
}

#' Density of gweibull2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dgweibull2 <- function(x, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * rmutil::dgweibull(x[x > 0], s = sh1, m = sh1b, f = sc1)
  result[x <= 0] <- (1 - p0) * rmutil::dgweibull(-x[x <= 0], s = sh2, m = sh2b, f = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for gweibull2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rgweibull2 <- function(n, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- rmutil::rgweibull(sum(mix), s = sh1, m = sh1b, f = sc1)
  result[mix == 0] <- -rmutil::rgweibull(sum(1 - mix), s = sh2, m = sh2b, f = sc2)
  return(result)
}

#' Distribution Function of ggamma2
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pggamma2 <- function(q, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * rmutil::pggamma(q[q > 0], s = sh1, m = sh1b, f = sc1)
  result[q <= 0] <- (1 - p0) * (1 - rmutil::pggamma(-q[q <= 0], s = sh2, m = sh2b, f = sc2))
  return(result)
}

#' Quantile Function of ggamma2
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qggamma2 <- function(p, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- rmutil::qggamma((p[p > (1 - p0)] - (1 - p0)) / p0, s = sh1, m = sh1b, f = sc1)
  result[p <= (1 - p0)] <- -rmutil::qggamma(1 - p[p <= (1 - p0)] / (1 - p0), s = sh2, m = sh2b, f = sc2)
  return(result)
}

#' Density of ggamma2
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dggamma2 <- function(x, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * rmutil::dggamma(x[x > 0], s = sh1, m = sh1b, f = sc1)
  result[x <= 0] <- (1 - p0) * rmutil::dggamma(-x[x <= 0], s = sh2, m = sh2b, f = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for ggamma2
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rggamma2 <- function(n, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- rmutil::rggamma(sum(mix), s = sh1, m = sh1b, f = sc1)
  result[mix == 0] <- -rmutil::rggamma(sum(1 - mix), s = sh2, m = sh2b, f = sc2)
  return(result)
}

#' Distribution Function of burrgamma
#'
#' @param q quantile
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
pburrgamma <- function(q, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- q
  result[q > 0] <- (1 - p0) + p0 * rmutil::pggamma(q[q > 0], s = sh1, m = sh1b, f = sc1)
  result[q <= 0] <- (1 - p0) * actuar::pburr(-q[q <= 0], shape1 = sh2, shape2 = sh2b, scale = sc2, lower.tail = F)
  return(result)
}

#' Quantile Function of burrgamma
#'
#' @param p probability
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
qburrgamma <- function(p, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- p
  result[p > (1 - p0)] <- rmutil::qggamma((p[p > (1 - p0)] - (1 - p0)) / p0, s = sh1, m = sh1b, f = sc1)
  result[p <= (1 - p0)] <- -actuar::qburr(1 - p[p <= (1 - p0)] / (1 - p0), shape1 = sh2, shape2 = sh2b, scale = sc2)
  return(result)
}

#' Density of burrgamma
#'
#' @param x variable
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
dburrgamma <- function(x, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0, log = F) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  result <- x
  result[x > 0] <- p0 * rmutil::dggamma(x[x > 0], s = sh1, m = sh1b, f = sc1)
  result[x <= 0] <- (1 - p0) * actuar::dburr(-x[x <= 0], shape1 = sh2, shape2 = sh2b, scale = sc2)
  if (log == T) {
    result <- log(result)
  }
  return(result)
}

#' Random number generating for burrgamma
#'
#' @param n size
#' @param p0 mixing parameter
#' @param sh1 shape parameter of right branch
#' @param sh2 shape parameter of left branch
#' @param sh1b second shape parameter of right branch
#' @param sh2b second shape parameter of left branch
#' @param sc1 scale parameter of right branch
#' @param sc2 scale parameter of left branch
#'
#' @return
#' @keywords internal
#'
rburrgamma <- function(n, sh1 = 1, sh1b = 1, sc1 = 1, sh2 = 1, sh2b = 1, sc2 = 1, p0 = 0) {
  sh1 <- abs(sh1)
  sh1b <- abs(sh1b)
  sc1 <- abs(sc1)
  sh2 <- abs(sh2)
  sh2b <- abs(sh2b)
  sc2 <- abs(sc2)
  p0 <- plogis(p0)
  mix <- rbinom(n, 1, p0)
  result <- mix
  result[mix == 1] <- rmutil::rggamma(sum(mix), s = sh1, m = sh1b, f = sc1)
  result[mix == 0] <- -actuar::rburr(sum(1 - mix), shape1 = sh2, shape2 = sh2b, scale = sc2)
  return(result)
}

#' Simulation Method for margin Class
#'
#' @param x an object of class \linkS4class{margin}.
#' @param n length of realization.
#'
#' @return
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
#' @return
#' @export
#'
#' @examples
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
  objective <- function(theta_free, theta_fix, dens, y) {
    theta <- c(theta_free, theta_fix)
    dx <- do.call(dens, append(theta, list(x = y, log = TRUE)))
    -sum(dx)
  }
  fit <- optim(x@pars[!(x@fixed)],
    fn = objective,
    theta_fix = x@pars[x@fixed],
    dens = dens,
    y = y,
    method = tsoptions$method,
    hessian = tsoptions$hessian,
    control = control
  )
  x@pars[!(x@fixed)] <- fit$par
  new("marginfit", margin = x, data = y, fit = fit)
})

#' Show Method for margin Class
#'
#' @param object an object of class \linkS4class{margin}.
#'
#' @return
#' @export
#'
#'
setMethod("show", "margin", function(object) {
  if (!is(object, "marginfit")) {
    cat("name: ", object@name, "\n", sep = "")
    cat("parameters: \n")
    print(coef(object))
  }
  if (is(object, "marginfit")) {
    cat("name: ", object@margin@name, "\n", sep = "")
    if (sum(object@margin@fixed) != 0) {
      cat("parameters: \n")
      print(coef(object@margin))
    }
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

#' Calculate Standardized Ranks of Data
#'
#' @param x a vector or time series of data.
#'
#' @return A vector or time series of standardized ranks in the interval (0,1)
#' @export
#'
#' @examples
#' strank(rnorm(100))
strank <- function(x) {
  U <- rank(x, na.last = "keep", ties.method = "random") / (length(x) + 1)
  if (inherits(x, "zoo")) {
    attributes(U) <- attributes(x)
  }
  U
}

#' logLik Method for marginfit Class
#'
#' @param object an object of class \linkS4class{marginfit}.
#'
#' @return an object of class logLik
#' @export
#'
setMethod("logLik", "marginfit", function(object) {
  ll <- -object@fit$value
  attr(ll, "nobs") <- length(object@data)
  attr(ll, "df") <- length(object@fit$par)
  class(ll) <- "logLik"
  ll
})

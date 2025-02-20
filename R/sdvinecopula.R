#' Stationary d-vine copula processes
#'
#' Class of objects for stationary d-vine copula processes. See \link{sdvinecopula} for more details.
#'
#' @slot name name of the stationary d-vine copula process.
#' @slot modelspec list containing the family, rotation, and name of KPACF
#' @slot pars list comprising of the parameters.
#'
#' @export
#'
setClass("sdvinecopula", contains = "tscopula", slots = list(
  name = "character",
  modelspec = "list",
  pars = "list"
))

#' Constructor function for sdvinecopula process
#'
#' This function sets up a stationary d-vine process of finite or infinite order. In general, models models consist of
#' base copulas and substituted copulas, which form a copula sequence. The copulas are parameterized using the
#' Kendall partial autocorrelation function (kpacf) specified by the \code{kpacf} argument.
#'
#' The default choice is the
#' kpacf of a standard ARMA process which is implemented in the function \code{\link{kpacf_arma}}. The parameters
#' of the kpacf should be set as a named list using the \code{pars} argument; the required parameters should usually
#' be clear from the documentation of the chosen kpacf function and must be correctly named.
#'
#' The base copula sequence is specified by the name of a single 1-parameter copula in \code{basefamily},
#' which can be Gauss, Frank, Gumbel, Clayton or Joe. For choices other than Gauss or Frank, the
#' user must specify the rotation that should be used for positive dependencies (0 or 180) and the rotation that
#' should be used for negative dependencies (90 or 270). These are specified by the arguments \code{baseposrot} and
#' \code{basenegrot} respectively.
#'
#' The number of copulas to be replaced is determined by the length of the \code{family} vector, which should be
#' a vector of character strings representing copula names. The default value is \code{NULL} meaning no substitutions.
#' If \code{family} has length one, the same copula name is repeated for all substitutions.
#' The substituted families can be Gauss, Gumbel, Clayton, Joe, Frank, t and BB1 copulas as implemented by the
#' \code{\link[rvinecopulib]{bicop_dist}} in the \code{rvinecopulib} package.
#'
#' For the substituted copulas (other than Gauss, t and Frank) the user should specify the rotation that should be used for
#' positive dependencies (0 or 180) and the rotation that should be used for negative dependencies (90 or 270). These are
#' specified by the arguments \code{posrot} and \code{negrot} respectively, which have default values 0 and 90.
#' These vectors should be the same length as \code{family} or length one, in which case the values are repeated.
#'
#' In practice, the sequence of copulas will be truncated at the last copula for which the kpacf exceeds \code{tautol}.
#' The \code{maxlag} parameter is typically used to force the truncation to take place at a lower lag (to increase speed).
#' This can also be achieved by increasing the value of \code{tautol}.
#'
#' If one or more of the substituted copulas are t or BB1 copulas the argument \code{auxpar} should be used to
#' specify the additional parameters. These are the degree-of-freedom parameter for t and the delta parameter for BB1;
#' the former must be greater or equal 2 and the latter greater or equal 1.
#'
#' @param pars a list containing the parameters of the model
#' @param kpacf a character string giving the name of the Kendall pacf (default is kpacf_arma)
#' @param family vector of family names for copula substitutions
#' @param posrot vector of rotations for substituted families under positive dependence (default is 0)
#' @param negrot vector of rotations for substituted families under negative dependence (default is 90)
#' @param basefamily scalar specifying base copula family (default is "gauss")
#' @param baseposrot scalar specifying copula rotation under positive dependence (default is 0)
#' @param basenegrot scalar specifying copula rotation under negative dependence (default is 0)
#' @param auxpar vector of additional parameters for two-parameter copulas
#' @param tautol scalar value at which kpacf is truncated
#' @param maxlag a scalar which can be used to force a given value for maximum lag
#'
#' @return An object of class \linkS4class{sdvinecopula}.
#' @export
#'
#' @examples
#' sdvinecopula(pars = list(ar = 0.95, ma = 0.85), kpacf = "kpacf_arma",
#' family = c("Gumbel", "clayton"), posrot = c(0, 180), negrot = c(90, 270),
#' tautol = 1e-04)
sdvinecopula <- function(pars = list(ar = 0, ma = 0),
                         kpacf = "kpacf_arma",
                         family = NULL,
                         posrot = 0,
                         negrot = 90,
                         basefamily = "gauss",
                         baseposrot = 0,
                         basenegrot = 0,
                         auxpar = NA,
                         tautol = 1e-04,
                         maxlag = Inf
                         ) {
  if (is.null(names(pars)))
    stop("parameters must be named")
  arflag <- FALSE
  if ((kpacf == "kpacf_arma") | (kpacf == "kpacf_sarma4") | (kpacf == "kpacf_sarma12")){
    period <- switch(kpacf,
                     "kpacf_arma" = 0,
                     "kpacf_sarma4" = 4,
                     "kpacf_sarma1" = 12)
    if ((length(pars$ma) == 0) & (length(pars$sma) == 0)){
      if ((length(pars$ar) > 0) | (length(pars$sar) > 0)){
        maxlag <- length(pars$ar) + period * length(pars$sar)
        arflag <- TRUE
      }
    }
  }
  basefamily <- tolower(basefamily)
  if ((basefamily == "frank") | (basefamily == "gauss")){
    baseposrot <- 0
    basenegrot <- 0 # radial symmetry
  }
  nreplace <- length(family)
  if (nreplace > 0){
    family <- tolower(family)
    if (arflag & (nreplace > maxlag)){
      stop("Too many copula substitutions for AR process")
    }
    if (nreplace > 1)
    {
      if (length(family) == 1)
        family <- rep(family, nreplace)
      if (length(posrot) == 1)
        posrot <- rep(posrot, nreplace)
      if (length(negrot) == 1)
        negrot <- rep(negrot, nreplace)
    }
    if ((length(posrot) != nreplace) | (length(negrot) != nreplace) | (length(family) != nreplace))
      stop("Length of family and rotations must be equal")
    posrot[(family == "frank") | (family == "t") | (family == "gauss")] <- 0 # radial symmetry
    negrot[(family == "frank") | (family == "t") | (family == "gauss")] <- 0 # radial symmetry
    twoparfamily <- family[(family == "t") | (family == "bb1")]
    ntwopar <- length(twoparfamily)
    if (ntwopar > 0){
      if (length(auxpar) != ntwopar)
        stop("Require second parameter for two-parameter copulas")
     pars$auxpar <- auxpar
    }
  }
  modelspec <- list(nreplace = nreplace,
                    family = family,
                    posrot = posrot,
                    negrot = negrot,
                    basefamily = basefamily,
                    baseposrot = baseposrot,
                    basenegrot = basenegrot,
                    kpacf = kpacf,
                    tautol = tautol,
                    maxlag = maxlag,
                    arflag = arflag,
                    npar = length(unlist(pars)))
  output <- new("sdvinecopula",
      name = paste("stationary d-vine"),
      modelspec = modelspec,
      pars = pars
  )
  if (nreplace > 0){ # check replacements only, or all copulas in AR process
    ncheck <- ifelse(arflag, maxlag, nreplace)
    check <- mklist_sdvine(output, ncheck)
  }
  output
}


#' @describeIn sdvinecopula Coef Method for sdvinecopula class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("coef", c(object = "sdvinecopula"), function(object) {
  unlist(object@pars)
})

#' Find effective maximum lag
#'
#' @param tau values of Kendall partial autocorrelation function.
#' @param modelspec sdvinecopula model specification
#'
#' @return a value for the effective maximum lag floored at one.
#' @keywords internal
#'
effective_kmax <- function(tau, modelspec){
  if (modelspec$arflag)
    return(modelspec$maxlag)
  else
    return(min(max(1,which(abs(tau) > modelspec$tautol)), modelspec$maxlag))
}

#' @describeIn sdvinecopula Calculate Kendall's tau values for pair copulas in type 3 d-vine copula
#'
#' @param object an object of the class.
#' @param lagmax maximum value of lag to be considered.
#'
#' @export
#'
setMethod("kendall", c(object = "sdvinecopula"), function(object, lagmax = 20) {
  kpacf <- eval(parse(text = object@modelspec$kpacf))
  tau <- kpacf(lagmax, object@pars)
  k <- effective_kmax(tau, object@modelspec)
  tau[which((1:lagmax) > k)] <- 0
  tau
}
)


#' @describeIn sdvinecopula Show method for sdvinecopula class
#'
#' @param object an object of the class.
#'
#' @export
#'
setMethod("show", c(object = "sdvinecopula"), function(object) {
  cat("object class: ", is(object)[[1]], "\n", sep = "")
  cat("name: ", object@name, "\n", sep = "")
  cat("kpacf: ", object@modelspec$kpacf,"\n", sep = "")
  nreplace <- object@modelspec$nreplace
  if (nreplace > 0){
    cat(nreplace, "explicit copula substitutions:\n")
    tau <- kendall(object, lagmax = nreplace)
    family <- object@modelspec$family
    rot <- as.character(ifelse(tau >= 0, object@modelspec$posrot,object@modelspec$negrot))
    rot[family %in% c("gauss", "frank", "t")] <- ""
    cat(" - families:", paste(family, rot, sep =""), "\n", sep = " ")
  }
  if (object@modelspec$maxlag > object@modelspec$nreplace){
    cat("base family:", object@modelspec$basefamily, "\n")
    if (!(object@modelspec$basefamily %in% c("gauss", "frank")))
      cat(" - positive and negative rotations:", object@modelspec$baseposrot, object@modelspec$basenegrot, "\n")
  }
  tau <- kendall(object, lagmax = 100) # hard coded limit
  EML <- effective_kmax(tau, object@modelspec)
  ML <- object@modelspec$maxlag
  if (ML <= EML)
    cat(" - maximum lag is", ML, "\n")
  else
    cat(" - effective maximum lag is", EML,
        "at tolerance", object@modelspec$tautol, "\n")
  cat("parameters: \n")
  print(coef(object))
})



#' Make list of pair copulas for sdvinecopula object
#'
#' @param x an object of class sdvinecopula
#' @param maxlag maximum possible lag to consider
#'
#' @return a list of pair copulas
#' @keywords internal
#'
mklist_sdvine <- function(x, maxlag){
  kpacf <- eval(parse(text = x@modelspec$kpacf))
  tauvals <- kpacf(maxlag, x@pars)
  k <- effective_kmax(tauvals, x@modelspec)
  pc_list <- vector("list", k)
  auxpar <- 1
  nreplace <- x@modelspec$nreplace
  for (i in 1:k) {
    if (i <= nreplace){
      fam <- x@modelspec$family[i]
      if (tauvals[i] >= 0)
        rot <- x@modelspec$posrot[i]
      if (tauvals[i] < 0)
        rot <- x@modelspec$negrot[i]
    }
    else{
      fam <- x@modelspec$basefamily
      if (tauvals[i] >= 0)
        rot <- x@modelspec$baseposrot
      if (tauvals[i] < 0)
        rot <- x@modelspec$basenegrot
      }
    coppars <- ktau_to_par(
      family = fam,
      tau = tauvals[i]
    )
    if (fam == "t"){
      coppars <- c(coppars, x@pars$auxpar[auxpar])
      auxpar <- auxpar + 1
    }
    if (fam == "bb1"){
      par2 <- x@pars$auxpar[auxpar]
      par1 <- (coppars+2)/par2 -2
      coppars <- c(par1, par2)
      auxpar <- auxpar + 1
    }
    pc_list[[i]] <- rvinecopulib::bicop_dist(
      family = fam,
      rotation = rot,
      parameters = coppars)
  }
  pc_list
}

#' @describeIn sdvinecopula Simulation method for sdvinecopula class
#'
#' @param object an object of the class.
#' @param n length of realization.
#'
#' @export
#'
setMethod("sim", c(object = "sdvinecopula"), function(object, n = 1000) {
  pc_list <- mklist_sdvine(object, n-1)
  simdvine(pc_list, n, innov = NA, start = NA)
})


#' Objective function for sdvinecopula process
#'
#' @param theta parameters of kpacf
#' @param modelspec list specifying model
#' @param u data
#'
#' @return Value of objective function at parameters.
#'
#' @keywords internal
#'
sdvinecopula_objective <- function(theta, modelspec, u) {
  n <- length(u)
  kpacf <- eval(parse(text = modelspec$kpacf))
  tauvals <- kpacf((n-1), theta)
  if (is.na(sum(tauvals)))
    return(NA)
  k <- effective_kmax(tauvals, modelspec)
  pc_list <- vector("list", k)
  tpar <- 1
  nreplace <- modelspec$nreplace
  for (i in 1:k) {
    if (i <= nreplace){
      fam <- modelspec$family[i]
      if (tauvals[i] >= 0)
        rot <- modelspec$posrot[i]
      if (tauvals[i] < 0)
        rot <- modelspec$negrot[i]
    }
    else{
      fam <- modelspec$basefamily
      if (tauvals[i] >= 0)
        rot <- modelspec$baseposrot
      if (tauvals[i] < 0)
        rot <- modelspec$basenegrot
    }
    coppars <- ktau_to_par(
      family = fam,
      tau = tauvals[i]
    )
    if (fam == "t"){
      auxpar <- theta[substring(names(theta),1,2) == "au"]
      coppars <- c(coppars, auxpar[tpar])
      tpar <- tpar + 1
    }
    if (fam == "bb1"){
      auxpar <- theta[substring(names(theta),1,2) == "au"]
      par2 <- auxpar[tpar]
      par1 <- (coppars+2)/par2 -2
      coppars <- c(par1, par2)
      tpar <- tpar + 1
    }
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

#' @describeIn sdvinecopula Prediction method for sdvinecopula class
#'
#' @param object an object of the class.
#' @param data vector of past data values.
#' @param x vector of arguments of prediction function.
#' @param type type of prediction function ("df" for density, "qf" for quantile function
#' or "dens" for density).
#'
#' @export
#'
setMethod("predict", c(object = "sdvinecopula"), function(object, data, x, type = "df") {
  pc_list <- mklist_sdvine(object, length(data)-1)
  switch(type,
         "df" = Rblatt(pc_list, data, x),
         "qf" = IRblatt(pc_list, data, x),
         "dens" = Rblattdens(pc_list, data, x))

})

#' Residual function for sdvinecopula object
#'
#' @param object a fitted sdvinecopula object.
#' @param data the data to which copula is fitted.
#' @param trace extract trace instead of residuals.
#'
#' @return vector of model residuals
#' @keywords internal
#'
resid_sdvinecopula <- function(object, data = NA, trace = FALSE){
  n <- length(data)
  pc_list <- mklist_sdvine(object, n-1)
  k <- length(pc_list)
  for (i in 1:k)
    if (pc_list[[i]]$rotation %in% c(90,270))
      pc_list[[i]]$rotation <- 360 - pc_list[[i]]$rotation
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



#' Generalized lagging for fitted sdvinecopula objects
#'
#' @param copula a sdvinecopula object
#' @param data the data to which copula is fitted
#' @param lagmax the maximum lag value.
#' @param glagplot logical value indicating generalized lag plot.
#'
#' @return If \code{glagplot} is \code{TRUE} a list of generalized lagged datasets
#' of maximum length 9 is returned to facilitate a generalized lagplot.
#' If \code{glagplot} is \code{FALSE} a vector of length \code{lagmax} containing
#' the Kendall rank correlations for the generalized lagged datasets is returned.
#' @keywords internal
glag_for_sdvinecopula <- function(copula, data, lagmax, glagplot = FALSE) {
  if (glagplot)
    lagmax <- min(lagmax, 9)
  pc_list <- mklist_sdvine(copula, lagmax)
  k <- length(pc_list)
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
  if (k >1){
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

#' Automatic estimation of stationary d-vine copula with Gaussian base family
#'
#' This function begins with a model based exclusively on Gaussian pair copulas and
#' systematically replaces a finite number of copulas with non-Gaussian alternatives that
#' lower a specified information criterion.
#'
#' The number of copulas to replace is specified by \code{nreplace} but the algorithm may stop
#' prematurely if the information criterion is barely lowered over \code{nstrike} consecutive steps.
#'
#' This is a slow function to execute. Increasing \code{tautol} or setting a value for
#' \code{maxlag} will speed it up at the cost of some accuracy in the final model. The default
#' choices of copula are all one-parameter but the t copula can be added if desired.
#' For copulas without radial symmetry different rotations are tried to accommodate negative
#' and positive dependencies.
#'
#' @param model an object of class \linkS4class{tscopulafit} such that \code{model@tscopula} is of class
#' \linkS4class{armacopula} or \linkS4class{sarmacopula}.
#' @param nreplace fixed number of copulas to replace
#' @param tautol tolerance for setting the effective maximum lag in list of copulas
#' @param maxlag fixed value for maximum lag
#' @param ICtol tolerance for deciding on a premature stop to copula replacement
#' @param nstrike number of negligible improvements before premature stop
#' @param criterion information criterion to be used (AIB, BIC or AICc)
#' @param choices vector of copula names to be used for replacements
#' @param verbose logical parameter specifying whether a verbose output is desired
#'
#' @return An object of class \linkS4class{tscopulafit} such that \code{model@tscopula} is of class
#' \linkS4class{sdvinecopula}.
#' @export
#'
auto_dvine <- function(model,
                       nreplace = 5,
                       tautol = 1e-04,
                       maxlag = Inf,
                       ICtol = 0.5,
                       nstrike = 3,
                       criterion = "AIC",
                       choices = c("gumbel", "clayton", "frank", "joe"),
                       verbose = TRUE) {
  if (is(model@tscopula, "armacopula"))
    kpacf <- "kpacf_arma"
  else if (is(model@tscopula, "sarmacopula"))
    kpacf <- "kpacf_sarma"
  else
    stop("Initial model should be an ARMA or SARMA copula")
  order <- as.numeric(model@tscopula@modelspec)
  if (kpacf == "kpacf_sarma") {
    kpacf <- paste(kpacf, order[length(order)], sep = "")
    order <- order[-length(order)]
  }
  IC <- switch(criterion,
               AIC = AIC,
               BIC = BIC,
               AICc = AICc)

  if (verbose) {
    cat("Model order : ", order, "\n")
    cat(
      "Initial Gaussian model: ",
      criterion,
      IC(model),
      "; convergence",
      model@fit$convergence,
      "\n"
    )
  }
  start <- model@tscopula@pars
  data <- model@data
  bestmod <- model
  improvement <- FALSE
  copulas <- rep("gauss", nreplace)
  posrots <- rep(NA, nreplace)
  negrots <- rep(NA, nreplace)
  auxpar <- NA # for t and BB1 copulas
  strike <- 0
  for (k in 1:nreplace) {
    lastmod <- bestmod
    if (verbose)
      cat("Determining copula", k, ": ")
    for (cop in choices) {
      if (verbose)
        cat(cop, " ")
      if ((cop == choices[length(choices)]) & verbose)
        cat("; selected ")
      notradial <- (cop != "frank") &
        (cop != "gauss") & (cop != "t")
      copulas[k] <- cop
      posrots[k] <- 0
      negrots[k] <- 90
      if (cop == "t") {
        if (is.na(auxpar[1]))
          tryauxpar <- 10
        else
          tryauxpar <- c(auxpar, 10)
      }
      else
        tryauxpar <- auxpar
      mod <- sdvinecopula(
        family = copulas[1:k],
        posrot = posrots[1:k],
        negrot = negrots[1:k],
        kpacf = kpacf,
        pars = start,
        auxpar = tryauxpar,
        tautol = tautol,
        maxlag = maxlag
      )
      modfit <- fit(mod, data)
      if (IC(modfit) < IC(bestmod)) {
        improvement <- TRUE
        bestmod <- modfit
      }
      if ((kendall(modfit)[k] > 0) & notradial)
      {
        posrots[k] <- 180
        mod <- sdvinecopula(
          family = copulas[1:k],
          posrot = posrots[1:k],
          negrot = negrots[1:k],
          kpacf = kpacf,
          pars = start,
          auxpar = tryauxpar,
          tautol = tautol,
          maxlag = maxlag
        )
        modfit <- fit(mod, data)
        if (IC(modfit) < IC(bestmod)) {
          bestmod <- modfit
          improvement <- TRUE
        }
      }
      if ((kendall(modfit)[k] < 0) & notradial)
      {
        negrots[k] <- 270
        mod <- sdvinecopula(
          family = copulas[1:k],
          posrot = posrots[1:k],
          negrot = negrots[1:k],
          kpacf = kpacf,
          pars = start,
          auxpar = tryauxpar,
          tautol = tautol,
          maxlag = maxlag
        )
        modfit <- fit(mod, data)
        if (IC(modfit) < IC(bestmod)) {
          improvement <- TRUE
          bestmod <- modfit
        }
      }
    }
    if (improvement) {
      copulas[k] <- bestmod@tscopula@modelspec$family[k]
      posrots[k] <- bestmod@tscopula@modelspec$posrot[k]
      negrots[k] <- bestmod@tscopula@modelspec$negrot[k]
      if (copulas[k] == "t") {
        auxpar <- bestmod@tscopula@pars$auxpar
      }
      start <- bestmod@tscopula@pars
      improvement <- FALSE
    }
    else
    {
      copulas[k] <- "gauss"
      posrots[k] <- 0
      negrots[k] <- 0
    }
    showrot <- ""
    if ((copulas[k] != "frank") &
        (copulas[k] != "gauss") & (copulas[k] != "t"))
      showrot <- as.character(ifelse(kendall(bestmod)[k] > 0, posrots[k], negrots[k]))
    if (verbose)
      cat(
        copulas[k],
        showrot,
        " ; ",
        criterion,
        IC(bestmod),
        "; convergence",
        bestmod@fit$convergence,
        "; K =",
        bestmod@fit$EML,
        "\n"
      )
    magnitude <- IC(lastmod) - IC(bestmod)
    if (magnitude < ICtol) {
      strike <- strike + 1
      if (verbose)
        cat("Strike ", strike, "\n")
    }
    else
      strike <- 0
    if (strike == nstrike)
      return(bestmod)
  }
  bestmod
}

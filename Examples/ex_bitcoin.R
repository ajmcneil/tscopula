data(bitcoin)
X <- (diff(log(bitcoin))[-1]) * 100 # log-returns (as percentages)
length(X)
plot(X)
U <- strank(X)

tsoptions <- list(hessian = TRUE, method = "Nelder-Mead", avoidzero= FALSE)

copmod_Gauss <- armacopula(pars = list(ar = 0.95, ma = -0.85))
mod_Gauss <- vtscopula(copmod_Gauss, Vlinear(0.45))
fit_Gauss <- fit(mod_Gauss, U, tsoptions = tsoptions)

copmod_Frank <- dvinecopula2(family = "frank",
                             pars = list(ar = 0.95, ma = -0.85),
                             maxlag = 30)
mod_Frank <- vtscopula(copmod_Frank, Vlinear(0.45))
fit_Frank <- fit(mod_Frank, U, tsoptions = tsoptions)

copmod_Gauss2 <- armacopula(list(ar = 0.96, ma = -0.84))
mod_Gauss2 <- vtscopula(copmod_Gauss2, V2p(delta=0.42, kappa=1))
fit_Gauss2 <- fit(mod_Gauss2, U, tsoptions = tsoptions)

copmod_Frank2 <- dvinecopula2(family = "frank",
                              pars = list(ar = 0.96, ma = -0.84),
                              maxlag = 30)
mod_Frank2 <- vtscopula(copmod_Frank2, V2p(delta=0.42, kappa=1))
fit_Frank2 <- fit(mod_Frank2, U, tsoptions = tsoptions)

copmod_Frank3 <- dvinecopula2(family = "frank",
                             pars = list(ar = 0.95, ma = -0.85),
                             maxlag = Inf)
mod_Frank3 <- vtscopula(copmod_Frank3, Vlinear(0.45))
fit_Frank3 <- fit(mod_Frank3, U, tsoptions = tsoptions)

AIC(fit_Gauss, fit_Gauss2, fit_Frank, fit_Frank2, fit_Frank3)
#
plot(fit_Frank3, plotoption = 1)
plot(fit_Frank3, plotoption = 2)
plot(fit_Frank3, plotoption = 3)
plot(fit_Frank, plottype = "vtransform")

marg_dweibull <- fit(margin("doubleweibull",
                            pars = c(mu = 0.2, shape =0.8, scale = 2.7)), X)
mod_dweibull <- tscm(fit_Frank3, margin = marg_dweibull)
mod_dweibull <- fit(mod_dweibull, as.numeric(X), tsoptions = tsoptions)
AIC(mod_dweibull)
plot(mod_dweibull, plotoption = 1)
plot(mod_dweibull, plotoption = 2)
plot(mod_dweibull, plotoption = 3)
plot(mod_dweibull, plottype = "margin")
plot(mod_dweibull, plotoption = 2, plottype = "margin")
plot(mod_dweibull, plottype = "vtransform")
plot(mod_dweibull, plottype = "volprofile")
plot(mod_dweibull, plottype = "volproxy")
plot(mod_dweibull, plotoption = 2, plottype = "volproxy")

library(stats4)
suppressMessages(library(quantmod))
BTCUSD <- suppressWarnings(getSymbols("BTCUSD=X", auto.assign = FALSE,
                                      from = "2015-12-31", to = "2019-12-31"))
BTCUSD <- na.omit(BTCUSD) # omit missing values
BTCUSD <- BTCUSD$`BTCUSD=X.Close` # closing values
plot(BTCUSD)

X <- (diff(log(BTCUSD))[-1]) * 100 # log-returns (as percentages)
length(X)
plot(X)

V <- strank(abs(X))


armamod_Gauss <- dvinecopula2(family = "gauss", pars = list(ar = 0.1, ma = 0.1), maxlag = 30)
fit_Gauss <-fit(armamod_Gauss, V)
armamod_Joe <- dvinecopula2(family = "joe", pars = list(ar = 0, ma = 0), maxlag = 30)
fit_Joe <- fit(armamod_Joe, V)
armamod_Gumbel <- dvinecopula2(family = "gumbel", pars = list(ar = 0, ma = 0), maxlag = 30)
fit_Gumbel <- fit(armamod_Gumbel, V)
armamod_Frank <- dvinecopula2(family = "frank", pars = list(ar = 0, ma = 0), maxlag = 30)
fit_Frank <-fit(armamod_Frank, V)
AIC(fit_Gauss, fit_Joe, fit_Gumbel, fit_Frank)

armamod_Gauss <- dvinecopula2(family = "gauss", kpacf = "kpacf_exp",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Gauss2 <-fit(armamod_Gauss, V)
armamod_Joe <- dvinecopula2(family = "joe", kpacf = "kpacf_exp",
                            pars = list(c(-1,0)), maxlag = 30)
fit_Joe2 <- fit(armamod_Joe, V)
armamod_Gumbel <- dvinecopula2(family = "gumbel", kpacf = "kpacf_exp",
                               pars = list(c(-1,0)), maxlag = 30)
fit_Gumbel2 <- fit(armamod_Gumbel, V)
armamod_Frank <- dvinecopula2(family = "frank", kpacf = "kpacf_exp",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Frank2 <-fit(armamod_Frank, V)
AIC(fit_Gauss2, fit_Joe2, fit_Gumbel2, fit_Frank2)

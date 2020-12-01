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

#

armamod_Gauss <- dvinecopula2(family = "gauss", kpacf = "kpacf_exp",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Gauss2 <-fit(armamod_Gauss, V)
armamod_Joe <- dvinecopula2(family = "joe", kpacf = "kpacf_exp",
                            pars = list(c(-1,0)), maxlag = 30)
fit_Joe2 <- fit(armamod_Joe, V)
armamod_Gumbel <- dvinecopula2(family = "gumbel", kpacf = "kpacf_exp",
                               pars = list(c(-1,-0.2)), maxlag = 30)
fit_Gumbel2 <- fit(armamod_Gumbel, V)
armamod_Frank <- dvinecopula2(family = "frank", kpacf = "kpacf_exp",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Frank2 <-fit(armamod_Frank, V)
AIC(fit_Gauss2, fit_Joe2, fit_Gumbel2, fit_Frank2)
#

armamod_Gauss <- dvinecopula2(family = "gauss", kpacf = "kpacf_pow",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Gauss3 <-fit(armamod_Gauss, V)
armamod_Joe <- dvinecopula2(family = "joe", kpacf = "kpacf_pow",
                            pars = list(c(-1,0)), maxlag = 30)
fit_Joe3 <- fit(armamod_Joe, V)
armamod_Gumbel <- dvinecopula2(family = "gumbel", kpacf = "kpacf_pow",
                               pars = list(c(-1,0)), maxlag = 30)
fit_Gumbel3 <- fit(armamod_Gumbel, V)
armamod_Frank <- dvinecopula2(family = "frank", kpacf = "kpacf_pow",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Frank3 <-fit(armamod_Frank, V)
AIC(fit_Gauss3, fit_Joe3, fit_Gumbel3, fit_Frank3)

#

armamod_Gauss <- dvinecopula2(family = "gauss", kpacf = "kpacf_fbn",
                              pars = list(c(1)), maxlag = 30)
fit_Gauss4 <-fit(armamod_Gauss, V)
armamod_Joe <- dvinecopula2(family = "joe", kpacf = "kpacf_fbn",
                            pars = list(c(1)), maxlag = 30)
fit_Joe4 <- fit(armamod_Joe, V)
armamod_Gumbel <- dvinecopula2(family = "gumbel", kpacf = "kpacf_fbn",
                               pars = list(c(1)), maxlag = 30)
fit_Gumbel4 <- fit(armamod_Gumbel, V)
armamod_Frank <- dvinecopula2(family = "frank", kpacf = "kpacf_fbn",
                              pars = list(c(1)), maxlag = 30)
fit_Frank4 <-fit(armamod_Frank, V)
AIC(fit_Gauss4, fit_Joe4, fit_Gumbel4, fit_Frank4)

#

armamod_Gauss <- dvinecopula2(family = "gauss", kpacf = "kpacf_arfima1",
                              pars = list(phi = 0.9, theta = 0.8, H = 0), maxlag = 30)
fit_Gauss5 <-fit(armamod_Gauss, V)
armamod_Joe <- dvinecopula2(family = "joe", kpacf = "kpacf_arfima1",
                            pars = list(phi = 0.9, theta = 0.8, H = 0), maxlag = 30)
fit_Joe5 <- fit(armamod_Joe, V)
armamod_Gumbel <- dvinecopula2(family = "gumbel", kpacf = "kpacf_arfima1",
                               pars = list(phi = 0.9, theta = 0.8, H = 0), maxlag = 30)
fit_Gumbel5 <- fit(armamod_Gumbel, V)
armamod_Frank <- dvinecopula2(family = "frank", kpacf = "kpacf_arfima1",
                              pars = list(phi = 0.9, theta = 0.8, H = 0), maxlag = 30)
fit_Frank5 <-fit(armamod_Frank, V)
AIC(fit_Gauss5, fit_Joe5, fit_Gumbel5, fit_Frank5)

fit_Frank
u <- sim(fit_Frank5)
acf(qnorm(u), 50)

#

armamod_Gauss <- dvinecopula2(family = "gauss", kpacf = "kpacf_exp2",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Gauss6 <-fit(armamod_Gauss, V)
armamod_Joe <- dvinecopula2(family = "joe", kpacf = "kpacf_exp2",
                            pars = list(c(-1,0)), maxlag = 30)
fit_Joe6 <- fit(armamod_Joe, V)
armamod_Gumbel <- dvinecopula2(family = "gumbel", kpacf = "kpacf_exp2",
                               pars = list(c(-1,0)), maxlag = 30)
fit_Gumbel6 <- fit(armamod_Gumbel, V)
armamod_Frank <- dvinecopula2(family = "frank", kpacf = "kpacf_exp2",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Frank6 <-fit(armamod_Frank, V)
AIC(fit_Gauss6, fit_Joe6, fit_Gumbel6, fit_Frank6)

#

armamod_Gauss <- dvinecopula2(family = "gauss", kpacf = "kpacf_pow2",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Gauss7 <-fit(armamod_Gauss, V)
armamod_Joe <- dvinecopula2(family = "joe", kpacf = "kpacf_pow2",
                            pars = list(c(-1,0)), maxlag = 30)
fit_Joe7 <- fit(armamod_Joe, V)
armamod_Gumbel <- dvinecopula2(family = "gumbel", kpacf = "kpacf_pow2",
                               pars = list(c(-1,0)), maxlag = 30)
fit_Gumbel7 <- fit(armamod_Gumbel, V)
armamod_Frank <- dvinecopula2(family = "frank", kpacf = "kpacf_pow2",
                              pars = list(c(-1,0)), maxlag = 30)
fit_Frank7 <-fit(armamod_Frank, V)
AIC(fit_Gauss7, fit_Joe7, fit_Gumbel7, fit_Frank7)


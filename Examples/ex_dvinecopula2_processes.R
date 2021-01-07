# A model with Joe copula
mod <- dvinecopula2(family = "joe",
                        pars = list(ar = 0.97, ma = -0.85),
                        maxlag = 30)
mod

set.seed(113)
data <- sim(mod, 1000)
ts.plot(data)


armamod_Gauss <- dvinecopula2(family = "gauss",
                              pars = list(ar = 0, ma = 0),
                              maxlag = 30)
fitGauss <- fit(armamod_Gauss, data)
fitGauss

armamod_Joe <- dvinecopula2(family = "joe",
                            pars = list(ar = 0, ma = 0),
                            maxlag = 30)
fitJoe <- fit(armamod_Joe, data)
fitJoe

plot(fitJoe, plotoption = 1)
plot(fitJoe, plotoption = 2)
plot(fitJoe, plotoption = 3)


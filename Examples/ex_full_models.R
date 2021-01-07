# simulate and fit Gaussian ARMA
set.seed(13)
data <- 0.5 + 2*arima.sim(list(ar =0.95, ma =-0.85), 1000)

armacop <- armacopula(pars = list(ar =0.01, ma =0.01))
fit0 <- fit(tscm(armacop, margin("norm")), data)
fit <- fit(fit0, data, method = "full")

plot(fit, plotoption = 1)
plot(fit, plotoption = 2)
plot(fit, plotoption = 3)
plot(fit, plottype = "margin")
plot(fit, plotoption = 2, plottype = "margin")





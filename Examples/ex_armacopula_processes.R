# create armacopula process specification
arma11 <- armacopula(list(ar = 0.95, ma = -0.85))
arma11

# simulate a realisation
set.seed(13)
data <- sim(arma11, 1000)
ts.plot(data)

# fit ARMA(1,1) specification to data
arma11spec <- armacopula(list(ar = 0.1, ma = 0.1))
modfit <- fit(arma11spec, data)
modfit
coef(modfit)

# plot fit
plot(modfit, plotoption = 1)
plot(modfit, plotoption = 2)
plot(modfit, plotoption = 3)
plot(modfit, plotoption = 4)
plot(modfit, plotoption = 5)

# resimulate fitted model
data2 <- sim(modfit, 1000)
ts.plot(data2)

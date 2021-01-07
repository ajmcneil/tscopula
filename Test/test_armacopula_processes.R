# create armacopula process specification
ar1 <- armacopula(list(ar = 0.7))
ar1

# simulate a realisation
data <- sim(ar1, 1000)
ts.plot(data)

# fit specification to data
ar1fit <- fit(ar1, data)
ar1fit

# create armacopula process specification
arma11 <- armacopula(list(ar = 0.95, ma = -0.85))
arma11

# simulate a realisation
data <- sim(arma11, 1000)
ts.plot(data)

# fit specification to data
tmp <- fit(arma11, data)
tmp
coef(tmp)

# plot fit
plot(tmp, type = 1)
plot(tmp, type = 2)
plot(tmp, type = 3)
plot(tmp, type = 4)
plot(tmp, type = 5)

# apply Kalman filter to data
head(kfilter(tmp@tscopula, data))
tail(kfilter(tmp@tscopula, data))

# resimulate fitted model
tmp2 <- sim(tmp, 1000)
ts.plot(tmp2)

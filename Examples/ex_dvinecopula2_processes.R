# set a high accuracy for fits
acc <- 10*(.Machine$double.eps)

# AR(3) example
set.seed(13)
truear <- c(0.5,0.2,-0.3)
ar3mod <- armacopula(pars = list(ar = truear))
u1 <- sim(ar3mod, n = 1000)
ar3modB <- dvinecopula2(pars = list(ar = c(0.2, 0.1, 0.05)))
fit11 <- fit(ar3mod, u1, control = list(reltol = acc))
fit11
fit12 <- fit(ar3modB, u1, control = list(reltol = acc))
fit12
coef(fit12)
#############################################
# ARMA(1,1) example
set.seed(13)
truear <- 0.95
truema <- -0.85
armamod <- armacopula(list(ar = truear, ma = truema))
u2 <- sim(armamod, 1000)
armamodB <- dvinecopula2(pars = list(ar = 0.1, ma = 0.1))
fit21 <- fit(armamod, u2, control = list(reltol = acc))
fit21
fit22 <- fit(armamodB, u2, control = list(reltol = acc))
fit22
coef(fit22)
#######################################
# MA(2) example
set.seed(13)
truema <- c(0.5, 0.2)
mamod <- armacopula(list(ma = truema))
u3 <- sim(mamod, 1000)
mamodB <- dvinecopula2(pars = list(ma = rep(0.1,2)))
fit31 <- fit(mamod, u3, control = list(reltol = acc))
fit31
fit32 <- fit(mamodB, u3, control = list(reltol = acc))
fit32
coef(fit32)
##########################################

pacf0 <- pacf(qnorm(u2), lag.max = 30, plot = FALSE)
plot(2/pi * asin(pacf0$acf[,,1]), type = "h", ylab = "KPACF")
lines(kpacf_arma(k = 1:30, theta = fit22@fit$par), col = "red")

# Faster fits with maxlag and default tolerance
armamodB2 <- dvinecopula2(pars = list(ar = 0.1, ma = 0.1), maxlag = 30)
armamodB2
fit22B <- fit(armamodB2, u2)
fit22B
fit22

# Simulation needs a maxlag or it takes forever
simmod <- dvinecopula2(pars = list(ar = 0.95, ma = -0.85), maxlag = 30)
simdata <- sim(simmod, 1000)
ts.plot(simdata)
fit(armamodB2, simdata)

# Now a model with Joe copula
simmod2 <- dvinecopula2(family = "joe", pars = list(ar = 0.95, ma = -0.85), maxlag = 30)
set.seed(13)
simdata2 <- sim(simmod2, 1000)
ts.plot(simdata2)
armamod_Gauss <- dvinecopula2(family = "gauss", pars = list(ar = 0.1, ma = 0.1), maxlag = 30)
fit(armamod_Gauss, simdata2)
armamod_Joe <- dvinecopula2(family = "joe", pars = list(ar = 0, ma = 0), maxlag = 30)
fit(armamod_Joe, simdata2)

armamod_Joe2 <- dvinecopula2(family = "joe", kpacf = "kpacf_exp",
                            pars = list(c(-1, 0)), maxlag = 30)
fit(armamod_Joe2, simdata2)

# This seems like a reasonable model
goodmod <- dvinecopula2(family = "joe", kpacf = "kpacf_exp",
                       pars = list(c(-2, -0.2)), maxlag = 30)
gooddata <- sim(goodmod, 2000)
ts.plot(gooddata)
hist(gooddata)

# This seems like a model that should be forbidden
oddmod <- dvinecopula2(family = "joe", kpacf = "kpacf_exp",
                       pars = list(c(-1, -0.1)), maxlag = 2)
odddata <- sim(oddmod, 10000)
ts.plot(odddata)
hist(odddata)

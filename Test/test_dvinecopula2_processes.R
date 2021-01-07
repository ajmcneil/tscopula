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



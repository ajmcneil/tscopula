# simulate Gaussian ARMA
set.seed(13)
data <- 0.5 + 2*arima.sim(list(ar =0.95, ma =-0.85), 1000)
# two copula-based models
armacop <- armacopula(pars = list(ar =0.01, ma =0.01))
dvinecop <- dvinecopula2(pars = list(ar =0.01, ma =0.01))

# 3 ways to fit a Gaussian ARMA
acc <- (.Machine$double.eps)*10
fit1 <- arima(data, order = c(1,0,1), method = "ML",
              optim.method = "Nelder-Mead",
              optim.control = list(reltol = acc))
# 2-stage copula fits
fit20 <- fit(tscm(armacop, margin("norm")), data)
fit30 <- fit(tscm(dvinecop, margin("norm")), data)
coef(fit1)
coef(fit20)
coef(fit30)
# full copula fits
fit2 <- fit(fit20, data, method = "full",
            control = list(reltol = acc))
fit3 <- fit(fit30, data, method = "full",
            control = list(reltol = acc))
coef(fit1)
coef(fit2)
coef(fit3)
logLik(fit1)
logLik(fit2)
logLik(fit3)


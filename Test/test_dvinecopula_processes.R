# AR1
set.seed(1)
ar1mod <- armacopula(pars = list(ar = 0.4))
ar1modB <- dvinecopula(family = "gauss", pars = list(0.4))
acc <- 10 * .Machine$double.eps
y <- sim(ar1mod, n = 1000)
fit1 <- fit(ar1mod, y, control = list(reltol = acc, warn.1d.NelderMead = F))
fit1
fit2 <- fit(ar1modB, y, control = list(factr = acc, warn.1d.NelderMead = F))
fit2

# AR2
set.seed(2)
ar2mod <- armacopula(pars = list(ar = c(0.5, 0.2)))
ar2modB <- dvinecopula(family = "gauss", pars = list(0.5, 0.2))
acc <- 10 * .Machine$double.eps
y <- sim(ar2mod, n = 1000)
fit1 <- fit(ar2mod, y, control = list(reltol = acc))
fit1
fit2 <- fit(ar2modB, y, control = list(factr = acc))
fit2
ARMAacf(ar = fit1@tscopula@pars$ar, pacf = TRUE)

# AR3
set.seed(3)
ar3mod <- armacopula(pars = list(ar = c(0.5, 0.2, -0.3)))
ar3modB <- dvinecopula(family = "gauss", pars = list(0.5, 0.2, -0.3))
acc <- 10 * .Machine$double.eps
y <- sim(ar3mod, n = 1000)
fit1 <- fit(ar3mod, y, control = list(reltol = acc))
fit1
fit2 <- fit(ar3modB, y, control = list(factr = acc))
fit2
ARMAacf(ar = fit1@tscopula@pars$ar, pacf = TRUE)

# Mixed model
mixmod <- dvinecopula(family = c("Gaussian", "Clayton", "Gumbel"), pars = list(0.3, 1.0, 1.5))
mixmod
y <- sim(mixmod, n = 1000)
ts.plot(y)
fit <- fit(mixmod, y, list(hessian = TRUE))
fit
coef(fit)





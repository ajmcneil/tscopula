mixmod <- dvinecopula(
  family = c("Gaussian", "Frank", "Clayton","t"),
  pars = list(0.3, 2, 1.2, c(0.4,5)),
  rotation = c(0, 0, 180,0)
)
mixmod

set.seed(17)
data <- sim(mixmod, n = 2000)
ts.plot(data)

modspec <- dvinecopula(
  family = c("Gaussian", "Frank", "Clayton","t"),
  pars = list(0, 1, 0.5, c(0,10)),
  rotation = c(0, 0, 180,0)
)
fitmod <- fit(modspec, data, tsoptions = list(hessian = TRUE))
fitmod
coef(fitmod)
# plot fit
plot(fitmod, plotoption = 1)
plot(fitmod, plotoption = 2)
plot(fitmod, plotoption = 3)

data2 <- sim(fitmod, n =1000)
ts.plot(data2)

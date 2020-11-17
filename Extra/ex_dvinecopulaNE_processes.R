# Non-exchangeable model
mixmod3 <- dvinecopulaNE(
  family = c("Gaussian", "Frank", "Gumbel"), pars = list(0.3, c(1.1, 1), c(1.2, 1)),
  rotation = c(0, 0, 180), exchangeable = c(T, F, F)
)
mixmod3
y3 <- sim(mixmod3, n = 1000)
ts.plot(y3)

fit3 <- fit(mixmod3, y3, list(hessian = TRUE))
fit3
coef(fit3)

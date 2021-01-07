# create skew Student distribution
model <- margin("sst", pars = c(df = 4, gamma = 0.8, mu = -0.1, sigma = 2))
model

# simulate data from distribution
set.seed(13)
data <- sim(model, n = 5000)

# fit distribution to data
fit <- fit(model, data)
fit


sdweibull_mod <- margin("sdoubleweibull",
                        pars = c(mu = -0.5, shape = 0.8, scale = 1.5, gamma = 0.9))
dweibull_mod <- margin("doubleweibull")
laplace_mod <- margin("laplace")
slaplace_mod <- margin("slaplace")

sdweibull_mod
set.seed(133)
data <- sim(sdweibull_mod, 2000)
hist(data)

fit(laplace_mod, data)
fit(slaplace_mod, data)
fit(dweibull_mod, data)
fit(sdweibull_mod, data)

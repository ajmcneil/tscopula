# EXAMPLES WITH ARMACOPULA

# create vtarmacopula process specification
vtar1 <- vtscopula(armacopula(list(ar = 0.7)),
  Vtransform = Vlinear()
)
vtar1

# simulate a realisation
set.seed(19)
data <- sim(vtar1, 2000)
ts.plot(data)

# fit specification to data
tmp <- fit(vtar1, data, list(fulcrum = 0.5))
tmp

# create vtarmacopula process specification
vtarma11 <- vtscopula(armacopula(list(ar = 0.95, ma = -0.85)),
  Vtransform = Vlinear(delta = 0.55)
)
vtarma11

# simulate a realisation
set.seed(19)
data <- sim(vtarma11, 2000)
ts.plot(data)

# fit specification to data
tmp <- fit(vtarma11, data)
tmp

# plot fit
plot(tmp, plotoption = 1)
plot(tmp, plotoption = 2)
plot(tmp, plotoption = 3)
plot(tmp, plotoption = 4)
plot(tmp, plotoption = 5)
plot(tmp, plottype = "vtransform")


# create vtarmacopula process specification
vtarma11 <- vtscopula(armacopula(list(ar = 0.95, ma = -0.85)),
  Vtransform = V2p(delta = 0.55, kappa = 0.8)
)
vtarma11

# simulate a realisation
set.seed(17)
data <- sim(vtarma11, n = 5000)
ts.plot(data)

# fit specification to data
tmp <- fit(vtarma11, data)
tmp

# create vtarmacopula process specification
vtar2 <- vtscopula(armacopula(list(ar = c(0.6, 0.3))),
  Vtransform = V2p(delta = 0.55, kappa = 0.8)
)
vtar2

# simulate a realisation
set.seed(13)
data <- sim(vtar2, 2000)
ts.plot(data)

# fit specification to data
tmp <- fit(vtar2, data, list(fulcrum = 0.55))
tmp

# profile likelihood
profilefulcrum(data, tscopula = vtar2, locations = seq(from = 0, to = 1, length = 41))
abline(v = 0.55)

# EXAMPLES WITH DVINECOPULA

# create vtdvinecopula specification
mod <- vtscopula(dvinecopula(family = c("Joe", "Gauss", "Clayton"), pars = list(3, -.7, 3)),
  Vtransform = V2p(delta = 0.55, kappa = 1.2)
)
mod
set.seed(1234)
y <- sim(mod)
ts.plot(y)

# fixed fulcrum
fit_c <- fit(mod, y, list(fulcrum = 0.55))
fit_c

# free fulcrum
fit_c2 <- fit(fit_c, y)
fit_c2

# second example
mod <- vtscopula(dvinecopula(family = c("Joe", "t", "Clayton"), pars = list(2, c(0.4, 6), 3)),
  Vtransform = V2p(delta = 0.55, kappa = 1.2)
)
mod
set.seed(1234)
y <- sim(mod)

fit_c <- fit(mod, y, list(fulcrum = 0.55))
fit_c

fit_c2 <- fit(fit_c, y, control = list(trace = TRUE))
fit_c2
coef(fit_c2)

# EXAMPLES WITH DVINECOPULA2

copmod <- dvinecopula2(family = "joe",
                       pars = list(ar = 0.95, ma = -0.85),
                       maxlag = 30)
copmod2 <- dvinecopula2(family = "gauss",
                       pars = list(ar = 0.95, ma = -0.85),
                       maxlag = 30)
mod <- vtscopula(copmod,
                 Vtransform = V2p(delta = 0.55, kappa = 1.2)
)
mod2 <- vtscopula(copmod2,
                 Vtransform = V2p(delta = 0.55, kappa = 1.2)
)
mod
set.seed(1234)
y <- sim(mod)
ts.plot(y)

fit_c <- fit(mod, y, list(fulcrum = 0.55))
fit_c

fit_c2 <- fit(fit_c, y)
fit_c2
coef(fit_c2)

fit_alt <- fit(mod2, y)
fit_alt

# FULL MODELS WITH V-TRANSFORMS

fullmod <- tscm(mod, margin("slaplace", pars = c(mu = 1, scale = 2, gamma = 0.8)))
fullmod
set.seed(14)
data <- sim(fullmod, n= 1000)
hist(data)
tmp <- fit(fullmod, data)
tmp2 <- fit(tmp, data, method = "full")
tmp2

plot(tmp2, plotoption = 1)
plot(tmp2, plotoption = 2)
plot(tmp2, plotoption = 3)
plot(tmp2, plottype = "margin")
plot(tmp2, plotoption = 2, plottype = "margin")
plot(tmp2, plottype = "vtransform")
plot(tmp2, plottype = "volprofile")
plot(tmp2, plottype = "volproxy")
plot(tmp2, plotoption = 2, plottype = "volproxy")


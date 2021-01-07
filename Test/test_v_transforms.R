
V0 <- Vsymmetric() # symmetric V-transform
V1 <- Vlinear(0.6) # linear V-transform
V2 <- V2p(delta = 0.45, kappa = 1.5) # 2-parameter V-transform
V3 <- V3p(delta = 0.45, kappa = 1.5, xi = 1.2) # 3-parameter V-transform
V3
V2b <- V2b(delta = 0.45, kappa = 1.1) # 2-parameter beta V-transform
V3b <- V3b(delta = 0.5, kappa = 3, xi = 1.2) # 3-parameter beta V-transform

# some examples of plots
plot(V0)
plot(V1)
plot(V3)
plot(V3, shading = FALSE)
plot(V0, type = "inverse")
plot(V1, type = "inverse")
plot(V3, type = "inverse")
plot(V0, type = "gradient")
plot(V1, type = "gradient")
plot(V3, type = "gradient")
plot(V0, type = "pdown")
plot(V1, type = "pdown")
plot(V3, type = "pdown")

plot(V2b)
plot(V3b)
plot(V2b, type = "inverse")
plot(V3b, type = "inverse")
plot(V2b, type = "gradient")
plot(V3b, type = "gradient")
plot(V2b, type = "pdown")
plot(V3b, type = "pdown")

# stochastically inverting a v-transform
set.seed(13)
n <- 10000
U <- runif(n)
V <- vtrans(V3, U)
U2 <- stochinverse(V3, V)
plot(U, U2)
hist(U2, prob = TRUE)
abline(h = 1)
sum(round(U, 4) == round(U2, 4)) / n # how often do U and U2 coincide?

# Now compute the true probability of coincidence
pcoincide(V3)


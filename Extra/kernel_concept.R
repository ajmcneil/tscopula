library(tscopulaplus)
suppressMessages(library(quantmod))
BTCUSD <- suppressWarnings(getSymbols("BTCUSD=X", auto.assign = FALSE,
                                      from = "2015-12-31", to = "2019-12-31"))
BTCUSD <- na.omit(BTCUSD) # omit missing values
BTCUSD <- BTCUSD$`BTCUSD=X.Close` # closing values
plot(BTCUSD)

X <- (diff(log(BTCUSD))[-1]) * 100 # log-returns (as percentages)
length(X)
plot(X)

V <- strank(abs(X))
######


emp_ACF <- acf(qnorm(V))
wt <- drop(emp_ACF$acf)
kern_ACF <- ksmooth(1:30, wt[-1], bandwidth = 5, kernel = "normal", x.points = 1:30)
lines(kern_ACF$x, kern_ACF$y, col = "blue")
kern_PACF <- drop(.Call(stats:::C_pacf1, c(1,kern_ACF$y), lag.max = max(1:30)))[1:30]


set.seed(12345)
smooth_GAUSS <- qnorm(sim(dvinecopula(family = "gauss", pars = kern_PACF)))
acf(smooth_GAUSS)

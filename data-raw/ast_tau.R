## code to prepare `ast_tau` dataset goes here

load("ast_tau.RData")
tauspline <- smooth.spline(ast_tau$tau, ast_tau$nu)
usethis::use_data(tauspline, internal = TRUE)

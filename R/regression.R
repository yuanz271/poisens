# Regression ----
pois_reg <- function(x, y, maxit = 1000, tol = 1e-8) {
  reg <- list(coefficient = fit_pois_reg(x, y, maxit, tol), family = "poisson")
  class(reg) <- "reg"
  reg
}

negbin_reg <- function(x, y, phi, maxit = 1000, tol = 1e-8) {
  reg <- list(coefficient = fit_negbin_reg(x, y, phi, maxit, tol), family = "negbin", phi = phi)
  class(reg) <- "reg"
  reg
}

negbin_reg2 <- function(x, y, maxit = 1000, tol = 1e-8) {
  fit0 <- pois_reg(x, y, maxit, tol)
  mu0 <- predict(fit0, x, FALSE)
  ph0 <- phi_ml(y, mu0, maxit, tol)
  for (i in 1:maxit) {
    fit1 <- negbin_reg(x, y, ph0)
    mu1 <- predict(fit1, x, FALSE)
    ph1 <- phi_ml(y, mu1, maxit, tol)
    if (is.nan(ph1) || abs(ph1 - ph0) < tol)
      break
    ph0 <- ph1
  }
  fit1$phi <- ph1
  fit1$phi_se <- sqrt(1/attr(ph1, "info"))
  fit1
}

predict.reg <- function(obj, newdata, log = FALSE) {
  stopifnot(is.matrix(newdata))
  if (log) newdata %*% obj$coefficient else exp(newdata %*% obj$coefficient)
}

# Elastic net ----
# Design matrix x must have intercept column at first.
pois_lambda <- function(x, y, alpha = 1, k = 10) {
  ybar <- mean(y)
  lmax <- max(abs(t(x[, -1]) %*% (y - ybar))) / alpha
  lmin <- 0.0001 * lmax
  exp(seq(from = log(lmax), to = log(lmin), length.out = k))
}

pois_enet <- function(x, y, alpha = 1, lambda, intercept = TRUE, maxit = 10000, tol = 1e-8) {
  reg <- list(coefficient = fit_pois_enet(x, y, alpha, lambda, intercept, maxit, tol),
              family = "poisson", alpha = alpha, lambda = lambda)
  class(reg) <- c("reg", "enet")
  reg
}

negbin_lambda <- function(x, y, phi, alpha = 1, k = 10) {
  ybar <- mean(y)
  lmax <- max(abs(t(x[, -1]) %*% (y - ybar))) / phi / (alpha * (ybar + 1/phi))
  lmin <- 0.0001 * lmax
  exp(seq(from = log(lmax), to = log(lmin), length.out = k))
}

negbin_enet <- function(x, y, phi, alpha = 1, lambda, intercept = TRUE, maxit = 10000, tol = 1e-8) {
  reg <- list(coefficient = fit_negbin_enet(x, y, phi, alpha, lambda, intercept, maxit, tol),
              family = "negbin", phi = phi, alpha = alpha, lambda = lambda)
  class(reg) <- c("reg", "enet")
  reg
}

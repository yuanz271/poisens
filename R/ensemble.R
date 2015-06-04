# Ensemble of regressions ----
pois_reg_ens <- function(x, y, nreg = 50, nsample = nrow(x), maxit = 1000, tol = 1e-8) {
  index <- replicate(nreg, sample(nrow(x), size = nsample, replace = TRUE), simplify = FALSE)
  ens <- list(coefficient = foreach(i = index, .combine = cbind, .multicombine = TRUE, .errorhandling = "remove") %do% {
    pois_reg(x[i, ], y[i], maxit = maxit, tol = tol)$coefficient
  }, family = "poisson")
  class(ens) <- c("reg_ens")
  ens
}

negbin_reg_ens <- function(x, y, phi, nreg = 50, nsample = nrow(x), maxit = 1000, tol = 1e-8) {
  index <- replicate(nreg, sample(nrow(x), size = nsample, replace = TRUE), simplify = FALSE)
  ens <- list(coefficient = foreach(i = index, .combine = cbind, .multicombine = TRUE, .errorhandling = "remove") %do% {
    negbin_reg(x[i, ], y[i], phi = phi, maxit = maxit, tol = tol)$coefficient
  }, family = "negbin", phi = phi)
  class(ens) <- c("reg_ens")
  ens
}

negbin_reg_ens1 <- function(x, y, nreg = 50, nsample = nrow(x), maxit = 1000, tol = 1e-8) {
  index <- replicate(nreg, sample(nrow(x), size = nsample, replace = TRUE), simplify = FALSE)
  fits <- foreach(i = index, .multicombine = TRUE, .errorhandling = "remove") %do% {
    negbin_reg2(x[i, ], y[i], maxit = maxit, tol = tol)
  }
  coef <- foreach(fit = fits, .combine = cbind, .multicombine = TRUE) %do% {
    fit$coefficient
  }
  phi <- foreach(fit = fits, .combine = c, .multicombine = TRUE) %do% {
    fit$phi
  }
  info <- foreach(fit = fits, .combine = c, .multicombine = TRUE) %do% {
    attr(fit$phi, "info")
  }
  ens <- list(coefficient = coef, phi = phi, info = info, family = "negbin")
  class(ens) <- c("reg_ens")
  ens
}

negbin_reg_ens2 <- function(x, y, nreg = 50, nsample = nrow(x), maxit = 1000, tol = 1e-8) {
  fit0 <- pois_reg(x, y, maxit, tol)
  mu0 <- predict(fit0, x, FALSE)
  ph0 <- phi_ml(y, mu0, maxit, tol)
  for (i in 1:maxit) {
    fit1 <- negbin_reg_ens(x, y, ph0, nreg, nsample)
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

predict.reg_ens <- function(obj, newdata, log = FALSE) {
  stopifnot(is.matrix(newdata))
  as.matrix(if (log) rowMeans(newdata %*% obj$coefficient) else rowMeans(exp(newdata %*% obj$coefficient)))
}

# Ensemble of enet ----
pois_enet_ens <- function(x, y, alpha = 1, lambda, intercept = TRUE, nreg = 50, nsample = nrow(x), maxit = 10000, tol = 1e-8) {
  index <- replicate(nreg, sample(nrow(x), size = nsample, replace = TRUE), simplify = FALSE)
  ens <- list(coefficient = foreach(i = index, .multicombine = TRUE, .errorhandling = "remove") %do% {
    pois_enet(x[i, ], y[i], alpha, lambda, intercept, maxit = maxit, tol = tol)$coefficient
  }, family = "poisson")
  class(ens) <- c("enet_ens")
  ens
}

negbin_enet_ens <- function(x, y, phi, alpha = 1, lambda, intercept = TRUE, nreg = 50, nsample = nrow(x), maxit = 10000, tol = 1e-8) {
  index <- replicate(nreg, sample(nrow(x), size = nsample, replace = TRUE), simplify = FALSE)
  ens <- list(coefficient = foreach(i = index, .multicombine = TRUE, .errorhandling = "remove") %do% {
    negbin_enet(x[i, ], y[i], phi, alpha, lambda, intercept, maxit = maxit, tol = tol)$coefficient
  }, family = "negbin", phi = phi)
  class(ens) <- c("enet_ens")
  ens
}

predict.enet_ens <- function(obj, newdata, log = FALSE) {
  stopifnot(is.matrix(newdata))
  Reduce(`+`, lapply(obj$coefficient, function(coef) if (log) newdata %*% coef else exp(newdata %*% coef))) / length(obj$coefficient)
}



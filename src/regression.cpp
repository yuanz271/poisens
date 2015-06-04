// [[Rcpp::depends(RcppArmadillo)]]
#include <limits>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "distribution.h"

using namespace arma;

template<class T>
arma::vec irls(const arma::mat& x, const arma::vec& y, T dist,
                    const unsigned int maxit, const double tol) {
  vec old_beta(x.n_cols, fill::zeros), beta(x.n_cols, fill::zeros), mu(x.n_rows);
  mu.fill(mean(y));
  vec z = dist.wresponse(y, mu);
  mat w = diagmat(ones<vec>(x.n_rows));
  for(unsigned int i = 0; i < maxit; ++i) {
    beta = solve(x.t() * w * x, x.t() * w * z);
    mu = exp(x * beta);
    w = diagmat(dist.weight(mu));
    z = dist.wresponse(y, mu);
    if (norm(old_beta - beta) < tol)
      break;
    old_beta = beta;
    if (i + 1 == maxit)
      Rcpp::stop("max iteration");
  }
  return beta;
}

// [[Rcpp::export]]
arma::vec fit_pois_reg(const arma::mat& x, const arma::vec& y,
                     const unsigned int maxit, const double tol) {
  return irls(x, y, Poisson(), maxit, tol);
}

// [[Rcpp::export]]
arma::vec fit_negbin_reg(const arma::mat& x, const arma::vec& y, const double phi,
                     const unsigned int maxit, const double tol) {
  return irls(x, y, Negbin(phi), maxit, tol);
}


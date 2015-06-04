// [[Rcpp::depends(RcppArmadillo)]]
#include <limits>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "distribution.h"

using namespace arma;

double soft(const double x, const double r) {
  if (x > r)
    return x - r;
  else if (x < -r)
    return x + r;
  else
    return 0.0;
}

template<class T>
vec coord_desc(const mat& x, const vec& y, const vec& b, T dist,
               const double alpha, const double lambda, const bool intercept,
               const unsigned int maxit, const double tol) {
  vec old_beta = b;
  vec etaj(x.n_rows);
  vec mu = exp(x * b);
  vec w = dist.weight(mu);
  vec z = dist.wresponse(y, mu);
  vec beta(x.n_cols, fill::zeros);
  vec xwx = diagvec(x.t() * diagmat(w) * x);

  for(unsigned int i = 0; i < maxit; ++i) {
    beta[0] = intercept ? accu(w % (z - x * beta + beta[0])) / xwx[0] : 0.0;
    for(unsigned int j = 1; j < x.n_cols; ++j) {
      etaj = x * beta - x.col(j) * beta[j];
      beta[j] = soft(accu(x.col(j) % w % (z - etaj)), lambda * alpha) / (xwx[j] + lambda * (1.0 - alpha));
    }
    if (norm(beta - old_beta) < tol)
      break;
    old_beta = beta;
    if (i + 1 == maxit)
      Rcpp::stop("not converged");
  }
  return beta;
}

template<class T>
mat fit_enet(const mat& x, const vec& y, T dist,
             const double alpha, const vec& lambda, const bool intercept,
             const unsigned int maxit, const double tol) {
  mat betas(x.n_cols, lambda.size());
  vec old_beta(x.n_cols, fill::zeros), beta;
  old_beta[0] = std::log(mean(y));

  for(unsigned int k = 0; k < lambda.size(); ++k) {
    for(unsigned int i = 0; i < maxit; ++i) {
      beta = coord_desc(x, y, old_beta, dist, alpha, lambda[k], intercept, maxit, tol);
      if (norm(beta - old_beta) < tol)
        break;
      old_beta = beta;
      if (i + 1 == maxit)
        Rcpp::stop("max iteration");
    }
    betas.col(k) = beta;
  }

  return betas;
}

// [[Rcpp::export]]
arma::mat fit_pois_enet(const arma::mat& x, const arma::vec& y,
              const double alpha, const arma::vec& lambda, const bool intercept,
              const unsigned int maxit, const double tol) {
  return fit_enet(x, y, Poisson(), alpha, lambda, intercept, maxit, tol);
}

// [[Rcpp::export]]
arma::mat fit_negbin_enet(const arma::mat& x, const arma::vec& y, const double phi,
              const double alpha, const arma::vec& lambda, const bool intercept,
              const unsigned int maxit, const double tol) {
  return fit_enet(x, y, Negbin(phi), alpha, lambda, intercept, maxit, tol);
}

#include <limits>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include "distribution.h"

using namespace arma;

double Poisson::likelihood(const vec& y, const vec& mu) {
  double lik = 0.0;
  for(unsigned int i = 0; i < y.n_elem; ++i) {
    lik += std::log(R::dpois(y[i], mu[i], 0));
  }
  return lik;
}

vec Poisson::weight(const vec& mu) {
  return mu;
}

vec Poisson::wresponse(const vec& y, const vec& mu) {
  return log(mu) + (y - mu) / mu;
}

void Poisson::print() {
  Rcpp::Rcout << "Poisson" << std::endl;
}


Negbin::Negbin(const double a) : phi(a) {}

double Negbin::likelihood(const vec& y, const vec& mu) {
  double lik = 0.0;
  if (std::abs(phi) < std::numeric_limits<double>::epsilon()) {
    for(unsigned int i = 0; i < y.n_elem; ++i) {
      lik += std::log(R::dnbinom_mu(y[i], 1 / phi, mu[i], 0));
    }
  } else {
    for(unsigned int i = 0; i < y.n_elem; ++i) {
      lik += std::log(R::dpois(y[i], mu[i], 0));
    }
  }
  return lik;
}

vec Negbin::weight(const vec& mu) {
  return mu / (1 + mu * phi);
}

vec Negbin::wresponse(const vec& y, const vec& mu) {
  return log(mu) + (y - mu) / mu;
}

void Negbin::print() {
  Rcpp::Rcout << phi << std::endl;
}

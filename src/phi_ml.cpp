#include <Rcpp.h>
#include <cmath>
#include <Rmath.h>

using namespace Rcpp;
using namespace R;

double score(const NumericVector& y, const NumericVector& mu, const double phi) {
  return - sum(digamma(y + 1/phi) - digamma(1/phi) - std::log(phi) - log(mu + 1/phi) - (y - mu) / (mu + 1/phi)) / std::pow(phi, 2);
}

double info(const NumericVector& y, const NumericVector& mu, const double phi) {
  return - sum(trigamma(y + 1/phi) - trigamma(1/phi) + phi - 2 / (1/phi + mu) + (y + 1/phi) / pow(mu + 1/phi, 2)) / std::pow(phi, 4);
}

double fisher(const NumericVector& mu, const double phi) {
  return - sum(trigamma(mu + 1/phi) - trigamma(1/phi) + phi - 1 / (1/phi + mu)) / std::pow(phi, 4);
}

// [[Rcpp::export]]
NumericVector phi_ml(const NumericVector& y, const NumericVector& mu, const unsigned int maxit, const double tol) {
  double delta = 0;
  double info;
  double ph = sum(pow(y / mu - 1, 2)) / y.size();

  for (unsigned int i = 0; i < maxit; ++i) {
    ph = std::abs(ph); // Rcpp::abs always return 1 for scalar argument. Use std::abs instead
    info =  fisher(mu, ph);
    delta = - score(y, mu, ph) / info;
    ph = ph - delta;
    if (std::abs(delta) < tol)
      break;
    if (i + 1 == maxit) {
      Function warning("warning");
      warning("max iteration");
    }
  }
  if (ph < 0.0)
    ph = 0.0;
  NumericVector phi = NumericVector::create(ph);
  phi.attr("info") = info + std::numeric_limits<double>::epsilon();
  return phi;
}


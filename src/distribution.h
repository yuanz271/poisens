template<class T>
class Distribution {
private:
  T& instance() {
    return static_cast<T&>(*this);
  }

public:
  double likelihood(const arma::vec& y, const arma::vec& mu) {
    return instance().likelihood(y, mu);
  }

  arma::vec weight(const arma::vec& mu) {
    return instance().weight(mu);
  }

  arma::vec wresponse(const arma::vec& y, const arma::vec& mu) {
    return instance().wresponse(y, mu);
  }

  void print() {
    instance().print();
  }

};

class Poisson : public Distribution<Poisson> {
public:
  double likelihood(const arma::vec& y, const arma::vec& mu);

  arma::vec weight(const arma::vec& mu);

  arma::vec wresponse(const arma::vec& y, const arma::vec& mu);

  void print();

};

class Negbin : public Distribution<Negbin> {
private:
  double phi;

public:
  Negbin();

  Negbin(const double);

  double likelihood(const arma::vec& y, const arma::vec& mu);

  arma::vec weight(const arma::vec& mu);

  arma::vec wresponse(const arma::vec& y, const arma::vec& mu);

  void print();

};

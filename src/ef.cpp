#define TMB_LIB_INIT R_init_ef
#include <TMB.hpp>

// log likelihood of a single multipass electrofishing sample
template<class Type>
Type ll_sample(vector<Type> p, vector<Type> y) {

  int npasses = y.size();
  vector<Type> pi(npasses);
  Type ll = Type(0);

  for (int i = 0; i < npasses; i++) {
    pi(i) = p(i);
    for (int j = 0; j < i; j++) {
      pi(i) *= 1. - p(j);
    }
    if (y(i) > 0) {
      ll += y(i) * log(pi(i));
    }
  }
  ll -= y.sum() * log(pi.sum());

  return ll;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(offset);
  DATA_FACTOR(sample_id);
  DATA_SPARSE_MATRIX(X);
  DATA_SPARSE_MATRIX(Z);
  DATA_INTEGER(sample_re);

  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);
  PARAMETER(log_sigma);

  int n_samples = sample_id.maxCoeff() + 1;

  vector<Type> eta = X * beta + offset;
  if (sample_re == 1) {
    // add sample random effect
    eta += Z * alpha;
  }

  vector<Type> p = invlogit(eta);
  vector<vector<Type>> ps = split(p, sample_id);
  vector<vector<Type>> ys = split(y, sample_id);

  Type nll = Type(0); // negative log likelihood

  if (sample_re == 1) {
    // add sample random effect
    nll -= dnorm(alpha, Type(0), exp(log_sigma), true).sum();
  }

  for (int i = 0; i < n_samples; i++) {
    nll += -ll_sample(ps(i), ys(i));
  }

  // delta method for sigma
  Type sigma = exp(log_sigma);
  ADREPORT(sigma);
  REPORT(sigma);

  REPORT(ps);

  return nll;
}

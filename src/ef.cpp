#include <TMB.hpp>

template<class Type>
Type nll_sample(vector<Type> p, vector<Type> y) {

  int npasses = y.size();
  vector<Type> pi(npasses);
  Type nll = 0;

  for (int i = 0; i < npasses; i++) {
    pi(i) = p(i);
    for (int j = 0; j < i; j++) {
      pi(i) *= 1. - p(j);
    }
    if (y(i) > 0) {
      nll += y(i) * log(pi(i));
    }
  }
  nll -= y.sum() * log(pi.sum());

  return nll;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(offset);
  DATA_FACTOR(sample_id);
  DATA_MATRIX(A);

  PARAMETER_VECTOR(alpha);

  int n_samples = sample_id.maxCoeff() + 1;

  // calculate all the within pass probs required
  vector<Type> eta = A * alpha + offset;
  vector<Type> p = invlogit(eta);
  vector<vector<Type>> ps = split(p, sample_id);
  vector<vector<Type>> ys = split(y, sample_id);

  Type nll = 0.0; // negative log likelihood

  for (int i = 0; i < n_samples; i++) {
    nll += nll_sample(ps(i), ys(i));
  }

  return nll;
}

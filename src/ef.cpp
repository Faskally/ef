#include <TMB.hpp>

// log likelihood of a single multipass electrofishing sample
template<class Type>
Type ll_sample(vector<Type> p, vector<Type> y) {

  int npasses = y.size();
  vector<Type> pi(npasses);
  Type ll = 0;

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
  DATA_MATRIX(A);

  PARAMETER_VECTOR(alpha);

  int n_samples = sample_id.maxCoeff() + 1;

  // calculate all the within pass probs required
  vector<Type> eta = A * alpha + offset;
  vector<Type> p = invlogit(eta);
  vector<vector<Type>> ps = split(p, sample_id);
  REPORT(ps)
  vector<vector<Type>> ys = split(y, sample_id);
  REPORT(ys)

  Type nll = 0.0; // negative log likelihood

  for (int i = 0; i < n_samples; i++) {
    nll -= ll_sample(ps(i), ys(i));
  }

  return nll;
}

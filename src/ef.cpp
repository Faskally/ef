#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER(N); // number of multipass observations
  DATA_INTEGER(K); // number of parameters
  DATA_INTEGER(S); // max number of passes
  DATA_IVECTOR(npasses); // number of passes per observation - vectors of S passes
  DATA_IMATRIX(y); // data - N vectors of S passes
  DATA_IVECTOR(yT); // data - N total catches
  DATA_IMATRIX(offset); // N offsets - S passes
  DATA_IARRAY(A); // the design matrices (N x K) - one for each of s passes

  PARAMETER_VECTOR(alpha); // K alpha params

  array<Type> p(N, S); // N x S
  array<Type> pi(N, S); // N x S
  vector<Type> piT(N); // N

  Type nll = 0.0; // negative log likelihood

  // loop over observations
  for (int i = 0; i < N; i++) {

    // calculate all the within pass probs required
    for (int j = 0; j < S; j++) {
      for (int k = 1; k < K; k++) {
        p(i,j) = 1.0 / (1.0 + exp(-1.0 * A(j,i,k) * alpha(k) - offset(j,i)));
      }
    }

    // calculate all the marginal(?) probs required
    pi(i,1) = p(i,1);
    for (int j = 1; j < S; j++) {
      pi(i,j) = p(i,j);
      for (int jj = 0; jj < (j + 1); jj++) {
        pi(i,j) *= 1 - p(i,jj);
      }
    }

    // calculate the probability of capture - needs number of passes
    piT(i) = pi(i,1);
    for (int j = 1; j < npasses(i); j++) {
        piT(i) += pi(i,j);
    }

    // calculate the likelihood of the data
    for (int j = 1; j < S; j++) {
      nll += y(j,i) * log(pi(i,j));
    }
    nll += -1.0 * yT(i) * log(piT(i));
  }

  ADREPORT(alpha);

  return nll;
}

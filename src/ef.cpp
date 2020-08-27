#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER(N); // number of multipass observations
  DATA_INTEGER(K); // number of parameters
  DATA_INTEGER(s); // max number of passes
  DATA_IVECTOR(npasses); // number of passes per observation - vectors of s passes
  DATA_IMATRIX(y); // data - N vectors of s passes
  DATA_IVECTOR(yT); // data - N total catches
  DATA_IMATRIX(s_offset); // N offsets - s passes
  DATA_IARRAY(A); // the design matrices (N x K) - one for each of s passes

  PARAMETER_VECTOR(alpha); // K alpha params

  PARAMATER_MATRIX(p); // s x N
  PARAMATER_MATRIX(pi); // s x N
  PARAMETER_VECTOR(piT); // N

  // calculate all the within pass probs required
  for (int pass = 0; pass < s; pass++) {
    for (int i = 0; i < N; i++) {
      p[pass,i] = 0;
      for (int j = 0; j < K; j++) {
        p[pass,i] += 1.0 / (1.0 + exp(-1.0 * A[pass,i,j] * alpha[j] - s_offset[pass,i]));
      }
    }
  }

  // calculate all the marginal(?) probs required
  for (int i = 0; i < N; i++) {
    pi[1,i] = p[1,i];
  }
  for (int passi = 1; passi < s; passi++) {
    for (int i = 0; i < N; i++) {
      pi[passi,i] = p[passi,i];
    }
    for (passj in 1:(passi-1)) {
      for (int i = 0; i < N; i++) {
        pi[passi,i] = pi[passi,i] .* (1-p[passj,i]);
      }
    }
  }

  // calculate the probability of capture - needs number of passes
  for (int i = 0; i < N; i++) {
    piT[i] = pi[1,i];
  }
  for (i in 1:N) {
    for (int passi = 0; passi < npasses[i]; passi++) {
      piT[i] = piT[i] + pi[passi,i];
    }
  }

  Type nll = 0.0;

  // calculate the probability of capture - needs number of passes
  for (int pass = 0; pass < s; pass++) {
    for (int i = 0; i < N; i++) {
      nll += y[pass,i] * log(pi[pass,i]);
    }
  }
  for (int i = 0; i < N; i++) {
    nll += -1.0 * yT[i] * log(piT[i]);
  }


  REPORT(alpha);

  REPORT(p);
  REPORT(pi);
  REPORT(piT);

  ADREPORT(alpha);

  return nll;
}

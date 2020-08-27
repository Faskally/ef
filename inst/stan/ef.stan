data {
  int<lower=0> N; // number of multipass observations
  int<lower=0> K; // number of parameters
  int<lower=0> s; // max number of passes
  int<lower=0> npasses[N]; // number of passes per observation - vectors of s passes
  row_vector[N] y[s]; // data - vectors of s passes
  row_vector[N] yT; // data - total catch
  vector[N] s_offset[s]; // offset - s passes
  matrix[N,K] A[s]; // the design matrices - one for each pass
}
parameters {
    vector[K] alpha;
}
transformed parameters {
  vector[N] p[s];
  vector[N] pi[s];
  vector[N] piT;
  // calculate all the within pass probs required
  for (i in 1:s) {
    p[i] = 1.0 ./ (1.0 + exp(-1.0 * A[i] * alpha - s_offset[i]));
  }
  // calculate all the marginal(?) probs required
  pi[1] = p[1];
  for (i in 2:s) {
    pi[i] = p[i];
    for (j in 1:(i-1)) {
      pi[i] = pi[i] .* (1-p[j]);
    }
  }
  // calculate the probability of capture - needs number of passes
  piT = pi[1];
  for (i in 1:N) {
    for (j in 2:npasses[i]) {
      piT[i] = piT[i] + pi[j,i];
    }
  }
}
model {
  // calculate the probability of capture - needs number of passes
  for (i in 1:s) {
    target += y[i] * log(pi[i]);
  }
  target += -1.0 * yT * log(piT);
}

functions {
  vector mfrmr_logits(real eta, row_vector steps) {
    int C = num_elements(steps) + 1;
    vector[C] logits;
    logits[1] = 0;
    for (k in 2:C)
      logits[k] = logits[k - 1] + eta - steps[k - 1];
    return logits;
  }
}
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=2> C;
  int<lower=0> P;
  int<lower=1> L;
  int<lower=1> G;
  int<lower=0> Q;
  array[N] int<lower=1, upper=J> person;
  array[N] int<lower=1, upper=G> owner;
  array[N] int<lower=1, upper=C> y;
  matrix[N, P] X;
  matrix[L, P] location_basis;
  matrix[C - 1, Q] step_basis;
  vector<lower=0>[P] location_sd;
  real<lower=0> step_sd;
}
parameters {
  vector[J] theta;
  vector[P] location_free;
  matrix[G, Q] step_free;
}
transformed parameters {
  vector[L] locations = location_basis * location_free;
  matrix[G, C - 1] steps = step_free * step_basis';
  vector[N] eta = theta[person] + X * location_free;
}
model {
  theta ~ std_normal();
  location_free ~ normal(0, location_sd);
  to_vector(step_free) ~ normal(0, step_sd);
  for (n in 1:N)
    y[n] ~ categorical_logit(mfrmr_logits(eta[n], steps[owner[n]]));
}
generated quantities {
  vector[N] log_lik;
  array[N] int y_rep;
  for (n in 1:N) {
    vector[C] logits = mfrmr_logits(eta[n], steps[owner[n]]);
    log_lik[n] = categorical_logit_lpmf(y[n] | logits);
    y_rep[n] = categorical_logit_rng(logits);
  }
}

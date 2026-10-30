functions {
  vector mfrmr_logits(real eta, row_vector steps) {
    int C = num_elements(steps) + 1;
    vector[C] logits;
    logits[1] = 0;
    for (k in 2:C) logits[k] = logits[k - 1] + eta - steps[k - 1];
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
  int<lower=0> T;
  int<lower=0> S;
  int<lower=0> A;
  array[N] int<lower=1, upper=J> person;
  array[N] int<lower=1, upper=G> owner;
  array[N] int<lower=1, upper=C> y;
  matrix[N, P] X;
  matrix[L, P] location_basis;
  vector[L] location_offset;
  vector[N] eta_offset;
  matrix[G * (C - 1), T] step_basis;
  vector[G * (C - 1)] step_offset;
  vector<lower=0>[P] location_sd;
  vector<lower=0>[T] step_sd;
  matrix[A, S] slope_basis;
  matrix[N, S] slope_X;
  vector<lower=0>[S] log_slope_sd;
  vector[2] population_mean_prior;
  vector[2] population_log_sd_prior;
}
parameters {
  vector[J] theta_raw;
  real population_mean;
  real population_log_sd;
  vector[P] location_free;
  vector[T] step_free;
  vector[S] log_slope_free;
}
transformed parameters {
  real<lower=0> population_sd = exp(population_log_sd);
  vector[J] theta = population_mean + population_sd * theta_raw;
  vector[L] locations = location_offset + location_basis * location_free;
  matrix[G, C - 1] steps = to_matrix(step_offset + step_basis * step_free, C - 1, G)';
  vector[N] eta = theta[person] + eta_offset + X * location_free;
  vector[A] log_component_slopes = slope_basis * log_slope_free;
  vector[A] component_slopes = exp(log_component_slopes);
  vector[N] effective_slopes = exp(slope_X * log_slope_free);
}
model {
  theta_raw ~ std_normal();
  population_mean ~ normal(population_mean_prior[1], population_mean_prior[2]);
  population_log_sd ~ normal(population_log_sd_prior[1], population_log_sd_prior[2]);
  location_free ~ normal(0, location_sd);
  step_free ~ normal(0, step_sd);
  log_slope_free ~ normal(0, log_slope_sd);
  for (n in 1:N)
    y[n] ~ categorical_logit(effective_slopes[n] * mfrmr_logits(eta[n], steps[owner[n]]));
}
generated quantities {
  vector[N] log_lik;
  array[N] int y_rep;
  for (n in 1:N) {
    vector[C] logits = effective_slopes[n] * mfrmr_logits(eta[n], steps[owner[n]]);
    log_lik[n] = categorical_logit_lpmf(y[n] | logits);
    y_rep[n] = categorical_logit_rng(logits);
  }
}

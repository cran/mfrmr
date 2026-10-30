// Joint rating-scale imputation example: known N(0,1) Person population.
// Each missing response uses the same Person draw as that person's observed
// responses. Unassigned events are absent from both observed and missing arrays.
functions {
  vector rsm_logits(real eta, vector step) {
    int K = num_elements(step) + 1;
    vector[K] logits;
    logits[1] = 0;
    for (k in 2:K) logits[k] = logits[k - 1] + eta - step[k - 1];
    return logits;
  }
}
data {
  int<lower=1> P;
  int<lower=2> R;
  int<lower=2> C;
  int<lower=3> K;
  int<lower=1> N_obs;
  int<lower=1> N_mis;
  array[N_obs] int<lower=1, upper=P> person_obs;
  array[N_obs] int<lower=1, upper=R> rater_obs;
  array[N_obs] int<lower=1, upper=C> criterion_obs;
  array[N_obs] int<lower=1, upper=K> y_obs;
  array[N_mis] int<lower=1, upper=P> person_mis;
  array[N_mis] int<lower=1, upper=R> rater_mis;
  array[N_mis] int<lower=1, upper=C> criterion_mis;
  matrix[R, R - 1] Q_rater;
  matrix[C, C - 1] Q_criterion;
  matrix[K - 1, K - 2] Q_step;
  real<lower=0> prior_sd;
}
parameters {
  vector[P] theta;
  vector[R - 1] rater_free;
  vector[C - 1] criterion_free;
  vector[K - 2] step_free;
}
transformed parameters {
  vector[R] severity = Q_rater * rater_free;
  vector[C] difficulty = Q_criterion * criterion_free;
  vector[K - 1] step = Q_step * step_free;
}
model {
  theta ~ std_normal();
  rater_free ~ normal(0, prior_sd);
  criterion_free ~ normal(0, prior_sd);
  step_free ~ normal(0, prior_sd);
  for (n in 1:N_obs)
    y_obs[n] ~ categorical_logit(rsm_logits(
      theta[person_obs[n]] - severity[rater_obs[n]] - difficulty[criterion_obs[n]], step));
}
generated quantities {
  array[N_mis] int y_mis;
  matrix[N_mis, K] probability_mis;
  real log_likelihood = 0;
  for (n in 1:N_obs)
    log_likelihood += categorical_logit_lpmf(y_obs[n] | rsm_logits(
      theta[person_obs[n]] - severity[rater_obs[n]] - difficulty[criterion_obs[n]], step));
  for (n in 1:N_mis) {
    vector[K] logits = rsm_logits(
      theta[person_mis[n]] - severity[rater_mis[n]] - difficulty[criterion_mis[n]], step);
    probability_mis[n] = softmax(logits)';
    y_mis[n] = categorical_logit_rng(logits);
  }
}

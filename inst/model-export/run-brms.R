# Run Rscript run.R after reviewing README.txt and priors.
# Run Rscript run.R --summarize to review fit.rds without compiling or sampling.
summarize_mfrm_brms <- function(fit) {
  source("model.R", local = TRUE)
  source("summarize.R", local = TRUE)
  parameters <- extract_mfrm_parameters(fit)
  saveRDS(parameters, "parameters.rds")
  summarize_mfrm_draws(mfrm_brms_draws(fit, parameters), bundle, mfrm_brms_sampler(fit), "brms")
}

run_mfrm_brms <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) && !identical(args, "--summarize")) stop("Usage: Rscript run.R [--summarize]")
  if (!requireNamespace("brms", quietly = TRUE)) stop("Install brms first.")
  if (identical(args, "--summarize")) {
    fit <- readRDS("fit.rds")
  } else {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("Install CmdStanR and CmdStan first.")
    source("model.R", local = TRUE)
    fit <- brms::brm(formula, data = model_data, family = family, prior = prior,
      stanvars = stanvars, backend = "cmdstanr", seed = 20261010,
      chains = 4, cores = 4, iter = 2000, warmup = 1000,
      control = list(adapt_delta = 0.95))
    # Preserve the fit before any reporting or predictive processing.
    saveRDS(fit, "fit.rds")
    saveRDS(list(y_rep = brms::posterior_predict(fit), log_lik = brms::log_lik(fit)),
      "predictive-draws.rds")
  }
  summarize_mfrm_brms(fit)
}

if (sys.nframe() == 0L) run_mfrm_brms()

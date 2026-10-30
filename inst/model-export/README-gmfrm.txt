Bayesian two-family GMFRM model exported by mfrmr

Run Rscript run.R from this directory after installing CmdStanR and CmdStan.
Export itself requires neither. The script requests four chains and retains
diagnostics. Review its resource settings and inspect divergences, R-hat,
effective sample sizes and posterior predictive fit before interpreting draws.

The summary screen and summary-decision.csv show the sampling review before
numerical tables. sampling-review.csv records each check, its criterion and
next action. parameter-summary.csv joins original labels to abilities, facet
locations, step offsets and slopes, with reference/direction notes. Its
value_role distinguishes estimated parameters, fixed anchors and values
determined by constraints. Fixed values have no estimated SD or interval
(NA); zero calibration uncertainty is not inferred. Internal coordinates
remain available for diagnostic review. posterior-summary.rds retains all
tables and interpretation notes.

The screen checks at least four chains, R-hat < 1.01, bulk/tail ESS >= 100
per chain, zero divergences/tree-depth hits and E-BFMI >= 0.3. These are
screening criteria, not a certificate of inference validity. Tree-depth hits
flag efficiency separately from divergences. Missing diagnostics require
review. Finite intervals remain visible when sampling needs review, with
interval_status beside each row. No flag does not establish model fit, prior
adequacy or calibrated coverage. Inspect endpoint MCSE for the intended use;
no universal acceptable fraction of interval width is imposed.

This Stan bundle supports Rscript run.R --summarize to regenerate the five report
files from fit.rds without compiling or sampling. Run it from this directory;
keep model-data.rds and summarize.R with the saved fit.
CmdStanR/posterior are needed for replay, but a CmdStan toolchain is not. Draws are unchanged.

The response likelihood retains the source specification:
  log(P(Y=k+1) / P(Y=k)) = a_first * a_second *
    (theta + signed facet locations - step[k]).
The product of slopes multiplies the whole adjacent-category predictor.
This is the logistic metric (L-metric), with scaling constant D = 1; no extra
1.7 or 1.702 factor is introduced. A one-unit ability increase changes each
adjacent-category log odds by a_first * a_second, holding other values fixed.
This response-link convention is distinct from the ability prior.
Abilities have independent normal(0,1) priors;
realized abilities are not centered across Persons. The first owner's log
slopes sum to zero, so its slopes have geometric mean one. Second-owner log
slopes are free. Locations retain their source signs and centering; each step
owner's offsets sum to zero and need not be ordered.

model-data.rds contains the design matrices, priors and mappings. levels.csv
maps Person and location labels; slope-levels.csv maps component_slopes and
log_component_slopes to owners and levels. Category 1 is the lowest source
score; categories.csv retains unused declared categories. steps[g,k] is the
offset for transition k to k+1; model-data.rds$step_owners identifies g.

The default marginal prior SD is 1 for locations, step offsets and component
log slopes. Centered blocks use orthonormal coordinates with SD multiplied by
sqrt(m/(m-1)); their constrained components have the declared marginal SD and
are dependent. Uncentered blocks use independent normal coordinates. These
are explicit proper Bayesian priors, not MML penalties or an empirical
coverage correction. Change location_sd, step_sd or log_slope_sd in
model-data.rds$stan_data for prior sensitivity; update $priors metadata too.

The model follows the product-slope response equation in Uto and Ueno (2020),
equation 9, DOI 10.1007/s41237-020-00115-7, while preserving the fitted mfrmr
constraints. It does not reproduce the article's Appendix program: that code
includes a 1.7 multiplier and applies densities to constrained components in
a different parameterization. The priors here are defined directly on free
normal coordinates; no lognormal density or Jacobian for an exponential
parameter change is needed. The exponentials are derived quantities.

Posterior credible intervals do not reproduce or validate MML confidence
intervals. This export does not establish sampling coverage, global
identification, convergence or the adequacy of the selected priors. Log
likelihood and replicated ratings condition on sampled abilities of existing
Persons; they are not new-Person predictions or leave-one-Person-out scores.
Repeated ratings retain the source conditional-independence assumption.

This format supports native, unanchored, unit-weight two-family GPCM MML
sources with fixed N(0,1), second-owner steps, and no population regression,
interactions, dummy facets or shrinkage. One-family GPCM uses a separate
affine template in Stan or brms; two-family brms export remains unsupported. Retained ratings are exported, not omitted input rows or
the original MML optimization. The bundle includes ratings and original labels;
review it before sharing.

Official documentation: https://mc-stan.org/cmdstanr/
https://mc-stan.org/learn-stan/diagnostics-warnings.html
https://mc-stan.org/posterior/reference/mcse_quantile.html

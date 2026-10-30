Bayesian RSM/PCM/GPCM model exported by mfrmr

Run Rscript run.R from this directory after installing CmdStanR and CmdStan;
the brms format also needs brms. Export itself runs no estimator. Review the
priors and runner settings before sampling. The bundle includes observed
ratings and original labels; review it before sharing.

The likelihood is adjacent-category logistic (D = 1):
  log(P(Y=k+1) / P(Y=k)) = effective_slope *
    (theta + signed facet locations - step[k]).
RSM/PCM use effective_slope = 1. A one-family GPCM uses its owner's slope,
whose log slopes sum to zero (geometric mean one). For two-family GPCM the
product of both slopes multiplies the entire predictor; only the first
family's log slopes sum to zero. Two-family GPCM is Stan-only.

Check model-data.rds$population to identify the population model. If absent,
abilities are independent normal(0,1). For estimated-normal one-family Stan
bundles, theta = population_mean + population_sd * theta_raw, where theta_raw
is standard normal. The mean and SD are sampled, not fixed MML estimates.
Realized Persons are not centered in either model. The default priors are
population_mean ~ normal(0,1) and log(population_sd) ~ normal(0,0.5), giving
the SD a median of 1 and central 95% prior interval about 0.38--2.66.
Choose these assumptions using population_prior when exporting and inspect
prior predictions and sensitivity. The two stan_data population_*_prior
vectors store c(location, scale); model-data.rds also records their meaning.
levels.csv reports the marginal Person prior SD, including mean/SD uncertainty.
Only the intercept-only estimated population is supported, and only in Stan.
Its schema-2 summaries require the newer reader; fixed-population summaries
retain schema 1. Retain the matching model-data.rds and summarize.R.

Fixed facet locations, group means and step anchors retain their source
values, units and centering convention. model-data.rds$constraints records
them. locations = location_offset + location_basis * location_free, and
the stacked step vector = step_offset + step_basis * step_free. Step rows
are ordered first by owner, then transition. These affine maps enforce the
constraints exactly, including parameters determined by the remaining sum.
Steps need not be ordered. There is no sampling of fixed coordinates.

Normal priors are specified on orthonormal free contrasts, centered at zero.
The affine prior mean shares the required remaining total equally among the
free levels in each centered block. Coordinate SDs are scaled so each varying
location, step or component log slope has the requested marginal prior SD.
Independent uncentered coordinates use that SD directly. Fixed or uniquely
determined levels instead have zero prior SD. These are declared priors on
the free subspace, not an automatic conditional prior from an earlier fit.
levels.csv and steps.csv give prior means, marginal SDs and Fixed indicators;
slope-levels.csv maps GPCM component slopes to their owners.

All inference conditions on the supplied anchors and group means. Their
zero posterior spread does not mean the earlier calibration was known without
error. This export does not propagate uncertain external calibration.
Person anchors and JML sources remain outside this MML exporter. Population
regression, non-unit weights, interactions, dummy facets, post-fit shrinkage,
testlets and random-rater models are not supported. This export does not add
a native step-anchor fitting argument or enable anchored two-family fitting.

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

Both formats support Rscript run.R --summarize to regenerate the five report
files from fit.rds without compiling or sampling. Run it from this directory;
keep model-data.rds, summarize.R and (for brms) model.R with the saved fit.
CmdStanR/posterior or brms are needed for the respective replay, but a CmdStan
toolchain is not. brms also regenerates parameters.rds. Draws are unchanged.

The brms script saves fit.rds, parameters.rds and predictive-draws.rds; its
custom log_lik(), posterior_predict() and posterior_epred() callbacks survive
saved replay. row_id identifies the retained rating design, including its
step/slope owners. Predictions condition on sampled abilities of the existing
Persons and retained design rows; they are not new-Person predictions or
leave-one-Person-out scores. Repeated ratings retain conditional independence.

These proper structural priors are Bayesian additions to the response model;
posterior summaries do not reproduce MML estimates or confidence intervals.
Neither exporting nor completing chains validates the source fit's inference.
For prior sensitivity, export again with the intended prior arguments into a
separate directory, keeping its model-data, fit and summaries together.

https://mc-stan.org/cmdstanr/
https://mc-stan.org/learn-stan/diagnostics-warnings.html
https://paulbuerkner.com/brms/articles/brms_customfamilies.html

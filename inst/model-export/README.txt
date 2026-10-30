Bayesian RSM/PCM model exported by mfrmr

Run Rscript run.R from this directory. Export itself requires neither Stan nor
brms. Estimation requires CmdStanR and CmdStan; the brms format also needs brms.
Review and adjust the run script for your computing resources. The scripts
request four chains and retain diagnostics. Inspect divergences, R-hat,
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

Both formats support Rscript run.R --summarize to regenerate the five report
files from fit.rds without compiling or sampling. Run it from this directory;
keep model-data.rds, summarize.R and (for brms) model.R with the saved fit.
CmdStanR/posterior or brms are needed for the respective replay, but a CmdStan
toolchain is not. brms also regenerates parameters.rds. Draws are unchanged.

model-data.rds is a base-R list containing data, design matrices, mappings and
prior settings. levels.csv contains original Person and facet labels. Category
1 is the lowest source score; categories.csv provides the complete ladder,
including any unobserved declared categories. steps[g, k] is the offset for
transition k to k+1 for owner g; model-data.rds$step_owners gives owner labels.
Locations include the original sign separately: higher values decrease the
adjacent log odds for a facet with Sign = -1. The brms helper
extract_mfrm_parameters(fit) returns parameters on this same scale.

Both formats implement:
  log(P(Y=k+1) / P(Y=k)) = theta + signed facet locations - step[k].
This is the logistic metric (L-metric), with scaling constant D = 1 and no
extra 1.7 or 1.702 normal-ogive approximation factor. Increasing ability by
one unit increases adjacent-category log odds by one, holding other values
fixed. This response-link convention is distinct from the ability prior.
Abilities independently follow normal(0,1), retaining the source MML population
scale without centering the observed Persons. Centered facet locations and
each owner's steps sum to zero. Steps need not be ordered. Centered components
use an orthonormal contrast basis with independent normal coordinates scaled
by sqrt(m/(m-1)); each of the m components then has the marginal prior SD
recorded in model-data.rds$priors. A centered singleton is fixed at zero.
Uncentered locations have independent normal priors with the recorded SD.

These structural priors are explicit Bayesian additions to MML. Posterior
credible intervals do not reproduce or validate MML confidence intervals.
Inspect prior sensitivity. To change priors coherently in either format, change
the location_sd and step_sd coordinate scales in model-data.rds$stan_data;
the priors field is descriptive metadata and should be updated as well.

The brms custom family uses the same adjacent-category likelihood, constraints
and priors as the direct Stan model. Its posterior_epred() returns category
probabilities. posterior_predict() and log_lik() are also supplied and survive
saving/reloading the fitted object. Predictions and pointwise log likelihoods
condition on sampled abilities for existing Persons; these are not new-Person
predictions or leave-one-Person-out validation. Repeated rows retain the
source model's conditional-independence assumption.

This bundle represents ordinary unanchored, unit-weight RSM/PCM MML sources
without population regression, interactions, dummy facets or shrinkage.
Anchored and one-family GPCM models use a separate affine template.
Two-family GPCM has a separate Stan-only export; this RSM/PCM bundle does not
implement it. JML, testlet and random-rater export remain unsupported. Export does not
certify source convergence or subsequent sampling. Only retained observations
are exported; omitted input rows and the original optimizer are not replayed.

The bundle contains ratings and original labels. Review before sharing.
Official documentation:
https://mc-stan.org/cmdstanr/
https://mc-stan.org/learn-stan/diagnostics-warnings.html
https://mc-stan.org/posterior/reference/mcse_quantile.html
https://paulbuerkner.com/brms/articles/brms_customfamilies.html

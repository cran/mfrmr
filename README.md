# mfrmr

[![GitHub](https://img.shields.io/badge/GitHub-mfrmr-181717?logo=github)](https://github.com/Ryuya-dot-com/R_package_mfrmr)
[![R-CMD-check](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/pkgdown.yaml)
[![test-coverage](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/test-coverage.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Native R package for many-facet Rasch model (MFRM) estimation, diagnostics,
and reporting workflows.

## Start here first

If you are new to `mfrmr`, use this route first and ignore the longer feature
lists below until it works end to end.

- Fit with `method = "MML"`
- Diagnose with `diagnostic_mode = "both"` for `RSM` / `PCM`
- Read `summary(fit)` and `summary(diag)` before branching into plots/reports

```r
library(mfrmr)
toy <- load_mfrmr_data("example_core")

fit <- fit_mfrm(
  toy,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM",
  quad_points = 15
)

diag <- diagnose_mfrm(
  fit,
  diagnostic_mode = "both",
  residual_pca = "none"
)

summary(fit)
summary(diag)
plot_qc_dashboard(fit, diagnostics = diag, preset = "publication")
chk <- reporting_checklist(fit, diagnostics = diag)
```

If that route works, the next natural step is:

- reporting: `build_apa_outputs()` / `apa_table()`
- misfit case review: `build_misfit_casebook()`
- weighting review: `compare_mfrm()` -> `build_weighting_audit()`
- strict follow-up: `plot_marginal_fit()` / `plot_marginal_pairwise()`
- operational linking review: `audit_mfrm_anchors()` -> `detect_anchor_drift()` -> `build_linking_review()`
- linking/design: `subset_connectivity_report()`

## What this package is for

`mfrmr` is designed around four package-native routes:

- Estimation and diagnostics: `fit_mfrm()` -> `diagnose_mfrm()`
- Reporting and manuscript preparation: `reporting_checklist()` -> `build_apa_outputs()`
- Misfit case review: `build_misfit_casebook()` -> `casebook$group_view_index` /
  `casebook$group_views` -> source-specific follow-up plots
- Linking, anchors, drift, and DFF:
  `audit_mfrm_anchors()` / `detect_anchor_drift()` -> `build_linking_review()`
  or `subset_connectivity_report()` -> `anchor_to_baseline()` / `analyze_dff()`
- Legacy-compatible export when required: `run_mfrm_facets()` and related compatibility helpers

If you want the shortest possible recommendation:

- Final estimation: prefer `method = "MML"`
- Fast exploratory pass only: use `method = "JML"`
- Preferred `RSM` / `PCM` fit screen: `diagnose_mfrm(..., diagnostic_mode = "both")`
- First visual screen: `plot_qc_dashboard(..., preset = "publication")`
- First reporting screen: `reporting_checklist()`
- First case-review screen: `build_misfit_casebook()` and then inspect
  `casebook$group_view_index`
- First weighting-policy screen: `build_weighting_audit()`
- First operational linking screen: `build_linking_review()`

## Minimum input contract

`mfrmr` expects long-format rating data: one row per observed rating.

- Required columns:
  - one person column
  - one ordered score column
  - one or more non-person facet columns supplied in `facets = c(...)`
- Score-column rules:
  - scores should be ordered integer category codes such as `0/1`, `1/2`, or `1:5`
  - binary two-category scores are supported under the same `RSM` / `PCM` interface
  - fractional scores are rejected; recode them explicitly before fitting
  - non-numeric score labels are dropped with a warning if coercion fails
  - when `keep_original = FALSE`, unused intermediate categories are collapsed to a contiguous internal scale and recorded in `fit$prep$score_map`
  - if the intended scale has unused boundary categories, such as a 1-5 scale with only 2-5 observed, set `rating_min = 1, rating_max = 5` so the zero-count boundary category remains in the fitted support
  - if the intended scale has unused intermediate categories, such as 1, 3, 5 observed on a 1-5 scale, also set `keep_original = TRUE`
  - `summary(describe_mfrm_data(...))` reports retained zero-count categories in `Notes`, the printed `Caveats` block, and `$caveats`; `summary(fit)` carries the full structured rows into printed `Caveats` and appendix/export role `analysis_caveats`, with `Key warnings` as a short triage subset
- Typical optional columns:
  - `Subset` for disconnected-form or linking work
  - `Weight` for weighted analyses
  - `Group` when downstream fairness or DFF workflows need grouping metadata

Minimal pattern:

```r
names(df)
# [1] "Person" "Rater" "Criterion" "Score"

fit <- fit_mfrm(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM"
)
```

## Main capabilities

Core analysis:

- estimation with `fit_mfrm()` under `MML` or `JML`
- fit diagnostics with `diagnose_mfrm()`, `plot_qc_dashboard()`, and residual PCA follow-up
- strict marginal follow-up for `RSM` / `PCM` via `diagnostic_mode = "both"`, `plot_marginal_fit()`, and `plot_marginal_pairwise()`
- package-native tables and summaries via `summary()`, `reporting_checklist()`, and `facet_statistics_report()`

Reporting and QA:

- APA/report drafting with `build_apa_outputs()`, `apa_table()`, and `build_summary_table_bundle()`
- visual/report routing with `build_visual_summaries()` and `reporting_checklist()`
- QC workflows with `run_qc_pipeline()` and `plot_qc_pipeline()`
- reproducible export helpers such as `export_mfrm_bundle()`, `build_mfrm_manifest()`, and `build_mfrm_replay_script()`

Linking, fairness, and advanced review:

- bias and DFF workflows through `estimate_bias()`, `estimate_all_bias()`, `analyze_dff()`, and `dif_report()`
- anchoring and linking via `anchor_to_baseline()`, `detect_anchor_drift()`, and `build_equating_chain()`
- precision/targeting views via `compute_information()`, `plot_information()`, and `plot_wright_unified()`
- equivalence and audit helpers such as `analyze_facet_equivalence()`, `describe_mfrm_data()`, and `audit_mfrm_anchors()`

Advanced or compatibility scope:

- legacy-compatible one-shot wrapper: `run_mfrm_facets()` / `mfrmRFacets()`
- simulation and planning helpers: `simulate_mfrm_data()`, `evaluate_mfrm_design()`, `build_mfrm_sim_spec()`, `extract_mfrm_sim_spec()`, `predict_mfrm_population()`
- future-unit posterior scoring: `predict_mfrm_units()` and `sample_mfrm_plausible_values()`

## Latent regression status

`mfrmr` now includes a first-version latent-regression branch inside
`fit_mfrm()`. Activate it with `method = "MML"`,
`population_formula = ~ ...`, and one-row-per-person `person_data`.

Current supported boundary:

- ordered-response `RSM` / `PCM`
- one latent dimension
- conditional-normal person population model
- person covariates supplied through an explicit person table and expanded
  through `stats::model.matrix()`, including numeric/logical predictors and
  factor/character categorical predictors
- posterior scoring and plausible-value draws that condition on the fitted
  population model

What to inspect after fitting:

- `summary(fit)$population_overview` shows the posterior basis, residual
  variance, and any omitted-person counts.
- `summary(fit)$population_coefficients` shows the latent-regression
  coefficients.
- `summary(fit)$population_coding` shows how categorical covariates were coded.
- `summary(fit)$key_warnings` and `summary(fit)$caveats` flag issues that
  should be reviewed before reporting or exporting results.

Beginner quick start:

```r
# response data: one row per rating event
# person data: one row per person, with the same Person IDs
person_tbl <- unique(dat[c("Person", "Grade", "Group")])

fit_pop <- fit_mfrm(
  data = dat,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM",
  population_formula = ~ Grade + Group,
  person_data = person_tbl,
  population_policy = "error"
)

s_pop <- summary(fit_pop)
s_pop$population_overview      # posterior basis, residual variance, omissions
s_pop$population_coefficients  # latent-regression coefficients
s_pop$population_coding        # categorical levels / contrasts / encoded columns
s_pop$caveats                  # complete-case and category-support warnings
```

Use `population_policy = "omit"` only when complete-case removal is intended,
then report the omitted-person and omitted-row counts. Coefficients in
`population_coefficients` are conditional-normal population-model parameters,
not a post hoc regression on EAP/MLE scores.

Reference checks for this branch:

```r
bench_pop <- reference_case_benchmark(
  cases = c("synthetic_latent_regression_omit", "synthetic_conquest_overlap_dry_run"),
  method = "MML",
  model = "RSM",
  quad_points = 5,
  maxit = 30
)

summary(bench_pop)
bench_pop$population_policy_checks  # complete-case omission check
bench_pop$conquest_overlap_checks   # package-side ConQuest preparation check
```

The ConQuest preparation case checks only package-side preparation. It does not run
ConQuest. When actual ConQuest output tables are available for the documented
overlap case, use the external-table comparison helpers:

```r
bundle <- build_conquest_overlap_bundle(fit_overlap, output_dir = "conquest_overlap")
normalized <- normalize_conquest_overlap_files(
  population_file = "conquest_population.csv",
  item_file = "conquest_items.csv",
  case_file = "conquest_cases.csv"
)
audit <- audit_conquest_overlap(bundle, normalized)
summary(audit)$summary
audit$attention_items
```

Treat this as a scoped comparison, not as full ConQuest numerical equivalence.
ConQuest must be run separately and the extracted tables must be reviewed.

Current non-goals for this branch:

- `JML` latent regression
- bounded `GPCM` latent regression
- multidimensional population models
- arbitrary imported design specifications
- the full ConQuest plausible-values workflow

This should be described as first-version overlap with the ConQuest
latent-regression framework, not as ConQuest numerical equivalence.

`predict_mfrm_population()` remains a simulation-based scenario-forecasting
helper. It should not be described as the latent-regression estimator itself.

## Bounded GPCM support

`GPCM` is now part of the supported core package scope, but only within a
bounded route. Use `gpcm_capability_matrix()` to see the current release
boundary in one place.

- Supported core: fitting, `summary()` / `print()`, posterior scoring,
  `compute_information()`, Wright/pathway/CCC plots, and category reports.
- Supported with caveat: `diagnose_mfrm()` and direct slope-aware simulation
  are exploratory; `reporting_checklist()` routes only the direct table/plot
  path.
- Not supported in this release: fair-average reporting, score-side exports,
  APA/report bundles, QC pass/fail pipelines, linking synthesis, planning /
  forecasting, posterior predictive computation, and `MCMC`.

The unsupported helpers depend on score-side or planning assumptions that are
validated for the Rasch-family route but not yet generalized to bounded `GPCM`.

## Equal weighting and when to prefer Rasch-MFRM

`mfrmr` treats `RSM` / `PCM` as the package's equal-weighting reference
models. In that Rasch-family route, category discrimination is fixed, so the
operational scoring contract does not let the psychometric model reweight some
item-facet combinations more heavily than others.

bounded `GPCM` serves a different purpose. It allows estimated slopes, so some
observed design cells become more influential than others through
discrimination-based reweighting. This often improves fit, but a better-fitting
`GPCM` does not automatically make it the preferred operational model.

The package therefore recommends:

- prefer `RSM` / `PCM` when equal contributions of items and raters are part of
  the substantive scoring argument
- use bounded `GPCM` when you explicitly want to inspect or allow
  discrimination-based reweighting and can defend that choice on validity
  grounds
- read `RSM` / `PCM` versus `GPCM` as a model-choice or sensitivity question,
  not as a contest in which fit alone decides the winner

One more distinction matters. The `weight =` argument in `fit_mfrm()` is for an
observation-weight column. That is different from the equal-weighting question
discussed above. Observation weights adjust how rating events enter estimation
and summaries; they do not turn a Rasch-family fit into a discrimination-based
model.

## Documentation map

The README is only the shortest map. The package now has guide-style help pages
for the main workflows.

- Workflow map:
  `help("mfrmr_workflow_methods", package = "mfrmr")`
- Visual diagnostics map:
  `help("mfrmr_visual_diagnostics", package = "mfrmr")`
- Reports and tables map:
  `help("mfrmr_reports_and_tables", package = "mfrmr")`
- Reporting and APA map:
  `help("mfrmr_reporting_and_apa", package = "mfrmr")`
- Linking and DFF map:
  `help("mfrmr_linking_and_dff", package = "mfrmr")`
- Compatibility layer map:
  `help("mfrmr_compatibility_layer", package = "mfrmr")`
- Bounded `GPCM` scope:
  `help("gpcm_capability_matrix", package = "mfrmr")`

Companion vignettes:

- `vignette("mfrmr-workflow", package = "mfrmr")`
- `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`
- `vignette("mfrmr-reporting-and-apa", package = "mfrmr")`
- `vignette("mfrmr-linking-and-dff", package = "mfrmr")`
- `vignette("mfrmr-mml-and-marginal-fit", package = "mfrmr")`

## Installation

```r
# GitHub (development version)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Ryuya-dot-com/R_package_mfrmr", build_vignettes = TRUE)

# CRAN (after release)
# install.packages("mfrmr")
```

If you install from GitHub without `build_vignettes = TRUE`, use the guide-style
help pages included in the package, for example:

- `help("mfrmr_workflow_methods", package = "mfrmr")`
- `help("mfrmr_reporting_and_apa", package = "mfrmr")`
- `help("mfrmr_linking_and_dff", package = "mfrmr")`

Installed vignettes:

```r
browseVignettes("mfrmr")
```

## Core workflow

```
fit_mfrm() --> diagnose_mfrm() --> reporting / advanced analysis
                    |
                    +--> analyze_residual_pca()
                    +--> estimate_bias()
                    +--> analyze_dff()
                    +--> compare_mfrm()
                    +--> run_qc_pipeline()
                    +--> anchor_to_baseline() / detect_anchor_drift()
```

1. Fit model: `fit_mfrm()`
2. Diagnostics: `diagnose_mfrm()`
3. Optional residual PCA: `analyze_residual_pca()`
4. Optional interaction bias: `estimate_bias()`
5. Differential-functioning analysis: `analyze_dff()`, `dif_report()`
6. Model comparison: `compare_mfrm()`
7. Reporting: `apa_table()`, `build_apa_outputs()`, `build_visual_summaries()`
8. Quality control: `run_qc_pipeline()`
9. Anchoring & linking: `anchor_to_baseline()`, `detect_anchor_drift()`, `build_equating_chain()`
10. Compatibility-contract audit when needed: `facets_parity_report()`;
    this audits package output contracts, not external FACETS numerical
    equivalence
11. Reproducible inspection: `summary()` and `plot(..., draw = FALSE)`

## Choose a route

Use the route that matches the question you are trying to answer.

| Question | Recommended route |
|---|---|
| Can I fit the model and get a first-pass diagnosis quickly? | `fit_mfrm()` -> `diagnose_mfrm()` -> `plot_qc_dashboard()` |
| Which reporting elements are draft-complete, and with what caveats? | `diagnose_mfrm()` -> `precision_audit_report()` -> `reporting_checklist()` |
| Which tables and prose should I adapt into a manuscript draft? | `reporting_checklist()` -> `build_apa_outputs()` -> `apa_table()` |
| Is the design connected well enough for a common scale? | `subset_connectivity_report()` -> `plot(..., type = "design_matrix")` |
| Do I need to place a new administration onto a baseline scale? | `make_anchor_table()` -> `anchor_to_baseline()` |
| Are common elements stable across separately fitted forms or waves? | fit each wave -> `detect_anchor_drift()` -> `build_equating_chain()` |
| Are some facet levels functioning differently across groups? | `subset_connectivity_report()` -> `analyze_dff()` -> `dif_report()` |
| Do I need old fixed-width or wrapper-style outputs? | `run_mfrm_facets()` or `build_fixed_reports()` only at the compatibility boundary |

## Additional routes

After the canonical `MML + both` route above, these are the next shortest
specialized routes.

Shared setup used by the snippets below:

```r
library(mfrmr)
toy <- load_mfrmr_data("example_core")
```

### 1. Quick first pass

```r
fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                method = "MML", model = "RSM", quad_points = 7)
diag <- diagnose_mfrm(fit, diagnostic_mode = "both", residual_pca = "none")
summary(diag)
plot_qc_dashboard(fit, diagnostics = diag, preset = "publication")
```

### 1b. Preferred MML + marginal-fit route

```r
fit_final <- fit_mfrm(
  toy,
  "Person",
  c("Rater", "Criterion"),
  "Score",
  method = "MML",
  model = "RSM",
  quad_points = 15
)

diag_final <- diagnose_mfrm(
  fit_final,
  diagnostic_mode = "both",
  residual_pca = "none"
)

summary(fit_final)
summary(diag_final)
```

For `RSM` / `PCM`, this is the recommended final-analysis route when you want
legacy continuity plus the newer strict marginal screening path.

### 2. Design and linking check

```r
diag <- diagnose_mfrm(fit, residual_pca = "none")
sc <- subset_connectivity_report(fit, diagnostics = diag)
summary(sc)
plot(sc, type = "design_matrix", preset = "publication")
plot_wright_unified(fit, preset = "publication", show_thresholds = TRUE)
```

### 3. Manuscript and reporting check

```r
# Add `bias_results = ...` if you want the bias/reporting layer included.
chk <- reporting_checklist(fit, diagnostics = diag)
apa <- build_apa_outputs(fit, diag)

chk$checklist[, c("Section", "Item", "DraftReady", "NextAction")]
cat(apa$report_text)
```

## Estimation choices

The package treats `MML` and `JML` differently on purpose.

- `MML` is the default and the preferred route for final estimation.
- `JML` is supported as a fast exploratory route.
- Downstream precision summaries distinguish `model_based`, `hybrid`, and `exploratory` tiers.
- Use `precision_audit_report()` when you need to decide how strongly to phrase SE, CI, or reliability claims.

Typical pattern:

```r
toy <- load_mfrmr_data("example_core")

fit_final <- fit_mfrm(
  toy, "Person", c("Rater", "Criterion"), "Score",
  method = "MML", model = "RSM", quad_points = 15
)

diag_final <- diagnose_mfrm(
  fit_final,
  diagnostic_mode = "both",
  residual_pca = "none"
)

precision_audit_report(fit_final, diagnostics = diag_final)
```

## Mathematical note for expert users

The package help and the `mfrmr-mml-and-marginal-fit` vignette carry the full
expert-facing definitions, but the two key ideas are:

1. `MML` optimizes the marginal likelihood

\[
L(\beta) = \prod_{n=1}^{N} \sum_{q=1}^{Q} w_q \, p(\mathbf{x}_n \mid \theta_q, \beta),
\]

with Gauss-Hermite nodes \(\theta_q\) and weights \(w_q\), and then forms
posterior summaries such as EAP scores from the same latent-integrated bundle.

2. The strict marginal branch compares observed summaries with
posterior-integrated expectations rather than `EAP` plug-in residuals. In
schematic form, for a flagged cell or grouped summary \(g\),

\[
\mathbb{E}_{\hat\beta}(N_{gc}) = \sum_{n=1}^{N}\sum_{q=1}^{Q}
\Pr(\theta_q \mid \mathbf{x}_n, \hat\beta) \, I(n \in g)\,
\Pr(X_n = c \mid \theta_q, \hat\beta).
\]

This is why `legacy` residual diagnostics and `marginal_fit` diagnostics are
reported separately in `mfrmr`.

In many-facet reporting, keep `facet-level reliability` and `inter-rater
agreement` separate. The former summarizes separation/precision on the latent
scale; the latter summarizes observed agreement across matched rater contexts.

Literature positioning:

- for MML, see Bock and Aitkin (1981)
- for many-facet severity/fit/agreement practice, see Linacre (1989) and
  Eckes (2005)
- for limited-information item-fit logic that motivates, but is not identical to,
  this many-facet screening layer, see Orlando and Thissen (2000, 2003)
- for generalized residual logic that informs, but does not fully calibrate, the
  current pairwise screening summaries, see Haberman and Sinharay (2013)
- for posterior-predictive follow-up as a separate model-checking family, see
  Sinharay et al. (2006)
- for the broader fit-assessment landscape and practical-significance framing,
  see Sinharay and Monroe (2025)

For the longer expert note, run:

```r
vignette("mfrmr-mml-and-marginal-fit", package = "mfrmr")
```

## Help-page navigation

Core analysis help pages include practical sections such as:

- `Interpreting output`
- `Typical workflow`

Recommended entry points:

- `?mfrmr-package` (package overview)
- `?fit_mfrm`, `?diagnose_mfrm`, `?run_mfrm_facets`
- `vignette("mfrmr-mml-and-marginal-fit", package = "mfrmr")`
- `?analyze_dff`, `?analyze_dif`, `?compare_mfrm`, `?run_qc_pipeline`
- `?anchor_to_baseline`, `?detect_anchor_drift`, `?build_equating_chain`
- `?build_apa_outputs`, `?build_visual_summaries`, `?apa_table`
- `?reporting_checklist`, `?facet_quality_dashboard`, `?estimate_all_bias`
- `?export_mfrm_bundle`, `?build_mfrm_manifest`, `?build_mfrm_replay_script`
- `?analyze_facet_equivalence`, `?plot_facet_equivalence`
- `?mfrmr_workflow_methods`, `?mfrmr_visual_diagnostics`
- `?mfrmr_reports_and_tables`, `?mfrmr_reporting_and_apa`
- `?mfrmr_linking_and_dff`, `?mfrmr_compatibility_layer`

Utility pages such as `?export_mfrm`, `?as.data.frame.mfrm_fit`, and
`?plot_bubble` also include lightweight export / plotting examples.

## Performance tips

- Start with `method = "JML"` only when you explicitly want a quick exploratory fit.
- Prefer `method = "MML"` for final estimation, but tune `quad_points` to match your workflow.
- `quad_points = 7` is a good fast iteration setting; `quad_points = 15` is a better final-analysis setting.
- Use `diagnose_mfrm(fit, residual_pca = "none")` for a quick first pass, then add residual PCA only when needed.
- Reuse diagnostics objects in downstream helpers such as `plot_bubble()` and `run_qc_pipeline()` to avoid repeated work.

## Documentation datasets

- `load_mfrmr_data("example_core")`: compact, approximately unidimensional example for fitting, diagnostics, plots, and reports.
- `load_mfrmr_data("example_bias")`: compact example with known `Group x Criterion` differential-functioning and `Rater x Criterion` interaction signals for bias-focused help pages.
- `load_mfrmr_data("study1")` / `load_mfrmr_data("study2")`: larger Eckes/Jin-inspired synthetic studies for more realistic end-to-end analyses.
- Direct dataset access also works with `data("mfrmr_example_core", package = "mfrmr")` and `data("mfrmr_example_bias", package = "mfrmr")`.

## Quick start

```r
library(mfrmr)

data("mfrmr_example_core", package = "mfrmr")
df <- mfrmr_example_core

# Fit
fit <- fit_mfrm(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",
  model = "RSM",
  quad_points = 7
)
summary(fit)

# Fast diagnostics first
diag <- diagnose_mfrm(fit, residual_pca = "none")
summary(diag)

# APA outputs
apa <- build_apa_outputs(fit, diag)
cat(apa$report_text)

# QC pipeline reuses the same diagnostics object
qc <- run_qc_pipeline(fit, diagnostics = diag)
summary(qc)
```

## Main objects you will reuse

Most package workflows reuse a small set of objects rather than recomputing
everything from scratch.

- `fit`: the fitted many-facet Rasch model returned by `fit_mfrm()`
- `diag`: diagnostic summaries returned by `diagnose_mfrm()`
- `chk`: reporting and manuscript-draft checks returned by `reporting_checklist()`
- `apa`: structured APA/report draft outputs returned by `build_apa_outputs()`
- `sc`: connectivity and linking summaries returned by `subset_connectivity_report()`
- `bias` / `dff`: interaction screening and differential-functioning results returned by
  `estimate_bias()` and `analyze_dff()`

Typical reuse pattern:

```r
toy <- load_mfrmr_data("example_core")

fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                method = "MML", model = "RSM", quad_points = 7)
diag <- diagnose_mfrm(fit, residual_pca = "none")
chk <- reporting_checklist(fit, diagnostics = diag)
apa <- build_apa_outputs(fit, diag)
sc <- subset_connectivity_report(fit, diagnostics = diag)
```

## Reporting and APA route

If your endpoint is a manuscript or internal report, use the package-native
reporting contract rather than composing text by hand.

```r
diag <- diagnose_mfrm(fit, residual_pca = "none")

# Add `bias_results = ...` to either helper when bias screening should
# appear in the checklist or draft text.
chk <- reporting_checklist(fit, diagnostics = diag)
chk$checklist[, c("Section", "Item", "DraftReady", "Priority", "NextAction")]

apa <- build_apa_outputs(
  fit,
  diag,
  context = list(
    assessment = "Writing assessment",
    setting = "Local scoring study",
    scale_desc = "0-4 rubric scale",
    rater_facet = "Rater"
  )
)

cat(apa$report_text)
apa$section_map[, c("SectionId", "Available", "Heading")]

tbl_fit <- apa_table(fit, which = "summary")
tbl_reliability <- apa_table(fit, which = "reliability", diagnostics = diag)
```

For a question-based map of the reporting API, see
`help("mfrmr_reporting_and_apa", package = "mfrmr")`.

## Visualization recipes

If you want a question-based map of the plotting API, see
`help("mfrmr_visual_diagnostics", package = "mfrmr")`.

```r
# Which plotting routes are supported for this run?
chk <- reporting_checklist(fit, diagnostics = diag)
chk$visual_scope

# Where should each figure go in a paper or appendix?
visual_reporting_template("manuscript")[, c("FigureFamily", "CaptionSkeleton")]
visual_reporting_template("surface")[, c("FigureFamily", "ResultsWording")]

# Wright map with shared targeting view
plot(fit, type = "wright", preset = "publication", show_ci = TRUE)

# Pathway map with dominant-category strips
plot(fit, type = "pathway", preset = "publication")

# Category probability curves
plot(fit, type = "ccc", preset = "publication")

# 3D-ready category probability surface payload.
# This is exploratory data for downstream rendering, not a default APA figure.
surface <- plot(fit, type = "ccc_surface", draw = FALSE)
head(surface$data$surface)
surface$data$category_support

# Linking design matrix
sc <- subset_connectivity_report(fit, diagnostics = diag)
plot(sc, type = "design_matrix", preset = "publication")

# Unexpected responses
plot_unexpected(fit, diagnostics = diag, preset = "publication")

# Displacement screening
plot_displacement(fit, diagnostics = diag, preset = "publication")

# Facet variability overview
plot_facets_chisq(fit, diagnostics = diag, preset = "publication")

# Residual PCA scree and loadings
pca <- analyze_residual_pca(diag, mode = "both")
plot_residual_pca(pca, mode = "overall", plot_type = "scree", preset = "publication")

# Bias screening profile
bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion")
plot_bias_interaction(bias, plot = "facet_profile", preset = "publication")

# One-page QC screen
plot_qc_dashboard(fit, diagnostics = diag, preset = "publication")
```

## Linking, anchors, and DFF route

Use this route when your design spans forms, waves, or subgroup comparisons.

```r
data("mfrmr_example_bias", package = "mfrmr")
df_bias <- mfrmr_example_bias
fit_bias <- fit_mfrm(df_bias, "Person", c("Rater", "Criterion"), "Score",
                     method = "MML", model = "RSM", quad_points = 7)
diag_bias <- diagnose_mfrm(fit_bias, residual_pca = "none")

# Connectivity and design coverage
sc <- subset_connectivity_report(fit_bias, diagnostics = diag_bias)
summary(sc)
plot(sc, type = "design_matrix", preset = "publication")

# Anchor export from a baseline fit
anchors <- make_anchor_table(fit_bias, facets = "Criterion")
head(anchors)

# Differential facet functioning
dff <- analyze_dff(
  fit_bias,
  diag_bias,
  facet = "Criterion",
  group = "Group",
  data = df_bias,
  method = "residual"
)
dff$summary
plot_dif_heatmap(dff)
```

For linking-specific guidance, see
`help("mfrmr_linking_and_dff", package = "mfrmr")`.

## DFF / DIF analysis

```r
data("mfrmr_example_bias", package = "mfrmr")
df_bias <- mfrmr_example_bias
fit_bias <- fit_mfrm(df_bias, "Person", c("Rater", "Criterion"), "Score",
                     method = "MML", model = "RSM", quad_points = 7)
diag_bias <- diagnose_mfrm(fit_bias, residual_pca = "none")

dff <- analyze_dff(fit_bias, diag_bias, facet = "Criterion",
                   group = "Group", data = df_bias, method = "residual")
dff$dif_table
dff$summary

# Cell-level interaction table
dit <- dif_interaction_table(fit_bias, diag_bias, facet = "Criterion",
                             group = "Group", data = df_bias)

# Visual, narrative, and bias reports
plot_dif_heatmap(dff)
dr <- dif_report(dff)
cat(dr$narrative)

# Refit-based contrasts can support ETS labels only when subgroup linking is adequate
dff_refit <- analyze_dff(fit_bias, diag_bias, facet = "Criterion",
                         group = "Group", data = df_bias, method = "refit")
dff_refit$summary

bias <- estimate_bias(fit_bias, diag_bias, facet_a = "Rater", facet_b = "Criterion")
summary(bias)

# App-style batch bias estimation across all modeled facet pairs
bias_all <- estimate_all_bias(fit_bias, diag_bias)
bias_all$summary
```

Interpretation rules:

- `residual` DFF is a screening route.
- `refit` DFF can support logit-scale contrasts only when subgroup linking is adequate.
- Check `ScaleLinkStatus`, `ContrastComparable`, and the reported classification system before treating a contrast as a strong interpretive claim.

## Model comparison

```r
fit_rsm <- fit_mfrm(df, "Person", c("Rater", "Criterion"), "Score",
                     method = "MML", model = "RSM")
fit_pcm <- fit_mfrm(df, "Person", c("Rater", "Criterion"), "Score",
                     method = "MML", model = "PCM", step_facet = "Criterion")
cmp <- compare_mfrm(RSM = fit_rsm, PCM = fit_pcm)
cmp$table

# Request nested tests only when models are truly nested and fit on the same basis
cmp_nested <- compare_mfrm(RSM = fit_rsm, PCM = fit_pcm, nested = TRUE)
cmp_nested$comparison_basis

# RSM design-weighted precision curves
info <- compute_information(fit_rsm)
plot_information(info)
```

## Design simulation

```r
spec <- build_mfrm_sim_spec(
  n_person = 50,
  n_rater = 4,
  n_criterion = 4,
  raters_per_person = 2,
  assignment = "rotating",
  model = "RSM"
)

sim_eval <- evaluate_mfrm_design(
  n_person = c(30, 50, 80),
  n_rater = 4,
  n_criterion = 4,
  raters_per_person = 2,
  reps = 2,
  maxit = 15,
  sim_spec = spec,
  seed = 123
)

s_sim <- summary(sim_eval)
s_sim$design_summary
s_sim$ademp

rec <- recommend_mfrm_design(sim_eval)
rec$recommended

plot(sim_eval, facet = "Rater", metric = "separation", x_var = "n_person")
plot(sim_eval, facet = "Criterion", metric = "severityrmse", x_var = "n_person")
```

Notes:

- Use `build_mfrm_sim_spec()` when you want one explicit, reusable data-generating mechanism.
- Use `extract_mfrm_sim_spec(fit)` when you want a fit-derived starting point for a later design study.
- Use `extract_mfrm_sim_spec(fit, latent_distribution = "empirical", assignment = "resampled")` when you want a more semi-parametric design study that reuses empirical fitted spreads and observed rater-assignment profiles.
- Use `extract_mfrm_sim_spec(fit, latent_distribution = "empirical", assignment = "skeleton")` when you want a more plasmode-style study that preserves the observed person-by-facet design skeleton and resimulates only the responses.
- `summary(sim_eval)$ademp` records the simulation-study contract: aims, DGM, estimands, methods, and performance measures.

## Population forecast

```r
spec_pop <- build_mfrm_sim_spec(
  n_person = 50,
  n_rater = 4,
  n_criterion = 4,
  raters_per_person = 2,
  assignment = "rotating",
  model = "RSM"
)

pred_pop <- predict_mfrm_population(
  sim_spec = spec_pop,
  n_person = 60,
  reps = 2,
  maxit = 15,
  seed = 123
)

s_pred <- summary(pred_pop)
s_pred$forecast[, c("Facet", "MeanSeparation", "McseSeparation")]
```

Notes:

- `predict_mfrm_population()` forecasts aggregate operating characteristics for one future design.
- It does not return deterministic future person or rater true values.

## Future-unit posterior scoring

```r
toy_pred <- load_mfrmr_data("example_core")
toy_fit <- fit_mfrm(
  toy_pred,
  "Person", c("Rater", "Criterion"), "Score",
  method = "MML",
  quad_points = 7
)

raters <- unique(toy_pred$Rater)[1:2]
criteria <- unique(toy_pred$Criterion)[1:2]

new_units <- data.frame(
  Person = c("NEW01", "NEW01", "NEW02", "NEW02"),
  Rater = c(raters[1], raters[2], raters[1], raters[2]),
  Criterion = c(criteria[1], criteria[2], criteria[1], criteria[2]),
  Score = c(2, 3, 2, 4)
)

pred_units <- predict_mfrm_units(toy_fit, new_units, n_draws = 0)
summary(pred_units)$estimates[, c("Person", "Estimate", "Lower", "Upper")]

pv_units <- sample_mfrm_plausible_values(
  toy_fit,
  new_units,
  n_draws = 3,
  seed = 123
)
summary(pv_units)$draw_summary[, c("Person", "Draws", "MeanValue")]
```

Notes:

- `predict_mfrm_units()` scores future or partially observed persons under the
  fitted scoring basis.
- For ordinary `MML` fits, that basis is the fitted marginal calibration.
- For latent-regression `MML` fits with covariates, supply one-row-per-person
  background data for the scored units and the posterior summaries will
  condition on the fitted population model.
- Intercept-only latent-regression fits (`population_formula = ~ 1`) can
  reconstruct that minimal scored-person table from the person IDs in
  `new_units`.
- For `JML` fits, the scoring layer remains a post hoc reference-prior
  approximation rather than a latent-regression fit.
- It returns posterior summaries, not deterministic future true values.
- `sample_mfrm_plausible_values()` exposes posterior draws under the same
  fitted scoring basis; the ordinary `MML` route is fixed-calibration, while
  active latent-regression fits use the fitted population model.
- Non-person facet levels in `new_units` must already exist in the fitted calibration.

## Prediction-aware bundle export

```r
bundle_pred <- export_mfrm_bundle(
  fit = toy_fit,
  population_prediction = pred_pop,
  unit_prediction = pred_units,
  plausible_values = pv_units,
  output_dir = tempdir(),
  prefix = "mfrmr_prediction_bundle",
  include = c("manifest", "predictions", "html"),
  overwrite = TRUE
)

bundle_pred$summary
```

Notes:

- `include = "predictions"` only writes prediction artifacts that you actually supply.
- Use `predict_mfrm_units()` and `sample_mfrm_plausible_values()` only with an
  existing fitted calibration. For latent-regression fits, keep the scoring
  `person_data` contract explicit when the fitted population model includes
  covariates rather than treating the scored outputs as ordinary
  fixed-calibration summaries.
- When a latent-regression fit is exported with `include = c("script", "html")`,
  the bundle writes a fit-level replay person-data sidecar for the replay
  script, while the HTML bundle exposes only an artifact index for that sidecar
  rather than embedding raw person-level rows.

## DIF / Bias screening simulation

```r
spec_sig <- build_mfrm_sim_spec(
  n_person = 50,
  n_rater = 4,
  n_criterion = 4,
  raters_per_person = 2,
  assignment = "rotating",
  group_levels = c("A", "B")
)

sig_eval <- evaluate_mfrm_signal_detection(
  n_person = c(30, 50, 80),
  n_rater = 4,
  n_criterion = 4,
  raters_per_person = 2,
  reps = 2,
  dif_effect = 0.8,
  bias_effect = -0.8,
  maxit = 15,
  sim_spec = spec_sig,
  seed = 123
)

s_sig <- summary(sig_eval)
s_sig$detection_summary
s_sig$ademp

plot(sig_eval, signal = "dif", metric = "power", x_var = "n_person")
plot(sig_eval, signal = "bias", metric = "false_positive", x_var = "n_person")
```

Notes:

- `DIFPower` is a conventional detection-power summary for the injected DIF target.
- `BiasScreenRate` and `BiasScreenFalsePositiveRate` summarize screening behavior from `estimate_bias()`.
- Bias-side `t`/`Prob.` values are screening metrics, not formal inferential p-values.

## Bundle export

```r
bundle <- export_mfrm_bundle(
  fit_bias,
  diagnostics = diag_bias,
  bias_results = bias_all,
  output_dir = tempdir(),
  prefix = "mfrmr_bundle",
  include = c("core_tables", "checklist", "manifest", "visual_summaries", "script", "html"),
  overwrite = TRUE
)

bundle$written_files

bundle_pred <- export_mfrm_bundle(
  toy_fit,
  output_dir = tempdir(),
  prefix = "mfrmr_prediction_bundle",
  include = c("manifest", "predictions", "html"),
  population_prediction = pred_pop,
  unit_prediction = pred_units,
  plausible_values = pv_units,
  overwrite = TRUE
)

bundle_pred$written_files

replay <- build_mfrm_replay_script(
  fit_bias,
  diagnostics = diag_bias,
  bias_results = bias_all,
  data_file = "your_data.csv"
)

replay$summary
```

## Anchoring and linking

```r
d1 <- load_mfrmr_data("study1")
d2 <- load_mfrmr_data("study2")
fit1 <- fit_mfrm(d1, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)
fit2 <- fit_mfrm(d2, "Person", c("Rater", "Criterion"), "Score", method = "JML", maxit = 25)

# Anchored calibration
res <- anchor_to_baseline(d2, fit1, "Person", c("Rater", "Criterion"), "Score")
summary(res)
res$drift

# Drift detection
drift <- detect_anchor_drift(list(Wave1 = fit1, Wave2 = fit2))
summary(drift)
plot_anchor_drift(drift, type = "drift")

# Screened linking chain
chain <- build_equating_chain(list(Form1 = fit1, Form2 = fit2))
summary(chain)
plot_anchor_drift(chain, type = "chain")
```

Notes:

- `detect_anchor_drift()` and `build_equating_chain()` remove the common-element link offset first, then report residual drift/link residuals.
- Treat `LinkSupportAdequate = FALSE` as a weak-link warning: at least one linking facet retained fewer than 5 common elements after screening.
- `build_equating_chain()` is a practical screened linking aid, not a full general-purpose equating framework.

## QC pipeline

```r
qc <- run_qc_pipeline(fit, threshold_profile = "standard")
qc$overall      # "Pass", "Warn", or "Fail"
qc$verdicts     # per-check verdicts
qc$recommendations

plot_qc_pipeline(qc, type = "traffic_light")
plot_qc_pipeline(qc, type = "detail")

# Threshold profiles: "strict", "standard", "lenient"
qc_strict <- run_qc_pipeline(fit, threshold_profile = "strict")
```

## Compatibility layer

Compatibility helpers are still available, but they are no longer the primary
route for new scripts.

- Use `run_mfrm_facets()` or `mfrmRFacets()` only when you need the one-shot wrapper.
- Use `build_fixed_reports()` and `facets_output_file_bundle()` only when a
  fixed-width or legacy export contract is required.
- For routine work, prefer package-native routes built from `fit_mfrm()`,
  `diagnose_mfrm()`, `reporting_checklist()`, and `build_apa_outputs()`.

For the full map, see
`help("mfrmr_compatibility_layer", package = "mfrmr")`.

External-software wording should stay conservative:

```r
chk <- reporting_checklist(fit, diagnostics = diag)
chk$software_scope
summary(chk)$software_scope
```

- `mfrmr native`: primary analysis surface.
- `FACETS`: compatibility-style wrappers and exports for handoff; results
  remain `mfrmr` estimates unless a separate external audit is performed.
- `ConQuest`: narrow external-table audit path for the documented latent-
  regression overlap; use scoped comparison wording.
- `SPSS`: CSV/data-frame/reporting handoff only; no native SPSS integration.

## Legacy-compatible one-shot wrapper

```r
run <- run_mfrm_facets(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "JML",
  model = "RSM"
)
summary(run)
plot(run, type = "fit", draw = FALSE)
```

## Public API map

Model and diagnostics:

- `fit_mfrm()`, `run_mfrm_facets()`, `mfrmRFacets()`
- `diagnose_mfrm()`, `analyze_residual_pca()`
- `estimate_bias()`, `bias_count_table()`

Differential functioning and model comparison:

- `analyze_dff()`, `analyze_dif()` (compatibility alias), `dif_interaction_table()`, `dif_report()`
- `compare_mfrm()`
- `compute_information()`, `plot_information()` for design-weighted precision curves

Anchoring and linking:

- `anchor_to_baseline()`, `detect_anchor_drift()`, `build_equating_chain()`
- `plot_anchor_drift()`

QC pipeline:

- `run_qc_pipeline()`, `plot_qc_pipeline()`

Table/report outputs:

- `specifications_report()`, `data_quality_report()`, `estimation_iteration_report()`
- `subset_connectivity_report()`, `facet_statistics_report()`
- `measurable_summary_table()`, `rating_scale_table()`
- `category_structure_report()`, `category_curves_report()`
- `unexpected_response_table()`, `unexpected_after_bias_table()`
- `fair_average_table()`, `displacement_table()`
- `interrater_agreement_table()`, `facets_chisq_table()`
- `facets_output_file_bundle()`, `facets_parity_report()` for compatibility-style output coverage checks
- `bias_interaction_report()`, `build_fixed_reports()`
- `apa_table()`, `build_apa_outputs()`, `build_visual_summaries()`

Output terminology:

- `ModelSE`: model-based standard error used for primary estimation summaries
- `RealSE`: fit-adjusted standard error, useful as a conservative companion
- `fair_average_table()` keeps historical display labels such as `Fair(M) Average`,
  and also exposes package-native aliases such as `AdjustedAverage`,
  `StandardizedAdjustedAverage`, `ModelBasedSE`, and `FitAdjustedSE`

Plots and dashboards:

- `plot_unexpected()`, `plot_fair_average()`, `plot_displacement()`
- `plot_interrater_agreement()`, `plot_facets_chisq()`
- `plot_bias_interaction()`, `plot_residual_pca()`, `plot_qc_dashboard()`
- `visual_reporting_template()` -- beginner-oriented figure placement,
  caption skeleton, results wording, interpretation, and overclaim-avoidance table
- `plot(fit, type = "ccc_surface", draw = FALSE)` -- 3D-ready category
  probability surface payload for exploratory downstream rendering
- `plot_bubble()` -- Rasch-convention bubble chart (Measure x Fit x SE)
- `plot_dif_heatmap()` -- differential-functioning heatmap across groups
- `plot_wright_unified()` -- unified Wright map across all facets
- `plot(fit, show_ci = TRUE)` -- approximate normal-interval whiskers on Wright map

Export and data utilities:

- `export_mfrm()` -- batch CSV export of all result tables
- `as.data.frame(fit)` -- tidy data.frame for `write.csv()` export
- `describe_mfrm_data()`, `audit_mfrm_anchors()`, `make_anchor_table()`
- `mfrm_threshold_profiles()`, `list_mfrmr_data()`, `load_mfrmr_data()`

Some legacy FACETS-style numbered names remain unexported.
`facets_parity_report()` keeps its historical function name, but the intended
use is package-output coverage checking, not a FACETS execution or
external software-equivalence claim.

## FACETS reference mapping

See:

- `inst/references/FACETS_manual_mapping.md`

## Packaged synthetic datasets

Installed at `system.file("extdata", package = "mfrmr")`:

- `eckes_jin_2021_study1_sim.csv`
- `eckes_jin_2021_study2_sim.csv`
- `eckes_jin_2021_combined_sim.csv`
- `eckes_jin_2021_study1_itercal_sim.csv`
- `eckes_jin_2021_study2_itercal_sim.csv`
- `eckes_jin_2021_combined_itercal_sim.csv`

The same datasets are also packaged in `data/` and can be loaded with:

```r
data("ej2021_study1", package = "mfrmr")
# or
df <- load_mfrmr_data("study1")
```

Current packaged dataset sizes:

- `study1`: 1842 rows, 307 persons, 18 raters, 3 criteria
- `study2`: 3287 rows, 206 persons, 12 raters, 9 criteria
- `combined`: 5129 rows, 307 persons, 18 raters, 12 criteria
- `study1_itercal`: 1842 rows, 307 persons, 18 raters, 3 criteria
- `study2_itercal`: 3341 rows, 206 persons, 12 raters, 9 criteria
- `combined_itercal`: 5183 rows, 307 persons, 18 raters, 12 criteria

## Citation

```r
citation("mfrmr")
```

## Acknowledgements

`mfrmr` has benefited from discussion and methodological input from
[Dr. Atsushi Mizumoto](https://mizumot.com/) and
[Dr. Taichi Yamashita](https://kugakujo.kansai-u.ac.jp/html/100000882_en.html).

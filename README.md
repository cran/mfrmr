# mfrmr

[![GitHub](https://img.shields.io/badge/GitHub-mfrmr-181717?logo=github)](https://github.com/Ryuya-dot-com/R_package_mfrmr)
[![R-CMD-check](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/pkgdown.yaml)
[![test-coverage](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/Ryuya-dot-com/R_package_mfrmr/actions/workflows/test-coverage.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Native R package for many-facet Rasch model (MFRM) estimation, diagnostics,
and reporting workflows.

## What this package is for

`mfrmr` is designed around four package-native routes:

- Estimation and diagnostics: `fit_mfrm()` -> `diagnose_mfrm()`
- Reporting and manuscript preparation: `reporting_checklist()` -> `build_apa_outputs()`
- Linking, anchors, drift, and DFF: `subset_connectivity_report()` -> `anchor_to_baseline()` / `analyze_dff()`
- Legacy-compatible export when required: `run_mfrm_facets()` and related compatibility helpers

If you want the shortest possible recommendation:

- Final estimation: prefer `method = "MML"`
- Fast exploratory pass: use `method = "JML"`
- First visual screen: `plot_qc_dashboard(..., preset = "publication")`
- First reporting screen: `reporting_checklist()`

## Features

- Flexible facet count (`facets = c(...)`)
- Estimation methods: `MML` (default) and `JML` (`JMLE` internally)
- Models: `RSM`, `PCM` (`PCM` uses step-facet-specific thresholds on a shared observed score scale)
- Legacy-compatible one-shot wrapper (`run_mfrm_facets()`, alias `mfrmRFacets()`)
- Bias/interaction iterative estimation with compatibility exports
- App-style multi-pair bias orchestration (`estimate_all_bias()`)
- Optional fixed-width text reports for console/log audits
- APA-style narrative output helpers (`build_apa_outputs()`)
- Reporting checklist and facet dashboard helpers (`reporting_checklist()`, `facet_quality_dashboard()`)
- Bundle export, manifest, and replay-script helpers (`export_mfrm_bundle()`, `build_mfrm_manifest()`, `build_mfrm_replay_script()`)
- Facet equivalence analysis with TOST / ROPE summaries (`analyze_facet_equivalence()`)
- Visual warning and summary maps (`build_visual_summaries()`)
- Exploratory residual PCA summaries (`overall` / `facet` / `both`)
- Differential facet functioning (DFF) analysis with method-appropriate classification (`analyze_dff()`, `analyze_dif()`, `dif_report()`)
- Model comparison with information criteria on a common likelihood basis (`compare_mfrm()`)
- RSM design-weighted precision curves (`compute_information()`, `plot_information()`)
- Unified Wright map across facets (`plot_wright_unified()`)
- Anchoring and linking across administrations (`anchor_to_baseline()`, `detect_anchor_drift()`, `build_equating_chain()`)
- Automated 10-check QC pipeline (`run_qc_pipeline()`, `plot_qc_pipeline()`)
- Simulation-based design planning (`simulate_mfrm_data()`, `evaluate_mfrm_design()`)
- Fit-derived semi-parametric / plasmode-style simulation specs via empirical latent supports, resampled rater profiles, or observed response skeleton reuse with optional `Weight` carry-over and optional person-level `Group` carry-over from `source_data` (`extract_mfrm_sim_spec(..., latent_distribution = "empirical", assignment = "resampled" / "skeleton")`)
- Scenario-level population forecasting from explicit simulation specs (`build_mfrm_sim_spec()`, `extract_mfrm_sim_spec()`, `predict_mfrm_population()`)
- Fixed-calibration posterior scoring for future or partially observed units (`predict_mfrm_units()`)
- Approximate plausible-value draws for future or partially observed units under a fixed MML calibration (`sample_mfrm_plausible_values()`)
- Data descriptives and anchor audit helpers (`describe_mfrm_data()`, `audit_mfrm_anchors()`)
- Anchor export for linking workflows (`make_anchor_table()`)

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

Companion vignettes:

- `vignette("mfrmr-workflow", package = "mfrmr")`
- `vignette("mfrmr-visual-diagnostics", package = "mfrmr")`
- `vignette("mfrmr-reporting-and-apa", package = "mfrmr")`
- `vignette("mfrmr-linking-and-dff", package = "mfrmr")`

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
10. Compatibility audit when needed: `facets_parity_report()`
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

## Start here

If you are new to the package, these are the three shortest useful routes.

Shared setup used by the snippets below:

```r
library(mfrmr)
toy <- load_mfrmr_data("example_core")
```

### 1. Quick first pass

```r
fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                method = "MML", model = "RSM", quad_points = 7)
diag <- diagnose_mfrm(fit, residual_pca = "none")
summary(diag)
plot_qc_dashboard(fit, diagnostics = diag, preset = "publication")
```

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

fit_fast <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                     method = "JML", model = "RSM", maxit = 50)
fit_final <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                      method = "MML", model = "RSM", quad_points = 15)

diag_final <- diagnose_mfrm(fit_final, residual_pca = "none")
precision_audit_report(fit_final, diagnostics = diag_final)
```

## Help-page navigation

Core analysis help pages include practical sections such as:

- `Interpreting output`
- `Typical workflow`

Recommended entry points:

- `?mfrmr-package` (package overview)
- `?fit_mfrm`, `?diagnose_mfrm`, `?run_mfrm_facets`
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

- Start with `method = "JML"` when you want a quick exploratory fit.
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
# Wright map with shared targeting view
plot(fit, type = "wright", preset = "publication", show_ci = TRUE)

# Pathway map with dominant-category strips
plot(fit, type = "pathway", preset = "publication")

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

- `predict_mfrm_units()` scores future or partially observed persons under a fixed MML calibration.
- It returns posterior summaries, not deterministic future true values.
- `sample_mfrm_plausible_values()` exposes fixed-calibration posterior draws as approximate plausible-value summaries.
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
- Use `predict_mfrm_units()` and `sample_mfrm_plausible_values()` only with an existing MML calibration.

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
- `facets_output_file_bundle()`, `facets_parity_report()`
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
- `plot_bubble()` -- Rasch-convention bubble chart (Measure x Fit x SE)
- `plot_dif_heatmap()` -- differential-functioning heatmap across groups
- `plot_wright_unified()` -- unified Wright map across all facets
- `plot(fit, show_ci = TRUE)` -- approximate normal-interval whiskers on Wright map

Export and data utilities:

- `export_mfrm()` -- batch CSV export of all result tables
- `as.data.frame(fit)` -- tidy data.frame for `write.csv()` export
- `describe_mfrm_data()`, `audit_mfrm_anchors()`, `make_anchor_table()`
- `mfrm_threshold_profiles()`, `list_mfrmr_data()`, `load_mfrmr_data()`

Legacy FACETS-style numbered names are internal and not exported.

## FACETS reference mapping

See:

- `inst/references/FACETS_manual_mapping.md`
- `inst/references/CODE_READING_GUIDE.md` (for developers/readers)

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

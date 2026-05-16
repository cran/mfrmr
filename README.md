# mfrmr

[![GitHub](https://img.shields.io/badge/GitHub-mfrmr-181717?logo=github)](https://github.com/Ryuya-dot-com/mfrmr)
[![R-CMD-check](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/pkgdown.yaml)
[![test-coverage](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/test-coverage.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Native R package for many-facet ordered-response measurement models: the
Rasch-family `RSM` / `PCM` route, plus the package's bounded `GPCM` extension
where explicitly documented.

## Start here first

If you are new to `mfrmr`, use this route first and ignore the longer feature
lists below until it works end to end.

- Fit with `method = "MML"`
- Diagnose with `diagnostic_mode = "both"` for `RSM` / `PCM`; for bounded
  `GPCM`, use the direct diagnostic route and read the caveats in
  `gpcm_capability_matrix()`
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
  model = "RSM"
  # quad_points defaults to 31 (publication tier); set 7 or 15 for
  # exploratory iteration.
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

- reporting: `build_apa_outputs()` / `apa_table()` for the full `RSM` /
  `PCM` manuscript route; `build_summary_table_bundle()` /
  `export_summary_appendix()` for direct appendix handoff
- misfit case review: `build_misfit_casebook()`
- weighting review: fit both an `RSM`/`PCM` and a bounded `GPCM` model with `fit_mfrm()`, then pass the two fits to `build_weighting_review(rsm_fit, gpcm_fit)`. `compare_mfrm()` is a complementary information-criterion summary over the same pair.
- confirmatory facet interaction review: fit an `RSM`/`PCM` model with
  explicit `facet_interactions = "FacetA:FacetB"`, inspect
  `interaction_effect_table(fit)`, and compare it to the additive fit on the
  same likelihood basis.
- strict follow-up: `plot_marginal_fit()` / `plot_marginal_pairwise()`
- operational linking review: `review_mfrm_anchors()` -> `detect_anchor_drift()` ->
  `build_linking_review()` for `RSM` / `PCM`
- linking/design: `subset_connectivity_report()`

## What this package is for

`mfrmr` is designed around five package-native routes:

- Estimation and diagnostics: `fit_mfrm()` -> `diagnose_mfrm()`
- Reporting and manuscript preparation: `reporting_checklist()` ->
  `build_apa_outputs()` for the full `RSM` / `PCM` manuscript route, or
  `build_summary_table_bundle()` -> `export_summary_appendix()` for direct
  appendix handoff
- Misfit case review: `build_misfit_casebook()` -> `casebook$group_view_index` /
  `casebook$group_views` -> source-specific follow-up plots
- Linking, anchors, drift, and DFF:
  `review_mfrm_anchors()` / `detect_anchor_drift()` -> `build_linking_review()`
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
- First weighting-policy screen: `build_weighting_review()`
- First operational linking screen: `build_linking_review()` for `RSM` /
  `PCM`; for bounded `GPCM`, keep anchor/drift helpers as direct exploratory
  support

## Minimum input contract

`mfrmr` expects long-format rating data: one row per observed rating.

- Required columns:
  - one person column
  - one ordered score column
  - one or more non-person facet columns supplied in `facets = c(...)`
- Score-column rules:
  - scores should be ordered integer category codes such as `0/1`, `1/2`, or `1:5`
  - binary two-category scores are supported under the same ordered-response
    interface
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
- equivalence and review helpers such as `analyze_facet_equivalence()`, `describe_mfrm_data()`, and `review_mfrm_anchors()`

Design-adequacy review and partial pooling:

- hierarchical-structure and sample-adequacy review with `detect_facet_nesting()`, `facet_small_sample_review()`, `compute_facet_icc()`, `compute_facet_design_effect()`, and the combined `analyze_hierarchical_structure()`
- empirical-Bayes / James-Stein shrinkage for small-N facets via `fit_mfrm(..., facet_shrinkage = "empirical_bayes")` or post-hoc `apply_empirical_bayes_shrinkage()`, with `shrinkage_report()` as the accessor
- missing-code pre-processing through `fit_mfrm(..., missing_codes = TRUE)` (FACETS / SPSS / SAS sentinels such as `99`, `999`, `-1`, `"N/A"`, `""` converted to `NA`) or the standalone `recode_missing_codes()` helper
- APA output adapters `as_kable.apa_table()` and `as_flextable.apa_table()` for RMarkdown / Quarto / Word / PowerPoint handoffs

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
review <- review_conquest_overlap(bundle, normalized)
summary(review)$summary
review$attention_items
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
  are exploratory; `evaluate_mfrm_recovery()` checks direct parameter recovery
  rather than design operating characteristics; `reporting_checklist()`,
  `build_summary_table_bundle()`, and `export_summary_appendix()` route only
  the direct table/plot path. `fair_average_table()` and `estimate_bias()` use
  the slope-aware element-conditional GPCM kernel. For fair averages, the
  historical SE columns remain scaled facet-measure SEs; use
  `fair_average_table(fair_se = TRUE)` to request structural delta-method
  fair-average SEs for non-person rows when the MML Hessian is available. For
  bias screening, the SE / `t` / `Prob.` columns are
  conditional plug-in screening quantities, and bounded-GPCM rows also carry
  conditional profile-likelihood columns for follow-up review.
- Not supported in this release: FACETS-style score-side exports, the APA
  writer, fit-based report/export bundles, QC pass/fail pipelines, linking
  synthesis, planning / forecasting, posterior predictive computation, and
  `MCMC`.

The unsupported helpers depend on FACETS-style score-side, narrative-export,
or planning assumptions that are validated for the Rasch-family route but not
yet for bounded `GPCM`.

For release review, the optional script
`system.file("validation", "recovery-validation.R", package = "mfrmr")`
defines core `RSM` / `PCM` / bounded-`GPCM` recovery cases, an extended
latent-regression case, structured release-review steps, and CSV/RDS/Markdown
summaries. It is intentionally separate from routine tests because the useful
settings are long-running Monte Carlo checks. The summary separates recovery
metric status from uncertainty status so unavailable coverage columns do not
look like failed parameter recovery. Printing the validation object or calling
`summary(validation)` shows the release-level status first.

For direct recovery checks, `plot(evaluate_mfrm_recovery(...), ...)` shows
recovery summaries, row-level errors, truth-estimate scatter, and replication
status. After `assess_mfrm_recovery()`, use
`plot(recovery_review, type = "status")` for checklist status counts and
`plot(recovery_review, type = "metrics", metric = "rmse")` for the
parameter-group metric review. The recommended reading order is:
`summary(recovery_review)`, then the status plot, then the metric plot, and
only then the row-level recovery table for the parameter groups that need
follow-up. The `draw = FALSE` plot data include `reading_order` and
`guidance` fields for this handoff.

Read the validation outputs in this order:

- `topline_release_decision`: the release-level recovery conclusion. Its
  `ReleaseRecoveryStatus` uses recovery metrics, convergence, and Monte Carlo
  precision as the primary evidence.
- `release_decision_table`: the same decision by validation case, with a short
  interpretation and any uncertainty limitation.
- `domain_decision_table`: the diagnostic split among recovery metrics,
  uncertainty, Monte Carlo precision, and the broader overall status.

In particular, do not treat `OverallStatus = "review"` as a release-level
recovery failure by itself. In the validation bundle, `UncertaintyStatus =
"review"` can mean that SE/coverage evidence is intentionally reported as a
separate limitation while recovery metrics remain acceptable.

For a source-grounded release review plan, read the packaged evidence map and
its structured checklist:

```r
file.show(system.file(
  "validation", "release-evidence-map-0.2.0.md",
  package = "mfrmr"
))

read.csv(system.file(
  "validation", "release-evidence-checklist-0.2.0.csv",
  package = "mfrmr"
))

file.show(system.file(
  "validation", "external-parameter-recovery-simulation-0.2.0.md",
  package = "mfrmr"
))
```

It links the 0.2.0 release checks to the ordered-response model literature,
FACETS/Winsteps fit conventions, and ADEMP-style simulation-study reporting.
The checklist classifies each item as a release blocker, caveat-managed item,
or post-release roadmap item.

The external parameter-recovery summary records a separate common-data
simulation workflow. It supports the distinction between recovery checks,
cross-engine agreement, and design endorsement: sparse stress designs can
converge and agree across engines while still showing recovery, coverage,
precision, or role-bias risk. The large generated datasets and engine outputs
are not bundled with the package; the validation bundle includes a sourceable
review helper for re-reading a local `Parameter_Recovery_Simulation` output
directory, checking expected CSV schemas, and recording file fingerprints when
that external workflow is refreshed.

## Equal weighting and when to prefer the Rasch-family route

`mfrmr` treats `RSM` / `PCM` as the package's equal-weighting reference
models. In that Rasch-family route, category discrimination is fixed, so the
operational scoring contract does not let the psychometric model reweight some
item-facet combinations more heavily than others.

Bounded `GPCM` serves a different purpose. It allows estimated slopes, so some
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

## Model selection guide and report wording

Use the model argument to match the score interpretation first, then use fit
statistics and diagnostics as checks on that interpretation.

| Choose | When it is the right starting point | Report wording |
|---|---|---|
| `RSM` | The rubric is intended to share the same category thresholds across items, criteria, or other step-facet levels. | "We fit a many-facet rating-scale Rasch model, treating category thresholds as common across the step facet." |
| `PCM` | Category thresholds may differ by item or criterion, but equal contribution of rating events remains part of the scoring argument. | "We fit a many-facet partial-credit Rasch model, allowing step thresholds to vary by the designated step facet." |
| bounded `GPCM` | You explicitly want a slope-aware sensitivity model and can defend discrimination-based reweighting. | "We fit a bounded generalized partial-credit many-facet model as a slope-aware sensitivity analysis." |

Avoid these shortcuts:

- do not describe the whole package surface as "a Rasch-only MFRM" now that
  bounded `GPCM` is implemented
- do not write that bounded `GPCM` is better for operational scoring solely
  because `AIC`, `BIC`, or log-likelihood improves
- do not use FACETS-style score-side, APA writer, QC pass/fail, or linking-
  synthesis language for bounded `GPCM` unless `gpcm_capability_matrix()`
  marks that route as supported

In a manuscript, a defensible model-choice sentence is:

> We treated `RSM`/`PCM` as the equal-weighting operational reference and used
> bounded `GPCM` to inspect whether allowing discrimination-based reweighting
> changed the substantive conclusions.

After fitting candidate models, use `build_model_choice_review()` to keep the
same guidance attached to the actual fit objects:

```r
review <- build_model_choice_review(RSM = fit_rsm, GPCM = fit_gpcm)
summary(review)

# Add the detailed reweighting review when an RSM/PCM reference and bounded
# GPCM sensitivity fit were estimated on the same response data.
review <- build_model_choice_review(RSM = fit_rsm, GPCM = fit_gpcm,
                                    run_weighting_review = TRUE)
```

## Documentation map

The README is only the shortest map. The package now has guide-style help pages
for the main workflows.

- Workflow map:
  `help("mfrmr_workflow_methods", package = "mfrmr")`
- Visual diagnostics map:
  `help("mfrmr_visual_diagnostics", package = "mfrmr")`
- Reports and tables map:
  `help("mfrmr_reports_and_tables", package = "mfrmr")`
- Output helper guide:
  `mfrmr_output_guide()`
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
- `vignette("mfrmr-gpcm-scope", package = "mfrmr")`
- `vignette("mfrmr-facets-migration", package = "mfrmr")`

A two-page landscape cheatsheet of the public API ships at
`system.file("cheatsheet", "mfrmr-cheatsheet.pdf", package = "mfrmr")`
(pre-rendered) and `system.file("cheatsheet", "mfrmr-cheatsheet.Rmd",
package = "mfrmr")` (source). Open the PDF directly for a quick
printable reference, or knit the `.Rmd` with `rmarkdown::render()`
when you want a customised version.

## Installation

```r
# GitHub
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Ryuya-dot-com/mfrmr", build_vignettes = TRUE)

# CRAN (when available)
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
                    +--> interaction_effect_table()
                    +--> analyze_dff()
                    +--> compare_mfrm()
                    +--> run_qc_pipeline()
                    +--> anchor_to_baseline() / detect_anchor_drift()
```

1. Fit model: `fit_mfrm()`
2. Diagnostics: `diagnose_mfrm()`
3. Optional residual-structure screen: `analyze_residual_pca()`
4. Optional interaction bias: `estimate_bias()`
5. Optional model-estimated facet interactions:
   `interaction_effect_table()`
6. Differential-functioning analysis: `analyze_dff()`, `dif_report()`
7. Model comparison: `compare_mfrm()`
8. Reporting: `apa_table()`, `build_apa_outputs()`, `build_visual_summaries()`
9. Quality control: `run_qc_pipeline()`
10. Anchoring & linking: `anchor_to_baseline()`, `detect_anchor_drift()`, `build_equating_chain()`
11. FACETS output-contract review when needed: `facets_output_contract_review()`;
    this checks package output contracts, not external FACETS numerical
    equivalence
12. Reproducible inspection: `summary()` and `plot(..., draw = FALSE)`

Dimensionality wording is deliberately conservative. Residual PCA and Q3-style
local-dependence screens are exploratory follow-up evidence, not standalone
proofs that unidimensionality has been established and not implementations of
DIMTEST/UNIDIM. For MFRM manuscripts, combine global residual fit, element fit,
residual PCA, and local-dependence checks, and use limited wording such as
"evidence consistent with essential unidimensionality under the specified facet
structure."

## Choose a route

Use the route that matches the question you are trying to answer.

| Question | Recommended route |
|---|---|
| Can I fit the model and get a first-pass diagnosis quickly? | `fit_mfrm()` -> `diagnose_mfrm()` -> `plot_qc_dashboard()` |
| Which reporting elements are draft-complete, and with what caveats? | `diagnose_mfrm()` -> `precision_review_report()` -> `reporting_checklist()` |
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

### 4. Hierarchical structure and sample-adequacy review

Use this when rater counts are small, raters may be nested in schools
or regions, or a reviewer asks for ICC / design-effect evidence that
the additive fixed-effects many-facet model cannot partition out on its own.

```r
review <- facet_small_sample_review(fit)
review$facet_summary         # worst level per facet + SampleCategory
summary(review)              # counts of sparse / marginal / standard / strong

nest <- detect_facet_nesting(toy, c("Rater", "Criterion"))
plot(nest)                   # nesting index heatmap

# Combined bundle (ICC uses lme4, connectivity uses igraph, both Suggests):
h <- analyze_hierarchical_structure(toy, c("Rater", "Criterion"), score = "Score",
                                    person = "Person")
summary(h)
```

`reporting_checklist(fit, hierarchical_structure = h)` then marks the
"Hierarchical structure review" item ready.

### 5. Empirical-Bayes shrinkage for small-N facets

When a facet has 3-10 levels, the fixed-effects many-facet model retains wide
per-level SEs. Empirical-Bayes partial pooling (Efron & Morris, 1973)
dominates the MLE under squared-error loss whenever `K >= 3`.

```r
# Integrated path: shrinkage applied as part of the fit.
fit_eb <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                   method = "MML", quad_points = 15,
                   facet_shrinkage = "empirical_bayes")
shrinkage_report(fit_eb)
plot(fit_eb, type = "shrinkage", show_ci = TRUE)

# Post-hoc path: apply to an existing fit.
fit_post <- apply_empirical_bayes_shrinkage(fit)
head(fit_post$facets$others[, c("Facet", "Level", "Estimate",
                                 "ShrunkEstimate", "ShrinkageFactor")])
```

### 6. Missing-code pre-processing

`fit_mfrm(..., missing_codes = TRUE)` converts the default
FACETS / SPSS / SAS sentinels (`"99"`, `"999"`, `"-1"`, `"N"`, `"NA"`,
`"n/a"`, `"."`, `""`) to `NA` on the `person`, `facets`, and `score`
columns before estimation. Replacement counts are kept in
`fit$prep$missing_recoding` and surfaced by
`build_mfrm_manifest()$missing_recoding`. The default
(`missing_codes = NULL`) is strictly backward-compatible.

```r
fit <- fit_mfrm(
  dirty_data, "Person", c("Rater", "Criterion"), "Score",
  missing_codes = TRUE           # or supply a custom character vector
)
fit$prep$missing_recoding
```

A standalone `recode_missing_codes()` helper is exported for users who
prefer to recode before calling `fit_mfrm()`.

## Estimation choices

The package treats `MML` and `JML` differently on purpose.

- `MML` is the default and the preferred route for final estimation.
- `JML` is supported as a fast exploratory route.
- Downstream precision summaries distinguish `model_based`, `hybrid`, and `exploratory` tiers.
- Use `precision_review_report()` when you need to decide how strongly to phrase SE, CI, or reliability claims.

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

precision_review_report(fit_final, diagnostics = diag_final)
```

## Mathematical note for expert users

Full marginal-likelihood and strict-marginal derivations, along with the
literature positioning (Bock & Aitkin, 1981; Linacre, 1989; Eckes, 2005;
Orlando & Thissen, 2000; Haberman & Sinharay, 2013; Sinharay & Monroe,
2025), are collected in the dedicated vignette:

```r
vignette("mfrmr-mml-and-marginal-fit", package = "mfrmr")
```

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
everything from scratch. The canonical list is kept up to date in
`summary(fit)` under "Next actions"; the items below are a short
orientation pointer.

- `fit`: the fitted model object returned by `fit_mfrm()`
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

A task-oriented index of the plotting surface lives at
`help("mfrmr_visual_diagnostics", package = "mfrmr")`, and worked
publication examples are collected in
`vignette("mfrmr-visual-diagnostics", package = "mfrmr")`.
The common starter patterns are:

```r
plot(fit, type = "wright", preset = "publication", show_ci = TRUE)
plot(fit, type = "pathway", preset = "publication")
plot(fit, type = "ccc", preset = "publication")
plot_qc_dashboard(fit, diagnostics = diag, preset = "publication")
```

A second-wave teaching / drift / agreement layer ships for follow-up
inspection; it is not a default reporting figure set:

```r
plot_guttman_scalogram(fit, diagnostics = diag)       # teaching ordering view
plot_residual_qq(fit, diagnostics = diag)             # residual tail follow-up
plot_rater_agreement_heatmap(fit, diagnostics = diag) # compact pairwise agreement
plot_rater_trajectory(list(T1 = fit_a, T2 = fit_b))   # requires anchor-linked waves
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
plot_dif_summary(dff)
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
plot_dif_summary(dff)

# Optional display controls for review meetings or appendices
plot_dif_heatmap(dff, metric = "t", flag_threshold = 2,
                 show_values = FALSE, scale_limit = 3)
plot_dif_summary(dff, ci_level = 0.90,
                 effect_thresholds = c(screen = 0.5))
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
- Residual-method classifications are screening labels, not ETS A/B/C severity categories.
- Check `ScaleLinkStatus`, `ContrastComparable`, and the reported classification system before treating a contrast as a strong interpretive claim.

## Model-estimated facet interactions

For confirmatory interaction hypotheses, `fit_mfrm()` can estimate explicit
two-way non-person facet interactions in the model likelihood.

```r
fit_add <- fit_mfrm(df, "Person", c("Rater", "Criterion"), "Score",
                    method = "MML", model = "RSM")

fit_rxcrit <- fit_mfrm(df, "Person", c("Rater", "Criterion"), "Score",
                       method = "MML", model = "RSM",
                       facet_interactions = "Rater:Criterion")

interaction_effect_table(fit_rxcrit)
compare_mfrm(Additive = fit_add, RaterCriterion = fit_rxcrit, nested = TRUE)
```

Rules for interpretation:

- Name the facet pair explicitly. The package does not add all possible
  interactions automatically.
- The current scope is two-way interactions between non-person facets for
  `RSM` and `PCM`; GPCM, person-involving, higher-order, and random-effect
  interaction terms are deferred.
- The interaction matrix uses zero marginal sums, so each estimate is a
  deviation from the additive main-effects MFRM. Positive values indicate
  higher-than-expected scores for that facet-level combination; negative
  values indicate lower-than-expected scores.
- `interaction_effect_table()` reports model-estimated fixed effects.
  `estimate_bias()` and `estimate_all_bias()` remain residual screening tools
  for exploratory bias review.
- Sparse cells matter. Use `min_obs_per_interaction` and inspect the `Sparse`
  column before reporting substantive interaction claims.

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
  maxit = 30,
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
- `evaluate_mfrm_design()` is a Monte Carlo design-evaluation helper. It can show how separation, reliability, strata, RMSE, and fit-screen rates change as facet counts vary; use `mfrm_generalizability()` plus `mfrm_d_study()` for observed G-study components and analytic D-study projections.

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
  maxit = 30,
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
  maxit = 30,
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
chk$facets_positioning
chk$software_scope
summary(chk)$software_scope
```

- `mfrmr native`: primary analysis surface.
- `FACETS`: FACETS-style reporting and handoff surfaces; results remain
  `mfrmr` estimates unless external FACETS output is supplied for explicit
  comparison.
- `ConQuest`: narrow external-table review path for the documented latent-
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

The full exported function index (with categories such as *Model and
diagnostics*, *Bias and DFF*, *Anchoring and linking*, *Reporting and
APA*, *Plots and dashboards*, *Simulation and design*, and *Export
utilities*) is generated from roxygen. Within R the same grouping is
available through the topic help pages
`?mfrmr_workflow_methods`, `?mfrmr_visual_diagnostics`,
`?mfrmr_reports_and_tables`, `?mfrmr_reporting_and_apa`,
`?mfrmr_linking_and_dff`, and `?mfrmr_compatibility_layer`.

Output-terminology note: `ModelSE` is the model-based standard error
used for primary summaries; `RealSE` is the fit-adjusted companion.
`fair_average_table()` keeps the historical display labels
(`Fair(M) Average`, `Fair(Z) Average`) alongside package-native
aliases `AdjustedAverage`, `StandardizedAdjustedAverage`,
`ModelBasedSE`, and `FitAdjustedSE`.

Reliability terminology note: `diagnostics$reliability` reports
Rasch/FACETS-style separation, strata, and separation reliability. These
indices answer whether persons, raters, criteria, or other facet elements are
distinguishable on the fitted logit scale. They are not intra-class
correlations. Use `compute_facet_icc()` only when you want a complementary
random-effects variance-share summary on the observed-score scale; for
non-person facets such as raters, a large ICC is systematic facet variance,
not better reliability.

Scope note: `mfrmr` does not estimate latent-class mixture models or
response-time / careless-rating adjustments. Use person fit, residual
matrices, Q3-style local-dependence screens, rater drift, and DFF diagnostics
as screening evidence, not as substitutes for an explicit mixture or
response-time model.

## FACETS reference mapping

A reference table mapping FACETS-program output tables (Table 1, Table 5,
Table 7, ...) to the `mfrmr` helper functions that produce substantively
corresponding or adjacent package-native reports ships with the installed
package. Open it with:

```r
file.show(system.file("references", "FACETS_manual_mapping.md", package = "mfrmr"))
```

The mapping is a package-output contract reference, not evidence that
FACETS was executed or that numerical FACETS equivalence has been
established for any given fit. The intended workflow is to estimate and
report from `mfrmr` objects, then use FACETS-style routes only for transition,
handoff, or explicit external-table review.

## Packaged synthetic datasets

Lazy-loaded under `data/` and accessed either by name or via the
canonical loader:

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

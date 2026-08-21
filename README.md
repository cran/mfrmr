# mfrmr <img src="man/figures/logo.png" align="right" height="160" alt="mfrmr hex logo" />

[![GitHub](https://img.shields.io/badge/GitHub-mfrmr-181717?logo=github)](https://github.com/Ryuya-dot-com/mfrmr)
[![R-CMD-check](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/pkgdown.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`mfrmr` fits unidimensional many-facet ordered-response models in R.
It supports rating-scale (`RSM`) and partial-credit (`PCM`) models, together
with a documented bounded `GPCM` extension. A facet can represent a rater,
item, task, criterion, form, occasion, or another observed role that affects
an ordered score.

For GPCM, *bounded* refers to the documented model and workflow scope; it does
not mean finite parameter box constraints. Unsupported combinations and
inference states are reported explicitly rather than silently treated as
ordinary estimates.

The package extends Rasch-family RSM/PCM work with MML, modern diagnostics,
reproducibility, network review, and reporting support. It is not a general
FACETS replacement: each `fit_mfrm()` call uses one response-model family and
one observed score scale, and the current public API does not provide mixed
response families, multiple independent rating scales, general threshold
anchoring, or an imported versioned frozen-calibration bundle for operational
scoring. Posterior scoring from an existing fitted object is a separate
analysis route.

The recommended workflow is:

```text
long-format data + describe_mfrm_data()
  -> fit_mfrm(method = "MML")
  -> print(fit) / summary(fit)$decision
  -> if ready: summary(profile = "facets")
     if review/blocked: follow summary(fit)$decision$NextAction
  -> native Wright map with uncertainty
  -> focused diagnostics
  -> report and export
```

The native Wright map is the required first fitted-scale figure in this
workflow. A separate FACETS-style renderer is available when a familiar
asterisk ruler and labelled category transitions are useful.

The three basic fitted-object methods have deliberately different jobs:
`print(fit)` is a compact triage view, `summary(fit)` is the canonical
structured review surface, and `plot(fit, draw = FALSE)` returns reusable plot
data. The printed `Decision` and structured `summary(fit)$decision` are the
common first-read for RSM, PCM, and bounded GPCM. All three methods carry the
same model/readiness basis; fitted-model plots also
retain a `scale_contract` table so the latent-coordinate and discrimination
scales are not inferred from axis labels alone.

Package website: <https://ryuya-dot-com.github.io/mfrmr/>

Source code: <https://github.com/Ryuya-dot-com/mfrmr>

Questions and bug reports:
<https://github.com/Ryuya-dot-com/mfrmr/issues>

## Installation

Install the CRAN package with:

```r
install.packages("mfrmr")
```

Install the GitHub version with:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "Ryuya-dot-com/mfrmr",
  build_vignettes = TRUE
)
```

Load the package, inspect the six-step beginner route, and list the installed
guides:

```r
library(mfrmr)
mfrmr_output_guide("beginner")[, c(
  "Question", "MainFunction", "NextStep"
)]
browseVignettes("mfrmr")
```

For the shortest end-to-end explanation, open
`help("mfrmr_workflow_methods", package = "mfrmr")` or
`vignette("mfrmr-workflow", package = "mfrmr")`. Use
`help("mfrmr_visual_diagnostics", package = "mfrmr")` when choosing a figure
and `help("mfrmr_reporting_and_apa", package = "mfrmr")` when moving from a
reviewed fit to tables and manuscript-draft output.

## Data format

`fit_mfrm()` expects one row per observed rating event. At minimum, the data
need:

- one person identifier;
- one or more non-person facet columns;
- one ordered numeric score column.

For example:

```r
head(data.frame(
  Person = c("P01", "P01", "P02", "P02"),
  Rater = c("R1", "R2", "R1", "R2"),
  Criterion = c("Content", "Content", "Content", "Content"),
  Score = c(3, 2, 2, 3)
))
```

Facet identifiers may be character or factor columns. Scores must represent
ordered categories. Use `NA` for missing responses, or document and recode
special missing-value codes before fitting. The design need not be fully
crossed, but it must contain enough links among persons and facet levels to
support the intended comparisons.

### Response type and frequencies

The current response likelihood is ordered categorical. Binary scores are the
two-category special case, and RSM, PCM, and bounded GPCM cover ordered
polytomous scores. Although their category probabilities form a vector that
sums to one, they are not unordered nominal-response or multinomial-logit
models. Poisson, negative-binomial, and grouped binomial-trial count responses
are also outside the current `fit_mfrm()` scope. Integer counts supplied as
`Score` are interpreted as ordered category codes, not as Poisson or related
counts.

A positive numeric `weight` weights one row's conditional ordered-category
likelihood and can represent a defensible row-replication weight. It is not a
general collapsed-person frequency table. Under MML, powering responses within
one Person's conditional pattern is not equivalent to replicating a complete
Person response pattern after marginalization. It also does not change the
response family or model dependence among repeated ratings. Non-unit
observation-weight fits are excluded from the common MML information-criterion
panel. [FACETS has separate `Bn` binomial-trial
and `P` Poisson response models](https://www.winsteps.com/facetman64/models.htm);
those are not reproduced by mfrmr's binary ordered-score route.

Each row should represent a distinguishable rating event. Exact duplicate
Person-by-facet combinations are retained but trigger a warning and a Data
review state because the package does not model within-cell dependence. For a
legitimate re-rating or replicated scoring study, include an event or occasion
facet that distinguishes the observations before fitting.

For conventional score sentinels such as `99`, `-1`, `N`, or `.`, set
`missing_codes = TRUE`. This convenience policy recodes the score column only;
it preserves person and facet identifiers because short labels such as `N` can
be legitimate IDs. Supplying an explicit character vector instead applies
that user-declared code set across the selected model columns, so inspect the
returned `missing_recoding` record before fitting or reporting exclusions.

The main workflow below uses a compact synthetic dataset with a connected
two-rater assignment, moderately unequal rater workloads, and six planned
criterion-level omissions represented by absent long-format rows, not `NA` or
sentinel scores. Its groups contain 24 persons each; three omissions per group
leave 141 observed rows in Group A and 141 in Group B:

```r
library(mfrmr)

dat <- load_mfrmr_data("example_operational")
data("mfrmr_example_operational_design", package = "mfrmr")
head(dat)
table(dat$Score)
list_mfrmr_data(details = TRUE)[, c("Key", "PrimaryUse", "Design", "CountBasis")]
```

All bundled examples are synthetic. `example_operational` is the applied
teaching dataset, not an empirical reference dataset. Use `example_core` only
when an idealized complete crossing is useful for a fast example, and use
`example_bias` for demonstrations with deliberately planted DFF/bias effects.

## Complete MML workflow

### 1. Check the data and intended score scale

Before fitting, state the complete rubric support and inspect retained,
missing, and zero-count categories:

```r
data_review <- describe_mfrm_data(
  data = dat,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  rating_min = 1,
  rating_max = 4,
  expected_design = mfrmr_example_operational_design
)
data_summary <- summary(data_review)
data_summary$structural_missingness
data_summary$design_connectivity
```

The assignment roster contains no scores. It tells `describe_mfrm_data()`
which Person x Rater x Criterion cells were planned. The review therefore
reports the six expected-but-unobserved cells separately from ordinary column
`NA` counts, while confirming that both observed Person-facet graphs remain
connected. Without `expected_design`, structural missingness is reported as
not assessed because an absent row may simply mean that the cell was never
assigned.

`rating_min` and `rating_max` retain unobserved boundary categories in the
data-support review. A boundary gap remains review evidence for the separate
element-boundary contract; it is not by itself an unsupported free-step
contrast. If an intended intermediate category is unobserved, use
`keep_original = TRUE` in both `describe_mfrm_data()` and `fit_mfrm()`. In a
polytomous fitted ladder, that retained internal gap creates an exact
adjacent-step recession direction, so fitting stops before optimization.
Otherwise non-consecutive observed scores such as `1, 2, 4, 5` are mapped to a
contiguous internal scale; always review the reported score map before
interpreting steps.

### 2. Fit the model

For a new analysis, start with marginal maximum likelihood:

```r
fit <- fit_mfrm(
  data = dat,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  rating_min = 1,
  rating_max = 4,
  method = "MML",
  model = "RSM"
)
```

`MML` integrates over the person distribution and returns posterior person
summaries. The default uses 31 quadrature points. Record that setting and
examine quadrature sensitivity when the application requires it. Eligible
fits below 15 points retain raw AIC/BIC/SABIC for screening, and fits at
15--30 points retain them for review, but automatic deltas, criterion weights,
preferred-model labels, evidence ratios, and LRT are disabled below 31 points.
Use 31--60 points as a comparison starting grid and 61 or more as a denser
sensitivity grid. A close or consequential comparison still requires a
prespecified common-grid sensitivity check; q>=31 alone is not evidence that
integration error is negligible.

After fitting a bounded GPCM-MML object as `fit_gpcm`, run the comparison
explicitly rather than making `summary()` refit the model in the background:

```r
q_review <- gpcm_mml_quadrature_sensitivity(
  fit_gpcm,
  data = dat,
  quad_points = c(31, 41)
)
summary(q_review)
apa_table(q_review)
```

The review reports changes in marginal likelihood per Person, relative slopes,
raw local-curvature SEs, population SD, and fitted probabilities. It does not
assign a universal stable/unstable cutoff, make raw slope SEs inferentially
eligible, or change the fit-readiness decision.

Use `model = "PCM", step_facet = "Criterion"` when category steps differ
across that facet. Choose the model from the scoring design and measurement
rationale, not from a single fit statistic.

### 3. Read the fit-only summary

The default summary is deliberately lightweight:

```r
fit_summary <- summary(
  fit,
  profile = "fit",
  detail = "brief"
)

fit_summary$decision
fit_summary$overview
fit_summary$status
fit_summary$readiness
fit_summary$data_review
fit_summary$settings_overview
fit_summary$facet_overview
fit_summary$person_overview
fit_summary$step_overview
```

Read `fit_summary$decision` first. It translates the stored readiness contract
into four practical questions: did the fit gates pass, has formal precision
support been evaluated, what evidence prevents formal use, and what should be
done next. A fit-only summary deliberately returns
`FormalInference = "No"` even when `InferenceReady = TRUE`, because convergence
and estimability do not by themselves validate standard errors, confidence
intervals, or reliability. Evaluate that separate contract with diagnostics:

```r
diag <- diagnose_mfrm(fit, residual_pca = "none")
summary(fit, diagnostics = diag)$decision
# Equivalent precision-aware decision:
summary(diag)$decision
```

The decision is a presentation of existing evidence, not a new statistical
test or an automatic model-selection rule. `FormalInference = "No"` still
permits explicitly labelled diagnostic inspection, but not formal SE/CI,
reliability, or significance claims. It does not mean that changing optimizer
settings until the answer becomes `"Yes"` is appropriate. Follow `NextAction`
and retain the original reason in reports.

Then review convergence and estimation settings. When optimizer code zero is
reached before the common terminal-gradient check passes, `fit_mfrm()` makes a
bounded sequence of warm-started polishing attempts when the requested setting
is at least as strict as the public default (`reltol <= 1e-9`). It retains the
best non-worsening stage under the recorded selection rule. The requested and
selected-stage settings, every attempted stage, terminal gradients, parameter
changes, and evaluation counts remain in `fit$opt$optimizer_polish`. For
L-BFGS-B, the native `factr` and `pgtol` controls are also recorded; do not
interpret `EffectiveReltol` as a native L-BFGS-B argument.

Treat `maxit` as a computational ceiling, not a convergence criterion or a
control to tune until preferred estimates appear. The package default is
`maxit = 400`; smaller values in executable examples exist only to shorten
checks. Prespecify the estimator and controls before inspecting results. If a
fit ends with `ConvergenceStatus = "iteration_limit"`, do not interpret or
compare its estimates. Refit the same data, model, method, anchors, optimizer,
tolerance, and quadrature rule using the next ceiling in a prespecified
sequence. Use a result only after `FitReadiness = ready`,
`InferenceReady = TRUE`, and `Numerical = pass`; do not select among runs by
coefficient size, fit statistics, significance, or agreement with an expected answer. Material
differences between separately ready runs indicate numerical instability that
requires review.

`InferenceReady` is deliberately a conservative fit-level first screen, not a
publication decision. It is `TRUE` only when Input, Estimability, Category,
Boundary, and Numerical components all pass. Then inspect the separate Design,
Stability, Diagnostics, and Reporting workflow rows in
`fit_summary$readiness`. A fit object saved before this versioned record existed
is labelled `legacy_unknown`; an older `InferenceReady = TRUE` value does not
make its summaries, results, or plots interpretation-ready. Refit it under the
current version to establish current readiness. Before optimization, mfrmr now
checks the estimator-specific constrained RSM/PCM free-coordinate design.
Exact rank deficiency stops with a structured `mfrmr_estimability_error`;
optimization cannot turn it into a usable fit. A disconnected design that is
full rank only under its declared constraints, an MML panel linked through a
common latent-population assumption rather than shared Persons, or a boundary-
constant facet level remains a reporting hold or review even when the optimizer
gradient is small. Inspect `fit$data_review$estimability` for the rank,
nullity, parameter blocks, and check scope. In an otherwise supported fit, a Reporting status such as
`ready_for_diagnostics_and_reporting_follow_up` means that fitting succeeded
and the next diagnostic stage is pending; it does not mean that optimization
failed or that the result is already manuscript-ready. Plots can still be
generated for diagnosis, but their
returned data carry `interpretation_status = "review_only"`, their subtitle is
marked `REVIEW ONLY`, the drawn title carries the same banner, and the plotting
call warns before substantive or cross-subset interpretation. The remaining
tables describe the fitted scale; they do not create universal acceptance
thresholds.

### 4. Request the comprehensive FACETS-organized summary

Use the `facets` profile for the main review:

```r
facets_summary <- summary(
  fit,
  profile = "facets",
  detail = "brief"
)

facets_summary
facets_summary$status
facets_summary$section_status
facets_summary$required_visual

res <- facets_summary$results
res$readiness
res$plot_map[, c(
  "Type", "Available", "InterpretationStatus", "InterpretationReady"
)]
```

This profile organizes model information, measures, uncertainty, fit evidence,
precision, category/step information, and plot routes in a reading order that
will be familiar to FACETS users. It computes the documented diagnostics when
they are needed and returns the resulting `mfrm_results` object in
`facets_summary$results`.

The profile name describes organization, not software execution:

- FACETS is not called;
- all estimates remain `mfrmr` estimates;
- sections that require a study-specific contrast, such as DIF or DFF, are not
  inferred automatically;
- residual PCA and multi-wave linking or drift analyses remain explicit
  follow-up decisions.

Pass a matching `diagnose_mfrm()` object through `diagnostics = diag` to
reuse work already completed. Use `compute = "never"` for a no-computation
review; dependent sections are then marked as not computed.

Person identifiers are omitted from the brief console view. Request detailed
person output only when the analysis purpose and data-handling plan require it.

### 5. Draw the required native Wright map

The native renderer is the primary targeting figure because it keeps
uncertainty visible:

```r
plot(
  res,
  type = "wright",
  renderer = "native",
  show_ci = TRUE,
  preset = "publication"
)
```

Read the vertical axis as the shared logit scale:

- the person distribution shows where the sample lies;
- facet levels show severity, difficulty, or the active signed orientation;
- step locations show the fitted category structure;
- intervals show uncertainty where an appropriate SE is available;
- gaps between the person distribution and the modeled locations may indicate
  weak targeting in that region.

Person and non-person uncertainty can have different statistical bases.
For example, an MML person interval is based on posterior uncertainty, whereas
a non-person interval may use observed-information uncertainty. Keep the
reported SE method with any table or figure interpretation.

The package defaults to a negative orientation for ordinary non-person facets:
higher measures indicate greater severity or difficulty. Facets named in
`positive_facets` have the reverse scoring direction. State the active
orientation in the figure caption.

Set `draw = FALSE` to obtain the fitted coordinates for a custom `ggplot2`,
Quarto, or accessibility-aware figure.

`top_n = Inf` retains and labels every fitted coordinate. The native text layer
keeps the fitted points fixed, moves only colliding labels, and connects them
with leader lines; `label_points` records both coordinates. Step thresholds
share a vertical ladder and are labelled with both the score transition and
the fitted logit. For genuinely dense maps, use the finite `top_n` compact view
or a larger output device rather than silently omitting interior labels.

When the fit review diagnoses boundary-separated facet levels and no explicit
`wright_range` is supplied, the native and FACETS-style renderers use the same
robust central range. Boundary levels appear at the appropriate ruler end with
triangles, while `OriginalEstimate`, `CI_Lower`, and `CI_Upper` retain the
untruncated values. `DisplayEstimate`, `DisplayCI_Lower`, `DisplayCI_Upper`, and
the `CIClipped*` / `CISuppressed` fields describe only what was drawn. This
prevents a huge separation interval from compressing the interpretable center
or being mistaken for a complete visible interval.

### 6. Add the optional FACETS-style Wright map

Use study-specific rubric labels rather than anonymous category numbers:

```r
rubric_labels <- c(
  "1" = "Beginning",
  "2" = "Developing",
  "3" = "Proficient",
  "4" = "Advanced"
)

plot(
  res,
  type = "wright",
  renderer = "facets",
  category_labels = rubric_labels,
  rows_per_logit = 2,
  show_ci = FALSE,
  preset = "publication"
)
```

The FACETS-style renderer provides:

- one common logit ruler;
- a person-frequency column drawn with `*`;
- signed facet headers and facet-level labels;
- horizontal step lines;
- adjacent-score labels such as
  `1 Beginning -> 2 Developing`;
- expected-score midpoint markers;
- an explicit `* = n person(s)` legend.

Set `persons_per_star` when the same star density must be used across several
figures. Otherwise, the renderer chooses a compact value and prints it below
the map.

The step display is designed to make the rating scale readable:

- a solid horizontal line is an estimated adjacent-category step location;
- its label names the lower and upper observed score categories and prints the
  fitted logit value in brackets;
- its height is the step location on the shared logit scale;
- a shorter dotted line is an expected-score midpoint crossing, not an
  additional model parameter.

For `RSM`, the step pattern is shared by the relevant observations. For
`PCM`, read the step column within its curve or step-facet group. Disordered or
closely spaced steps are prompts to inspect category use, sample support, and
the scoring design; they are not automatic instructions to collapse
categories.

Use the actual rubric wording from the instrument. The names of
`category_labels` must match the original scores. A two-column data frame with
columns `Score` and `Label` is also accepted.

With `draw = FALSE`, the returned `facets_style` component exposes
`category_labels`, `step_ruler`, `score_transitions`, and `settings` tables.
For line-printer reconstruction, `RulerValue` records the nearest discrete
ruler row; `DrawValue` records the exact coordinate used for step and midpoint
lines in the current renderer.

For a closer FACETS visual comparison, leave `show_ci = FALSE`. Setting
`show_ci = TRUE` adds `mfrmr` uncertainty whiskers to the asterisk ruler. That
hybrid is useful analytically, but it is intentionally an extension of the
FACETS-style layout; its footer identifies the interval level and the dot at
each whisker marks the corresponding fitted facet location.

Boundary-separated levels otherwise make an estimate-derived ruler very wide.
The automatic boundary-aware range described above is used by default; set an
explicit, reported range such as `wright_range = c(-4, 4)` when the application
requires a prespecified display scale. Out-of-range levels remain visible in
parentheses at the ruler ends. With `show_ci = TRUE`, endpoint triangles and
the footer identify intervals that extend beyond or are omitted from the
displayed ruler; exact values remain in the returned data.

### Visual correspondence is not numerical equivalence

The FACETS renderer reproduces the main reading grammar of a FACETS Table
6-style variable map. It does not reproduce a FACETS run.

| Question | What `mfrmr` provides |
| --- | --- |
| Familiar visual layout | Asterisk person counts, shared ruler, signed facet columns, step lines, and rubric-labelled transitions |
| Estimate source | The fitted `mfrmr` object |
| Pixel-identical output | Not guaranteed across graphics devices, fonts, or FACETS versions |
| Numerical equivalence | Not established by selecting `renderer = "facets"` |
| External comparison | Available after the user supplies output from a separately run FACETS analysis |

A defensible numerical comparison requires the same response records, score
coding, model, estimator, identification constraints, anchors, facet
orientation, extreme-score handling, and convergence criteria. FACETS commonly
uses JMLE, whereas the starter `mfrmr` workflow above uses MML with EAP person
summaries when its population-model assumptions are suitable. Those choices
target related but different calculations.

When a JMLE-oriented comparison is required, refit with `method = "JML"` and
still document every remaining setting. Matching the estimator family alone
does not establish equivalence.

For an exported FACETS fit table, use `read_facets_fit_table()` followed by
`facets_fit_review(fit, diagnostics = res$diagnostics, facets_fit = ...)`.
This reads supplied output; it does not launch or automate FACETS.

### 7. Review the Infit pathway, including persons

The fit pathway uses Infit on the horizontal axis and measure on the vertical
axis:

```r
plot(
  res,
  type = "fit_pathway",
  fit_stat = "Infit",
  fit_scale = "mnsq",
  include_person = TRUE,
  top_n_person = 12,
  person_labels = "none",
  facet_labels = "flagged",
  show_ci = TRUE,
  preset = "publication"
)
```

The vertical axis remains the fitted logit measure. The horizontal reference
at Infit MnSq = 1 represents model expectation; the outer lines are review
guides. Selected persons use a different point shape from non-person facet
levels. Completing this follow-up can change the Reporting readiness row while
leaving the already-passed Numerical row unchanged.

`top_n_person` limits the displayed person layer so a large study remains
readable. The selected persons are those with the largest fit distances;
non-person rows remain available. Use `person_subset` for a prespecified case
list, or `top_n_person = Inf` only when displaying every person serves a clear
purpose.

Fit statistics are evidence for review, not automatic exclusion rules.
Investigate response patterns, design cells, score support, and practical
consequences before changing data or operational decisions.

The separate `type = "pathway"` route displays expected scores and
dominant-category regions across theta; `type = "fit_pathway"` displays Infit
or Outfit against the fitted measure.

### 8. Build a report and export the results

Start with the brief result summary:

```r
results_summary <- summary(res, view = "brief")
results_summary$overview
results_summary$triage
results_summary$next_actions
results_summary$plot_map
```

Build a report-oriented object from the same fitted results:

```r
report <- mfrm_report(
  res,
  style = "qc"
)

summary(report, view = "reader")
report$first_screen
report$report_index
```

`mfrm_report()` organizes existing evidence. It does not refit the model or
turn diagnostic thresholds into a validity decision.

Export a reader-oriented, controlled analysis archive:

```r
exported <- export_mfrm_results(
  res,
  output_dir = "mfrmr-results",
  prefix = "analysis01",
  preset = "starter"
)

exported$summary
exported$written_files
```

The `starter` preset includes the result summary, report tables, replay code,
manifest, native Wright map, and focused plot routes. Existing files are not
overwritten unless `overwrite = TRUE` is requested explicitly.

This preset is a controlled **analysis archive**, not a deidentified or
automatically shareable deliverable. It includes a complete `.rds` result
object, and its tables, HTML, plots, replay code, and manifest can retain direct
person identifiers, person-level estimates, original labels, or local paths.
Review and transform every file under the study's data-handling policy before
sharing it. After making that assessment, `acknowledge_sensitive = TRUE` can
suppress the export warning; it does not remove or pseudonymize any data.

## A practical reading order

For most studies, review the output in this order:

1. Confirm the data roles, score range, row retention, and model.
2. Confirm convergence and the estimation settings.
3. Inspect the native Wright map with uncertainty.
4. Review targeting, facet measures, person summaries, and category steps.
5. Inspect Infit/Outfit and the person-inclusive Infit pathway.
6. Add residual, bias/DIF, linking, or interaction analyses only when the
   design and research question require them.
7. Review report caveats before exporting or drafting substantive claims.

Do not reduce this sequence to a single pass/fail index. Fit, precision,
targeting, category function, fairness, and validity answer different
questions.

## FACETS users

The closest translation of common FACETS concepts is:

| FACETS concept | `mfrmr` route |
| --- | --- |
| Data/specification roles | Explicit `person`, `facets`, and `score` arguments |
| JMLE-oriented fit | `fit_mfrm(method = "JML")` |
| MML analysis | `fit_mfrm(method = "MML")` |
| Measures and SEs | `summary(fit, profile = "facets")` and its result tables |
| Variable/Wright map | Native `renderer = "native"` or optional `renderer = "facets"` |
| Infit/Outfit review | `diagnose_mfrm()` and `type = "fit_pathway"` |
| Fair average | `fair_average_table()` |
| Bias/interaction screen | `estimate_bias()` and related review functions |
| Anchors | `anchors` and `group_anchors` in `fit_mfrm()` |
| External fit comparison | `read_facets_fit_table()` then `facets_fit_review()` |

For an existing FACETS-oriented script, `run_mfrm_facets()` provides a
one-call wrapper around package fitting and diagnostics. For new work, the
explicit `fit_mfrm()` workflow is easier to review and is recommended.

Use these guides for the full migration details:

- [Migrating from FACETS to mfrmr](vignettes/mfrmr-facets-migration.Rmd)
- [Visual diagnostics](vignettes/mfrmr-visual-diagnostics.Rmd)

The installed `references/FACETS_manual_mapping.md` maps concepts and output
routes. It is not evidence that FACETS was executed.

## ConQuest MML comparison

`mfrmr` can prepare and review a narrow external-table comparison for an MML
latent-regression case. The supported overlap is:

- `RSM` or `PCM`;
- binary responses;
- exactly one non-person facet, treated as the item facet;
- a unidimensional latent-regression MML fit;
- one numeric person covariate in addition to the intercept;
- complete rectangular person-by-item response data.

Within that scope, comparison targets include the population regression
slope, residual variance, centered item estimates, and case-level EAP
estimates.

The handoff is explicit: mfrmr prepares the analysis files, the user runs
ConQuest separately, and mfrmr then normalizes the requested exports.

```r
# fit_lr must satisfy the documented overlap conditions.
bundle <- build_conquest_overlap_bundle(
  fit = fit_lr,
  output_dir = "conquest-overlap"
)

# Run the generated .cqc file in ConQuest separately, then normalize its
# parameter, regression, covariance, and case-EAP CSV files.
conquest_tables <- normalize_conquest_overlap_exports(
  bundle,
  parameter_file = "conquest-overlap/conquest_overlap_conquest_parameters.csv",
  regression_file = "conquest-overlap/conquest_overlap_conquest_reg_coefficients.csv",
  covariance_file = "conquest-overlap/conquest_overlap_conquest_covariance.csv",
  case_file = "conquest-overlap/conquest_overlap_conquest_cases_eap.csv",
  conquest_version = "5.47.5",
  conquest_edition = "demo/free",
  run_date = Sys.Date()
)

conquest_review <- review_conquest_overlap(
  bundle,
  conquest_tables
)

summary(conquest_review)
conquest_review$attention_items
```

`mfrmr` does not execute or control ConQuest, and it does not parse arbitrary
raw ConQuest reports. The user runs ConQuest separately; mfrmr can normalize
the native comparison CSV exports requested by its generated command. The
review reports coordinate-level differences; it does not declare general
software equivalence from a single dataset or tolerance.

The generated ConQuest command uses the fitted mfrmr quadrature-point count;
the bundle records both values. Record the actual ConQuest version, edition,
and run date during normalization so the external comparison remains
reproducible.

This public bundle route does not cover multidimensional models, arbitrary
imported design matrices, bounded `GPCM` latent regression, JML latent
regression, or the full ConQuest plausible-values workflow. Separately, the
item-only bounded-GPCM parameterization can be compared only after the response
kernel, slope grouping, threshold coordinates, latent-scale identification,
retained rows, and category map have been matched. A result in that restricted
overlap does not extend automatically to a multifacet ConQuest generalized-item
design.

The ConQuest overlap bundle is also a controlled analysis bundle. Its long and
wide response files contain person identifiers and responses; the person-data
file contains identifiers and the covariate; and both mfrmr and ConQuest
case-EAP files contain identifiers and person-level estimates. The helper warns
when writing these files and creates `*_privacy_notice.csv`. Store the bundle in
an approved restricted location and pseudonymize or redact it as required
before sharing.

## Model and interpretation boundaries

| Area | Supported route | Important boundary |
| --- | --- | --- |
| Latent structure | One latent dimension | No multidimensional or Q-matrix engine |
| Facets | Multiple observed facet roles | The design must remain connected for the intended contrasts |
| `RSM` | Shared step structure | The common rating-scale assumption must be substantively defensible |
| `PCM` | Step structure associated with `step_facet` | Specify the step facet explicitly when the default is not intended |
| Bounded `GPCM` | Documented slope-aware core with `slope_facet == step_facet`; MML estimates the common scale by default | Not an unrestricted many-facet GPCM implementation |
| Estimation | `MML` and `JML`/`JMLE` | Estimator choice changes person summaries and residual-fit basis |
| Latent regression | Conditional-normal, unidimensional MML population model | Not arbitrary ConQuest design-matrix or multidimensional population modeling |
| Diagnostics | Residual/EAP and strict marginal screening routes | A flag is not a deletion, fairness, or validity decision |
| FACETS and ConQuest | Exported-table review within documented overlap | Neither external program is executed by `mfrmr` |

For bounded `GPCM`, inspect the capability table before choosing a downstream
helper:

```r
gpcm_capability_matrix()
vignette("mfrmr-gpcm-scope", package = "mfrmr")
```

Free-slope GPCM fits can be estimation-converged while parameter-level
inference remains review-only. Read `print(fit)`, `summary(fit)$decision`, and
the slope table's `ParameterStatus` and `PrimaryEstimate` before interpreting
the finite optimizer trace in `OptimizerEstimate`. `Optimizer*SE` and
`Optimizer*CI` are diagnostic quantities; ordinary slope SEs and confidence
intervals remain unavailable until the parameter-specific readiness checks
pass. The GPCM scope vignette explains each status and the appropriate next
action.

For MML, the default `gpcm_mml_identification = "free_population"` estimates
an intercept-only population distribution while relative slopes satisfy a
geometric-mean-one constraint. The corresponding fixed-latent-SD optimizer
coordinate is retained in `FixedLatentSDOptimizerEstimate`. Use
`gpcm_mml_identification = "fixed_standard_normal"` only when a deliberately
matched legacy or external comparison requires that identification.

The bounded GPCM is an **aligned single-owner relative-slope GPCM**:

$$
\log\frac{P(Y_o=k)}{P(Y_o=k-1)}
 = \alpha_{g(o)}\{\eta_o-\tau_{g(o),k}\},
\qquad \prod_g\alpha_g=1.
$$

Exactly one facet owns both slopes and steps, because
`slope_facet == step_facet`. Setting all slopes to one recovers the
equal-discrimination PCM kernel. This is narrower than the generalized MFRM of
[Uto and Ueno (2020)](https://doi.org/10.1007/s41237-020-00115-7), whose task
and rater slopes enter multiplicatively as `alpha_i * alpha_r` and whose step
owner is stated separately. The current
criterion-owned and rater-owned fits are therefore separate restricted
many-facet GPCM strata, not interchangeable fits of the full Uto--Ueno model.

The selected facet receives one slope per level, not one common slope for the
whole fit. For example, `slope_facet = "Criterion"` estimates a relative slope
for every criterion; `slope_facet = "Rater"` estimates one for every rater.
The other facets retain additive location effects but no slope block. Those
effects are still inside the complete adjacent-category predictor multiplied
by the selected slope; the current kernel is not a loading-only model with
unscaled facet intercepts. Criterion and rater slopes cannot be estimated
simultaneously in this model.

`plot(fit_gpcm, type = "ccc")` and `category_curves_report(fit_gpcm)` retain
the estimated slope for each step-facet level. These are reference-profile
curves, however: additive facet main effects and fitted interactions are fixed
at zero. The native plot facets multiple curve groups instead of overlaying
unlabelled traces, and both native and ggplot displays state the reference
profile condition. Inspect `CurveBasis`, `PredictorOffset`, and
`settings$curve_basis` before interpreting a curve as if it represented a
particular observed Person-by-facet cell.

For unpenalized JML, all-minimum or all-maximum Person patterns can have an
unbounded primary ability estimate. A finite adjusted display, when supplied,
is kept separate from the likelihood-based primary status. The same principle
applies to any GPCM slope whose boundary status has not been resolved: a finite
optimizer iterate is not automatically a finite maximum suitable for ordinary
inference.

PCM and bounded GPCM can be reviewed on the same data with
`compare_mfrm(fit_pcm, fit_gpcm)`, `build_weighting_review(fit_pcm, fit_gpcm)`,
or `build_model_choice_review(..., run_weighting_review = TRUE)`. Selectable
information-criterion ranking requires comparable, inference-ready MML fits
on a common quadrature grid with at least 31 points. Although PCM is the
all-unit-slope reduction of the aligned GPCM kernel, the current automatic
nesting contract withholds the PCM-versus-GPCM chi-square LRT and records
`PCM_in_GPCM_ic_only` instead.

The practical comparison is available without reconstructing the model
matrices manually. After fitting the same data as `fit_pcm` and `fit_gpcm`
(see the GPCM scope vignette), use:

```r
choice <- build_model_choice_review(PCM = fit_pcm, GPCM = fit_gpcm)
choice$model_roles[, c(
  "Model", "StepCoordinates", "FreeStepParameters",
  "SlopeCoordinates", "FreeSlopeParameters",
  "FitReadiness", "FormalInference", "Interpretation"
)]
```

The coordinate columns count reported values; the free-parameter columns also
apply the fitted constraints. `build_weighting_review()` keeps JML and MML
evidence separate. Comparable, inference-ready MML fits can contribute
information criteria, whereas an unpenalized JML likelihood difference is not
turned into an automatic PCM-versus-GPCM choice.

Cross-software slope values are not automatically matched estimands. FACETS
does not jointly fit Muraki's free-slope polytomous GPCM. [FACETS' reported
element discrimination](https://www.winsteps.com/facetman64/t7menu.htm) is a
post-fit diagnostic computed after the Rasch measures and does not feed back
into the other estimates; it must not be treated as a free-GPCM slope estimate
from `mfrmr`. TAM can estimate GPCM
slopes through its 2PL/GPCM MML route, but its many-facet fitting route does
not estimate those slopes. The current `immer` estimation routes provide PCM-
design and hierarchical-rater references rather than a matched free-GPCM fit.
Consequently, FACETS and `immer` are appropriate only for documented
equal-discrimination overlaps or deliberately different-model sensitivity
checks. A TAM GPCM comparison is numerical only after the response kernel,
slope grouping, threshold parameterization, latent-scale identification,
retained rows, category map, and covariance information have been matched.

For strict MML diagnostics, keep the two evidence bases distinct:

```r
diag_both <- diagnose_mfrm(
  fit,
  residual_pca = "none",
  diagnostic_mode = "both"
)

summary(diag_both)$diagnostic_basis
```

The legacy route evaluates residual-oriented summaries using person estimates.
The strict marginal route uses latent-integrated expected values. They are
complementary screens and should not be collapsed into one decision rule.

## Documentation

The package includes the following vignettes:

- [End-to-end workflow](vignettes/mfrmr-workflow.Rmd)
- [MML estimation and marginal-fit diagnostics](vignettes/mfrmr-mml-and-marginal-fit.Rmd)
- [Migrating from FACETS](vignettes/mfrmr-facets-migration.Rmd)
- [Visual diagnostics](vignettes/mfrmr-visual-diagnostics.Rmd)
- [Reporting and APA-oriented output](vignettes/mfrmr-reporting-and-apa.Rmd)
- [Linking and DFF](vignettes/mfrmr-linking-and-dff.Rmd)
- [Bounded GPCM scope](vignettes/mfrmr-gpcm-scope.Rmd)

Open a guide with `vignette("mfrmr-workflow", package = "mfrmr")` or list
them with `browseVignettes("mfrmr")`. A printable API overview is installed
at `cheatsheet/mfrmr-cheatsheet.pdf`. Function-level help starts with
`?fit_mfrm`, `?summary.mfrm_fit`, and `?plot.mfrm_fit`.

## Citation

Use the package citation supplied with the installed version:

```r
citation("mfrmr")
```

The underlying ordered-response models are described by Andrich (1978),
Masters (1982), and Muraki (1992). See the package help and vignettes for the
references relevant to each model and diagnostic route.

## License

`mfrmr` is licensed under the
[MIT License](https://opensource.org/licenses/MIT).

## Acknowledgements

`mfrmr` has benefited from discussion and methodological input from
[Dr. Atsushi Mizumoto](https://mizumot.com/) and
[Dr. Taichi Yamashita](https://kugakujo.kansai-u.ac.jp/html/100000882_en.html).

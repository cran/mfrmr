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

The package extends Rasch-family RSM/PCM work with MML, modern diagnostics,
reproducibility, network review, and reporting support. It is not a general
FACETS replacement: each `fit_mfrm()` call uses one response-model family and
one observed score scale, and the current public API does not provide mixed
response families, multiple independent rating scales, general threshold
anchoring, or fixed-calibration operational scoring.

The recommended workflow is:

```text
long-format data + describe_mfrm_data()
  -> fit_mfrm(method = "MML")
  -> summary(profile = "facets")
  -> native Wright map with uncertainty
  -> focused diagnostics
  -> report and export
```

The native Wright map is the required first fitted-scale figure in this
workflow. A separate FACETS-style renderer is available when a familiar
asterisk ruler and labelled category transitions are useful.

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

Load the package and list its installed guides:

```r
library(mfrmr)
browseVignettes("mfrmr")
```

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

`rating_min` and `rating_max` retain unobserved boundary categories. If an
intended intermediate category is also unobserved, use `keep_original = TRUE`
in both `describe_mfrm_data()` and `fit_mfrm()`. Otherwise non-consecutive
observed scores such as `1, 2, 4, 5` are mapped to a contiguous internal scale;
always review the reported score map before interpreting steps.

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
examine quadrature sensitivity when the application requires it.

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

fit_summary$overview
fit_summary$status
fit_summary$readiness
fit_summary$data_review
fit_summary$settings_overview
fit_summary$facet_overview
fit_summary$person_overview
fit_summary$step_overview
```

Start with convergence and estimation settings. When optimizer code zero is
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
sequence. Use a result only after `Converged = TRUE`, `InferenceReady = TRUE`,
and `Numerical = pass`; do not select among runs by coefficient size, fit
statistics, significance, or agreement with an expected answer. Material
differences between separately ready runs indicate numerical instability that
requires review.

`InferenceReady` is deliberately a numerical status, not a publication
decision. Require `Numerical = pass`, then inspect the separate Data, Design,
and Stability rows in `fit_summary$readiness`. A disconnected design or a
boundary-constant facet level remains a reporting hold even when the optimizer
gradient is small. In an otherwise supported fit, a Reporting status such as
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

`top_n = Inf` retains every fitted coordinate. The native text layer remains
collision-aware, so a dense map may leave some retained points unlabeled;
the returned plot data and retention table remain the complete record.

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

The handoff is:

```r
# fit_lr must satisfy the documented overlap conditions.
bundle <- build_conquest_overlap_bundle(
  fit = fit_lr,
  output_dir = "conquest-overlap"
)

# Run the generated .cqc file in ConQuest separately, then normalize
# the four CSV files requested by that command.
conquest_tables <- normalize_conquest_overlap_exports(
  bundle,
  parameter_file = "conquest-overlap/conquest_overlap_conquest_parameters.csv",
  regression_file = "conquest-overlap/conquest_overlap_conquest_reg_coefficients.csv",
  covariance_file = "conquest-overlap/conquest_overlap_conquest_covariance.csv",
  case_file = "conquest-overlap/conquest_overlap_conquest_cases_eap.csv",
  conquest_version = "5.47.5",
  conquest_edition = "demo/free",
  run_date = "2026-07-23"
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
the four native CSV exports requested by its generated command. The comparison
reports differences and does not declare software equivalence from a fixed
tolerance.

The generated ConQuest command uses the fitted mfrmr quadrature-point count;
the bundle records both values. Record the actual ConQuest version, edition,
and run date during normalization so the external comparison remains
reproducible.

A public
[aggregate comparison record](https://github.com/Ryuya-dot-com/mfrmr/blob/main/inst/validation/conquest-mml-overlap-0.2.2.md)
of a matched 31-node check with ConQuest 5.47.5 is available in the source
repository. The record is excluded from
the installed CRAN package and supports only the overlap case stated above;
identifier-bearing response and case-level files are not included in the
package.

This route does not cover multidimensional models, arbitrary imported design
matrices, bounded `GPCM` latent regression, JML latent regression, or the full
ConQuest plausible-values workflow.

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
| Bounded `GPCM` | Documented slope-aware core with `slope_facet == step_facet` | Not an unrestricted many-facet GPCM implementation |
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

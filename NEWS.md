# mfrmr 0.3.2

- Loading mfrmr defers data-manipulation dependencies until they are needed.
  Examples for extraction and display reuse supplied calibrations; estimation
  examples continue to fit the supplied rating data.


- Expand method abbreviations in the package description and correct the
  reference identifier for Raju (1990).

- Clarify the drift model's normal quadrature grid and the available D-Study
  result accessors. `lpSolve` is now required for MML structural-boundary
  audits, keeping inference-eligibility checks available in standard installs.
- Posterior-bundle import and Stan/brms model-generation tools are not included
  in 0.3.2. Read saved bundles with their producer version; BayesianMGMFRM.jl
  interoperability is planned for a later release.
- Print common-reference Fair Average scores and intervals, pairwise bias
  contrasts and joint-model LRTs. Word fallback tables use the document default
  style. Reports reuse reliability results only for the matching calibration.
- Cut-score inversion follows the same finite facet means and reference
  discrimination as Fair(M); excluded persons no longer hide a common cut logit.
- Dataset loading covers all 13 bundled tables and provides
  `load_mfrmr_data("list", details = TRUE)` for their different column roles.
- Package citations include the CRAN DOI.


- Estimation, diagnostics, and G-Study results use `mfrm_fit`,
  `mfrm_diagnostics`, and `mfrm_gtheory`. Refit saved calibration objects from
  earlier versions using their original data and model specifications; changing
  a class label does not convert a saved object.
- Profile-adjusted JML controls and comparison-only fit objects are no longer
  supported. `JML_jeffreys` is a distinct estimator and does not reproduce
  profile-adjusted JML results. Compare fitted models with `anova()` or
  `compare_mfrm()`.
- Bias and DFF analyses require a current fitted calibration. Select the test
  explicitly with `method`; former `residual` and `refit` switches no longer
  silently choose another procedure.
- MML automatic quadrature compares two orders and records numerical sensitivity;
  unstable checks withhold inferential reporting.
- ANOVA, MINQUE(0), and Henderson G-Study estimates use consistent result
  fields while retaining estimator-specific assumptions. Balanced ANOVA keeps
  MLS intervals; Henderson intervals with Satterthwaite degrees of freedom
  below one are withheld.
- Fair Average `confint()` returns conditional ability-based score intervals.
  Legacy calibration-bootstrap and joint-delta options are not supported.

- Person-fit tables use consistent columns across supported estimators,
  including `$table`, conditional JML ability correction, and matching
  summary/plot screening flags and unavailable-state reasons.
- Explicit composite weights retain their supplied scale and named order;
  omitted weights still define an equal-weight average.
- Fits with non-unit observation weights or facet locations linked only through
  population assumptions require review before inference.
- Rater cards report central-tendency correlations by criterion and explain
  why raw residual correlations alone do not establish halo.
- Markdown and HTML report tables preserve literal pipe/backslash characters
  and HTML-special characters in identifiers and diagnostic text.

- Tidy extractors use the shared `generics` interface, including when broom
  or tidymodels is attached.
- G-Study composite weights require explicit `method = "urgenova"`; they
  no longer override a requested ANOVA or MINQUE(0) estimator.
- Henderson component intervals with Satterthwaite degrees of freedom below
  one are withheld with a reason in `confint()`, tidy tables, and reports.
- Workflows accept all supported fitted calibrations. Diagnostics are reported
  only when they match the calibration attached to the workflow.

- Calibration and new-data prediction reject duplicate identifiers in the
  person or facet covariate tables they use.
- Workflow score bounds and extreme-score flags use the fitted calibration's
  category support, including when an earlier saved workflow is reopened.
- A G-Study above four non-person facets reports its limit and leaves the
  calibration and diagnostics available in the one-call workflow.

- The unified calibration and design check reject factor-valued score
  columns before R can substitute their internal codes. Saved unified fits
  with such inputs are marked for review.
- Unified RSM and PCM MML detect rank deficiency in the joint design
  for population means and free facet locations before estimation. Saved
  fits are reviewed for the same issue when used for inference.

- Re-estimating a workflow clears diagnostics and estimator comparisons from
  the previous calibration. Saved diagnostics that disagree with the attached
  calibration are withheld until diagnostics are run again.
- New-data WLE checks person IDs and category support before scoring. EAP
  prediction with supplied ratings now explains that this scoring path is
  unavailable. Omit `newdata` to retrieve fitted EAP scores.

- Fitted models with discrimination estimates display as Generalized
  Many-Facet Rasch Model (GMFRM), with the facets carrying discrimination.
  Following Uto and Ueno (2020), this name denotes an extension of the
  constant-discrimination many-facet Rasch model. Use `model = "GMFRM"` or its retained code `"GPCM"`
  in code. The display label is available as `summary(fit)$model_label`.

- `mfrm_analysis()` connects calibration, diagnostics, G-Study and D-Study,
  plots, and Markdown, HTML, and optional Word reporting in `mfrmr`.
  The workflow steps, `tidy()`, `glance()`, `augment()`, and
  `render_mfrm_report()` are available from the same package. Saved workflows
  retain their class, and `tidymfrm()` is a compatibility alias.

- Assessment reports share estimator notes, facet separation, agreement with
  model expectations, extreme-score proportions, anchor displacement, and
  multivariate G-Study covariance and effective weights across formats.
  A common reference plan can include conditional original-scale intervals.
  Select a rater for a focused feedback report. Connectivity audits can also
  report algebraic connectivity and panel effective resistance.
- Baseline drift screening uses declared stable reference anchors for linked
  comparisons. Unqualified comparisons with separately centered estimates are
  reported as unavailable. Anchor displacement is reported separately as a
  conditional one-step diagnostic.

- Rating plans now apply their counts to information curves for fitted calibrations.
  Observed exposure can also be reported per person. The original
  one-per-context reference curve remains available by default.
- The design check reports Person-by-rater components and individual linking
  persons separately from the multi-facet graph.
- Pathway plots can show Infit and Outfit side by side. Residual PCA is
  available on request for fitted calibrations, with a first-contrast plot.

- Fits from a shared-design score matrix retain the ratings from their own
  replication for residuals, diagnostics, and reuse. The batch fitter checks
  every column against the common design and rubric.
- An unidentified urGENOVA variance decomposition or too few complete
  observations signals `mfrm_gstudy_unavailable` for callers that retain
  results from other analyses.

- Local GPCM bias calculations retain fitted discrimination products and slope
  formulas. Joint bias refits preserve the source model specification and use
  the common nested-likelihood checks. Unsupported split-DFF specifications are
  rejected; supported refits retain weights, population covariates and controls.
- Bias, DFF and diagnostic outputs distinguish unavailable inference from a
  negative screen, including summaries, plots and reasons for missing values.
  Split DFF checks both refits; fixed-calibration methods check the source.
- G-Study and D-Study return NA with a reason for zero total variance or invalid
  projected variance, preserving a defined zero coefficient where appropriate.
  Saved G-Study displays recalculate these states from stored variances.

- GPCM Wald slope intervals now use log-slope uncertainty and exponentiation,
  matching their documentation and retaining positive endpoints. Model and
  person-clustered sandwich intervals both use this construction.
  Slope plots use these intervals too. Unified person, facet and slope plots
  withhold interval bars when calibration is not ready, retaining descriptive
  points and the reasons in plot notes.
- The measurement tutorial uses matching 31-point integration for nested
  RSM/PCM/GPCM comparisons. The G-Study tutorial explains equal-budget planning,
  relative versus absolute decisions, and the limitations of interval and
  design projections.

- Correct uncertainty and validation in the unified engine:
  - Reject unidentified G-Study designs; use the actual assignment pattern in
    Gaussian plug-in component covariances, including unbalanced observations.
    Honor `missing = "error"` and require positive integer D-Study counts.
  - Require matching observations, weights, likelihood settings and eligible
    fits for nested LRTs. Record optimizer convergence separately from inference
    readiness; detected structural MML boundary directions prevent intervals
    and ordinary classification.
  - Use posterior variance decomposition for EAP reliability and recompute
    diagnostics and fair scores from the selected primary EAP/WLE score.
  - Include full structural covariance, steps, slopes and interactions in
    conditional prediction SEs; person-score uncertainty remains excluded.
  - Apply `missing_codes` before unified estimation, honor weighted common
    `reference_plan` scoring, and reject duplicate unit-rater ratings in
    agreement calculations. Common-plan unified scores are plug-in predictions;
    request a conditional interval with `fair_se = TRUE`. The interval holds
    the calibration and reference plan fixed. The point prediction is not a
    posterior mean of the transformed score.
  - Document category support as `Level`, `Min`, `Max`, and add explicit S3
    method signatures, uncertainty definitions and eligibility conditions.

- Consolidate the public API around the analysis workflow and standard R S3
  methods (`predict()`, `coef()`, `vcov()`, `confint()`, `anova()`, `plot()`,
  `fitted()`, `residuals()`, `logLik()`, `AIC()`, `BIC()`, `nobs()`, `summary()`,
  and `print()`), removing external-software wrappers (`TAM`, `ConQuest`, and
  `FACETS` importers/crosswalk guides).
- Upgrade `fit_mfrm()` (unified C++ estimation),
  `check_mfrm_connectivity()`, and `mfrm_batch_fit()`:
  - Native C++ solvers for Joint Maximum Likelihood (`"JML"`), Jeffreys/Firth
    penalized JML (`"JML_jeffreys"` with exact Block-Schur inverse-information
    SEs), Marginal Maximum Likelihood (`"MML"` with observed and person-clustered
    sandwich covariances, `EAP`, and Warm's 1989 `WLE`), Conditional Maximum
    Likelihood (`"CML"` for `"RSM"` and `"PCM"` via elementary symmetric
    polynomials with full conditional covariance and `WLE`), and GMFRM
    (`model = "GPCM"`, `"MML"` for one or two discrimination facets and
    `"JML"` for one discrimination facet, with `"full_predictor"` and `"ability_only"`
    slope parameterizations, `slope_prior_sd` log-normal regularization, Warm
    GPCM `WLE`, and 4-level nested affine LRTs).
  - Full C++ support for `group_anchors`, `facet_shrinkage` (`"empirical_bayes"`
    and `"laplace"`), multi-facet collinearity rank checks, and MML latent-regression
    boundary separation checks (`lpSolve`).
  - Automatic pre-estimation connectivity and interior-category support warnings
    (`warn_connectivity = TRUE`).
  - Extended `step_constraints` (fixed and shared thresholds across `MML` with
    fixed or estimated normal populations, `CML`, `JML`, and `JML_jeffreys`,
    with Wald, sandwich, and 1D profile-likelihood `confint(..., se_type = "profile")`
    intervals and affine-subspace `anova()` / `compare_mfrm(nested = TRUE)` LRTs).
  - Declarative `category_support` for mixed binary/polytomous criteria with
    exact zero probability outside supported categories.
  - Explanatory `RSM`/`PCM` (`facet_formulas` + `facet_data`) for linear
    covariate decomposition of facet locations, including `predict()` on unseen
    facet levels and nested LRTs.
- Accelerate `estimate_bias()` and `analyze_dff()` in C++:
  - `estimate_bias()` supports single-pair and `interaction_facets = "all"` C++
    bias tables, defaulting to the FACETS 4.5.2/4.5.3 `RajuSA-Welch`
    (von Davier, 2025; Raju, 1990) integer-truncated degrees of freedom
    (`df_method = "facets"`, with `"classical"` Welch-Satterthwaite also
    available), plus joint C++ interaction refits (`method = "joint"`).
  - `analyze_dff()` adds warm-started parallel C++ split-element `"mml_lrt"`
    (with latent-regression impact control and sandwich SEs), `"cml_lrt"`,
    `"raju_welch"`, `"mh"`, and `"logistic"`.
- Unify `mfrm_g_study()` (adding non-iterative Henderson Method 1
  `method = "urgenova"` for crossed and nested multi-facet designs, multivariate
  composite `weights`, and Satterthwaite/delta-method `confint()`) and
  `interrater_agreement()` (exact/adjacent agreement, weighted kappa, Gwet's
  `AC2`, Krippendorff's `alpha`, and `ICC`).
- Unify all `plot()` S3 methods on `ggplot2` (including `"wright"`, `"pathway"`,
  `"steps"`, `"category_curves"`, `"persons"`, `"raters"`, `"person_fit"`,
  `"residual_qq"`, `"slopes"`, and `"qc"`) with color-vision-safe `"okabe_ito"`,
  `"viridis"`, and strict `R == G == B` `"monochrome"` palettes paired with
  redundant point `shape` and `linetype` encodings.

# mfrmr 0.3.1

- Add `fit_mfrm_drift()` (`"mfrm_drift"`) for non-Bayesian Markov
  rater-severity drift modeling across ordered scoring waves, sessions, or days
  without MCMC sampling (adapted and generalized from Uto, 2023). Examinee
  abilities are integrated by normal quadrature, inner parameters use the exact
  analytical Louis (1982) observed information matrix with optional sum-to-zero
  fixed facets (`facets`) and shared (`"RSM"`) or rater-specific (`"PCM"`) steps,
  and rater-specific log-drift SDs are estimated by Laplace-marginalized
  Empirical Bayes with Kass-Steffey (1989) hyperparameter covariance propagation
  (`TotalSE` vs `ConditionalSE`), cross-wave bridge-person checks, reference-wave
  `FairM` prediction (`predict()`), and CVD-safe trajectory, contrast, and
  drift-SD plots (`plot()`, `plot_data()`, `as_ggplot()`).

- Supplying either rating bound now preserves the rubric by default. Use
  `category_policy = "collapse"` to reproduce earlier fits with gaps. Facet
  interval tables from MML fits carry the source fit and identification status.

- ZSTD comparisons retain unavailable decisions. Design simulations count
  exact values of -2 and 2 in the ZSTD screening rate. `fit_measures_table()`
  shows both strict and inclusive ZSTD flags at the requested cutoff and
  preserves the df convention of saved diagnostics. Rebuild older fit tables
  from their diagnostics to update the companion columns.
- Dashboard MnSq screening includes both endpoints in the review band. Rebuild
  saved dashboards from the fit and diagnostics to update their flags.
  Diagnostic summaries and misfit casebooks accept `lower` and `upper`.

- Reject additive RSM/PCM MML fits when the observed ratings and latent
  regression permit population or facet locations to diverge. Designs requiring
  a numerical boundary check use the suggested package `lpSolve`. Recalculating
  adjusted scores from an affected saved fit also reports the boundary.

- Direct RSM and PCM MML fits support fixed and partially shared centered
  step offsets through `spec$step_constraints`, including estimated normal
  `population_formula` designs, affine step-subspace nested likelihood-ratio
  tests in `compare_mfrm(..., nested = TRUE)`, and reuse in a new R session.
- `diagnose_mfrm()` computes MML `Person` separation and reliability from the
  Bayesian EAP posterior-variance decomposition
  `Var(EAP) / (Var(EAP) + mean(SE^2))` rather than unshrunk variance
  subtraction. Recreate saved MML diagnostics to update `Person` reliability.
- `fair_average_table()` supports `reference_plan` on native additive RSM/PCM
  MML fits (returning strictly monotone EAP plug-in scores plus conditional
  quadrature `PosteriorExpectedScore` and `PosteriorScoreSD`) and records
  `pcm_exposure_audit` when incomplete PCM step-scale exposure causes
  measure-vs-fair rank inversions.
- `compute_information()` adds `exposure_basis = "per_person_average"` to scale
  observed context weights and `ConditionalSEM` to one active Person without
  constructing a manual rating plan.
- `plot(x, type = "fit_pathway", fit_scale = "zstd")` keeps dimensionless
  FACETS fourth-moment df primary in its plot payload when
  `zstd_method = "facets"`.
- Upgrade `ggplot2` diagnostics and color-vision accessibility across presets:
  - `standard`, `publication`, and `compact` use Okabe-Ito and Paul Tol
    color-vision-safe palettes paired with redundant point shapes, line types,
    and direct status markers; `monochrome` uses pure neutral grayscale with a
    monotonic luminance ramp for diverging correlation and residual heatmaps.
  - `plot(fit, type = "wright", renderer = "facets")` renders the FACETS Table
    6 star, signed-facet, and scale-threshold layout directly as a `ggplot`
    object, and native Wright maps show facet orientation polarity (`(+)` vs
    `(-)`).
  - Dedicated `ggplot2` converters replace generic tabular fallbacks for
    `"facet"`, `"person"`, `"step"`, `"person_fit"` (both `"meansquare"`
    Infit-vs-Outfit bubble scatters with marginal histograms and `"loglik"`
    ranked `lz*`/`lz` lollipops), `"marginal_fit"`, `"marginal_pairwise"`,
    `"unexpected"`, `"facets_chisq"`, `"guttman"`, `"residual_matrix"`,
    `"trajectory"`, `"reliability"`, `"shrinkage"` (`view = "funnel"`), and
    named-facet design grids.
  - `plot(fit, type = "qc")` renders all nine documented diagnostic panels;
    `plot(fit, type = "apa_summary")` adds a structured model-and-readiness
    summary card; `plot(fit, type = "pathway")` places step markers on the
    expected-score curve; `plot(fit, type = "thresholds")` offsets disordered
    step crossings with dashed brackets and triangle markers;
    `plot(fit, type = "residual_qq")` adds a 95% normal order-statistic
    pointwise reference band; category-curve and `mfrm_d_study()` projections
    pair colors with distinct line types and integer facet-count ticks;
    `preset` is now accepted across `facet_equivalence`, `facet_dashboard`,
    `run_qc_pipeline()`, `compute_information()`, and report bundles; and
    composite dashboards suppress outer coordinate axes when restyled with
    `ggplot2` themes.
- D-study planning tables now show the facet counts for each scenario,
  including comparisons saved with 0.3.0.
- Peer-review design summaries count reciprocal pairs between distinct students.
  For earlier summaries, rerun `review_peer_design()` on the saved simulation.
- Casebook person views include person fit flags. Warnings and source counts
  include element fit. Unexpected-response tables retain their columns when
  no evaluated rows are flagged.
- Data preparation shows why ratings were dropped. Fixed-facet reports show
  the recorded reasons when intervals are unavailable.
- Fit summaries point to identification and interval eligibility. `mfrm_report(...,
  style = "rater")` now supports native additive one-family GPCM fits alongside
  RSM/PCM and two-family GPCM, withholding working values when primary slope
  estimates are unavailable.
- Replace internal `stringr` helpers with base R equivalents and remove `stringr`
  from `Imports`, avoiding runtime failures when optional `stringi` binaries are
  unavailable.

# mfrmr 0.3.0

- Consolidate related S3 method and `mfrm_report()` component help pages into
  shared reference topics while preserving all help aliases, formulas, scale
  and identification notes, unavailable-result states and runnable examples.
  Extended historical and benchmark discussions move to website articles.

- Reuse matching PCA-free diagnostics in optional diagnostic-screening reports.
  Report audit fields and optimization replay remain unchanged; report work
  remains outside the replicate's `ElapsedSec` value.

- Return a typed, zero-row component index for empty saved plot data instead
  of failing in `plot_data_components()`. Existing plot data remain unchanged.

- Use `score_comparison_table(..., screen = "williamson2012")` instead of the
  former `automated_scoring_table()` interface. The generic default still
  applies no thresholds. Score/reference column names replace the earlier
  automated/human names; see `mfrmr_api_migration_030` for argument and output
  mappings. Earlier saved result lists remain readable.

- Consolidate facet statistics, precision, bias iteration/pairwise, DFF,
  checklist and reference-review reports in `mfrm_report(x, type = ...)`.
  Their dedicated entry functions are removed; specialized help and old help
  aliases remain. Result classes, calculation defaults and saved-result methods
  are retained. Use complete named component options; unsupported or duplicate
  options are rejected, including options formerly ignored by `dif_report()`.
  See `mfrmr_api_migration_030` for the seven replacement calls.

- Use `mfrm_report(..., type = "iterations")` for replayed estimation diagnostics
  and `type = "bias_interaction"` for ranked interaction-bias reports. The former
  dedicated entries are removed. Replay still runs optimization; saved bias
  estimates can be reported without recalibration. Result classes and methods
  remain available.

- Compare D-study plans through `compare_mfrm(d, assumption = "normal", ...)`.
  The former `mfrm_multivariate_d_compare()` entry is removed; result classes,
  paired covariance, inference scope and saved-result methods are retained.

- Prevent integer overflow when comparing large future rating counts.
  Integer- and double-stored counts now yield the same comparison calculations.

- Consolidate specification, data-quality, category-structure and category-curve
  reports in `mfrm_report(fit, type = ...)`. Their dedicated entry functions are
  removed; see `mfrmr_api_migration_030` for replacements. Named component options,
  defaults, result classes, summaries and plots retain their meanings.
  Collected-result reports continue to use `mfrm_report(res)`.

- Make `mfrmr_output_guide()` show the five central analysis routes by default.
  Use `scope = "all"` for the full specialist catalogue. Explicit scopes,
  including the six-row `"public"` and `"beginner"` routes, keep their meanings.

- Share standard summary-table assembly and printing while preserving table order, conversion rules, reporting roles, descriptions and console output.

- Retain both series and their colours in appendix-role plot legends and
  correct their ggplot conversion, including older saved plot data. Keep
  appendix preset names separate from visual style settings. Older plot data
  with only one recorded colour use the default for the missing series;
  recreate the plot from its original bundle and palette for both custom colours.

- Keep posterior contrasts with different parameter quantities or scales out
  of a shared plot axis. Mixed summaries and reports remain available;
  select a common quantity and scale for each figure. Label slope contrasts
  and log-slope contrasts separately from location contrasts. Recreate plots
  from saved posterior results; no new sampling is required.

- Correct unavailable-plot status and missing uncertainty reasons when
  corrected JML finds no supported point solution. Valid points remain
  available when only RootSE is unavailable. Recreate results and reports
  from a saved fit to update the status; refitting is unnecessary.

- Speed up two-family GPCM EM fitting by omitting unused reporting tables
  during likelihood evaluations. Numerical controls and returned results
  are unchanged.

- Update imputed-feature validation tests for the logical-class restoration
  in `mice` 3.20.0. Clustering behavior and input validation are unchanged.

- Connect rating-plan comparisons with rater feedback and saved results in
  the visual tour. Clarify uncertainty and draft-completeness messages while
  preserving estimates, availability decisions and saved field names.
  Decision text is consistent in tables, print and HTML, including redisplayed
  older saved decisions. Rating-plan links open the comparison directly.
  The G-theory saving example retains resource reviews, report context and
  edited ggplot figures alongside the corresponding studies.

- Allow `compute_information()` to evaluate an explicit per-Person rating
  plan for supported additive MML fits, using fitted facet levels and rating
  counts. Saved plots retain the plan, its fixed-calibration interpretation
  and any use of exploratory GPCM coordinates. The default observed-exposure
  calculation is unchanged.

- Add separate limiting-likelihood point estimates for supported one- and
  two-family GPCM MML fits with one globally unused endpoint category.
  Original data and full-model boundary status are retained. Available slopes,
  original category labels, reduced-coordinate tables and numerical failures
  are carried through summaries and reports. Limiting curves include the
  unused endpoint at probability zero; ordinary intervals remain unavailable.
  Displaying an old saved fit does not trigger a refit.

- Correct MML boundary reporting for unpenalized RSM/PCM and GPCM fits with provable
  improving likelihood paths from extreme free facet levels or unused endpoint categories,
  including centered facets with a free compensating location or population mean.
  Finite Person EAPs no longer label these full models as having a finite maximum.
  Saved-fit summaries explain retained optimizer coordinates, and ordinary
  structural intervals are withheld on boundary evidence. Likelihood calculations
  and category recoding defaults are unchanged.
  Two-family GPCM no longer reports endpoint-boundary optimizer values as
  primary point estimates. Rebuild affected saved result collections from
  their saved fits; model refitting is unnecessary for this reporting correction.
  Boundary optimizer values are not used for fitted curves or rater estimates.

- Label empty exact-F confidence sets explicitly in D-study interval plots,
  distinguishing them from numerical unavailability and an estimate of zero.
  Estimates and interval calculations are unchanged.

- Clarify the inputs, returned objects and first tables to read for model
  summaries, diagnostics, result collection, reports and bias screening.
  Consolidate repeated explanations while retaining uncertainty and
  availability conditions. Explain estimate direction, posterior SD versus
  structural SE, fit-adjusted reliability and GPCM working versus primary
  values. Clarify one-step displacement, diagnostic averaging, cumulative-logit
  coefficient direction and observed-score covariance units. Numerical
  calculations are unchanged.

- Shorten the `fit_mfrm()` entry help around inputs, results and a first
  analysis. Detailed argument definitions, equations and estimator-specific
  conditions remain available in the installed `fit_mfrm_technical_details`
  reference. Fitting behavior and defaults are unchanged.

- Refine the shared plotting theme with a white canvas, lighter grids and
  more spacing. Axis and legend text now scale with ordinary ggplot theme
  edits. Existing presets and palettes remain available; saved statistical
  plot data are unchanged. The visual guide explains customization, export
  dimensions and saving edited ggplot objects.

- Wright-map group distributions now use a common, adjustable density bandwidth
  and retain Person-level rug marks, empirical CDF values, counts and exclusion
  reasons. Conflicting group labels are rejected, and the main map's Person
  exclusions also apply to the group display. These are descriptive distributions
  of fitted estimates, not latent population densities or DIF/DFF tests.

- Correct unified report text and exported figure descriptions that labelled
  saved profile-likelihood intervals as normal intervals. Mixed-method reports
  retain each interval's method and distinguish source estimates from
  reoptimized profile centers; estimates, bounds and availability are unchanged.

- `ordinal_logistic_table()` now accepts two or more groups, retaining the
  existing two-group results. It separates omnibus uniform, nonuniform and
  joint tests from reference/focal point odds ratios, and records each group's
  counts and ranges. All groups must share a positive-width matching range;
  repeated Person-target ratings remain unsupported. The extension does not
  establish small-sample error rates or account for matching-score uncertainty.

- `mantel_haenszel_table()` now offers explicit `method = "gmh"` for matched
  multi-category and multi-group association tests. It retains declared
  categories, group counts, excluded responses and unavailable tests; binary
  MH defaults are unchanged. GMH requires asymptotic inference and independent
  Persons within each target, and does not estimate a severity difference or
  a common odds ratio.

- `mfrm_facet_intervals()` now supports native-scale main effects and contrasts
  in RSM/PCM MML models with one fitted two-facet interaction. Intervals use
  full joint covariance, retain dependent levels and numerical refusals, and
  distinguish equal-level average severity differences from criterion-specific
  differences. This does not establish small-sample coverage.

- Simplify package and reference help around analysis questions, inputs and
  results. Keep worked procedures in existing tutorials and calculation/uncertainty
  specifications in specialist help. Consolidate data review, screening, plotting,
  intervals, reporting and saved-analysis guidance, preserving model conditions.
  Remove development-stage language from comparison guidance; clarify CML scope,
  diagnostic screening, conditional uncertainty and displacement calculations.
  Repair PCA and workflow plot-data examples. Align plotting and saved-posterior
  instructions with the current interface. Measurement, observed-score planning
  and rater feedback remain optional analysis branches.

- Reject diagnostics from a differently sized fitted dataset with the intended
  mismatch message, without unrelated factor-comparison warnings.

- Speed up continuous Person intervals in fitted-object and portable-calibration
  scoring by reusing integrated probability masses to bracket the endpoints.
  Integration tolerances, tail checks and scoring qualification are unchanged.

- Speed up `score_mfrm_persons()` for ordinary RSM fits by reusing the native
  response-likelihood calculation. Continuous integration, interval checks,
  cut-score probabilities and calibration-conditional uncertainty are retained.

- Refresh the bundled example diagnostics with the affected-response counts
  used by the ZSTD safeguards, retaining their fitted models.

- Consolidate plotting documentation around public `plot()` methods and the
  migration guide. Internal renderer topics no longer appear in installed help;
  view controls, interpretation and available result tables remain documented.
- Use reproducible saved synthetic fits in centrality and Person-scoring examples
  so readers can focus on the reported quantities without recalibrating models.

- Direct RSM/PCM MML with an intercept-only normal population and direct
  anchors on the noncentered facet can now estimate population SD. Fixed or
  adaptive integration is supported. Model-based facet/population intervals
  can include supplied independent-source anchor covariance, including its
  contribution to estimated log variance. Shared Persons, drift, source-to-new
  population mean contrasts and profile/sandwich extensions remain outside
  this route; numerical checks do not guarantee interval coverage.

## Migration from 0.2.4

- `fit_mfrm()` groups advanced model choices in `spec` and numerical settings
  in `control`. Existing flat calls and `do.call(fit_mfrm, ...)` remain
  supported when arguments use their complete names. Only the first four
  positional arguments (`data`, `person`, `facets`, `score`) retain their old
  positions; name later arguments. Abbreviations routed through `...`, such
  as `quad_point`, must become `quad_points` or
  `control = list(quad_points = ...)`. See `?mfrmr_api_migration_030`.
  A named `spec` can be reused across `describe_mfrm_data()`,
  `review_mfrm_anchors()`, `fit_mfrm()` and `export_mfrm_model()` within their
  respective supported scopes; a valid native specification need not be
  exportable. Detailed equations, identification rules and return-object
  provenance are in `?fit_mfrm_technical_details`.

Because this development series removes 35 overlapping public entry points in
favor of arguments and standard S3 methods, this interface-changing release is versioned as
`0.3.0`. Update scripts using the former names before upgrading:

- Use `mfrm_cluster(x, k, method = "pam")`, choosing `"pam"`,
  `"hierarchical"`, or `"kmeans"` as the method.
  The default remains PAM; method-specific geometry, settings and result
  classes are unchanged. The three method-specific names are no longer exported.
- Use `update(calibration, status = "validated")`, choosing `"validated"`,
  `"frozen"`, `"superseded"`, or `"retired"` as the target state.
  The four former transition functions are no longer exported. Validation
  requirements, permitted transitions and immutable frozen records are retained.
- Use `normalize_conquest_overlap_tables(..., input = "files")` for extracted
  table files; `input = "tables"` keeps the existing data-frame route. The
  former `normalize_conquest_overlap_files()` is no longer exported.
- Use `mml_quadrature_sensitivity()` for supported MML families and
  `read_facets_results()` for FACETS Scorefiles and fit tables (pass the
  returned `mfrm_facets_results` object or file path directly to
  `facets_fit_review()`). The former `gpcm_mml_quadrature_sensitivity()`,
  `import_facets_fit_table()`, and `read_facets_fit_table()` names are no
  longer exported.
- Use `import_irt_fit(fit, model = ..., source = "auto")` to import supported
  `mirt`, `TAM`, and `eRm` fits; `import_mirt_fit()`, `import_tam_fit()`, and
  `import_erm_fit()` are no longer exported separately.
- Use `mfrmr_output_guide()` with `scope = "intervals"`, `"visual_templates"`,
  `"gpcm_capabilities"`, or `"calibration_capabilities"` for the static route
  and capability matrices formerly exported as `mfrmr_interval_guide()`,
  `visual_reporting_template()`, `gpcm_capability_matrix()`, and
  `mfrm_calibration_capabilities()`.
- Read `fit$config$anchor_review` (or call `review_mfrm_anchors()`) and
  `diagnostics$precision_review` (or call `precision_review_report()`) directly;
  the `anchor_review()` and `precision_review()` wrapper accessors are no longer
  exported.
- Use `analyze_dff()` for the former `analyze_dif()` alias, and
  `review_mfrm_imputations(..., impute_ids = ids)` for the former
  `mfrm_response_imputations(..., impute = ids)` alias. These aliases are no
  longer exported; the latter reviews supplied completions and does not impute.

See `help("mfrmr_api_migration_030")` for exact replacements. Existing saved
result classes (including legacy `mfrm_design_evaluation_plot_data` and
`mfrm_signal_detection_plot_data` objects) retain their `summary()`, `plot()`,
and `as_ggplot()` methods. These interface changes do not change the underlying
statistical targets or inference eligibility.

## New capabilities and fixes

- Withhold standardized Infit/Outfit (ZSTD) when a diagnostic scope contains
  finite JML optimizer proxies for infinite primary Person estimates, or when
  the selected standardization has fewer than one degree of freedom. Such
  values could otherwise become spuriously large and positive. Mean squares,
  row counts and point estimates retain their existing definitions; affected
  summaries explain the missing ZSTD. Recompute diagnostics from saved fits
  without refitting before reusing affected older reports.

- Adaptive quadrature uses the original R calculation if the native posterior-mode
  search cannot bracket a root or converge. An unfinished search no longer returns
  a finite endpoint as a successful mode; non-finite evaluations remain explicit errors.

- Numerical regression checks now allow explicitly bounded floating-point
  differences when adaptive moments or sampling covariances are recomputed
  with alternative BLAS/LAPACK implementations. Diagnostic decisions and
  saved-object identity remain checked separately; estimation algorithms and
  convergence criteria are unchanged by this portability fix.

- A website-only ICNALE GRA essay-rating workflow connects licensed-data
  preprocessing, RSM/PCM and GPCM comparisons, local rater feedback and a
  fixed-CAF observed-score assessment plan. It includes common-scale location
  and descriptive Pathway plots, and bounded TAM/ConQuest comparisons.
  The source repository supplies code;
  corpus data and individual results are not redistributed or run during
  package checks or website builds.

- Direct RSM/PCM and GPCM MML can now use a bounded curvature restart for
  up to 200 free parameters, including adaptive integration and estimated
  normal populations. This addresses relative-likelihood stopping before
  the gradient criterion is met. Positive curvature, objective and gradient
  checks remain required; numerical accuracy thresholds are unchanged.

- RSM/PCM location intervals no longer reject every target solely because
  individual raters leave categories unused. All categories must remain
  supported in each fitted step scope, and the full joint information,
  stationarity, boundary and finer-quadrature checks must still pass.
  Local sparse-category warnings remain visible; available Wald intervals
  are model-based approximations, not a finite-sample coverage guarantee.

- Individual two-family GPCM rater sheets now distinguish saved limiting
  estimates from full-model working values when no finite maximum exists.
  Unidentified or unevaluated targets are shown as unavailable; binary step
  offsets fixed by constraint remain zero. Sheets retain the limiting location
  reference, distinguish original and retained exposure, and do not transfer
  full-model intervals or residuals to the limiting analysis. Recreate older
  sheets from their saved results to obtain the corrected reporting.

- D-Study plan differences from `confint(..., method = "mean_square_region",
  reference = ...)` now retain shared error-component cancellations. Bounds
  can be narrower under the same normal-model simultaneous confidence region;
  point estimates and individual-plan intervals are unchanged. Diagnostics
  retain the refinement and any reason for using the preceding enclosure.

- Joint RSM/PCM interaction-cell intervals now support direct MML with an
  estimated normal population (including latent regression) and adaptive
  integration. Full-information covariance includes estimated population
  parameters. Existing fixed-standard-normal EM/hybrid support is retained;
  numerical eligibility does not establish small-sample coverage or rater bias.

- `mfrm_facet_intervals(..., method = "profile")` adds an experimental,
  explicitly requested pairwise contrast for direct RSM/PCM MML with an
  estimated normal ability population. Residual variance may reach zero and
  is reoptimized throughout the interval. Source and profile estimates remain
  separate; quadrature retries, refusals and unqualified coverage follow
  saved plots and reports. Adaptive integration handles concentrated Person
  likelihoods, while weighted matrix products reduce repeated curvature work;
  the zero-variance calculation and numerical acceptance criteria are retained. This does not replace the default Wald procedure
  or supply a general small-sample correction.

- Two-family GPCM `confint(fit, parm = "steps")` now supplies experimental
  pointwise intervals for sum-zero step offsets using full joint marginal
  information. Dependent steps retain cross-covariances; binary fixed steps
  have zero SE by constraint and no inferential interval. Numerical refusals,
  coverage qualifications and step identities persist in plots, reports and
  exports. These are not location-inclusive category boundaries or ordering tests.

- Two-family GPCM summaries and reports now distinguish numerical interval
  availability from unestablished coverage, and identify Person posterior SDs
  as conditional on fitted calibration. The scope guide explains target-specific
  uncertainty and the constraint transformation used for step intervals. This
  clarification does not change existing estimates or intervals.

- Facet count reviews now separate distinct persons from rating rows and show
  per-person exposure, marginal facet crossings and category counts, including
  unused categories. Numeric count bands replace qualitative console/plot
  labels; existing `N`, `SampleCategory` and threshold names remain supported.
  Two-family GPCM summaries and reports expose the same descriptive counts,
  separately for full-model and limiting-fit data. Counts do not certify
  precision or a sufficient sample size.
- Two-family GPCM reports point estimates separately from SE and interval
  eligibility. Converged finite solutions populate `PrimaryEstimate` with
  their numerical status; automatic integration checks no longer require
  a full covariance matrix and do not certify local curvature. All-low or
  all-high levels of any facet identify a sufficient condition for no finite
  joint maximum, independently of location centering. Supported boundary cases
  retain separately fitted limiting estimates for the remaining targets,
  including rater-owned steps and an extreme criterion without slopes.
  The reduction requires one affected facet outside the reference slope family,
  at least two remaining levels and retention of all other facet levels.
  When the affected facet owns steps or has no slopes, full-model optimization is skipped;
  unidentified slopes/steps of removed levels and binary steps fixed at zero have distinct
  status in summaries and saved reports. Deliberately skipped optimization is
  reported as `not_run`, separately from numerical failure. Boundary summaries
  direct curves to an available limiting fit and retain unavailable intervals. Numerical
  failures retain optimizer traces. Summaries, result tables and reports
  preserve these distinctions; older saved fits are not silently upgraded.
- Experimental two-family GPCM MML now accepts two slope facets among
  additional fixed additive facets. In `spec`, `slope_reference` names the
  geometric-mean-one family, independently of `step_facet` and
  `noncenter_facet`. The original full-predictor response equation, fixed
  N(0,1) population and experimental inference limits remain unchanged.
  The supported path connects fitting, checked intervals, curves, response
  diagnostics, reports and conditional new-Person scoring. Extended portable
  calibrations use format 7; original two-facet format 6 artifacts and legacy
  fits retain their meaning. Stan export from an existing fit now preserves
  extended roles through posterior reload, contrasts, intervals, plots and
  reports. Log-slope priors remain tied to the source slope reference;
  changing that reference with fresh priors changes the Bayesian model.
  Two-family simulation, JML, direct data/spec export, brms and new-Person
  posterior scoring remain unsupported. This extension is included in the
  planned `0.3.0` release.
  Numerical convergence does not establish a unique finite optimum or
  adequate estimation precision for every rating design.
- Unanchored two-family Stan bundles omit per-rating `eta` and
  `effective_slopes` from saved output. These are deterministic work vectors;
  parameter draws, `log_lik`, `y_rep`, diagnostics and posterior replay remain
  available. Previously saved bundles retain their replay path.
- Reloaded two-family fits now load the sparse-matrix arithmetic methods before
  producing a rater report in a fresh R session. Location constraints in
  reports and experimental intervals follow `noncenter_facet`, independently
  of which facet owns slopes or category steps.

- Estimated-normal-population RSM/PCM fits now record a conditional analytic
  identification check in
  `fit$data_review$estimability$normal_population_identification`. It uses
  the constrained design of Persons with repeated ratings and states its
  assumptions and unsupported cases explicitly. This separates identification
  evidence from numerical audit limits without changing interval eligibility
  or formal-inference restrictions; see `?fit_mfrm_technical_details`.

- Plotting help and supplementary examples now describe the `ggplot` return
  value consistently and use `plot_data()` to retrieve tables and notes.
- Threshold ladders now draw the stored step thresholds instead of an empty
  display. Missing steps remain unconnected, and `highlight_disorder = FALSE`
  is preserved when saving and redrawing a plot.

- Information plots now draw only the requested information or conditional
  SEM series. `type = "both"` uses separate panels and vertical scales;
  numerical values and the complete `plot_data()` payload are unchanged.
  The visual tour explains how these exposure-weighted curves differ from
  per-Person information, separation reliability, agreement and G/Phi.

- Five core tutorials remain available through `browseVignettes("mfrmr")`:
  workflow, model choice, G theory, visual diagnostics and reporting.
  Specialized tutorials are provided on the documentation website for the
  corresponding version; their source and worked calculations are retained.
  Help pages and output guides link to these articles directly instead of
  suggesting unavailable local `vignette()` calls. The workflow help explains
  how to distinguish published website documentation from development help.

- `recommend_mfrm_design()` now accepts a saved paired D-study comparison
  with scenario-keyed costs, resource limits and G/Phi/SEM point targets.
  `summary(result, planning = TRUE)` separates feasibility, point attainment
  and paired interval direction, retaining every plan and unavailable value.
  Reports, plots and appendix exports keep the review and its assumptions.
  The observed-score two-crossed-facet scope is unchanged; no winner or
  selection-adjusted guarantee is supplied. Simulation recommendations retain
  their existing behavior and do not use these new planning arguments.

- Posterior scoring now checks predictive category probabilities and
  calibration variation in expected scores at both quadrature orders,
  including Persons with no observed ratings. Predictive rows carry their
  own `IntegrationMaxDiff`, `Status`, and `Reason`; inaccurate predictions
  are withheld while valid ability summaries remain available. Increasing
  `quad_points` can resolve a failed predictive check.
- The introductory RSM example explicitly estimates population mean and SD
  and reports the current inference restriction for that route. Help explains
  the substantive location/spread constraints in the unchanged fixed N(0,1)
  defaults, CML BIC's sensitivity to extreme-score Person counts, and the
  distinction between numerical readiness and finite-sample validity.
- Optional CmdStanR execution and test dependencies (`cmdstanr`, `jsonlite`,
  `posterior`) are declared, including the CmdStanR repository. A required
  release verification command records actual live Stan execution and rejects
  skipped tests; see `CONTRIBUTING.md`.


- `export_mfrm_model(..., format = "stan")` now accepts optional
  `run_cmdstan = TRUE` and `cmdstan_args = list(...)` when `cmdstanr` and a
  configured CmdStan toolchain are available, compiling and sampling via
  `inst/model-export/run-stan.R`, saving `fit.rds`, `predictive-draws.rds`,
  `posterior-summary.rds`, and `posterior-draws.rds`, and returning the
  validated `mfrm_posterior_summary` on `$posterior` (failing closed before
  writing files when `cmdstanr` or CmdStan is absent). In addition,
  `predict(object, newdata, ...)` (`predict.mfrm_posterior_summary()`) scores
  out-of-sample Persons for calibrated ordinary RSM/PCM and one-family GPCM
  bundles (fixed standard-normal populations; also estimated intercept-only
  normal populations for one-family GPCM in Stan), separating
  calibrated-posterior marginalization (`conditioning = "posterior"`, with
  law-of-total-variance `WithinSD`, `BetweenSD`, `SD`, and
  `CalibrationVarShare`) from fixed posterior-mean plug-in scoring
  (`conditioning = "plugin_mean"`), latent ability intervals from future rating
  posterior predictive distributions (`type = c("ability", "predictive",
  "both")`), and optional cut-score classification probabilities (`cut_scores`),
  with downstream `summary()`, `confint()`, `as.data.frame()`, `plot()`,
  `mfrm_report()`, `apa_table()`, and `export_summary_appendix()` support.

- Exported Bayesian summaries (`summarize.R`) now save a checksum-bound
  chain-by-iteration draw array (`posterior-draws.rds`) alongside
  `posterior-summary.rds` without requiring the `posterior` package during
  replay. `read_mfrm_posterior(path, draws = c("auto", "require", "none"))`
  validates free-to-constrained affine/Helmert basis identities and connects
  saved joint draws to shared-calibration facet, step-transition, and slope
  contrasts, arbitrary equal-tailed credible levels, empirical draw views, and
  reports across `summary(x, contrasts = ..., level = ...)`, `confint()`,
  `as.data.frame()`, `plot(..., quantity = "contrast", view = c("intervals",
  "draws"))`, `mfrm_report()`, `apa_table(..., which = "contrasts")`, and
  `export_summary_appendix()`, while keeping lightweight summary-only
  `posterior-summary.rds` workflows unchanged. Operations requiring joint draws fail explicitly
  when those draws are absent.

- All `plot()` S3 methods now render with `ggplot2` (`ggplot2 (>= 3.4.0)` in
  `Imports`) under a unified scientific theme (`preset = "standard"`,
  `"publication"`, `"compact"`, or `"monochrome"`) and invisibly return a
  `ggplot` object (`draw = FALSE` builds the `ggplot` without opening a
  graphics device). Extract underlying tables, matrices, and figure cautions
  with `plot_data(p)` and `plot_data_components(p)` instead of `p$data$...`;
  `plot_data()` and `as_ggplot()` preserve backward compatibility with saved
  `0.2.x` `mfrm_plot_data`, `mfrm_design_evaluation_plot_data`, and
  `mfrm_signal_detection_plot_data` RDS objects. Two-continuous-variable
  scatterplots—including the new `plot(dif, type = "scatter")` route after
  `analyze_dff()`, paired and extended model comparisons, fit pathway, fair
  average, bias, PCA scores, displacement, person-fit, and recovery
  scatterplots—support aligned X/Y marginal histograms via
  `marginal = c("auto", "histogram", "none")` computed strictly over the finite
  plotted pairs (`"auto"` enables histograms for single-panel views with at
  least eight finite pairs while keeping `unavailable` / fixed-anchor exclusions
  explicit).

- `export_mfrm_model()` now accepts long-format rating `data` and a model
  specification list `spec` directly, without running native MML optimization,
  Hessian evaluation or EAP scoring first. Ordinary RSM/PCM require a fixed
  standard-normal ability population in both Stan and brms. One-family GPCM
  additionally supports an estimated intercept-only normal population in Stan;
  its fixed-standard-normal specification supports Stan and brms. Supported
  fixed non-Person anchors and group means are retained. Population regression
  and non-Person facet interactions remain outside the export scope. Existing `mfrm_fit`
  export remains unchanged, and `write_json = TRUE` (`format = "stan"`) adds
  CmdStanR-compatible `model-data.json` alongside `model-data.rds` and CSV
  label/category maps without requiring `cmdstanr` or `jsonlite`.

- Compiled C++ (`cpp11`) evaluation kernels now accelerate fixed-grid and
  adaptive Gauss-Hermite MML estimation and information/covariance routines,
  2D Gauss-Hermite quadrature in `fit_mfrm_testlet()`, conditional JML, CML
  elementary symmetric functions, and marginal fit expectations, while sparse
  co-rating aggregation removes $O(R^2)$ table expansion in
  `describe_mfrm_data()` and `diagnose_mfrm()`. Estimands, quadrature rules,
  convergence criteria, CML existence checks, GPCM profile-likelihood curvature
  checks, and pure-R reference evaluation
  (`options(mfrmr.use_cpp11_backend = FALSE)`) are unchanged.

- FACETS terminology, coverage, visual and fit-df guides now use
  `facets_guide(topic = ...)`. Their explanations and return structures are
  retained. The legacy `mfrmRFacets()` name is replaced by `run_mfrm_facets()`;
  use the same arguments. See `?mfrmr_api_migration_030` for script migration.

- Reports now recognize the historical `JMLE` label as the same estimator as
  `JML` when matching a saved fit to freshly computed diagnostics. Estimates
  are unchanged, and diagnostics from a different fitted analysis still fail.


# mfrmr 0.2.4

**Published on CRAN: 2026-10-30.** This release connects rating analysis, assessment
planning and reporting. You can move from a fitted model or rating plan to
figures, tables and an editable report while keeping its assumptions visible.

## What is easier now?

- **Turn an analysis into a report.** `mfrm_report()` now brings generalizability
  and decision studies (G/D studies),
  anchor/linking reviews, individual or batch bias screens and model-choice reviews into
  the shared reporting workflow. Editable templates help you write the Methods
  and Results, and exports keep the underlying tables and analysis. You can
  record what high scores mean alongside the results.
  See [reporting examples](https://github.com/Ryuya-dot-com/mfrmr/blob/v0.2.5/vignettes/mfrmr-reporting-and-apa.Rmd).
- **Compare plans before collecting more ratings.** Examine how specified
  rater/task counts affect the precision of a score, weighted composite or
  difference. G/D-study reports connect sources of variation to plan tables,
  figures and supported intervals. Simulation recommendations now show their
  target, criteria and candidate outcomes in summaries, plots and reports.
  See [assessment planning](https://github.com/Ryuya-dot-com/mfrmr/blob/v0.2.5/vignettes/mfrmr-gtheory.Rmd).
- **Find rating patterns that need attention.** Check assignment links before
  fitting, examine observed agreement, and prepare supported rater feedback.
  Network and residual analyses connect to tables, plots and reports; their
  exploratory findings prompt investigation, not automatic rater exclusion.
  See [agreement](https://github.com/Ryuya-dot-com/mfrmr/blob/v0.3.2/website-archive/articles/mfrmr-rating-agreement.Rmd)
  and [networks and reporting](https://github.com/Ryuya-dot-com/mfrmr/blob/v0.3.2/website-archive/articles/mfrmr-networks-and-reporting.Rmd).

Continue a default single-family GPCM analysis in Stan while retaining its
estimated normal population and supported fixed anchors. `export_mfrm_model()`
writes the code, data and parameter maps with explicit priors.
RSM/PCM and fixed-population single-family GPCM also support brms;
two-family GPCM supports Stan. `read_mfrm_posterior()` brings saved summaries,
including population mean/SD, back to figures and reports with credible
intervals and sampling review. Population regression remains unsupported.
See [model export](https://github.com/Ryuya-dot-com/mfrmr/blob/v0.2.5/man/export_mfrm_model.Rd).

Supported MML models can show fitting progress and automatically increase
numerical integration accuracy up to a user-set maximum. Two-family GPCM
retains experimental inference limits.
The [GPCM guide](https://github.com/Ryuya-dot-com/mfrmr/blob/v0.3.2/website-archive/articles/mfrmr-gpcm-scope.Rmd)
and [detailed changes](https://github.com/Ryuya-dot-com/mfrmr/blob/main/inst/NEWS-0.2.4-details.md) explain these features,
calibration reuse, external-software comparisons and numerical improvements.

The [visual gallery](https://github.com/Ryuya-dot-com/mfrmr/blob/v0.3.2/website-archive/articles/visual-gallery.Rmd)
adds ggplot2 examples with interactive Plotly views, static figures and the
underlying values for rater feedback, assessment planning and score probabilities.

Source installation requires cpp11 0.5.4 or later for compatible C++ headers.

Complete worked examples remain in the help pages. Short primary workflows run
normally; supplementary examples show how to execute them explicitly with
`example("topic", run.dontrun = TRUE)`, including their printed results. Related short summary and plot methods are
documented with their generating function; their existing help aliases remain.

Model setup, residual PCA, Q3, residual plots, halo analyses and agreement
summaries avoid calculations that their requested results do not use. The same
input checks, numerical settings and interpretation limits apply.
Simulation specifications and design reviews also avoid repeated calculations,
with unchanged design tables, recommendation rules and numerical settings.
Testlet fitting and scoring avoid redundant calculations. Scoring also avoids
work for design patterns unused by the requested Persons. Quadrature, continuous
posterior integration and accuracy checks are unchanged.
RSM and PCM marginal maximum likelihood calculations reuse response
probabilities for repeated rating conditions on a shared quadrature grid.
Integration and convergence settings are unchanged.
Model comparisons reject stale or nonstationary MML solutions before expensive
uncertainty calculations, while retaining the existing eligibility criteria.

## Before reusing an older analysis

**Installing 0.2.4 does not update saved calculations.** Corrections affect
some fit statistics, bias/fair-score calculations, uncertainty and comparisons.
Keep your original inputs, fitted objects and settings, then use the
[updating guide](https://github.com/Ryuya-dot-com/mfrmr/blob/main/inst/NEWS-0.2.4-details.md#updating-saved-analyses) to identify
which results need rebuilding. `help("mfrmr_workflow_methods")` explains the
applicable checks; not every correction requires refitting.

Simulation files also check the package implementation and numerical environment
before loading or resuming a run. An implementation change can make an older
file incompatible even when its stored numerical values remain readable.
Retain the original package installation and environment to continue that run;
use a new file path when recomputing after an update. The default
`file_refit = "never"` rejects mismatches without replacing the saved file.

Use `plot()` on result objects instead of retired standalone drawing helpers;
`help("mfrmr_plot_migration")` lists replacements. Model comparisons now reject
likelihood-ratio tests when the larger model omits terms from the smaller one;
recreate affected comparisons from the existing fits.

The detailed change record is also installed at
`system.file("NEWS-0.2.4-details.md", package = "mfrmr")`.

# mfrmr 0.2.3.1

mfrmr 0.2.3.1 was a focused CRAN maintenance release. It did not change the
public R API or fitted-model contracts introduced in 0.2.3.

* Corrected the compiled-header configuration used under link-time
  optimization so all C++ translation units use R's configured `Rboolean`
  definition.
* Removed expired FACETS/Winsteps hyperlinks while retaining the substantive
  model distinctions and bibliographic references.

# mfrmr 0.2.3

mfrmr 0.2.3 improves the reliability and interpretation of the existing
many-facet workflows, with particular emphasis on the bounded GPCM introduced
before this release. No existing exported function has been removed.

## Bounded GPCM

* Clarified the supported GPCM as an aligned, single-owner relative-slope
  model. A slope multiplies the complete adjacent-category predictor, including
  the ability, facet, interaction, and owned-step terms. This is distinct from
  loading-only GPCM variants and multiplicative generalized MFRM families.
* Corrected the default scale identification for GPCM marginal maximum
  likelihood fits. When no population formula is supplied, an intercept-only
  population model now provides the latent location and scale while relative
  slopes use their documented geometric-mean constraint.
* Added parameter-level status for free GPCM slopes. Standard errors and
  confidence intervals are reported as ordinary inferential results only when
  the fitted solution and the parameter both satisfy the relevant readiness
  checks. Optimizer- and Hessian-based quantities remain available as clearly
  labelled diagnostic evidence when formal inference is not supported.
* Improved handling of extreme response patterns. Free JML Person estimates
  remain infinite for all-minimum or all-maximum responses; finite optimizer
  values are computational traces. Optional adjusted-score tables and plot
  placements do not replace the fitted estimates. Prior-regularized marginal
  EAP estimates remain available but do not override a blocked source fit.
* Added `gpcm_mml_quadrature_sensitivity()` for explicit same-data comparison of
  a fitted GPCM-MML result across quadrature grids. It reports changes in the
  marginal likelihood, relative slopes, raw observed-information quantities,
  population scale, and fitted probabilities without silently refitting during
  `summary()` or `print()`. Results support `summary()`, `print()`,
  `as.data.frame()`, and `apa_table()`.

## Estimation and readiness

* Added pre-fit checks for category and step support and strengthened
  constrained-estimability checks for sparse many-facet designs. Unsupported
  step contrasts and structurally unidentified comparisons now stop or fail
  closed before they can be presented as ordinary estimates.
* Strengthened numerical and boundary reporting for JML and GPCM fits.
  Optimizer convergence, terminal gradients, local rank, curvature, boundary
  evidence, and formal inference readiness are retained as separate concepts.
* Corrected the alignment of MML Person EAP estimates and posterior standard
  deviations when observations or Person rows are filtered during fitting.
* Fit-level readiness is propagated to Person results, plots, model-choice
  reviews, reports, and exported tables. A finite value from a blocked fit is
  no longer presented as an unrestricted inferential estimate.

## Summaries, diagnostics, and plots

* `print()` and `summary()` for fitted models now begin with a plain-language
  interpretation status, the availability of formal inference, the main reason
  for any hold, and a suggested next action.
* Formal-inference reporting now distinguishes satisfied fit-readiness
  requirements from the
  separate precision contract. A fit-only summary does not claim formal
  inference until matching diagnostics support the standard-error,
  confidence-interval, and reliability basis.
* GPCM uncertainty labels are consistent across fitted objects, summaries, and
  plots. Intervals calculated for observation-table Wright or facet displays
  are labelled `screening_only` and are not described as confidence intervals.
* PCM/GPCM model-choice reviews retain comparison warnings and keep information
  criteria, automatic preferences, and likelihood-ratio results unavailable
  when the compared fits do not share an eligible basis. Model-choice warnings
  can also be included in exported summary appendices.
* Differential-facet-functioning refits now replay the baseline response
  family, scoring range, estimator, weighting, and numerical settings. They
  fail closed for model structures that do not yet have a valid subgroup replay
  and linking contract. `min_obs` remains a computability guard rather than a
  claim about power or adequate sample size.
* Residual-PCA results now carry machine-readable exploratory-use guards.
  Returned components and warnings are descriptive diagnostics and do not
  automatically support dimensionality or subscore decisions.
* Corrected Person-involving bias-screen collection for explicitly requested
  facet pairs and expanded machine-readable FACETS feature-coverage guidance.
* FACETS table import retains the reported numeric text and displayed decimal
  precision so rounded output is not silently treated as hidden full-precision
  data.
* The Infit-versus-measure pathway now places Infit on the horizontal axis and
  the fitted measure on the vertical axis by default. GPCM category curves
  retain estimated step-facet slopes and now state
  explicitly that additive facet effects and interactions are fixed at zero.
  Multiple GPCM curve groups are shown in separate panels with category
  legends instead of as unlabelled overlaid traces.
* Native Wright maps label every retained facet level using displaced text and
  leader lines. Step thresholds are displayed as a vertical ladder with the
  score transition and fitted logit in each label.

## Reporting and reproducibility

* Reproducibility manifests describe inputs through semantic summaries such as
  class, dimensions, fields, and missingness instead of environment-sensitive
  serialization hashes. The `digest` package is no longer required by the
  distributed package.
* Reporting, APA tables, model summaries, and plots use the same readiness and
  estimate-use vocabulary. Added fields may affect code that assumes an exact
  number of rows or columns in a summary component; selecting fields by name is
  recommended.
* APA output is described as an APA/JARS-informed drafting aid rather than a
  compliance certificate. The reporting guide now lists study-level fields
  that the fitted model cannot supply, and weighted-method prose distinguishes
  row-level likelihood weighting from replicated Person response patterns.
* External programs are not required to install or use mfrmr. Comparisons with
  FACETS, ConQuest, TAM, or immer remain model- and estimator-specific and are
  not general interchangeability claims.

## Scope

* This release does not add unrestricted GPCM, loading-only or multiplicative
  generalized MFRM estimation, native multidimensional estimation, or new
  public multivariate G-theory and D-study functionality.

# mfrmr 0.2.2

* Standardized the package's canonical joint-maximum-likelihood label as
  `"JML"` across fitted objects, engine state, manifests, and replay scripts.
  `method = "JMLE"` remains accepted only as a backward-compatible input alias
  and now resolves immediately to `"JML"`.

* Revised first-contact guides and result guidance to use reader-facing
  wording while retaining documented API and status vocabulary.

* Clarified that `maxit` is a prespecified computational ceiling rather than a
  result-selection control. Iteration-limited fits now direct users to keep the
  specification fixed, follow a prespecified ceiling sequence, and withhold
  interpretation until the numerical-readiness criteria are satisfied.

* Replaced blanket `\dontrun{}` and `@examplesIf interactive()` guards with
  checkable examples or `\donttest{}` blocks. Only the two workflows that need
  separately generated ConQuest files remain `\dontrun{}`, and only the local
  Shiny viewer remains interactive-only. This makes most examples executable
  while keeping workflows that require separately generated files clearly
  marked.
- Corrected bounded-GPCM score-side delta-method uncertainty to use the
  expected-score derivative `ScoreSlope * Var`. `ScoreSideLogitSE` remains the
  logit-side component SE, while `ScoreSideSE` and its interval columns now
  apply `ScoreSlope * Var * ScoreSideLogitSE` on the expected-score scale.
- Refit DFF/DIF contrasts are now explicitly exploratory: separate-subgroup
  plug-in standard errors are labeled as conditional on baseline anchors and
  as omitting baseline-anchor uncertainty and cross-refit covariance. Refit
  rows no longer receive ETS A/B/C, formal-inference, or primary-reporting
  eligibility.
- Bias summaries and multi-pair bias collections now use `ScreenPositive` as
  the primary label and expose explicit screening-only eligibility metadata.
  Historical `Significant` names remain as compatibility aliases.
- `import_erm_fit()` now reads the current `eRm` `Person Parameter` /
  `Std.Error` schema as well as historical estimate labels, preserves usable
  person IDs, and rejects ambiguous or misaligned schemas instead of silently
  returning empty or recycled person rows.
- `q3_statistic()` and its print method now identify the result as mfrmr's
  standardized, Person-by-level aggregated-residual Q3-style screen. Legacy
  `YenFlag` names remain for compatibility, while fixed 0.20/0.30 rules are
  explicitly described as uncalibrated heuristics rather than raw-residual
  Yen Q3 critical values.
- `as_kable.apa_table(format = "pipe")` now appends an APA note once after the
  complete Markdown table. Previously the vectorized append could repeat the
  same note after every rendered table line.
- Added the package hex sticker.
- Tightened the FACETS positioning contract against the current 64-bit 4.5.1
  software target: coverage rows describe package-native surfaces, not
  external numerical equivalence, and mixed models, multiple scales,
  threshold anchoring, and fixed-calibration scoring remain outside 0.2.2.
- Corrected the `interrater_agreement_table()` documentation: `ExpectedExact`
  is computed from fitted category-probability vectors, not marginal-frequency
  chance agreement.
- Clarified that exact Person-by-facet duplicate rows are retained but place
  Data readiness under review; legitimate repeated ratings should carry a
  distinguishing event or occasion facet.

- Design, signal-detection, and population-prediction summaries now expose a
  deterministic named-facet review as `structural_design_review`. The review
  reports design balance, coverage, connectivity, and readiness without
  implying Monte Carlo performance or arbitrary-facet simulation support.

## Estimation performance

- A code-zero solution whose terminal gradient still requires review now
  triggers a bounded warm-started polish ladder when the portable tolerance
  setting is at least as strict as the public default. Each stage records its
  optimizer, portable setting, native L-BFGS-B controls when applicable,
  objective, terminal gradient, maximum parameter change, evaluations, and
  elapsed time; the best non-worsening stage is retained rather than assuming
  that stricter controls improve every fit monotonically.
- Direct, hybrid, and EM MML engines now apply the same terminal-gradient check
  to `InferenceReady`. EM relative log-likelihood convergence remains visible
  as an engine-specific stopping condition but no longer overrides the common
  numerical-readiness contract.
- `fit_mfrm()` now shares likelihood and analytical-gradient work at an
  identical parameter vector. MML direct and EM paths reuse quadrature
  probabilities and posterior quantities, while JML reuses category
  probabilities and stable observed log probabilities.
- The compiled cpp11 probability kernels are now the default for supported
  RSM/PCM MML work, with automatic pure-R fallback. Set
  `options(mfrmr.use_cpp11_backend = FALSE)` for an explicit reference-path
  comparison; GPCM continues to use its validated R kernel.
- `optimizer = "auto"` selects limited-memory L-BFGS-B for MML and for large
  JML parameter vectors; `"BFGS"` and `"L-BFGS-B"` remain explicit choices.
  The requested and actual methods are recorded for summaries, exports, and
  replay. The portable `reltol` setting is mapped to L-BFGS-B `factr` and
  `pgtol`; actual stage controls are recorded alongside the requested and
  selected-stage settings.
- Repeated likelihood calculations reuse work within each fit without
  enlarging saved fit objects or sharing cached values across parallel fits.
- Measurement-graph component detection now avoids repeated row-wise lookups
  while preserving the established subset labels and ordering. This reduces
  first-fit overhead for larger long-format rating designs.

## Summary workflow

- Fit summaries now separate Numerical, Data, Design, Stability, Diagnostics,
  and Reporting readiness. Disconnected measurement graphs and
  boundary-constant or single-level facet support remain explicit reporting
  holds even when numerical optimization succeeds.
- Wright, FACETS-style Wright, pathway, and related fit plots carry additive
  fit-readiness metadata. Review-only displays remain available for diagnosis
  but warn, mark their returned subtitle and drawn title, and do not silently
  promote availability to interpretability.
- `plot_apa_figure_one()` now emits one consolidated readiness warning per
  call, retains the readiness table and interpretation note on the composite,
  and visibly labels a non-ready result as a manuscript-oriented draft for
  review rather than a finished publication figure.
- Native and FACETS-style Wright maps now share a robust automatic range when
  boundary-separated facet levels are diagnosed. Exact estimates and CI bounds
  remain in the returned tables; ruler-end triangles, clipping metadata, and
  plot footers prevent truncated intervals from being read as complete, while
  the native returned legend uses the same keys as the rendered legend.
- `summary(fit)` now supports `profile = "fit"`, `"facets"`, and
  `"reporting"`. The default fit profile remains fast and does not compute
  diagnostics. The opt-in FACETS profile organizes fitted measures, fit,
  precision, categories, steps, and plot routes in a familiar reading order;
  it does not imply that FACETS was run or that its estimates are numerically
  equivalent.
- FACETS and reporting profiles can reuse a matching `mfrm_diagnostics`
  object. The returned summary records provenance and section availability,
  and `compute = "never"` prevents automatic diagnostic computation.
- Bias/DIF, residual PCA, and anchor-drift or linking analyses remain explicit
  follow-up decisions because their interpretation depends on the study
  design.
- `detail = "brief"` gives a selective console view without person
  identifiers. Full structured results remain available through the returned
  object.
- The concise summary presents the visual workflow in order: the required
  native Wright map with facet uncertainty and labelled step locations, the
  optional FACETS-style Wright ruler, and the optional Infit pathway. Person
  rows in the pathway remain opt-in.

## Examples and teaching data

- `example_operational` adds a reproducible 48-person teaching dataset with a
  connected two-rater assignment, moderate workload imbalance, and six
  planned omissions. It is the primary applied tutorial dataset;
  `example_core` remains an explicitly idealized complete-crossing
  example, and `example_bias` remains the planted-effect diagnostic example.
- `mfrmr_example_operational_design` declares the 288 planned assignment cells
  separately from the 282 observed scores. `describe_mfrm_data()` can compare
  an explicit `expected_design` with observed cells, report planned omissions
  and unexpected observations, review Person-facet graph components, summarize
  sparse links and duplicate cells, and keep person labels out of its default
  compact output. Without a roster, structural missingness is reported as not
  assessed rather than inferred from a hypothetical complete crossing.
- `list_mfrmr_data(details = TRUE)` now explains the design and intended role
  of every bundled synthetic dataset. The compact examples can be regenerated
  with fixed seeds. Combined-study
  objects now explain that relabeling prevents identifier collisions but does
  not establish a common scale without an explicit anchor/linking design.
- Precomputed vignette tables now follow the same successful operational MML
  route as the displayed workflow and record their source dataset, schema,
  MD5 checksum, and package version.

## Safer first analyses

- The public default remains `reltol = 1e-9` for the initial optimizer stage;
  bounded polishing is invoked only when `reltol <= 1e-9` and code zero
  precedes the terminal-gradient check. The fitted object records requested and
  selected-stage controls for replay. Model specification, design,
  identification, and inferential assumptions remain separate review
  questions.
- Non-finite scores or weights, blank person/facet identifiers, and fractional
  `maxit` or `quad_points` values now fail before expensive optimization with a
  focused correction. Duplicate Person-by-facet cells warn once per fit,
  report both affected rows and duplicate cells, and propagate a Data review
  state downstream.
- `missing_codes = TRUE` now applies the conventional sentinel set to scores
  while preserving person and facet IDs. An explicit character vector remains
  an explicit request to apply those codes across all selected model columns;
  the review records the scope used for each column.
- On-the-fly ConQuest overlap examples now use the same `1e-9` tolerance.
  Their bundle summaries, settings, written README files, and compact console
  summaries report the actual mfrmr fit controls, MML engine, terminal
  gradient, convergence state, and inference readiness. A fit requiring
  convergence review is clearly withheld from the external comparison step.
- `fit_mfrm()` now gives focused guidance for common undeclared missing-value
  codes and records score-category recoding in the fitted object. It also
  distinguishes an explicitly silent anchor policy from policies that report
  anchor review information.
- Partial-credit fits infer the step facet only when a familiar item-like role
  is unambiguous. Otherwise, the warning shows how to set `step_facet`
  explicitly. Rating-scale fits report that `step_facet` and `slope_facet` are
  not used.
- Direct data-frame input to `mfrm_results()` is limited to data with
  recognizable measurement roles. Ambiguous columns now lead to an explicit
  `fit_mfrm(..., method = "MML")` instruction instead of a guessed analysis.
- `describe_mfrm_data()` computes agreement automatically only when a
  rater-like facet is present. Agreement output names the facet actually used
  and avoids presenting a generic facet as a rater.
- Latent regression rejects a non-person-centered parameterization that would
  confound the population intercept with the measurement scale.
## Interpretation and compatibility boundaries

- Optimizer code zero is no longer treated as sufficient evidence of a clean
  solution when the terminal gradient remains large. Summaries label this
  state as requiring review and explain the diagnostic basis.
- `facets_feature_coverage()` and `gpcm_capability_matrix()` now present concise
  user-facing capability, limitation, and recommended-route information.
  Only documented user-facing columns are returned.
- FACETS-style plots reproduce a reading convention, not FACETS numerical
  estimation. ConQuest comparison helpers cover documented unidimensional MML
  overlap and do not automate ConQuest or claim general numerical equivalence.
- The generated ConQuest overlap command now states quadrature MML explicitly
  and requests four machine-readable CSV outputs.
  `normalize_conquest_overlap_exports()` reads those files, reconstructs the
  sum-constrained item location, trims fixed-width person identifiers, and
  prepares them for `review_conquest_overlap()`.
- The documented binary, item-only, one-covariate MML handoff was compared in
  a matched 31-node run with ConQuest 5.47.5 Demonstration Version. The result
  supports that narrow handoff and is not a claim of general numerical
  equivalence.
- `export_mfrm_results()` now labels every preset as a potentially identifying
  analysis archive, warns before writing unless the risk is explicitly
  acknowledged, and records privacy status in its summary, HTML index, and
  written-files manifest. Fit-level `export_mfrm_bundle()` archives follow the
  same warning and metadata contract, and the lower-level `export_mfrm_tables()` CSV
  writer now records per-file handling metadata. ConQuest overlap bundles
  likewise warn on file export and include an artifact-level privacy inventory
  for response, covariate, and case-EAP files.

## First-use workflow

- The recommended workflow is now data -> `fit_mfrm()` -> fit summary ->
  required Wright map -> focused diagnostics -> `mfrm_report()` or
  `export_mfrm_results()`.
- `summary(mfrm_results(...), view = "brief")` and
  `summary(mfrm_report(...), view = "reader")` provide stable, concise views
  over the corresponding structured objects.
- `export_mfrm_results(preset = "starter")` writes a reader-first result
  folder with an index, required Wright-map image, selected tables, report
  files, replay code, and a reproducibility manifest.
- The README, workflow vignette, and help pages now begin with the same compact
  analysis route and direct specialist questions to focused follow-up
  helpers.

## Wright maps and fit pathways

- Wright maps retain the native renderer as the default. Native maps can show
  facet SE or confidence-interval whiskers alongside fitted step locations with
  `show_ci = TRUE`, while fitted coordinates remain unchanged.
- `renderer = "facets"` adds an opt-in FACETS Table 6-style visual grammar:
  a shared logit ruler, person-frequency asterisks, signed facet columns, all
  fitted facet levels, horizontal score-transition lines, and optional rubric
  labels. The renderer reproduces a display convention, not FACETS estimation
  or numerical output.
- Both Wright renderers return tidy draw-free data for custom graphics. The
  native `top_n` display remains compact; the FACETS-style data retain every
  fitted location.
- `plot(..., type = "fit_pathway")` adds a separate fit-oriented display with
  Infit or Outfit on the x-axis and measure logits on the y-axis. Screening
  bands, measure intervals, and optional ZSTD companions are explicit.
- Person rows can be added to the fit pathway with bounded selection and
  independent person/facet label controls. The existing expected-score
  `type = "pathway"` is unchanged.

## Reporting and migration support

- `facets_term_crosswalk()` and `facets_visual_contract()` document the
  correspondence between FACETS terminology and mfrmr outputs while keeping
  visual compatibility separate from numerical equivalence.
- `plot_data()`, `plot_data_components()`, and `as_ggplot()` make plot
  coordinates, annotations, reference lines, and guidance available for
  custom R graphics.
- Plot helpers consistently support `preset = "monochrome"` for
  print-friendly figures.
- `export_mfrm_bundle(..., include = "html")` provides a fit-level
  HTML/CSV/replay bundle without first creating an `mfrm_results` object.
- Model-comparison output can be routed through
  `build_model_choice_review()` and `build_summary_table_bundle()`, with
  explicit guidance for equal-weighting RSM/PCM models, bounded GPCM
  sensitivity analyses, and latent-regression reporting.

# mfrmr 0.2.1

## Results, reports, and export

- `mfrm_results()` adds a comprehensive first-screen object for an existing
  fit, a `run_mfrm_facets()` result, or a long-format data frame. It gathers
  diagnostics, available tables, plot routes, status information, next
  actions, and reproducible code without replacing the lower-level helpers.
- `mfrm_results(include = ...)` supports purpose presets for publication,
  FACETS migration, validation, bias, local misfit, linking, network review,
  and bounded GPCM review.
- `mfrm_report()` converts an `mfrm_results` object into a navigable reporting
  plan. Its first screen, report index, template index, evidence boundaries,
  cautious wording, and next actions keep detailed tables available without
  turning diagnostics into pass/fail decisions.
- `export_mfrm_results()` writes selected result tables, report files,
  draw-free plot data, images, replay code, RDS output, and a written-files
  manifest. `export_mfrm_bundle()` remains the broader fit-centered archive.
- `launch_mfrmr_viewer()` provides an optional Shiny reader over an existing
  `mfrm_results` object. It displays stored results and does not refit the
  model or change diagnostics.
- `mfrmr_output_guide("public")` maps the shortest fit, results, report,
  viewer, export, and specialist routes. Additional guides cover FACETS,
  ConQuest, binary data, simulation, linking, response time, and R-first
  visualization.

## Interpretation and reporting accuracy

- APA output now describes mean-square fit relative to the selected screening
  band instead of labeling overall fit “acceptable” or “elevated.” Band
  position is presented as a review signal, not a validity decision.
- MML reports state that person measures are EAP estimates and that
  residual-based fit statistics are evaluated at those measures. Comparisons
  intended to match JMLE-based FACETS output should use `method = "JML"` and
  aligned settings.
- Small-df ZSTD values are withheld when the transformation is unstable.
  FACETS/Winsteps output may still show a value under different sparse-cell
  conventions; such pairs are labeled as availability or standardization
  differences rather than automatically as fit differences.
- MML person separation and reliability are based on EAP measures and
  posterior SDs. They are kept distinct from JMLE-based FACETS reliability and
  from observed inter-rater agreement.
- APA tables and narratives report the measure-CI basis and fitted sign
  convention when available. Separation reliability, agreement, fit, and
  validity remain separate reporting claims.
- `precision_review_report()`, `fit_measures_table()`, and
  `facets_fit_review()` expose the fit, ZSTD, df, separation, and uncertainty
  bases needed before drafting technical conclusions.

## Focused review and planning

- Result and report objects can carry explicit bias, local-misfit/pathway,
  linking/anchor, precision, network, and response-time sections. Missing or
  unrequested sections remain visible as such.
- Recovery summaries expose `reading_order`, `condition_review`, and
  fit/separation operating characteristics. Bounded-GPCM `slope_regime`
  labels and extended sensitivity evidence remain separate from recovery
  metrics, convergence, and uncertainty availability; they are not automatic
  adequacy decisions.
- Resampling and simulation tools add person-clustered subsampling/bootstrap,
  sparse linked rating designs, connectedness summaries, and peer-assessment
  assignment checks. These are stability, design, or operating-characteristic
  diagnostics rather than calibrated tests or automatic decisions.

# mfrmr 0.2.0

## Scope and compatibility

- This version strengthens mathematical identification, uncertainty
  reporting, diagnostic tables, recovery tools, and draw-free visual output
  for RSM, PCM, and the documented bounded GPCM implementation.
- Breaking change: former exported `*_audit*` helper names, compatibility
  classes, and duplicate output fields were removed in favor of the canonical
  `*_review*` names. Stable accessors include `anchor_review()` and
  `precision_review()`.
- `facets_positioning_guide()`, `facets_feature_coverage()`, and
  `facets_output_contract_review()` describe supported FACETS-style tables,
  migration routes, and known differences. mfrmr estimates remain
  package-native unless external FACETS output is supplied for comparison.
- `mfrmr_output_guide("facets")`, `mfrmr_output_guide("conquest")`, and
  `mfrmr_output_guide("r")` provide focused entry points for users moving from
  FACETS or ConQuest and for users who want reusable R plot data.
- `write_mfrm_residual_file()` and `write_mfrm_subset_file()` add standalone
  residual and connected-subset files for external review.

## Estimation and fit statistics

- RSM, PCM, and bounded GPCM step profiles now use the correct sum-to-zero
  parameter count. MML structural covariance output provides uncertainty for
  non-person facets, steps, and bounded-GPCM slopes when the observed
  information is available.
- Measure tables record confidence level, interval method, eligibility, and
  interpretation basis. `compare_mfrm()` records the BIC sample-size basis,
  including weighted fits, and withholds unsupported likelihood-ratio tests
  with an explicit reason.
- Bounded-GPCM simulation and fitting use the same geometric-mean-one
  relative-slope identification. Expected scores, information, category
  curves, fair averages, and bias screening use slope-aware probabilities.
- `fair_average_table(fair_se = TRUE)` adds structural delta-method
  uncertainty where supported. `estimate_bias()` uses slope-aware information
  and can report conditional profile-likelihood screening quantities.
- `diagnose_mfrm(fit_df_method = "engine" | "facets" | "both")` exposes the
  package and FACETS-style df/ZSTD conventions separately.
  `facets_fit_review()` and `read_facets_fit_table()` support row-aligned
  comparison with existing FACETS tables without treating convention
  differences as estimation errors.
- `compute_person_fit_indices()` computes polytomous `lz` from observed
  category probabilities. Snijders-corrected `lz_star` is reported for
  compatible JML/fixed-effect person estimates and remains unavailable for
  MML/EAP scores. The incorrectly named `ECI4` output was removed; use
  `OutfitZSTD` for the corresponding standardized chi-square quantity.

## Diagnostics and visualization

- `fit_measures_table()` adds FACETS-style element fit tables, configurable
  threshold profiles, measure intervals, df-sensitivity summaries, and
  draw-free fit plots.
- `data_quality_report()` reports row retention, score-support gaps,
  zero/sparse category use by facet level, restricted response patterns,
  quality flags, original-to-internal score mapping, and dashboard plot data.
- `analyze_residual_pca(parallel = TRUE)` adds residual-permutation parallel
  analysis and dedicated plots. It remains exploratory dimensionality
  evidence.
- `category_curves_report()` adds category probabilities, cumulative
  probabilities, total information, category-specific information, boundary
  summaries, and overview/focused plots.
- `plot_data()` and `plot_data_components()` expose long-form data,
  annotations, styles, and settings from supported `draw = FALSE` plots;
  monochrome and interval guides support print-oriented reporting.
- Response-time QC, design connectedness, rater-effect networks, and halo
  screening receive dedicated summaries and plots. They remain descriptive
  evidence, not speed parameters, logit estimates, or automatic exclusions.
- `mfrm_d_study()` extends observed-score generalizability output to planned
  rater/facet-count comparisons. Its residual-scaling assumptions are
  reported explicitly; it is not a substitute for an unidentified
  interaction decomposition.

## Recovery, model choice, and reporting

- `evaluate_mfrm_recovery()` and `assess_mfrm_recovery()` report parameter
  recovery, convergence, coverage, Monte Carlo precision, uncertainty
  availability, score support, and user-specified practical thresholds in
  separate summaries and plots.
- `build_model_choice_review()` combines fitted-model comparisons, model-role
  guidance, downstream support, cautious wording, and optional weighting
  review for RSM, PCM, and bounded GPCM candidates.
- `build_summary_table_bundle()` and `export_summary_appendix()` accept a
  broader set of fit, recovery, person-fit, precision, and comparison objects
  for report and appendix handoff.
- DIF plots add comparable scales, value labels, flag thresholds, confidence
  intervals, and interpretation metadata. Input validation for DFF/DIF
  helpers now fails earlier with clearer messages.
- Citation and interpretation corrections clarify mean-square screening
  ranges, Q3 residual conventions, sample-size guidance, ICC bands,
  shrinkage uncertainty, and the limits of pairwise bias SE approximations.

## Bounded GPCM boundary

- `gpcm_capability_matrix()` is the authoritative support map. Supported,
  caveated, and unavailable routes include a recommended alternative and the
  evidence needed for broader use.
- Direct fitting, posterior scoring, information, category plots, recovery,
  fair averages, conditional bias screening, and selected reporting/planning
  helpers are available where marked.
- Full unrestricted discrimination structures, full FACETS score-side
  equivalence, posterior-predictive checks, and heavy Bayesian backends are
  outside this version's supported scope.
- Structured `mfrmr_gpcm_scope_error` conditions identify the unsupported
  area and recommended route instead of returning a partial result.

## Defaults and performance

- No defaults changed from 0.1.6:
  `quad_points = 31`, `diagnostic_mode = "both"`,
  `plot(fit)` showing the Wright map, and `keep_original = FALSE`.
- Users upgrading directly from 0.1.5 should note that these defaults were
  introduced in 0.1.6.
- The cpp11 MML backend is used by default for supported RSM and PCM work;
  `options(mfrmr.use_cpp11_backend = FALSE)` selects the pure-R reference
  path. Unsupported kernels fall back automatically.

## 0.1.6

- Changed the default diagnostic mode from legacy-only to both legacy and
  strict-marginal diagnostics, increased MML quadrature points from 15 to 31,
  and made the Wright map the default `plot(fit)` output. The former overview
  remains available with `type = "bundle"`.
- Added estimated facet interactions, empirical-Bayes shrinkage,
  hierarchical/sample-adequacy review, missing-code preprocessing, APA output
  adapters, confidence intervals across major plots, Q3 diagnostics, expanded
  person-fit indices, observed-score generalizability helpers, import adapters
  for mirt/TAM/eRm, resumable MML fits, and additional diagnostic plots.
- Improved fit summaries, replay scripts, input validation, examples,
  large-design diagnostics, and the printable cheatsheet.

## 0.1.5

- Simplified the first-use fit, diagnostic, and reporting workflow.
- Added MML latent regression with EAP scoring, the first bounded-GPCM
  fitting route, binary and non-consecutive score support, strict-marginal
  follow-up plots, report/appendix helpers, and clearer uncertainty and
  support boundaries.
- Added focused overlap and handoff guidance for FACETS, ConQuest, mirt, TAM,
  and eRm.

## 0.1.4 to 0.1.1

- Improved metadata, references, help-page examples, output documentation,
  and cross-platform portability while preserving the public analysis
  workflow.

## 0.1.0

- Introduced package-native many-facet RSM/PCM estimation with MML and JML,
  arbitrary facet counts, FACETS-style bias and fixed-width reports,
  APA-oriented summaries, residual-PCA diagnostics, visual summaries,
  anchoring helpers, and synthetic example data.

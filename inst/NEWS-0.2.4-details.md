# Detailed changes in mfrmr 0.2.4

This per-function record complements the concise release highlights in NEWS.
For saved analyses, start with `help("mfrmr_workflow_methods")`; for reporting,
use `vignette("mfrmr-reporting-and-apa")`.

Complete worked examples remain in each help topic. Short primary workflows run
with `example("topic")`; supplementary workflows show an explicit
`example("topic", run.dontrun = TRUE)` command. Long fitting and external-software
examples retain their existing opt-in status, inputs and numerical settings.
Executing help examples displays their tables and summaries through R's standard
example interface.

Fixed-quadrature MML fitting and observed-information calculations now aggregate
observation scores before projecting them onto facet, interaction, step and
slope parameters. This avoids repeated grouping at each integration point in
RSM, PCM and single-family GPCM, with the same likelihood, quadrature and
convergence settings. Floating-point summation order can differ; this does not
extend the supported models or their inference guarantees.

## Updating saved analyses

Keep the original objects and settings. Use the affected step below and the
central migration guide, `help("mfrmr_workflow_methods")`, before regenerating
results. A wording or display correction alone does not require refitting.

| Affected work | What changes / where to start |
| --- | --- |
| Standalone drawing functions | Use `plot()` on the fit or result. `help("mfrmr_plot_migration")` lists replacements for the 33 retired exports. |
| Fit, precision, fair scores, bias and reports | Corrections address fit df/ZSTD, PTMEA, anchor handling, fair-score targets and uncertainty/availability labels. Recompute affected diagnostics and derived results from a compatible fit; older fits lacking current estimation checks require refitting. |
| Model comparisons | RSM/PCM LRTs require every RSM interaction to be retained in the PCM. Recreate affected comparisons with `compare_mfrm()` from the existing eligible fits; information criteria and LRT eligibility are distinct. |
| G/D studies, ICC and design recommendations | Corrections retain missing-data accounting, covariance and metric-specific availability, negative-component policies, failed-run denominators, and full-precision recommendations. Follow the migration guide for whether to recompute the study or only its summaries. |
| Person scores, intervals and portable calibrations | Preserve the recorded prior, scoring algorithm and interval meaning. Updating labels does not replace old grid-based intervals; changing the algorithm requires rescoring or a new reviewed artifact. |
| External imports | Re-import affected source files/objects and rebuild dependent output. FACETS group flags, delimiters and uncertainty placeholders are retained; imported values are not converted into proof of native-model equivalence. |
| Renamed entry points | `simulate_mfrm_recovery()` replaces `evaluate_mfrm_recovery()`; `review_peer_design()` replaces `build_peer_review_design_review()`; `gpcm_capability_matrix()` consolidates the earlier GPCM scope helpers. |

## Changes by function and workflow

- Structural design summaries now count appendix tables separately for each
  preset (`"all"`, `"recommended"`, `"compact"` and manuscript sections).
  Previously, some counts and fractions pooled rows from all presets and
  repeated the pooled result under each label. This also affected their
  selection plots. Correct filtering removes redundant grouping work.
  Recreate affected summaries, plots and report bundles from the original
  design or simulation result; fitted estimates and simulation runs do not
  need to be repeated. A single-preset appendix export was already filtered.

- `build_conquest_overlap_bundle(fit = ...)` now records the quadrature-point
  count from the supplied fit, rather than the default used without a fit.
  Rebuild an affected bundle to correct its settings table; this does not
  require refitting or change the fitted estimates.

- ERGM results now support `coef()`, `vcov()` and `confint()` alongside
  `summary()`, `plot()` and reporting. Coefficient intervals retain unavailable
  states and conditional-network limitations. GPCM bootstrap and curve-interval
  results have dedicated summaries that retain unresolved trials and cautions.
  Facet bar plots reserve space for labels and display finite confidence bounds
  within the axis range, while restoring the caller's graphical margins.

- Diagnostic summaries print readable reliability-use labels while preserving
  the structured values and inference restrictions. `estimate_bias()` checks
  missing or invalid diagnostics before inspecting extrapolation settings,
  so incomplete inputs receive the relevant `diagnose_mfrm()` guidance.
  These corrections do not require refitting or recalculating estimates.

- Workflow and reporting guides connect assessment questions to design,
  measurement, agreement, G/D studies, linking, comparison and exploration.
  The editable GPCM report now correctly describes its estimated normal
  population mean/variance, exports those settings, and explains when to add
  pairwise kappa. Its previous fixed-standard-normal Methods wording was wrong;
  fitted estimates and algorithms are unchanged.
  The template now takes its inference decision from the assembled report,
  including the collected precision diagnostics. Stan export now retains the
  default one-family GPCM's estimated intercept-only normal population, with
  explicit mean/log-SD priors and connected posterior population summaries.
  Model-export errors explain
  when an estimated population is unsupported; switching to a fixed N(0,1)
  population changes the model rather than merely enabling a file format.

- The README and release highlights now lead with assessment tasks and link
  to specialist guides. An editable one-family GPCM report template combines
  model summaries, fit, figures and a separate observed ordinal-alpha analysis.
  Multipanel probability plots reserve additional axis/title space. Curve
  export help now correctly identifies the `probabilities` table.

- Marginal pairwise agreement diagnostics validate observation probabilities
  once and evaluate quadrature-node products together within each response
  context. Integration rules, diagnostic thresholds and unavailable-result
  handling are unchanged; no model refitting or migration of saved results
  is required.

- RSM-versus-PCM likelihood-ratio tests now require the PCM to retain every
  interaction term in the RSM. Previously, an RSM interaction omitted from the
  PCM could be overlooked by the nesting review. Such pairs retain eligible
  information-criterion comparisons, but no chi-square test is reported.
  Recreate saved comparisons with `compare_mfrm()` to apply the corrected
  review; the fitted models need no refitting for this correction.

- The GPCM capability guide now consistently describes both fixed-grid EM and
  adaptive direct fitting for two-family models, including opt-in automatic
  integration with a maximum point count. An obsolete EM-only statement has
  been removed. Fitting defaults and experimental interval qualifications
  are unchanged.

- New two-family adaptive direct fits cap the ordinary optimizer-gradient
  review tolerance at `N * 1e-6`, retaining any stricter existing tolerance.
  This lets the existing bounded polishing mechanism address an unfinished
  per-Person score before interval qualification. Candidate stages still
  require a non-worsening objective; starts, quadrature and subsequent
  information/precision criteria are unchanged. The rule is recorded in the
  fit, and older saved fits retain their original convergence review.

- After a successful two-family automatic-integration fit, slope and location
  intervals reuse its complete numerical qualification for the unchanged fit
  within the same R session. Only the latest qualification is retained, within
  the information-memory budget; changed fits, budgets or numerical backends
  require fresh checks.
  Warnings and refusals are not cached. Saved model objects are unchanged and a
  new R session recalculates the checks. No admission criterion is relaxed.

- `fit_mfrm(..., progress = TRUE)` displays fitting phases and elapsed time,
  with direct objective counts and detailed two-family EM, initialization and
  automatic-quadrature stages. It reports no estimated completion percentage
  or remaining time; the indicator is cleared on return, error or interruption.
  The default remains silent and saved calibration settings are unaffected.

- Two-family EM reuses response probabilities at exactly identical parameter
  vectors across marginal/complete likelihood and gradient calls. E-step
  counts remain fresh; numerical targets, starts, stopping and interval
  qualification criteria are unchanged. The GPCM guide now distinguishes
  paired numerical delivery, coverage, computation cost, small-sample stopping
  differences, and the scope of bias/anchor workflows.

- Automatic MML integration accepts `mml_accuracy_control=list(max_points=...)`
  as a hard ceiling including reference rules. Two-family GPCM now supports
  automatic refinement while preserving fixed-grid EM or adaptive direct
  fitting and its initialization policy. Person likelihood/moment precision
  and existing local-information checks remain explicit; cap exhaustion and
  other refusals retain the candidate and diagnostic history. Success does
  not establish global optimality or nominal interval coverage.
  Two-family automatic histories use `check_refused` when a prerequisite
  check stops refinement. The former label `nonintegration_refusal`
  did not identify the cause; a failed higher-order information check also
  stops refinement. Inspect the retained check table. Existing saved histories
  keep their original labels; the stopping rules and numerical results are unchanged.

- Two-family GPCM quadrature sensitivity now reuses the full information
  calculation for diagnostic standard errors and component intervals within
  each evaluated fit. Numerical checks, unavailable intervals, warnings and
  interval definitions are unchanged; no fitting or integration default changes.

- D-study interval tables now distinguish approximate F, which has documented
  small-sample undercoverage, from the experimental numerical projection of
  a conservative normal-model joint mean-square region. The method choice,
  point estimates, interval endpoints and availability are unchanged.

- MML help now separates overall fit readiness from target-specific interval
  eligibility and explains how integration choice affects downstream outputs.
  Supported native-location and RSM/PCM population intervals retain their own
  checks and coverage cautions. Fixed integration remains the default; its
  default order is not a guarantee of numerical accuracy.

- Automatic MML integration checks reuse each rule's likelihood and posterior
  calculations when evaluating gradients. Identical response/design/weight/
  population-covariate patterns retain separate Person checks. Accuracy
  tolerances, refinement orders, defaults and inference eligibility are unchanged.

- One-family GPCM results now accept matching saved `predict_mfrm_units()`
  scores through `mfrm_results(..., scores = ...)`. Reports and exports retain
  EAP, posterior SD, continuous conditional intervals, per-Person priors,
  numerical checks and omissions without rescoring. Calibration/population
  estimation uncertainty is excluded, and review-only results remain labelled.
  Attached scores are separate from ordinary diagnostic normal bands.

- Imported FACETS element plots accept a `label_map` keyed by facet and exported
  element ID, allowing local display names while preserving source labels,
  values and status. Unknown or duplicate keys are rejected. Help and the
  migration guide explain Japanese PNG/PDF output with an appropriate device
  and installed font; this does not require Japanese input to FACETS.

- FACETS Scorefile import now recognizes semicolon-separated output, including
  its facet heading line and `score.N.txt` filenames. Previously this layout
  could be misread as fixed-width text or fail to identify its facet/element
  columns. Separator detection ignores punctuation inside quoted headings.
  Reported numeric tokens, element IDs and duplicate labels are preserved;
  existing imported results need no refit.

- FACETS result plots now label source `Status=1` as a direct/group anchor
  flag, without implying that each flagged measure is individually fixed.
  `SourceStatus = "anchored"`, raw group codes, numerical values and
  availability are preserved. The summary count is now `AnchorFlagged`
  (previously `Anchored`). Replot saved imports to refresh their legend and
  interpretation notes; source files need not be reread. Scorefile group codes
  and SEs do not identify the anchor subtype; consult the original
  specification or Table 7. Extreme group members can instead carry Status
  -2/-3 and are not included in `AnchorFlagged`. Their finite reported
  locations do not identify an adjustment rule; their fit placeholders remain
  unavailable for plotting.
  Measure and category-curve plots also keep their legends separate from axis
  labels during saved-plot replay.

- Ordinary additive RSM/PCM JML now accepts group means with explicit
  `centering_scope = "anchor_origin"`. Direct or group anchors supply their
  facet's origin; other facets retain the declared centering convention.
  Default whole-facet and legacy remainder centering are unchanged. FACETS
  per-element group contributions must first be converted to a repeated mean
  target. The estimator retains constrained maximum likelihood; matching
  FACETS group totals does not guarantee matching every numerical result,
  extreme convention or uncertainty measure.

- `estimate_bias()` now retains directly fixed finite main-effect anchors in
  ordinary additive RSM/PCM JML when their observed scores are all minimum or
  maximum, including with the default `xtreme = 0`. Previously, requesting
  extreme-cell extrapolation could change whether a Person-by-facet cell was
  present. Unadjusted boundary cells still have unavailable screening SEs and
  tails. Recreate affected bias tables and collections from the saved fit and
  matching diagnostics; no main calibration refit is required. Extreme-score
  reporting now explains that the inward correction can reverse a bias sign
  relative to the observed score departure.

- Jeffreys JML now attempts bounded local finishing when BFGS stops before
  satisfying the existing penalized score criterion. It retains the objective,
  gradient tolerance, category support and anchor constraints, with finishing
  diagnostics available in `checks$polish`.
- Experimental Jeffreys delta intervals now document demonstrated small-sample
  undercoverage. In a six-response binary case with true probability 0.5,
  exact coverage of nominal 95% limits is 78.125%, despite correct local
  covariance and an unbiased point at that probability. Numerical availability
  does not qualify inference; `CIEligible` remains false. Exact basic-refit
  quantiles also under-cover in this special case (45.2709% at true probability
  0.1); increasing the number of refits cannot repair that interval limitation.

- `fit_mfrm(..., method = "JML", jml_penalty = "jeffreys")` now offers
  explicit canonical information-penalized estimation for identified additive
  RSM/PCM models. It preserves declared rubric categories, including unobserved
  categories and constant outcomes. A separate result exposes Person/facet/step
  points and numerical checks without inheriting ordinary SEs or FACETS
  diagnostics. Common-plan scores and experimental refit intervals use this
  estimator throughout calibration, anchor extraction and scoring, preserving
  shared response events. Experimental `method = "delta"` intervals use
  analytic penalized score and anchor derivatives, retaining shared-event
  covariance. Inverse penalized curvature is not treated as sampling variance;
  failed derivatives or conditioning retain points with unavailable uncertainty.
  Mixed-estimator pipelines are refused. Ordinary defaults and the generator-only option are
  preserved. Finiteness and numerical stationarity do not establish score-bias
  reduction, global optimality, incidental-parameter consistency or coverage.
- Ordinary additive RSM/PCM JML can finish stalled numerical fits while holding
  extreme-Person optimizer coordinates fixed. Sparse Newton updates retain
  every observation, the original likelihood and the full-gradient check.
  Primary extreme estimates remain infinite; this does not remove structural
  nonidentification, correct statistical bias or qualify interval coverage.
- Common-plan refit intervals now report Monte Carlo enclosures for both
  interval endpoints, with a separate `mc_level`. Binomial order-statistic
  bounds retain unresolved refits in the requested denominator and expose
  infinite limits when simulation precision cannot be bounded. The original
  empirical limits are unchanged. `MCFinite` describes finite enclosures,
  not adequate precision, nominal coverage or bias correction; plotting and
  saved results retain the diagnostic and its interpretation.
- Common-plan refit intervals offer explicit `generator_method = "jeffreys"`.
  This experimental option fits a canonical information-penalized generating
  model while retaining the original ordinary-JML points and every ordinary
  calibration/anchor/scoring refit. It can handle extreme generating Persons
  and original limiting expected scores without substituting finite reporting
  measures. Generating parameters, probabilities, targets and numerical checks
  are saved. Declared rubric support and qualified source fits remain required;
  failed generators and unresolved refits are not silently replaced. The default
  ordinary generator and delta/bounded methods are unchanged. Finite generating
  estimates do not establish score-bias reduction or interval coverage.
- Common-plan expected-score intervals now offer explicit `method = "bounded"`.
  This conservative alternative uses independent bounded scoring responses
  directly and retains the entire score range for unobserved plan cells.
  Extreme Persons do not require finite reporting replacements or a fitted
  generating model. The finite-sample confidence statement requires fixed
  exogenous assignments, plans and contrasts, and stable within-cell means.
  Results preserve the original JML points, identify uninformative bounds and
  show observed plan support. A separate `DesignEstimate` is available for
  complete plans. These bounds can be wide, do not correct JML bias and do not
  provide simultaneous coverage for data-selected comparisons. Refit remains
  the default; its and delta's experimental qualifications are unchanged.
  Both rubric endpoints must be explicitly declared when fitting, with unchanged
  full category coding. Observed extrema, missing range provenance and
  data-dependent category collapse are rejected; the known-support requirement
  cannot be satisfied by labelling sample extrema as declared bounds.
- Ordinary additive RSM/PCM JML can now finish a stalled optimizer with bounded
  sparse Newton steps using the original likelihood and constraints. Each
  accepted step improves stationarity without resolved likelihood loss;
  the original gradient criterion still applies. The optimization history
  records the steps, solve residual and failures. Extreme-person boundaries
  remain distinct from numerical stopping, and this repair does not qualify
  confidence intervals or correct estimator bias.
- `read_facets_results()` imports supported FACETS score exports for standalone
  measure, point-only Wright and fit plots through `plot()`. It retains element
  IDs, labels, source precision, anchors/extremes, unavailable rows and file
  identity. Shared-axis Wright maps require an explicit compatible-scale
  declaration. Extreme fit placeholders are withheld from fit plots. Import
  preserves source values and does not refit or supply new scoring/intervals.
  Named CSV/TSV supports declared encodings including Japanese labels; validated
  fixed-width imports require ASCII. A small synthetic FACETS export is included.
  An optional `graph_file` reads named CSV/TSV or validated fixed-width
  Graphfiles for category-probability (`type = "ccc"`) and expected-score
  (`type = "expected"`) plots. Source scale/category codes, raw coordinates
  and precision are retained; missing rows break lines and probability-sum
  diagnostics respect printed precision. Graphfile axes are separate from
  element measures, with no automatic category recoding or scale alignment.
- Ordinary additive RSM/PCM JML fair-average results accept `reference_plan`
  for Person expected scores under one common weighted rating/scale plan.
  `confint()` adds an experimental fixed-Person refit procedure for these
  scores and paired contrasts. Estimated anchors can be recalibrated in every
  replication, with explicit response-event maps sharing observations across
  calibration and scoring. Saved intervals retain unresolved trials,
  numerical checks, anchor draws and the uncertainty basis through summary,
  plotting and export. Existing FairM/FairZ targets are unchanged. The new
  procedure is model based and does not yet qualify finite-sample coverage,
  correct JML bias, predict future observed scores or establish fairness.
  `method = "delta"` adds an experimental deterministic alternative through
  the same estimation and centering maps. It retains target covariance,
  calibration propagation and shared-event cross terms, using sparse solves
  for the selected targets. It requires finite stage estimates and positive
  information; failures preserve points with missing uncertainty and reasons.
  Normal limits remain unqualified for coverage and are not clipped to the
  score support. The default remains `method = "refit"`.
  If an unselected Person has an extreme primary estimate in the generating
  fit, delta results now preserve finite selected points with unavailable
  uncertainty and a reason. Such Persons are not silently removed or assigned
  finite reporting measures to obtain an interval.
- Ordinary additive RSM/PCM JML now accepts
  `centering_scope = "anchor_origin"`: direct anchors set their facet's origin
  without an additional zero-sum constraint. Choose an unanchored, active
  `noncenter_facet`. Group-mean origins are also supported as described above. Sparse
  interior RSM/PCM cases match FACETS 4.5.0 under these constraints. Anchor
  review, saved settings and replay retain the convention. The default
  whole-facet and legacy remainder conventions are unchanged; supplied anchor
  values remain fixed rather than carrying external calibration uncertainty.
  Detected joint additive recession directions affecting calibration stop this
  route; an identified origin alone does not guarantee a finite maximum.
  When every retained response of the Person anchors is at the same declared
  minimum or maximum, this route now stops before optimization without an
  optional LP dependency. FACETS' finite output for this case is not reproduced
  as an ordinary JML estimate. Interior or opposite-boundary support from
  additional Person anchors is not rejected by this specific rule.
- Ordinary additive RSM/PCM JML fair-score tables now retain finite direct
  anchors in extreme-row filtering under both display conventions. With
  `totalscore = FALSE`, an all-minimum/all-maximum anchored element no longer
  incorrectly removes its observations from other facet totals or PCM scale
  mixtures. Fixed values and unanchored extreme-score adjustments are preserved.
  Recreate affected tables and plots using the saved fit and matching diagnostics.
- FACETS migration guidance now explains comparisons with TAM JML: explicit
  category support, separate bias/extreme-score corrections, Person MLE versus
  WLE, conditional SEs, diagnostic trimming and weight denominators. Matched
  sparse and partially anchored RSM/PCM regressions preserve ordinary-JML
  behavior; they do not add a TAM JML importer or certify general coverage.
  Guidance also distinguishes a TAM 4.3-25 free-intercept anchoring path that
  moves the supplied Person anchor or fails to converge from a correctly
  translated unanchored reference solution.
- Ordinary additive RSM/PCM JML diagnostics now report `PTMEAExpected`
  beside observed `PTMEA`, with its unadjusted moment reference, frequency
  count and availability reasons. The finite-count factor `(N-2)/N` follows
  Linacre's published correction and matches scoped FACETS 4.5.0 PCM output,
  including unequal replication weights. Fair-average tables retain these
  columns when supplied by matching diagnostics. Below frequency count two,
  or with unbounded calibration or invalid moments, the corrected reference
  is unavailable. It is a conditional reporting approximation, not an exact
  sample-correlation expectation or a test. Recreate diagnostics from the
  saved fit to obtain the new columns.
- Ordinary additive RSM/PCM JML reports now offer
  `fair_average_table(..., xtreme = .3, xtreme_method = "facets")` and the
  corresponding `plot(fit, type = "fair_average", ...)` options. This solves
  a corrected weighted-total equation using the administered calibration,
  retains primary estimates and fixed anchors, and labels conditional reporting
  SEs separately. The default average-score display convention is unchanged.
  `estimate_bias(..., xtreme = .5)` similarly supplies finite reporting values
  for extreme cells when the corrected root is inside the search limits.
  These routes exclude unanchored extreme calibration elements; group-mean
  anchors and fitted interactions are unsupported. Extrapolation is not
  original-MLE inference or a calibration-inclusive uncertainty method.
  Broader extreme-score anchoring and exclusion equivalence remains unverified.
- Bias tables retain continuous count-minus-one df for Student screening tails
  and add `DFDisplay` for the nearest-integer display comparison. FACETS 4.5.0
  fractional-count tails in checked examples match an invalid fractional
  extension of an integer-df finite sum, which can change a 0.05 screen.
  mfrmr retains the continuous Student distribution. New `ProbabilityReference`,
  `ProbabilityStatus` and `ProbabilityReason` columns record its basis and
  availability. Nonpositive count-minus-one values remain in `DFUnrounded`;
  reference df, display df and tail area are unavailable, while a finite bias
  estimate and conditional information SE can remain usable for description.
- `mfrm_d_study()` now offers `negative_components = "genova_algorithm"`
  or `"genova_ems"` for complete balanced univariate ANOVA sources with
  supported crossed or nested facets. The default `"raw"` is unchanged.
  Results retain original and effective components, and report the selected
  convention in tables and plot notes. EMS recalculates components using the
  source design before projecting future counts. `confint()` and interval-aware
  summaries/plots pair adjusted points with original-source population intervals,
  retaining `IntervalEstimate`, `IntervalBasis` and unchanged bounds. These are
  not recentered adjusted-estimator intervals; sampling diagnostics retain their
  original basis. Saved intervals must match the selected treatment. Supported
  population plan differences use `confint(..., reference = ...)`; adjusted-
  estimator delta comparisons remain unavailable.
- `mfrm_facet_intervals(interaction = "FacetA:FacetB")` now reports jointly
  fitted interaction-cell effects and prespecified contrasts, using the full
  marginal-information covariance. The bounded route supports one two-facet
  interaction in fixed-population direct RSM/PCM MML. Finite points, fixed or
  unavailable intervals, cell constraints and interpretation notes are retained.
  These pointwise intervals are distinct from `estimate_bias()` screens and
  do not imply simultaneous inference or a fairness decision.

- Fixed-SD RSM/PCM anchor-origin facet and population intervals can now include
  a supplied full anchor-source covariance with an explicit independent-sample
  declaration. Point estimates and conditional intervals remain visible. Direct
  anchor contributions and cross-covariance are included; saved population
  results reuse the propagation without refitting. This first-order method
  does not correct drift, dependent samples or uncertainty in the fixed SD.

- Direct RSM/PCM MML now accepts `population_sd` to fix the residual population
  SD explicitly while estimating population coefficients. With an intercept-only
  population and direct anchors on `noncenter_facet`, this allows a cohort mean
  and new rater severities to be estimated without centering that rater panel.
  Existing defaults are unchanged. The full constrained design is checked for
  identification; GPCM, group-anchor origins and interactions are outside this
  extension. Default population and facet intervals condition on the supplied
  SD and anchors; source-error propagation is a separate opt-in described above.
  Neither option establishes anchor invariance.
  Saved replay and quadrature refits preserve the fixed-SD specification.

- Fit-measure tables, summaries and reports retain the recorded calculation
  basis, distinguishing EAP plug-in residuals, JML optimizer coordinates and
  imported source statistics. Missing definitions remain unrecorded; finite
  values and screening flags are preserved.
- Reports now show separation, reliability and strata by facet, with model and
  fit-adjusted values, input uncertainty and missing-value reasons. Existing
  cross-facet ranges remain compatibility summaries, not overall reliability.
  Person separation based on EAP/posterior SD is not renamed as EAP reliability.
  Rebuild reports from saved results to obtain these explanations; no refit is
  required.

- `category_curves_report(fit, newdata = contexts)` now provides full-context
  curves and cumulative 50% thresholds for native one- and two-family GPCM
  MML. Fitted locations, interactions and effective discriminations are retained.
  Thresholds use analytic binary roots or stable bracketed root finding, preserve numerical precision, and
  distinguish adjacent intersections, centered steps and unavailable roots.
  Saved summaries and plots retain source readiness, original category labels
  and the absence of calibration intervals. Two-family values remain provisional.
- Cumulative threshold calculations now check local accuracy in ability units,
  including very small positive discriminations. Upper-tail cumulative plots
  retain the complementary 50% boundary lines and correct category labels.
- MML accuracy errors now keep the last optimizer candidate with its own
  integration settings, free-coordinate names and diagnostics, separately from
  the last attempted stage. Initial values are identified when optimization
  never returned. For errors with `provenance_version = 2`, `last_integration`
  and `last_quad_points` describe `last_parameters`; use `last_attempt` for the
  latest failed stage. Retained candidates have not passed numerical accuracy
  checks and are not fitted-result objects.
- Category-average ordering warnings now tolerate floating-point roundoff in
  equal means while retaining substantive decreases and missing-value cautions.

- PCM fair-average plot intervals now use the derivative of the same weighted
  administered-scale expectation as the reported fair average. They remain
  conditional delta-method intervals; uncertainty in other calibration
  parameters and the reference population is not included. The public
  `plot(..., type = "fair_average")` route accepts the table's score-selection
  and reporting-scale controls.
- ConQuest imports interpret factor-valued estimates by their labels and
  preserve nonnumeric case-estimate flags when reading native CSV exports.
  Recreate affected comparisons from the original tables or CSV files.
- Unexpected-response plots now share the table's residual-only cutoff of 3,
  preserve saved screening thresholds and show a valid empty result when all
  responses were evaluated and none were flagged.
- Model metadata, connectivity review and portable calibrations retain Unicode
  facet names in non-UTF-8 locales, including corrected JML. Result CSV exports use UTF-8; import them with
  `read.csv(..., encoding = "UTF-8", check.names = FALSE)`.


Unreleased release candidate. This version combines native many-facet RSM/PCM
conditional estimation, reusable calibration and new-Person scoring,
external-feature clustering and multivariate G/D-studies, with corrections to
uncertainty, subgroup comparisons and design planning.

* FACETS-style diagnostics, fit-measure tables and ZSTD fit-pathway plots
  now default to fourth-moment fit df. Refreshing table/plot companions
  preserves the supplied diagnostics' recorded WHEXACT choice. Both WHEXACT
  choices use the Wilson--Hilferty cube root (`FALSE` = YES, `TRUE` = NO);
  `"legacy_linear"` retains the previous TRUE calculation, and
  `fit_df_method = "engine"` retains the previous df. PTMEA uses the oriented
  measure sum, displacement follows the facet sign, and chi-square df counts
  only eligible estimates and SEs. PCM fair scores average expectations under
  administered scales; non-Person references exclude nonfinite Person estimates
  and report the reference counts. Fixed anchors retain their displayed value,
  and malformed anchor schemas always stop fitting.

* Bias expected scores and departures now refer to the baseline calibration,
  with adjusted residuals in separate columns. The summary tests all cells
  against zero with df equal to their count. Cell t df remains count minus 1:
  this matches the executed FACETS 4.5.0 example despite the online Table 13
  description of count minus 2. Unexpected-response screening now defaults
  to absolute standardized residual >= 3; the former probability OR rule is
  available explicitly. Tied modes no longer yield asymmetric rankings.
  Recompute diagnostics and bias analyses from saved fits before rebuilding
  reports; older calculation bundles are refused when paired with a fit.
  No main calibration refit is required. A matched binary JML comparison with
  FACETS 4.5.0 checks these report definitions within printed precision; this
  is not general cross-software or small-sample inferential equivalence.

* Category median boundaries now solve cumulative 50% probabilities and retain
  adjacent Andrich steps separately. Out-of-range half-score crossings stay
  unavailable. Recreate old category-structure reports. Category generating
  measures are distinguished from Person means and expected response scores;
  outcome-selected residuals no longer receive unit-baseline ZSTD or fit flags.
  Chi-square calculations are stable under a common location shift, and their
  summaries preserve unavailable decisions before display truncation. Weighted
  totals retain decimal precision and conditional SE fallbacks scale once.
  Bias-count reports separate logit shifts from baseline and remaining score
  residuals. Fixed-field imports validate headings and separate facet suffixes;
  ambiguous layouts are rejected. Coverage descriptions now distinguish native
  observation exports, graph components and partial tables from FACETS formats.

* G/D-study help now distinguishes executed GENOVA 3.1 negative-component
  conventions from mGENOVA 2.1. Raw ANOVA components are retained by mfrmr;
  matching software defaults requires matching the component treatment as
  well as the design. urGENOVA's unbalanced estimator and absence of D-study
  output are explicitly distinguished from native MINQUE0 and D-studies.

* Ordinary direct MML now offers `mml_integration = "auto"`. It checks each
  Person's log marginal likelihood, EAP and posterior SD, plus marginal-gradient
  stability and stationarity, against refined adaptive integration. Failed
  checks trigger warm-started adaptive refits; exhausted or failed attempts
  raise an error with the numerical history. The final rule consistently
  determines the likelihood, information and posterior calculations. Controls,
  per-Person checks and replay settings are retained. This is opt-in numerical
  accuracy checking for RSM, PCM and one-family GPCM, including supported
  population regressions; it does not change the default, certify integration
  error bounds or extend formal inference eligibility.

* Adaptive MML reuses likelihood and gradient calculations for Persons with
  identical responses, rating designs, observation weights and population
  covariates. Independent Persons retain their separate contributions to the
  marginal likelihood. This reduces repeated computation when patterns recur;
  distinct covariates or observation patterns can leave little work to reuse.
  Integration orders, derivatives and inference requirements are unchanged.

* ConQuest comparison bundles now reject unsupported anchors, orientation,
  constraints, shrinkage, population-intercept and category settings before
  export. Generated commands reset prior analysis settings and explicitly set
  item centering, error termination, binary codes and posterior scoring controls.
  Native CSV import validates coefficient indices, dimensions and nonnegative
  population variance. Item comparisons center both programs on the same
  matched finite items, eliminating artificial differences when external items
  are missing; fewer than two matched items leave centered metrics unavailable.
  `MfrmrCentered` records this comparison coordinate while `CenteredEstimate`
  retains the original bundle coordinate. Recompute older reviews from their
  saved bundle and external tables; no calibration refit is needed for that
  correction. Diagnostic explanations distinguish native row residuals from
  ConQuest posterior-averaged parameter fit. This does not add ConQuest fit,
  threshold or full-covariance replication to the comparison bundle.

* ConQuest bundles now use reversible person and covariate aliases; CSV import
  preserves leading zeros, long numeric IDs and literal `NA` identifiers.
  `conquest_control` sets external calibration and posterior scoring controls
  independently of native quadrature (31 external calibration nodes by default).
  Native export normalization checks available item labels and accepts the
  generated review text as `report_file` to record stopping/model metadata.
  Reviews flag extra IDs and failed stopping, and distinguish unknown stopping
  status from a met criterion. EAP errors and the population intercept now
  require a verified common origin or an explicit `conquest_location_shift`;
  raw differences remain available. `Facility` is no longer an automatic
  difficulty-column alias. Recreate reviews from saved exports and matching
  bundles, retaining old bundles for old-ID exports. Neither a met stopping
  criterion nor default controls certify numerical precision.

* RSM/PCM joint interaction models now accept `mml_engine = "em"` and
  `"hybrid"` with a fixed standard-normal population and fixed quadrature.
  Hybrid fitting uses EM initialization followed by direct optimization;
  direct remains the default. Ordinary EM checks auxiliary and marginal
  likelihood ascent before accepting numerical M steps. Joint-cell intervals
  use the same full information and readiness checks across these engines.
  EM alone does not establish stationarity, global optimality or adequate
  integration, and need not be faster. Existing weight/anchor restrictions
  on joint-cell intervals remain in force. A remaining direct-gradient stall
  can trigger one polish with the likelihood constant removed from relative
  stopping; likelihoods retain their original scale in all outputs. Saved
  facet-interval objects also support `confint()` extraction at their stored
  confidence level, retaining unavailable limits.

* Simulation design recommendations now connect to `summary()`, `plot()`,
  `as_ggplot()`, `plot_data()` and `mfrm_report()`. Reports retain targets,
  threshold sources, replication aggregates, failed/unavailable criteria and
  workload constraints. High facet separation reliability describes distinct
  locations, not agreement; `facets = "Person"` targets Person differentiation.
  The underlying ranking and defaults are unchanged. Main-effects D-study
  approximations also enter shared reports, retaining their residual assumptions
  and newly saved source components/counts without pretending to estimate omitted
  interactions or confidence intervals. Missing historical metadata stay explicit.
  Fit summaries and measurement reports add `numerical_review` with separate
  optimization, curvature, starts, integration and reproducible-refit actions;
  this guidance does not alter numerical or inference acceptance criteria.

* `fit_mfrm(..., mml_starts = c(0, 0.5, -0.5))` optionally reviews specified
  starting values for ordinary MML. Summary and report tables retain every
  attempt's likelihood, gradient, status, warnings and failures; long tables
  preserve the named initial and terminal coordinates. The default start is
  always included and a worse converged point never replaces a better
  unfinished result. Existing readiness checks still apply. This checks
  specified starts, not global optimality or confidence-interval validity;
  integration sensitivity is a separate review. Default fitting is unchanged.

* `mfrm_report()` now accepts `estimate_all_bias()` collections directly.
  Reports and appendix exports retain the full requested search, all stored
  signed cell estimates, unavailable screens, failed pairs, settings and
  calculation history. No recalibration or multiplicity adjustment is added.
  The Linking and DFF guide connects conditional exploration to explicit joint
  interaction models and compatible model comparisons, distinguishing
  same-data exploratory selection from prespecified inference.

  Bias-cell Infit/Outfit summaries now use one grouped aggregation per
  interaction instead of repeating the same setup for every cell. The cell
  likelihood, conditional uncertainty, search boundaries and complete-cell
  availability rules are unchanged; this reduces reporting overhead without
  lowering numerical precision.

  Marginal pairwise diagnostics reuse facet-label ordering and validate an
  exactly shared posterior once per context. Unequal posteriors retain the
  original pairwise tolerance checks, and unavailable observations stay visible.
  This reduces repeated preparation without changing integrated probabilities
  or screening rules.

* Pairwise bias reports now orient local measures using the target-facet sign
  and condition local SEs on the baseline calibration. The default
  `df_method = "facets"` uses modified Rasch-Welch df rounded down for the tail
  area, matching an unequal-count FACETS 4.5.0 example. `DFUnrounded` retains
  the continuous value; `df_method = "classical"` reproduces the earlier Welch
  reference. Weighted counts govern either rule, with missing inputs retained
  as unavailable. Recreate older reports to obtain `CalculationVersion = 3`;
  first recreate old bias bundles from saved fits and current diagnostics.
  Summaries cover the full family before display truncation. These conditional
  screens do not propagate calibration uncertainty or establish fairness.

* `load_mfrmr_data()` now uses R's lazy-loaded package data directly, avoiding
  repeated dataset discovery and deserialization. It returns the same tables;
  editing a returned table does not alter the packaged data. Short help examples
  now run in ordinary checks without redundant `donttest` wrappers.

* Linking reviews now retain low-support warnings when no drift comparisons
  are flagged. An empty drift table previously caused a column-type error
  while combining these warnings.
  Diagnostics and manifest export also recognize the legacy `JMLE` label in
  saved fits without changing the original object or its joint estimates.
  Diagnostic summaries wrap overall-fit explanations below the numeric table,
  preserving the full inference caution at ordinary console widths.
  Supplying another model's diagnostics now consistently reports the
  fit/diagnostics mismatch, including models with different optional constraint
  entries; the source-identity checks remain unchanged.
  Point-only testlet calibration tables explain missing standard errors
  without an incidental missing-column warning or fabricated intervals.

* Signal-detection simulations now record seeds for complete trials, retain
  generation failures with their stage/reason, and orient DIF contrasts as
  reference minus focal. Refit contrast means require documented common-scale
  linking. Summaries expose point, classification and non-target denominators,
  and correct the Monte Carlo SE of non-target rates (previously always missing).
  Missing adjusted screening statistics no longer fall back to raw values.
  Save new results separately: stochastic fitting/diagnostics now use the trial
  stream, and old results lacking direction/linking metadata cannot reconstruct
  the corresponding contrast means. Raw saved results remain available.
  DIF and bias outputs retain their existing descriptive/screening limits.
  Subgroup refits also preserve an explicit repeated-rating independence
  declaration from the baseline fit; they do not infer independence for legacy
  or undeclared repetitions.

* Population and assignment templates now support joint Person-level sampling
  through `population_assignment = "matched"` in simulation-spec construction
  and extraction. This default preserves covariate/group/assignment associations;
  use `"independent"` for a deliberately independent scenario. Template IDs and
  shared group labels are checked, and source template IDs are retained in truth
  metadata. A population `Group` now also identifies response groups instead of
  allowing a second unrelated group draw. Old specifications without the new
  setting retain independent sampling; conflicting definitions of `Group` error.
  Preserve old generated data and use a new study file for the corrected scenario.

* Observed simulation skeletons retain repeated rating rows and their weights.
  Each row generates a separate conditionally independent response, and study
  runners declare this known generation assumption when refitting. This does
  not model dependence among repeated source ratings or turn weights into
  additional generated observations. Remove accidental source duplicates first.

* Latent-regression simulation now preserves the template/fitted basis for
  `scale()`, polynomials and splines, together with categorical coding, through
  generation and refitting. Previously each sample could use a different basis
  with the same coefficients. Rebuild old transformed specifications from their
  original template/fit and use a new study file; retain previous results with
  their source environment. Fit-derived GPCM specifications now align slopes
  to their actual owner facet, including custom names, and refresh planning
  metadata. Forecast help and returned notes distinguish Monte Carlo precision
  from uncertainty in the source fit or generating assumptions.

* Normal simulation generators no longer sort rater/criterion effects before
  assigning them to level IDs. Each ID now receives an independent normal draw,
  avoiding an unintended association with ordered IDs and ID-specific steps or
  slopes. This changes seeded generated datasets. Preserve earlier data/results
  and their software environment; new runs use the corrected generator.

* Observed-data resampling now keeps distinct Person and stratum identities,
  refuses input columns that would be overwritten by metadata, and avoids
  redundant coverage top-ups. Stratum manifests distinguish base, added and
  final cluster counts. `sampling_scheme` and saved interpretation distinguish
  same-size Person-cluster bootstrap draws from reduced-size and coverage-augmented
  stability designs. Person resampling conditions on the observed non-Person
  facets; it does not resample new raters/tasks or qualify interval coverage.
  Old specifications are revalidated before drawing; corrected labels/top-up
  can change draws under the same seed. Preserve old results separately.

* Design evaluation now applies the identified PCM/GPCM step-location truth
  mapping used by recovery studies, and retains representable large-error RMSEs.
  New summaries withhold legacy PCM/GPCM recovery metrics whose coordinates
  cannot be reconstructed; other diagnostics remain. Rebuild old summaries
  before recommending a design, and use a new study/file for corrected recovery.
  Unavailable recommendation extrema return `NA` instead of misleading infinities.

* Recovery summaries now calculate Monte Carlo SEs over independent trials,
  including unequal available-target counts, instead of treating parameter rows
  within a trial as independent. Centered structural facet coverage uses the
  matching covariance through `mfrm_facet_intervals()`; source uncertainty is
  retained separately as `RawSE`. JML and Person-score point recovery remain
  available, with explicit reasons for unavailable structural coverage.
  Summaries, assessments and plot data retain expected/available counts and
  distinguish conditional coverage from interval availability. Mean-aligned
  bias and complete identified log-slope group bias are no longer treated as
  adequacy evidence when their average is zero by construction. Parameter-specific
  assessment thresholds now retain their names and apply to the intended groups.
  PCM/GPCM facet truth now absorbs generating step means; fitted population
  intercept truth follows the same identification. `summary()`, assessment and
  summary/coverage plotting refresh old aggregates without refitting, but withhold
  legacy intervals whose covariance is missing. Start a new study/file for the
  corrected truth mapping and interval evidence; see `simulate_mfrm_recovery()`.

* `simulate_mfrm_recovery()` and `evaluate_mfrm_design()` now accept `file`
  and `file_refit` for checked result reuse and interrupted-study resumption.
  Saved trials retain failures and preallocated seeds. Reuse checks inputs,
  fitting controls, requested replications, RNG kind and implementation/software
  identity; incompatible files error by default. Serial checkpoints retain each
  completed trial; future design evaluations retain each completed design row.
  Completed RDS files contain ordinary result objects with reproducibility
  metadata. Recovery trials now control the RNG across fitting and diagnostics
  as well as data generation. See the function help for replacement policies
  and the limits of reproducibility across systems.

* The output guide now connects exploratory all-pair bias analysis to its
  summary, plots and saved availability information. G-theory guidance starts
  from components, plan coefficients, measurement SEM and error-source plots,
  with interval details in the method help. GPCM capability guidance distinguishes
  external two-family Stan code export from unavailable native MCMC estimation.
  These changes clarify existing routes without changing calculations or
  statistical qualifications.

* mirt imports now handle mixtures of binary and partial-credit items, including
  item-specific category counts and rating-scale offsets. Parameter extraction
  distinguishes inapplicable coefficient-table cells from missing parameters.
  Existing source slopes, ability coordinates and uncertainty conventions are
  preserved. Re-import the source mirt fit; no model re-estimation is needed.

* Facet-equivalence ROPE summaries now preserve symmetric tail probabilities
  and very narrow normal-reference bands. Unresolvable display intervals retain
  estimates and SEs with `CIReason`, `DeviationCIReason` or `CI90Reason`;
  mean-deviation status is unavailable when its interval cannot be represented.
  The fixed 5% TOST rule is unchanged and remains separate from display intervals.
  Recompute equivalence summaries from saved eligible fits before exporting them.

* Model-comparison evidence ratios now use AIC differences directly, retaining
  finite ratios between candidates whose normalized Akaike weights both round
  to zero. `LogEvidenceRatio` and `RatioReason` preserve the log ratio and explain
  numeric overflow/underflow. Recreate comparisons from saved fits; estimation,
  criterion values, candidate weights and model-selection restrictions are unchanged.

* External-feature group profiles now preserve representable sample SDs when
  squaring very large or very small values would overflow or underflow.
  Group memberships and distances are unchanged. Recreate affected clustering
  results from the original features or retained imputations to update profiles;
  no MFRM fit or imputation model needs to be refitted. Help clarifies original-unit
  profile SDs, singleton silhouettes and the equal-spacing convention for ordered
  factors in Gower dissimilarities.

* Multiple-imputation pooling now preserves extreme confidence tails and
  small positive within-imputation fractions in the finite-data degrees of
  freedom. `IntervalReason` explains numerical interval limits while keeping
  point estimates and variances; saved notes state the independent-completion
  assumption. Re-pool existing valid analyses to update tables and plots.
  Input review rejects colliding numeric IDs and detects exact changes to
  observed numeric predictors. Pooling also revalidates the original roster
  and its saved review; affected inputs must be corrected before analysis.

* Model-expected inter-rater agreement no longer uses finite optimizer proxies
  when the corresponding primary estimates are unavailable. Observed
  agreement is retained; expectation counts and `ExpectedReason` explain
  missing values, including repeated cells and unavailable probabilities.
  Help and reports clarify the conditional point-score calculation, including
  MML EAP plug-in probabilities. Recompute diagnostics from saved fits and
  recreate agreement reports/exports; calibration refitting is unnecessary.
  `fit_rater_accuracy()` also rejects distinct numeric identifiers that would
  collapse to the same text label. Supply unique character IDs before fitting.

* Halo networks now distinguish shared contexts from independent Persons.
  Repeated Persons across tasks, or context definitions omitting Person,
  retain descriptive correlations but withhold ordinary correlation p-values,
  with explicit counts and reasons. `alpha = 1` disables the p-value filter
  as documented, including for finite correlations without a test. Other
  filters still apply. Recreate `rater_halo_network_analysis()` results from
  saved fits and matching diagnostics, then regenerate summaries, figures
  and exports; no model refit is needed. Review flags remain descriptive
  screens, without a validated halo-detection interpretation.

* Mantel--Haenszel calculations now avoid integer overflow in large strata.
  The defect could remove an asymptotic interval or make an otherwise valid
  odds ratio and test unavailable. Asymptotic inference uses the same
  count-moment, log-odds variance and continuity-correction formulas with
  double-precision products; exact conditional inference is unchanged.
  Recompute affected `mantel_haenszel_table()` results from the original
  response data, then regenerate exports. No MFRM refit is involved.

* Linking offsets and paired SEs now avoid overflow/underflow in squared SEs,
  inverse weights and weighted sums. Drift/chain results expose unavailable
  calculations, fallback screening and actual contributing elements.
  Singleton residual spread is missing, not zero; unavailable flagged counts
  and cumulative offsets retain their status and reason. Recompute
  `detect_anchor_drift()` / `build_equating_chain()` from saved fitted waves;
  no model refit is needed. Regenerate plots and exports. ggplot conversion
  now preserves drift, per-wave estimates and stored intervals, cumulative
  offsets and common-element counts instead of selecting the first numeric
  column. Unavailable chain links are not bridged, and all-zero/unavailable
  drift views render without inventing a result. These remain conditional
  review screens, without full linking uncertainty or invariance tests.

* Refit DFF no longer substitutes baseline estimates when a subgroup fit
  fails or is skipped. Absent target levels remain as unavailable rows with
  reasons. Sparsity and Welch-screen counts now use each target-by-group
  cell's retained responses; subgroup totals and Person counts are separate.
  Missing SEs or shared Persons withhold the independent-subgroup screen
  while retaining actual point estimates. Direction labels describe location
  differences, including exact equality, without implying an ordering of
  response curves with different steps or slopes. Reports and plots retain
  the availability information. Rerun `analyze_dff(method = "refit")` from
  saved baseline fits and original data, then regenerate tables and figures;
  subgroup fitting is required. Earlier refit objects are rejected by
  summary/report/plot methods. Formal DFF inference remains unavailable.

* Bias screens now retain search-boundary point estimates with explicit
  statuses while withholding their SEs and screening tail areas. Conditional
  information includes every row's positive variance; Infit/Outfit and the
  heterogeneity screen no longer silently use an available subset. Cell-level
  counts and reasons explain missing results. Summaries distinguish missing
  screens from zero flags and retain unavailable cells in Bonferroni/Holm
  family sizes, without claiming formal multiplicity control. Facet labels
  containing `||` no longer merge distinct bias cells. Recompute
  `estimate_bias()` / `estimate_all_bias()` from saved fits and matching
  diagnostics; calibration refitting is unnecessary.

* Explicit two-crossed-facet G score intervals now distinguish a numerically
  unresolved test cutoff from an empty confidence set. Very small confidence
  levels previously could report a completed empty set or inadequately
  checked endpoints. Cutoff checks and endpoint searches now scale their
  accuracy to the critical value. Points are retained, with missing bounds,
  explanatory statuses and saved precision diagnostics when unresolved.
  Recompute `confint(d, method = "score", assumption = "normal")` from saved
  D-study summaries; no model refit is needed. The method retains its
  experimental calibration and finite-grid limitations.

* Explicit two-facet approximate-F and mean-square-region intervals now
  distinguish numerically unresolved quantiles/regions from zero uncertainty.
  SEM intervals also report limits lost through square-root or unit
  transformations. Missing limits retain points and explanatory statuses;
  diagnostics save the relevant quantiles or normalized region bounds.
  Recompute `confint(d, ...)` from saved D-study summaries; no model refit
  is needed. The approximate-F method retains its finite-sample limitations,
  and the region method retains its conservative normal-model scope.
  Combined base R interval plots now size their left margin to the plan
  labels, including numbering used to distinguish duplicate plans.

* One-facet spectral F intervals now reject numerically unresolved F quantiles
  and ratio bounds, retaining the original point and an explicit reason.
  Shared confidence-set mappings for G/Phi plan differences also withhold
  intervals whose widths or nonzero endpoints disappear through rounding,
  and use more accurate extrema for nearly equal future counts. Identical
  counts still give exactly zero when the base set is available. Recompute
  `confint(d, ...)` from saved D-study results to apply these checks; spectral
  F uses retained ratings without model fitting. Other methods retain their
  own computational and coverage qualifications.

* REML profiles now preserve a composite when score units and weights are
  changed by reciprocal factors. Previously, extreme differences in units
  could incorrectly make the composite appear constant or erase its component
  estimate. Profile endpoints also use a cutoff-scaled accuracy check:
  numerically unresolved cutoffs and intervals lost when restoring score units
  return missing bounds with reasons, retaining points and representable REML
  centers. Exact identical-plan direct contrasts remain zero. Recompute
  `confint(..., method = "reml_profile", assumption = "normal")` from a saved
  G/D-study retaining its used ratings; the separate REML model is fitted again.
  These repairs do not change the experimental coverage qualifications.

* Between-score covariance intervals now scale the two scores separately.
  This fixes SEs that became zero, or intervals that disappeared, when the
  scores used very different units. SEs remain calculable when only their
  squares fall outside floating-point range. Collapsed covariance intervals,
  including residual-contrast intervals, return missing limits with reasons
  while retaining point estimates. Generalized intervals distinguish
  unresolved quantile probabilities and empirical ranks. Recompute
  `confint(g, ...)` from a saved G-study; no model refit is needed.
  The returned `mean_products$score_scale` records each score's divisor;
  `scale` restores covariance units, not both marginal variance columns.

* Variance-component and measurement-SEM intervals now share the scalar
  mean-square calculation and reject positive intervals whose width rounds
  away, including residual-contrast intervals and SEM square-root images.
  Points and calculable SEs remain available with explicit interval statuses.
  This prevents extremely small confidence levels from producing misleading
  zero-width intervals. Recompute `confint()` from saved G/D-study results;
  no source-model refit is needed.

* Two-crossed-facet MLS G/Phi intervals now resolve narrow sets near one
  using shifted and scaled polynomial coordinates. This corrects endpoint
  errors previously reported as completed calculations. Root separation
  outside numerical precision returns missing bounds and an explicit status,
  while retaining the point estimate. Recompute `confint(d, ...)` from a saved
  D-study to apply the repair; no G-study refit is needed.

* Main-effects G/D studies now distinguish a zero variance estimate from
  components that the observed design cannot separate. Fitted allocations
  carry identification flags; identified error sums remain usable, while
  projections depending on arbitrary allocations return missing values with
  reasons. Zero object variance now gives G/Phi of zero when the corresponding
  error is positive, and no longer hides calculable error variances.
  Scaled coefficient sums and successive residual division avoid overflow;
  quantities outside numerical range have explicit reasons. Recreate an
  older G-study from its ratings for the component-separation check.
  Saved components can be reused for the numerical projection repairs.
  Selecting a non-Person object facet now retains Person among the default
  random measurement conditions, instead of silently dropping it.

* `plot(d, type = "error_components")` now explains multivariate D-study
  error variance by source for each future plan. Crossed and nested counts
  use their corresponding divisors; composite weights include cross-score
  covariance. Negative contributions remain visible, with original metric
  availability and totals retained in `plot_data()`. Base and ggplot views
  use saved components without refitting. Contributions are point projections
  in squared score units, not SEM, confidence intervals or an optimal plan.

* `estimate_all_bias()` now supports `summary()`, concise printing and `plot()`.
  The summary combines retained cell tables with their pair identifiers,
  observed counts and availability; the overview plot shows signed shifts
  across pairs, and `pair` selects an individual heatmap or ranked view.
  Failed pairs and cells without usable screening statistics are distinguished
  from evaluated cells with no flags. Repeated facet combinations are rejected.
  Saved collections can be summarized and plotted without refitting.
  These remain conditional exploratory screens, without collection-wide
  multiplicity control.

* One-family GPCM MML location and contrast intervals are now described as
  supported pointwise Wald approximations under explicit regularity conditions,
  consistently with full joint MML information and delta-method uncertainty.
  This also clarifies the estimated-normal RSM/PCM native-location route.
  Estimates, numerical eligibility checks and interval calculations are unchanged;
  finite-sample coverage is not guaranteed. Corrected JML, optional slope
  profiles and two-family extensions retain their separate qualifications.

* D-study G intervals now withhold limits when the base F-ratio width is
  lost to numerical rounding, including intervals for plan differences.
  Delta G/Phi intervals also withhold limits whose positive width rounds
  away, while retaining the coefficient and its SE. This can occur at an
  extremely small requested confidence level. Returned statuses explain
  the missing limits. Recompute `confint(d, ...)` from a saved D-study to
  apply the checks; no G-study refit is needed.

* Shared-rater population-SD profiles now verify that endpoints attain the
  likelihood-ratio cutoff, preventing spuriously coincident limits for narrow
  profiles. Known-ability-SD profiles also check nuisance curvature at both
  integration orders, so stationary nuisance maxima cannot pass as optima.
  Numerical failures retain candidate likelihood ratios, information checks
  and explanations in the error's `profile` field. Recompute an older interval
  with `confint(fit, parm = "rater_sd")` to apply these checks; saved intervals
  do not change automatically. These repairs do not establish coverage.

* MML uncertainty and inferential model comparison now reject known marginal
  response-dimension deficits and complete-pattern first-order rank failures.
  A numerically positive Hessian can no longer override these findings in GPCM
  slope/curve intervals or local model-comparison checks. Descriptive points
  and unavailable-result reasons are preserved. The dimension check distinguishes
  rating rosters and population covariates; it is not a minimum Person count.

* Full RSM/PCM calibration can now accept one Person for JML, CML and
  fixed-population MML, subject to identification, category and numerical checks.
  One-Person JML/MML uses `lpSolve` to check additive recession before fitting;
  a detected or unresolved boundary returns no finite estimate. Data/anchor
  review also accepts one Person. Results retain the single-Person caution;
  available model-information intervals do not establish finite-sample coverage
  or between-Person variation, and Person-cluster sandwich inference remains
  unavailable. Estimated-population MML, GPCM, corrected JML and comparison-only
  calibration retain their at-least-two-Person support limit.
  The Wright map uses a visible display bin for a lone Person estimate near
  zero and retains a note distinguishing its width from uncertainty.

* Added a central model-and-estimator choice guide, linked from the README,
  first-analysis tutorial and help guides. It starts from assessment questions,
  compares JML, corrected JML, MML and CML and the available model families,
  and explains the benefits and costs of additional complexity, population
  scale choices and valid model comparisons. The GPCM explanation distinguishes
  ordinary total-score conditioning from other possible sufficient statistics,
  and shows why two-family scale constraints and connected crossings alone
  do not establish marginal identification.

* For eligible ordinary fixed-population RSM/PCM MML fits, failed joint
  information no longer stops facet-interval reporting. Finite points remain
  with missing model/sandwich uncertainty and saved reasons. `vcov()` returns
  a named `NA` matrix when its storage fits the information budget; default
  `confint()` retains points and reasons in diagnostics. Constraint-fixed
  targets remain distinct. Source qualification and cluster-input checks
  still apply; no regularized inverse or new coverage guarantee is supplied.

* CML now retains finite numerical point estimates when identification,
  finite-MLE and stationarity checks pass but conditional information is
  unavailable. Summaries, facet intervals and default `confint()` retain
  missing uncertainty with reasons; `vcov()` returns an `NA` matrix with
  diagnostic attributes. Comparison-only CML similarly separates retained
  points from unavailable intervals. Information warnings still limit point
  accuracy in weak directions; identification, boundary and convergence
  failures do not gain estimates or inference through this change.

* Removed the blanket ten-response-row minimum. Small supported datasets now
  reach the estimator's model-specific checks. RSM/PCM model-information
  intervals can also retain singleton category counts after fresh category,
  convergence, full-information and finer-quadrature checks; empty categories
  and identification failures remain restricted. The source fit retains its
  category review. Interval tables, `vcov()` and `confint()` preserve the
  small-sample caution; numerical availability does not guarantee coverage.

* **Breaking API names:** `simulate_mfrm_recovery()` replaces
  `evaluate_mfrm_recovery()` for repeated data generation and refitting;
  `assess_mfrm_recovery()` continues to review saved simulation results and
  now rejects unused options. `review_peer_design()` replaces
  `build_peer_review_design_review()` for existing peer-assessment assignment
  metadata. These renames preserve calculations and saved-result classes.
  GPCM availability guidance now uses `gpcm_capability_matrix()` throughout:
  `view = "scores"` replaces `gpcm_score_side_contract()`, and
  `view = "unavailable"` replaces `gpcm_runtime_guard_coverage()`. Each view
  retains its status meanings and columns; the compatibility guide maps calls.

* **Breaking CSV export rename:** use `export_mfrm_tables()` instead of
  `export_mfrm()`. Arguments, CSV filenames/values and the returned file table
  are unchanged. Enriched facet CSVs now keep SE methods, inference/fit
  availability and reasons, in the source facet-table order.
  `export_mfrm_results()` remains the results/RDS archive route;
  `export_mfrm_bundle()` remains the fit-level route with an optional refitting
  script. The reports guide and compatibility guide explain these choices.
  ZIP replacement now contains only the current export's files, checks the
  command outcome and archive contents, and retains an older ZIP if creation
  fails without reporting it as newly written. Bundle file lists no longer
  pick up preexisting tables that were not written by the current call.
  Results manifests and file counts now include plot-error CSVs.

* **Breaking G/D-study API change:** use `mfrm_g_study()` for one or several
  scores and `mfrm_d_study()` for plans. These replace the public
  `mfrm_multivariate_gstudy()` and `mfrm_multivariate_d_study()` names.
  D-study dispatch preserves the source model, estimator and saved-result
  classes; existing RDS objects remain usable. Both routes accept facet labels
  or `n_<facet>` plan columns, retaining their existing output-column names.
  Main-effects D-study `summary()` now returns its coefficient table.
  Previously ignored extra arguments, including unsupported composite weights,
  are rejected. See the compatibility guide for migration and model-specific limits.
  D-study curve legends now recalculate text spacing when recorded plots are
  redrawn at a different size, preventing overlapping plan labels in saved
  figures and tutorials.

* **Breaking change:** `confint(fit)` now targets the free coefficients from
  `coef(fit)` for supported native MML/CML fits. Coefficient names or positions
  select subsets, with covariance qualification and availability reasons
  retained. CML comparison fits preserve unidentified targets as missing rows.
  For the former GPCM slope interface, explicitly use
  `confint(fit, parm = "slopes", ...)`, including for standardized slopes,
  ratios/differences, sandwich/profile intervals and slope plots/reports.
  Saved slope-interval objects and other result classes' defaults are unchanged.
  See `?confint.mfrm_fit` and the compatibility guide.

* MML response diagnostics now match saved Person estimates by identifier,
  correcting misassigned expected scores/probabilities after reordering the
  Person table. Legacy diagnostic rows and score/residual exports explain
  JML optimizer coordinates versus primary estimates, and distinguish a
  positive variance below the diagnostic cutoff from a zero variance.
  Overall/element mean squares no longer count unavailable residuals as zero
  misfit: incomplete groups retain `NA`, row coverage and a reason. Regenerate
  saved diagnostics and exports from the existing fit to obtain the corrections.

* Ordinary native fits now provide `fitted()`, `residuals()` and source-rating
  `predict()` methods. These use conditional probabilities at saved primary
  Person estimates, with MML EAP plug-in values distinguished from posterior
  mean predictions. Extreme JML Persons use limiting probabilities rather than
  finite optimizer proxies; zero-variance Pearson residuals are unavailable
  with a reason. Row identities, modeled-category mappings and interpretation
  limits accompany the values. New-data ability scoring retains its dedicated
  interfaces. See `?mfrm_fit_responses` and the compatibility guide.

* Native fits now provide standard `coef()`, `vcov()`, `logLik()` and `nobs()`
  methods with explicit estimator and coordinate meanings. Default coefficients
  align with free-coordinate covariance; `coef(type = "locations")` and other
  component choices retain labelled primary estimates. MML BIC uses Persons,
  while `nobs(type = "ratings")` counts retained response rows. Unsupported
  information criteria remain `NA`, and multiple-model AIC/BIC calls direct
  users to `compare_mfrm()` so likelihood bases and source compatibility are
  checked. See `?mfrm_fit_extractors` and the compatibility guide for migration.

* Individual rater sheets now retain a recipient-safe interval `Reason` in
  their uncertainty table and explain failed numerical or profile checks.
  Unavailable intervals omit both endpoints, including partially computed
  profiles; the analyst's original result is unchanged. Printed sheets also
  show uncertainty and its availability. Recreate sheets from saved results
  to obtain the new explanations; existing saved sheets are not rewritten.

* Experimental G/Phi REML profiles now explain when their coefficient center
  rounds to one and prevents interval calculation. The source optimizer's
  numerical success is retained separately from the unavailable interval;
  a rounded center is not reported as perfect reliability or zero uncertainty.

* Testlet and shared-rater normal intervals now retain row-specific availability
  reasons and source checks in `attr(ci, "diagnostics")` and `attr(ci, "checks")`.
  Calibration/summary tables distinguish unrequested intervals, unsupported
  targets, variance boundaries and numerical failures. Facet intervals retain
  plain-language explanations alongside status codes, including insufficient
  independent-cluster rank. Person scoring now explains prior-only rows in
  `Reason`. The interval guide explains where to find and how to preserve these
  qualifications. Recreate summaries/interval displays from saved fits to obtain
  the new fields; reading an older saved table does not update it.

* D-study exact-F G intervals now withhold endpoints when a positive interval
  width rounds away on the coefficient scale. The retained log-ratio interval
  can still support representable plan differences. Phi generalized intervals
  also withhold rounded-away widths and refuse difference intervals when their
  base pivotal interval has lost numeric resolution. A failed base no longer
  yields a spurious zero-width difference interval. Raw pivot diagnostics remain
  available; recreate affected intervals from the saved D-study, without refitting.

* `classification_table()` now directly integrates narrow normal-probability
  bands when CDF subtraction would lose precision, retaining the band width
  before centering at the supplied ability. Recompute affected tables from the
  same calibration and plan. This improves numerical evaluation of the Rudner
  approximation; it does not change its statistical assumptions or qualifications.

* D-study count columns now resolve one-to-one to random facets, including
  facet names starting with `n_`. The simplified `mfrm_d_study()` now refuses
  duplicate, missing, or extra count columns instead of silently reusing or
  ignoring them. Supply exactly one raw or `n_<facet>` count column per facet.
  Literal facet names, including spaces, are preserved in default/list plans
  and returned count columns.
  Recompute affected projections from their saved G-study; no source refit
  is needed. Multivariate crossed/nested projections and connected interval
  calculations also avoid integer overflow when multiplying large planned
  counts. Integer and double storage now give the same results. This does not
  change the source information or establish interval coverage for large plans.

* The simplified G/D-study APIs now retain full-precision G/Phi values,
  consistent with the values used for their reference-status labels. Previously,
  rounding the stored coefficients to four decimals could move them across a
  reference threshold or erase small plan differences. Round only for display.
  Recreate G-study summaries from the existing MFRM fit, or D-study projections
  from the saved G-study components, to recover full precision.

* Optimizer polishing no longer selects a materially worse likelihood merely
  because a candidate has a smaller gradient or passes its convergence check.
  The acceptance comparison now allows only floating-point rounding error.
  Incomplete or invalid MML solution-information summaries remain ineligible
  with a reason instead of raising an internal length-zero error.

* Two-family GPCM/GMFRM EM now preserves quadrature weights and normalized
  expected counts at extreme finite likelihoods. Its log-slope derivatives
  and the ordinary GPCM JML/MML derivatives also avoid cancellation from large,
  nearly tied category logits. GPCM
  profile intervals now check positive nuisance curvature at both integration
  orders; agreement between starts and small gradients alone cannot qualify a
  stationary maximum or saddle of the negative log likelihood. Failed checks
  retain missing bounds and reasons. Refit affected calibrations and recompute
  dependent intervals; saved profile intervals require recomputation to apply
  the new curvature checks.

* JML/MML probability and likelihood calculations now normalize after removing
  shared offsets, preventing cancellation when finite category logits
  are large and nearly tied. R and C++ likelihoods, adaptive integration and
  posterior probabilities use the same stable ordering. Optimizer code zero
  no longer passes the numerical-readiness check when the terminal gradient
  is missing, incomplete or contains nonfinite values. Affected saved fits
  require refitting/review; changing their readiness flags is not a repair.

* New-Person and portable-calibration scoring now share the stable posterior
  normalization used by model fitting. Large common log-likelihood offsets
  no longer distort prior weights, EAPs or posterior SDs. Continuous posterior
  intervals and cut-score probabilities also subtract the likelihood offset
  before adding the prior density. Recompute affected saved scores and draws
  from their eligible calibration; scoring-only corrections do not require a
  model refit. These changes do not add calibration uncertainty to conditional
  posterior intervals.

* Ordinary source-Person, testlet and shared-rater scoring also preserve the
  normal prior when centering large log likelihoods. Testlet nested weights,
  gradients and posterior SDs, shared-rater likelihoods/derivatives, and the
  connected response diagnostics use stable normalization. A testlet derivative
  calculation with one latent point now retains its matrix dimensions.
  Unresolved response-normalization checks remain unavailable. Recompute
  affected scores/diagnostics; refit first if the calibration likelihood was
  affected. Conditional intervals still exclude calibration uncertainty and
  shared-rater scoring retains its Laplace-approximation limits.

* Shared model preparation now rounds accepted near-integer scores to the
  nearest category instead of truncating toward zero, and preserves the
  supplied precision of numeric weights. Grouped/rowwise tibbles are processed
  without inherited dplyr grouping. Reserved facet names (`Person`, `Score`,
  `Weight`, `score_k`) are rejected before they can collide with internal
  columns; integer category and weight-total overflows also fail early.
  Saved fits affected by incorrect category coding or column collisions must
  be refitted from the original corrected data. See `?mfrmr_compatibility_layer`.

* Parameter bookkeeping also rejects facets named `theta`, `steps`,
  `interactions`, `log_slopes`, `beta` or `log_sigma2`, with renaming/refitting
  guidance. Similar names such as `log_sigma2_group` remain supported: exact
  lookup prevents them from being mistaken for absent population or slope
  blocks, including in adaptive MML gradients and parameter maps. Incorrect
  free-vector lengths now fail instead of silently dropping trailing values.

* Diagnostic consumers now check that supplied diagnostics belong to the fitted
  analysis, including tables, plots, bias screens, reviews and file exports.
  Matching labels and location estimates alone no longer suffice: the stored
  calibration, thresholds, population specification and observation identity
  must also match. Exports reject mismatches before writing files. Older native
  diagnostics without calibration metadata must be regenerated from their saved
  fit; model refitting is unnecessary. See `?mfrmr_compatibility_layer` for the
  0.2.4 migration instructions.

* Supplied bias results, unexpected-response tables and displacement tables now
  retain their fitted-analysis identity. Reports, dashboards, QC, casebooks and
  archive/replay helpers reject results from another fit, including mismatched
  members of a bias collection. Older results without this identity must be
  recreated from the saved fit and matching diagnostics using the original
  screening settings; model refitting is unnecessary. Standalone summaries and
  plots retain their existing input rules.

* Prediction exports and replay now verify the source fit of unit scores,
  plausible values and fit-derived design forecasts, including fitted
  population coding and scoring status. Mismatches are rejected before archive
  files are created. Older unit scores/draws can be regenerated from their saved
  calibration; older fit-derived forecasts require rerunning their simulations.
  Explicit simulation-specification forecasts remain separate scenarios, with
  their source recorded in the manifest. See `?mfrmr_compatibility_layer`.

* `detect_facet_nesting()` now calculates entropy directly from weighted
  frequencies, preserving fractional weights and avoiding weight-sized row
  replication. It uses the same complete, positive-weight rows for all
  entropies in a pair, correcting negative indices caused by inconsistent
  missing-value handling. `UsedRows`/`ExcludedRows` identify each pair's scope.
  Invalid or missing weight columns now raise errors; zero-weight rows are
  excluded. Recalculate saved weighted or incomplete-data nesting reviews.

* Structural reviews supplied to `reporting_checklist()` now verify the
  reviewed columns and retained rating rows, including scores when ICC was
  computed. Matching raw-data reviews and subsets of fitted facets remain
  supported, with their scope named explicitly. Older reviews must be
  regenerated with their original review settings; no MFRM refit is needed.
  Standalone reviews remain available. See `?mfrmr_compatibility_layer`.

* `facet_small_sample_review()` now reports the lowest count band as
  `WorstCategory`, correcting an aggregation that selected the highest band.
  It also validates custom thresholds. Its help, saved notes, plot and generated
  manuscript text identify counts as retained rating rows and avoid unsupported
  claims about precision, interval coverage, targeting or publication adequacy.
  The reporting-checklist item is now named `Facet rating-count screening`;
  update scripts that select its former `Facet sample-size adequacy` label.

* The G-study sampling-covariance help now combines a directly computed balanced
  example with saved incomplete-design results and their complete regeneration
  recipe. Both variance and between-score covariance targets remain illustrated;
  the estimator, assumptions and numerical checks are unchanged.

* Help now introduces estimator and output choices before specialist fitting
  details, uses a common results/report workflow, and supplies complete GPCM
  interval/curve and native/standardized facet-contrast examples, including the
  person-level data required for the fitted population SD. CML reporting
  instructions use its own tables, plots and saved interval objects; sending a
  CML fit to `mfrm_results()` now explains that supported alternative directly.

* Normal-reference intervals now evaluate upper-tail probabilities directly,
  preserving finite critical values for valid confidence levels extremely close
  to one, including Bonferroni-adjusted GPCM slopes and curve grids. The same
  correction applies to normal score/rater intervals and their tables and
  plots. It prevents spurious missing, infinite or full-support limits caused
  by rounding the complementary probability to one; inference qualifications
  and genuinely unrepresentable bounds retain their existing meanings.

* CML Wald intervals now retain finite normal critical values for valid confidence
  levels extremely close to one. Non-fixed contrasts whose variance underflows
  to zero are rejected with rescaling guidance instead of returning zero-width
  intervals. Confidence levels follow the same scalar-input contract across
  full and comparison-only CML interval routes.

* Output, interval and plot guides now distinguish CML structural/comparison
  intervals from full MML interaction intervals and population-SD standardization.
  CML guidance uses dedicated plots, tables and RDS saving; ordinary Person
  diagnostics, Wright maps and consolidated reports remain unavailable.

* Prospective interaction reviews now report `IdentificationBasis` and
  `IdentificationDetail` for each contrast. They distinguish unrestricted-Person
  response-model identification, dependence on declared population restrictions,
  and CML conditional support, while preserving unidentified, unevaluated and
  unsupported states. MML conditional results remain separate from global
  marginal identification. Fitted comparisons and interval tables retain this
  interpretation alongside numerical and uncertainty qualifications.

* `mfrm_facet_intervals(comparison_review = ...)` now returns model-information
  intervals for declared interaction contrasts from matching full RSM/PCM MML
  fits. `scale = "standardized"` uses the same fit's estimated intercept-only
  normal population SD and full joint delta-method covariance. Native intervals
  also support origin-specific population means; regression residual SD is not
  relabeled as marginal population SD. Saved plots/reports retain target weights,
  group order, scale and numerical qualifications. Global marginal identification
  and finite-sample coverage remain separate from numerical interval availability.

* Native RSM/PCM CML now supports `fit_scope = "comparisons"` with a CML
  comparison review. Identification uses conditional support at the observed
  Person totals and rosters, including incomplete assignments. Identified
  contrasts receive conditional-information SEs, pointwise 95% Wald intervals
  and joint covariance after finite-MLE and numerical checks; unidentified
  targets stay missing. The route estimates no population SD. These intervals
  are asymptotic, without a general finite-sample coverage guarantee.
  `mfrm_facet_intervals(fit)` reuses these saved comparisons at a requested
  confidence level, with existing interval plots, APA tables and RDS reuse.
  Missing estimates retain labeled rows without plotted points or zero-filling.

* `fit_mfrm(method = "CML", category_policy = "preserve")` now fits full-rank
  many-facet RSM/PCM conditional models, including supported facet interactions.
  A generalized ESF recursion conditions on Person totals, with analytic
  conditional information, identification and finite-MLE checks. The separate
  `mfrm_cml_fit` supports summary, joint covariance, conditional log likelihood
  and native-scale `mfrm_facet_intervals()` with existing interval plots.
  CML estimates no Person abilities or population SD; population-SD-standardized
  contrasts remain a full joint delta-method MML calculation. Anchors, GPCM,
  and ordinary fit diagnostics/reports are outside this CML scope.
  Conditional-information intervals are asymptotic. Help and returned notes
  distinguish increasing Person count from information growth in every
  identified direction, and explain the assumptions for response-dependent
  assignment, stopping and missingness. Comparison help gives the covariance
  transformation in the identified subspace.

* `fit_mfrm(fit_scope = "comparisons", comparison_review = review)` estimates
  identified RSM/PCM interaction contrasts when the full conditional model is
  rank deficient. The separate `mfrm_comparison_fit` retains the original
  likelihood, removes only redundant directions, and leaves unidentified
  targets missing. JML recession, numerical convergence/curvature and MML
  integration checks determine availability of descriptive estimates. Ordinary full-fit guards
  remain active; global MML identification and formal uncertainty are not
  certified. `summary()` and `interaction_effect_table()` retrieve the saved
  comparison table.

* `interaction_effect_table(fit, comparison_review = review)` and
  `summary(fit, comparison_review = review)` return explicitly weighted
  RSM/PCM joint-fit interaction contrasts. They check matching responses,
  origins, constraints and population design, full conditional rank,
  readiness, the retained likelihood/gradient and finer MML integration.
  Unavailable rows retain reasons and separate optimizer traces. Descriptive
  values preserve the MML identification limitations; no partial fit, SE,
  confidence interval or formal DFF test is supplied.

* `describe_mfrm_data(comparisons = ...)` reviews explicitly weighted
  facet-by-origin interaction comparisons before estimation. RSM/PCM reviews
  distinguish full, partial and target nonidentification under the declared
  conditional model; positive MML conditional ranks remain separate from
  marginal identification. GPCM interactions are explicitly unsupported.
  An optional `comparison_review` in `analyze_dff()` checks response/origin
  alignment and preserves this separate question in summaries, plots and
  reports. It adds no estimates or formal DFF inference and does not relax
  existing fit guards.

* Experimental score-based G intervals resolve their optimizer and reference
  distribution explicitly, including sessions where `stats` is not attached.

* DFF summary plots retain every source comparison, with display status,
  counts and comparison conditions. Both plot views retain source settings
  and residual response exclusions when available. Returned notes and the
  linking tutorial explain group direction, native units, sparse cells and
  the different linking screens. Existing estimates and inference rules are
  unchanged; a displayed contrast does not establish formal DFF inference.
  Base-R summaries reserve space for full comparison labels instead of
  clipping longer facet names at the default figure width.

* Saved conditional-normal population intervals now support `plot(ci)` for
  coefficients and `plot(ci, type = "contrasts")` for supplied linear contrasts,
  with `as_ggplot()` conversion. Unavailable intervals retain their points and
  labeled rows; predictor coding, contrast weights, numerical checks and
  experimental qualifications accompany the plot. Residual variance remains
  in the original table, on its own scale. No inference is recalculated and
  the statistical scope and coverage qualifications are unchanged.

* Wright maps accept `highlight = list(Rater = c("R01", "R02"))` through
  `plot()`. A `> ` label identifies selected facet levels in native and
  FACETS-style views and saved ggplot conversion. Selected levels are added
  to the usual `top_n` selection; counts and notes disclose the resulting
  display. Automatic native ranges now use locations before compaction.
  Fitted values and uncertainty are unchanged; highlighting is not a test.

* Rater-location profiles now omit magnitude bands by default. Explicit
  `show_bands = TRUE` draws neutral gray distance guides; numerical ranges
  replace `gentle`/`moderate`/`strict` in the returned `Band` column. These
  ranges do not establish homogeneity or practical equivalence. Old saved
  profiles retain their stored display when converted with `as_ggplot()`;
  recreate the plot from the saved fit and matching diagnostics for the new
  default. Estimates, uncertainty and equivalence calculations are unchanged.

* Planning help distinguishes MFRM simulation from observed-score G/D studies
  and explains why neither automatically requires the other. The G-theory
  tutorial adds an optional integer-count grid, and the linking guide explains
  the assumptions of connecting two schools through common recordings.
  Simulation help no longer describes higher rater separation reliability as
  universally better. These documentation changes add no estimator or plot API.

* Measurement plot notes now explain coded-score direction using the recorded
  facet signs, with Step and fit definitions in the relevant views. The rater
  profile axis now says facet location without assuming favorable score coding. Bias plots
  retain their interaction definition, sign and magnitude-only profile caveats.
  Help and the visual tutorial explain reverse-coded comprehensibility and
  calibrated judges without treating agreement as overfit or a bias flag as
  proof of unfairness. Plot notes stay outside figures by default; no new
  plotting argument, recoding, estimator or screening rule is introduced.

* Plot printing now counts legend entries correctly and points to public
  value-extraction and supported redraw methods. The visual tutorial uses a
  focused fit screen with an identifying table and shows how to save a figure
  with its values and interpretation notes. Estimation and fit rules are unchanged.

* The methods guide now connects construct claims to measurement, G-theory
  and rater-feedback evidence, including conflicting or inconclusive results.
  The G-theory tutorial adds a bundled-data covariance figure, its uncertainty,
  a reporting example and saved reuse. The visual guide explains what each
  figure can establish. No estimator or inference rule changes.

* Joint REML covariance profiles now inspect already computed equivalent
  optimizer candidates when the first selection fails numerical checks.
  The original convergence, gradient, covariance-cone and target constraints
  remain unchanged; worse local solutions are not substituted. Candidate
  selection is saved with diagnostics. This avoids unnecessary unavailable
  results without adding fits or claiming better statistical calibration.
  The methods guide now compares interval choices by source design and
  target, including incomplete data, and separates sampling SE from SEM.

* `confint(d, method = "delta", assumption = "normal")` adds G/Phi sampling
  SEs and explicit normal-reference intervals for supported ANOVA/MINQUE
  sources, including nested and incomplete designs. It uses the original
  component estimator and its full sampling covariance, without refitting.
  Defaults are unchanged. Summaries, forests and ribbons retain wide bounds
  outside `[0,1]`, with boundary and approximation qualifications. Undefined
  points or negative/nonfinite variance estimates remain unavailable; zero
  first-order variance retains SE zero without a degenerate interval.
  Balanced scalar component covariance now uses its algebraically equivalent
  independent-mean-square formula to avoid cancellation at zero boundaries.

* One-facet `confint(d, method = "spectral_f", assumption = "normal")`
  supplies G intervals for retained complete or incomplete source ratings,
  with paired plan differences through `reference`. A spectral F pivot
  removes facet effects without likelihood fitting or simulation. It requires
  independent normal Person effects and residuals, positive residual variance
  and remaining residual contrasts. Original ANOVA/MINQUE points, empty sets
  and numerical failures remain explicit in summaries and saved curve ribbons.
  Sparse sources can yield very wide intervals. Defaults are unchanged;
  this method does not supply Phi, SEM, nested or multifacet intervals.

* One-facet `confint(d, parm = "RelativeSEM", method = "residual_contrast",
  assumption = "normal")` now transforms the residual chi-square interval
  for complete or incomplete source ratings. Adding `reference` returns
  paired SEM differences from the same interval, preserving dependence.
  Original ANOVA/MINQUE points, separate projected centers and unavailable
  states remain visible in summaries, forests and saved curve ribbons.
  AbsoluteSEM, G/Phi and multifacet SEM require other methods. Normal residuals,
  positive residual variance and omissions preserving the model are required.

* `confint(g, method = "residual_contrast", assumption = "normal")` supplies
  residual-only intervals from retained G-study ratings, including identifiable
  incomplete sources. Variance limits use the exact normal-model chi-square
  law; signed covariance limits are conservative enclosures that can be wide.
  Original ANOVA/MINQUE points and the separate `IntervalCenter` remain
  distinct. Other components and designs without residual contrasts retain
  explicit unavailable limits. `plot(g, intervals = TRUE, ...)` supports the
  same method and saved reuse. G/Phi and other variance components require
  separate methods; normal-model coverage is not robustness to nonnormal residuals.

* G-study, D-study, paired-comparison and interval-summary printing now honor
  `digits` and `row.names` without changing stored values. Interval-summary
  rows identify their confidence levels, and nested plans explicitly label
  future child counts per parent. Completed approximate intervals retain their
  qualifications without being mislabeled as rows requiring diagnostic review.
  Ordinary G-study/comparison summaries reject
  unused arguments rather than silently ignoring them. The G-theory tutorial
  now leads with a bundled-data assessment question, interpretation and saved
  reuse; optional design/method examples follow, with formulas and detailed
  coverage evidence linked to the relevant help pages.

* `summary(d, intervals = TRUE, assumption = "normal")` collects supported
  G/Phi and measurement-SEM intervals with plan counts, methods, availability,
  source design and score weights. Select metrics with `parm`, or use the
  existing `reference` argument for paired plan differences. Saved interval
  objects can be supplied without recalculation and are checked against the
  D-study. Original estimates, separate REML centers and method qualifications
  remain visible. Ordinary `summary(d)` still returns the coefficient table.

* Individual G/Phi intervals from a balanced two-crossed-random-facet source
  now default to approximate MLS null inversion under `assumption = "normal"`.
  The existing `method = "mls"` name selects this explicitly for coefficients;
  SEM retains its distinct error-variance calculation. Original points, empty
  sets and numerical failures remain visible. Undefined squared MLS widths
  are conservatively retained in the interval, flagged by `FormulaComplete`
  and saved diagnostic ranges. Existing curve ribbons and method-comparison
  plots reuse these results. Coverage is approximate and small-source intervals
  can be almost `[0,1]`; this does not extend to incomplete sources or plan
  differences. Other explicitly selected methods remain available.

* `vcov(g, assumption = "normal")` now supplies sampling covariance and
  component SEs for balanced ANOVA and identifiable MINQUE(0), including
  incomplete and supported nested sources. One score, a fixed composite,
  or two-score covariance components can be selected. Gaussian product-moment
  correction removes plug-in bias in identifiable sampling-covariance entries;
  negative variance estimates, unavailable SEs and unidentified entries remain
  explicit. These describe the original estimator, not a separate REML fit,
  and do not validate normal-reference confidence intervals.

* Balanced two-crossed-facet G now accepts explicit experimental
  `confint(d, method = "score", assumption = "normal")`. The extended-nuisance
  score calculation reuses mean squares and preserves ANOVA points, including
  unavailable points. Numerical diagnostics distinguish empty confidence
  sets from failed calculations; existing forests and curve ribbons retain
  these states and support saved replay. Calibration is approximate, with
  limited finite-sample evidence. This is not automatic and does not supply
  Phi, SEM, incomplete/nested-source or plan-difference intervals.

* Two-score covariance components now also accept explicit experimental
  `method = "reml_profile"`, including identifiable incomplete sources.
  The joint normal model preserves signed bounds, original ANOVA/MINQUE
  estimates and separate REML centers. Numerical and boundary diagnostics,
  unavailable results and saved plot reuse remain inspectable; sampling SEs
  are not supplied. Chi-square calibration can under-cover with few facet
  levels and incomplete-source coverage is not established. This is never
  an automatic choice or a replacement for the balanced covariance methods.

* G-study scalar variance components now accept explicit experimental
  `confint(g, method = "reml_profile", assumption = "normal")`, including
  identifiable incomplete sources. Original ANOVA/MINQUE points, separate
  REML centers, accepted and candidate bounds, and numerical diagnostics are
  retained in squared score units. Sampling SEs are not supplied by this
  inversion. Existing component plots show these intervals while preserving
  negative points; saved plot objects retain their full interval evidence.
  The method is never automatic or a guarantee of nominal coverage.
  Two-score covariance uses the separate joint profile described above.

* D-study G/Phi and measurement-SEM intervals now offer explicit experimental
  `method = "reml_profile"` for balanced and identifiable incomplete or unequal
  sources under normal effects. It fits a separate REML target likelihood;
  original ANOVA/MINQUE points remain in `Estimate`, with REML estimates in
  `IntervalCenter`. The chi-square calibration is approximate, with limited
  finite-sample and boundary evaluation, and is never an automatic default.
  Saved intervals work with existing comparison plots and curve ribbons.
  Numerical diagnostics and exploratory failed bounds remain inspectable;
  no automatic fallback or matrix repair is supplied. Adding `reference`
  supplies G/Phi/SEM differences for one-facet plans or nested plans with
  a common future child count, including incomplete sources. This maps one
  shared profile interval, retaining interior extrema, original points,
  separate REML centers and approximate calibration. Other crossed/nested
  comparisons now constrain the difference directly in the shared likelihood,
  using optional nloptr. They retain detected accepted pieces, both zero-
  coefficient branches, failed candidates and score units for SEM differences.
  Their chi-square calibration is experimental, particularly at singular
  zero contrasts; finite-sample coverage is not established. These intervals
  are pointwise, not simultaneous, and are never chosen automatically.
  Shared-region G/Phi difference enclosures now intersect with the analytic
  parameter range implied by the future error weights, preserving their
  conservative coverage justification while removing impossible differences.
  SEM profiles constrain the error variance before taking square roots, retain
  score units, and search for the interval around the REML center; remote
  disconnected accepted regions are not enumerated. Boundary checks showed
  below-nominal inclusion, so numerical completion is not a coverage guarantee.

* Balanced ANOVA G-studies now provide signed two-score covariance-component
  intervals through `confint(g, score = c("V", "W"), ...)`, with explicit
  `gci` or `regression_region` choices under normal effects. Generalized
  intervals retain numerical endpoint bands and approximate coverage;
  conservative enclosures can be very wide. SE and interval availability are
  separate, including one-df sources. `plot(g, type = "covariance",
  intervals = TRUE, ...)` shows one interval plot for the selected pair;
  the point-only covariance heatmaps remain available. No matrix repair,
  clipping of signed bounds, or automatic method fallback is performed.

* Balanced nested plans with changed child counts now support explicitly
  selected experimental approximate-F and mean-square-region G/Phi intervals.
  The F statistic is derived from expected mean squares; help identifies
  discrepancies in the printed nested Phi formulas and reports boundary
  limitations found in a small evaluation. Automatic reduction methods keep
  their unchanged-child restriction.
  For balanced crossed/nested two-facet plans, `mean_square_region` also
  supplies SEM intervals and `reference` difference enclosures for G/Phi/SEM.
  These preserve a common simultaneous source-region event and known error
  ordering, but need not be the tight image of a difference and can be wide.
  Individual projections, reference counts, failed calculations and negative
  source estimates remain visible; no automatic fallback is introduced.

* Balanced nested D-study plans now have G/Phi intervals and `reference`
  comparisons when every future child count per parent equals the source
  count. Averaging within Person-by-parent cells yields the one-facet
  sampling law: G uses exact normal-model F inference; Phi retains its
  generalized interval and approximate coverage. Both random facets and
  their counts remain in the source and result. Nested SEM comparisons
  require a common future child count, which may differ from the source;
  they preserve the base chi-square/MLS qualification. Curve ribbons and
  saved comparison plots use the same `plot()` entry points.

* D-study `confint()` now accepts `parm = "RelativeSEM"` and `"AbsoluteSEM"`
  for complete balanced crossed or supported nested sources under explicit
  normal random effects. Intervals combine mean squares before taking square
  roots: chi-square for a single term, Graybill-Wang MLS for positive sums,
  Ting MLS for signed sums. Methods, raw variance bounds, SEs and unavailable
  states remain visible. MLS is approximate and can under-cover or be
  conservative; help reports a limited normal-model evaluation.
  One-facet sources also accept `reference` for SEM differences using their
  common error-SD interval. SEM curves support opt-in ribbons and saved reuse;
  forests and ggplot conversions retain score units and full interval ranges.

* G-study component plots now add pointwise intervals with
  `plot(g, intervals = TRUE, assumption = "normal")`, using the existing
  balanced ANOVA calculations for one score, multiple score panels or a
  fixed composite. Base and ggplot preserve negative estimates, all
  interactions/residuals and square-marked points without intervals.
  Returned data retain SEs, raw endpoints, methods and unavailable reasons;
  saved plots convert to ggplot without recalculation. These are approximate
  Ting MLS intervals (chi-square for residuals), not simultaneous intervals,
  covariance-element intervals or measurement SEM.

* One-facet Phi plan changes are now available with
  `confint(d, parm = "Phi", reference = 1, assumption = "normal")`.
  The full single-rating generalized interval is mapped to the difference,
  including interior extrema. The method retains approximate frequentist
  coverage, reference context, unavailable points and zero-boundary diagnostics.
  Joint numerical endpoint bands quantify Monte Carlo error across the returned
  comparisons. `plot()` and ggplot conversion use the difference scale.

* Multivariate D-study curves can now calculate and display coefficient
  ribbons with `plot(d, metric = "Phi", intervals = TRUE, assumption = "normal")`.
  Saved `confint()` results can also be supplied through `intervals` to avoid
  recalculation, using a list for several scores/composites. Curves use one
  interval method and confidence level.
  Source summaries, weights, plans and point estimates are checked before
  plotting. Base and ggplot ribbons preserve missing-interval gaps, open
  points and isolated interval segments. Notes and complete diagnostics
  remain available through `plot_data()` and saved plot objects.
  Forest plots remain the separate view for comparing interval methods.

* D-study `confint()` now accepts an explicit `method`. Existing one-facet
  defaults are unchanged. Complete balanced two-random-facet ANOVA sources
  can select experimental `"approximate_f"` and `"mean_square_region"`
  intervals for G or Phi, including a vector of methods for comparison.
  Neither is selected automatically. Results retain negative-component flags,
  availability, validation limits and numerical diagnostics. `plot()` on the
  returned table groups future plans on a single forest axis by default, with
  adjacent method intervals and a shared legend. Color, shape and line type
  distinguish methods; monochrome and optional separate panels are supported.
  Unavailable rows are retained, and ggplot conversion is available.
  Notes remain outside the figure. The tutorial compares methods on bundled
  data and explains why agreement alone does not establish robust coverage.

* G-theory interval help and the tutorial now distinguish available methods
  from two-facet coefficient methods still under evaluation. They explain
  the reasons for deferral, approximation and boundary limits, and the
  distinction between coefficient SEs and interval coverage, with references.
  The tutorial reports the conditions and denominators of a small method
  evaluation. This documentation adds no new interval estimator.

* One-facet D-study G changes can now be obtained with
  `confint(d, parm = "G", reference = 1, assumption = "normal")`.
  The interval is the range of the difference over the common F confidence
  set, including its interior extremum. Under the stated model it has at
  least nominal coverage and can be conservative. It preserves unavailable
  point/interval states and reference context; it does not add a delta-method
  extension or nested or two-facet coefficient intervals.

* `confint()` for a D-study now supplies G and Phi intervals from balanced one-facet
  ANOVA sources under explicit independent normal random effects. An exact F
  pivot uses the source sample sizes and transforms to each future count.
  Ratio bounds are intersected with the admissible parameter space before
  transformation; empty intervals and zero residual mean squares are explicit.
  `parm = "Phi"` uses the corrected agreement-ICC generalized pivot, then
  transforms its bounds to future counts. It preserves raw point estimates
  and returns Monte Carlo endpoint bands and boundary diagnostics. Phi's
  coverage is approximate; the caller's RNG seed is preserved. Two-facet
  coefficient intervals require the separate methods described above.

* `confint()` for a balanced ANOVA G-study now returns sampling SEs and
  pointwise Ting MLS intervals for one score or a prespecified fixed composite,
  with an explicit normal-random-effects assumption. One or two crossed
  facets and the supported balanced nested design are covered. Raw negative
  estimates and interval endpoints remain visible alongside the intersection
  with nonnegative variances. Residual intervals use chi-square quantiles.
  Covariance-element, incomplete-source and individual D-plan coefficient
  intervals are not supplied by this method. Component plots can show these
  intervals with the explicit option described above.

* New two-family scoring records can be reattached to their unchanged fit
  after switching R versions: their source identifier no longer includes
  the version of R that serialized it. Older identifiers still require a
  matching runtime; re-score with the retained fit if an older attachment
  is refused. Re-estimation is unnecessary.

* Shared random-rater models now flag fixed-facet levels containing only
  minimum or maximum scores. A converged optimizer and tiny positive
  curvature no longer admit SEs, intervals or prediction for this known
  nonfinite-calibration case. The check does not rule out all other forms
  of separation.

* Help for saved-result `plot()` methods now explains display controls,
  returned tables and interpretation, with examples on the same pages. Residual
  heatmaps are described in score units, QC detail panels show saved check
  explanations, and anchor forest whiskers are distinguished from intervals
  for changes between waves.

* Adaptive integration evaluates each Person's quadrature points together,
  including category moments, reducing repeated setup during fitting, scoring
  and uncertainty checks. Quadrature reviews avoid computing unused category
  moments. Integration orders and numerical acceptance criteria are unchanged.

* Design-evaluation summaries reuse resolved planning metadata when preparing
  their structural design review, avoiding repeated construction for older
  saved evaluations. Summary values, recommendations and notes are unchanged.

* Result bundles and expected-score pathway plots now reuse their stored
  diagnostics, including during plot export. If `compute = "never"` was used
  without diagnostics, these plots keep the curves and omit fit-statistic
  tables rather than starting a new diagnostic calculation. Supply matching
  diagnostics or rebuild with `compute = "auto"` to add those tables.

* Long plotting examples now wrap across lines so their arguments remain
  visible in the PDF reference manual.

* The linking-chain example now uses the default optimization allowance and
  shows fit summaries before linking. It explains the retained thin-support
  warning instead of proceeding from fits truncated at 30 iterations.

* Rater-effect help now distinguishes severity, central-category use,
  halo/horn impressions and local dependence. It explains the difference
  between the available descriptive indices and the unimplemented FM-SC
  and MRFM-H models, with links to their original references. Horn effects
  are not inferred from negative correlations or high severity. The methods
  guide also corrects an obsolete partial-anchor centering description.

* Hierarchical ICC plots now retain the interval method and unavailable-bound
  explanation in `plot_data(p, "notes")`, where `p` is the returned plot
  object. Figures remain free of embedded titles and footnotes, and saved
  plot returns preserve the notes without changing ICC estimates.

* Fit summaries, specification tables and methods text now report the
  centering convention. Native fits without convention metadata retain legacy
  remainder centering; missing external specifications are not inferred.
  The linking guide distinguishes reopening saved results, repeating the
  original analysis, changing constraints by refitting, and refreshing derived
  output. It includes a runnable migration example and explains why an old
  portable calibration cannot be converted by editing its labels. Recreate
  summaries from the saved fit to obtain the new descriptions; loading an
  existing summary or interval object does not update it automatically.

* One-family GPCM slope and curve intervals now retain an explicit statement
  that fixed anchors and group means are treated as known. The qualification
  follows interval tables, printing and returned plot notes without appearing
  inside the figure. Model/sandwich covariance calculations are unchanged.
  Help distinguishes these conditional intervals from the separate restrictions
  on anchored location/profile intervals, portable calibration and source
  readiness for fitted-object scoring.

* Refit DFF preserves baseline direct anchors on the target facet and JML
  Persons. Fixed target rows keep their estimates but report no group
  contrast or test; equal supplied anchors are not evidence of invariance.
  Rerun affected subgroup analyses from the baseline fit; loading saved DFF
  results does not repair them. Facets containing fixed values have no
  aggregate reliability, separation or strata, including when their fixed
  variance is represented as zero.
  Bias-screening help now describes the implemented bounded conditional
  likelihood, the positive-score direction of `Bias Size`, its conditional
  SE, and the post-shift meaning of `Obs-Exp Average`. The bias estimator
  itself is unchanged.

* Portable calibrations now record whole-facet, legacy remainder and uncentered
  facet constraints explicitly and check them against stored coordinates and
  free dimensions. Uncentered facets are no longer mislabeled as sum-to-zero.
  Existing files retain their original labels and scores; re-extract from the
  source fit to obtain the more precise record. Extraction also accepts a
  quadrature review whose retained fits have a different list order, while
  still checking each named grid, source fit, data and fitting settings.

* Constraint-fixed facet sums and group means no longer acquire tiny nonzero
  SEs from subtracting anchor offsets when constructing their Jacobian. The
  calculation now uses the analytic constraint derivative. Fixed targets keep
  SE zero and missing inferential bounds; free comparisons retain full joint
  model or sandwich covariance. Anchor-review messages also clarify that
  overlapping direct and group anchors retain both constraints.
  Recompute affected intervals from the unchanged fit; this derivative repair
  does not require refitting or change a saved result merely on loading it.
  Both covariance methods now use the fixed-target symbol in base-R interval
  plots, consistently with the ggplot view and returned explanation.

* New ordinary fits use `centering_scope = "facet"`: centered facets include
  direct anchors and group totals in their zero sum. Use `"remainder"` to
  reproduce earlier partial-anchor constraints. Contradictory whole-facet or
  group constraints stop before optimization. Anchor review, same-model refits
  and refitting exports retain the convention; older saved fits keep their
  original meaning. Population and uncentered-facet settings remain explicit.
  `make_anchor_table()` now preserves numeric precision by default
  (`digits = NULL`), avoiding rounding-induced conflicts when a whole facet
  is fixed. Explicit rounding remains available and requires constraint review.
  Unanchored native intervals accept the new constraint metadata; saved
  legacy inference retains its source identity when the equations are unchanged.

* MML now rejects Person direct/group anchors, which were previously retained
  in the configuration but ignored by estimation. Refit affected analyses
  using a supported specification; reloading an older fit does not repair it.
  JML remains the supported route for fixed Person values; the centering
  convention is now explicit for supported location constraints.

* Anchor help and the linking tutorial now start from packaged data and explain
  fixed values, group means, scale constraints and observed overlap separately.
  They document partial-anchor centering, settings that need to be resupplied,
  conditional uncertainty and the distinction from complete calibration reuse.
  They distinguish available conditional fixed-population intervals from the
  still unsupported anchored estimated-population interval configuration.
  Anchored-fit summaries identify fixed-value consistency checks explicitly;
  zero differences do not demonstrate stability across administrations.

* Beginner examples explain loading packaged data, column roles and result
  interpretation. Existing crossed/nested G-theory and agreement teaching
  data are now available through `data()`, so these workflows no longer
  require writing data-generation loops. The classification example starts
  from the packaged rating data and its planned assignment roster.

* The testlet and shared-rater tutorials reuse their packaged synthetic
  calibration results. Shared-rater scoring and bootstrap displays also reuse
  saved results, while retaining complete recomputation calls and interpretation
  limits. Rebuilding these displays no longer repeats the corresponding fits.

* `export_mfrm_results()` now generates saved-result replay scripts for every
  supported model. Ordinary-model scripts no longer depend on `fit` or
  `diagnostics` variables from the original session. Requesting `"replay"`
  includes the complete RDS; use `export_mfrm_bundle()` for refitting.

* The individual-rater tutorial distinguishes supported RSM/PCM sheets from
  experimental two-family GPCM sheets and unavailable model combinations.
  It also shows how to retain interpretation notes when sharing selected
  tables as CSV files, separately from the full analyst archive.

* Feature comparisons accept unchanged imputation results reopened from
  separate RDS files. Recreated formula environments no longer cause a false
  model mismatch. Retained data, settings, formula expressions and completions
  are still checked; comparison does not evaluate or refit the formulas.

* The external-feature tutorial uses a larger imputation-diagnostic figure
  so variable and summary labels remain readable. Its saving example retains
  both compared analyses and distinguishes a summary CSV from the full RDS.

* Workflow help now correctly lists ggplot conversion for saved pooled
  multiple-imputation intervals. The assigned-score tutorial keeps the
  probability-figure description outside the plotting area.

* Model guides distinguish ability scoring from response-probability prediction.
  Shared-rater help explains that a supplied scoring table replaces source
  responses, even when rater IDs are reused. Testlet help clarifies that a
  column used both for local membership and as a fixed facet requires known
  fixed levels when scoring new data.

* Shared-rater and testlet calculations reuse repeated likelihood and scoring
  operations to reduce computation time. Shared-rater calculations retain
  each Person's rating counts and shared-rater dependence; testlet scoring
  reuses exact density values within each scoring call. Integration settings,
  numerical checks and uncertainty interpretations are unchanged.

* Testlet and shared-rater reports distinguish response diagnostics that were
  not supplied from diagnostics that could not produce usable results. The
  report points to the separate diagnostic call without computing it. Ordinary
  residual diagnostics remain unsupported for these models; the available
  posterior summaries are descriptive and have no ordinary-model cutoffs.

* Calibration help explains validation and freezing in terms of the saved
  model and settings.

* Reports and numerical-check notes give plain-language review actions and
  limitations. The reporting checklist describes current plot support,
  including the absence of a native 3D information surface.
  Manifest help distinguishes optimizer convergence from numerical checks for
  inference and from the statistical adequacy of an interval method.

* Portable score plots now use their dedicated ggplot conversion even with
  explicit `component = "data"`, preserving interval/precision/edge-mass axes,
  stored bounds and review states. Other components and unsupported conversion
  arguments are rejected instead of entering generic graphics. Help distinguishes
  saving an analysis, a reusable calibration and an already scored batch, and
  explains how to handle older files with incomplete source or score records.
  Summary-appendix exports also retain analysis tables named `overview` or
  `reporting_map`, which were previously omitted as reserved metadata names.
  Conflicting analysis filenames receive a `table_` prefix and, when needed,
  a numeric suffix; existing metadata filenames remain unchanged. Inspect
  `written_files` for the actual CSV paths. HTML retains the same analysis tables.

* Severity profiles now check supplied diagnostics against the fitted analysis
  and retain SE sources, interval-use restrictions and readiness in saved plot
  data. Fixed values, points without usable SEs and unavailable estimates remain
  distinct; finite points without intervals are included in the axis range.
  Dedicated ggplot conversion preserves stored bounds, rows and guide bands.
  Help clarifies that these are normal uncertainty bands, not likelihood-profile
  intervals or pairwise tests; exploratory JML uncertainty remains exploratory.

* Bias scatter and ranked ggplot views now preserve their selected axes,
  saved screening thresholds, point colors and conditional interval status.
  Saved reports reject conflicting plotting-time thresholds, and their notes
  use the stored rules. Help corrects the scatter-axis description and explains
  that repeated views do not provide independent evidence of bias. Generic
  conversion of the other bias views is explicitly unavailable; their base
  and custom-data routes remain available.

* Halo ggplot views now show the selected correlation distribution, rater
  metric or signed retained-edge matrix instead of unrelated count columns.
  Figures retain settings, pair denominators and interpretation caveats;
  unassessed values remain distinct from no configured flag. Matrix scales
  remain -1 to 1 with grey missing/filtered entries. Network figures omit
  automatic titles and use base/custom graphics. Help distinguishes observed
  score associations, residual dependence and indirect rater centrality.

* Threshold ladders now draw one set per step-owner level when GPCM step and
  slope owners differ. Missing steps remain explicit, and ladder notes explain
  the centered-step target and limits of disorder screening. Rating-scale
  threshold figures omit an automatic title and size their bottom margin from
  the labels, preserving a readable panel on smaller devices.

* Category-count and rating-scale threshold ggplot conversions now preserve
  the selected values, expected counts, owner/step order, missing values and
  palettes. Previously generic conversion could select category codes instead
  of counts or steps. Automatic conversion of category-curve bundles, boundary
  views and threshold ladders is now explicitly unavailable because generic
  conversion lost their targets or grouping; use their existing base renderer
  or `plot_data()` for custom graphics. Fitted CCC conversion remains available.

* Ordinary-model reports and APA drafts now distinguish observed agreement
  from the fitted-model expected baseline, show the stored matched-opportunity
  denominator and available-pair count, and describe unavailable results
  explicitly. Aggregate reports retain the numerical agreement summary and
  its interpretation in table/CSV output as well as Markdown/HTML. Reporting
  reuses stored diagnostics without recalculating agreement or fitting;
  existing saved APA drafts need regeneration to acquire the revised wording.

* Agreement figures converted with `as_ggplot()` now draw the selected exact
  agreement, correlation or score-difference metric. The generic conversion
  previously selected the first numeric column, which could plot matched
  counts instead. Dedicated views preserve pair order, model-expected values,
  missing comparisons, colors and returned notes. Agreement heatmaps retain
  the facet labels and fixed metric scale; monochrome now works in both
  renderers, and base heatmap margins accommodate the level labels.

* Plots of saved agreement results now use their saved review thresholds
  instead of default reference lines. Conflicting plotting-time threshold
  overrides are rejected so the reference and stored flags cannot disagree;
  recreate the agreement result from the existing fit to change the rules.
  Help and the visual tutorial explain saved plotting and the combined flag.
  Agreement calculations, pair counts and result-table values are unchanged.

* Residual-correlation heatmaps converted with `as_ggplot()` now label both
  axes with the selected facet and use a fixed -1 to 1 correlation scale.
  Unavailable cells remain distinct from zero correlations. Q3-style results
  and heatmap payloads retain aggregation and interpretation notes when saved;
  figures still omit automatic titles and captions. Residual-PCA reference
  legends now wrap to two rows to fit smaller figures.

* Residual-PCA and Q3-style help clarifies analysis units, supported saving
  and plotting routes, reference limits and literature links, and distinguishes
  the loading-column limit from the full eigenvalue spectrum. The visual
  tutorial computes each PCA scope once instead of twice. Correlations,
  eigenvalues, flags and reference calculations are unchanged.

* Agreement guidance now distinguishes chance agreement, fitted model-expected
  agreement, pairwise ICC and residual association, with explicit matching,
  missingness and saving rules. `interrater_agreement_table()` help clarifies
  repeated-cell averaging, pair-context denominators, weighted mean correlation,
  fractional expected counts and the summary's scope when `top_n` is used.
  Calculations, defaults and returned values are unchanged.

* G-theory help now separates component sampling uncertainty (not supplied),
  projected measurement SEM, and approximate intervals for prespecified plan
  differences. The plan-comparison help gives its quadratic-form covariance,
  paired delta calculation and methodological references. Main-effects G/D
  help points to the separate interaction-capable route, including for one
  score. These documentation changes do not change numerical calculations
  or extend interval support.

* **Breaking plotting API change:** the 33 standalone drawing `plot_*()`
  helpers are no longer exported. Use `plot()` on the fitted model or saved
  analysis result; `?mfrmr_plot_migration` maps every former helper to its
  replacement. Only `plot_data()` and `plot_data_components()` retain the
  public `plot_` prefix. Fit-level views use `type`, applicable `view`/`metric`
  controls, and `comparison` for paired models or ordered trajectories.
  Existing renderers, numerical results, notes and model restrictions remain;
  drawing a comparison does not establish linked scales or run an IC/LRT test.
  Automatic residual Q-Q conversion with `as_ggplot()` is now refused: its
  generic fallback produced a bar chart rather than a Q-Q plot. Use the base
  renderer or specify both quantile axes from `plot_data()` for a custom graph.

* Older reference material now distinguishes the one- and two-family MML
  engines, the multivariate G/D route and the separate main-effects adapter.
  The printable quick reference is synchronized with its source and includes
  G/D planning and interpretation boundaries. Marginal-fit documentation now
  explains that its diagonal grouped-count variance is a screening approximation:
  conditional independence given ability does not remove covariance induced by
  integrating shared ability. Calculations and existing result schemas are unchanged.

* The beginner workflow now explicitly selects either the practice analysis
  or the user's CSV before diagnostics, reporting and saving. The CSV route
  does not inherit the practice assignment roster. Model-family guidance also
  distinguishes ordinary RSM/PCM diagnostics from two-family GPCM and corrected
  JML results that share the same fit class.

* The output/workflow guides now distinguish ordinary RSM/PCM from the
  one-family, two-family and corrected-JML routes sharing `mfrm_fit`.
  They align G-study component/covariance views, D-study curves and existing
  PCA/cluster ggplot conversions with the implemented methods, and correctly
  describe two-family profiles under fitted fixed or adaptive integration.
  These are guidance corrections; no estimator or supported design changes.

* `fit_measures_table()` now enforces the ordinary-model restriction even
  when diagnostics are supplied: two-family GPCM, corrected JML, testlet and
  shared-rater fits cannot acquire ordinary fit indices from another result.
  Existing ordinary and detached-table outputs retain their values and
  interpretation. Help clarifies that interval plots use the returned table
  bundle and keep their title and notes outside the figure by default.
  Converting these plots with `as_ggplot()` now preserves each facet/level,
  the saved normal bands and fixed/unavailable symbols; the generic fallback
  previously dropped the bands and collapsed levels into facet bars.

* `mfrmr_methods_reference` now explains how to review rating exposure,
  confounded facet effects despite connected graphs, planned assignments,
  owner-level category support, missingness and anchor uncertainty before
  fitting. Flexible input designs and available estimates do not establish
  general robustness or extend the requirements of downstream outputs.

* `plot(fit, type = "reliability")` now rejects corrected-JML and two-family
  GPCM fits even when diagnostics are supplied. Plotted values retain existing
  precision/use flags in their returned table, and ordinary-JML notes identify
  exploratory screening. Reliability, separation and strata calculations are
  unchanged; corrected-JML `RootSE` is not substituted for their input SEs.

* The interval guide now includes the experimental RSM/PCM population route.
  Its help explains numerical availability, unavailable bounds, log-variance
  display limits and the separate interpretation of conditional Person scores.
  Estimation, eligibility checks and reported numerical results are unchanged.

* Multivariate G-study plots now accept multiple original `score` names for
  component panels, or `type = "covariance"` for component matrix heatmaps.
  `ncol` controls panel columns. Base and ggplot preserve shared scales, signed
  values, full-matrix diagnostics and the combined highest interaction/residual.
  Heatmaps show numeric cell labels and support color or monochrome without
  standardizing scores, repairing matrices or estimating sampling intervals.
  Existing single-score/composite calls and saved component plots remain usable.

* Multivariate D-study curves now compare explicitly selected original scores
  or named composites in one figure, with `metric` selecting G, Phi, relative
  SEM or absolute SEM. Each line keeps the other facet count fixed; identities,
  weights and unavailable points remain in the plot data. Optional `reference`
  values are user-chosen guides, recorded with their interpretation. The
  main-effects D-study plot now also omits automatic 0.70/0.80 lines; request
  `reference = c(0.7, 0.8)` to restore them. Surface plots reject such guides.
  Base and ggplot figures return titles and cautions separately. Estimation
  and D-study projection formulas are unchanged. `panel_by` can now separate
  the non-horizontal facet count into panels; `ncol` controls their wrapping.
  Base and ggplot renderers retain shared axes, consistent score styles,
  nested count units and panels with unavailable estimates.

* `plot()` now displays raw variance components from a saved multivariate
  G-study for selected scores or one prespecified weighted score. Crossed, nested and
  incomplete sources retain their component identities, negative estimates
  and matrix-diagnostic cautions. Base and optional ggplot renderers preserve
  values and return explanations outside the figure. The G-theory tutorial
  connects this view to D-study projections and paired plan differences;
  percentage/error-share interpretations are not added. Component intervals
  are opt-in for eligible balanced ANOVA results, as described above.
  Labels explicitly name the highest-order interaction combined with residual
  variation, while `Source` identifiers and numeric values stay unchanged.
  The tutorial shows all crossed interactions and both SEMs, distinguishing
  residual variation, relative/absolute measurement error and sampling uncertainty.

* The G-theory planning tutorial now connects data/design review, G-study
  estimation, D-study projections, prespecified paired comparisons and saved
  reopening, with explicit incomplete-source and unavailable-result examples.
  G-theory functions have their own reference-index group. Returned plot notes
  no longer repeat captions that differ only in outer whitespace or recombine
  separate component warnings during D-study ggplot conversion. Saved prose,
  numerical values and estimation methods are unchanged.

* `classification_table()` now retains `SourceStatus` in category- and
  sum-score probabilities, cut-point precision and the population classification
  matrix, as well as its existing index tables. Numerically available CSV rows
  therefore retain a source-review flag. Returned notes and help distinguish
  joint population errors from error rates conditional on a true category.
  Classification probabilities, decision rules and eligibility are unchanged.

* Individual rater sheets now include an `interpretation_notes` table in
  their object and `output = "tables"` results, retaining the overall review,
  audience-specific guidance and section limits when tables are saved or
  shared. Only the sheet's recipient-safe explanations are copied; private
  source notes remain excluded. Numerical tables and HTML/Markdown presentation
  are unchanged.

* Experimental-method limitations now travel as interpretation-note tables in
  two-family GPCM and population-interval result/report exports, including
  GPCM results with no interval attachment. Corrected-JML notes also make
  incomplete method qualification explicit. Reporting help distinguishes this
  from a partial interval with a missing endpoint. Numerical results, interval
  eligibility and the supported model scope are unchanged.

* `help("mfrmr_methods_reference")` now maps assessment questions to
  calculations, interpretation limits and methodological references, with
  DOI links for the foundational response models and MML framework. README,
  package/workflow help and the introductory vignette provide shorter reading
  paths to the existing examples, data checks and result explanations.
  This changes documentation only, not calculations or supported inference.

* Exported Stan and brms bundles now put the sampling review and next action
  before numerical summaries. Parameter tables combine original labels,
  reference/direction notes, fixed/estimated status, 95% credible bounds and
  diagnostics. Fixed values have no estimated SD or interval (`NA`);
  `interval_computed` now refers to varying parameters only. Bounds for varying
  parameters remain visible with diagnostic concerns and explicit interval
  status. `Rscript run.R --summarize` rebuilds both formats' reports from saved
  draws without compiling or sampling. Keep the matching bundle and new
  `summarize.R` with the fit. Finite endpoints and a clean sampling screen do
  not certify inference validity or frequentist coverage. Priors, likelihoods
  and sampler settings are unchanged.

* `read_mfrm_posterior()` connects the matching exported summary/model-data
  pair to `summary()`, `plot()`, `confint()` and `mfrm_report()`, including
  ggplot conversion, long plot tables and appendix exports. Fixed values,
  credible-interval semantics, priors and sampling review survive saved reuse.
  The reader requires a model-data checksum from the current summary runner;
  regenerate older development summaries with the current `run.R` and
  `summarize.R` using `--summarize`. No new sampling is needed. This is a
  read-only marginal-summary object, not a native fit or a joint-draw interface;
  new interval levels, contrasts, predictions and native MML diagnostics are
  not reconstructed. A small synthetic example is included.

* `mfrmr_interval_guide("gpcm")` now includes one-family location intervals,
  experimental two-family component intervals, slope profiles and corrected-JML
  root uncertainty. The latter has no default confidence level or structural
  bounds. Help and the GPCM tutorial distinguish numerical availability,
  coverage evidence, partial endpoints and posterior uncertainty; no estimator,
  interval calculation or admission rule changes.

* Category reports now retain raw and weighted category counts within the
  fitted step-owning facet, with explicit denominators, original score labels
  and zero-exposure states. Stored step SEs/intervals are joined by owner and
  step with estimate/source checks; unverified legacy sources remain review
  outputs. Bare-fit category reports no longer request full diagnostics or
  covariance implicitly; supply `diagnose_mfrm()` output for those companions.
  Tables, curve-profile maps and saved plot data preserve model roles. Ordinary
  PCM/GPCM fitting now rejects multiple or missing explicit step-owner names.
  No threshold model, interval method or two-family diagnostic scope is added.

* Saved slope intervals distinguish computed, partial, review-only and missing
  bounds in print, tables and plots. Partial endpoints remain visible without
  implying an infinite bound or a calibrated one-sided interval. Bias plots
  preserve profile search-limit status with dashed whiskers and retain failed
  profile bounds without substituting SE intervals. Returned notes state that
  bias uncertainty conditions on the fitted calibration.

* Reliability, separation and Strata retain their existing numerical formulas
  but are labelled approximate summaries. Source SE eligibility is recorded
  separately from formal inference for these derived indices. Reports explain
  the omitted covariance, EAP distinction and hypothetical meaning of Strata;
  diagonal facet heterogeneity statistics are labelled screening summaries.

* A new beginner vignette, "Rating agreement and score comparisons", works
  through composite-key matching, pairwise QWK, correlation and mean shifts,
  missing-rating audits, and comparisons with a second reference. It explains
  how panel alpha/AC2 and model-adjusted severity answer different questions,
  and shows how to save and reopen the full results. Custom ggplot examples
  display matched-score counts and pairwise QWK with sample sizes, preserving
  negative estimates and unavailable rows. Titles and notes stay outside the
  drawing and are saved with its data. The README, ordinary workflow, reporting
  and visual-diagnostics guides link to this executable, model-free lesson.

* Results export now saves every member of a multi-plot `fit` bundle as a
  separate PNG with its matching title/notes CSV. The starter index links
  each image to its notes. Drawing runs once per route; inconsistent page
  counts are reported as failures, and failed rendering preserves existing
  files. Single-figure filenames and figure annotation defaults are unchanged.
  Native Wright maps also reserve space for full axis labels when shortening
  them would make different groups indistinguishable.

* `mfrm_facet_intervals(..., scale = "standardized")` now reports explicit
  zero-sum severity contrasts in population SD units. The full joint delta
  method includes estimated SD uncertainty and its covariance with severity;
  it does not divide native interval limits by an estimated SD. Supported
  model-based routes use fixed N(0,1) or an estimated intercept-only normal
  population. Latent regression and standardized sandwich intervals are
  excluded. Saved tables, reports, exports and plot axes retain the scale and
  unqualified sampling-coverage status; plot explanations remain returned notes.
  Help and returned cautions explain the positive true-SD and marginal-model
  assumptions. The scalar delta variance includes the cross-covariance term,
  so accounting for estimated SD need not increase the SE. MML help also
  distinguishes finite-quadrature sensitivity checks from continuous-model
  inference and explains why increasing sample size alone does not remove
  integration error, misspecification or lack of identification.

* Two-family adaptive MML now accepts a named list of native starting vectors
  in `gpcm_mml_start`. Every supplied start is evaluated and optimized under
  the requested adaptive objective; the existing selection, convergence and
  output checks still apply. Summaries identify explicit initialization, and
  saved replay/quadrature refits retain the exact roster. This permits values
  from a separate numerical calculation to enter the standard fit, interval,
  scoring and portable-calibration workflows without importing its qualification.

* Experimental two-family GMFRM interval tables now retain
  `CoverageQualified = FALSE` separately from numerical `CIEligible`.
  This distinction persists for component Wald/profile and location/contrast
  intervals in saved results, reports and exports; it does not estimate coverage
  or change interval admission. Nonrepresentable location/contrast uncertainty
  is now missing in auxiliary model columns as well as the primary columns,
  preserving finite point estimates and the failure status.

* Package, fitting and external-model help now explicitly name the existing
  logistic metric (L-metric, D = 1). Guides and exported READMEs distinguish
  response-link scaling from normal ability distributions, explain slope-aware
  logit changes and caution against comparing unconverted parameters from
  D = 1.7 software. Estimation and model code are unchanged.

* `export_mfrm_model(..., format = "stan")` now supports native two-family
  GMFRM MML specifications. Exported likelihoods preserve product slopes,
  the fixed N(0,1) population, source location constraints and step ownership.
  `log_slope_prior_sd` sets explicit marginal log-slope priors;
  `slope-levels.csv` maps both families to Stan parameters. No fitting is done
  during export. Bayesian intervals do not qualify the existing MML intervals;
  two-family brms export remains unsupported.

* Figures now omit titles, subtitles and explanatory footnotes by default.
  `plot_data(p, "notes")` returns descriptions and cautions for base plot
  payloads and mfrmr ggplots, including after saving and reloading. Suggested
  titles remain in `plot_data(p, "title")`; printing a base plot payload
  displays its notes. Axes, legends, panel identifiers and unresolved-value
  markers remain. Existing `show_title`/`show_notes` switches now default to
  `FALSE`; `title`/`main` text is retained as metadata. Explicit switches and
  ggplot `labs()` support optional in-figure text. Numerical warnings and
  inference eligibility are unchanged.
  `export_mfrm_results(include = "plots")` also saves each figure's suggested
  title and notes in a matching CSV, linked from the starter index and listed
  in the export manifest, so interpretation cautions travel with the images.

* Experimental two-family component profile intervals now preserve adaptive
  integration for direct MML fits. Every constrained trial updates the
  integration nodes, scales and Jacobians; nuisance reoptimization, two starts
  and finer-order checks remain required. Saved outputs identify the integration
  method and engine. This enables an explicit sensitivity comparison with
  log-Wald intervals; finite-sample coverage improvement is not established.

* `export_mfrm_model()` writes self-contained Stan or brms model bundles for
  unit-weight RSM/PCM and one-family GPCM MML specifications on the fixed
  N(0,1) ability scale, including supported fixed facet/group anchors and
  retained step constraints. Both retain
  the response likelihood, category mapping, facet signs, centering and step
  ownership, with explicit matching Bayesian priors. Bundles include data,
  level mappings, executable R scripts and conditional posterior predictions.
  One-family slopes have geometric mean one. Affine bundles record prior means,
  marginal SDs and fixed coordinates in `levels.csv` and `steps.csv`; anchors
  are exact inputs and do not propagate external calibration uncertainty.
  No new native step-anchor or two-family-anchor fitting route is added.
  Unsupported extensions are rejected before writing. Posterior credible
  intervals are distinct from MML confidence intervals; no estimator is run
  and no Stan dependency is installed by export. Executable help and workflow
  examples show both formats and their mappings/priors; the output guide lists
  this external-model route separately from saved-analysis archives.


* Experimental two-family GMFRM intervals retain deficient empirical Person-score
  rank as a diagnostic and warning instead of rejecting solely on that rank.
  The model-based covariance is the inverse full observed Hessian, not an
  empirical score covariance. Positive unregularized information, derivative,
  stationarity and quadrature checks still apply; this change does not qualify
  small-sample coverage or add sandwich inference.
* `mfrmr_output_guide()` now locates observed agreement, score comparisons,
  centrality, criterion-referenced rater accuracy, fixed-plan classification
  and MH/ordinal DIF sensitivity in its existing scopes and psychometric map.
  Each route states its input, interpretation limits and saved-output workflow;
  agreement, severity separation and the different accuracy targets stay distinct.

* Residual `analyze_dff()` and `dif_interaction_table()` now match the supplied
  response multiset by Person, all facets, score and weight. Repeated rows no
  longer inherit the first matching group, and subsets no longer include
  unselected repetitions. A `cases` audit preserves exclusions; excess copies
  are rejected. Counts of distinct Persons are separate from response counts.
  Earlier saved residual results must be recomputed from their fit and data;
  no model refit is needed. Repeated-rating inference remains unsupported.

* `ordinal_logistic_table()` adds independent-Person cumulative-logit group
  comparisons through the optional `ordinal` package: uniform, nonuniform
  and joint LR tests, Nagelkerke R-squared changes, common-range effect
  contrasts and a separate proportional-odds diagnostic. Declared category
  order, matching ranges, exclusions and numerical failures remain auditable.
  Inference is asymptotic; repeated ratings, uncertain matching and unequal
  rating-plan exposure are not handled. Existing DFF inference is unchanged.

* `mantel_haenszel_table()` adds binary conditional group comparisons with
  explicit matching strata and one response per Person-target cell. It
  retains MH and conditional-ML estimates, exact or asymptotic inference,
  pointwise intervals, adjusted p-values, stratum counts and a row audit.
  Missing/no-overlap/constant strata and unavailable results remain visible;
  repeated Person-target ratings are rejected. Matching validity, common
  exposure and independence remain substantive requirements. This does not
  qualify residual/refit DFF inference.

* DFF documentation now explains the matching, rating-plan and independence
  conditions needed for MH/ordinal-logistic comparisons, including
  nonuniform effects and the different cumulative-logit assumptions. The
  typical workflow correctly describes residual contrasts as having no flags.
  A runnable shared-roster comparison now aligns effect directions, displays
  method-specific counts/exclusions and numerical status, and retains complete
  results for saving. Native effect units and multiplicity families remain
  separate; no combined robustness classification is assigned.

* `classification_table()` evaluates classification accuracy, consistency
  and ability-scale conditional SEM for an explicit fixed RSM/PCM rating
  plan. Lee-style summed-score recursion and the Rudner ML normal
  approximation retain separate operational cuts and true ability cuts.
  Population indices require supplied ability weights. Category/sum
  probabilities, joint classification cells, exact-cut precision, plan,
  category mapping and source readiness are retained for review and export.
  Zero-information normal results remain unavailable. Calibration uncertainty,
  random plans, weighted/repeated cells and automatic cut selection are not
  included; the two methods need not describe the same decision rule.

* `pairwise_kappa_table()` adds quadratic (QWK), linear and unweighted Cohen
  agreement on each pair's common complete units. `gwet_ac2_table()` adds
  panel AC2 and unweighted AC1, using equally weighted unit proportions that
  include singletons. Both require a declared category universe and retain
  weights, component quantities, matching audits and unavailable-result reasons.
  Numeric categories use their spacing; labels use the declared order.
  Negative estimates remain visible. The distinct chance/missingness models,
  influence of unused categories on AC2 and constant-rating interpretation
  are explicit. No interval, quality cutoff or model fitting is added.
  Existing score-comparison QWK uses the same distance calculation, fixing
  undefined results when finite opposite-sign category endpoints have an
  overflowing span; the statistical definition is unchanged.

* `krippendorff_alpha_table()` computes pooled agreement across raters for
  nominal, ordinal, interval and ratio scores. Only units with at least two
  nonmissing ratings contribute to observed and expected disagreement;
  singleton ratings, missing scores and absent assignments remain auditable.
  Ordinal scores require a declared category order and use pairable-frequency
  distances, rather than QWK weights. Repeated cells are rejected, negative
  alpha is retained, and undefined coefficients have reasons. No confidence
  interval, quality cutoff or reference-accuracy claim is added. Alpha and
  pairwise ICC share input matching; the ICC formulas and outputs are unchanged.

* `pairwise_icc_table()` computes two-way observed-score ICCs on each rater
  pair's common complete units, without fitting an MFRM. It distinguishes
  absolute agreement from consistency, single ratings from the mean of two,
  and selected-rater from random-rater interpretations. Pair-specific samples,
  scaled ANOVA ingredients, missing scores and absent assignments are retained.
  Negative ICCs remain visible; undefined cases retain reasons. Repeated cells
  are rejected. No interval, quality cutoff or panel ICC is inferred.

* `fit_rater_accuracy()` models exact agreement with supplied reference ratings
  using a binary additive Rasch model, with raters as the objects of measurement.
  It returns accuracy logits, benchmark-facet difficulties, match/exclusion
  counts and source review status. JML retains infinite estimates for extreme
  raters; MML explicitly estimates a normal rater population and returns EAPs
  and conditional posterior SDs. Conflicting references and duplicate cells
  are rejected. References are fixed; no accuracy cutoff, confidence interval
  or rater-by-domain interaction test is added.

* Reliability reports, APA text and diagnostic summaries now explicitly
  distinguish high rater separation reliability (distinguishable severities)
  from high Person separation reliability (distinguishable abilities). Low
  rater values can also reflect imprecision and do not establish agreement.
  `plot(fit, type = "reliability")` retains this interpretation in its help,
  caption and saved notes, and removes the universal .70/.90 quality lines.
  Numerical reliability estimates are unchanged.

* `rater_centrality_table()` places RSM raw-residual versus expected-score
  correlations and Infit/Outfit alongside rater-specific PCM threshold SDs
  and variances, following the separate model bases in Eckes and Jin (2021).
  It retains category counts, exposure, step order, source review status and
  unavailable-index reasons. Two supplied fits must share the retained rating
  roster and category/facet constraints. The function computes descriptive
  indices from existing unit-weight ordinary JML/MML fits without refitting;
  it does not classify raters, estimate Bayesian centrality or adjust abilities.

* `score_comparison_table()` compares aligned scores and reference ratings
  from human raters, automated systems or other scoring procedures. It returns
  QWK, correlation, baseline degradation and both pooled- and reference-SD
  mean differences. No thresholds apply by default; `screen = "williamson2012"`
  selects the automated-scoring examples, with optional custom limits.
  Optional saved RSM/PCM severity and Infit/Outfit retain source scope and
  review status. The earlier `automated_scoring_table()` interface preserves
  its arguments, output names and default screen through the same calculation.
  No model fitting, recomputed diagnostics or deployment approval is implied.

* Conditional ability-probability tables and scoring guides now explain that
  numerical checks do not validate the scoring prior for a new population.
  Estimated-prior uncertainty remains excluded; a larger calibration sample
  does not correct population mismatch, and a fixed reference prior is not a
  universal remedy. Summaries, table notes and exported table descriptions
  share this interpretation without changing saved probabilities or rescoring.

* Experimental RSM/PCM population and native facet intervals now refuse an
  interior fit when setting residual variance to zero, with the other estimates
  retained, materially improves its likelihood. The comparison accounts for
  the observed integration-order difference and also runs when population
  intervals reuse saved facet covariance. Points and refusal reasons remain
  available. This is a candidate check, not a full boundary optimization or a
  sampling-coverage qualification.

* `mfrm_population_intervals()` accepts named `contrasts` of fitted population
  coefficients. Pointwise intervals use the coefficient block of the full
  nuisance-adjusted covariance. Reports and CSV exports retain contrast weights,
  covariance and unavailable targets without refitting. `vcov()` retains its
  coefficient/log-variance coordinates. These remain experimental intervals
  with unqualified sampling coverage; contrasts concern conditional population
  means, not individual ability differences or causal effects.

* Ordinary RSM/PCM scoring accepts `cut_scores` in `predict_mfrm_units()` and
  `score_mfrm_calibration()`; ordinary RSM `score_mfrm_persons()` also supports
  them. Saved `decisions` tables contain continuous probabilities above each
  cut and within ability bands, with cut labels, priors and availability.
  Summaries and `apa_table(..., which = "decisions")` reuse these probabilities;
  prediction exports retain them. Portable score batches and their saved
  summaries now also support `build_summary_table_bundle()` and
  `export_summary_appendix()`, retaining unavailable Person dispositions.
  Scoring help pages share the interpretation and link to executed examples.
  Calibration, population and cut parameters
  are held fixed. These are not tests of Person differences, future pass rates
  or automatic CEFR classifications. Calls without cuts retain their behavior.

* Corrected JML reuses transition probabilities across correction steps at
  the same parameter vector. The optional cache is limited to one quarter
  of `mfrmr.max_information_bytes`; larger calculations retain the existing
  blockwise computation. State limits, equations and numerical checks are
  unchanged, and cached values are not shared across parameter evaluations.

* Fixed-grid RSM/PCM MML with an estimated population and fixed-grid GPCM MML
  can use the existing guarded curvature restart after ordinary polishing
  stalls. The step uses all free coordinates, including population parameters,
  and retains the 64-parameter limit and all convergence, integration and
  information checks. Previously these fits could remain stalled with optimizer
  code zero and a gradient too large for inference.

* Corrected JML reduces repeated work in probability evaluations, with
  unchanged profile equations, numerical tolerances and handling of nonfinite
  values.

* Corrected JML checks for singleton assignment rosters before differentiating
  its covariance equations. Fixed-roster covariance remains unavailable in
  that case, with the same reason; point estimation and random-roster
  covariance rules are unchanged.

* Adaptive integration reuses cumulative category thresholds while evaluating
  each Person at a fixed parameter vector. This reduces repeated work in
  fitting, information calculations and conditional scoring without changing
  integration rules, numerical tolerances or source-readiness checks.

* `predict_mfrm_units()` admits conditional scoring for full-rank RSM/PCM
  latent regressions after source-consistency, local-solution and integration
  checks. The fitted predictor coding and residual population variance are
  retained. Conditional scores exclude calibration/population estimation
  uncertainty; global inference readiness is unchanged. GPCM regression still
  requires review.

* `import_tam_fit()` now identifies all-low and all-high observed responses
  from the fitted category support, preserving source EAPs and posterior SDs.
  Previously every Person was marked nonextreme. Missing responses and
  response masks are respected; entirely missing or unverifiable rows retain
  unavailable flags. Re-import the source TAM fit to update older imports.

* `fit_mfrm()`, `describe_mfrm_data()` and `review_mfrm_anchors()` accept
  `repeated_ratings = "independent"` to record an explicit assumption that
  repeated Person-by-facet rows are distinct, conditionally independent rating
  events. This removes only the duplicate-cell input-review reason; it keeps
  every rating, the same likelihood and all other source/inference checks.
  Default duplicate review is unchanged. Preparation notes, saved fits,
  replay settings and fitted-object scoring retain the declaration. Corrected
  JML accepts the same explicit declaration under its existing independence
  assumption. Two-family GPCM retains its existing input restrictions.
  The warning no longer incorrectly states that the model mathematically
  requires one observation per cell. Independence is an analyst assumption,
  not something verified by an event ID or this option.

* Two-family GPCM interval checks now ignore cached sparse QR factors when
  comparing the saved model with its reconstruction. Clearing or changing
  that numerical cache no longer invalidates an unchanged model; design
  entries, constraints and model metadata still require exact agreement.
  EM result assembly also permits eigensolver roundoff in a saved normal
  quadrature rule, with a relative check on every weight including the tails.
  The supplied nodes and weights, likelihood and score checks are preserved.

* Multivariate D-study comparison plots now fit standard landscape graphics
  devices, including recorded plots redrawn at a different size.
  More compact margin spacing preserves both panels, labels, reference points
  and intervals without changing the plotted values.

* Saved two-family GPCM results now load their required sparse-matrix methods
  when the package is attached in a fresh R session. Reports and score exports
  no longer depend on Matrix having been loaded by an earlier fit or another
  package; replay still performs no fitting or scoring.

* Latent-regression scoring now preserves the fitted basis for `scale()`,
  `poly()` and prediction-aware spline terms. Previously the same Person's
  conditional prior could change with the other Persons in the scoring batch,
  and a single-Person batch could fail. Saved fits retain the training terms;
  older fits reconstruct them from matching retained training data or stop
  explicitly when that basis cannot be recovered. Scoring batches need not
  have a full-rank covariate matrix because they do not re-estimate coefficients.
  Estimation, population-parameter inference and source-readiness rules are
  unchanged.

* `fit_mfrm()` help clarifies that extreme Persons are retained: ordinary
  JML distinguishes infinite Person estimates from optimizer traces, corrected
  JML uses profile-boundary contributions, and MML integrates their responses
  and reports posterior EAPs. The `has_exclusions` Person-boundary label does
  not mean input deletion. New-Person posterior scoring is distinguished from
  calibration. This is a documentation clarification, with no estimator change.

* `predict_mfrm_units()` now checks EAP and posterior-SD integration for every
  native scoring batch, including ordinary RSM/PCM with the retained prior.
  Previously that route could return a coarse-grid score even when supplying
  the identical prior explicitly required a finer grid. Failed checks stop
  under the default policy; increasing `scoring_quad_points` changes only
  scoring, and explicit review output retains the failure labels. New saved
  scores must retain their numerical checks for summary/export. Legacy scores
  keep their recorded procedure. Calibration and posterior interval definitions
  are unchanged; these intervals still exclude calibration uncertainty.

* `mfrm_facet_intervals()` adds native-scale model-information Wald intervals for
  estimated normal RSM/PCM populations, including full-rank conditional-normal
  latent regression, and intercept-only normal one-family GPCM MML, including
  the latter's explicit fixed-N(0,1) restriction. Direct fixed/adaptive fitting,
  unit weights and centered additive facets/steps are required. Full joint
  information and same-point finer-grid checks retain nuisance uncertainty;
  numerical refusals preserve points with unavailable bounds. Results, plots,
  reports and exports retain these checks and the native scale. The fixed-normal
  GPCM branch now recognizes its actual saved population metadata instead of
  rejecting a valid fit. Regression output retains fitted coefficients and
  residual variance, with their full cross-information included in structural
  uncertainty. These are not population-coefficient or variance intervals.
  Finite-sample coverage is not guaranteed. Use `scale = "standardized"` for
  explicit zero-sum contrasts with an intercept-only normal population;
  sandwich extensions and formal JML structural intervals remain unavailable.

* `mfrm_population_intervals()` reports experimental pointwise coefficient
  and residual-variance intervals for estimated normal RSM/PCM populations.
  It uses the full joint information inverse; variance bounds are formed on
  the log scale. Saved native-location intervals can supply the already
  qualified covariance without recomputing it. Source checks, numerical
  refusals and unavailable targets persist through results, reports and
  export/reopen. Population summaries remain point-only; GPCM population
  intervals and repeated-sample coverage are not qualified by this addition.

* Two-family GPCM fitted-curve plots now state once that calibration intervals
  are unavailable, instead of covering every curve point with a cross. Facet
  labels use separate lines so rating contexts are easier to read. All saved
  values, contexts and unavailable bounds remain unchanged; one-family curves
  still mark individual unavailable intervals and preserve gaps in ribbons.

* Help and output guides now consistently describe experimental two-family
  GPCM Person scoring, saved-score attachments, location/contrast intervals
  and posterior residuals. Portable-calibration guides distinguish corrected
  JML reference-prior scoring from two-family MML's fixed N(0,1) prior, which
  cannot be overridden. Obsolete statements that these routes were
  unavailable have been corrected. Conditional scoring still excludes
  calibration uncertainty, and interval coverage remains unqualified.

* Saved response diagnostics now distinguish complete, partial and unavailable
  residual output in `mfrm_results()`. The `response_overview` table reports
  selected, available, unresolved, missing-score, zero-variance and unselected
  source rows. Two-family GPCM, testlet and shared-rater reports show these counts
  on their first screen; incomplete numerical results no longer receive an
  unconditional available label.
  Saved estimates and the numerical acceptance checks are unchanged.
  Older saved `mfrm_results` objects without the overview recover it from their
  retained diagnostic rows when summarized, reported, exported or opened in a
  supported results viewer, without re-estimation.

* `mfrm_facet_intervals()` now supports experimental two-family GPCM MML
  location estimates and within-facet contrasts with `method = "model"`.
  It transforms the inverse full marginal information under the fitted location
  constraints, retaining slope/step nuisance uncertainty. Source, stationarity,
  information and quadrature checks must pass; otherwise points and reasons
  remain with missing intervals. Plots, result attachments, reports and exports
  preserve the fixed N(0,1) scale and experimental status. Sampling coverage is
  unqualified, and location differences do not imply uniform rating differences
  when slopes or steps differ. No sandwich, step/curve or rater-quality inference
  is added; individual rater sheets continue to display slope intervals only.

* Multivariate D-study printing now separates observed source counts and row
  usage from complete future plans and displays composite weights. Saved paired
  plan comparisons retain this source context without raw ratings or identifiers;
  older comparisons explicitly report when it is absent. Comparison plots now
  default to G/Phi or SEM on their original scale for every plan, including a
  marked reference, connected in supplied plan order. Use `view = "differences"`
  for the previous paired difference intervals. Both views identify future
  complete crossed plans; point projections do not acquire confidence intervals
  from the difference intervals. Coefficients, intervals and `summary()`
  tables are unchanged. Save full objects with `saveRDS()`;
  a coefficient CSV alone does not retain the design and interpretation context.

* `mfrm_report(..., style = "rater")` now creates experimental individual
  feedback sheets for two-family GPCM MML. The sheet separates location,
  component slope, step offsets, category use and assignment exposure, and
  reuses matching saved slope intervals and posterior residuals when supplied.
  Unresolved and unselected rows remain counted; location/step intervals and
  diagnostic cutoffs remain unavailable. Object, tables, Markdown and standalone
  HTML omit source identifiers and free-form source notes. No fitting or new
  inference occurs, and no rater-quality or training-effect claim is added.

* Saved two-family scores now reject inconsistent interval levels or algorithm
  labels before summary and result collection. Result attachments also compare
  their saved calibration values and coding with the supplied fit, so a retained
  source identity cannot make a changed calibration match. These checks require
  no fitting, integration or rescoring.

* Two-family `mfrm_results(fit, scores = scores)` now accepts saved native
  new-Person predictions and format-6 portable scores with matching source
  identity. Results, reports, CSV and replay retain unrounded EAP/SD/intervals,
  source/batch checks, omissions and review-only labels without rescoring.
  Older outputs without source identity remain standalone-readable but need
  regeneration before attachment; older portable calibrations must also be
  re-extracted. No refitting is required. These conditional intervals exclude
  calibration uncertainty and do not qualify coverage or population transport.

* The existing portable calibration workflow now supports experimental
  two-family GPCM MML in file format 6. It preserves ordered slope components,
  location/step constraints, category codes, fixed N(0,1) prior and source
  checks without training responses or Person estimates. Saved artifacts score
  new batches in a fresh session with separate integration checks. Unobserved
  combinations of known levels are flagged; prior overrides and repeated-event
  extensions are unsupported. Formats 1--5 retain their original meaning.
  Validation/freezing does not qualify coverage or population transport, and
  posterior intervals still exclude calibration uncertainty.

* `predict_mfrm_units()` now supports experimental conditional new-Person
  EAP, posterior SD and continuous posterior intervals for two-family GPCM
  MML. Both slope owners and original category steps are retained, with known
  levels, unit weights and the fixed N(0,1) prior. Source identity/convergence,
  finer source integration and each scoring batch are checked separately.
  Review-only results remain flagged; invalid calibrations cannot be enabled
  by review. Save/reopen and summaries retain the calibration and checks.
  Intervals exclude calibration uncertainty; coverage and population transport
  remain unqualified. Portable calibration follows the separately validated
  format 6 route above.

* Two-family results now collect all fitted non-Person locations, category
  steps and component slopes for tables and CSV export. Reports include
  location and step values alongside slopes, preserving facet identities,
  centering references and unavailable intervals without recalculation.
  Category steps are identified as centered offsets, not standalone thresholds.

* Two-family result summaries no longer recommend unsupported ordinary
  diagnostics or flag a Wright map as a missing required figure. `plot(res)`
  uses the sole available saved curve, slope-interval or response-diagnostic
  plot; multiple saved plots require an explicit `type`. Missing plots are
  explained without calculating new results, and unavailable intervals stay
  visible in the selected saved display.

* `mfrm_results(fit)` now accepts two-family GPCM fits with its default
  arguments, collecting saved estimates and numerical status without new
  fitting or diagnostic calculations. Saved curves, intervals and response
  diagnostics can still be attached explicitly; unsupported section requests
  remain errors. Replay instructions preserve the complete saved result even
  when no additional inference is attached.

* Adaptive GPCM MML can now use the existing bounded curvature restart when
  ordinary optimizer polishing leaves a failed gradient review. This uses
  the same adaptive likelihood with a fixed population and at most 64 free
  parameters. A proposal needs positive, well-conditioned curvature, a smaller
  gradient and no objective deterioration beyond roundoff; failures remain in
  the stage history. Information and interval checks are unchanged.
  Weak-information refusals now show the curvature-change, inverse-residual
  and standardized-displacement values with their limits.

* Two-family adaptive direct MML now compares neutral and same-data EM-derived
  starts using the same adaptive likelihood. This repairs an observed failure
  where a neutral start reached a poor, nearly flat solution despite a small
  gradient. `gpcm_mml_start = "neutral"` retains the earlier procedure.
  Selection preserves a better unfinished candidate's convergence status;
  starts, failures, seed and optimizer histories and costs remain saved.
  Summaries, results, reports and quadrature refits retain the initialization
  policy. Interval checks remain separate and sampling coverage is unqualified.

* `mfrm_response_diagnostics()` now accepts two-family adaptive direct MML
  fits. Source likelihood and gradient checks use the adaptive objective;
  prediction integrates each Person's complete observed record with both
  slope families and calibration fixed. Higher-order row checks preserve
  unresolved probabilities and groups. Selecting output rows does not shorten
  the conditioning record. Existing plots, reports, exports and saved replay
  retain the integration method and descriptive interpretation, without
  expectation-one references, fit cutoffs or calibration-uncertainty propagation.

* Two-family adaptive direct MML now supports experimental component log-Wald
  intervals through `confint(fit, parm = "slopes")`. Full observed information differentiates
  the moving-node likelihood, including all nuisance parameters. Fresh source,
  convergence, score-derivative, stationarity and higher-order adaptive-integration
  checks retain estimates with missing bounds when unsuccessful. Saved interval
  settings, plots, reports and exports retain the actual integration method.
  These numerical checks do not establish sampling coverage. Adaptive profile
  intervals remain unavailable.

* Two-family GPCM now permits `mml_engine = "direct"` with
  `mml_integration = "adaptive"` in `fit_mfrm()`. This uses the existing
  moving-node likelihood and gradient while retaining the fixed N(0,1)
  population and both slope owners. Summaries, conditional curves, saved
  reports and quadrature-order refits retain the actual engine and integration
  method. Fixed-grid EM remains available. Adaptive two-family profile
  intervals are unavailable; numerical
  convergence does not establish identification or interval coverage.

* `mml_quadrature_sensitivity()` and its GPCM alias now accept
  `adaptive_quad_points` for two slope families. The existing adaptive
  integration calculation uses both slope components and compares integrals
  at unchanged calibration parameters. Per-Person results distinguish
  integration error from changes introduced by refitting. Summaries and saved
  results retain these diagnostic quantities without adding Person scoring
  or adaptive EM; the fixed-calibration review alone does not change interval eligibility.

* The GPCM guide adds a judged-sport example based on official Olympic score
  protocols. It distinguishes competition totals from latent ability, preserves
  judge identities and advancement, and explains the additional uncertainty
  needed for winner selection. GPCM and multivariate D-study help clarify that
  individual intervals and reliability coefficients do not supply winner
  probabilities or simultaneous rank confidence sets. No new ranking estimator
  or changes to fitted scores are introduced.

* The GPCM guide adds an empirical writing-data workflow using the separately
  available `sirt::data.ratings1` dataset. It explains uneven assignment,
  unused factor levels, absent categories and unknown planned assignments,
  and shows why numerical convergence alone is insufficient for rater
  feedback. It also distinguishes the supported model from published OSCE
  extensions and clarifies experimental two-family interval availability.

* Corrected-JML summaries, reports and exports now explain which repeated
  sampling the reported `RootSE` describes. Fixed assignment counts still
  allow the composition of a new cohort to vary. This uncertainty does not
  describe reassessing the same Persons at fixed abilities. Help clarifies
  the sampling assumptions and small-sample limits; estimates and covariance
  calculations are unchanged.

* Adaptive MML fitting and calibration checks perform fewer repeated matrix
  calculations for RSM, PCM and GPCM. The likelihood, integration rules and
  numerical accuracy requirements are unchanged.

* Numerical integration sensitivity checks for RSM and PCM no longer emit
  an unnecessary warning about a missing slope column.

* `recommend_mfrm_design()` avoids building unrelated planning tables when
  selecting a design. Recommendations still use unrounded results and retain
  failed replications in the convergence and workload checks.

* Two-family GPCM MML--EM now supports an explicitly selected component profile
  through `confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t1"))`, using
  the fitted facet and level names. Other slopes, locations and steps are
  reoptimized while the population and identification constraints remain fixed.
  Saved results retain owner identity, numerical checks, unresolved endpoints,
  the same-target Wald comparison and profile plots/reports/exports. This is
  an experimental local likelihood-ratio interval; better coverage is not
  established. No default estimator or interval method changes.
  Profile calculations are faster and more reliable when starting estimates
  differ substantially. After an inaccurate outward trial, the search can
  find an endpoint closer to the fitted slope, while retaining the failed
  trial in tables and plots. Failures inside the interval search still leave
  the affected endpoint unavailable. Numerical accuracy requirements are
  unchanged; profiles remain optional sensitivity analyses, not a routine
  replacement for Wald intervals.
  A failed optimization can also retry with rescaled search coordinates when
  a small slope makes other coefficients difficult to estimate. The same
  model, parameter constraints and numerical checks apply; unsuccessful
  searches still leave the affected endpoint unavailable.

* Two-family MML--EM more reliably reaches the accuracy requested by
  `em_score_tol` when an M step stops too early. Accepted steps still preserve
  likelihood ascent, and `fit$opt$em_trace` records convergence and numerical
  difficulties. Improved numerical convergence does not guarantee reliable
  confidence intervals, particularly when slopes are poorly determined.

* Two-family GPCM MML--EM fits now connect to `mfrm_response_diagnostics()`.
  The existing response model integrates each Person's ability posterior while
  holding both slope families, locations and steps fixed. Original observed
  assignments and row identities are retained; selecting rows does not shorten
  the conditioning record. Category probabilities, expected scores, mixture
  variances and descriptive Infit/Outfit receive row-level integration checks.
  Unresolved rows remain visible and are not silently excluded from groups.
  Saved results support paired/scatter plots, ggplot conversion, reports and
  exports without refitting. There is no expectation-one reference, automatic
  fit cutoff, formal fit test or calibration-uncertainty propagation.

* Corrected shared-owner GPCM JML now connects to `predict_mfrm_units()` and
  the portable calibration lifecycle. New-Person EAP uses the corrected point
  calibration and a separate normal reference prior. File format 5 preserves
  the explicit correction order, assignment-sampling assumption and fresh
  adjusted-equation/root checks, without training responses or Person
  estimates. It does not apply ordinary JML likelihood checks to a corrected
  root. Covariance failure does not suppress a valid calibration. Score tables,
  summaries and portable score plots retain the experimental estimator and
  conditional uncertainty; residual calibration bias may remain. Structural
  intervals and corrected Person ML/WLE are not added. Formats 1–4 retain their
  existing meaning, and every score batch still receives integration checks.

* `mfrm_response_diagnostics()` now supports corrected GPCM JML. It uses
  saved corrected calibration and reprofiled Person abilities to compute
  conditional probabilities, expected scores, variances and descriptive
  residuals on observed rows. It performs no posterior integration and
  excludes calibration/Person estimation uncertainty. Extreme-score
  probabilities remain available when zero variance makes standardization
  undefined. Infit is retained when its aggregate ratio is defined, even if
  Outfit is unavailable; no rows are silently dropped. Paired/scatter plots,
  ggplot conversion, results, reports and saved exports preserve this
  distinction without cutoffs or rater-quality flags. Saved attachments
  must match the slopes, Person profiles and exact fitted roster.

* `fit_mfrm()` adds an explicit, experimental `jml_correction_order` for
  shared-owner GPCM JML with centered fixed facets and observed unit-weight
  ratings. It solves a profile-score-adjusted equation and retains numerical
  attempts, unresolved roots, point estimates and covariance failures
  separately. `RootSE` describes local adjusted-equation variation; residual
  structural bias may remain, and no confidence limits are supplied. Summary,
  point/distribution plots, ggplot conversion, results, reports and saved
  exports preserve this interpretation. Unsupported ordinary diagnostics,
  likelihood comparison are refused; conditional post-hoc EAP and portable
  scoring use the distinct corrected-equation checks described above. The default MML
  and uncorrected JML calls are unchanged; correction order is not selected
  automatically. Help and the scope guide explain the sampling assumptions,
  supported assignments and remaining inferential limits.

* `mml_quadrature_sensitivity()` and its GPCM entry point now support two
  slope families with fixed-grid MML--EM. Refits preserve owner order,
  the score ladder and EM controls; comparisons use both slopes in category
  probabilities and full joint information for raw diagnostic SEs. Saved
  component intervals retain their numerical checks and connect to existing
  plots and tables. Person-score comparisons are explicitly unavailable. Adaptive
  direct refits and fixed-calibration integration comparisons are now supported
  as described above. No fixed
  number of integration points is sufficient for every dataset; compare
  estimates and intervals across grids before interpreting sensitive results.

* Two-family slope intervals now retain a warning for a small optimization
  residual instead of automatically withholding every interval. The
  standardized Newton displacement must be at most 0.01; values above 1e-4
  carry the additional caution and remain visible in saved numerical checks,
  plots and reports. The fit, estimates and EM stopping rule are unchanged.
  Mean-score, information, derivative and quadrature checks still apply.
  This numerical approximation does not establish sampling coverage.
  Sparse assignments can leave slopes poorly determined even when numerical
  checks pass. Bias and underestimated uncertainty remain possible.

* `fit_mfrm()` now accepts an ordered pair of slope facets for a provisional
  two-family GPCM MML--EM route. It requires explicit fixed-standard-normal
  identification, exactly two facets, second-owner steps, unit weights and no
  anchors or population covariates. The first family's slopes have geometric
  mean one; the second is free on the fixed ability scale. Original facet
  names and owner-level identities are retained. `em_score_tol` controls the
  marginal-score stopping rule per Person; unsupported options are rejected.
  Summary and conditional probability/information curves connect to plots,
  saved results, reports and exports. These curves have no calibration
  intervals. Ordinary diagnostics, Wright/Pathway maps, other parameter intervals,
  and model ranking remain unavailable. Conditional new-Person and portable
  scoring are provided by the later extensions described above. The
  capability table and workflow guide distinguish this scope from one-family
  GPCM. Existing single-family defaults and EM-to-direct fallback are unchanged.

* `confint()` now offers explicitly experimental component-slope intervals for
  the two-family fixed-standard-normal MML--EM model. Full joint marginal
  information includes nuisance parameters and cross-family covariance.
  Source, score-derivative, stationarity and higher-order quadrature checks
  are retained through plots, tables, reports and exports. Failed checks keep
  estimates with missing bounds and reasons. Both families retain their own
  scale reference; no unit-slope p-values or automatic rater-quality decisions
  are supplied. The default scale for two families is `"standardized"` (SD fixed
  at one). One-family defaults are unchanged. These numerical checks do not
  establish global identification, boundary absence or sampling coverage.

* Summary help now distinguishes slope owners, primary estimates and numerical
  optimizer values, and explains why a scale constraint alone does not establish
  inferential reliability. It identifies the fields for the actual MML engine,
  iteration count and numerical stopping rule. Corrected JML now has the
  explicit experimental workflow above; its formal structural intervals remain
  unavailable.

* `as_ggplot()` now converts residual-PCA scree, reference, excess and loading
  plots explicitly. Previously, reference plots could become unrelated bars
  through generic column guessing. Conversions retain the saved numerical
  values, styles and replay data; incomplete payloads fail clearly. Base
  reference plots now show a legend. Summaries and reference plots explain
  that scanning componentwise cutoffs does not control the chance of any false
  flag, and bootstrap status distinguishes unavailable observed PCA from a
  failed refit.

* Residual PCA from a fitted object now skips unrelated marginal-fit tables,
  reducing repeated diagnostic work during the model bootstrap while preserving
  observation residuals and PCA results. The visual guide now distinguishes
  shared-Person links from fixed parameter anchors, explains why a connected
  sparse design can still lack usable residual correlations, and describes
  when testlet and correlated-trait interpretations yield the same ratings.

* `analyze_residual_pca()` adds `parallel_method = "model_bootstrap"` for
  RSM/PCM MML fits with fixed standard-normal population, fixed quadrature,
  additive facets, unit weights and no anchors or shrinkage. Ratings are
  regenerated on the retained assignment and the model is refitted once per
  replicate for all requested PCA scopes. Every attempt and warning is saved;
  an affected comparison is unavailable if any replicate fails. Existing
  summary, scree/excess plots and ggplot conversion distinguish this reference
  from the unchanged default residual permutation. Both remain exploratory,
  without calibrated dimensionality decisions or multiplicity adjustment.

* The visual diagnostics guide and residual-PCA/Q3 help now explain how to
  examine a one-ability assumption before fitting a multidimensional model.
  They distinguish current exploratory screens from model-generated reference
  tests and residual-network/EGA methods, clarify sparse-overlap and
  community-count interpretation, and identify external comparison requirements.
  No calibrated dimensionality test or residual-network estimator is added.

* Testlet result/report exports now retain the complete saved scoring roster,
  including explicitly omitted responses and original column names, alongside
  observed block counts. Testlet scores and random-rater probability predictions
  also record effect-sharing and integration meanings in their settings;
  observed-rater conditional distributions remain distinct from new-rater
  population predictions. Older saved scores without a roster remain readable.

* GPCM, testlet and random-rater help now distinguish discrimination from local
  dependence and identify which ratings share a latent effect. The guide
  explains performance versus rating-occasion blocks, observed/new-effect
  prediction and why local variance does not by itself establish halo or
  equal task weights. Existing model availability is unchanged.

* The GPCM guide and workflow help now explain proposed two-slope applications
  for task review and assessor feedback, with speaking, music and clinical
  assessment examples. They distinguish response models from MML/JML estimation
  and EM/direct optimization, and explain why the former "bounded GPCM" label
  should be replaced by explicit structural and output restrictions. The
  provisional two-family workflow above has its own narrower output scope.

* The GPCM guide clarifies ability-SD changes, scale identification and the
  difference between common rated performances and fixed parameter anchors.
  Multivariate D-study plot help explains when to use 2D comparisons, why SEM
  is not a confidence band, and why a higher G/Phi in a more heterogeneous
  population need not mean lower measurement error. Existing plot methods
  and model availability are unchanged.

* GPCM predictions and probability-based diagnostics now share the conditional
  probability calculation used by likelihood evaluation. Full category log
  probabilities are retained directly, including when an extremely small
  probability rounds to zero.
  Curve interval results now retain a `contexts` table identifying predictions
  for combinations of known facet levels not observed together in the retained
  rating data. Printed output, default plot captions and report tables preserve
  this distinction; it does not certify interval reliability.

* The GPCM guide explains how MML--EM can estimate task and rater slopes
  together and why the response equations of `sirt::rm.facets()` and the
  Uto--Ueno model generally differ. The provisional two-family `fit_mfrm()`
  route uses the whole-predictor slope product; this does not establish
  free-estimator equivalence with sirt.
  For the two-family formulation, it distinguishes the task-family geometric
  reference from free rater slopes: a rater slope of one does not mean the
  average rater or provide a competence threshold.
  Facet names are distinguished from model roles: example labels such as
  Task/Rater are not mandatory column names, and labels alone do not determine
  the slope reference or step owner.
  The guide also connects Wang and Liu's latent-regression formulation to
  `population_formula` and distinguishes ability-dependent response sensitivity,
  information and external scoring accuracy when interpreting rater feedback.
  It explains why complete task-rater crossings alone do not establish
  identification of two slope families when persons contribute too few ratings.
  The GPCM and workflow guides also distinguish multivariate observed-score
  G/D studies from multiple GPCM slopes, with literature on complementary and
  joint IRT/G-theory approaches. No automatic conversion or joint fitter is added.

* The linking guide now starts with assignment checks before fitting and
  distinguishes common-person links from fixed parameter anchors. Runnable
  examples cover common persons, rotating subsets and random subsets, with
  task-specific skeletons explained separately. The workflow guide clarifies
  that MML is a default rather than an automatically selected estimator.
  The GPCM explanation uses separate slope-owner and step-owner indices to
  match the supported MML model. Estimation behavior is unchanged.

* `analyze_facet_equivalence()` now marks fixed locations explicitly and omits
  their marginal intervals. Contrasts with estimated levels and deviations
  from the estimated facet mean retain their covariance-based uncertainty.
  Matching diagnostics with fixed-value labels now give the same results as
  the direct fit route. This also fixes an MML path that could still display
  a zero-width interval for an anchor after the general diagnostic repair.

* `fit_measures_table()` now carries SE methods and interval interpretation
  through its report tables, summaries and table exports. JML normal bands
  retain their exploratory status; missing provenance in older inputs is
  explicit. The `measure_ci` plot preserves these labels when changing the
  interval level and displays finite fixed/no-interval estimates instead of
  dropping them. Captions can be hidden with `show_notes = FALSE` while their
  explanation remains in saved plot data. Estimation and SE formulas are unchanged.

* Diagnostic tables now identify location values fixed by anchors or
  identification constraints. Such values no longer receive exploratory JML
  SEs or zero-width MML confidence intervals. Sampling SEs and intervals are
  marked not applicable, including attached diagnostics, reports and Wright
  maps. Facet reliability/separation and variability tests involving fixed
  values are unavailable; free-level SEs and response-fit statistics remain
  available. Supplied-anchor uncertainty is not estimated.

* Category variance now uses centered moments in information curves, curve
  intervals, response diagnostics, posterior category summaries and the
  shared-rater/testlet derivative calculations. The R and compiled posterior
  paths agree. This fixes cancellation that could turn small positive upper-tail
  information into zero and disagree with summed category information.
  Existing diagnostic variance thresholds and inference eligibility rules remain
  unchanged; this numerical repair is not a new interval-coverage guarantee.

* JML estimation help corrects the displayed facet/location SE formula to
  include observation weights and squared GPCM slopes, matching the existing
  calculation. It distinguishes exploratory SEs, local curvature checks for
  portable scoring, and formal inference, with a binary Rasch reference for
  the incidental-parameter limitation. Estimates and SE calculations are
  unchanged; formal JML slope intervals remain unavailable. The help also
  explains why applying a common item-bias multiplier to all relative GPCM
  slopes cannot correct their ratios.

* Calibration and prediction help now distinguish the literature basis for
  fixed-calibration EAP from checks of numerical accuracy and
  JML estimation. Fixed-score agreement with external software does not imply
  matching free estimators, calibration bias or interval coverage. The
  `quadrature_review` help now specifies that its RSM/PCM requirement applies
  to MML, not JML.

* Shared-owner GPCM JML now supports conditional fitted-object and portable
  EAP scoring with a post-hoc reference prior (file format 4). Fresh joint
  likelihood, gradient and unregularized full-curvature checks are separate
  from the retained global audits. Known Person/additive/slope boundary
  certificates prevent automatic scoring even when optimizer traces and local
  curvature are finite. Incomplete global audits stay labelled as incomplete
  through saved artifacts and score summaries. Unit weights and no anchors or
  interactions are required. This does not supply formal JML slope intervals,
  global-optimum or calibration-bias guarantees. Frozen scoring reuses the
  existing kernel and does not repeat source curvature calculations.

* RSM/PCM JML now supports the existing portable calibration lifecycle, using
  file format 3 and a post-hoc standard-normal reference prior. Extraction
  requires a finite identified source, unit weights, no anchors/interactions,
  and passing fresh joint-likelihood/gradient checks. Each new batch checks
  EAP/SD integration; summaries retain source evidence and distinguish the
  reference prior from a distribution estimated by JML. Explicit scoring-prior
  sensitivity remains available. JML calibration bias and uncertainty are not
  included in posterior intervals. Formal JML slope inference remains
  unavailable. Existing MML file formats are unchanged.

* Portable score summaries distinguish checks on the reported EAP/SD from
  comparisons of alternative integration grids. With adaptive scoring, the
  fixed grid in the comparison is explicitly identified as unused for the
  reported scores. The GPCM tutorial now connects calibration, a fresh scoring
  session, prior sensitivity and saved score/review tables, including people
  without valid responses. Score plots state how many people have no score;
  their zero reference is labelled as the scale origin, rather than assuming
  it equals an estimated or explicitly changed prior mean.

* `confint()` now offers experimental `method = "profile"` for one explicitly
  named GPCM relative `slope`. Nuisance calibration and population parameters
  are reoptimized under the sum-zero log-slope constraint. Two-start, gradient
  and integration checks preserve unresolved endpoints; a finite search limit
  is not labelled an infinite confidence bound. Shared and separate owners
  are supported within the estimated-normal, unit-weight scope without anchors
  or interactions. Saved curves, Wald comparisons and numerical records connect
  to `plot(..., type = "profile")`, ggplot conversion, tables and reports without
  refitting. Model/sandwich defaults are unchanged. Profile coverage, global
  optimality and superiority over Wald are not established. A small-sample
  saved-fit comparison found substantially lower profile interval availability;
  the method remains for explicitly requested sensitivity analysis, not a
  routine replacement for Wald intervals.

* Portable GPCM MML calibration now uses the existing extract/validate/freeze/
  save/load/score API. It retains both slope/step owners, one relative-slope
  family, an estimated intercept-only normal prior and passing conditional
  source-check evidence. Unit weights and no anchors/interactions are required.
  Shared and separate owners are supported. Replay validates the stored source
  evidence without refitting; each batch must pass numerical scoring checks.
  Failed integration stops instead of returning review-only numerical scores.
  `score_mfrm_calibration()` also accepts an explicit common normal scoring
  prior for current artifacts, preserving original and actual priors through
  summaries and plot data. Existing RSM/PCM default artifact semantics remain
  unchanged. These are conditional posterior scores, not a global-maximum,
  population-transport or interval-coverage guarantee.

* `predict_mfrm_units()` and `sample_mfrm_plausible_values()` now accept
  `scoring_prior = list(mean = ..., sd = ...)` for sensitivity to a common
  normal scoring prior. For MML the default retains the fitted prior; JML uses
  a post-hoc standard-normal reference prior. The calibration
  stays fixed, its source checks remain in force, and numerical checks use the
  specified prior. Original and scoring prior values remain unrounded through
  estimates, summaries, draws and exports. Background-covariate population
  overrides are refused. Help includes a plot comparing identical response
  patterns under two priors; this does not establish a new cohort's distribution
  or estimate its distribution.

* Fitted-object scoring now distinguishes local calibration eligibility from
  global inference readiness. GPCM MML and intercept-only normal population
  models can return conditional scores after fresh likelihood, gradient,
  unregularized-information and integration checks. Per-person checks compare
  reported EAP/SD with adaptive references; unresolved scores require explicit
  review and retain their labels through summaries, draws and exports. Invalid
  calibrations or inconsistent stored priors cannot be enabled by review.
  Unevaluated boundary audits are labelled as unevaluated. These changes do not
  validate a prior for a new cohort.
  The scoring help explains the retained-prior default for a later cohort,
  distinguishes prior sensitivity from recalibration, and clarifies why a
  narrower interval from a narrower prior is not additional rating evidence.

* Adaptive MML fitting now rejects numerically unevaluable trial values,
  including overflow with finite GPCM slopes, so the optimizer can shorten
  its step. Failed evaluations cannot reuse an earlier likelihood from the
  cache. Invalid starting values still fail explicitly, and final convergence
  checks remain unchanged; this does not impose an upper discrimination limit
  or establish that a fitted model has a finite global maximum.

* Ordinary fitted-model plots now accept named `level` as an alternative to
  `ci_level`, also through results and ggplot entry points. Previously `level`
  could be silently ignored. Supplying both names is an error; the default
  remains 0.95 and old positional calls retain their meanings. Saved interval
  results retain their original confidence level.

* Conversion of saved core plot payloads now rejects unsupported settings
  instead of silently ignoring them, including attempted title or interval
  changes. All saved plot-data views also reject a new `type`, including
  extended-model and D-study converters. Create a new payload from the fit/result to change plot settings,
  or use ggplot label controls after conversion. Complete CCC views retain
  their documented conversion-specific appearance options.

* The reporting guide now starts from saved results, figures, reports and
  analyst archives, with specialist manuscript tools described separately.
  The fixed-rater tutorial explains report-purpose and output-format defaults,
  explicitly selects a saved interval for recipient sheets, and demonstrates
  exporting and reopening an analysis without recalculating its intervals.
  Export presets select files; plotting presets select appearance.

* GPCM MML now allows different slope and step facets, for example criterion
  discrimination with rater-specific category steps. One slope family retains
  geometric-mean-one identification and multiplies the entire adjacent-category
  predictor. Fitted-object scoring, information, eligible slope/curve intervals,
  matched PCM comparison, same-design bootstrap, and saved results preserve both
  owners. CCC and expected-score pathway plots show labeled step/slope pairs;
  single-facet fit flags are not assigned to those joint profiles. JML, weighting
  reviews and simulation/design workflows still require a shared owner. This
  extension does not establish finite-sample coverage or provide simultaneous
  criterion/rater slope families.

* Random-rater interval help now explains why increasing bootstrap size alone
  cannot resolve unbounded limits: at 95% with 499 planned draws, 13 unresolved
  studentized errors for a rater make both endpoints infinite. This describes
  the existing failure-preserving rule, not a new acceptable-failure threshold
  or a coverage guarantee.

* Random-rater bootstrap trial records now retain the optimizer code,
  numerical/information checks, Person quadrature orders and discrepancies,
  and separate variance-boundary checks. Missing check values remain missing
  when a refit fails. These records carry through existing result tables,
  reports and exports, helping distinguish numerical problems from interval
  performance. Interval formulas, eligibility rules and saved older results
  are unchanged; additional historical checks cannot be reconstructed by
  reopening an old result.

* The feedback and external-feature tutorials now show how to save and reopen
  analyses and retain a reviewed HTML sheet. They distinguish the recipient's
  standalone sheet from analyst archives, and a single-completion PCA view
  from the full imputed analysis. Saved results reuse their original estimates
  and display choices; reopening does not update them for new data.

* `as_ggplot()` now converts pooled fixed-facet multiple-imputation intervals.
  It uses saved t-interval endpoints, preserving target order, contrasts,
  degrees of freedom and complete-data information cautions. Fixed targets
  have no inferential interval; missing intervals use crosses and infinite
  bounds use arrows. Default and explicit table conversion retain the full
  view, without refitting, repooling or recalculating confidence limits.

* `as_ggplot()` now converts silhouette and numeric/categorical feature-profile
  views. Silhouettes preserve negative widths, saved order and the overall
  mean line. Numeric profiles retain original-unit means, medians and counts,
  with separate symbols and small vertical offsets. Categorical profiles retain
  group counts, unused levels and original level order on a fixed proportion
  scale. Default and explicit table conversion preserve the complete view;
  source data remain extractable. No groups or intervals are estimated.

* `as_ggplot()` now converts imputation co-membership heatmaps with the saved
  ID order, all-imputation fractions and a fixed zero-to-one scale. Grey fill
  plus crosses distinguish unavailable pairs from zero. Default and explicit
  `component = "matrix"` retain the same view; imputation count, exclusions,
  selected IDs and source values remain available through `plot_data()`.
  Selecting IDs or hiding labels does not renormalize or sample values.

* `as_ggplot()` now converts saved external-feature dendrograms. Merge heights,
  leaf order, requested group boxes and label settings are retained, including
  tied-height cuts and excluded-ID records. The default and `component = "tree"`
  use the same complete view. No clustering or partition is recomputed;
  heights do not measure significance or branch support.

* `as_ggplot()` now converts external-feature PCA scree, scores and loadings
  views using saved components, group colours/shapes and label settings.
  Scree plots distinguish retained components, scores keep equal axis units,
  and loadings keep feature order and the zero reference. The default and
  `component = "table"` use the same dedicated view; transformation metadata
  and excluded IDs remain available through `plot_data()`. No PCA or grouping
  is refitted.

* `options(mfrmr.plot_preset = "publication")` now sets the session default
  for plots using the common preset controls. Explicit call settings take
  precedence; removing the option restores `"standard"`. Supported conversions
  of saved plot payloads retain their saved preset when the option changes.
  Plots with separate palette controls keep their own settings. Estimates,
  intervals, screening thresholds and the global ggplot theme are unchanged.

* Bubble-chart ggplot conversion now uses saved circle radii, facet colours
  and facet order instead of equal-sized points and a new default palette.
  The native bubble plot also honours the monochrome preset; explicit palette
  overrides remain supported. Relative radii and reference lines are retained,
  while physical sizes can differ between base graphics and ggplot.

* Eight plotting functions now accept `title` alongside the supported
  legacy `main` spelling: marginal fit/pairwise, unexpected responses,
  interrater agreement, facet chi-square, bubble, bias interaction and facet
  dashboard plots. `title = NULL` suppresses the heading; omission and legacy
  `main = NULL` retain the default. Supplying both names is an error. Existing
  positional calls, numerical data, screening settings and interpretation
  notes remain unchanged. See `?mfrmr_visual_diagnostics` for the exact routes.

* `mfrm_report(..., style = "rater")` now prepares an individual feedback
  sheet from saved native additive RSM/PCM results. HTML includes category-use
  bars, numerical tables and print styling; plain-language and researcher
  presentations share the same selected values. Missing diagnostics and
  individual intervals are explained without new calculations. All four
  output formats omit the source fit and source identifiers. Rater selection
  and ambiguous interval choices are explicit; fit statistics do not trigger
  automatic quality labels or exclusion.

* The visual-diagnostics tutorial now offers six linked figure previews with
  accessible labels, executed examples and data-extraction code. The gallery
  uses the plot guide for purposes and ggplot support, including a new subset
  coverage route. Coverage help distinguishes relative observed facet-level
  counts from completion of planned assignments. External-feature trees and
  composite D-study planning have compact examples with interpretation limits.

* `mfrmr_output_guide("plots")` now maps selected figure purposes to result
  classes, plot calls, reusable data and ggplot routes. It distinguishes
  dedicated converters, plots that already return ggplot, generic table views
  and unavailable conversions with alternatives. Existing guide scopes retain
  their table format. This adds guidance; estimators and renderers are unchanged.

* REML profile calculations reject invalid variance trials before evaluation
  and rebuild the numerical likelihood state after a failed factorization.
  A failed trial no longer corrupts later evaluations at valid parameter
  values; numerical failures remain distinct from accepted interval bounds.

## Getting started

* Results now reject unsupported prediction attachments before routing interval
  inputs. An incomplete fit object receives an explanatory input error instead
  of an internal missing-logical-value error. Valid fitted analyses and their
  numerical results are unchanged.

* GPCM bootstrap output now preserves recorded selected-refit history in
  interval cautions, printed output and reporting tables. An assembled
  reanalysis is no longer displayed as though it were a complete rerun under
  one procedure. Saved draws, interval values and acceptance decisions are
  unchanged. The GPCM guide distinguishes the evidence for each inference
  target and retains the adverse probability-coverage findings; bootstrap
  coverage and a general finite-sample guarantee remain unestablished.

* The GPCM guide now reports a reanalysis of 800 saved fits for standardized
  slope differences, category probabilities and per-rating information.
  Probability intervals showed undercoverage in small incomplete designs,
  including Bonferroni-adjusted grid families. Refitting all 100 datasets in
  one affected condition with the updated optimizer left every parameter and
  interval endpoint unchanged; these are matched reanalyses, not additional
  independent replications. Curve printing and default figure subtitles now
  explicitly identify approximate intervals. Help distinguishes numerical
  availability from nominal coverage and historical fitted estimates from a
  new evaluation of the current estimator. These findings do not change the
  interval formulas or remove the APIs.

* Ordinary result summaries and report decisions now use the saved diagnostic
  precision assessment instead of describing it as unevaluated. Reports also
  retain that assessment when a separate precision-review section was not
  requested. Missing assessments remain unreviewed, and a positive assessment
  cannot override restrictions on the source fit. No estimates, intervals or
  eligibility rules change.
  The introductory workflow now uses explicit rubric bounds and
  `category_policy = "preserve"`, demonstrates reopening saved results, explains
  session-dependent screening bands, and links rater differences, missing-score
  imputation, external-feature groups and observed-score G/D planning to their
  dedicated workflows.

* GPCM curve plots now mark estimates with unavailable intervals using crosses
  and report their count in a caption. Ribbons stop at unavailable grid points.
  Long numerical cautions wrap within the default figure. `caption = NULL`
  removes the caption while retaining markers and saved failure reasons.
  The plotting guide explains how to use a logarithmic axis for very wide
  positive slope/ratio intervals; interval calculations are unchanged.
  Markdown report summaries now retain GPCM interval counts (finite,
  unbounded and unavailable) and saved inference cautions/reasons, matching
  the information already retained in result tables.

* MML information calculations now distinguish poor conditioning from an
  unusable inverse. Positive ill-conditioned information receives two levels of
  numerical refinement and checks of curvature stability, unregularized
  inversion and the curvature-scaled gradient. A verified inverse can support
  GPCM intervals and local IC/bootstrap checks with an explicit caution;
  failed verification leaves inference unavailable. No eigenvalue floor is
  used for the accepted inverse. This common MML calculation also applies to
  RSM/PCM. Weak information can still produce extremely wide or unreliable
  intervals. Cautions follow slope/curve intervals, plots, tables and reports;
  bootstrap checks retain refinement diagnostics. The existing information
  workspace budget also checks the additional refinement allocation.

* Fixed-facet intervals, practical-equivalence analyses and pooled imputation
  results now retain weak-information cautions when a numerically verified
  inverse is used. Warnings remain with tables, saved results and plot data;
  fixed-facet reports also retain the numerical review. Pooling identifies the
  affected imputations and never discards a failed completion. These warnings
  do not change point estimates, covariance formulas or eligibility criteria.
  Pooled plots now accept optional `title` and `subtitle`, including `NULL`;
  hiding plot text does not remove the saved caution. Quadrature reviews also
  retain the covariance calculation's explanation alongside its status.

* Small fixed-grid GPCM MML fits now review negative numerical curvature even
  after a small terminal gradient. Up to three BFGS restarts can rescale the
  search coordinates to recover a better finite solution; failed recovery
  retains the original estimate with a numerical warning. The complete stage
  history records curvature and any review error. This does not regularize the
  information matrix. Retained solutions undergo separate interval,
  information-criterion and bootstrap checks. Quadrature accuracy remains a
  separate requirement.

* GPCM bootstrap checks now retain the population SD and the minimum/maximum
  population-SD-standardized optimizer slopes, without additional fitting or
  information calculations. These diagnostics help separate changes in scale,
  near-zero discrimination and category-support problems. Help explains why
  stable slopes do not establish finite thresholds or valid Wald intervals.
  APA check tables preserve small scale, slope, gradient and information
  values using significant digits, including in HTML/LaTeX table conversion,
  instead of displaying a small positive value as zero. Raw checks stay numeric.
  The additional diagnostics do not change replicate admission or certify
  a boundary solution.

* GPCM bootstrap results now distinguish refit errors from returned estimates
  withheld by eligibility checks. `trials` records the last stage and which
  model fits returned; `checks` retains category support, numerical status and
  already-computed information diagnostics. `refit_draws` saves returned
  alternative-model parameters for diagnosis, including rejected estimates;
  these are never substituted for accepted `draws` in confidence intervals.
  Saved checks connect to APA tables, reports and exports.
  A confirmed singleton-only category warning no longer excludes a GPCM
  point estimate from basic slope bootstrapping when the remaining identity,
  convergence and unregularized-information checks pass. The source fit uses
  the same rule. `BootstrapEligible` and `WaldEligible` remain separate;
  cautions persist in intervals, printed output, default plot subtitles and
  reports. Wald, IC and LRT rules are unchanged. This does not establish
  finite-sample coverage or admit singular-information estimates.

* Results with attached inference now show complete-result save/reload code
  in summaries and HTML. This prevents the displayed reproduction route from
  omitting saved RSM/PCM intervals, GPCM inference or posterior diagnostics.
  Exported replay continues to reload the saved RDS. Starter export indexes
  now show saved RSM/PCM and GPCM inference figures with descriptions, keeping
  their targets separate from ordinary Wright/Pathway uncertainty.

* Saved RSM/PCM fixed-facet intervals now connect to `apa_table()`,
  `mfrm_results(intervals = ...)`, reports and exports. Named results have
  `facet_` plot routes. Source matching checks the fitted data, parameters,
  constraints, population and integration settings; replay reloads saved
  results without refitting or recomputing covariance. Method, confidence
  level, contrast coefficients, cluster mapping and unavailable outcomes
  remain available. `as_ggplot()` supports the same interval comparisons,
  with line types and offsets as well as color. Titles, subtitles, captions,
  reference lines and legends are optional; noninteger confidence levels are
  displayed without rounding to a different level. The estimator and interval
  calculations are unchanged, and ordinary Wright/Pathway displays retain
  their own uncertainty and fit meanings.

* Fixed-grid GPCM response-pattern checks now reuse cumulative category
  probabilities within each rating design and parameter vector. All response
  patterns, posterior weights and information checks are retained; the point
  estimator and default interval calculation are unchanged. Adaptive nodes
  are excluded from this reuse.

* Saved GPCM slope/curve intervals and bootstrap results now connect to
  `apa_table()`, `mfrm_results(intervals = ...)`, reports and exports.
  `plot()`, `as_ggplot()` and `plot_data()` preserve the selected target,
  method, level and multiplicity adjustment. Plots retain unavailable and
  unbounded outcomes, use non-color cues, and allow titles/annotations to be
  removed. Named saved results provide matching result-plot routes. Exact
  fitted-source matching prevents mixing analyses; export replay reloads
  the saved results without refitting. These displays supplement location
  and fit diagnostics and do not classify raters or select scoring weights.

* GPCM help now relates the fitted model and inference targets to ConQuest and
  TAM, including TAM's documented many-facet slope construction. It explains
  why shared additive rater effects need not match mfrmr's slope-scaled rater
  effects, and why standardized-slope uncertainty needs a joint covariance
  transformation. Existing item-only numerical comparisons are distinguished
  from unverified interval comparisons and public import/export support.
* `confint(fit, parm = "slopes", level = 0.95)` supplies approximate pointwise
  relative-slope intervals for eligible GPCM MML solutions. `diagnose_mfrm()`
  uses the same joint-information log-Wald calculation; attached diagnostics
  carry its decision into the fit and weighting-review slope tables. Current
  likelihood, gradient, category, integration and unregularized-information
  checks replace the blanket interval restriction. Missing bounds retain
  reasons. The default retains geometric-mean-one relative slopes. Explicit
  options provide population-SD-standardized slopes, named ratios/differences
  with Wald tests, Bonferroni adjustment, and Person or larger independent-cluster
  sandwich covariance. Standardized intervals include scale cross-covariances;
  sandwich intervals do not repair bias or informative rating assignment.
* `bootstrap_mfrm_gpcm()` simulates independent Persons and scores under the
  fitted population on the analyzed assignment. Saved draws support basic
  bootstrap slope intervals; a matched PCM null instead supplies a bootstrap
  LRT. Every planned refit, warning and failure remains in the result. Unresolved
  draws widen interval limits and bound p-values rather than being discarded.
  Small-sample accuracy and misspecification robustness are not guaranteed.
* `mfrm_curve_intervals()` supplies logit/log delta intervals for category
  probabilities and per-rating information at specified ability values. Plots
  use lines and ribbons, color plus line types, and optional titles/notes.
  Bonferroni adjustment covers the requested grid, not a continuous band.
* IC/LRT and slope intervals now share the joint information calculation and
  a configurable dense-matrix workspace budget, replacing the IC/LRT-specific
  80-coordinate cutoff. `mml_quadrature_sensitivity()` preserves explicit
  population formulas, covariates and coding, and reports changes in slope
  interval endpoints and availability. A workspace estimate is not a total
  memory or capacity guarantee.

* `compare_mfrm(..., nested = TRUE)` now provides a PCM/GPCM equal-slope
  likelihood-ratio test for eligible MML pairs. It verifies the same population
  design, step/facet constraints and interactions, G-1 relative slope contrasts,
  and regular numerical solutions. The ordinary chi-square reference is
  asymptotic. Incompatible settings, non-unit weights, unstable solutions and
  negative likelihood gains retain explicit reasons instead of a p-value.
  `build_weighting_review(..., nested = TRUE)` carries the same test and status;
  its default remains `FALSE`. Slope intervals use their own checks.
* `compare_mfrm()` now separates GPCM MML information-criterion comparison
  from slope-interval eligibility. It reevaluates the retained likelihood,
  gradient and unregularized local information, while preserving the existing
  data, likelihood, parameter-count and integration checks. `ICFitEligible`
  and `ICFitReview` explain the solution decision; `ICComparable` controls
  deltas, weights and preferred candidates. Weighting/model-choice reviews
  retain that decision. These local checks do not prove global optimality or
  enable slope intervals or establish the nesting required by a likelihood-ratio test.
* User-facing model names, messages and help now use "GPCM" instead of
  "bounded GPCM". The single slope/step facet and each unavailable operation
  are described separately. The GPCM guide distinguishes slope intervals,
  information-criterion comparison and the PCM/GPCM likelihood-ratio test,
  including their distinct eligibility conditions. Unidimensionality is not a reason
  to prohibit these methods, and unit slopes are not a boundary null.
  The terminology change itself does not alter statistical eligibility.
* `mfrmr_output_guide("feedback")` now connects rater-feedback questions to
  fixed-facet intervals, ordinary or extended-model residual review,
  shared-rater uncertainty and known-truth screening evaluation. Each route
  identifies its required inputs, follow-up tables/plots and interpretation
  limits. This adds guidance for existing functions, not a new estimator or
  a rater-quality classification rule.
* The package description and MML tutorial now state the GPCM slope/step
  structure and its separate inference limits.
* GPCM diagnostics no longer treat missing or outdated slope-eligibility
  records as permission to report ordinary SEs or confidence intervals.
  Local covariance calculations remain in explicitly named `Optimizer*`
  diagnostic columns; regularized covariance does not authorize primary
  slope uncertainty. SE and interval eligibility are checked separately.
  Re-run `diagnose_mfrm()` or `confint()` on saved fits to refresh these tables;
  available slope intervals do not resolve other warnings about the fitted model.
* Use `mfrm_cluster_pam()` for Gower/PAM grouping and
  `review_mfrm_imputations(..., impute_ids = ...)` to check supplied score
  completions. The older `mfrm_cluster()` and `mfrm_response_imputations()`
  calls remain supported with their original arguments. Computations and
  saved-result classes are unchanged; guides, examples and the compatibility
  table identify the recommended names.
* `fit_mfrm()`, `describe_mfrm_data()` and `review_mfrm_anchors()` now accept
  `category_policy = "preserve"` or `"collapse"`. These choices correspond to
  `keep_original = TRUE` and `FALSE`; the default behavior is unchanged.
  Conflicting explicit old/new choices stop. A preserved, unobserved internal
  category still prevents fitting an unsupported step; this option does not
  estimate missing category information. Saved fits retain the existing replay
  representation.
  Data-review overviews and fit settings now report `CategoryPolicy` and
  `ScoreRecoded`, distinguishing the selected rule from actual recoding.
  Older saved results with missing records report an unknown policy or mapping
  rather than infer the choice from the observed categories.
* Fitting help groups arguments by data, rubric, model, computation and output.
  GPCM availability guidance separates its one-facet slope/step structure
  from limits on formal uncertainty and model comparison, without changing
  estimation or statistical eligibility.
* The GPCM tutorial and comparison help describe the separate checks for
  information-criterion ranking, the matched PCM/GPCM likelihood-ratio test
  and relative-slope intervals. A larger quadrature grid alone does not
  qualify any of these outputs. The weighting-review examples retain the
  corresponding output-specific decisions.
* Linking, reporting, workflow and visual-diagnostic tutorial figures now have
  descriptive alternative text, including both response-time views. Source
  archives include executed output and figures for all fifteen tutorials.
* Extended-model help examples now reuse packaged synthetic fits, scores,
  diagnostics and bootstrap results. Summaries, figures and interval-level
  changes run immediately; expensive recomputation calls remain visible as
  comments, with a complete regeneration script. Data, quadrature settings
  and numerical checks are preserved. The 19 bootstrap trials illustrate
  the workflow, not accurate tail quantiles or a coverage guarantee.
* JML help examples that previously stopped after 30 iterations now allow
  sufficient iterations for their example data to converge. Their estimator,
  data and model specification are unchanged; convergence warnings are not
  suppressed.
* ConQuest export help now explains the local command-file workflow and the
  separate controls used for posterior EAP calculations. Bundle guidance no
  longer treats every inference-readiness restriction as optimizer failure;
  numerical comparison remains distinct from validation of standard errors.
* The workflow guide now distinguishes defaults that change the statistical
  analysis from calculation and display controls. The opening example declares
  its rubric bounds, preserved categories and fixed ability population.
  Help clarifies that testlet membership does not automatically add fixed
  difficulty effects, and that an explicit dashboard `misfit_warn` replaces
  both review bounds. Existing defaults are unchanged.
* The workflow guide and README now explain function names and show how to
  choose an analysis from a concrete question. Extended-model scoring guidance
  points to `score_mfrm_persons()` for people already in the fitted model and
  distinguishes the model-specific meanings of `predict()`.
* Help introductions clarify which functions review supplied imputations,
  evaluate simulations, summarize attributes or plan G/D studies. The beginner output guide
  uses plainer questions and no longer describes a particular plot as mandatory.

## Design decisions and saved output

* Plot conversion now refuses unsupported PCA, clustering, pooled-MI,
  fixed-facet interval, screening-performance and D-study difference views
  even when `component` is specified. Previously, that option could silently
  replace the intended axes or intervals with a generic bar chart.
  Multivariate D-study scenario conversion retains its dedicated panels with
  `component = "series"`; use `plot_data()` for other tables and custom graphics.
  Non-supported views give that guidance even when ggplot2 is not installed.
* PCA score plots distinguish groups by both colour and point shape and retain
  the encoding in `plot_data()`. Scree plots distinguish retained and omitted
  components with filled and open points, including monochrome output.
  Tutorial figures now include alternative text describing the plotted target
  and the meaning of intervals, missing results and descriptive groups.
* `mfrmr_output_guide()` now includes `"features"`, `"imputation"` and
  `"gtheory"` routes, distinguishing their dedicated tables, plots and RDS
  saving from the fitted-model `mfrm_results()` reporting workflow.
* Extended-model interval plots and testlet Person scoring now also work in
  R sessions where the stats package is installed but not attached.
* Help for `fit_mfrm_imputed()` and `pool_mfrm_imputed()` now renders function
  links and code formatting correctly.

* Design-evaluation summaries retain full-precision metrics; `digits` now
  controls printing only. Rounding can no longer turn a value just below a
  minimum, or above a maximum, into a passing design recommendation.
  Rebuild older saved summaries from the original evaluation before requesting
  recommendations; no new simulation or model fitting is required.
* The default fit-level archive works without optional prediction objects.
  Predictions are included by default when supplied; explicitly requesting
  them without supplying an object still gives an explanatory error.
* Design-planning help distinguishes generating assumptions from analysis
  settings and explains numerical failures and separation reliability.

## Rater feedback

* Shared-rater and testlet calibration now show point estimates and approximate
  SEs with missing bounds by default. Explicit pointwise normal approximations
  remain available through `confint(fit, parm = "calibration", level = ...)`,
  `summary(fit, calibration_intervals = "normal", level = ...)`, and
  `mfrm_results(fit, calibration_intervals = "normal", calibration_level = ...)`.
  Testlet plots accept `intervals = "normal", level = ...`; common results
  preserve this choice in plots, reports and saved replay. Numerical failures
  and estimated variance boundaries still withhold bounds. This change does
  not improve or qualify finite-sample coverage. Rebuild outputs from older
  saved fits to apply it without refitting; existing result bundles keep their
  stored tables. Population-SD and conditional Person intervals retain their
  separate methods and interpretation.

* Shared-rater fitting help now includes a bounded independent comparison of
  local calibration-likelihood changes using saved joint posterior samples.
  The tested changes met the stated numerical tolerance; this does not qualify
  interval coverage, full population-SD profiles or variance boundaries.
  Extended-model reports and testlet summaries/figures now identify fixed-facet
  and step intervals explicitly as observed-information normal approximations.

* Testlet fitting now selects a converged optimization start when a failed
  start has an indistinguishable objective within floating-point precision.
  Gradient, integration, search-bound and information checks are unchanged.
  A 480-dataset comparison included fourteen fits affected by this
  start-selection issue. Refit affected saved models and regenerate their
  scores; printing an older unavailable result does not repair it.
* The testlet tutorial now reports an estimated-population comparison with
  ordinary RSM, including original failures, their separate numerical repair,
  Monte Carlo uncertainty and known-calibration references. Corrected Person
  coverage was 94.4% versus 84.9% with 120 Persons and positive local dependence,
  and 91.7% versus 83.3% with 24 Persons. Testlet EAP mean squared error was
  higher in those two balanced conditions. Conditional intervals retain their
  stated interpretation; neither a universal 95% coverage guarantee nor a
  uniformly more accurate model is established.
* Shared-rater fits now report `checks$PersonQuadratureStable` separately
  from overall numerical readiness. When Person integration is insufficient,
  the warning identifies it and explains refitting with more quadrature
  points; the tutorial gives the corresponding workflow. Integration
  tolerances, default order, estimators and interval restrictions are unchanged.
* Added a testlet application tutorial connecting task/criterion allocation,
  unequal rubric lengths, possible halo and task-specific dependence to the
  existing fit and scoring APIs. A complete four-rating enumeration compares
  model-conditional posterior precision and marginal EAP reliability; a
  five-versus-two-criterion example compares the influence of one-point changes.
  Saved scoring results provide plots, tables and report connections. Help
  distinguishes dependence from halo causation, task difficulty from local
  variance, and likelihood adjustment from equal task weighting. The current
  model still estimates one common local variance; no task-specific variance
  model, general planning API or new coverage claim is introduced.
* Individual shared-rater prediction intervals are no longer supplied
  automatically in fits, summaries, plots or report tables. Point estimates,
  conditional SDs and first-order `PredictionSE` remain available.
  `confint(fit, parm = "raters")` and `plot(fit, intervals = "normal")`
  explicitly inspect the normal approximation without refitting. The
  interval-width view and sorting require that explicit selection. New
  summaries/reports of earlier saved fits follow this default; raw saved
  objects and historical plot payloads retain their values. This restriction
  does not correct coverage or promote bootstrap intervals as a replacement.
  In a prespecified 800-dataset estimated-population study, conditional rater
  coverage was 91.6--91.8% with six raters and 94.0--94.1% with 24; finite
  interval availability was 99% and 87--88%, respectively. None met the joint
  coverage/availability criterion. All optimizers and information checks
  passed, but 53 fits failed higher-order Person quadrature checks, and one
  additional fit estimated zero rater variance. The tutorial reports Monte
  Carlo uncertainty and paired ordinary-MFRM point comparisons. These results
  do not qualify Person/testlet intervals or a minimum rater count.
* `score_mfrm_persons()` scores the complete source roster under ordinary
  RSM MML, testlet or shared-rater calibration, with optional Person selection.
  Ordinary scores use continuous posterior moments and equal-tail endpoints,
  retaining the fitted normal mean/variance. Missing-only Persons remain
  prior-only; zero ability variance does not produce zero-width intervals.
  `compare_mfrm(..., person_scores = ...)` checks matching rosters and compares
  EAPs after population-origin alignment, retaining conditional endpoints and
  unavailable rows without model-difference intervals or automatic ranking.
* Saved extension results support model-aware Wright maps and descriptive
  fit pathways through `type = "wright"` and `"fit_pathway"`. Facet/rater
  positions include the mean step at a zero-other-effect reference; only
  Person conditional intervals are drawn. The fit pathway uses posterior
  residual summaries, with no ordinary fit cutoffs. English labels, shapes,
  monochrome, text alternatives, panel/Person selection and optional
  annotations are available in native and ggplot views. Maps and their tables
  join static reports and saved replay without new scoring. Older testlet
  scores without roster metadata need rescoring from the existing fit.

* `mfrm_response_diagnostics()` now integrates latent uncertainty for ordinary
  RSM MML, testlet and shared-rater fits before computing category probabilities, predictive
  means/variances and descriptive Infit/Outfit. Output selection preserves the
  full conditioning roster; missing or unresolved rows remain explicit.
  Shared-rater predictions use normalized category-specific joint Laplace
  integrals and retain the normalization defect; numerical agreement does not
  certify approximation accuracy. Paired/scatter displays support ggplot,
  distinct symbols, monochrome, text alternatives and optional annotations.
  Saved diagnostics attach to `mfrm_results()` for reports and replay.
  These same-data summaries have no expectation-one reference, cutoffs, ZSTD
  or p-values; ordinary plug-in fit indices are not directly comparable.

* `compare_mfrm()` accepts matching saved `response_diagnostics` for an
  ordinary RSM versus a testlet/shared-rater model. Category probabilities,
  means, full variances and descriptive Infit/Outfit share one definition.
  Selected events are matched by contents and multiplicity, including missing
  scores, rather than row numbers. Paired/difference plots offer metric and
  category selection, monochrome and optional annotations. Saved comparison
  tables and figures connect to reports and replay without recomputation.
  The results are same-data descriptions, with no model preference or fit
  cutoffs. Ordinary plug-in diagnostics retain their existing definition.

* `fit_measures_table()` now uses mean squares for its default directional
  screen. ZSTD values and df comparisons remain visible; `ZSTDOnly` identifies
  additional combined-rule flags. Use `flag_basis = "mnsq_or_zstd"` for the
  earlier combined rule. Missing indices no longer establish a negative
  screen; `ScreenComplete` identifies incomplete rows. Low mean squares are
  described as low residual variability, not poor rater quality or grounds
  for automatic exclusion. Console follow-up prioritizes high mean squares.
  Saved reports identify the actual rule; earlier stored tables retain their
  combined-rule interpretation.
* `mfrm_screening_sensitivity()` compares declared threshold bands on saved
  simulation statistics, separating underfit, overfit and their union. It
  preserves planned trials, unknown outcomes and replication-level Monte
  Carlo intervals. Labelled tile and curve views support ggplot, monochrome
  output, text alternatives and optional titles/notes. ZSTD is an explicit
  alternative rule; no best threshold or universal error guarantee is selected.

* `compare_mfrm()` now accepts one ordinary RSM MML fit and one testlet or
  shared-rater fit for descriptive, matched-event facet comparisons. It checks
  category/population/constraint compatibility and repeated or omitted events,
  centers effects over the same facet levels, and retains raw values and failed
  checks. Paired and mean-versus-difference plots support ggplot conversion,
  custom/hidden annotations and text alternatives. Matching comparisons attach
  to `mfrm_results(..., comparison = ...)` for reports and saved exports.
  No automatic ranking or difference intervals are implied. Person-score
  comparison requires the saved scoring outputs described above. Predictive comparisons require the separately supplied diagnostics
  described above. Ordinary fits now retain omitted events and input-row
  indices; older fits with unrecorded omissions require refitting for this route.

* Shared-rater and testlet estimate plots now offer interval, precision
  (estimate versus interval width), and empirical cumulative-distribution
  views. Shared display controls cover sorting, monochrome output, point/text
  size, labels, reference lines and custom or hidden titles/captions. Colours
  are supplemented by symbols; plot data retain text alternatives, complete
  tables and exclusion reasons. Cumulative distributions exclude prior-only
  scores and do not estimate a latent population distribution. Shared-rater
  severity and bootstrap intervals now also convert to ggplot; bootstrap
  comparisons retain dashed ordinary intervals and arrows for infinite ends.
  These presentation changes reuse saved estimates without refitting.

* `score_mfrm_random_rater()` supplies conditional Person EAPs, posterior SDs
  and continuous equal-tail intervals with joint shared-rater integration.
  Selecting Person IDs retains the full scoring roster; supplied new data
  replace the roster. Prior-only and unavailable rows stay explicit. The
  conditional Laplace approximation holds calibration fixed and does not
  establish general approximation accuracy or interval coverage. Saved scores
  support plots, ggplot conversion and `mfrm_results(fit, scores = scores)`
  alongside response predictions and rater bootstrap intervals. Numerical
  checks, omitted rows and the complete scoring roster remain in the result.
  A fixed-calibration joint-posterior comparison now supports the stated
  numerical tolerances for 48 selected scores on twelve linked and reduced
  three-category RSM rosters. Help separates approximation accuracy from
  interval coverage, retains the reference-precision limitations, and explains
  how to select a few output Persons while keeping the full scoring data.
* `fit_mfrm_random_rater()` now estimates normal ability variance by default.
  Use `person_sd = 1` for a known standard normal population, or another
  positive known SD. Calibration-adjusted rater uncertainty includes the
  ability-SD coordinate; rater-SD profiles re-estimate it, and bootstrap
  generation/refits preserve the fitted population and known/estimated choice.
  Both estimated zero variances and the ability search bound are checked.
  Up to 241 Person quadrature points can be requested when lower orders fail
  the unchanged accuracy checks.
  An estimated ability boundary withholds regular intervals and rater-SD
  profiling. Reports distinguish both variance components; earlier saved fits
  retain their known N(0,1) population. Existing fixed-population coverage
  evidence does not qualify this estimated-population workflow.
* `fit_mfrm_testlet()` now estimates normal ability variance by default,
  alongside the common local variance. Use `person_sd = 1` for the earlier
  fixed N(0,1) model, or specify another known positive SD. Both estimated
  zero-variance boundaries are checked explicitly. Zero ability variance
  withholds individual scores; search-bound failures retain their checks.
  Conditional scoring, prior-only rows, summaries, reports and saved results
  use the matching population. Earlier saved fits retain N(0,1); changing
  that assumption requires refitting and regenerating scores. Person intervals
  exclude uncertainty in estimating both variances; earlier fixed-population
  coverage studies do not qualify this new fit.
* Random-rater and testlet help now connects the effect-sharing assumptions
  to measurement-model research, distinguishes the fixed ability-variance
  restriction from scale identification, and explains which ordinary-MFRM
  comparisons preserve the data and statistical target. Descriptive matched-facet
  comparison uses `compare_mfrm()`; predictive diagnostics and their cross-model
  comparisons use the matching response-probability definition.
* Testlet and random-rater fits now connect to `mfrm_results()`,
  `mfrm_report()` and `export_mfrm_results()`. Reports retain calibration,
  numerical checks, omitted-score accounting and interval meanings, plus
  explicitly supplied matching predictions and random-rater bootstrap results.
  CSV/RDS exports preserve prior-only and unavailable rows, failed bootstrap
  trials and infinite interval endpoints. Replay reloads saved results without
  fitting, scoring or resampling. Older predictions must be regenerated from
  their saved fit to attach source metadata; no model refit is needed.
  Numerical checks are not model-fit diagnostics or general coverage evidence.
* Testlet results now connect to `as_ggplot()` and the common plot-data
  accessors, preserving interval labels, prior-only symbols and unavailable
  rows. `mfrmr_output_guide("models")` and the workflow help compare fixed
  facets, shared random raters and testlets, including their different
  prediction targets. The tutorial covers figure and table export. Ordinary
  diagnostics and the Shiny viewer remain unsupported for the extended
  models and provide targeted guidance. Replot or convert saved results
  to use the display changes; no refit or rescoring is required.
* Added `fit_mfrm_testlet()` for an RSM with explicit, non-overlapping
  Person-specific testlet memberships, additive fixed facets and a common
  normal local-effect variance. Nested quadrature, analytic gradients,
  multiple starts, an exact zero-variance submodel and information checks
  support calibration; unequal observed block sizes are allowed. `predict()`
  provides continuous conditional Person intervals, retaining unavailable
  and prior-only rows. Saved fits, English plots, an interval-guide entry and
  a complete tutorial cover fitting, feedback, missing assigned scores and
  reuse. Estimated zero variance withholds regular calibration intervals.
  Person intervals hold calibration fixed and do not propagate its estimation
  uncertainty. This local effect is distinct from a rater shared across
  persons. The separate result class does not support ordinary-model
  diagnostics, response-MI pooling or portable calibration conversion.
  Earlier models require a new fit with explicit membership. General coverage,
  calibration-aware Person intervals and correlated/heterogeneous testlets
  remain unfinished.
* Added `mfrm_random_rater_intervals()` to compare model-based bootstrap
  prediction intervals with ordinary intervals for observed random raters.
  It generates shared rater effects, persons and scores on the analyzed
  assignment and refits the same model. Saved prediction errors support
  studentized and unscaled intervals and different levels without refitting.
  Failed refits and unavailable studentizers remain in the planned count;
  uncertainty about their roots widens limits, possibly to infinity.
  Plots distinguish unbounded endpoints from finite ordinary intervals.
  Positive source SD and prediction SEs are required. This candidate does
  not replace population-SD profiles or establish general coverage; boundary
  sources, simultaneous classification and missingness-model uncertainty
  remain outside its scope. In a separate 24-dataset comparison (12 per
  condition, 99 refits each), six-rater coverage was 91.7% for ordinary and
  100% for studentized intervals, with 40% greater mean width. The improvement
  occurred in one dataset; its paired change had MCSE 8.3 percentage points.
  With 24 raters, coverage was 95.5% and 95.1%. All intervals were finite.
  This small comparison does not qualify 95% coverage or justify making
  bootstrap intervals the default. The interval guide and tutorial distinguish these
  realized-rater intervals from population-SD profiles.
* Added `fit_mfrm_random_rater()` for an RSM with one normal rater effect
  shared across all persons rated by that rater. Frequentist approximate MML
  uses Person quadrature and a joint rater Laplace integral through optional
  RTMB (>= 2.0). The distinct result class retains numerical checks, full
  covariance, observed-rater feedback, plots and saved-result support.
  `confint()` profiles the population SD, including zero. Estimated variance
  boundaries withhold ordinary calibration and rater intervals. `predict()`
  distinguishes observed raters from replacement raters at specified abilities;
  it returns marginal score probabilities conditional on calibration, not
  joint predictions, person scoring or parameter-uncertainty intervals.
  A complete tutorial explains population assumptions, missing assigned scores,
  numerical failures and interpretation. PCM, anchors, testlets, latent
  regression and multidimensional abilities are outside this route. Earlier
  fixed-rater analyses require a new fit for this model; changing a saved
  object's class is not a migration. A bounded 160-dataset normal-population
  pilot returned all population-SD profile intervals. At true SD 0.7,
  individual-rater interval coverage averaged 84.2% with six raters and 94.4%
  with 24 raters (40 datasets each). Few-rater interval calibration remains
  unresolved; these results do not establish general coverage.
* Added `mfrm_screening_performance()` and its plot for evaluating declared
  screening rules against known simulation truth. A complete planned roster
  retains failed and incomplete trials. Per-target and any-target family rates
  keep separate denominators, exact Monte Carlo intervals and bounds for
  unresolved outcomes. A tutorial connects simulated ratings to rater flags.
  In a bounded matched-budget study, the specified Infit/Outfit union detected
  only 6/100 and 2/100 contaminated-rater cases under two sparse assignments;
  absence of a warning cannot certify rater quality.
* Corrected `evaluate_mfrm_signal_detection()` so unavailable bias statistics
  and descriptive-only DIF classifications remain unavailable rather than
  counting as negative screens. Target summaries retain full precision,
  planned/available counts and Monte Carlo bounds; non-target results retain
  cell availability. Rebuild saved summaries and plots from the original
  evaluation for target corrections. Earlier non-target averages lack the
  required cell counts and are withheld; recovering them requires rerunning
  the simulation. These are conditional screening rates, not general accuracy
  guarantees or familywise error control.
* Added `mfrm_facet_intervals()` for pointwise intervals on fixed facet estimates
  and named linear contrasts, including rater differences. For inference-ready
  fixed-standard-normal RSM/PCM MML fits with unit weights, compare ordinary
  observed-information covariance with a one-way sandwich using whole persons
  or explicitly declared larger independent clusters. The full covariance
  preserves fitted constraints and covariance between facet levels.
  A plot and complete tutorial compare methods on the same estimates and
  retain unavailable intervals. The sandwich describes variability around the
  working model's target; it does not correct misspecification bias or provide
  general coverage, small-cluster, random-rater or diagnostic-accuracy guarantees.
  `mfrmr_interval_guide()` includes this route and assigned-score MI pooling.
* The introductory CSV examples preserve literal `NA` identifiers and recode
  missing-score markers only in the score column. If an earlier import
  removed valid IDs, reimport the original file and rerun the analysis.
* The workflow tutorial connects coverage and category review to rater
  severity, uncertainty and screening flags, with focused tables and guidance
  on interpreting individual ratings separately from rater-level summaries.
* Reduced overhead when assembling strict pairwise diagnostic summaries for
  large rating datasets, preserving statistics and unavailable-result flags.

* Facet dashboards identify unavailable diagnostics and preserve fit-readiness
  restrictions in summaries and plots. Zero observed flags is not a complete
  diagnostic pass. Rater severity profiles label restricted fits `REVIEW ONLY`.
  A nonfinite Infit or Outfit value cannot hide a flag from the other,
  finite index or create a flag on its own.
  Rebuild saved dashboards from existing fits and diagnostics to obtain the
  added availability and readiness fields; a model refit is unnecessary.
* Dashboard plots now use the same default misfit band as dashboard tables and
  retain a saved dashboard's thresholds for their guide lines. Recreate plots
  to correct earlier defaults or guide lines; no model refit is needed.
* Dashboard exports include screening settings and interpretation notes in
  CSV/text and HTML. Re-export existing fits to retain this context.
  HTML also shows unflagged levels with unavailable diagnostics and uses
  focused columns; the CSV files retain the full diagnostic detail.

## Multivariate G-theory and assessment planning

* Added a fixed-task-set example using existing multivariate score components
  and a random-rater facet. It connects task-specific scores, prespecified
  weights, rater-count scenarios and plots, and checks the result against a
  directly weighted score. Help distinguishes this target from sampling new
  tasks, explains workload and common-rater requirements, and retains whole-row
  omission accounting for incomplete score vectors. No general fixed-facet
  estimator is added.
* Added `mfrm_g_study()` and `mfrm_d_study()` for
  numeric observed scores with fixed score components and one or two random
  measurement facets. A single score is also supported. Select common tasks,
  raters or other named facets through `rater`/`task` or `facets`. ANOVA handles
  complete balanced data; explicit `method = "minque0"` also handles incomplete
  designs whose moment equations separate the covariance components.
* Specify `nesting = c(Rater = "Task")` for persons crossed with task-specific
  rater teams. This estimates five components, with child identities local to
  each parent and shared across persons and scores. MINQUE(0) also permits
  unequal observed team sizes. D-studies preserve the fitted nesting and
  project complete balanced plans; `Raters` means raters per task. Nesting
  within persons, partly shared children, score-specific identities and unequal
  future allocations remain unsupported.
* D-studies report original scores and named weighted composites, including
  signed differences. Vector or matrix weights are used without normalization;
  between-score covariances enter composite projections. Incomplete or unequal
  source designs require an explicit future grid. Observed pool counts are not
  future per-person replication or the reliability of a sparse roster.
* Raw negative or indefinite covariance estimates are retained. Each score,
  composite and metric has its own availability status: a non-PSD component
  no longer suppresses every G/Phi/SEM. `ComponentPSD` separately flags the
  component matrices for review. Zero universe variance gives zero G/Phi when
  the corresponding error is positive; zero total variance is undefined.
  Recompute saved D-studies from their G-studies to adopt these rules.
* Added `mfrm_multivariate_d_compare()` for prespecified differences in G, Phi
  and SEM between future plans from the same two-crossed-facet G-study.
  Explicit `assumption = "normal"` requests approximate pointwise paired-delta
  intervals, preserving dependence between plans. Complete ANOVA and incomplete
  MINQUE(0) sources are supported. Nonnormal-robust, one-facet, nested and
  simultaneous intervals are not provided. Unavailable intervals retain their
  reasons; a returned point estimate need not have an available interval.
* G/D-study help provides runnable examples, including task-only, named-facet,
  nested-team and composite planning. It distinguishes ranking from absolute
  decisions, measurement error from sampling uncertainty, planned assignments
  from missing scores, and component identification from precision. G/Phi are
  not pass/fail accuracy or latent MFRM reliability. Omitting missing rows does
  not impute scores or correct informative assignment. Guidance explains raw
  negative-component conventions when comparing GENOVA-family results.
  Updating only future counts or composite weights reuses the G-study;
  changing its source data or model requires a new G-study. Saved-result
  guidance distinguishes these actions from replotting stored values.
* D-study base plots and ggplot exports preserve selected scores, weights,
  planned counts, metric-specific omissions and component diagnostics. Nested
  child counts are labeled per parent. Plan-comparison plots show differences
  and intervals against zero; `plot_data()` supplies their exact values.
  Automatic ggplot conversion is unsupported for plan-comparison plots.
* Corrected single-score MINQUE(0) matrix handling and kept distinct interaction
  IDs containing periods separate. Rerun affected G-studies. Main-effects
  D-study curves now distinguish every other facet count and retain missing
  estimates; surface plots require other conditions to be fixed within panels.
  Their automatic ggplot fallback is disabled because it could misrepresent
  the selected quantities.
* `mfrm_generalizability()` now requires explicit `missing = "omit"` for missing
  scores or selected facet values, refuses malformed scores, and retains input,
  used and excluded row accounting. Older saved results require a new G-study
  from original data to obtain that accounting; an MFRM refit is unnecessary.

## Missing scores on assigned ratings

* A paired repeated-sampling assessment of joint-RSM imputations now informs
  the help. In 200 datasets with 80 Persons, three raters and two criteria,
  MI coverage of a fixed-rater contrast was 96.5% under observed-score-dependent
  MAR missingness (198 available intervals; 95% Monte Carlo bounds 92.9--98.6%).
  Low-score-dependent MNAR missingness reduced coverage to 26.0% and introduced
  +0.574-logit bias despite a similar overall missing fraction near 14%.
  Direct MML and Bayesian inference showed the same vulnerability. Two MAR
  posteriors missed the sampling-diagnostic threshold and remain counted as
  unavailable. These are bounded results for the specified imputer and design;
  the generic supplied-imputation and pooling calculations are unchanged.
* The assigned-score tutorial now uses a joint RSM imputation example with
  uncertainty in calibration and shared Person abilities. Forty supplied
  posterior predictive completions run without a Bayesian toolchain; an
  accompanying R/Stan script regenerates them and retains sampling diagnostics.
  A probability display, direct observed-score comparison and separate
  lower-score sensitivity analysis connect missing-score assumptions to a
  fixed-rater contrast. The example preserves unassigned cells and observed
  scores. Agreement in one example does not establish repeated-sampling
  coverage or exact compatibility between Bayesian imputation and MML analysis.
* Pooling help distinguishes uncertainty in an imputation model from bootstrap
  inference for the entire analysis. Imputer resampling or a common response
  likelihood alone does not establish valid Rubin intervals. The pooling
  calculation is unchanged. Direct observed-score MML can estimate the same
  target under the same likelihood and ignorable missingness assumptions;
  imputation creates no additional observed information.
* Added `mfrm_response_imputations()` to review supplied score completions from
  a retained imputation model. Explicit event IDs and eligible missing-score
  selections preserve observed scores, person/facet identities and assignment.
  Unassigned rows remain unassigned; no complete rating grid is constructed.
  Unselected assigned nonresponse requires an explicit omission decision.
* `fit_mfrm_imputed()` fits every completion with a common RSM/PCM MML
  specification, category ladder and fixed-standard-normal person scale.
  All fits, failures and warnings remain available; incomplete or ineligible
  analyses cannot be silently omitted from subsequent pooling.
* `pool_mfrm_imputed()` combines eligible non-person facet estimates or
  prespecified linear contrasts with full observed-information covariance and
  Rubin's rules. Results retain within/between/total covariance, reference
  degrees of freedom and Monte Carlo SEs. A plot displays the stored pointwise
  intervals and identifies fixed targets without inferential intervals.
  Person EAPs and posterior SDs are not eligible for this pooling route.

## External features and exploratory groups

* Added `mfrm_pca()` and `mfrm_cluster_kmeans()` for explicitly selected
  numeric external attributes. They retain scaling, feature weights, PCA
  component choices, initialization settings, original-unit profiles, omitted
  IDs and missingness reasons. PCA input to k-means is used without whitening
  or rescaling; factors and arbitrary category codes are not converted.
* PCA provides scree, score and loading plots with exact plotted tables;
  matching group results can annotate the score view. Existing profiles,
  silhouettes, saved results and group comparisons support k-means. Comparison
  summaries identify each analysis's distance, fitted space, scaling and
  retained component count. Recreate comparisons from saved analyses to add
  these columns; the existing analyses need no refit. Different geometries do
  not make silhouettes a common-scale criterion for choosing a feature set.
* `mfrm_cluster_imputed(method = "kmeans")` supports direct numeric features or
  explicit PCA reduction within each retained completion. It reuses the same
  reviewed `mice` model and preserves all completions and their transformations.
  No PCA basis, component scores or group labels are pooled. Updated examples
  connect numeric preprocessing, paired comparisons and original-unit profiles.
* K-means can omit silhouettes explicitly to avoid their pairwise distance
  allocation. Uncomputed values remain unavailable in summaries and plots.
  Imputation co-membership matrices still require quadratic memory. Seeds
  preserve the caller's random state; optimization failures are reported.
* The external-feature tutorial now uses `mice` chain plots to review
  imputation behavior, so the basic example does not require `rstan`.
* Added `mfrm_features()` to prepare one row per Person, rater or task from
  selected external attributes, preserving types, IDs and missingness reasons.
  `mfrm_cluster()` uses Gower distances and PAM; `mfrm_cluster_hierarchical()`
  uses average or complete linkage and retains a dendrogram. Users choose
  features, weights and group counts. Results include profiles, silhouettes
  and, for PAM, representative entities. Groups do not establish latent classes,
  ability levels or rater quality.
* Missing features stop clustering by default; explicit complete-case selection
  preserves excluded IDs with unavailable memberships. Numeric-range overflow
  and relative-weight underflow are refused rather than silently changing
  feature contributions.
* `mfrm_cluster_imputed()` uses user-fitted `mice` completions and explicitly
  selected eligible missing cells. IDs, observed values, remaining missingness,
  reasons and model diagnostics are preserved. Pairwise co-membership proportions
  describe sensitivity to those completions, not membership probabilities or
  pooled inference. A failed analysis stops the comparison; hierarchical trees
  remain specific to each completion.
* `mfrm_cluster_compare()` compares fitted groupings across group counts,
  weights, feature selections and methods without refitting. It reports
  label-invariant pair changes, adjusted Rand indices, sizes and silhouettes.
  Analyses must contain the same entities, and shared features must retain their
  values and types. Imputation comparisons must reuse the same fitted `mice`
  object and preserve completion pairing; no setting is selected automatically.
* Plots show silhouettes, feature profiles, stored hierarchies and imputation
  co-membership heatmaps. `plot_data()` retains plotted values and exclusions,
  distinguishing unavailable pairs from zero co-membership. Automatic
  `as_ggplot()` conversion is refused because its generic fallback could display
  cluster numbers or mislabel the intended quantities.
* The external-feature tutorial uses fictional experience, workload, specialty,
  training and certification attributes. It covers separate Person/rater/task
  groupings linked to assignments, missingness and imputation choices, setting
  comparisons, sparse designs, and the roles and limits of PCA and k-means.

## Updating saved analyses

Keep the original objects and analysis settings. Installing this version does
not update saved numerical results or previously exported reports. See
"Updating saved analyses for 0.2.4" in
`help("mfrmr_workflow_methods", package = "mfrmr")` for the starting object
and action needed for each affected workflow.

* Reprinting updates display wording. Re-summarizing can update summaries
  derived from retained scores, draws or simulation runs; it cannot replace
  an old diagnostic calculation or grid-based interval.
* Most corrections below can reuse a fit with current estimation checks.
  Older native fits lacking those checks must be refitted from their original
  data and settings before inferential reuse. Recomputing diagnostics alone
  cannot supply missing fit checks.
* Recreate affected results before exporting again. Valid portable calibrations
  retain their recorded scoring algorithm; adopting continuous intervals
  requires a newly created calibration. External results can be re-imported
  from saved source-package fits without re-estimation.

## Portable calibration and scoring

* The separate-session CSV scoring example now preserves text identifiers
  such as `001`, `1` and literal `NA`. Earlier automatic CSV conversion could
  merge or replace identifiers before scoring. Reimport affected response
  files with the corrected column types and rescore; the calibration does
  not need refitting.

* Added portable calibration for one-scale RSM/PCM MML fits with a fixed
  standard-normal Person distribution. A reviewed calibration can be saved,
  transferred and used to score compatible new Persons without the source fit
  or training responses. Extraction, validation, freezing and scoring preserve
  categories, facet labels, anchors, supported two-way facet interactions and
  the scoring basis. Score intervals
  condition on the saved calibration and prior; they exclude uncertainty in
  the estimated calibration parameters.

* Portable calibration extraction now requires the reviewed object's exact
  highest-grid fit together with its quadrature-sensitivity result. Users
  choose the evaluated grids and judge the observed movement; response-linked
  fits remain outside the portable artifact and should be archived separately
  when needed. Extraction checks each retained fit's actual grid, convergence,
  score map, model and fitting settings, including anchors and integration
  mode; inconsistent review objects are rejected.

* Portable calibration artifacts now record an explicit scoring algorithm and
  a scoring grid independent of the source fit's integration grid.

* Portable score tables and summaries retain the estimate and uncertainty
  basis, calibration identifiers, scoring algorithm and requested interval
  level when exported to CSV. The interval level remains unrounded. Older
  score results recover the algorithm and level from their stored settings
  when re-summarized. Score estimates, SDs, and interval values are unchanged.

* Portable score `print()` and `summary()` output now stays compact at ordinary
  console widths, direct calibration-method help is available, and base and
  ggplot2 score displays both distinguish review states by colour and shape.

* Calibration displays describe the scoring prior, and score summaries explain
  reasons for review in plain language. Printed scores and interval plots state
  the requested interval level and distinguish continuous posterior quantiles
  from saved grid endpoints, whose posterior mass may differ from that level.
  Integration-review summaries explain that numerical comparisons do not
  automatically classify stability or change readiness, without internal codes.

* Portable calibration score batches now have concise `print()` and structured
  `summary()` methods plus interval, precision, and quadrature-edge review
  plots. The displays retain review and not-scored dispositions and state that
  posterior uncertainty excludes calibration-parameter uncertainty.

* FACETS-facing scope guidance now distinguishes the supported native
  portable-calibration workflow from unsupported FACETS or third-party
  calibration-file import, and removes stale version-number wording from
  current capability statements.

## Changes affecting existing analyses

* Normal-model shrinkage no longer loses positive plug-in SDs or retained
  weights when shrinkage factors round to one. The moment calculation avoids
  overflowing sums of squares; unrepresentable variances leave the original
  estimates unchanged with unavailable factors and an explanatory note.
  `shrink_person = TRUE` now requires ordinary JML: EAP posterior SDs are not
  sampling SEs for another normal-model update. Reapply shrinkage to saved fits
  with the original prior settings (`shrink_person = FALSE` for other
  estimators), then regenerate reports, plots and exports. Updated reports
  carry `CalculationVersion = 2`; no MFRM refit is needed.
  `as_ggplot()` now preserves the per-level original/shrunk comparison and
  saved plug-in bands instead of converting this view to facet-level bars.

* Reapplying empirical-Bayes shrinkage now replaces previous Person adjustments:
  `shrink_person = FALSE` removes stale Person shrinkage columns. Replay scripts
  retain the latest post-fit adjustment as a separate step after the original
  fit, including when diagnostic SEs were attached before shrinkage. Refresh
  older shrinkage results with the original prior and supported Person settings before
  regenerating replay scripts. Automatic replay selects fit mode for workflow
  objects with post-fit shrinkage.
* ICC calculations no longer reject positive variances using a fixed cutoff
  tied to score units, and returned variances are no longer rounded to six
  decimal places. Constant retained scores yield unavailable variances and
  ICCs. Constant-response bootstrap refits count as unavailable draws and
  withhold intervals; rerun saved bootstrap results to apply this check.
* Design-effect help and printed summaries now identify the per-facet,
  average-cluster-size approximation. `EffectiveN` is a descriptive equivalent
  row count, not a count of independent Persons or a precision estimate for
  the full crossed, nested, weighted, or unequal-cluster-size design. The
  design-effect formula itself is unchanged.
* ICC analyses now preserve numeric character/factor score labels and reject
  malformed scores. Missing scores or selected grouping values stop the
  analysis unless `missing = "omit"` is explicit. Results retain input, used,
  and excluded row counts, excluded row positions, and missing columns.
  Design effects use the ICC model's retained observation and grouping-level
  counts, so excluded rows no longer inflate cluster sizes or effective sample
  sizes. Rerun older ICC results before calculating design effects. Omission
  does not impute scores or correct missing-data bias.
* Withdrew `ci_method = "profile"` from `compute_facet_icc()` and
  `analyze_hierarchical_structure()`: transforming separate variance-component
  intervals while fixing the other components did not produce a
  profile-likelihood interval for the ICC ratio. Choose `"boot"` explicitly
  for parametric percentile intervals, subject to the fitted Gaussian model.
  The interval guide now also points to this explicit bootstrap route.
* ICC bootstrap results now retain all draws, requested and unavailable counts,
  convergence/singularity diagnostics, and refit warnings/errors. Any missing
  or nonconverged replicate, or a fit warning, withholds intervals; successful
  draws are not silently selected for quantiles. Converged boundary fits remain
  in the distribution. Invalid replicate counts, CPU counts, and seeds are
  rejected rather than truncated. Snow workers now load lme4 before refitting.
  Saved interval results must be rerun from
  their original data/settings before printing, summarizing, or plotting with
  this version. ICC plots explain unavailable intervals and no longer apply
  reliability-band reference lines to every facet's variance share.

* Selecting rows or columns from a D-study table now preserves its calculation
  and interpretation information. Printing selected columns no longer mistakes
  a newly computed result for an older saved result. Coefficients are unchanged;
  older results still require regeneration as described in the help.

* Interrater summaries now retain unavailable comparisons and report the number
  of classified pairs. Expected agreement is withheld when category probabilities
  are incomplete or repeated ratings have been averaged within a context.
  Context identifiers containing separators no longer merge distinct contexts.
  Plots retain unavailable pairs and omit unestimated self-agreement values.
  Agreement and network reports require diagnostics matching the supplied fit.

* Rater networks now exclude zero-weight edges and leave directional indices
  unavailable when a rater has no retained directional comparisons. Distances
  refer to reachable pairs. Design reviews distinguish selected subsets from
  the full observed design and withhold a complete assessment when checks are
  missing. Halo reviews retain unavailable comparisons and no longer provide a
  Welch test of dependent pair correlations; its compatibility columns are
  missing. Recreate older agreement and rater-network results from the existing
  fit, matching diagnostics and original settings; no model refit is needed.

* Response-time reviews now read numeric factor labels correctly, retain groups
  with no valid times, and report valid and excluded row counts. Quantiles are
  computed correctly when an explicitly supplied cutoff argument evaluates to
  `NULL`; reversed rapid/slow cutoffs are rejected. Printed flags explain the
  descriptive rules without internal codes. Rates use valid times only and do
  not establish rapid guessing or low effort. Recreate older reviews from the
  original timed-event data and settings; no MFRM refit is needed.

* Model summaries now give one consistent formal-inference decision. The full
  display no longer reports inference as ready solely because the fit checks
  passed when precision support remains unassessed. Optimizer success is
  reported from the stored return code.

* Fit and summary displays explain population assumptions, GPCM discrimination,
  diagnostic availability and information-criterion restrictions without
  internal implementation/status labels. JML summaries state that no normal
  population was estimated. Structured results and numerical values are
  unchanged; reprint saved summaries to use the revised display.

* Fitted-object Person scores and posterior draws now retain their prior,
  calibration method, weighting and conditional uncertainty meaning in exported
  tables; score tables also retain the interval level. Summaries and HTML output explain source-fit
  restrictions without internal status codes. JML-based scoring is labelled
  as posterior EAP with a standard normal prior, rather than a new JML estimate.
  HTML bundles show basic model details and the public interpretation decision
  in place of raw manifest dumps;
  detailed manifest files remain available for reproducibility.

* Plausible-value summaries now use the requested interval level for empirical
  draw quantiles, instead of always using 95%. These finite-draw limits are
  distinguished from the companion continuous posterior interval. Existing
  scores, posterior SDs, posterior interval bounds and draws are unchanged.
  Re-summarize saved results and recreate exports; no refit or resampling is
  needed. To retain numerical priors absent from older estimated-population
  results, re-score using the existing fit.

* Imported eRm item difficulties now reverse source easiness signs and convert
  cumulative category coefficients to adjacent thresholds. TAM imports derive
  item difficulties and thresholds from category logits instead of coefficient
  names. Multi-facet TAM imports report combined response-condition difficulties
  and retain the original coefficient table; they do not reconstruct separate
  facet measures. Transformed SEs remain missing when joint covariance is needed.

* mirt imports now include rating-scale offsets and require a supported
  unidimensional Rasch or partial-credit model with positive slopes. Additional
  factor scores can no longer be mistaken for SEs. EAP estimates and conditional
  posterior SDs are labelled explicitly. For TAM, source person-fit statistics
  use WLE scores while imported person estimates remain EAP.

* Imported summaries describe the source scale without assuming native mfrmr
  population or slope constraints. Wright maps show point estimates only, and
  source posterior SDs no longer produce separation reliability. Re-import
  existing source-package fits to update saved bundles and derived displays;
  no model re-estimation is needed. Native model curves, response-level QC and
  portable calibration require a native fit.

* Marginal category and agreement diagnostics retain unavailable probabilities,
  counts and classifications. Missing context pairs no longer disappear from
  agreement summaries, and complete-scope residuals/RMSDs are withheld when
  contributing inputs are missing. Known cutoff crossings remain flagged when
  a companion rule is unavailable. Summaries, reports and exported tables show
  classified and unclassified counts.

* GPCM marginal diagnostics now summarize each declared threshold family
  separately. Marginal plots select from all rows within the requested facet
  and metric, retain unavailable rows, and display availability counts. Supplied
  diagnostics must match the fit. Ordinary marginal output explains its basis
  without internal status or literature codes.

* Marginal expectations condition on the same responses through Person
  posteriors, with calibration held fixed. Help now explains that residual
  scales omit cross-response covariance and calibration uncertainty; the
  cutoffs remain descriptive review rules. Recreate older marginal diagnostics,
  summaries, plots and exports from the existing fit and original diagnostic
  settings. No model refit is required.

* Rating-scale tables, QC and report text retain unavailable category counts
  and threshold comparisons. Adjacent steps are compared within each threshold
  family without skipping missing estimates; missing families require review.
  A binary scale has no applicable adjacent-threshold comparison. Equal
  thresholds count as nondecreasing, a description of point estimates rather
  than evidence that categories function adequately. Threshold plots label
  families and connect only available adjacent steps within a family.

* Missing category fit statistics and expected counts no longer become zero
  flags or counts. Usage totals cover all declared categories even when unused
  rows are hidden. Printed category summaries explain availability in plain
  language. Recreate older rating-scale/category-structure results, QC and
  reports with the original fit, diagnostics and settings; no model refit is
  required.

* Mean-square and ZSTD summaries now retain classified and unclassified
  element counts consistently in diagnostic summaries, report text and output
  tables. Missing results no longer imply zero flags or an overall fit failure.
  A known threshold crossing remains flagged when the companion statistic is
  missing; overall flag rates require complete classification.

* Fit rankings no longer mix ZSTD values and mean-square deviations under one
  label. Rankings use complete paired ZSTD values, falling back to mean-square
  deviations only when no paired ZSTD is available. Report text describes
  interval availability without internal column names. Recreate older
  diagnostic summaries, reports and exports from existing diagnostics; no
  model refit is required.

* Person-fit calculations use small positive probabilities without a lower
  floor. Missing or invalid response information makes the affected person's
  statistic and flags unavailable, with response counts retained. The ability
  correction requires numerical convergence, a finite person estimate, unit
  observation weights and complete derivative information. Otherwise, the
  uncorrected index is labelled explicitly when it can be calculated.

* Person-fit summaries explain limitations in plain language, and plots retain
  persons whose Infit/Outfit values are missing when their likelihood-based
  index is available. Normal-reference cutoffs are screening aids; they do not
  guarantee individual or multiple-person false-positive rates. Recreate older
  saved person-fit results with the matching fit and diagnostics; no MFRM
  re-estimation is needed.

* Unexpected-response summaries retain evaluated and unclassified response
  counts. Overall rates and before/after reductions remain unavailable when
  classification is incomplete, and QC requests review. Older saved screening
  results and QC summaries require regeneration with the original settings and
  existing fit. No model refit is required.

* Q3-style screens now count each unordered pair once and retain unavailable
  pairs with their shared-person counts and reasons. Missing correlations and
  flags remain unavailable. Printed summaries distinguish available pairs
  from all candidate pairs, and supplied diagnostics must match the fit.

* Residual PCA no longer replaces undefined correlations with zero or smooths
  invalid correlation matrices. Affected analyses remain unavailable with an
  explanation, which the visual-diagnostics tutorial now checks before plotting.
  Residual-permutation comparisons also remain unavailable if
  any requested permutation fails, while retaining the successful count.
  Recreate older PCA and Q3 summaries with the original settings and existing
  fit or diagnostics; the MFRM fit does not need to be re-estimated.

* Residual-PCA plots and summaries now describe exploratory evidence without
  fixed-value claims about a second dimension. Permutation references explain
  that fitted-model uncertainty is omitted. Facet-equivalence summaries use
  plain-language decisions; plots clearly describe deviations from the facet
  mean, separately from pairwise equivalence tests. Equivalence calculations
  are unchanged.

* G-study variance components now retain full precision, so changing score
  units does not erase small components before G/Phi are calculated. Facet-level counts
  reflect rows actually used by the mixed model; numerical warnings remain
  visible as reasons to review the result.

* D-study projections no longer substitute zero for missing or invalid error
  components. Affected coefficients remain unavailable. Planned facet counts
  must be positive integers, including when supplied as numeric factor labels.
  G/D summaries and plots explain that the coefficients are conditional
  observed-score planning values with variance-estimation uncertainty omitted.
  Recreate older G/D results from the existing MFRM fit; the MFRM itself does
  not need to be re-estimated.

* Shrinkage reports calculate mean shrinkage over eligible estimate/SE pairs,
  preserve unavailable factors, and validate supplied prior SDs. A zero prior
  variance means full pooling at zero; the previous help example stated the
  opposite. Zero plug-in SE after full pooling does not mean perfect precision.
  Plot whiskers are labelled descriptive plug-in bands, with their omitted
  uncertainty retained in returned tables. Reapply shrinkage with the original
  settings and regenerate saved reports, plots and exports; no MFRM refit is
  needed. Original estimates and predictions are unchanged.

* Reporting no longer recommends shrinkage solely because facet counts are
  small, or calls estimates stable solely because counts are large. Methods
  text explains the conditional adjustment and omits internal status codes.

* QC reliability/separation checks now require `separation_facets` to name
  facets whose levels should be distinguished. They are not assessed by
  default: similar rater severity is not automatically a quality failure.
  Unrequested checks do not affect the overall verdict. QC flags remain
  screening rules, not evidence of statistical validity.

* QC no longer replaces unavailable fit, element-misfit, unexpected-response
  or connectivity results with passing values. It verifies that supplied
  diagnostics match the fit, checks threshold order within each ladder, and
  avoids assuming that the first facet is a rater. Recreate saved QC results
  with `run_qc_pipeline()` using existing fits and current diagnostics; no
  model refit is required.

* Imported measurement diagnostics now retain finite-estimate/SE counts and
  withhold reliability when a required SE is missing. They no longer create
  joint facet chi-square statistics or p-values from imported marginal SEs
  alone. Source uncertainty conventions are retained. Re-import
  the existing source-package fit with `compute_fit = TRUE` to update saved
  diagnostics. Native QC requires a native fit with response-level diagnostics.

* Printed drift and linking reviews explain omitted source-fit, offset and
  cross-fit uncertainty, and identify `Offset_SD` as residual spread rather
  than an offset SE. They present review guidance without internal status
  codes. Linking calculations and flags are unchanged.

* Reliability and separation tables now show how many finite estimates and
  corresponding SEs were used. Non-finite estimates are excluded from both
  components; indices are unavailable when a finite estimate lacks its SE.
  QC no longer ignores an unavailable facet when taking the minimum index,
  and precision checks request review when no comparable values are available.
  Fit-adjusted values are described as companion indices, not confidence bounds.
  Recreate saved diagnostics with `diagnose_mfrm(fit)` and regenerate dependent
  reports and exports; no model refit is needed.

* Weighting and model-choice summaries distinguish numerical convergence from
  eligibility for inference and explain restrictions in plain language.
  GPCM ranking and PCM/GPCM tests now use separate eligibility checks. Rebuild saved reviews from the existing fits to obtain the
  revised explanations; older saved reviews request this recreation instead of
  displaying outdated guidance. Machine-readable status fields remain in the
  returned objects.

* JML fit summaries now state that fitting applies neither an extreme-score
  adjustment nor a finite-item bias correction. Free extreme Person estimates
  remain infinite; fixed anchors retain their values and JML uncertainty
  remains exploratory.

* Adjusted-score tables retain the original `PrimaryMeasure`, a readable
  `MeasureBasis`, and the exact `ExtremeAdjustment` amount alongside the
  displayed `Measure`, including in CSV exports. `xtreme` changes display
  values only. Measure SEs are now unavailable on replaced rows because the
  original SE does not describe the replacement. Fitting is unchanged.
  A display replacement is unavailable when its required Person-mean
  reference is unbounded.

* Non-Person FairM is unavailable when a free JML Person estimate is infinite.
  Its mean Person reference previously used finite optimizer traces, which
  do not define the mean of the primary estimates. FairZ uses a zero reference
  and is unchanged. Tables and CSV exports record `FairMReference`; summaries
  and plots explain unavailable values. Recompute older diagnostics and
  recreate saved fair-score tables from the existing fit and original options
  before summary or plotting; no model refit is needed.

* Facet SEs from a regularized Hessian or an observation-table fallback are
  now distinguished from ordinary model-based SEs and cannot authorize
  confidence-interval reporting. Their numerical bands remain available for
  diagnostic review. Recompute older diagnostics with `diagnose_mfrm(fit)` and
  recreate precision reports from the existing fit; no model refit is needed.

* Precision reviews distinguish numerical convergence from support for
  inference. JML guidance no longer implies that changing to MML alone makes
  uncertainty valid. Model comparisons explain withheld rankings and tests in
  plain language, with detailed status codes retained in the returned tables.
  GPCM rankings and PCM/GPCM tests use the separate comparison checks
  described above.

* Residual comparisons from `analyze_dff()` / `analyze_dif()` now report group
  differences in observed-minus-expected scores without SEs, p-values,
  confidence intervals or positive/negative classifications. These differences
  do not isolate differential functioning. Existing test columns remain as
  `NA` for compatibility. `plot(dif, type = "summary")` refuses residual confidence
  intervals. Recreate older residual results, summaries and reports using
  the existing fit and original data; no model refit is needed.

* `dif_interaction_table()` retains observed/expected scores and descriptive
  scaled residuals, with `t`, `df`, p-values and `flag_t` set to `NA`.
  `flag_bias` compares absolute residual means in score units, not logits.
  The retained `p_adjust` and `abs_t_warn` arguments have no effect on residual
  outputs. Group-comparison plots and reports no longer imply a DFF test.

* Residual DIF detection and false-positive rates from
  `evaluate_mfrm_signal_detection()` are unavailable rather than counted as
  zero. Re-summarizing older evaluation objects applies this restriction
  without repeating their simulations; target residual contrasts remain
  available.

* New fitted-object Person intervals and newly created portable calibrations
  invert the continuous posterior CDF. Selecting interval endpoints directly
  from quadrature nodes could give incorrect posterior mass even with accurate
  EAP and SD. Scoring algorithm v2 identifies this change; saved v1 calibrations
  retain their original grid-interval calculation and an explicit note. EAP,
  SD, fitting and grid-based plausible-value draws are unchanged. Intervals
  remain conditional on the point calibration and exclude its uncertainty.

* Estimated-population individual scoring now requires explicit
  `readiness_policy = "review"`; an active population model no longer bypasses
  incomplete source-identification checks. The same guard applies to plausible
  values, whose notes now retain the source review restriction. Earlier
  population predictions claiming scoring readiness must be regenerated before
  summary, printing or structured export. Numerical calculations are unchanged;
  local rank and optimizer convergence do not establish population-scoring validity.

* Non-unit observation weights now prevent ordinary-inference approval in the
  shared fit-readiness record. Diagnostic SE/normal-band calculations remain
  available, but do not license formal inference or facet-equivalence decisions.
  Explicit all-unit weights preserve the unweighted estimation/inference path.
  Older saved fits and diagnostics must be refitted or rechecked, and old
  equivalence results recomputed, before using their inferential conclusions.

* Facet equivalence now requires an inference-ready MML fit, matching supplied
  diagnostics, and unregularized covariance with estimable joint contrasts.
  Contrast rank is checked from the model constraints so covariance roundoff
  cannot admit a singular anchored comparison.
  Pairwise TOST SEs include covariance between facet effects. The heterogeneity
  summary uses a joint Wald test; the former unsupported Bayes-factor heuristic
  is unavailable. Forest plots show deviations from the equally weighted facet
  mean with uncertainty in that mean included. Old equivalence bundles must be
  recomputed before summary, printing, or plotting. These corrections do not
  establish finite-sample coverage or extend GPCM/JML inference support.

* Anchor candidate export and baseline reuse now fail closed unless the source
  fit is inference-ready under the current readiness contract. Explicit
  `readiness_policy = "review"` extraction remains available for inspection,
  while anchor application and anchor CSV export retain the strict default.

* Design recommendations now require connected Person-rater and
  Person-criterion assignments by default (`require_connected = TRUE`).
  Evaluation records component counts before fitting, and the screen includes
  failed runs. Disconnected or unassessed designs cannot pass unless the
  requirement is explicitly disabled; saved objects without component counts
  are unassessed. Sparse overlap reviews are shown separately and now include
  failed runs. Indirect connections remain valid: a pair without common
  persons is not automatically a disconnected design, and connectedness alone
  does not establish model identification or sufficient precision.

* Fit summaries, precision reviews, APA tables, Wright/pathway plots, visual
  reports and fit-level export helpers now check that supplied diagnostics
  belong to the fit and its current readiness state. Previously, MML diagnostics
  could incorrectly give a JML summary a positive formal-inference decision.
  Custom table captions/notes cannot bypass the check. Matching saved
  diagnostics remain reusable; recompute diagnostics after refitting.

## Estimation and numerical accuracy

* Fixed-grid, fixed-standard-normal RSM/PCM MML fits with up to 64 free
  parameters can use a guarded curvature restart when ordinary optimizer
  polishing stalls. Convergence and gradient tolerances are unchanged;
  failed proposals retain the previous solution and their review history.
  Refit an affected saved analysis to recompute its numerical status.

* `fit_mfrm(mml_integration = "adaptive")` uses posterior-mode/curvature
  Gauss-Hermite integration with moving-node gradients for direct MML.
  Likelihood, covariance, person scoring, refit sensitivity and replay retain
  the integration mode; supported portable calibrations store an explicit
  adaptive scoring identity. Fixed integration remains the default. Adaptive
  integration is unavailable for EM estimation and its checkpoints. Some
  identification and boundary checks are also unavailable with adaptive
  integration; numerical convergence alone does not establish valid inference.

* Quadrature sensitivity, fitted-object prediction, and portable-calibration
  scoring accept optional `adaptive_quad_points`, comparing the fixed grid
  with posterior-mode/curvature-adapted grids at unchanged parameters and
  prior. Unrounded per-Person log-marginal, EAP and posterior-SD differences,
  adaptive-order changes, and failure reasons remain available for review.
  This diagnostic leaves estimates, intervals, draws, and readiness unchanged.

* Gauss-Hermite weights now use a scaled Hermite recurrence, preserving tiny
  positive weights lost by the eigensolver at high orders. Fitting, posterior
  scoring, and portable calibration share the corrected rule. New artifacts
  identify the revised algorithm; existing valid calibration grids remain
  readable. Orders whose weights cannot all be represented as positive finite
  doubles fail explicitly instead of silently dropping quadrature nodes.

* Added `mml_quadrature_sensitivity()` for explicit same-data RSM, PCM, and
  bounded-GPCM refits across user-selected integration grids. It reports
  continuous changes in likelihood, measurement coordinates, probabilities,
  EAP, and posterior SD without assigning a stability cutoff; the existing
  GPCM-specific function remains available.

* Fitted-object scoring now uses an explicit scoring quadrature grid rather
  than inheriting the fit-time grid. It refuses a one-point grid, rejects
  invalid weights, and fails closed when the source is not scoring-ready unless
  the user explicitly requests a labelled review-only calculation.

* Continuous RSM/PCM Person intervals reuse the existing compiled likelihood
  kernel when enabled. CDF evaluations omit unused category probabilities;
  integration tolerances and the moment/gradient calculations are unchanged.

* Direct MML optimization now rejects trial points whose population variance
  overflows or underflows after log-scale transformation. This fixes an abort
  in a five-category GPCM fit with all-maximum responses. L-BFGS-B can also
  request gradients at rejected slope/variance trials; those receive the
  constant penalty's gradient. Readiness still uses the actual terminal
  gradient, and invalid starting or retained parameters remain errors.

* MML `print()` and `summary()` output now states the optimization engine,
  fixed or adaptive Gauss--Hermite rule, quadrature order, one-dimensional
  latent structure, and the fitted population-scale identification.

* Automatic nesting review now checks that the population specification is
  shared. Adding a fixed interaction no longer passes this review when the
  population design also changes or cannot be verified. Person and column
  permutations are aligned, and nuisance estimates may differ between fits.

* Model-comparison warnings report the affected fits and recorded inference
  readiness reasons, distinguishing estimability review from numerical
  convergence. Likelihoods, comparison eligibility and LRT calculations are
  unchanged.

* MML EM checkpoints now bind the data, model, parameter layout, quadrature,
  package version, and engine stage. Incompatible checkpoints and completed
  pure-EM re-entry are refused, hybrid warm-start checkpoints are honored, and
  writes use checked same-directory replacement.

## Diagnostics, linking and design

* Design recommendations now require results for every requested facet in
  each candidate design. Previously, a design missing a requested facet
  could pass on its remaining facet alone. `FacetsMissing` identifies absent
  facets, and `FacetsRequired` counts the request rather than available rows.
  Default requests use stored facet names when available, so missing results
  no longer silently narrow the requested scope.

* `recommend_mfrm_design()` accepts `max_ratings` and
  `max_ratings_per_rater` to select designs within total and individual rater
  workload limits. Counts use actual generated rating rows, including linking
  assignments and failed fits, and limits apply to the maximum across all
  recorded replications. Missing workload records cannot pass an active limit;
  omitting the limits preserves the existing ranking behavior.

* Design-evaluation summaries, plots, and recommendations now include failed
  fits and diagnostic runs in the convergence-rate denominator. Previously,
  one successful run and nine failures could be reported as 100% convergence
  and pass a design recommendation. `Reps` now counts all recorded runs;
  `AvailableReps` counts returned facet results. Re-summarize an existing
  evaluation object to apply the correction without refitting.

* Design-evaluation summaries now include the appropriate convergence and
  performance notes without missing-column warnings. Previously, these notes
  could be omitted even when the run summaries contained the required values.

* `evaluate_mfrm_design(parallel = "future")` now dispatches replications
  through the active future plan while preserving the same preallocated
  design-replication seeds and stochastic inputs as serial execution.

* Fair-average tables and plots now also reject diagnostics from a different
  fitted analysis. This prevents mixing one fit's thresholds with another
  fit's measures and standard errors when computing FairM/FairZ.

* Fair-average summaries now use FairZ values and SEs when `reference = "zero"`,
  with an explicit `FairMetric`. Requested fair-score intervals carry separate
  diagnostic-only eligibility metadata through tables, summaries and plots.
  GPCM fair-score SEs remain unavailable if any gradient component is missing,
  rather than replacing the missing uncertainty with zero.

* `plot(fit, type = "fair_average", view = "measure")` connects person/facet measures
  to Fair Scores alongside the observed-score and gap views. Base/ggplot
  share colours and six point shapes, grayscale, and title/note controls;
  full notes and excluded rows remain available in the returned data.
  FairZ documentation now correctly describes a zero-reference expected
  score, not a z-score. Conditional RSM/PCM interval derivatives now use the
  actual FairM/FairZ profile and undo measure-unit scaling. GPCM bundle
  intervals honor the requested confidence level. These approximate intervals
  do not establish full-refit coverage; gap whiskers hold observed means fixed.

* Unexpected-response counts and percentages now use all flagged observations
  before applying the table's `top_n` limit. This also corrects default
  diagnostic summaries and bias-adjustment comparisons, where truncation
  could previously understate prevalence or overstate the reduction in flags.

* Generalizability-study and D-study results now state explicitly that their
  G/Phi coefficients are estimated on the observed numeric score scale, not on
  the fitted MFRM latent scale. Ordered scores continue to use the documented
  Gaussian linear mixed-model approximation in this complementary analysis.

* Corrected the G/Phi interpretation example: under the implemented
  nonnegative main-effect decomposition, `Phi <= G`; `Phi < G`, not `G < Phi`,
  indicates that absolute decisions carry additional facet-main-effect error.

* Linking-chain plots add `type = "offset_sensitivity"`. Each common-element
  deletion reruns preliminary offsets, screening, weighted/unweighted offsets
  and cumulative offsets with source estimates/SEs fixed. Returned tables expose
  signed changes, rescreened elements, actual weight contributors, screening
  fallback and unavailable comparisons. Missing offsets propagate rather than
  becoming zero changes; plots distinguish complete, partial and unavailable
  comparisons. Inconsistent source offsets require rebuilding the chain. No
  source model is refitted and no new uncertainty interval is calculated.

* Linking results record contributors and whether screening falls back to all
  finite rows. Non-finite preliminary weighted offsets now return an unavailable
  result rather than failing during the subsequent screening test.

* Linking chains add `type = "links"` and `type = "anchor_removal"` plot views.
  They summarize recorded retained connections and remove each common-element
  identity across all comparisons, holding screening fixed. Returned tables
  distinguish lost direct links from newly disconnected wave pairs, preserve
  existing disconnections and unknown-retention information, and retain all
  waves. Monochrome line/point encodings and title/note controls are supported.
  These are graph calculations; offsets, estimates and SEs are not recomputed,
  and zero new disconnections does not establish negligible statistical impact.

* The screened linking-chain graph retains isolated waves and handles no-common-
  element cases. Elements are shown separately for each reviewed link so that
  retained/excluded decisions cannot be merged across comparisons. Returned
  node/edge tables preserve source IDs, screening flags and interpretation notes;
  explicit IDs handle repeated names and punctuation. Four line types distinguish
  retention states in colour and monochrome. `show_title` and `show_notes` allow
  clean figures while notes remain printable. Drawing now uses base graphics
  without requiring `igraph`. This graph does not establish anchor invariance,
  precision or common-scale comparability.

* Facet-level category avoidance is now documented separately from global
  score support and from GPCM model choice. The legacy facet-dashboard
  `CentralTendencyFlag` is disabled by default because proximity of a severity
  estimate to the fitted origin does not diagnose restricted category use;
  `data_quality_report()` remains the response-based screening route.

* Adjusted-score summaries and plot notes explain the interval assumptions
  in plain language. Interval values and uncertainty metadata in extracted
  tables are unchanged. These approximate intervals remain unsuitable for
  confidence-interval decisions.

## Plots and reporting

* Report narratives now describe numerical convergence separately from formal
  inference readiness. A converged fit awaiting statistical review no longer
  has its convergence described as unknown.

* `plot(reference, comparison = comparison)` adds paired Wright distributions/locations and CCCs,
  plus matched-location and probability-difference views. Existing source
  coordinates and probability calculations are reused, with explicit group
  selection, original score labels, source readiness and returned notes.
  Monochrome/many-category CCCs use category panels by default. Incompatible
  recorded comparison bases fail without automatic alignment. Drawing uses
  optional ggplot2; draw-free data need no renderer. No difference SEs or CIs
  are calculated. Wright maps also handle a single person without applying the
  undefined one-observation FD histogram rule.

* Rater-severity profiles reserve space for complete labels and explain when
  the device is too narrow. FACETS-style maps warn about crowded facet labels;
  visual help explains font selection for non-Latin text.

* Rater-severity profiles omit band descriptions and legend entries when
  `show_bands = FALSE`. Band guidance is descriptive and no longer implies
  operational interchangeability or a training decision.

* Category-count plots in rating-scale, category-structure, and QC displays
  now include expected counts in the vertical range, preventing the expected
  line from being clipped when it exceeds the tallest observed bar.

* Fit-family Wright, pathway, and CCC plots now share CUD-informed series
  colours across base and ggplot. Bright yellow is omitted from line defaults;
  small data labels use neutral dark text. Expected-score/CCC curves vary line
  type even in colour, CCC overlays vary point shape, and Wright subgroup
  densities vary line type. ggplot conversion preserves category order, custom
  palettes, and monochrome presets. Large series sets still need suitable
  labels, panels, and canvas size; a palette alone cannot ensure readability.

* Fit plots accept `show_title` and `show_notes` (both default to `TRUE`).
  Users can omit figure titles and explanatory annotations while retaining
  axes, legends, and structural panel headings. Interpretation, curve-profile,
  and retention/display notes remain in `data$notes` and are printed by
  `print()`; ggplot conversions honor the flags and retain a `mfrmr_notes`
  attribute. Readiness and R warnings are unchanged by these display flags.

* Clean native pathways reserve enough bottom space for their axis title.
  FACETS-style Wright maps retain every column heading, split scale headings
  over two lines, and wrap facet cells and footnotes to the available width. Crowded headings and
  overflowing frequency stars warn about the needed width or star-count setting.

* Expected-score pathways now display the same reference-profile condition as
  CCC plots: fitted steps/slopes are retained and additive facet main effects
  and interactions are fixed at zero. Their returned `curve_basis` tables
  match, including in plot bundles; numerical curves are unchanged.

* ggplot conversions wrap titles, subtitles, and captions at 72 text columns
  and omit missing labels, preventing ordinary-width subtitle clipping and
  spurious `NA` captions. Users can override labels with `ggplot2::labs()`.
  Multi-group native CCC plots now share a margin legend even with five or
  fewer categories, keeping repeated keys off the fitted curves.

* Native Wright `top_n` now limits non-person facet locations independently
  of step thresholds. Many PCM thresholds no longer consume the limit and
  silently remove every item/rater location; returned locations and retention
  counts reflect this corrected selection.

* Plot label abbreviations now respect display width and distinguish long
  common-prefix names, including Japanese. Crowded native Wright/pathway
  labels prioritize fitted-point visibility and warn when more space is
  needed. Wider Wright maps keep their legend inside the page.

* Native CCC plots use a shared margin legend beyond five categories and
  distinct default colours beyond eight categories. Long panel names remain
  distinguishable, and single-group legends stay within the caller's panel.

* Native Wright and expected-score pathway labels now account for text width,
  text height, nearby labels, and fitted points when choosing their positions.
  Leader lines retain the link to each unchanged estimate; off-screen pathway
  thresholds stay off-screen. Drawing moves text without changing fitted
  coordinates; initial annotation positions remain available to custom renderers.

* Plot styling now restores only the graphical settings it changes, allowing
  single-panel plots to advance through a caller's grid. Native Wright maps
  and multi-group CCC plots clean up their own layouts before subsequent plots.
  Wright titles, pathway endpoint labels, and plot footnotes have more room;
  subgroup densities now draw in the Wright map's person panel.
  `plot(fit, type = "apa_summary")` uses a single-panel FACETS-style Wright map without
  CI whiskers (also reflected in its returned Wright payload), so all four
  panels share one page; its summary text wraps within the panel.

* Weighted APA reports distinguish retained rating rows from the sum of
  observation weights. Visual count summaries use retained rows as well.

* APA design text now reports each facet's own level count. APA and visual
  reporting reuse stored residual PCA results without silently computing
  omitted overall or facet-specific analyses.

## Reproducibility and documentation

* Replay CSV imports preserve literal identifiers (including leading zeros and
  `"NA"`), non-syntactic column names, and latent-regression factor levels and
  contrasts. Background-data paths also work with spaces in command-line runs.

* Diagnostics retain their replay settings, including fit standardization,
  interaction selection, and PCA limits. FACETS-workflow replay uses supported
  arguments and preserves the requested MML engine.

* Identifier text is normalized to UTF-8 during data preparation so native-encoded
  non-Latin labels do not fail in marginal diagnostic sorting.

* Results and summary-table archives omit tables with no column definitions instead of writing
  empty CSVs that `read.csv()` cannot read. Defined zero-row tables retain
  their column headers.

* Replay scripts with person-data CSV files now locate those files beside the
  replay script when invoked through nested `source()` or `sys.source()`.
  They no longer mistake the outer Rscript driver for the replay file.

* Replay scripts now preserve facet-interaction specifications and their
  support policy.

* Manuscript examples now pass the same diagnostics to tables and Wright maps,
  save full-precision CSVs with separate captions/notes and a PNG, and include
  input data for fit replay. Results-archive guidance distinguishes rebuilding
  results from refitting the original analysis.

* Introductory help examples, README, and the workflow vignette now start with
  loading data, fitting with default MML settings, `plot(fit)`, and
  `summary(fit)`. Overview tables are distinguished from individual estimates;
  diagnostic and reporting routes follow the basic example.

* Beginner examples now show the input rows, explain saved summaries, and
  select person, rater, and criterion estimates explicitly. Follow-up examples
  use the same default MML setup, draw the intended figures, distinguish basic
  summaries from comprehensive results, and show where exported files are saved.

* Follow-up examples now connect default MML fits to model comparison,
  diagnostic plots, score summaries, APA tables, and CSV export. Reporting
  checklist examples use logical flags and describe draft availability rather
  than permission to paste text unchanged. Agreement-plot guidance now matches
  the observed-score difference direction and model-expected agreement baseline.

* README and the workflow vignette now connect the packaged example to a user's
  own CSV: rating-row layout, identifier preservation, column mapping, declared
  score categories, missing-data review, and interpretation of individual
  estimates. The vignette adds a runnable CSV round trip, wide-to-long example,
  and input troubleshooting before the advanced workflow.

* Data-check and missing-code help now distinguish input missingness,
  replacement counts, retained rows, and full versus summary table names.
  Examples show which data to reuse after cleaning. Guidance on extreme scores
  now distinguishes MML posterior person scoring from fixed non-person facets.

* The workflow and reporting vignettes reuse one fit for actual diagnostic
  tables and figures. A manuscript coverage map connects reporting topics to
  package outputs and author-supplied study information. Worked reporting
  examples show estimates with uncertainty and distinguish severity, fit,
  separation reliability, and agreement; the example rubric is correctly 1–4.

* Reporting examples now include a short Methods and Results write-up,
  residual denominators, MML step SEs and intervals from existing diagnostics,
  and observed-versus-fair scores with an explicit reference. Guidance based
  on published reporting examples distinguishes separation, strata, and
  posterior-variance EAP reliability and identifies the relevant score units.

* Front-door and first diagnostic help examples now follow one explicit,
  inference-ready workflow: the applied synthetic data, named column roles and
  score support, RSM-MML fitting, fit review, Wright map, diagnostic summary,
  data/rating-scale/precision review, results, and report/export handoff. The
  historical `profile = "facets"` name no longer assumes experience with
  FACETS, TAM, or sirt, and examples no longer choose JML merely for speed.



## Connected assessment-planning reports and residual exploration

- `mfrm_report()` accepts G/D-study, paired D-study comparison and residual
  network objects directly. `build_summary_table_bundle()`, `apa_table()` and
  `export_summary_appendix()` retain their full tables, settings and cautions.
  Saved G/D intervals are checked against their source; no intervals are
  calculated while reporting. New component intervals retain source context;
  recreate older context-free intervals before attaching them.
- `templates/gtheory-report.Rmd` connects observed design and fixed weights to
  components, G/Phi, error variances, SEM, saved intervals, paired comparisons,
  three figures and CSV/HTML/RDS. Rubric direction, sampling universe, decision
  criteria and costs remain explicit manuscript prompts.
- `mfrm_residual_network()` constructs signed Pearson correlations of
  Person-by-facet-level standardized-residual means. Optional average/complete
  hierarchical clustering uses common-person distances and a supplied group
  count. All pairs and overlap reasons remain visible. `plot()` supplies the
  network, full matrix and available tree; the matrix supports `as_ggplot()`.
  This is descriptive exploration, not EGA, a fitted residual-network model,
  stable latent classes or rater-quality classification. Help distinguishes
  unsupervised grouping from supervised random-forest follow-up and explains
  leakage, calibration uncertainty and independent validation.

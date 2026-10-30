# FACETS Manual Mapping

This package implements native R estimation and FACETS-style reports.
Formula correspondence, numerical agreement under matched settings, and
inferential validity are separate claims. A local FACETS 4.5.0 comparison of
complete binary data with ordinary JML, unit weights and aligned facet signs
reproduces location estimates, SEs, MnSq, PTMEA and both WHEXACT conventions
within the external report's precision. Bias-cell estimates, baseline expected
scores and the zero-bias summary also agree in that example. This does not
establish equivalence for other models, versions, anchors, extreme treatments,
missingness patterns or replication weights.

A separate five-facet PCM comparison checks Criterion-specific 0:3 scales,
unit-weight nonextreme observations, fair averages and core fit statistics.
Independent response-probability calculations verify the native fair-score
aggregation and conditional plot derivatives. Small calibration differences
can move a fair average across a two-decimal rounding boundary, so identical
printed values are not the numerical acceptance rule. Expected correlations
agree after the published finite-count correction described below; four small
Outfit-df discrepancies remain unresolved. These results do not establish
every executable diagnostic as equivalent.
Scorefile `S.E.` belongs to the measure; it is not a fair-average SE. The
native plot's conditional interval holds the other calibration parameters
fixed and remains diagnostic-only, without established coverage.

Reviewed against the 0.2.4 development routes on October 24, 2026. Use
`facets_feature_coverage()` for route limits and
`help("mfrmr_methods_reference")` for calculation assumptions and literature.
These routes are not automatically available for two-family GPCM or corrected
JML merely because they share `mfrm_fit`. Consult `gpcm_capability_matrix()`.
Plotting uses `plot()` on fits or analysis results; see `mfrmr_plot_migration`.

## Corrected calculations and remaining differences

- **ZSTD:** diagnostics now default to the FACETS fourth-moment df convention.
  `whexact = FALSE` implements [WHEXACT=YES](https://www.winsteps.com/facetman64/whexact.htm);
  `TRUE` implements NO. Both use the cube root; NO omits the additive correction.
  `whexact = "legacy_linear"` reproduces the previous TRUE transformation.
  `fit_df_method = "engine"` retains the former df convention, and `"both"`
  exposes sensitivity. Identical MnSq can cross a ZSTD cutoff under different
  transformations or df. Category residual summaries are descriptive; their
  ZSTD and unit-baseline fit flags are unavailable.
- **PTMEA:** scores now correlate with the oriented observation-level measure
  sum, including Person rows, as defined in the
  [point-measure documentation](https://www.winsteps.com/facetman64/gstats.htm).
- **Expected PTMEA:** ordinary additive RSM/PCM JML diagnostics now expose
  `PTMEAExpected`, the covariance of the oriented predictor with the conditional
  expected score divided by the square root of the product of predictor
  variance and the sum of expected-score variance and `(N-2)/N` times mean
  conditional variance.
  All moments use frequency weights and `N` is their sum. This follows
  [Linacre's finite-count note](https://www.rasch.org/rmt/rmt221e.htm), which
  supplements the earlier unadjusted attenuation formula. It explains the
  retained PCM expected-correlation discrepancy; fresh unit, repeated,
  fractional and unequal-weight checks agree within printed precision where
  `N >= 2`. `PTMEAExpectedAsymptotic` retains the unadjusted ratio and status
  columns preserve unavailable scopes. At `N < 2`, mfrmr withholds the corrected
  reference rather than copying native boundary values. Group-mean anchors,
  fitted interactions and other estimation routes are outside this addition.
  The corrected quantity is a reporting approximation, not an exact expectation
  of a finite sample correlation or a calibrated inferential test.
- **Displacement:** the proposed update now follows the fitted facet sign and
  score derivative. For a negative facet it subtracts the predictor shift.
  It remains a conditional one-step diagnostic, not a complete released refit.
- **PCM fair averages:** expectations are evaluated under each administered
  scale and then averaged using response weights, following the
  [fair-average definition](https://www.winsteps.com/facetman64/fairaverage.htm).
  Non-Person reference means use finite Person estimates; eligible/excluded
  counts and the score basis are reported. An empty reference stays unavailable.
- **Facet chi-square:** estimates and positive finite SEs define one eligible
  set for both the statistic and df. Excluded levels are counted explicitly.
- **Bias scores and null:** `Expctd Score` and `Obs-Exp Average` now describe
  the baseline main analysis. Separate after-bias columns retain adjusted
  residuals. The summary tests all bias cells against zero with df equal to
  the eligible cell count, rather than centering on a fitted common bias.
- **Bias t df:** the [Table 13 manual](https://winsteps.com/facetman64/table13.htm)
  says count minus 2, but the executed FACETS 4.5.0 example reports df 7 for
  cells with count 8. mfrmr retains its count-minus-1 approximation, matching
  that executable example. Do not infer another version's behavior from the
  web text alone; these conditional tail areas do not establish calibrated
  small-sample fairness inference. In fixed-calibration fractional-count
  examples, native tail areas differ from ordinary Student tails beyond
  print precision, count-minus-two and adjacent-integer df alternatives.
  The checked native fractional tails are reproduced by using integer-df
  finite sums with fractional df in their arguments, while truncating the
  recurrence length and parity. This is not the continuous Student distribution;
  at count 13.99, a checked cell has native tail .0589 versus Student .04905446.
  This reproduces the observed behavior, not an inspection of the executable's
  source. mfrmr retains the continuous Student reference, independently checked
  by integrating its density. `ProbabilityReference`, `ProbabilityStatus` and
  `ProbabilityReason` preserve its basis and unavailable states. `DFDisplay`
  rounds positive df to the nearest integer (halves upward) for comparison.
  Nonpositive count-minus-one values stay in `DFUnrounded`, but cannot define a
  reference df or probability. These fields do not relabel a separate GPCM LR tail.
- **Anchors and extremes:** a fixed anchor is no longer replaced by an extreme
  display adjustment. Malformed nonempty anchor schemas now always stop fitting.
  The default `xtreme_method = "display"` retains the average-score display
  target. Ordinary additive RSM/PCM JML also supports
  `fair_average_table(..., xtreme = .3, xtreme_method = "facets")`:
  [FACETS Xtreme](https://www.winsteps.com/facetman64/xtremescore.htm)-style
  weighted-total extrapolation, using the element's minimum positive eligible
  weight and actual administered calibration. Primary estimates stay intact;
  corrected locations enter that element's fair scores and conditional
  information supplies a distinctly labeled reporting SE. Fair-score intervals
  remain unavailable for positive `xtreme`. `estimate_bias(..., xtreme = .5)`
  separately handles extreme cells with a root inside the specified search
  limits. Extreme main effects are excluded from its calibration rows.
  Binary unit/replication-weight examples and a fixed-scale PCM kernel agree
  with FACETS 4.5.0 to small numerical tolerances. Public PCM results also satisfy
  their own administered-scale score equations. This does not establish full
  PCM fitting or whole-program equivalence. Matching the origin constraints
  resolves the checked interior partial-anchor discrepancies, as described
  below. With only same-boundary Person-anchor responses in matched sparse
  RSM/PCM inputs, FACETS 4.5.0 retains all observations and the anchor value,
  but reports a finite calibration whose expected anchor total is about
  0.299 or 17.701 for observed totals 0 or 18 (`Xtreme=0.3`). A feasible
  location shift strictly improves the ordinary likelihood, so this is not
  a finite ordinary JML solution. mfrmr's `anchor_origin` route rejects this
  configuration before optimization without requiring `lpSolve`; additional
  interior or opposite-boundary Person-anchor support does not trigger that
  specific rule. With additional interior Person-anchor support, native
  fair scores and conditional Rater-by-Task bias calculations agree with
  independent calculations at the reported calibration. All responses,
  including those of the fixed extreme Person, remain eligible in these
  examples. Main calibrations still differ: free-Person scores are nonzero
  and maximum-score examples reach the main iteration limit. Conditional
  reporting agreement therefore does not establish calibration equivalence.
  Ordinary additive RSM/PCM JML now retains direct fixed anchors in both
  fair-score display conventions, including `totalscore = FALSE` and the
  conditional PCM plot derivative. The exact native adjustment and broader
  anchored exclusion designs remain unverified. Do not extend the manual's general
  extreme-exclusion rule to fixed anchors without checking the actual output.
  Group-mean anchors and fitted interactions are unsupported.
  Reporting SEs and screening tails are not uncertainty for the original
  infinite MLE, calibration-inclusive inference or fairness decisions.
- **Pairwise bias:** local estimates use `target + bias/TargetSign`, with
  conditional cell SEs. The default now uses the modified Rasch-Welch df in
  the [t-statistics manual](https://www.winsteps.com/facetman64/t-statistics.htm),
  rounded down for the tail area as confirmed in an unequal-count FACETS 4.5.0
  comparison. `DFUnrounded` retains the continuous value. Explicit
  `df_method = "classical"` retains the earlier continuous Welch reference.
  Recreate pairwise reports to obtain `CalculationVersion = 3`. Missing inputs
  withhold screens; calibration uncertainty and cross-cell covariance are omitted.
- **Unexpected responses:** the default is now `abs(z) >= 3` with
  `rule = "residual"`, consistent with the [FACETS threshold](https://www.winsteps.com/facetman64/unexpected.htm).
  The former exploratory rule is available through `abs_z_min = 2,
  prob_max = .30, rule = "either"`. Tied modes no longer arbitrarily favor the
  first category: `MostLikely` is missing, `MostLikelyCount` records the tie,
  and category distance uses the nearest mode. Screening is not a fairness test.

Recompute `diagnose_mfrm()` from the saved fit, then recreate `estimate_bias()`
and derived reports. Corrected native bundles carry `calculation_version = 2`;
older bundles are refused when combined with a fit. No main calibration refit
is required for these diagnostic corrections. Retain old outputs as historical
results rather than relabeling their values with the new definitions.

Category `AvgMeasure` now averages the oriented generating predictor;
`ExpectedMeasure` uses category-probability weights over retained observations.
`AvgPersonMeasure` and `ExpectedAverage` remain separately named Person-only
and predicted-score summaries. Grouping residuals by their realized category
does not preserve a universal unit MnSq baseline, so category ZSTD and their
fit flags are unavailable. This is a statistical limit even though
[Table 8.1](https://www.winsteps.com/facetman64/table8_1ratingscale.htm) displays
category Outfit. Counts are native retained-row counts, not automatically the
FACETS Total/Used extreme-selection counts.

Category `MedianThreshold` now denotes the cumulative 50% boundary, with
`AndrichThreshold` retaining the adjacent step. Half-score crossings outside
the requested grid stay unavailable; mode peaks refer to that finite grid.
Recreate old category-structure results before summary or plotting. Fractional
weighted totals retain requested decimal precision, conditional SE fallbacks
scale once, and chi-square summaries preserve unavailable decisions and the
full facet family before display truncation. Centered computations avoid
cancellation under a common location shift. The random chi-square is a spread
index; with equal SEs it cannot distinguish distributional shapes and cannot
establish normality or random sampling.

Fixed-field score imports require validated all-field headings with element
number after the statistics. Appended facet fields no longer enter element
labels; ambiguous, headerless or element-first layouts are rejected. Use named
delimited output for those layouts. Native score exports are observation-level,
whereas FACETS Scorefile imports are element-level. Co-observation graph
components are also distinct from identification subsets: inspect the model's
rank audit even when the graph has one component. Native Tables 5/11 cover
counts and selected summaries, not all FACETS response moments or variance
reports. The historical `conditional_probability` plot alias is ordinary P(k),
not the FACETS adjacent-conditional P(k)/(P(k-1)+P(k)); prefer `ccc` for that
native curve. Native information curves depend on the specified design.

Native MML uses EAP plug-in residual fit, not posterior-averaged TAM/ConQuest
parameter fit. Reporting defaults to all rows (`totalscore = TRUE`), whereas
[FACETS Totalscore](https://www.winsteps.com/facetman64/totalscore.htm) defaults
to excluding responses involving other extreme elements. Raw total/count
columns and weighted totals are distinct; `ObservedAverage` uses the latter.
Match [replication weights](https://www.winsteps.com/facetman64/weighting.htm),
observations and extremes before comparing. See `mfrmr-facets-migration` for
migration and decision boundaries. Agreement on rounded numbers does not
establish agreement at a cutoff or inferential calibration.

## Manual References Used
- Output table index (64-bit): <https://www.winsteps.com/facetman64/outputtableindex.htm>
- Output table index (32-bit archive): <https://www.winsteps.com/facetman/outputtableindex.htm>
- Table 5 measurable data summary: <https://www.winsteps.com/facetman/table5.htm>
- Table 7 measurement report: <https://www.winsteps.com/facetman/table7.htm>
- Table 7 reliability/chi-square notes: <https://www.winsteps.com/facetman/table7summarystatistics.htm>
- Table 7 agreement notes: <https://www.winsteps.com/facetman/table7agreementstatistics.htm>
- Table 8.1 rating scale report: <https://www.winsteps.com/facetman/table8_1ratingscale.htm>
- Table 8 bar-chart report: <https://www.winsteps.com/facetman/table8barchart.htm>
- Table 8 probability-curves report: <https://www.winsteps.com/facetman/table8curves.htm>
- Graph output (`Graphfile=`): <https://www.winsteps.com/facetman/graphoutputfile.htm>
- Table 9 bias-iteration report: <https://www.winsteps.com/facetman/table9.htm>
- Table 10 unexpected-after-bias report: <https://www.winsteps.com/facetman/table10.htm>
- Table 11 bias-calculation counts report: <https://www.winsteps.com/facetman/table11.htm>
- Table 12 DIF/bias summary report: <https://www.winsteps.com/facetman/table12.htm>
- Table 13 DIF/bias detail report: <https://www.winsteps.com/facetman/table13.htm>
- Table 14 pairwise bias report: <https://www.winsteps.com/facetman/table14.htm>
- Score file fixed-field columns: <https://www.winsteps.com/facetman64/scorefileinvisible.htm>

## Implemented (Direct or Close Compatibility Surface)
- Core multifacet estimation (RSM/PCM, MML/JML): `fit_mfrm()` / `mfrm_estimate()`
- Diagnostics core bundle (obs, fit, reliability, interactions, subsets): `diagnose_mfrm()`
- Table 1-style specification summary: `specifications_report()`
- Table 2-style data summary report: `data_quality_report()`
- Table 3-style iteration report (replayed): `estimation_iteration_report()`
- Table 4-style unexpected responses: `unexpected_response_table()`, `plot(fit, type = "unexpected")`
- Table 5-style measurable summary bundle: `measurable_summary_table()`
- Table 6.0.0-style subset/disjoint listing: `subset_connectivity_report()`
- Table 6.2-style facet-statistics graphic summary: `facet_statistics_report()`
- Table 7-style facet/person measures and fit summary: `diagnose_mfrm()` + `summary.mfrm_fit()`
- Table 7 reliability + facet chi-square style summaries: `diagnose_mfrm()$reliability`, `diagnose_mfrm()$facets_chisq`, `facets_chisq_table()`, `plot(fit, type = "facets_chisq")`
- Table 7 agreement style summaries: `diagnose_mfrm()$interrater`, `interrater_agreement_table()`, `rater_network_analysis()`, `rater_halo_network_analysis()`, `plot(fit, type = "agreement")`
- Table 6-style Wright/variable-map visual display:
  `plot(fit, type = "wright", view = "unified", renderer = "facets")` provides the shared
  logit ruler, person-frequency asterisks, signed facet columns, and labeled
  score-transition lines; `renderer = "native"` retains the package-native
  facet SE/CI display. Both routes expose editable draw-free data.
- Fit-oriented pathway display:
  `plot(fit, type = "fit_pathway", fit_stat = "Infit", include_person = TRUE, top_n_person = 12)`
  places Infit/Outfit on the x-axis and measure logits on the y-axis without
  changing the existing expected-score `type = "pathway"` contract.
- Table 8.1-style rating scale bundle: `rating_scale_table()`
- Table 8-style bar-chart and curves exporters: `category_structure_report()`,
  `category_curves_report()`, including cumulative probability, total
  information, and category-specific information curve data
- Output-file emulation (`GRAPH=` / `SCORE=` style): `facets_output_file_bundle()`
- Standalone residual and subset handoff files:
  `write_mfrm_residual_file()` and `write_mfrm_subset_file()`
- FACETS fit/score-file import: `read_facets_fit_table()` /
  `import_facets_fit_table()` for delimited and fixed-field score-file extracts
- Fair-average measurement summaries: `fair_average_table()`, `plot(fit, type = "fair_average")`.
  The column contract now identifies this as Table 7, retaining only the
  legacy object key `t12`. FACETS Table 12 is the
  [bias summary report](https://www.winsteps.com/facetman64/table12.htm).
- Displacement diagnostics (FACETS-style anchor drift check): `displacement_table()`, `plot(fit, type = "displacement")`
- Bias re-estimation iteration (Table 9 workflow): `estimate_bias()` / `estimate_bias_interaction()`
- Table 10-style unexpected-after-bias output: `unexpected_after_bias_table()`
- Table 11-style bias-count report: `bias_count_table()`
- Table 12/13/14 style bias outputs:
  - summary/detail: `estimate_bias()` result (`summary`, `table`, `chi_sq`)
  - higher-order mode (3+ facets): `estimate_bias(..., interaction_facets = c(...))`
  - pairwise: `bias_pairwise_report()` for structured review and
    `build_fixed_reports()` (`pairwise_table`, `pairwise_fixed`) for 2-way
    fixed-report handoff
  - Table 13 plot export: `bias_interaction_report()` → `plot(result, type = "bias")` (including `view = "facet_profile"`)
- Fixed-width report generation: `build_fixed_reports()`
- APA narrative/table helpers: `build_apa_outputs()`, `apa_table()`
- Residual PCA checks (overall + by facet): `analyze_residual_pca()` → `plot(result, type = "scree")` (use `mode = "facet"` and `facet` for a facet-specific view)
- QC dashboard (base graphics): `plot(fit, type = "qc")`
- Anchor workflow: `review_mfrm_anchors()`, `make_anchor_table()`
- Data packaging/loading helpers: `list_mfrmr_data()`, `load_mfrmr_data()`
- FACETS feature coverage boundary: `facets_feature_coverage()`
- FACETS output-contract review:
  `facets_output_contract_review()` checks the documented package columns and
  metrics. A numerical comparison still requires output from a separately run,
  documented FACETS version under aligned settings.

## Capabilities, Limitations, and Alternatives

### Wright and variable maps

- Capability: `plot(fit, type = "wright", view = "unified", renderer = "facets")` provides the shared
  ruler, person-frequency asterisks, signed facet columns, and labelled score
  transitions. The native renderer adds SE/CI information.
- Limitation: visual correspondence does not establish pixel identity or
  numerical equivalence across software versions.
- Alternative: retain the native uncertainty-aware map for analysis. When
  numerical comparison is required, supply aligned external FACETS output to
  `facets_fit_review()`.

### Specification, data, iteration, and measurement tables

- Capability: `specifications_report()`, `data_quality_report()`,
  `estimation_iteration_report()`, `measurable_summary_table()`, and
  `fit_measures_table()` provide structured R tables.
- Limitation: column order, fixed-width layout, and the complete FACETS
  iteration trace are not reproduced.
- Alternative: use the structured tables for analysis and export; use FACETS
  externally when its exact report layout is required.

### Rating-scale tables and curves

- Capability: `rating_scale_table()`, `category_structure_report()`, and
  `category_curves_report()` provide category, threshold, expected-score, and
  information-curve data.
- Limitation: FACETS printer artwork, graph-window behavior, and every legacy
  column option are not reproduced.
- Alternative: use package plot data for custom R graphics, or use FACETS
  externally for its program-specific display.

### Graph, score, residual, and subset handoff

- Capability: `facets_output_file_bundle()`,
  `write_mfrm_residual_file()`, and `write_mfrm_subset_file()` provide
  documented CSV/TSV or fixed-width handoff routes.
- Limitation: these are package-native contracts rather than complete FACETS
  command-file or fixed-field compatibility.
- Alternative: use the package files for reproducible R workflows; use FACETS
  externally when a downstream system requires its exact syntax.

### Bias and pairwise contrasts

- Capability: two-way bias runs provide structured pairwise output through
  `bias_pairwise_report()` and `build_fixed_reports()`.
- Limitation: FACETS layout and options are broader, and higher-order runs omit
  the pairwise section.
- Alternative: report the available structured contrasts and their screening
  boundary; redesign the contrast or use FACETS externally if a missing
  FACETS-specific report is required.

Package bias screens condition on the supplied calibration. Re-estimating a
bias cell does not propagate the joint uncertainty of all calibration
parameters. The accompanying intervals and pairwise summaries are not
full-covariance structural tests.

### Separation, agreement and assessment planning

High rater separation reliability indicates distinguishable severities,
not high observed agreement. Low separation can also result from imprecise
measures. Strata describes hypothetical distinguishable levels, not discovered
groups. See the [FACETS interpretation](https://www.winsteps.com/facetman64/reliability.htm)
and the estimator-specific formulas in `help("diagnose_mfrm")`.

In FACETS 4.5.0, the [G-theory menu](https://www.winsteps.com/facetman64/gtheory.htm)
generates an R program using the external package `gtheory`. The
[version history](https://www.winsteps.com/facgood.htm) records replacement of
that interface in 4.5.1; check the installed version separately from online
help. The menu defaults to all facet main effects, no interactions, observed
scores (`Obs`), and the first facet as the object of measurement. Select the
relevant interactions explicitly before comparing a full G-theory design.
This R route models numeric observed scores with Gaussian random effects;
it does not analyze Rasch parameter uncertainty.

The `gtheory` 0.1.2 calculation uses `lme4` REML with nonnegative variance
components. This differs from unconstrained ANOVA and from simply truncating
its negative estimates. The generated D-study passes the observed data and
score column to calculate averaging denominators. Check that D-study's `n`
column, rather than the G-study table's initial `n = 1`: for a complete
four-rater, four-task design with one rating per Person-by-Rater-by-Task cell,
Person-by-Rater and Person-by-Task use 4, and Rater-by-Task and the residual
use 16. Changing a score variable, interaction or future count changes the
question. Preserve boundary warnings and messages and the R environment
with the results. A message that R was launched does not establish
that its analysis completed; retain the generated script, its data and actual
R output. The archived dependency and temporary input paths also need to be
available when replaying an older FACETS program.

mfrmr's own
`mfrm_g_study()` / `mfrm_d_study()` route handles
supported univariate or multivariate observed-score designs, interactions,
composites and future complete plans. Its component/covariance plots and
D-study curves retain negative components and unavailable metrics by default.
For complete balanced univariate ANOVA sources, the D-study's explicit
`negative_components = "genova_algorithm"` or `"genova_ems"` option reproduces
GENOVA zero substitution or source-design EMS recalculation before future-plan
scaling. Original G-study estimates are retained separately from effective
components; summaries and plot notes record the choice. These options do not
implement REML or multivariate covariance repair. `confint(d)` retains
population-parameter intervals from the original source alongside adjusted
points, with explicit `IntervalEstimate` and `IntervalBasis`; the bounds are
not shifted to center on the adjusted point. Native GENOVA's component SEs
continue to describe the original estimator, even when D-study components are
adjusted. They are not SEs of the truncated/EMS estimator. Supported population
plan differences use `confint(d, reference = ...)`; the separate delta-estimator
comparison requires `"raw"`.
For balanced
ANOVA, `confint(g, assumption = "normal")` supplies component SEs and intervals
for one score or a fixed composite. Two-score covariance intervals now have
explicit generalized and conservative choices with signed bounds. Scalar
component `reml_profile` intervals also accept identifiable MINQUE(0) sources,
retaining original points and separate REML centers but not supplying sampling
SEs. Explicit joint `reml_profile` also supplies signed two-score covariance
intervals for identifiable incomplete sources. This is experimental, has no
sampling SE and does not establish finite-sample coverage or a default method.
Separately, `vcov(g, assumption = "normal")` provides sampling covariance and
SEs for the original ANOVA/MINQUE component estimates, including incomplete
sources. Product-moment correction retains negative or unidentified sampling
variances; it does not describe the separate REML estimator or validate
normal-reference confidence intervals.
For individual G/Phi in balanced two-crossed-facet sources, the default is
approximate MLS null inversion. Undefined MLS widths are conservatively
retained and flagged separately from numerical availability; small-source
intervals may be very wide. This is normal-model approximate inference,
not a transfer of FACETS intervals or a uniform-coverage guarantee.
Balanced crossed G also offers an explicit experimental `score` interval,
using expected information and extended nuisance mean-square expectations.
It preserves original points and distinguishes empty numerical confidence
sets from calculation failures; this does not supply Phi or paired differences.
For G/Phi and SEM, explicit experimental
`reml_profile` accepts balanced or identifiable incomplete sources under
normal effects. Its separate REML fit preserves the original points and uses
approximate likelihood-ratio calibration, without a finite-sample or boundary
coverage guarantee. SEM profiles constrain error variance directly and report
the interval around the REML center without enumerating remote disconnected
accepted regions. With `reference`, one-facet plans or nested plans with a
common future child count have shared-profile G/Phi/SEM difference intervals,
including incomplete sources. These preserve interior extrema, original
points and the base approximation; nonproportional future errors remain
outside this transformation. `mfrm_multivariate_d_compare()` has a separate
normal-theory paired-plan scope. The older `mfrm_generalizability()` /
`mfrm_d_study()` adapter uses a simpler main-effects model and an explicit
residual-scaling rule. These decompositions are not interchangeable. Follow
`vignette("mfrmr-gtheory", package = "mfrmr")`; no MFRM fit is needed for the
multivariate route, and none of these coefficients is MFRM separation reliability.

## Unsupported or External-Only FACETS Surfaces

The public `facets_feature_coverage()` table is the canonical capability
reference. Its `Capability`, `Limitation`, and `Alternative` columns describe
the appropriate route for each surface. External-only areas include:

- FACETS binomial-trial and Poisson-specific reports; use the package's
  two-category ordered-score route only when that model is appropriate.
- Complete FACETS line-printer report emulation; use structured package tables
  or run FACETS for the exact report.
- Arbitrary FACETS R/Web/Excel menu plots; use `plot_data()` with an R graphics
  system or run FACETS for those menu outputs.
- The FACETS/gtheory interface itself; use the package-native observed-score
  G/D routes above for their documented designs.
- FACETS network-menu output; use `subset_connectivity_report()`,
  `mfrm_network_analysis()`, `rater_network_analysis()`, or
  `rater_halo_network_analysis()` for reusable R node and edge tables.
- Winsteps control/data-file export; use Winsteps externally when that format
  is required.
- Arbitrary FACETS command-file and raw-report parsing; provide supported
  delimited or fixed-field fit extracts to `read_facets_fit_table()` instead.

## Anchoring Rules Encoded
- Direct anchors (`Facet`, `Level`, `Anchor`) are fixed.
- Group anchors (`Facet`, `Level`, `Group`, `GroupValue`) constrain group means.
- A level in both tables retains its fixed value and contributes to the group
  mean; incompatible constraints stop fitting.
- Default `centering_scope = "facet"` includes fixed values in the zero sum;
  `"remainder"` retains the earlier remainder constraint. Ordinary additive
  RSM/PCM JML also supports `"anchor_origin"`: direct or group anchors replace their
  facet's centering, matching the FACETS
  [centering specification](https://www.winsteps.com/facetman64/centerfacet.htm).
  Declare an unanchored, active noncenter facet. Convert native per-element
  group contributions to one repeated arithmetic-mean `GroupValue`; native
  group 0 requires direct anchors. mfrmr retains the ordinary constrained
  likelihood, not all native multi-group numerical updates or `X` extreme
  exclusions. For a Person anchor, Noncenter=Rater and centered Task
  resolve the checked sparse RSM/PCM discrepancy in native 4.5.0. Keeping
  Rater and Task centered instead imposes an additional restriction in mfrmr
  and is not a matched FACETS maximum when native score residuals persist.
- Missing `GroupValue` is treated as 0.
- Default review thresholds (not universal adequacy or identification rules):
  - common anchors per linking facet: `>= 5`
  - observations per element: `>= 30`
  - observations per score category: `>= 10`

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
  interpretation until the numerical-readiness gate passes.

* Replaced blanket `\dontrun{}` and `@examplesIf interactive()` guards with
  checkable examples or `\donttest{}` blocks. Only the two workflows that need
  separately generated ConQuest files remain `\dontrun{}`, and only the local
  Shiny viewer remains interactive-only. The release-readiness review now
  enforces that allowlist and flags CRAN-side package workload above ten
  minutes, based on ordinary examples, `donttest` examples, tests, and vignette
  rebuilding. Other top-level check components remain visible as diagnostics
  but do not inflate that package-controlled threshold.

- Added one authoritative repository roadmap and aligned release metadata and
  validation notes with the accepted 0.2.2 boundary. External numerical
  comparison and calibrated MML joint-stationarity gates are explicitly 0.2.3
  work rather than retroactive 0.2.2 requirements.
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
- Added the package hex sticker to the README and pkgdown-standard
  `man/figures/logo.png` location, while retaining the editable SVG source.
- Tightened the FACETS positioning contract against the current 64-bit 4.5.1
  software target: coverage rows describe package-native surfaces, not
  external numerical equivalence, and mixed models, multiple scales,
  threshold anchoring, and fixed-calibration scoring remain outside 0.2.2.
- Corrected the `interrater_agreement_table()` documentation: `ExpectedExact`
  is computed from fitted category-probability vectors, not marginal-frequency
  chance agreement. A focused regression test now guards that definition.
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
- Direct, hybrid, and EM MML engines now apply the same terminal-gradient gate
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
- Per-fit workspaces are local to one optimization and are discarded with the
  fit evaluator. They are not global, are not shared across parallel fits,
  and are not stored as large probability arrays in the returned fit object.
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
  of every bundled synthetic dataset. Fixed-seed generators for the compact
  examples are tracked in the public source repository. Combined-study
  objects now explain that relabeling prevents identifier collisions but does
  not establish a common scale without an explicit anchor/linking design.
- Precomputed vignette tables now follow the same successful operational MML
  route as the displayed workflow and record their source dataset, schema,
  MD5 checksum, and package version.

## Safer first analyses

- The public default remains `reltol = 1e-9` for the initial optimizer stage;
  bounded polishing is invoked only when `reltol <= 1e-9` and code zero
  precedes the terminal-gradient gate. The fitted object records requested and
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
- CRAN checks now exercise the complete introductory workflow once and
  use the exact README/default MML controls rather than a reduced quadrature
  setting. They retain lightweight compatibility/backend/artifact contracts.
  Repeated
  estimation, detailed plotting, simulation, and broad regression coverage
  remain in the complete local and GitHub Actions suite.
- A repository-level first-use workflow stress protocol covers linked, sparse,
  disconnected, shared-link, PCM, bounded-GPCM, extreme-score, separation,
  missing-code, and weighted scenarios across deterministic seeds. It keeps
  expectation matching separate from actual report readiness and is excluded
  from routine CRAN checks.

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
- A matched 31-node run with ConQuest 5.47.5 Demonstration Version is recorded
  in the public source repository's validation record (excluded from the
  installed CRAN package) for the documented binary, item-only, one-covariate
  MML overlap case. The result supports that narrow handoff and is not a claim
  of general numerical equivalence.
- `export_mfrm_results()` now labels every preset as a potentially identifying
  analysis archive, warns before writing unless the risk is explicitly
  acknowledged, and records privacy status in its summary, HTML index, and
  written-files manifest. Fit-level `export_mfrm_bundle()` archives follow the
  same warning and metadata contract, and the lower-level `export_mfrm()` CSV
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

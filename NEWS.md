# mfrmr 0.2.1

This release focuses on a clearer public workflow, a more readable reporting
surface, and source-aware review outputs. The main user-facing entry points are
`mfrm_results()`, `mfrm_report()`, `export_mfrm_results()`,
`launch_mfrmr_viewer()`, and `mfrmr_output_guide("public")`.

Highlights:

- A shorter public workflow now starts from `fit_mfrm()` -> `mfrm_results()` ->
  `mfrm_report()` -> `export_mfrm_results()`.
- `mfrm_report()` and `summary(mfrm_report(...))` provide first-screen report
  readiness tables, report routes, and cautious interpretation boundaries.
- Simulation and recovery summaries expose operating-characteristic,
  sparse-design, peer-review, and recovery-review tables for appendix and
  export workflows.
- Linking, anchor, bias, misfit/pathway, and network review outputs are routed
  through explicit result components and scoped follow-up helpers.
- Release-review artifacts now distinguish CRAN-time lightweight checks from
  full non-CRAN regression evidence, and the readiness gate checks that the
  local `R CMD check` log matches the target package version.

Detailed changes:

- The generated APA narrative no longer labels overall fit as
  "acceptable"/"elevated". It now reports whether mean-square fit fell
  within or outside the active screening band, names the band, and states
  that band position is screening evidence rather than a model-validity
  decision. Element-level wording uses the same screening-band language.
- The APA narrative now separates separation reliability from observed
  inter-rater agreement explicitly: reliability sentences state that
  separation indices are not inter-rater agreement, and agreement sentences
  are introduced as a separate quantity.
- APA Methods text for `MML` fits now states the estimation basis: person
  measures are EAP estimates, and residual-based fit statistics are
  evaluated at those EAP measures rather than at JMLE estimates. A matching
  fit-basis caution routes external FACETS comparisons to
  `method = "JML"`.
- The APA narrative and Table 1 note now report measure-level confidence
  intervals when available: the Results fit/precision section states the CI
  level, method, and `CIEligible` counts, and the facet-summary table note
  instructs reporting `CI_Lower` / `CI_Upper` for eligible rows.
- APA Wright-map and facet-summary notes now state the fitted sign
  convention explicitly (higher person values = higher ability; higher
  non-person facet values = greater severity/difficulty under the default
  negative orientation), including any `positive_facets` exceptions, instead
  of the previous "depending on facet orientation" wording. The pathway-map
  note and `plot.mfrm_fit()` help now state that the expected-score pathway
  display is distinct from the Bond-and-Fox-style measure-versus-fit bubble
  chart available via `plot_bubble()`.
- The MML-vs-FACETS residual basis is now documented across the fit
  surfaces: `diagnose_mfrm()`, `facets_fit_df_guide()`,
  `facets_fit_review()` (including a new residual-basis guidance row and
  interpretation-guide rows), the package statistical-background help, and
  the FACETS migration vignette all state that MML fit statistics are
  computed at shrunken EAP person measures and route JMLE-style comparisons
  to `method = "JML"`.
- Small-df ZSTD availability is now documented as an explicit boundary:
  mfrmr withholds ZSTD as `NA` when the applicable df falls below 1
  (Wilson-Hilferty instability), while FACETS/Winsteps under `WHEXACT` can
  report a value on the same sparse cell. The fit guides and
  `facets_fit_review()` now describe NA-vs-finite ZSTD pairs as availability
  differences rather than fit differences.
- The MML person-row separation/reliability convention is now documented:
  person rows apply the separation formulas to EAP measures with posterior
  SDs, which yields a conservative summary that is lower than the IRT
  empirical-reliability convention and is not numerically comparable to
  JMLE-based FACETS person reliability. The APA reliability sentence carries
  the same note when a Person row is present under `MML`.
- `preset = "monochrome"` is now accepted by all preset-based plot helpers,
  including `plot_fair_average()`, `plot_displacement()`,
  `plot_unexpected()`, `plot_interrater_agreement()`, `plot_marginal_fit()`,
  `plot_marginal_pairwise()`, `plot_facets_chisq()`, `plot_qc_dashboard()`,
  `plot_threshold_ladder()`, and the network/secondary plot methods, not
  only the fit-family plots.
- `assess_mfrm_recovery()` now reports uncertainty evidence separately from
  recovery metrics. This makes unavailable SE/coverage evidence visible without
  turning it into an implied RMSE or bias failure.
- `inst/validation/release-readiness.R` now follows the target package version
  from `DESCRIPTION`, selects versioned evidence-map/checklist files when
  available, and reports stale `R CMD check` logs as package-check review
  items when the check-log package version does not match the target release.
- The release-readiness GPCM scope review now also works against an installed
  package: when the `R/` capability-matrix source file is not present (for
  example inside `R CMD check` runs or CI artifact reviews), it reads
  `gpcm_capability_matrix()` and `gpcm_runtime_guard_coverage()` from the
  installed namespace instead of reporting a missing-source concern.
- `precision_review_report()` now includes a source-grounded fit/separation
  basis table. This keeps mean-square fit, ZSTD standardization,
  Rasch/FACETS-style separation, and package QC thresholds in separate
  reporting lanes.
- New `mfrm_results()` gives users a FACETS-style first-screen entry point:
  it accepts an existing fit, a `run_mfrm_facets()` object, or a standard
  long-format data frame, runs diagnostics automatically, gathers available
  tables/reviews/plot routes, and can emit a lightweight temporary HTML
  report. Its summary includes next-action routes and a replay-code scaffold.
  Existing table, report, analysis, and review helpers remain supported as
  detailed components. `mfrm_results_interactive()` adds an explicit opt-in
  column-selection wizard for interactive sessions only.
- `summary(mfrm_results(...))` now includes a `triage` table that orders
  unavailable, review, information, and OK signals across diagnostics, plots,
  tables, precision/reliability, reporting, model scope, and network review
  surfaces. `next_actions` uses that triage layer for first-screen routing.
- `mfrmr_output_guide("public")` now gives the shortest top-level public API
  map: explicit fit, comprehensive results, report readiness, optional viewer,
  download export, scoped guide, and opt-in interactive routes. The guide also
  carries `APILayer`, `ObjectRole`, and `DecisionBoundary` columns so
  top-level public surfaces, specialist follow-ups, advanced design-review
  rows, and migration/integration routes are separated before users scan the
  broader namespace, and so users can see whether a route estimates,
  summarizes, displays, exports, or only points to the next helper.
- `mfrmr_output_guide("entry")` still gives the first-screen creation routes,
  including explicit-fit, comprehensive-results, purpose-specific guide, and
  opt-in interactive routes. The guide labels lifecycle, user level,
  recommended-entry status, and includes advanced simulation/network route
  rows.
- The README now starts from a shorter public-surface workflow:
  `fit_mfrm()` -> `mfrm_results()` -> `mfrm_report()` ->
  `export_mfrm_results()`. Detailed table, report, review, and export helpers
  are presented as scoped follow-ups rather than as functions to memorize
  before the first analysis.
- `mfrmr_output_guide("viewer")` now maps local Shiny viewer workflows back
  to the `mfrm_results(include = ...)` object that should be created first,
  including publication, validation, bias-screen, pathway/misfit, and combined
  review routes.
- `plot.mfrm_diagnostic_screening()` now gives diagnostic-screening simulation
  output an integrated visualization route. The default overview combines
  legacy ZSTD, strict marginal, strict pairwise, strict combined, and optional
  report-review rates; focused views cover report signals, scenario contrasts,
  and runtime summaries. The draw-free return is an `mfrm_plot_data` object,
  so `plot_data(diag_eval, type = "overview", component = "plot_long")` can be
  used directly for ggplot2, plotly, Quarto, or custom export workflows.
  Draw-free diagnostic-screening plot objects also carry `overview`,
  `reading_order`, `next_actions`, `reporting_notes`, and `figure_recipes`,
  keeping custom figure/report handoffs aligned with the same interpretation
  boundaries used by `summary(diag_eval)`.
- `summary(evaluate_mfrm_diagnostic_screening(...))` now provides an explicit
  diagnostic-screening report surface, and `build_summary_table_bundle()` /
  `export_summary_appendix()` can export its scenario, performance,
  report-signal, contrast, and draw-free plot-data tables.
- Diagnostic-screening summaries now include `reading_order`, `next_actions`,
  `reporting_notes`, and `figure_recipes` tables so users can separate
  first-read tables, follow-up actions, figure/caption planning,
  appendix-only tables, operating-characteristic claims, and validation-gate
  boundaries before writing report text.
- `mfrmr_output_guide("simulation")` now routes users to diagnostic-screening
  evaluation and appendix export alongside data generation, design/recovery
  evaluation, network review, and peer-review design review.
- `mfrm_results(include = ...)` now accepts purpose presets:
  `"publication"`, `"validation"`, `"facets"`, `"bias"`, `"misfit_review"`,
  `"network"`, and `"gpcm_review"` in addition to `"standard"` and `"all"`.
  The bias preset surfaces facet-level bias-screen guidance while leaving
  interaction-bias facet-pair selection explicit; the misfit preset bundles
  unexpected responses, displacement review, and pathway-map fit annotations.
- `mfrmr_output_guide("binary")` and the README now make the ordinary
  person-item dichotomous route explicit: pass the item column as the single
  non-person facet and use `model = "RSM"`; with exactly two ordered
  categories this is the usual binary Rasch logit up to centering and
  threshold-identification conventions.
- New `launch_mfrmr_viewer()` provides an optional local Shiny reader for
  existing `mfrm_results` objects. It does not estimate models, read external
  web applications, or change diagnostics; it displays the already-built
  overview, triage, status, QC evidence, APA-style report text when available,
  bias screens, pathway/misfit review, tables, plots, and replay code. The QC,
  Report, Bias, and Pathway/Misfit tabs now include local section-status
  tables so unavailable or not-requested sections are visible at the point of
  inspection.
- New `export_mfrm_results()` writes a lightweight download folder from an
  existing `mfrm_results` object: summary CSVs, collected tables, HTML, RDS,
  replay code, and a written-files manifest, with optional PNG plot export and
  best-effort zip creation. This is the compact result-object handoff route;
  `export_mfrm_bundle()` remains the broader fit-centered analysis archive.
- `export_mfrm_results()` now accepts `include = "report"` to write
  `mfrm_report()` artifacts into the same download folder: report-table CSVs
  such as `report_index`, evidence summaries, and reporting templates, plus
  report Markdown and report HTML.
- New `mfrm_report()` turns an existing `mfrm_results` object into
  report-ready QC, APA, validation, reviewer, or technical section plans with
  claim-readiness, report-gap, evidence-boundary, and next-action tables. It
  is a reporting layer over existing results, not a new estimator, diagnostic,
  or acceptance rule.
- `mfrm_report()` now includes `first_screen`, a FACETS-like entry table with
  an overall row and one row per major evidence area. It reports status,
  readiness, the main issue, next action, and primary route before users open
  the detailed evidence and template tables.
- `print(mfrm_report(...))`, the help-page examples, and the README now follow
  the same short reading order: read `summary(report)` and `report$first_screen`
  first, use `report$report_index` and `report$template_index` as table
  indexes, and open detailed evidence tables only when those indexes point to
  them.
- `summary(mfrm_report(...))` now provides a shorter reader-facing report
  summary: overview status, first-screen rows, immediate actions, optional
  not-requested sections, claim-readiness counts, report gaps, boundary rows,
  and standard routes. This keeps the comprehensive report object easier to
  read without creating a new diagnostic decision rule.
- `mfrm_report(output = "html")` now starts from the same reader-facing
  guidance and report-summary tables before the full Markdown report. This
  keeps browser output aligned with the `summary(report)` and `first_screen`
  workflow.
- `mfrm_report()` now includes a compact `report_index` table that lists
  the major evidence areas, evidence status, readiness label, review-signal
  count, and the primary/template tables to inspect next. This keeps the
  expanded report surface navigable without hiding the detailed tables.
- `mfrm_report()` `report_index` now also includes explicit evidence,
  template, plot, export, and `mfrm_results(include = ...)` routes. These
  columns keep report drafting, figure review, and download handoff connected
  to the same evidence surface without adding a new diagnostic decision rule.
- `mfrm_report()` now also exposes fit-specific `fit_criteria`,
  `zstd_conventions`, and `fit_decision_policy` tables. These make the
  selected MnSq threshold band, alternative published fit bands, and
  engine-vs-FACETS-style df/ZSTD conventions visible before users write fit,
  separation, or reliability claims.
- `mfrm_report()` now adds result-specific fit evidence tables from the stored
  `mfrm_results()` fit-measures component: observed fit-status counts,
  threshold-profile sensitivity, df/ZSTD sensitivity counts, and row-level
  df-sensitive prompts. This keeps FACETS-style ZSTD differences visible
  without turning them into a new acceptance rule or a different MnSq signal.
- `mfrm_report()` now includes `fit_reporting_templates`, a cautious wording
  scaffold for APA, QC, validation, reviewer, and technical reports. The
  templates summarize observed fit counts, threshold-profile sensitivity,
  ZSTD convention, and df sensitivity in separate sentences so fit,
  separation, and reliability are not collapsed into one pass/fail claim.
- All `mfrm_report()` reporting-template tables now share evidence and claim
  boundary columns: `EvidenceTable`, `EvidenceRoute`, `BoundaryType`,
  `ClaimStrength`, and `RecommendedUse`. This makes template wording easier to
  trace back to its source table and helps keep descriptive, caveated, and
  follow-up-only claims separate.
- `mfrm_report()` now includes `template_index`, a stacked template index
  over all fit, precision, bias, misfit/pathway, and linking/anchor reporting
  templates. It lets users review boundary type, claim strength, recommended
  use, and evidence route before opening full template text.
- `mfrm_report()` now also includes precision-specific reporting surfaces:
  `precision_evidence_summary`, `precision_basis`, and
  `precision_reporting_templates`. These summarize separation, reliability,
  strata, precision tier, and review/warn checks while keeping
  Rasch/FACETS-style separation reliability distinct from inter-rater
  agreement, model fit, and standalone validity evidence.
- `mfrm_report()` now includes bias-specific reporting surfaces when the
  source result was built with `include = "bias"`: `bias_evidence_summary`
  and `bias_reporting_templates`. These keep facet-level bias screens,
  interaction-bias contrasts, DFF follow-up, and fairness conclusions
  separated so screen-positive rows are not reported as final fairness or
  invariance decisions.
- `mfrm_report()` now includes misfit/pathway reporting surfaces when the
  source result was built with `include = "misfit_review"`:
  `misfit_evidence_summary` and `misfit_reporting_templates`. These keep
  unexpected-response rows, displacement review, pathway-map evidence, and
  case-review wording separate so local misfit prompts are not reported as
  automatic exclusion, fairness, or validity decisions.
- `mfrm_report()` now includes linking/anchor reporting surfaces when the
  source result was built with `include = "linking"`:
  `linking_evidence_summary` and `linking_reporting_templates`. These keep
  anchor readiness, drift review, screened equating-chain review, and GPCM
  support boundaries separate so anchor evidence is not reported as automatic
  drift absence, completed equating, DFF support, or validity proof.
- `mfrm_results(include = "linking")` now adds the fitted object's stored
  anchor-review evidence and an operational linking-readiness surface to the
  comprehensive results object. It exposes `plot(res, type = "anchors")` for
  anchor-readiness visualization and routes drift/equating follow-up to
  explicit multi-fit calls such as `detect_anchor_drift()` and
  `build_equating_chain()`.
- Recovery review can now retain fit/separation operating characteristics as
  diagnostic context. These summaries help users inspect MnSq, ZSTD,
  separation, and reliability behavior without making them top-line recovery
  gates. Recovery assessment and validation summaries also expose
  `diagnostic_reporting_notes` so zero separation/reliability and ZSTD
  sensitivity are routed into report caveats rather than recovery or release
  decisions.
- `assess_mfrm_recovery()` now includes `condition_review` and
  `condition_reporting_notes` for recovery simulations. For bounded `GPCM`,
  these tables separate slope-regime context and generated score-category
  support from recovery metrics before users interpret stress cases.
- Bounded-`GPCM` simulation specifications now carry `slope_regime` metadata.
  These labels are documented as package recovery-review labels, not
  literature-derived fit or adequacy cut points.
- The optional recovery-validation summary now separates the core release
  recovery decision from extended sensitivity evidence. This keeps stress
  cases visible without treating them as top-line release failures by default.
- `build_summary_table_bundle()` can now convert recovery-validation summaries
  into appendix-ready tables, including top-line decision, case, condition, and
  diagnostic reporting tables. `export_summary_appendix()` accepts the same
  recovery-validation summaries for CSV/HTML appendix handoff, and
  `export_mfrm_bundle(summary_tables = ...)` can co-locate those tables with a
  fit-based export bundle.
- `export_summary_appendix()` and `export_mfrm_bundle(summary_tables = ...)`
  now accept person-fit summary objects consistently with
  `build_summary_table_bundle()`.
- `precision_review_report()` can now be sent through
  `build_summary_table_bundle()` and appendix/export helpers. Its
  `fit_separation_basis` table remains a precision-review surface so fit,
  ZSTD, separation/reliability/strata, and QC thresholds are not mistaken for
  release or recovery success gates.
- `fit_measures_table()` and `facets_fit_review()` summaries can now be sent
  through the same appendix/export helpers. Their df/ZSTD sensitivity and
  optional external FACETS matching tables stay separate from MnSq fit status
  and top-line validation decisions.
- `reporting_checklist()` now includes a Global Fit row for the
  fit/separation reporting boundary, pointing users to
  `precision_review_report()`, `fit_measures_table()`, and
  `facets_fit_review()` before drafting fit, ZSTD, separation, or reliability
  claims.
- New observed-data resampling helpers, `build_mfrm_resampling_spec()` and
  `draw_mfrm_resamples()`, create person-clustered stratified subsample or
  bootstrap inputs with manifest tables for stratum representation and
  rater/facet coverage. These helpers are explicitly framed as stability or
  reproducibility inputs against a full-data reference, not true-parameter
  recovery evidence.
- `build_mfrm_sim_spec()` and `simulate_mfrm_data()` now support
  `assignment = "sparse_linked"` for planned-missing sparse rating designs.
  The generated data retain sparse-design metadata for design density,
  planned missingness, rater coverage, and rater-pair common-person links.
- `evaluate_mfrm_design()` and `evaluate_mfrm_recovery()` now carry sparse
  linked generators directly, including run-level and summary columns for
  planned missingness and rater-pair linking diagnostics.
- `build_summary_table_bundle()` now separates those sparse linked diagnostics
  into appendix-ready `sparse_design` tables for design-evaluation and
  recovery-simulation outputs, so planned missingness and rater linkage are
  not buried inside performance or recovery metrics. The table also labels
  zero common-person rater pairs and requested-link target shortfalls as
  design-review issues, not recovery failures.
- `summary(evaluate_mfrm_design(...))`, `summary(evaluate_mfrm_recovery(...))`,
  and `build_summary_table_bundle()` now include a compact `sparse_review`
  table when sparse linked designs are active. `plot.mfrm_design_evaluation()`
  also accepts sparse-design metrics such as `plannedmissingrate`,
  `mincommonpersons`, `zerocommonpairs`, and `pairsshorttarget`.
- New `build_mfrm_network_review()` synthesizes the existing design-network
  analysis into a reportable review surface: connectedness, articulation
  points, bridge edges, facet-level vulnerability, optional sparse-linking
  diagnostics, and a reporting map that keeps network evidence separate from
  MFRM estimates, fit, separation, and recovery gates.
- New `build_peer_review_sim_spec()` creates peer-review / peer-assessment
  simulation specifications where submissions and reviewers share the same ID
  universe, self-review can be structurally excluded, and common-link anchor
  submissions can be assigned many reviewers. Generated data carry
  peer-review design metadata for assignment density, reviewer load,
  reciprocal review pairs, and common submissions per reviewer pair.
  `build_peer_review_design_review()` turns the same metadata into
  appendix-ready assignment diagnostics, and `build_mfrm_network_review()` can
  include the metadata alongside graph connectedness review.
- README, help pages, and vignettes now show the recommended reading order for
  recovery review: `summary(recovery_review)`, condition notes/review,
  diagnostic notes/review, status plot, metric plot, then row-level recovery
  rows that need follow-up. Recovery assessment and validation summaries now
  expose this as a `reading_order` table so users can find the next table to
  inspect without opening plot data.

# mfrmr 0.2.0

Documentation accuracy pass plus research-grounded visualization and
GPCM bias-inference refinements. Documentation, citations, and band
attributions are corrected against primary sources, with mathematical
screening-SE corrections, Snijders-corrected person-fit reporting where
the assumptions are met, and clearer plot data for review.

This release keeps the 0.1.6 defaults, but it is not only an
infrastructure polish release. Public review helpers have been
consolidated on the `*_review*` names documented below, and the former
`*_audit*` public spellings, S3 compatibility classes, and duplicate
top-level fields have been removed as a deliberate breaking cleanup.

## Release overview

For most users, the main changes in 0.2.0 are:

- **More defensible mathematics and inference**: RSM/PCM/GPCM reductions,
  GPCM slope handling, fit df/ZSTD conventions, Snijders-corrected person fit,
  information curves, recovery simulations, and score-support edge cases are
  now covered by explicit regression tests.
- **Clearer FACETS relationship**: `mfrmr` is positioned as a package-native
  MFRM workflow with FACETS-style tables, review helpers, and migration routes,
  not as a promise to numerically reproduce every FACETS estimate.
- **Better user-facing diagnostics**: fit-measure tables, data-quality reports,
  category curves, residual-PCA follow-up, person-fit summaries, recovery
  checks, and reporting bundles now expose structured tables before asking
  users to interpret plots or console output.
- **R-first visualization access**: plot helpers increasingly return reusable
  `draw = FALSE` data, long-form plot tables, annotations, and style metadata
  so users can rebuild figures in ggplot2, plotly, Quarto, or other reporting
  workflows.
- **Quieter, more reproducible workflows**: routine preparation and
  rating-range messages are stored in fit/data objects and are opt-in at the
  console; long-running design evaluation shows progress only in interactive
  runs by default.

Breaking changes in 0.2.0 are intentional and concentrated around public naming
clarity: former exported `*_audit*` helper names and their compatibility S3
classes were removed in favour of the canonical `*_review*` surface. Model
defaults from 0.1.6 are retained.

The detailed notes below are organized as follows:

- user pathways, output contracts, visualization, and naming changes;
- mathematical/statistical corrections and regression-test coverage;
- recovery simulation and validation workflow;
- citation/documentation corrections;
- smaller feature additions, bug fixes, build hygiene, and deferred work.

## User pathways, output contracts, and terminology

- **FACETS positioning guide**: new `facets_positioning_guide()` makes the
  package boundary explicit for reports and migration notes. `mfrmr` is not
  presented as a FACETS numerical clone: estimates remain package-native
  unless external FACETS output is supplied for comparison, while FACETS-style
  wrappers, coverage tables, and output files serve transition, handoff, and
  report-organization purposes.
- **Report-ready FACETS relationship wording**: `reporting_checklist()` and
  `summary(reporting_checklist(...))` now carry a `facets_positioning` table
  so Quarto, appendix, and handoff workflows can quote the same boundary
  language used by `facets_positioning_guide()`. `build_summary_table_bundle()`
  includes this table as `facets_relationship_wording`.
- **FACETS output-contract review naming**:
  `facets_output_contract_review()` is now the sole public helper for
  FACETS-style output-contract review. The returned bundle class is
  `mfrm_facets_contract_review`; the helper checks package output columns and
  derived metric consistency against the FACETS-style contract, and does not
  claim numerical FACETS equivalence. Public result components use
  `column_review` and `metric_checks`, not older bookkeeping labels.
- **FACETS pathway for bias and Wright maps**: `mfrmr_output_guide("facets")`
  now explicitly routes FACETS users to Table 14-style bias/interactions via
  `estimate_bias()`, `bias_interaction_report()`, `bias_pairwise_report()`,
  and `plot_bias_interaction()`, and to variable-map review via
  `plot(fit, type = "wright")`, `plot_wright_unified()`, and
  `plot_data(type = "wright")`.
- **FACETS pathway for anchors and category outputs**:
  `mfrmr_output_guide("facets")` now also exposes direct/group anchor routes
  through `review_mfrm_anchors()`, `make_anchor_table()`, and
  `fit_mfrm(anchors = ..., group_anchors = ...)`, drift/linking review through
  `detect_anchor_drift()` and `plot_anchor_drift()`, and FACETS-style category
  and fair-average routes through `rating_scale_table()`,
  `category_structure_report()`, `category_curves_report()`, and
  `fair_average_table()`.
- **Standalone residual and subset file writers**: new
  `write_mfrm_residual_file()` writes observation-level observed, expected,
  residual, standardized residual, score-information, and optional category
  probability columns to CSV/TSV. New `write_mfrm_subset_file()` writes
  connected-subset summaries plus node-membership files for external linking
  review, without forcing users through the legacy graph/score output bundle.
- **Category-specific information curves**: `category_curves_report()` now
  includes a `category_information` table and
  `plot(..., type = "category_information")`. Category contributions use
  `a^2 P_k(theta) (k - E[X | theta])^2` and sum to the total curve
  information at each theta value; for PCM/RSM this reduces to the unit-slope
  form.
- **Cumulative probability curves**: `category_curves_report()` now also
  includes `cumulative_probabilities` and `cumulative_boundaries`, with
  `plot(..., type = "cumulative")` for the FACETS/Winsteps-style
  accumulated category-probability view. Both `P(X <= k)` and flipped
  `P(X >= k)` directions are returned; boundary rows report approximate
  theta values where `P(X <= k) = .5`, with crossing status columns so
  out-of-range or multiple-crossing boundaries are not over-interpreted.
- **Category-curve overview plot**: `plot(category_curves_report(...))`
  now defaults to an overview panel that shows category probabilities,
  cumulative probabilities, total information, and category-specific
  information together. Existing focused views remain available through
  `type = "ogive"`, `"ccc"`, `"cumulative"`, `"information"`, and
  `"category_information"`. The plot now also supports
  `preset = "monochrome"` with line-type separation and explicit cumulative
  `.5` boundary-line control through `boundary_status = "in_range"`,
  `"all"`, or `"none"`. `plot_data(..., component = "plot_long")` returns a
  ggplot2/plotly-friendly long table spanning ogive, category-probability,
  cumulative-probability, total-information, and category-specific-information
  series, with `curve_style` carrying the resolved color/line-type mapping.
- **Quieter rating-range provenance**: `fit_mfrm()` no longer emits an
  informational message for routine observed-score range inference by
  default. The same provenance is retained in `fit$prep`,
  `summary(fit)$settings_overview`, and `describe_mfrm_data()` so users can
  still tell whether `rating_min` / `rating_max` were declared or inferred.
  Set `options(mfrmr.show_inferred_rating_range = TRUE)` to restore the
  one-time message during interactive checks.
- **Structured data-preparation notes**: `prepare_mfrm_data()` now stores
  row retention and preparation notes in `fit$prep$row_retention` and
  `fit$prep$preparation_notes`. Row drops, whitespace trimming, duplicate
  person-by-facet cells, and single-level facets are therefore available to
  `summary(fit)` and `summary(describe_mfrm_data(...))` instead of existing
  only as transient console messages. Routine row-drop, trim, and
  single-level-facet messages are quiet by default; set
  `options(mfrmr.show_preparation_messages = TRUE)` to show them during
  interactive checks.
- **Fit-annotated pathway plot data**: `plot(fit, type = "pathway",
  draw = FALSE)` now returns R-friendly `pathway_long` and
  `pathway_annotations` tables alongside `fit_measures`, `fit_status`,
  `curve_fit_status`, and `fit_measure_status`. This lets users rebuild
  FACETS-style pathway maps in ggplot2, plotly, or Quarto while retaining the
  same underfit/overfit labels used by `fit_measures_table()`.
- **R-first plot-data contracts for bias and information plots**:
  `plot_bias_interaction(..., draw = FALSE)` now exposes `plot_long`,
  `plot_annotations`, `flag_summary`, and `plot_settings` across scatter,
  ranked, heatmap, and facet-profile views. `compute_information()` now stores
  `conditional_sem`, `information_long`, and a precision/SEM summary, and
  `plot_information(type = "sem")` / `"csem"` are supported aliases for the
  conditional standard-error-of-measurement curve.
- **Category-probability plot aliases and annotations**:
  `plot(category_curves_report(...), type = "category_probability")` and
  `type = "conditional_probability"` now route to the same category-probability
  curves as `type = "ccc"`, matching FACETS/Winsteps terminology without
  changing the underlying probability data. Draw-free plot data now include
  `plot_annotations` and `curve_summary` alongside `plot_long`,
  `curve_style`, `boundary_lines`, and `plot_settings`.
- **FACETS feature coverage matrix**: new `facets_feature_coverage()` gives a
  public, release-scoped map from the FACETS 64-bit output index to current
  `mfrmr` routes, separating `implemented`, `partial`, `not_implemented`, and
  `not_targeted` surfaces. This makes unsupported FACETS-specific outputs such
  as Winsteps control-file export, raw FACETS report parsing, and arbitrary
  Web/Excel menu plots explicit rather than implicit.
- **G-study / D-study planning route**: new `mfrm_d_study()` extends
  `mfrm_generalizability()` from observed variance-component review to planned
  design comparison. It reports projected `G` and `Phi` under alternative
  numbers of raters, criteria, or other random measurement facets, and exposes
  residual-scaling sensitivity assumptions so simplified G-study residuals are
  not silently over-interpreted. `plot(mfrm_d_study(...), draw = FALSE)` and
  `plot_data()` expose reusable coefficient/error-variance series for custom
  design-planning visuals. `plot.mfrm_d_study()` now supports line plots with
  `group_var`, ggplot2-like `panel_by` / `panel_grid` small multiples, and
  two-axis `heatmap` / `contour` views for rater-by-task design grids. An
  optional `surface3d` view is available for exploration, while heatmap/contour
  remain the recommended reporting displays. Plot data labels these
  coefficients as `MetricFamily = "G-theory"` so they are not confused with
  IRT or classical-test-theory reliability coefficients.
- **Connectivity network visualization**: `subset_connectivity_report()` now
  includes reusable node/edge tables, and `plot(..., type = "network")`
  provides an igraph-based co-observation graph when `igraph` is installed.
  With `draw = FALSE`, the returned plot data supports custom R visualization
  without depending on the base plotting default.
- **MFRM design-network analysis**: new `mfrm_network_analysis()` treats the
  person/facet-level observation design as an undirected weighted
  co-observation graph and returns graph-level connectedness, node degree and
  strength, betweenness/closeness, articulation points, bridge edges, and
  facet-level vulnerability summaries. These are explicitly framed as design
  linking diagnostics, not as person ability, rater quality, or model-fit
  statistics. `plot(..., type = "centrality")`, `plot(..., type =
  "facet_summary")`, and `plot(..., type = "network")` provide immediate
  visual checks with draw-free plot data.
- **Rater-effect network analysis**: new `rater_network_analysis()` adds a
  Lamprianou-style pairwise rater network route separate from design
  connectedness. It supports agreement, disagreement, and directed
  severity-direction networks, returns rater-level in/out strength,
  betweenness/closeness, a finite network severity index, retained edge
  tables, and all pairwise metrics used before thresholding. The help page
  states that these indices are descriptive network diagnostics rather than
  Rasch logit estimates or formal fit statistics. `plot(..., type =
  "network")`, `"severity"`, `"centrality"`, and `"matrix"` provide immediate
  visual checks and reusable plot data.
- **Halo-effect network screening**: new `rater_halo_network_analysis()`
  reshapes observed ratings into rater-by-criterion nodes, computes
  Spearman/Pearson/Kendall node-pair correlations, labels same-rater
  cross-criterion edges as `halo`, and contrasts halo-edge weights with
  non-halo edges. The default Bonferroni-adjusted edge filter follows the
  conservative network-screening strategy used in Lamprianou's halo example.
  The returned bundle includes `summary`, `node_metrics`, `edge_metrics`,
  `pair_metrics`, `halo_summary_by_rater`, and caveats that the Welch
  halo/non-halo comparison is descriptive because network edges are dependent.
  `halo_summary_by_rater` now includes `ReviewStatus` and `ReviewReason`
  based on same-rater cross-criterion mean weight, incident non-halo
  comparison weight, and retained halo-edge count, with labels framed as
  screening priorities rather than causal halo diagnoses.
  `plot(..., type = "edge_distribution")`, `"halo_summary"`, `"matrix"`, and
  `"network"` provide immediate visual review and draw-free plot data.
- **Fit-measures review table**: new `fit_measures_table()` gives a direct
  FACETS-style fit-measure view for raters, criteria, or other facet elements.
  It returns both R-friendly columns and a `facets_table` with labels such as
  `Infit MnSq`, `Outfit ZStd`, `Fit Status`, and `Review Reason`; `underfit`,
  `overfit`, and `mixed` subsets are included for immediate review.
- **Fit-threshold sensitivity summaries**: `fit_measures_table()` now returns
  `profile_summary_by_facet` and `profile_summary_overall`, reporting
  underfit, overfit, mixed, and any-flag rates under multiple literature-based
  MnSq bands from Linacre, Bond & Fox, and Wright & Linacre. The main table
  still uses the active review band, while `threshold_profiles` controls
  whether literature, active, all, or no profile summaries are returned.
- **Fit-measure CI display**: `fit_measures_table(ci_level = ...)` now adds
  approximate measure confidence intervals to both the R-friendly table and
  `facets_table`. `plot(fit_measures, type = "measure_ci",
  ci_level = ..., preset = "monochrome")` draws an interval plot with the
  requested confidence level.
- **FACETS df/ZSTD guide**: new `facets_fit_df_guide()` documents the
  engine-vs-FACETS-style degrees-of-freedom distinction and the MnSq-to-ZSTD
  transformation workflow. `fit_measures_table(fit_df_method = "both")` now
  exposes primary df plus FACETS-style companion df/ZSTD columns, and
  adds `df_sensitivity`, `df_sensitive`, and `df_sensitivity_summary` so users
  can identify rows whose |ZSTD| flag status or interpretation is
  convention-sensitive. The df-sensitivity screen exposes explicit
  `df_zstd_tolerance`, `df_zstd_large_shift`, and `df_ratio_tolerance`
  settings so FACETS-style reviews can be reproduced under stricter or more
  permissive rules. `plot(fit_measures, type = "df_sensitivity")` visualizes
  the largest engine-vs-FACETS-style ZSTD shifts. `facets_fit_review()` now
  uses the same row-level df-sensitivity engine and returns
  `df_sensitivity`, `df_sensitive`, and `df_sensitivity_summary`; the former
  `internal_comparison` component has been removed to keep the public output
  vocabulary consistent. External FACETS ZSTD tolerance is now named
  `external_zstd_tolerance`, separating external-table comparison from
  engine-vs-FACETS-style df sensitivity. The same
  `plot(..., type = "df_sensitivity")` route is available for
  `facets_fit_review()` bundles.
- **Data quality score-support review**: `data_quality_report()` now keeps the
  full fitted score support in `category_counts`, adds
  `score_support_review`, and surfaces `caveats` for zero-frequency categories.
  Intermediate gaps such as a declared 1-5 scale with observed `1, 2, 4, 5`
  are flagged either as retained zero-count categories (`keep_original = TRUE`)
  or as original-label gaps hidden by internal recoding (`keep_original =
  FALSE`). `plot(data_quality_report(...), type = "score_support")` highlights
  those categories, with `preset = "monochrome"` available.
- **Facet-level category-usage QC**: `data_quality_report()` now adds
  `category_usage_by_facet` and `category_usage_summary`, covering every fitted
  facet level crossed with the retained score support. This flags local zero
  and sparse category use, such as a rater who never uses the middle category
  even when the category appears elsewhere in the data.
  `plot(data_quality_report(...), type = "facet_category_usage")` provides a
  quick view of affected facet levels.
- **Data quality dashboard**: `plot(data_quality_report(...), type =
  "dashboard")` now combines row review, score-support category use,
  facet-level category-use issues, and missing/invalid row counts in one
  base-R view. With `draw = FALSE`, the returned plot data contains all four
  panel tables for report handoff.
- **Data quality flags**: `data_quality_report()` now includes
  `quality_flags`, a prioritized QC table that summarizes row exclusions,
  unknown design levels, score-support gaps, facet-level category-use cautions,
  and restricted facet response patterns with counts and next actions.
  `summary(data_quality_report(...))` previews this table when any priority
  flag is present.
- **Facet response-pattern QC**: `data_quality_report()` now adds
  `facet_response_patterns`, which flags facet levels that use only one score
  category, assign one category to at least the configured dominant-category
  cutoff, or use only boundary categories. This catches cases such as a rater
  assigning score 1 to all responses on a 1-5 scale.
- **Data quality overview and score map plots**: `data_quality_report()` now
  adds `quality_overview`, a compact area-level status table for rows,
  score support, facet-level category use, facet response patterns, and design
  matching. `plot(..., type = "quality_flags")` summarizes priority QC flags by
  area, `plot(..., type = "facet_response_patterns")` shows dominant local
  category use by facet level, and
  `plot(..., type = "score_map")` shows original-to-internal score mappings
  when labels have been recoded.
- **User-pathway and plot-data access**: `mfrmr_output_guide("facets")`,
  `mfrmr_output_guide("conquest")`, and `mfrmr_output_guide("r")` now give
  focused starting points for users arriving from FACETS, ConQuest, or
  R-first visualization workflows. New `plot_data()` extracts the full
  reusable plot-data list, or one named component, from any `mfrm_plot_data` object or
  mfrmr plot helper that supports `draw = FALSE`. New
  `plot_data_components()` lists each reusable plot-data component, its shape,
  role, accessor call, and custom-graphics notes so users can discover
  `plot_long`, annotation, settings, style, and review tables without reading
  the underlying list structure.
- **Monochrome plot preset**: plot helpers that use the package visual preset
  system now accept `preset = "monochrome"`. Color remains the default
  (`"standard"`), while monochrome supports print-oriented figures and
  color-independent review.
- **Interval-aware visualization guide**: new `mfrmr_interval_guide()` maps
  public 95% CI / uncertainty routes across fit-measure tables, Wright maps,
  fair averages, bias screens, displacement, DFF/DIF summaries, anchor drift,
  rater severity profiles, rater trajectories, manuscript Figure 1 composites,
  shrinkage, and ICC review. The guide records the interval basis and
  interpretation boundary so CI displays remain precision or screening
  evidence rather than automatic fit, fairness, or validity decisions.
- **Searchable help concepts**: high-level route guides and CI-capable plot /
  table helpers now carry Rd concept tags such as `confidence intervals`,
  `visual diagnostics`, `reporting workflow`, `route selection`, and
  `GPCM boundaries`. This makes `help.search()` useful for finding the right
  help page before users know the function name.
- **Shrinkage CI plot data**: `plot_shrinkage_funnel(show_ci = TRUE,
  ci_level = ...)` now draws approximate raw and shrunken estimate whiskers and
  returns `RawCI_Lower`, `RawCI_Upper`, `ShrunkCI_Lower`, `ShrunkCI_Upper`,
  and `CI_Level` for downstream graphics. `plot(fit, type = "shrinkage")` now
  uses the requested `ci_level` for CI whiskers instead of a fixed 95% level.
- **Shrinkage figure guidance**: `visual_reporting_template()` and the visual
  diagnostics vignette now include an empirical-Bayes shrinkage funnel route,
  caption skeleton, beginner check, and interpretation guardrails so users can
  report shrinkage movement without treating it as automatic rater-quality,
  bias, or validity evidence.
- **Response-time diagnostic layer**: new `response_time_review()` and
  `plot_response_time_review()` summarize response-time metadata outside the
  fitted MFRM likelihood. The review returns rapid/slow thresholds,
  event-level flags, person/facet/score summaries, and grouped plot data so
  timing patterns can be reviewed as descriptive QC rather than joint
  speed-accuracy parameters or automatic exclusion rules.
  `mfrmr_output_guide("response_time")` and the R-user pathway now expose this
  route alongside other review and reusable plot-data helpers.
  `mfrm_results(include = "response_time", response_time = ...,
  response_time_data = ...)` can now carry the same descriptive review into
  the first-screen result object, `summary(res)$next_actions`,
  `plot(res, type = "response_time")`, the local viewer payload, and
  `export_mfrm_results()` table exports without changing fitted MFRM
  estimates.
- **Bounded GPCM boundary text**: README, vignettes, help pages, and
  unsupported-path messages now state the current `GPCM` scope consistently:
  direct data generation, parameter recovery, fair averages, bias screening,
  summary-table bundles, appendix export, caveated APA/QC/export bundles, and
  exploratory linking review are available within documented caveats. Role-based
  design evaluation and population forecasting are also available as caveated
  bounded-`GPCM` sensitivity evidence when the requested design preserves the
  simulation specification's slope structure. Role-based diagnostic and
  signal-detection design screening is available as caveated slope-aware
  operating-characteristic evidence with `gpcm_boundary` output. Full
  FACETS-style score-side contract review or score-side equivalence, posterior
  predictive checks, and heavy backend routes remain outside the validated
  route. The validation artifacts now include a bounded-`GPCM` roadmap so
  those `supported_with_caveat`, `blocked`, and `deferred` rows are tracked as
  explicit release-scope decisions instead of implicit support gaps.
- **GPCM score-side boundary checks**: blocked bounded-`GPCM` score-side
  helpers now stop before producing partial outputs. The
  `facets_output_file_bundle()` and
  `facets_output_contract_review()` help pages also state that graph and
  package-native scorefile output are caveated bounded-`GPCM` routes while
  full FACETS output-contract review remains outside the bounded-`GPCM`
  boundary. `gpcm_score_side_contract()` records
  the estimand, uncertainty, reduction-test, schema, guard, and release-wording
  requirements that separate caveated scorefile output from full FACETS-style
  score-side review.
- **GPCM route guidance**: `gpcm_capability_matrix()` now includes
  `RecommendedRoute` and `NextValidationStep` columns, so each supported,
  caveated, blocked, or deferred helper family states both the current
  substitute workflow and the validation evidence needed before the boundary
  can move.
- **GPCM out-of-scope route guidance**: blocked and deferred bounded-`GPCM`
  routes now report the matching capability-matrix row, recommended substitute
  route, and next validation step before returning an error. The errors carry
  class `mfrmr_gpcm_scope_error` with helper, area, status, recommended-route,
  and next-validation-step fields for programmatic handling.
- **GPCM route-boundary coverage**: the capability-matrix tests and validation
  review now verify that every blocked or deferred bounded-`GPCM` row is
  represented by a structured stop condition or an explicit future-scope entry,
  and that those messages stay synchronized with the matrix's area, status,
  route, and next-validation-step fields. `gpcm_runtime_guard_coverage()` is
  exported as the public table for this route-boundary check.
  `mfrmr_output_guide("gpcm")` now points users to both
  `gpcm_capability_matrix()` and this coverage table.
- **GPCM release-readiness alignment**: `inst/validation/release-readiness.R`
  now checks that blocked and deferred `gpcm_capability_matrix()` rows have
  non-empty route guidance, are represented in the installed bounded-`GPCM`
  scope notes, and are covered by future-scope rows in the release-evidence
  checklist.
- **User-facing wording**: public visualization and reporting documentation now
  uses "plot data", "surface data", or "data handoff" instead of implementation
  terminology where possible, while retaining actual field names such as
  `plot_payloads` only where users need to access them.
- **RSM/PCM wording review**: package-level, fitting, information, bias, and
  export docs now distinguish the equal-weighting `RSM` / `PCM` reference route
  from the broader ordered-response model surface. Stale "Rasch-only" and
  "legacy-compatible" labels were narrowed where they described helpers that
  now also serve bounded `GPCM` or direct summary-table workflows.
- **Model-choice user guide**: README and the GPCM/MML vignettes now give
  user-facing guidance for choosing `RSM`, `PCM`, or bounded `GPCM`, including
  report wording templates and a warning that better `GPCM` fit is sensitivity
  evidence rather than an automatic operational-scoring decision. A
  documentation terminology regression test now guards against reintroducing the
  stale Rasch-only phrasing removed in this pass.
- **Fit-level model-choice review**: new `build_model_choice_review()` bundles
  `compare_mfrm()`, model-role guidance, downstream route availability, report
  wording templates, the bounded-`GPCM` support matrix, and an optional
  `build_weighting_review()` run so users can review `RSM` / `PCM` /
  bounded-`GPCM` candidates from the actual fitted objects.
- **Review-name migration completed as a breaking cleanup**:
  `review_mfrm_anchors()`, `precision_review_report()`,
  `review_conquest_overlap()`, `facet_small_sample_review()`, and
  `build_weighting_review()` are now the only exported review-name
  implementations. The former public `*_audit*` function spellings, their S3
  compatibility classes, and old public component names such as
  `orientation_audit`, `nesting_audit`, `hierarchical_audit`, and
  `shrinkage_audit` have been removed or renamed rather than shown as
  user-facing migration artifacts.
- **Compatibility registry narrowed**:
  `compatibility_alias_table()` now lists only retained compatibility names that
  remain part of the public package surface, such as `mfrmRFacets`,
  `analyze_dif`, `JMLE`, and long-standing output column aliases. It no longer
  advertises removed review-name migration artifacts.
- **Public output review wording**: linking-review tables and plot routes now
  use `anchor_review` labels in user-facing source metadata, model-choice
  raw objects and summaries expose only `weighting_review_status` /
  `weighting_review`, and diagnostics summary-table bundles export
  `precision_review` without the old duplicate table. Bias reports now expose
  `orientation_review`; model-comparison output uses `nesting_review`; and
  reproducibility manifests use `hierarchical_review` and `shrinkage_review`.
  DFF subgroup refit rows now record anchor-review notes in `LinkingReview`.
- **Review component accessors**: new `anchor_review()` and
  `precision_review()` helpers provide a stable route to package-native
  review components. They intentionally read only canonical `*_review` fields.
- **Reference review naming**: `reference_case_review()` is now the canonical
  package-native report-completeness helper, and reporting-checklist /
  cheatsheet wording now uses review labels for hierarchical and complete-case
  follow-up items.
- **Review-wording guardrail**: current public guides and generated help now
  avoid exposing `audit` as a user-facing package concept. Data-quality
  row-status output uses `row_review`; prediction provenance uses
  `row_review` / `population_review`; bias, model-comparison, and manifest
  components use `*_review` names; and ordinary user-facing guidance uses
  review/check/traceability terminology.
- **Output helper guide**: new `mfrmr_output_guide()` gives users a compact
  purpose-to-helper map for choosing among `*_table`, `*_report`,
  `*_review`, `*_bundle`, `export_*`, and compatibility routes. The
  compatibility guide now also states that old `*_audit` helper and component
  names are not part of the 0.2.0 public surface.

## Mathematical and inferential corrections

- **Identified step/threshold parameterization**: `RSM`, `PCM`, and
  bounded `GPCM` now optimize step/threshold profiles with the correct
  sum-to-zero degrees of freedom (`steps - 1` per profile). Earlier
  pre-release implementations centered the step vector after optimization
  but still left the centered-away null direction in the optimizer,
  AIC/BIC parameter count, and Hessian. The point-estimate scale is
  unchanged; the likelihood parameter count and observed-information
  basis are now aligned with the stated identification constraint.
- **MML joint covariance layer for structural parameters**:
  `diagnose_mfrm()` now reuses one observed-information covariance for
  non-person facet SEs and exposes the same covariance basis for
  step/threshold and bounded-`GPCM` slope uncertainty in
  `diagnostics$parameter_uncertainty`. Step rows get `SE`, normal
  `CI_Lower` / `CI_Upper`, and covariance status metadata. GPCM slope
  rows get log-slope SEs plus positive-scale delta-method SEs and
  log-normal confidence limits. `fit_mfrm(..., attach_diagnostics = TRUE)`
  attaches those structural SE columns to `fit$steps` and `fit$slopes`
  when the MML Hessian is available.
- **Measure-level CI contract**: `diagnose_mfrm()$measures` now records
  `CI_Level = 0.95` and `CI_Method = "Normal approximation"` alongside
  `CI_Lower` / `CI_Upper`. The interval calculation uses
  `qnorm(0.975)` rather than a rounded multiplier, while row-level
  `CIEligible`, `CIBasis`, and `CIUse` continue to distinguish primary
  reporting intervals from review or screening approximations.
- **Weighted BIC transparency**: `compare_mfrm()` now reports `WeightedN`,
  `ICSampleSize`, and `ICSampleSizeBasis` in its comparison table. This makes
  the BIC penalty basis explicit: ordinary fits use row count, while weighted
  fits use the sum of weights already used by the fitted model summary.
- **Residual-PCA boundary handling**: exploratory residual-PCA helpers now
  capture non-fatal PCA-engine warnings inside the returned PCA bundle instead
  of emitting them as loose warnings during diagnostics or plotting. Degenerate
  residual-correlation conditions therefore remain visible for review without
  looking like confirmatory test failures.
- **Residual-PCA parallel analysis**: `analyze_residual_pca()` now supports
  `parallel = TRUE` for residual-permutation parallel analysis. The null
  comparison permutes standardized residuals within residual columns, preserving
  column distributions and missingness while breaking residual association.
  PCA tables gain `ParallelMean`, `ParallelCutoff`,
  `ExcessOverParallelCutoff`, and `ExceedsParallelCutoff`, and
  `parallel_status` records availability and successful permutation counts.
  `plot_residual_pca()` adds `parallel_scree` and `parallel_excess` views.
  This is reported as exploratory follow-up evidence for dimensionality review,
  not as a standalone proof of unidimensionality or multidimensionality.
- **GPCM fair-average structural SEs**: `fair_average_table(fair_se = TRUE)`
  now adds opt-in structural delta-method SE and CI columns for bounded
  `GPCM` fair averages when the MML observed-information covariance is
  available. The original `SE` / `Model S.E.` / `Real S.E.` columns keep their
  measure-SE meaning; fair-average uncertainty is exposed in distinct columns
  such as `Fair(M) S.E.`, `AdjustedAverageSE`, and
  `AdjustedAverageCI_Lower` / `AdjustedAverageCI_Upper`. Person rows remain
  unavailable because MML person EAP estimates are conditioned on rather than
  included in the structural Hessian. `summary(fair_average_table(...))` and
  its print method now surface whether fair-average SEs were requested, how many
  rows are available, and the resulting status mix. `plot_fair_average(show_ci =
  TRUE)` uses these columns automatically for bounded-`GPCM` fit objects.
- **GPCM expected-score consistency**: the internal `expected_score_table()`
  route now uses the same response-probability bundle as diagnostics and
  category-count calculations, so bounded `GPCM` expected scores respect the
  fitted slope parameters instead of falling through to the PCM kernel.
  Fair-average documentation now also states that non-slope-facet rows use an
  identification-based reporting convention, not a FACETS score-side
  equivalence claim.
- **GPCM invalid-slope guard**: the low-level GPCM expected-score helper no
  longer treats non-finite, zero, or negative slopes as slope = 1. It returns
  unavailable expected scores instead, so a malformed bounded-`GPCM` object
  cannot silently become a PCM-style calculation in fair-average internals.
  The internal iteration-state replay helper also now has an explicit GPCM
  probability-kernel branch, avoiding a latent PCM fallback if that route is
  later moved inside the supported GPCM boundary.
- **RSM-to-PCM reduction checks**: the test suite now pins the `RSM` special
  case as common-threshold `PCM`. Under identical common thresholds, `RSM` and
  `PCM` must agree for category probabilities, unweighted and weighted
  log-likelihoods, response-bundle diagnostics, generated simulation data, and
  reconstructed simulation probabilities. The public `compare_mfrm(...,
  nested = TRUE)` path now also has a regression check that the reported LRT
  degrees of freedom equal the identified RSM-to-PCM step-structure difference.
- **Boundary-safe LRT reporting**: `compare_mfrm(..., nested = TRUE)` now
  records `comparison_basis$lrt_status` and `comparison_basis$lrt_reason`.
  Non-finite log-likelihoods, equal parameter counts, unsupported nesting, or
  negative likelihood-ratio statistics no longer fail silently or imply a
  model-choice conclusion; the LRT is withheld and the print/summary path states
  why it was not reported.
- **GPCM slope-scale consistency**: GPCM simulation specifications now treat
  supplied slopes as relative discriminations and normalize them to the same
  geometric-mean-one log-slope identification used by `fit_mfrm()`. Recovery
  summaries compare identified log slopes without an additional mean-alignment
  step, so absolute slope-scale bias is not hidden by the recovery table.
- **PCM reduction check for GPCM simulation**: the simulation tests now pin the
  special case in which bounded `GPCM` has unit slopes. With the same
  simulation specification, seed, and step-facet thresholds, the generated
  visible data and reconstructed category-probability matrix must match `PCM`
  to numerical tolerance. This guards the intended mathematical reduction
  without implying that a freely estimated `GPCM` fit should equal a `PCM` fit.
- **PCM reduction check for downstream diagnostics**: the unit-slope `GPCM`
  reduction is now also tested at the response-probability bundle layer. The
  bundle used by expected scores, variance-based fit diagnostics, and
  information calculations must match the `PCM` bundle for probabilities,
  expected category scores, score variances, fourth central moments, and score
  information. A companion sensitivity check verifies that non-unit slopes move
  the same quantities away from the `PCM` values, so the GPCM path is neither a
  hidden PCM fallback nor an unconstrained divergence.
- **Mathematical consistency regression tests**: the test suite now pins
  probability, expectation, variance, fourth-moment, information, and
  conditional-SEM identities across the low-level bounded-`GPCM` response
  bundle, `RSM` / `PCM` / bounded-`GPCM` category-curve reports, draw-free
  CCC/pathway plot data, and `compute_information()`. The checks also require
  facet-level information-contribution curves to aggregate back to the total
  information curve. This guards user-visible visualization and reporting
  tables against drifting away from the probability kernels.
- **Fit-measure consistency regression tests**: the test suite now verifies
  that `fit_measures_table()` preserves the documented df/ZSTD formulas,
  confidence-interval formulas, active fit-status labels, threshold-profile
  counts and rates, df-sensitivity status taxonomy, and draw-free
  `measure_ci` / `df_sensitivity` plot data. This pins the FACETS-style
  reporting surface to the same row-level calculations users see in the
  returned tables.
- **Data-quality consistency regression tests**: the test suite now verifies
  that `data_quality_report()` summary counts are recomputable from returned
  detail tables, that `quality_flags` and `quality_overview` summarize the
  same QC evidence, and that draw-free `quality_flags`,
  `facet_category_usage`, `facet_response_patterns`, and `score_map` plot data
  preserve score-support gaps, facet-level category-use issues, restricted
  rater response patterns, and original-label gaps hidden by score recoding.
- **GPCM bias SE (`estimate_bias()`)**: the conditional plug-in SE for
  the additive bias shift now uses the correct GPCM information
  \(\sum_i a_i^2 \mathrm{Var}(X_i)\). The previous pre-release
  implementation optimized the point estimate with the slope-aware GPCM kernel
  but used the PCM information \(\sum_i \mathrm{Var}(X_i)\) for
  `S.E.` / `t` / `Prob.`. The review label remains `"screening"`.
- **FACETS-style fit ZSTD df layer**: `diagnose_mfrm()` now accepts
  `fit_df_method = "engine"`, `"facets"`, or `"both"`. The default keeps the
  existing package-native df convention
  (`DF_Infit = sum(Var * Weight)`, `DF_Outfit = sum(Weight)`). The FACETS path
  adds the Wright-Masters/FACETS fourth-moment df approximation
  (`df = 2 / q^2`) and caps FACETS-style ZSTD values at +/-9. Use
  `fit_df_method = "both"` when comparing mfrmr fit flags with FACETS output:
  it preserves the engine ZSTD columns and adds `DF_Infit_FACETS`,
  `DF_Outfit_FACETS`, `InfitZSTD_FACETS`, and `OutfitZSTD_FACETS`.
- **FACETS fit review helper**: new `facets_fit_review()` separates
  engine-level fit-standardization differences from optional external
  FACETS table comparisons. The engine-vs-FACETS review compares engine and
  FACETS-style df/ZSTD values row-by-row and flags cases where the df
  convention changes the usual `|ZSTD| >= 2` screen. When a FACETS-like table
  is supplied, the external review matches rows by `Facet` / `Level` (or
  person labels for person-only
  tables) and classifies differences as `same`, `rounding`,
  `df_or_whexact_difference`, `mnsq_or_measure_difference`, or
  `needs_review`. This makes FACETS comparisons reproducible without treating
  mfrmr's package-native df convention as an error.
- **FACETS fit table import**: new `read_facets_fit_table()` /
  `import_facets_fit_table()` reads existing FACETS output into the
  `Facet` / `Level` / `Infit` / `Outfit` / `ZSTD` / df schema expected by
  `facets_fit_review()`. It supports already harmonized CSV/TSV-style tables,
  partial FACETS extracts with ZSTD and `TCount` but no MnSq/df columns, and
  FACETS `score.N.txt` files, including fixed-field score files using the
  FACETS manual column positions, with an optional `facet_map` for assigning
  score-file numbers to user-facing facet names. `facets_fit_review()` now
  returns `external_table_quality` so users can see duplicate
  `Facet` x `Level` rows and whether MnSq, ZSTD, df, and count columns were
  available in the supplied external table.
- **GPCM fair-average CI display**: `plot_fair_average(show_ci = TRUE)`
  no longer fabricates CIs from measure-level SEs for bounded `GPCM` fits.
  It now uses the opt-in structural fair-average SE columns when a fit object
  is supplied and records an unavailable-CI note for precomputed fair-average
  bundles that lack those columns.
- **GPCM bias likelihood checks**: `estimate_bias()` now adds
  conditional profile-likelihood columns for bounded `GPCM` rows:
  `LR ChiSq`, `LR d.f.`, `LR Prob.`, `Profile CI Lower`,
  `Profile CI Upper`, `Profile CI Level`, and `Profile CI Status`.
  These compare the fitted additive bias shift with zero while holding
  theta, steps, slopes, and other facet estimates fixed. They strengthen
  the GPCM screening evidence without turning it into standalone
  confirmatory fairness inference. `summary(estimate_bias(...))` now
  surfaces the profile-LR screen-positive count and carries these columns
  through the top-row review table when they are available.

## Research-grounded visualization refinements

- **Category-curve information output**: `category_curves_report()`
  now carries per-curve `ScoreVariance`, `Slope`, and `Information`
  columns. For `GPCM`, `Information` is computed as
  \(a^2 \mathrm{Var}(X \mid \theta)\), matching the Muraki/Samejima
  polytomous information identity; for `RSM` / `PCM`, this reduces to
  the usual score variance. `plot(category_curves_report(fit),
  type = "information")` returns the corresponding curve-level
  information plot data.
- **Bias heatmap review data**: `plot_bias_interaction(plot =
  "heatmap", draw = FALSE)` now returns `heatmap_cells`,
  `heatmap_matrix`, flag/count matrices, interpretation guidance, and
  reference notes. The display is documented as a FACETS Table 13-style
  screening follow-up, not confirmatory evidence.

## Recovery simulation workflow

- **`evaluate_mfrm_recovery()`** adds a dedicated parameter-recovery
  simulation route. It repeatedly simulates from a known MFRM generating
  setup, refits the requested model, and returns row-level truth/estimate
  comparisons plus summaries by parameter type. Location-like parameters
  are mean-aligned within replication before reporting recovery `Bias`,
  `RMSE`, `MAE`, correlation, and 95% coverage where standard errors are
  available; bounded-`GPCM` slopes are compared on the identified log-slope
  scale after the generator and fitter have both imposed the geometric-mean-one
  slope convention. The output also carries ADEMP-style simulation-study
  metadata so recovery checks are separated from broader design-evaluation
  claims.
- **`plot(evaluate_mfrm_recovery(...))`** adds review plots for recovery
  summaries, coverage, row-level error distributions, truth-estimate scatter,
  and replication status. `draw = FALSE` returns an `mfrm_plot_data` object
  with reusable plot tables and notes.
- **`assess_mfrm_recovery()`** adds a user-facing adequacy checklist for
  recovery simulations. It separates run completion, convergence, uncertainty
  availability, coverage, Monte Carlo precision, and optional practical
  RMSE/Bias thresholds into `ok` / `review` / `concern` style statuses with
  next-action text. `plot(assess_mfrm_recovery(...))` now provides checklist
  status-count and parameter-metric review plots so users can see which part of
  the assessment needs attention before reading the full tables. The
  `draw = FALSE` plot data include `reading_order`, `guidance`, and
  user-facing handoff tables such as `section_status` so follow-up starts with
  review/concern rows rather than raw row-level output.
- **Simulation refit score support**: simulation-based refit helpers now pass
  the generator's declared `1:score_levels` score support into `fit_mfrm()`.
  This keeps zero-count boundary categories in the fitted support during
  recovery, design-evaluation, diagnostic-screening, and bias-screening runs,
  and avoids repeated rating-range inference messages in release-validation
  logs.
- **Compact step-threshold specifications**: `build_mfrm_sim_spec()` and
  `simulate_mfrm_data()` now accept step-facet-specific thresholds as a named
  list or row-named numeric matrix, in addition to the existing long
  `StepFacet` / `StepIndex` / `Estimate` table.
- **Design-evaluation progress control**:
  `evaluate_mfrm_design(progress = interactive())` now shows the progress bar
  only in interactive sessions by default. Non-interactive tests, Quarto
  rendering, and batch scripts stay quiet unless users set `progress = TRUE`;
  users can also set `progress = FALSE` for fully silent exploratory runs.
- **Release recovery-validation protocol**: `inst/validation/recovery-validation.R`
  provides an optional long-run validation script for release review. It defines
  structured review steps, core `RSM` / `PCM` / bounded-`GPCM` recovery cases,
  an extended latent-regression case, practical thresholds, and a summary
  writer that produces top-line release-decision, case-level release-decision,
  review-step, case-plan, case-summary, metric-summary, overall decision-table,
  domain decision-table, run-note, RDS, and Markdown outputs without adding
  heavy Monte Carlo runs to routine package tests. The release decision uses
  recovery metrics, convergence, and Monte Carlo precision as the primary
  evidence; the domain decision table separately reports uncertainty status so
  missing JML coverage columns are not mistaken for recovery failure. Printing
  the validation object or calling `summary(validation)` now shows the
  release-level decision and case-level statuses before the full tables.
- **Recovery reporting handoff**: `build_summary_table_bundle()` now accepts
  `evaluate_mfrm_recovery()` and `assess_mfrm_recovery()` outputs directly,
  including ADEMP-style methods metadata, replication status, checklist rows,
  metric review rows, thresholds, notes, and appendix-preset roles.
- **Recovery appendix export**: `export_summary_appendix()` now recognizes
  recovery simulation and recovery assessment objects as direct inputs. The
  workflow vignette shows the full sequence from simulation specification to
  recovery plots, adequacy assessment, summary-table bundle, and appendix
  export.
- **Research-grounded release evidence map**:
  `inst/validation/release-evidence-map-0.2.0.md` gives a source-based
  review plan for 0.2.0. It links the release checks to Andrich's `RSM`,
  Masters' `PCM`, Muraki's `GPCM` and information-function work,
  FACETS/Winsteps fit conventions, and ADEMP-style simulation-study reporting,
  then separates release-gate checks from future-scope items. The
  companion `release-evidence-checklist-0.2.0.csv` provides a structured
  required / caveat / future-scope checklist for release review.
  `external-parameter-recovery-simulation-0.2.0.md` summarizes the separate
  common-data recovery and cross-engine agreement workflow and its limits without bundling the
  generated simulation datasets; the validation bundle also includes a sourceable
  helper for re-reading a local external-output directory when that workflow is refreshed.
- **Release-readiness protocol**:
  `inst/validation/release-readiness.R` turns the evidence map into a
  reproducible review object. It records eight review steps, parses an
  `R CMD check` log, checks the 0.2.0 version contract, verifies the CI
  workflow contract for warning failures and retained check artifacts, scans
  public docs for disallowed removed-helper wording, confirms evidence
  artifacts, and reports a top-line `ok` / `review` / `concern` gate summary
  without adding exported user-facing API.

## Citation and attribution corrections

- **Muraki DOI consistency**: `DESCRIPTION` now cites Muraki (1992) using
  the GPCM article DOI, `10.1177/014662169201600206`. The Muraki (1993)
  Applied Psychological Measurement reference in GPCM information help pages
  now uses `10.1177/014662169301700403`.
- **Wright (1998) page**: `R/api-shrinkage.R` references corrected from
  *Rasch Measurement Transactions*, 12(2), **638** to **632-633**
  (page 638 in the same RMT issue is a different paper; verified at
  <https://www.rasch.org/rmt/rmt122.htm>).
- **Linacre (1989, "2004") in `reporting_checklist()`**: the bare
  "Linacre (1989, 2004)" tag in `R/api-reporting-checklist.R` is now
  **Linacre (1989, 2002)**. The 2002 paper is "Optimizing rating scale
  category effectiveness," *JAM*, 3(1), 85-106 -- the canonical Linacre
  reference for rating-scale guidance. (No bibliographic entry existed
  for "Linacre (2004)".)
- **Eckes (2011) full reference**: the inline `(cf. Eckes, 2011; ...)`
  caveat in `dif_report()` now has a complete `@references` entry
  pointing to *Introduction to Many-Facet Rasch Measurement* (1st ed.,
  Peter Lang). McNamara & Knoch (2012) is also fully cited.
- **Mean-square fit ranges** in `?mfrmr-package` previously attributed
  the context-specific bands (high-stakes / clinical / survey) to
  Linacre (2002). The actual source is **Wright & Linacre (1994)**,
  *RMT* 8(3), 370. The band assignments were also swapped: high-stakes
  MCQ is **0.8-1.2** (not 0.6-1.4), survey is 0.6-1.4, clinical
  observation is 0.5-1.7. Corrected.
- **Yen Q3 (`q3_statistic()`)**: previously stated mfrmr's Q3 uses
  standardized residuals as if matching Yen (1984). Yen's eq. 7 (p. 127)
  uses **raw** residuals; mfrmr's standardized-residual choice is now
  documented as a deliberate departure. The `|Q3| > 0.20` cutoff was
  attributed to Yen but is from **Chen & Thissen (1997)**, *JEBS*,
  22(3), 265-289. Re-attributed.
- **Christensen et al. (2017) in `q3_statistic()`**: the central
  finding of Christensen et al. is that **no single critical value is
  appropriate** across designs and that a parametric bootstrap should
  be used. Documentation now states this clearly; the fixed
  `relative_offset = 0.20` is described as a screening default rather
  than as a re-implementation of `Q3_*`.
- **Morris (1983) posterior-SE correction formula** in
  `?apply_empirical_bayes_shrinkage` was dimensionally wrong:
  previously written `2 B^2 (tau^2 + SE^2)^2 / (K - 3)`, which is
  SE^4-units. The actual Morris (1983, eq. 4.1-4.2, p. 51) correction
  is `(2 / (K - r - 2)) * B^2 * delta^2`. Corrected, with re-derived
  magnitude examples (SE understated by ~73% at K=3, ~29% at K=5,
  ~7% at K=15).
- **Koo & Li (2016) ICC band boundary**: `compute_facet_icc()`
  previously placed ICC = 0.9 in **Excellent** (`>= 0.9`). Koo & Li
  (2016, p. 161) write "values **greater than 0.90** indicate excellent
  reliability" -- strict `>`. Code at `R/api-hierarchical-audit.R`
  now uses `> 0.9` for Excellent; ICC = 0.9 reads as Good.

## Documentation refinements

- **Linacre FACETS / Winsteps manuals**: cited years updated from 2023
  / 2024 to **2026** (current FACETS 4.5.0 = April 2026, Winsteps 5.11.0
  = March 2026 per <https://www.winsteps.com/index.htm>).
- **Bock & Aitkin (1981) clarification**: `?mfrmr-package` now notes
  that the default `mml_engine = "direct"` optimises the marginal
  log-likelihood by gradient methods (BFGS / L-BFGS-B), not by Bock &
  Aitkin's signature EM. The `"em"` and `"hybrid"` engines follow the
  EM template but with a BFGS M-step (rather than B&A's probit IRLS),
  because the target is the polytomous Rasch family rather than 2PL.
- **Linacre (1994) sample-size bands**: `mfrm_core.R` and `reporting.R`
  now describe the bands as "adapted from Linacre (1994)" rather than
  "follow Linacre (1994)". Only the 30-examinee floor is Linacre's;
  the `< 10 sparse` and `< 50 standard` watermarks are mfrmr-specific
  screening choices.
- **Snijders (2001) lz\\* correction**: `compute_person_fit_indices()`
  now computes the Snijders weight-projection correction for
  JML/fixed-effect person estimates, conditional on the fitted
  non-person calibration. The implementation uses the polytomous form
  `w_tilde_k = log(P_k) - c_n d log(P_k) / d theta`, with
  `c_n = Cov(log P, score) / I(theta)`. MML/EAP person scores keep
  `lz_star = NA` with `lz_star_status = "not_applicable_eap"` because
  EAP does not satisfy the ML/MAP/WLE estimating-equation setup.
- **Report-ready person-fit output**: `compute_person_fit_indices()` now
  adds practical 5% / 1% flag columns and compact `ReportIndex`,
  `ReportValue`, `ReviewStatus`, `ReviewReason`, and `ReportCaveat`
  columns. `ReportIndex` uses `lz_star` only when the Snijders correction
  was actually computed; otherwise it falls back to uncorrected `lz` with
  the status caveat left visible. `plot_person_fit()` now carries these
  person-fit indices in its draw-free plot data, adds reusable `plot_long`
  and `flag_summary` tables, and supports `fit_index = "loglik"` plus
  `preset = "monochrome"` for report-focused person-fit displays.
- **Person-fit summary and table-bundle handoff**:
  `compute_person_fit_indices()` now returns an `mfrm_person_fit_indices`
  data-frame subclass. `summary(person_fit)` gives overview counts,
  `ReviewStatus` / `ReportIndex` / `lz_star_status` summaries, top review
  rows, thresholds, caveats, and a reporting map. The same summary is now
  accepted by `build_summary_table_bundle()`, so person-fit review rows and
  Snijders-availability caveats can move into appendix/report workflows
  without custom table wrangling.
- **Marais (2013) `|Q3| > 0.30`**: documented as a community convention
  Marais cites, not as her own recommendation; her actual recommendation
  is the relative-to-mean comparison.

## Default changes

No defaults change between 0.1.6 and 0.2.0. The 0.1.6 defaults
(`quad_points = 31`, `diagnostic_mode = "both"`,
`plot.mfrm_fit(type = "wright")`, `keep_original = FALSE`) are retained.

Note for users upgrading directly from CRAN 0.1.5 to 0.2.0 (skipping
intermediate 0.1.6 builds): three defaults were flipped in 0.1.6
and remain on those values in 0.2.0 -- `diagnose_mfrm(diagnostic_mode)`
went from `"legacy"` to `"both"`, `plot(fit)` returns the Wright map
alone instead of a three-plot overview (the overview is still
available via `plot(fit, type = "bundle")`), and `fit_mfrm(quad_points)`
went from `15` to `31`. See the "mfrmr 0.1.6" section below for the
full description and revert paths.

## New features

### Continuous integration

New GitHub Actions workflows added alongside the existing
`pkgdown.yaml`: `R-CMD-check.yaml` runs the matrix on Ubuntu
(release / devel / oldrel-1) plus macos-latest and windows-latest
(release), and `test-coverage.yaml` runs `covr` with artifact
upload (no external service contacted).

### Differential-functioning display controls

`plot_dif_heatmap()` gains display controls for cell labels
(`show_values`, `value_digits`), absolute flag thresholds
(`flag_threshold`, `flag_color`), and shared symmetric color limits
(`scale_limit`) so several heatmaps can be drawn on a comparable scale.

`plot_dif_summary()` gains optional normal-approximation confidence
intervals, effect-threshold guide lines, method-aware axis labels, and
interpretation-guide data that downstream code can render
alongside the figure.

### Plot data printing

`print.mfrm_plot_data()` is now defined, so the headline `draw = FALSE`
return value renders as a compact summary (name, title, reusable data
shapes, legend / reference-line counts) instead of a raw list dump.

### Bounded GPCM fair-average and bias unblock (slope-aware)

`fair_average_table()` and `estimate_bias()` no longer hard-stop on
`GPCM` fits. Both helpers now use the slope-aware element-conditional
GPCM construction:

- **`fair_average_table()`**: for slope-facet element rows, the
  fair-average uses that element's own discrimination `a_{j*}` and
  threshold structure: `FA_{p,j*} = sum_k k * P_GPCM(X = k | theta_p,
  a_{j*}, delta_{j*})`. For non-slope facets (Person, Rater, ...), the
  fair-average uses the geometric-mean-one slope by GPCM
  identification, so the construction is continuous with the PCM
  Linacre fair-average and reduces to it exactly when all slopes
  equal one (regression-tested at machine precision).

- **`estimate_bias()`**: the per-cell bias parameter is the additive
  shift on the linear predictor that maximises the per-cell GPCM
  log-likelihood. The dispatch routes the inner `nll` and the
  per-iteration `category_prob` calls through the GPCM kernel instead
  of the PCM kernel; SE / t / Prob columns retain the screening-tier
  semantics documented in `?estimate_bias`.

Both helpers gain `method = "GPCM-slope-aware"` and a `caveat`
field that names the slope convention. For fair averages, the original
SE columns remain measure-level SEs, while `fair_se = TRUE` adds
structural delta-method fair-average SEs for non-person rows when the
MML Hessian is available. For bias values, the SE / t / Prob columns
retain their conditional screening interpretation.
See `?fair_average_table`, `?estimate_bias`, and
`gpcm_capability_matrix()` for the full support contract.

`build_apa_outputs()`, `facets_output_contract_review()`, and
`facets_output_file_bundle(include = "score")` remain blocked under
GPCM in 0.2.0; they require the same SE infrastructure to ship as
publication-quality outputs.

## Bug fixes

- `compute_person_fit_indices()` now computes `lz` from the model
  category probability of the observed category directly (true
  Drasgow, Levine & Williams (1985) polytomous form), via three new
  intermediate columns `PrObserved`, `ItemEntropy`, and `ItemVarLogP`
  on `compute_obs_table()`. The previous Gaussian-residual
  approximation overstated `Var[log P]` by roughly a factor of five
  on a 4-category fixture and pulled `lz` toward zero.
- The `ECI4` column is removed from `compute_person_fit_indices()`.
  The previous implementation was the standardized chi-square
  `(sum StdSq - n) / sqrt(2 * n)`, which is the linear (Smith)
  approximation to `OutfitZSTD`, not the Tatsuoka & Tatsuoka (1983)
  extended-caution index. Users who want the equivalent statistic
  should use `OutfitZSTD` directly. `lz_star` now uses the Snijders
  (2001) weight-projection correction where its estimating-equation
  assumptions are met, and otherwise stays `NA` with an explicit status.
- `displacement_table()$summary` now returns `NA_real_` for
  `MaxAbsDisplacement` and `MaxAbsDisplacementT` when every flagged
  level has zero information (so every `Displacement` is `NA`).
  Previously the helper called `max(..., na.rm = TRUE)` on an
  all-`NA` vector, which returned `-Inf` and emitted a "no
  non-missing arguments to max; returning -Inf" warning. The
  guarded version is regression-tested in `test-core-coverage-gaps.R`.
- `analyze_dff()` and `dif_interaction_table()` now reject invalid
  `p_adjust`, non-integer `min_obs`, invalid `focal` groups, and
  all-missing group columns up front, instead of failing later inside
  the contrast computation. Missing or empty group rows are dropped
  with a `message()`.

## Documentation

- `?analyze_dff`, `?plot_dif_summary`, `?mfrmr_linking_and_dff`,
  and `?mfrmr_visual_diagnostics` now distinguish residual-method
  screening labels from refit-method ETS A/B/C classifications more
  explicitly and route users to both `plot_dif_heatmap()` and
  `plot_dif_summary()`.

- `?compute_person_fit_indices` now describes when `lz_star` is computed
  and when it is intentionally left `NA`: JML/fixed-effect person scores
  receive the Snijders (2001) correction, whereas MML/EAP scores remain
  uncorrected because EAP is outside the Snijders estimating-equation
  setup.

- `?mfrm_generalizability` now discloses that the lme4 random-effects
  model is main-effects only (`Score ~ 1 + (1|Person) + (1|Facet) +
  ... + Residual`, no explicit `(1|Person:Facet)` interaction terms),
  which folds two-way interaction variance into Residual and can
  bias `G` downward. The companion `mfrm_d_study()` projects `G` and
  `Phi` under planned facet counts, but reports the residual-scaling
  assumption explicitly; users who need a full p x r x i decomposition
  should treat these projections as planning evidence, not as a
  substitute for separately estimated interaction components.

- `?q3_statistic` now discloses that, when the chosen facet has
  multiple residual rows per (Person, Level) cell because of
  additional facets in the design, the standardized residuals are
  mean-aggregated to one value per cell before the Pearson
  correlation. Yen's (1984) original definition takes the
  correlation over per-(Person, Item) residuals without aggregation,
  so the published `|Q3| > 0.20` threshold and the Christensen et
  al. (2017) critical values were derived for the original
  formulation; the values returned here should be treated as a
  screening summary rather than a direct substitute for those
  thresholds.

- `?bias_pairwise_report` now discloses that the contrast SE uses
  the independence approximation `sqrt(SE_i^2 + SE_j^2)`. For
  same-facet bias values that share a sum-to-zero identification
  the true `Cov(b_i, b_j) < 0`, so the reported SE is an
  over-estimate and the t-statistic / p-value are conservative
  (the true significance is higher than reported). For across-facet
  contrasts the covariance term is approximately zero and the
  approximation is appropriate.

- Two new vignettes ship in the `Migration and Scope` section of the
  pkgdown article navigation: `vignette("mfrmr-facets-migration")`
  walks Facets users through the corresponding `mfrmr` workflow and
  numeric contract checks, and `vignette("mfrmr-gpcm-scope")` documents
  which downstream helpers the bounded `GPCM` route currently
  supports versus restricts and what to use as a substitute when a
  helper is restricted.

## Build hygiene

`.Rbuildignore` tightened the `inst/references/` source-package boundary.
The two runtime / user-facing files in that directory --
`facets_column_contract.csv` (read at runtime by
`facets_output_contract_review()`) and `FACETS_manual_mapping.md` (the
FACETS Table to `mfrmr` helper mapping cited in the README) -- are
preserved.

## Performance note

The cpp11 MML backend (`src/mml_backend.cpp`, RSM and PCM only) is
opt-in via `options(mfrmr.use_cpp11_backend = TRUE)` for this release.
It is validated against the pure-R reference at `tolerance = 1e-12`
on a fixed regression fixture. The default flip to ON is planned for
a follow-up release after a cycle of community testing.

## Deferred to a follow-up release

Considered for 0.2.0 but not shipped in 0.2.0; carried over to a
later release:

- User-facing GPCM unblock for `build_apa_outputs()`,
  `facets_output_contract_review()`, and `facets_output_file_bundle(include =
  "score")`. (`fair_average_table()` and `estimate_bias()` are
  unblocked above.)
- A classical-DIF helper (working title `analyze_dif_classical()`)
  covering Mantel-Haenszel, logistic regression, and SIBTEST.
- Five additional Rasch / IRT classic plots (KIDMAP, TCC, expected
  score curve, cumulative ICC, information surface).
- A native classical-DIF vignette (the migration and bounded-GPCM-scope
  vignettes ship in this release; see the Documentation section above).

These are scheduled for a follow-up release.

# mfrmr 0.1.6

This release adds empirical-Bayes shrinkage for small-N facets, a
hierarchical-structure and sample-adequacy review layer, integrated
missing-code pre-processing, APA output adapters for Word / HTML,
model-estimated two-way non-person facet interactions, confidence-interval
propagation through the plot surface and the ICC
reporting family, and expanded reproducibility manifests. Six bug
fixes close issues that affected bias statistics, ZSTD sign, input
validation, and graphical state hygiene.

## Default changes (three breaking flips)

Three default values change in this release. Scripts that explicitly
pass the old value are unaffected; scripts that rely on the default
should be reviewed.

- `diagnose_mfrm(diagnostic_mode = ...)` default flips from `"legacy"`
  to `"both"`. Strict marginal screens are produced automatically for
  `RSM` / `PCM` fits without the caller having to request them.
  Pass `diagnostic_mode = "legacy"` to restore the earlier behaviour.
- `plot(fit)` default output is now the Wright map alone, returned as
  an `mfrm_plot_data` object. The previous three-plot overview
  (Wright + pathway + CCC) remains available via
  `plot(fit, type = "bundle")`, which returns an `mfrm_plot_bundle`
  with the same three slots.
- `fit_mfrm(quad_points = ...)` default increases from `15` to `31`
  so a default MML fit is stable enough for direct manuscript
  reporting. Pass `quad_points = 15` (or `7`) to restore the earlier
  iteration speed for exploratory scans.

## New features

### Model-estimated facet interactions

`fit_mfrm()` gains `facet_interactions` for confirmatory two-way interactions
between non-person facets in `RSM` and `PCM` fits, for example
`facet_interactions = "Rater:Criterion"`. These terms are estimated
simultaneously with the main MFRM parameters as fixed effects under zero
marginal-sum constraints, contributing `(A - 1) * (B - 1)` free parameters for
an `A x B` interaction block.

New supporting pieces:

- `interaction_effect_table(fit)` returns one row per interaction cell, with
  estimates, weighted counts, sparse-cell flags, and the identification note.
- `summary(fit)` reports a compact interaction overview when interaction terms
  are present.
- `compare_mfrm(..., nested = TRUE)` now recognizes same-family additive-vs-
  interaction comparisons as nested when all other structural settings match
  and the smaller model's interaction set is a subset of the larger model's
  set.

The feature is intentionally narrow for the initial CRAN-facing release:
person-involving interactions, higher-order interactions, GPCM interactions,
and random-effect facet interactions are deferred. Residual bias screening via
`estimate_bias()` and `estimate_all_bias()` remains separate from these
model-estimated fixed effects.

### Empirical-Bayes facet shrinkage

`fit_mfrm(..., facet_shrinkage = "empirical_bayes")` applies
James-Stein / empirical-Bayes shrinkage to each non-person facet's
fixed-effect estimates. `fit$facets$others` gains `ShrunkEstimate`,
`ShrunkSE`, and `ShrinkageFactor` columns, and `fit$shrinkage_report`
summarises the per-facet prior variance, mean shrinkage, and
effective degrees of freedom.

The estimator is the classical method-of-moments form (Efron & Morris,
1973):

- `tau_hat^2 = max(0, mean(delta_hat_j^2) - mean(SE_j^2))`, using the
  raw second moment under mfrmr's sum-to-zero identification (the
  facet mean is exactly 0 by construction, so no degree of freedom is
  consumed).
- `B_j = SE_j^2 / (tau_hat^2 + SE_j^2)` (shrinkage factor).
- `delta_hat_j^EB = (1 - B_j) * delta_hat_j` and
  `SE_j^EB = sqrt((1 - B_j) * SE_j^2)` (posterior mean / SE; the
  posterior SE treats `tau_hat^2` as known, omitting the Morris
  (1983) correction for `tau_hat^2` uncertainty).

Two post-hoc helpers make shrinkage available to existing fits:

- `apply_empirical_bayes_shrinkage(fit, facet_prior_sd = NULL,
  shrink_person = FALSE)` augments an existing `mfrm_fit`.
- `shrinkage_report(fit)` returns the per-facet summary table.

The `"laplace"` alias currently routes to the empirical-Bayes path
and is reserved for a future penalised-MML implementation.

Integration: `summary(fit)` exposes `FacetShrinkage` and
`FacetShrinkageTau2Mean`; `build_apa_outputs()` adds a Method-section
sentence naming the mode, mean `tau_hat^2`, and mean shrinkage
with a Efron & Morris (1973) citation; `build_mfrm_manifest()` gains
a `shrinkage_audit` table; `reporting_checklist()` gains an
"Empirical-Bayes shrinkage" item.

### Hierarchical structure and sample-adequacy review

Five new exported functions describe the observed design, flag
small-N facet levels, and quantify ICC / design effect. Estimation
remains fixed-effects MFRM; these helpers are purely descriptive and
do not alter the fit.

- `detect_facet_nesting(data, facets, person)` classifies every
  ordered pair of facets (plus Person, optionally) as *Fully nested*,
  *Near-perfectly nested*, *Partially nested*, or *Crossed* using the
  conditional-entropy index `1 - H(B|A)/H(B)`.
- `facet_small_sample_review(fit)` returns per-level
  `N / Estimate / SE / Infit / Outfit / SampleCategory` for every
  facet. `SampleCategory` is one of `"sparse"` (< 10), `"marginal"`
  (< 30), `"standard"` (< 50), `"strong"` (>= 50). Thresholds follow
  Linacre (1994) and are configurable.
- `compute_facet_icc(data, facets, score, person)` fits
  `lme4::lmer(Score ~ 1 + (1|Person) + (1|Facet1) + ...)` and reports
  the variance-component share per facet. Person uses the Koo & Li
  (2016) reliability bands; other facets use a "variance share" label
  (Trivial / Small / Moderate / Large).
- `compute_facet_design_effect(data, facets, icc_table)` computes the
  Kish (1965) `Deff = 1 + (m - 1) * rho` and effective N per facet.
- `analyze_hierarchical_structure(data, facets, ...)` bundles the
  four helpers above and (when `igraph` is available) a bipartite
  connectivity summary over Person * facet-level edges.

Fit- and reporting-stack integration:

- `fit$summary` carries `FacetSampleSizeFlag`, `FacetMinLevelN`, and
  `FacetSparseCount`.
- `reporting_checklist()` gains two items: "Facet sample-size
  adequacy" (auto-ready when the flag is `"standard"` / `"strong"`)
  and "Hierarchical structure review" (ready when the user passes
  `hierarchical_structure = analyze_hierarchical_structure(...)`).
- `build_apa_outputs()` adds a Method sentence naming the
  sample-adequacy band and linking to `facet_small_sample_review()`.
- `build_mfrm_manifest()` gains a `hierarchical_audit` table.
- `recommend_mfrm_design()$caveats` now points users at the three
  post-fit audit functions.

Optional dependencies `igraph` and `lme4` move to `Suggests`; when
either is absent the relevant report is omitted with a clear
`message()`.

### Missing-code pre-processing in the fit call

`fit_mfrm()` now accepts `missing_codes = NULL | TRUE | "default" |
<character vector>`, forwarded to `prepare_mfrm_data()`,
`review_mfrm_anchors()`, and `describe_mfrm_data()`. When active, the
standard FACETS / SPSS / SAS sentinels (`"99"`, `"999"`, `"-1"`,
`"N"`, `"NA"`, `"n/a"`, `"."`, `""` by default, or any caller-
supplied set) are converted to `NA` on the `person`, `facets`, and
`score` columns before any downstream processing. Replacement counts
are recorded in `fit$prep$missing_recoding` and surfaced through
`build_mfrm_manifest()$missing_recoding`. The default
(`missing_codes = NULL`) is strictly backward-compatible.

A standalone `recode_missing_codes()` helper is also exported for
users who prefer to recode before calling `fit_mfrm()`.

### APA output adapters

- `as_kable.apa_table()` converts an `apa_table` into a
  `knitr::kable()` object with the caption above and the note below.
  When `kableExtra` is installed the note becomes a proper table
  footnote; otherwise it is appended as Markdown.
- `as_flextable.apa_table()` produces a `flextable::flextable()`
  with caption and note pre-wired, suitable for `officer` / Word /
  PowerPoint exports.
- Two generics, `as_kable()` and `as_flextable()`, are exported so
  other mfrmr classes (or third-party wrappers) can register
  compatible methods.
- `build_apa_outputs(..., context = list(output_mode = "reflow"))`
  now returns the Method / Results paragraphs as single long lines
  per sentence-joined paragraph, which is the format Word / Quarto /
  RMarkdown prefer. The default `"wrapped"` keeps the 92-column
  layout for console readability.

`kableExtra` and `flextable` join `Suggests`.

### Shrinkage and review visualisations

- `plot(fit, type = "shrinkage")` renders a horizontal forest-style
  dotplot of original and shrunk facet-level estimates, with arrows
  indicating shrinkage direction, optional 95 % CI error bars
  (`show_ci = TRUE`), and a reference line at zero. When shrinkage
  is not applied the plot shows an unavailable-state message inviting
  the user to re-fit with `facet_shrinkage = "empirical_bayes"`.
- `plot.mfrm_facet_sample_review()` draws a horizontal bar chart of
  per-level observation counts coloured by Linacre band, with dashed
  vertical lines at the thresholds.
- `plot.mfrm_facet_nesting()` renders the pairwise nesting index as
  a heatmap with numeric cell labels.

All three methods follow the existing
`preset = c("standard", "publication", "compact")` convention and
use base-R graphics.

### Confidence intervals across the plot surface

- `plot_bias_interaction(show_ci = TRUE, ci_level = 0.95)` draws
  `BiasSize +/- z * SE` whiskers on the scatter and ranked views.
- `plot_displacement(show_ci = TRUE)` draws
  `Displacement +/- z * DisplacementSE` whiskers in the lollipop
  view.
- `plot_fair_average(show_ci = TRUE)` draws fair-average CI whiskers
  on the observed-score scale using a delta-method propagation
  `SE_fair = Var(X | Measure) * ModelSE` from the logit `Measure`
  error. Rows near a rating boundary (where the implied score
  variance is effectively zero) are excluded from the whiskers,
  drawn as open circles, and counted in the subtitle.
- `compute_facet_icc(ci_method = "profile" | "boot")` returns ICC
  confidence intervals in new `ICC_CI_Lower` / `ICC_CI_Upper` /
  `ICC_CI_Level` / `ICC_CI_Method` columns, propagated through
  `analyze_hierarchical_structure()` and drawn as whiskers on
  `plot.mfrm_hierarchical_structure(type = "icc")`. The default
  `ci_method = "none"` keeps the point-estimate-only behaviour.

### Additional visualisations

Fourteen additions across the plot surface, all base-R / additive
(default behaviours unchanged):

- **`plot_threshold_ladder()`** (new) — vertical ladder of
  Rasch-Andrich thresholds for RSM and PCM, with disordered-step
  crossings highlighted in the preset's `fail` colour. The returned
  object includes per-step `Group / Step / Threshold / Disordered`
  rows.
- **`plot(fit, type = "ccc_overlay")`** (new branch on
  `plot.mfrm_fit`) — observed category proportions binned by person
  measure overlaid on the model CCC curves, for an at-a-glance
  model-data fit visual.
- **`plot_person_fit()`** (new) — FACETS Table 6 style per-person
  Infit / Outfit bubble with the standard 0.5-1.5 acceptance band
  (Linacre, 2002).
- **`plot_bias_interaction(plot = "heatmap")`** (new mode) — diverging
  Rater x Criterion grid coloured by bias size, with flagged cells
  outlined for emphasis.
- **`plot(fit, type = "wright", group = ..., group_data = ...)`**
  (new option) — overlays per-group person-density curves on the
  Wright map's left density column, useful for DIF / DFF screening.
- **`plot_rater_severity_profile()`** (new) — per-rater severity
  ranking with CI whiskers and optional `+/-0.5` (gentle) and
  `+/-1.0` (strict) guidance bands for rater-training feedback.
- **`plot_anchor_drift(type = "forest")`** (new mode) — per-wave
  anchor-element CI forest with point estimate + `z * SE` whiskers.
- **`plot.mfrm_equating_chain()`** (new S3 method) — `type =
  "common_anchors"` (default bar chart of pairwise common-anchor
  counts) and `type = "graph"` (bipartite Wave x anchor element
  graph via `igraph`).
- **`plot_apa_figure_one()`** (new) — 2x2 publication composite
  bundling Wright map, rater severity profile, threshold ladder, and
  a one-panel summary block.
- **`plot_dif_summary()`** (new) — compact effect-size summary for
  [analyze_dff()] / [analyze_dif()] with ETS A / B / C colour coding.
- **`plot_guttman_scalogram()`** (new) — Person x facet-level
  observed-category matrix, ordered by person measure and location
  measure, with unexpected cells highlighted.
- **`plot_residual_qq()`** (new) — normal Q-Q plot of person-level
  standardized residuals for distributional misfit diagnostics.
- **`plot_rater_trajectory()`** (new) — per-rater severity
  trajectory across an ordered wave / session variable with CI
  whiskers. Accepts a named list of fits.
- **`plot_rater_agreement_heatmap()`** (new) — symmetric rater x
  rater agreement matrix colored by exact agreement (default) or the
  Pearson-style `Corr` column from `interrater_agreement_table()`.
  Quadratic-weighted kappa is not currently computed by that helper
  and is therefore not exposed as a `metric` option.

`igraph` is already in `Suggests`; the equating-graph view falls
back to the bar chart when `igraph` is not installed.

### Expanded test coverage

Direct regression tests for the 0.1.6 additions:

- `test-attach-diagnostics.R` — 18 assertions covering the
  `attach_diagnostics = TRUE` merge, type validation, idempotence,
  and MML / JML agreement checks.
- `test-icc-ci-method.R` — 25 assertions covering
  `compute_facet_icc(ci_method = "profile" / "boot")`, bootstrap
  seed reproducibility, range validation, deprecated
  `icc_ci_method` alias, and `plot.mfrm_hierarchical_structure(type
  = "icc")` integration.
- `test-ci-api-consistency.R` — 21 assertions covering the
  `lifecycle::deprecate_warn()` path for `conf_level`, `show_ci` /
  `ci_level` on `plot_fair_average` / `plot_displacement` /
  `plot_bias_interaction`, and CI column schema.
- `test-messaging-and-guards.R` — assertions covering the quiet-by-default
  inferred rating-range message, the opt-in one-time message,
  `analyze_dff(method = "refit")` `missing(diagnostics)` guard, and
  `missing_codes` integration.
- `test-lme4-confint-helper.R` — 17 assertions covering
  `.lme4_confint_components()` across terse and verbose lme4
  row-name conventions.
- `test-plotting-extras.R` + `test-plotting-screening.R` — 78
  assertions covering all 14 new plot helpers.

### Internal architecture

`row_max_fast()` and the three `category_prob_*` polytomous-response
kernels are now in `R/core-category-probabilities.R` instead of
inline in `R/mfrm_core.R`. Pure file-level reorganization; no
behaviour change. The remaining structural split of `mfrm_core.R`
(likelihood / optimizer / EM / gradients / prep / report tables) is
scheduled for a future release.

### Package-level MnSq misfit threshold

`mfrm_misfit_thresholds()` returns the lower / upper Linacre acceptance
band that mfrmr screens use when flagging element-level Infit / Outfit
MnSq misfit. Defaults are `c(lower = 0.5, upper = 1.5)` and can be
overridden globally via R options:

- `options(mfrmr.misfit_lower = 0.7)`
- `options(mfrmr.misfit_upper = 1.3)`

Helpers that consume the band include `summary(diagnose_mfrm(...))`
(`misfit_flagged` block + `key_warnings` auto-flag),
`build_misfit_casebook()` (the new `element_fit` source family),
the bias / misfit narrative inside `build_apa_outputs()`, and
`facet_quality_dashboard()` when `misfit_warn = NULL`. Setting the
options once at the top of an analysis script therefore changes
every downstream screen at once.

### Additional secondary plots

Four new public helpers extend the diagnostic plot family:

- `plot_local_dependence_heatmap(fit)` -- N x N Q3-style
  pairwise residual correlation heatmap between facet levels.
  Complements `plot_marginal_pairwise()` by showing every pair on a
  shared color scale rather than a top-N bar list.
- `plot_reliability_snapshot(fit)` -- compact facet x reliability /
  separation / strata bar overview built from
  `diagnostics$reliability`. Useful as a single small figure for "are
  persons / raters / criteria distinguishable?".
- `plot_residual_matrix(fit)` -- person x facet-level standardized
  residual heatmap. Complements `plot_guttman_scalogram()` by showing
  residual sign and magnitude rather than the raw response code.
- `plot_shrinkage_funnel(fit)` -- empirical-Bayes shrinkage
  caterpillar / funnel for fits augmented via
  `apply_empirical_bayes_shrinkage()`.

`plot_bubble()` gains a `view = c("measure", "infit_outfit")`
argument. The default `"measure"` keeps the historical Measure
(logit) x MnSq bubble layout; `view = "infit_outfit"` switches to the
Winsteps Table 30 layout (Infit MnSq on x, Outfit MnSq on y, bubble
size defaults to `N`). Both views return the same `mfrm_plot_data`
contract.

`plot_dif_heatmap(draw = FALSE)` now returns an `mfrm_plot_data` object
whose `data$matrix` is the metric matrix (was previously the bare
matrix only).

`plot_information(..., draw = FALSE)` outputs now include a
`series` field listing which curves the legend describes
(`"Information"`, `"SE"`, or both for `type = "both"`), so downstream
ggplot2 re-renderers can map the right column without inspecting
`type` manually.

### Reporting surface enrichments

- `summary(diagnose_mfrm(...))` now prints the **fixed-effect chi-square
  block** ("are all elements equal?") directly from
  `diag$facets_chisq` (`Facet`, `Levels`, `MeanMeasure`, `SD`,
  `FixedChiSq`, `FixedDF`, `FixedProb`, plus the random-effect
  counterparts when present) and the **inter-rater agreement summary**
  (Exact / Expected / Adjacent agreement, MeanAbsDiff, MeanCorr,
  RaterSeparation, RaterReliability) instead of leaving them in the
  diagnostics object only. The new `summary(diag)$facets_chisq` and
  `summary(diag)$interrater` slots also expose the same tables for
  programmatic use.
- `summary(diagnose_mfrm(...))$key_warnings` now names the worst
  MnSq-misfit elements (e.g. `MnSq misfit: Person:P023 (Infit=1.70,
  Outfit=2.40; outside 0.5-1.5).`) and prints a dedicated
  `MnSq misfit` block showing every flagged element. Threshold pair
  is exposed at `summary(diag)$misfit_thresholds` and is steered by
  `mfrm_misfit_thresholds()` (see above).
- `summary(diagnose_mfrm(...))` now prints a **category usage block**
  (one row per observed score with `Count`, `AvgMeasure`, and a
  `Disordering` flag when the average measure decreases across
  adjacent categories). Exposed programmatically at
  `summary(diag)$category_usage`.
- `summary(fit_mfrm(...))` now prints a **targeting block**
  (`Person mean - Facet mean`, plus `PersonSD` / `FacetSD` /
  `SpreadRatio`) for every non-person facet. Under the package's
  sum-to-zero identification this collapses to the person mean by
  construction; the row labels make that explicit and the spread
  ratio surfaces whether persons or facets dominate the test scale.
- `summary(estimate_bias(...))` now reports **Bonferroni** and
  **Holm** significant-cell counts alongside the raw screen-positive
  count. Both are exposed in `summary(bias)$overview` as
  `BonferroniSignificant` and `HolmSignificant`.
- `print(fit)` and `print(summary(fit))` now show an **"Attached
  diagnostics"** line when `fit_mfrm(..., attach_diagnostics = TRUE)`
  has merged per-element fit columns onto `fit$facets`. The
  attach-diagnostics path now extends to the person-facet table,
  so per-person `Infit`, `Outfit`, `InfitZSTD`, `OutfitZSTD`, and
  `PtMeaCorr` columns are visible in `summary(fit)$person_high`
  and `summary(fit)$person_low`.

### Internal architecture: file split

To improve navigability of the core estimation engine, four
self-contained sections moved out of `R/mfrm_core.R` into focused
files. All functions remain internal and the public API is
unchanged.

- `R/core-likelihood.R` -- polytomous Rasch likelihoods and
  cumulative response-probability helpers.
- `R/core-data-prep.R` -- data validation, indexing, and small
  formatting utilities.
- `R/core-anchor-review.R` -- anchor-table reading, normalization,
  and connectivity / overlap audit.
- `R/core-optimizer.R` -- optim() / EM dispatch and MML-EM
  scaffolding.

`R/api-simulation.R` similarly grew an
`R/api-simulation-future-branch.R` companion file holding the
future-branch design-schema layer. Public simulation entry points
(`simulate_mfrm_data`, `evaluate_mfrm_design`,
`evaluate_mfrm_diagnostic_screening`,
`evaluate_mfrm_signal_detection`) remain in `R/api-simulation.R`.
- `evaluate_mfrm_diagnostic_screening(include_report = TRUE)` can now retain
  `mfrm_results()` / `mfrm_report()` `report_index` signals at the replicate
  level and summarize them in `report_signal_summary`, keeping report-layer
  readiness separate from diagnostic-screening Type I and sensitivity proxies.
- `plot.mfrm_diagnostic_screening()` can now turn diagnostic-screening
  summaries into integrated plot-data bundles for overview rates, report
  signals, scenario contrasts, and runtime checks. The draw-free return follows
  the `mfrm_plot_data` / `plot_data()` contract.

`R/api-plotting-extras2.R` was renamed to
`R/api-plotting-screening.R` to drop the numerical suffix in favour
of a functional name; tests follow the same rename.

A new `tests/testthat/helper-fixtures.R` exposes
`make_toy_fit()` / `make_toy_diagnostics()` / `local_toy_fit()`
helpers so future tests can reuse the standard `example_core` fit
without retyping the `load_mfrmr_data()` + `fit_mfrm()` +
`diagnose_mfrm()` chain.

### Replay-script overhaul

`export_mfrm_bundle()` and `build_mfrm_replay_script()` now write a
self-contained replay package:

- The generated `replay.R` includes every argument that affected the
  original `fit_mfrm()` call. Earlier 0.1.x scripts silently dropped
  `missing_codes`, `mml_engine`, `slope_facet`, `anchor_policy`,
  `min_common_anchors`, `min_obs_per_*`, `facet_shrinkage`,
  `facet_prior_sd`, `shrink_person`, and `attach_diagnostics`, so
  fits that depended on those arguments did not actually replay.
- `fit_mfrm()` now records its inputs in
  `fit$config$replay_inputs` (post `match.arg`) so the bundle
  generator has a single source of truth.
- The replay script begins with a `utils::packageVersion("mfrmr")`
  guard that warns when the installed version differs from the
  recorded one.
- `export_mfrm_bundle(..., data = ...)` accepts the original analysis
  data; when supplied, the data is written into the bundle as
  `<prefix>_replay_data.csv` and the replay script reads from that
  co-located file. The recorded input hash is now computed against
  the user's original data (not the package's internal `prep$data`,
  which carries synthesised columns), so users can verify their CSV
  matches the recorded fingerprint.
- A new `tests/testthat/test-replay-roundtrip.R` actually sources the
  generated replay script in a fresh environment and compares the
  reproduced log-likelihood and person estimates to the original.

### Performance: `diagnose_mfrm()` on large designs

`calc_interrater_agreement()` (the inter-rater agreement helper that
`diagnose_mfrm()` calls when `Person` is part of `facet_cols`)
previously used a `list()` for the per-context probability lookup
and `c(exp_vals, ...)` accumulation inside a per-row loop. This
gave near-quadratic scaling: 6,400 observations took ~2 s, but
72,000 observations took ~141 s. The lookup is now an
`environment` (hash-backed for character keys) and `exp_vals` is
preallocated and filled by index, so the helper now scales linearly
in the number of observations. On the 72,000-observation benchmark
in the review, `diagnose_mfrm()` drops from ~141 s to ~15 s.

The `make_union_find()` helper used by the connectivity audit was
also rewritten with an iterative `find_root` (with path
compression) instead of the previous recursive form. Designs whose
union chain depth exceeded `options(expressions)` (default 5,000)
no longer error out with "evaluation is too deeply nested".

### Input validation: degenerate inputs surface earlier

`prepare_mfrm_data()` now:

- records how many rows it dropped due to missing values or non-positive
  weights, instead of dropping them silently;
- trims leading/trailing whitespace from `Person` and facet IDs
  and records the row count so " P01 " and "P01" do not silently become
  two persons;
- `warning()`s when the input contains duplicate Person x facet
  rows (which violate MFRM's conditional-independence assumption)
  but lets the fit continue rather than refusing it outright.

`fit_mfrm()` now treats `NaN` / `Inf` for `maxit`, `reltol`, and
`quad_points` as invalid input with a localised English error,
instead of falling through to R's locale-dependent
"missing value where TRUE/FALSE needed" message.

### Pre-rendered cheatsheet PDF

The two-page landscape cheatsheet now ships in pre-rendered form at
`system.file("cheatsheet", "mfrmr-cheatsheet.pdf", package = "mfrmr")`
alongside the existing `.Rmd` source. Users without a working LaTeX
toolchain can open the PDF directly; users who want to customize it
can still knit the `.Rmd` with `rmarkdown::render()`. The README and
`?mfrmr` package help now point at both files.

### Help-page examples: "what to look for" annotations

The most-visited help pages now embed concrete interpretation
comments inside their `@examples` blocks. Each shipped example
shows what value ranges or patterns indicate "good", what threshold
or rule of thumb applies, and what follow-up to run if the value
is off. Coverage in 0.1.6 includes:

- `?fit_mfrm` (convergence, person SD, targeting bands).
- `?diagnose_mfrm` (key_warnings, MnSq misfit lines, facets_chisq,
  inter-rater agreement minus expected).
- `?summary.mfrm_fit` and `?summary.mfrm_diagnostics` (overview,
  person distribution, top_fit ZSTD bands, facets_chisq, targeting).
- `?estimate_bias`, `?analyze_dff`, `?compute_facet_icc`,
  `?apply_empirical_bayes_shrinkage` (effect-size bands, Penfield
  classification, Koo & Li 2016 reliability bands, shrinkage factor
  interpretation).
- `?build_apa_outputs`, `?reporting_checklist`, `?plot_qc_dashboard`,
  `?plot.mfrm_fit` (manuscript-readiness signals, dashboard panel
  status, Wright / pathway / CCC interpretation).
- `?plot_bubble`, `?plot_dif_heatmap`, `?plot_local_dependence_heatmap`,
  `?plot_reliability_snapshot`, `?plot_residual_matrix`,
  `?plot_shrinkage_funnel`, `?plot_guttman_scalogram`,
  `?plot_residual_qq`, `?plot_rater_trajectory`,
  `?plot_rater_agreement_heatmap` (cell / band thresholds,
  reference-line interpretation).

### Help-page examples: lighter-weight `\donttest{}`

Several main entry points now expose a small fast-path block (a
`JML` fit on `example_core` plus a single diagnostic / plot call)
before the heavier `\donttest{}` block. The fast path is below
R CMD check's example-time budget and provides a regression net
that runs every check, while the full `\donttest{}` block
continues to showcase the larger MML / publication-route examples.
A final CRAN-timing pass keeps the active examples at lightweight
`maxit = 30` smoke fits so the printed examples demonstrate converged
objects without returning to the original long-running example surface.
Affected pages: `?fit_mfrm`, `?diagnose_mfrm`, `?plot_qc_dashboard`,
`?reporting_checklist`, `?build_apa_outputs`.

### Documentation

- `?mfrmr_visual_diagnostics` adds a "Cross-reference to FACETS /
  Winsteps tables" section that lists the closest mfrmr helper for
  each canonical Rasch / MFRM table or figure family (Wright map,
  pathway / probability curves, test information, misfit, bias /
  interaction, DIF / DRF, inter-rater agreement, anchoring /
  linking).
- `?mfrmr_visual_diagnostics` and the visual reporting template now
  enumerate the 4 secondary plot helpers and the 4 screening
  helpers added in 0.1.6.
- `?diagnose_mfrm` cites Wright & Masters (1982) at the
  separation / strata / reliability section and reproduces the
  formulae (G = TrueSD / RMSE, R = G^2 / (1 + G^2),
  H = (4G + 1) / 3) so the reliability outputs are traceable to
  source.
- `?fit_mfrm` example block now flags the `quad_points = 7` opening
  fit as an exploratory speed setting.
- The README and `?mfrmr` package help now point at the public
  cheatsheet (`system.file("cheatsheet", "mfrmr-cheatsheet.Rmd",
  package = "mfrmr")`).
- The bias / misfit APA narrative now spells out `|ZSTD|` (or
  `|MnSq - 1|` when ZSTD is unavailable) instead of the generic
  `|metric|` label.
- `build_misfit_casebook()` now also draws element-level Infit /
  Outfit MnSq misfit cases from `diagnostics$fit` (in addition to
  marginal cells, pairwise screens, unexpected responses, and
  displacement). The casebook therefore matches what its name
  implies.

### Yen Q3 local-dependence statistic

`q3_statistic(fit, diagnostics)` returns the Yen (1984) Q3 index
between every facet-level pair, with three published reporting
thresholds (Yen 0.20, Marais 0.30, Christensen et al. relative
0.20) and a textual `Interpretation` column that names which
flag(s) each pair triggered. The helper reuses the standardized-
residual pivot that `plot_local_dependence_heatmap()` already
draws, so the table and the heatmap stay numerically consistent.

### Extended person-fit indices

`compute_person_fit_indices(diagnostics, fit)` adds person-level
fit detail on top of the Infit / Outfit / ZSTD
columns that `diagnose_mfrm()` already exposes:

- **lz** (Drasgow, Levine & Williams, 1985): standardized
  log-likelihood under the fitted model.
- **lz\\*** (Snijders, 2001): estimated-ability correction computed for
  JML/fixed-effect person estimates conditional on the fitted non-person
  calibration; returned as `NA` for MML/EAP scores with an explanatory
  status.

The reported `lz` statistic is asymptotically standard normal under
the conditional-independence assumption; |lz| > 1.96 / 2.58 are the
5% / 1% reporting flags.

### Generalizability-theory adapter

`mfrm_generalizability(fit)` re-fits the rating data as a crossed
random-effects model `Score ~ 1 + (1 | Person) + (1 | Facet1) + ...`
via `lme4::lmer` and returns the canonical G / Phi coefficients
plus per-source variance components. Useful when a reviewer asks
for a generalizability-theory complement to the Rasch-style
separation / reliability statistics that `diagnose_mfrm()`
already emits.

### Import adapters: mirt / TAM / eRm

Three thin importers expose external fit objects via the same
`mfrm_fit` interface that the mfrmr plot and table helpers
consume:

- `import_mirt_fit(fit, model)` accepts a `mirt::mirt()` result.
- `import_tam_fit(fit, model)` accepts `TAM::tam.mml()` /
  `TAM::tam.jml()`.
- `import_erm_fit(fit, model)` accepts `eRm::PCM()` /
  `eRm::RM()` / `eRm::RSM()`.

The imported objects carry the `mfrm_imported_fit` class and
populate measurement-side slots (`facets$person`,
`facets$others`, `steps`, `summary`) only. Bias / DIF / anchor /
replay slots are explicitly not populated; full bundle import is
planned for a future release.

### Parallel parametric-bootstrap ICC

`compute_facet_icc(boot = "boot")` gains `ci_boot_parallel`
(`"no"` / `"multicore"` / `"snow"`) and `ci_boot_ncpus` arguments
that are forwarded to `lme4::bootMer()`. The per-replicate `cli`
progress bar is suppressed under parallel execution because
worker processes hold their own copy of the progress state.

### Parallel evaluate_mfrm_design (scaffold)

`evaluate_mfrm_design()` accepts a `parallel = c("no", "future")`
argument. When `"future"` is requested and the `future.apply`
Suggests package is installed, the rep loop within each design
row honours whatever `future::plan()` is currently active;
cross-design-row parallelism is planned for a future release. Without
`future.apply` the call falls back to serial execution with an
explicit message.

### Resumable MML EM fits

`fit_mfrm()` accepts a `checkpoint = list(file = ..., every_iter = ...)`
argument. When supplied to a `mml_engine = "em"` (or hybrid)
fit, the EM scaffolding writes its state to `file` every
`every_iter` outer iterations using `saveRDS()`. If the file
exists when a subsequent call starts, the engine resumes from the
recorded iteration. The direct `optim()` engine ignores the
checkpoint; non-EM fits run unaffected.

### GPCM verification tests

A new `tests/testthat/test-gpcm-verification.R` exercises every
`"supported"` and `"supported_with_caveat"` row of
`gpcm_capability_matrix()` on a toy dataset and asserts the
documented helper returns the expected shape. `"blocked"` and
`"deferred"` rows have negative tests that confirm the helper
either refuses to run or returns an explicit caveat. These tests
make the GPCM scope a contract that future commits cannot
silently shrink.

### Optional FACETS Table 7 style fit output on fit$facets$others

`fit_mfrm(attach_diagnostics = TRUE)` runs `diagnose_mfrm()` once
after the fit and merges the per-level `SE`, `Infit`, `Outfit`, and
`PtMeaCorr` columns onto `fit$facets$others`. This makes the facet
table look like a FACETS Table 7 summary without a separate call.
The default `FALSE` preserves the minimal `Facet` / `Level` /
`Estimate` layout from 0.1.5.

## Reproducibility

`build_mfrm_manifest()` gains several new tables so replay bundles
carry everything a deterministic re-run needs:

- `environment` now records `RNGKind`, `RNGSeedDigest`, `Locale`, and
  a UTC ISO-8601 timestamp in addition to the existing package and
  platform fields.
- `dependencies` (new) records the installed version of every
  `Imports` and `Suggests` dependency, with a `Role` column.
- `input_hash` (new) hashes the input data, anchors, group anchors,
  and `score_map` with SHA-256 (via `digest`, now in `Suggests`) or
  an MD5-of-RDS fallback. The hash is deterministic across sessions.
- `session_info` (new) unrolls `utils::sessionInfo()` into a long
  data frame (`Scope` / `Package` / `Version`).
- `hierarchical_audit`, `missing_recoding`, and `shrinkage_audit`
  (new) surface the three new audit layers in one place.

`digest` is added to `Suggests`.

## Bug fixes

- **Bias / interaction NA.** `estimate_bias()` and
  `estimate_all_bias()` previously returned `NA` for every cell's
  `S.E.`, `t`, `Prob.`, `Obs-Exp Average`, `Infit`, and `Outfit`,
  and `Significant` counts collapsed to zero. Root cause was an
  `nzchar(NA_character_)` call (which returns `TRUE`) in an internal
  predicate. Downstream helpers such as `bias_interaction_report()`
  and `plot_bias_interaction()` are now populated again.
- **`estimate_bias()` silent failure on typo'd facet names.** A
  mis-spelled `facet_a` / `facet_b` (e.g. `"Raters"` with trailing
  `s`) previously returned an empty `list()` with no warning. It now
  raises an informative error naming the available facets. Missing
  `diagnostics` argument likewise raises an explicit mfrmr error
  rather than falling through to R's locale-dependent missing-
  argument message.
- **ZSTD sign.** `zstd_from_mnsq()` was numerically unstable for very
  small degrees of freedom and could return large positive ZSTD when
  `MnSq` was close to zero, flipping the sign relative to the
  companion Outfit ZSTD for the same element. A `df >= 1` guard
  returns `NA` in degenerate cells.
- **Score out of range.** `prepare_mfrm_data()` now stops when any
  observed `Score` falls outside the declared
  `[rating_min, rating_max]` range. Previously negative `score_k`
  values passed through `m[cbind(i, 0)]`, silently dropping those
  rows from the likelihood while `n_obs` kept its original value.
- **Silent facet-name mismatches.** `sanitize_noncenter_facet()`,
  `sanitize_dummy_facets()`, and `build_facet_signs()` now emit a
  warning when supplied facet names are not part of the fitted
  model. Previously typos such as `positive_facets = "rater"`
  (lowercase) or `noncenter_facet = "Raters"` could silently flip
  the sign convention of facet measures.
- **Graphical state hygiene.** `apply_plot_preset()` and
  `.draw_shrinkage_plot()` now restore the user's `par()` on exit,
  per "Writing R Extensions" 2.1. All plot methods that relied on
  `apply_plot_preset()` inherit this automatically.
- **DFF contrast sign flip.** `analyze_dff()` adds a
  `ContrastDirection` column to the residual and refit branches.
  The two methods use opposite sign conventions by design, so the
  new column spells out which interpretation applies.
- **`compute_facet_icc()` singular fit.** Total variance below
  `sqrt(.Machine$double.eps)` is now reported as `ICC = NA` with
  `Interpretation = "Non-identifiable"` instead of a falsely
  meaningful value. The first `lme4` convergence diagnostic surfaces
  as a `message()` rather than being silently suppressed.
- **Extreme-person flag persistence.** `as.data.frame.mfrm_fit()`
  now carries the new `Extreme` column through to ggplot2 / CSV
  pipelines instead of dropping it.
- **`as_kable(format = "pipe")` output.** Previously silently
  returned HTML when `kableExtra` was installed and the `apa_table`
  carried a non-empty `note`. `"pipe"` now consistently returns the
  Markdown table with an appended `Note.` line.
- **`review_mfrm_anchors()` false positives.** Overlap-adequacy risk
  flags are skipped when no anchors or group anchors were supplied,
  so single-wave analyses no longer emit "high severity" warnings
  because `OverlapLevels == 0` everywhere.
- **Fractional-score tolerance.** Tightened from
  `sqrt(.Machine$double.eps)` (~`1.5e-8`) to `1e-6`, so integer
  codes like `1.0000001` that round-trip through CSV floats are now
  accepted. Genuinely fractional scores (`1.5`, `2.75`) are still
  caught.
- **Rating-range inference output.** When `rating_min` / `rating_max`
  are inferred from the observed scores, the provenance is now retained
  in fit summaries and data-description output rather than emitted as a
  routine message. Users who prefer the interactive reminder can set
  `options(mfrmr.show_inferred_rating_range = TRUE)`; `fit_mfrm()` still
  limits that opt-in message to one per fit.
- **Locale-independent error for `plot(fit, type = ...)`.** Passing an
  unknown `type` previously raised R's locale-dependent
  `match.arg()` error. It now raises an English mfrmr-style error
  listing the valid choices.
- **`plot_dif_heatmap(draw = FALSE)` return contract.** The helper
  documented an `mfrm_plot_data` object but invisibly returned the
  bare `matrix`, breaking the documented contract used by sibling
  `plot_*` helpers. It now returns an `mfrm_plot_data` whose `data`
  slot bundles `matrix`, `pairs`, `metric`, and `value_column`. Code
  that relied on the old shape should switch from `dim(heat)` to
  `dim(heat$data$matrix)`.
- **Approximate 95% CI whiskers on bias and displacement plots.**
  `plot_bias_interaction(show_ci = TRUE, ci_level = 0.95)` (scatter
  and ranked modes) now draws `BiasSize \u00b1 z \u00b7 SE` whiskers
  using the per-cell SE from [estimate_bias()]. `plot_displacement(
  show_ci = TRUE)` (lollipop mode) draws
  `Displacement \u00b1 z \u00b7 DisplacementSE` whiskers from the
  audit-table standard error. Both functions now populate
  `CI_Lower` / `CI_Upper` / `CI_Level` columns on the returned
  plot-data element so downstream pipelines can reuse the bounds.
  `plot_fair_average()` CI support (now implemented later in this
  release) uses a delta-method propagation because the
  fair-average SE lives on the logit scale while the plot uses the
  observed-score scale, which requires a delta-method transformation.

## Messaging improvements

- `fit_mfrm()` emits a one-time `message()` when called with
  `anchor_policy = "silent"` while the anchor review flags issues.
- `prepare_mfrm_data()` records whether `rating_min` / `rating_max`
  were inferred from the observed scores or supplied explicitly, and
  fit/data summaries surface that provenance. The former informational
  message is opt-in through `options(mfrmr.show_inferred_rating_range =
  TRUE)`. Row drops, ID trimming, and facets with only one observed
  level are recorded in `preparation_notes`; routine preparation
  messages are opt-in through `options(mfrmr.show_preparation_messages =
  TRUE)`.
- Non-numeric score labels (`"low"`, `"medium"`, `"high"`) now raise
  a targeted error up front instead of surfacing as the opaque
  "No valid observations remain" message.
- `detect_anchor_drift()` and `build_equating_chain()` thin-linking
  warnings now list per-facet retained-vs-threshold counts
  (e.g. `"Rater (3/5)"`).
- `bias_interaction_report()$summary` carries a `FlagStatus` column
  so empty `ranked_table` rows are no longer ambiguous between
  "nothing flagged" and "nothing computed".
- Latent-regression fits warn when the design matrix is
  near-singular (`rcond(mm) < 1e-8`), catching numerically collinear
  covariates rather than only exact rank deficiency.

## Documentation and citations

- `apply_empirical_bayes_shrinkage()` docstring and R comments now
  document the shrinkage-variance formula as
  `max(0, K^{-1} * sum(delta^2) - mean(SE^2))`, matching the
  implementation under sum-to-zero identification.
- `?fit_mfrm` "Input requirements" now states the MFRM conditional
  independence assumption (Linacre, 1989) and points at
  `diagnose_mfrm(..., diagnostic_mode = "both")` /
  `strict_pairwise_local_dependence` as the exploratory follow-up.
- `?fit_mfrm` example block now flags the `quad_points = 7` opening
  fit as an exploratory speed setting and reminds readers that the
  package default `quad_points = 31` is the publication tier, so the
  example no longer reads as a recommendation against the new default.
- `?fit_mfrm` now presents the recommended `quad_points` tiers as a
  `\tabular{}` block (`7` fast scan, `15` default, `31+` publication)
  so readers do not have to re-extract the recommendation from prose.
  The "adapted from Linacre (1994)" wording for the sample-size bands
  and the wall-clock cost of `diagnostic_mode = "both"` are also
  spelled out.
- `?fit_mfrm` `missing_codes` docstring is now an itemised list of
  the three branches (`NULL` / `TRUE` / custom vector) instead of
  dense prose.
- `?run_mfrm_facets` now notes that `method = "JML"` is the default
  for legacy FACETS-style output continuity and points users at
  `fit_mfrm(..., method = "MML")` for new analysis scripts.
- `?apply_empirical_bayes_shrinkage` now states that
  `EffectiveDF = Σ(1 − B_j)` matches the "effective number of
  parameters" from Efron & Morris (1973).
- `?compute_facet_icc` notes that Koo & Li (2016) recommend applying
  the reliability bands to the 95% confidence interval of the ICC,
  while the current implementation bands the point estimate only.
- `?analyze_facet_equivalence` adds Kass & Raftery (1995) as the
  reference for the BIC-based Bayes-factor approximation
  `BF_{01} ≈ exp((BIC_{H1} − BIC_{H0}) / 2)`.
- ZSTD docstring in `?mfrmr-package` now cites Wilson & Hilferty
  (1931) explicitly alongside Wright & Linacre (1994).
- `calc_displacement_table()` now cites the Winsteps user guide for
  the combined `|Displacement| > 0.5 logit` and `|t| > 2` flagging
  rule.
- `analyze_facet_equivalence()` docstring frames the
  `equivalence_bound = 0.5` default as a starting point.
- Added `print()` S3 methods for 13 classes that previously fell
  back to the default list printer (`mfrm_apa_outputs`, `mfrm_bias`,
  `mfrm_bundle`, `mfrm_design_evaluation`, `mfrm_diagnostics`,
  `mfrm_facet_dashboard`, `mfrm_future_branch_active_branch`,
  `mfrm_plausible_values`, `mfrm_population_prediction`,
  `mfrm_reporting_checklist`, `mfrm_signal_detection`,
  `mfrm_threshold_profiles`, `mfrm_unit_prediction`). Each delegates
  to the existing `summary()` method.

Reference citations corrected:

- Efron & Morris (1973) page range `379-402` (was `379-421`).
- McEwen (2018) BYU dissertation year (was 2017).
- Wright (1998) RMT 12(2) on extreme scores (was Wright 1988,
  which does not exist).
- Jones & Wind (2018) JAM 19(2), 148-161 (was "Wind & Jones, 2018,
  JAM 19(1), 1-19", which does not exist).
- Linacre (2023) *A User's Guide to Facets, Version 4.5* applied
  uniformly where comments had used 2024.

## Plot polish

- `plot_qc_dashboard()` plots a signed ZSTD distribution combining
  Infit and Outfit ZSTD, with reference lines at
  `-3 / -2 / 0 / 2 / 3`. The previous absolute-value histogram
  collapsed over-fit and under-fit tails.
- `plot_residual_pca(..., plot_type = "scree")` draws Rasch
  secondary-dimension reference lines at `1.0 / 1.4 / 2.0 / 3.0`,
  consistent with the Winsteps user guide. The legend and returned
  `reference_lines` record the new entries.
- The Residual-PCA `content_checks` entry uses a case-insensitive
  regex so the check passes whether the APA contract uses
  `"Residual PCA"` or the longer `"Exploratory residual PCA"`.

## Other additions

- `fit$facets$person` exposes `PosteriorSD` and `SE` aliases
  alongside the legacy `SD`. MML fits populate all three with the
  posterior SD under the Gauss-Hermite prior; JML fits set
  `SE = NA_real_` and note that per-person SEs should be pulled
  from `diagnose_mfrm()$measures`.
- `analyze_dff(method = "refit")` subgroup fits return a
  `LinkingReview` column that captures the anchor-review messages
  emitted during the refit, replacing the previous
  hard-coded `anchor_policy = "silent"` silence.
- `detect_anchor_drift()` returns `common_vs_reference` and
  `n_common_all_waves` alongside the existing pairwise
  `common_elements` table, for 3+ wave linking reviews.
- `analyze_residual_pca(..., pca_max_factors = "auto")` caps the
  factor count at `min(10, ncol - 1, nrow - 1)` per matrix; this
  value was previously silently coerced to `NA`.
- `describe_mfrm_data()` returns two new components:
  `missing_rate_summary` (per-column missing / non-missing counts)
  and `facet_crosstabs` (long-format pairwise observation-count
  tables, suitable for heatmap plotting).
- DESCRIPTION removes the duplicated `Author` / `Maintainer` fields
  auto-generated by CRAN; `Authors@R` remains the single source of
  truth. The `Description:` field is now three sentences (was 10
  lines of prose), improving CRAN web readability while retaining
  the two rating-scale / partial-credit DOI references.
- `inst/CITATION` now tracks `meta$Version` and `meta$Title`
  dynamically, so `citation("mfrmr")` prints the current installed
  version rather than a hard-coded string.

## Test suite

6,380+ tests pass (up from 6,343 in 0.1.5), with 0 failures and
0 errors. New test files:

- `test-shrinkage.R` (40 tests) covers the closed-form math, edge
  cases (`K < 3`, `tau^2 <= 0`, user-supplied prior), `fit_mfrm`
  integration, reporting and manifest trails, and the three new
  plot methods.
- `test-missing-codes-integration.R` (17 tests) covers `fit_mfrm`,
  `describe_mfrm_data`, `review_mfrm_anchors`, and manifest paths.
- `test-hierarchical-audit.R` (10 tests) covers the five new
  hierarchical-audit helpers and their integration points.

Pre-existing test-harness errors unrelated to 0.1.5 behaviour have
also been cleaned up (S3 dispatch, GPCM scope wording, internal-helper
prefixing with `mfrmr:::`).

# mfrmr 0.1.5

## Maintenance release

### First-use workflow

- Reworked `print(fit)`, `summary(fit)`, and
  `summary(diagnose_mfrm(...))` so results start with `Status`,
  `Key warnings`, and `Next actions`.
- Added a clearer recommended workflow in the README and help pages: fit with
  `MML`, review diagnostics with `diagnostic_mode = "both"`, then move to
  reporting helpers.
- Improved ordered-score handling and guidance, including binary
  two-category use, rejection of fractional score values, non-consecutive
  score-code mapping through `score_map`, and clearer warnings for retained
  zero-count categories.

### Estimation and scoring

- Added the first public latent-regression `MML` branch for ordered `RSM` /
  `PCM` fits with person covariates, including simulation and scoring support
  for the fitted population model.
- Added bounded `GPCM` support for the documented direct route, including core
  summaries, diagnostics, plots, posterior scoring, and information checks,
  while keeping unsupported downstream routes explicit.
- Extended ordered-response support and documentation for binary `RSM` / `PCM`
  use, fixed-calibration scoring after `JML`, and `PCM` information curves.

### Diagnostics, reporting, and visualization

- Added strict marginal follow-up plots through `plot_marginal_fit()` and
  `plot_marginal_pairwise()`.
- Strengthened the reporting surface with `reporting_checklist()`,
  `build_summary_table_bundle()`, `export_summary_appendix()`, and
  `visual_reporting_template()` for manuscript-oriented tables, appendix
  artifacts, and figure-placement guidance.
- Added structured caveats in summaries and appendix tables for retained
  zero-count score categories and latent-regression population-model
  omission/design issues.
- Added exploratory `plot(fit, type = "ccc_surface", draw = FALSE)` output
  for advanced visualization while keeping 2D Wright/pathway/category plots as
  the default reporting route.

### External-software scope

- Added scoped ConQuest overlap helpers and concise software-scope summaries
  for FACETS, ConQuest, and SPSS handoffs.
- Clarified latent-regression reporting outputs so coefficient reporting is
  kept separate from post hoc score regression.

# mfrmr 0.1.4

## CRAN resubmission

- Replaced a misencoded author name in documentation references so the PDF
  manual builds cleanly under CRAN's LaTeX checks.
- Revised DESCRIPTION references again to avoid incoming spell-check notes
  while preserving the requested author-year-doi citation format.

# mfrmr 0.1.3

## CRAN resubmission

- Revised `DESCRIPTION` references to use the requested `authors (year) <doi:...>`
  format.
- Added a documented return-value section for `facet_quality_dashboard()`,
  including the output class, structure, and interpretation.
- Replaced `\dontrun{}` with `\donttest{}` for executable examples so CRAN can
  exercise those examples during checks.

# mfrmr 0.1.2

## CRAN resubmission

- Further reduced CRAN check time by trimming the CRAN-only test subset to
  lightweight smoke tests after the incoming pretest still reported a Windows
  overall-checktime NOTE for version 0.1.1.

# mfrmr 0.1.1

## CRAN resubmission

- Revised `DESCRIPTION` metadata to avoid CRAN incoming spell-check notes on
  cited proper names.
- Reduced CRAN check time by skipping long integration and coverage-expansion
  test files during CRAN checks while keeping the full local test suite.

# mfrmr 0.1.0

## Initial release

- Native R implementation of many-facet Rasch model (MFRM) estimation without TAM/sirt backends.
- Supports arbitrary facet counts with `fit_mfrm()` and method selection (`MML` default, `JML`).
- Includes FACETS-style bias/interaction iterative estimation via `estimate_bias()`.
- Provides fixed-width report helpers (`build_fixed_reports()`).
- Adds APA-style narrative output (`build_apa_outputs()`).
- Adds visual warning summaries (`build_visual_summaries()`) with configurable threshold profiles.
- Implements residual PCA diagnostics and visualization (`analyze_residual_pca()`, `plot_residual_pca()`).
- Bundles Eckes & Jin (2021)-inspired synthetic Study 1/2 datasets in both `data/` and `inst/extdata/`.

## Package operations and publication readiness

- Added GitHub Actions CI for cross-platform `R CMD check`.
- Added `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, and `SECURITY.md`.
- Added citation metadata (`inst/CITATION`, `CITATION.cff`).
- Expanded README with explicit installation and citation instructions.

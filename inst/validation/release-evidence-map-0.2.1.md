# mfrmr 0.2.1 evidence map

This note tracks the 0.2.1 development focus after the 0.2.0 release:
bounded `GPCM` refinement, slope-aware recovery review, and uncertainty
interpretation.

## Literature basis checked

- Rasch (1960/1980), Probabilistic Models for Some Intelligence and
  Attainment Tests, supports the dichotomous person-item model as the
  two-category reference case for Rasch-family measurement.
- Andrich (1978), "A Rating Formulation for Ordered Response Categories",
  Psychometrika, 43, 561-573, doi:10.1007/BF02293814. This supports the
  rating-scale ordered-response parameterization; with one category boundary,
  the adjacent-category logit is the ordinary binary Rasch logit up to
  centering and threshold-identification conventions.
- Muraki (1992), "A Generalized Partial Credit Model: Application of an EM
  Algorithm", Applied Psychological Measurement, 16(2), 159-176,
  doi:10.1177/014662169201600206. This is the direct source for the GPCM
  model basis used by the bounded `GPCM` route.
- Muraki (1993), "Information Functions of the Generalized Partial Credit
  Model", Applied Psychological Measurement, 17(4), 351-363,
  doi:10.1177/014662169301700403. This supports the package's slope-aware
  information-function interpretation for polytomous GPCM-style items.
- Morris, White, and Crowther (2019), "Using simulation studies to evaluate
  statistical methods", Statistics in Medicine, 38(11), 2074-2102,
  doi:10.1002/sim.8086. This supports keeping aims, data-generating
  mechanisms, estimands, methods, and performance measures explicit in
  recovery validation outputs.
- Wind, Jones, and Grajeda (2023), "Does Sparseness Matter? Examining the
  Use of Generalizability Theory and Many-Facet Rasch Measurement in Sparse
  Rating Designs", Applied Psychological Measurement, 47(5-6), 351-364,
  doi:10.1177/01466216231182148. This supports treating sparse rating
  designs as planned missingness designs where rater commonality and linkage
  need to be inspected explicitly.
- Wind and Jones (2018), "The Stabilizing Influences of Linking Set Size and
  Model-Data Fit in Sparse Rater-Mediated Assessment Networks", Educational
  and Psychological Measurement, doi:10.1177/0013164417703733. This supports
  exposing linking-set size and rater-pair common-person counts as design
  diagnostics rather than fit statistics.
- DeMars, Shapovalov, and Hathcoat (2023), "Many-Facet Rasch Designs: How
  Should Raters be Assigned to Examinees?", NCME presentation. This supports
  the pragmatic common-linking-set option when many examinees receive only a
  small number of ratings.
- Farrokhi, Esfandiari, and Schaefer (2012), "A Many-Facet Rasch Measurement
  of Differential Rater Severity/Leniency in Three Types of Assessment",
  JALT Journal, 34(1), 79-102, doi:10.37546/JALTJJ34.1-3. This supports
  treating peer-assessor severity/leniency as an MFRM reporting context, while
  leaving the assignment-design choices to the package helper and design
  diagnostics.
- Uto and Ueno (2020), "A generalized many-facet Rasch model and its Bayesian
  estimation using Hamiltonian Monte Carlo", Behaviormetrika, 47, 469-496,
  doi:10.1007/s41237-020-00115-7. This supports the broader peer-assessment
  and generalized-MFRM context. The 0.2.1 helper does not implement the
  Bayesian Uto/Ueno model; it only supplies a fixed peer-review design
  generator for the package's current ordered-response simulation route.
- Wright and Linacre (1994), "Reasonable mean-square fit values", Rasch
  Measurement Transactions, 8(3), 370, and Linacre (2002), "What do Infit and
  Outfit, Mean-square and Standardized mean?", Rasch Measurement Transactions,
  16(2), 878. These support the package distinction between mean-square fit
  size and standardized ZSTD interpretation.
- Wright and Masters (1982), Rating Scale Analysis, and Wright and Masters
  (2002), "Number of Person or Item Strata: (4*Separation + 1)/3", Rasch
  Measurement Transactions, 16(3), 888. These support the G/R/H separation,
  reliability, and strata convention used in `diagnostics$reliability`.
- The FACETS user guide's WHEXACT / fit-standardization documentation supports
  treating ZSTD differences as convention-sensitive when mean-square fit agrees
  but df or Wilson-Hilferty-style standardization differs.

Zotero was checked locally for these sources. The local library contains
Muraki (1992), Muraki (1993), and Morris et al. (2019). Muraki (1993) cites
Samejima (1974), but the current 0.2.1 bounded-`GPCM` changes do not add a
Samejima normal-ogive or graded-response-model implementation.

## Fit and separation reporting boundary

Fit and separation are source-grounded reporting diagnostics, but they are not
single-number validation gates.

- Mean-square fit (`Infit`, `Outfit`) is the primary fit-size diagnostic.
  Wright and Linacre (1994) and Linacre (2002) support reading values relative
  to the expected value of 1 and treating broad bands as practical review
  conventions rather than universal validity criteria.
- ZSTD is a standardized version of mean-square fit. Linacre (2002) and the
  FACETS WHEXACT documentation support the 0.2.x design choice to compare MnSq
  first and to label df-driven ZSTD changes as convention-sensitive through
  `fit_df_method = "both"` and `facets_fit_review()`.
- Separation, reliability, and strata follow the Wright/Masters G/R/H
  convention: spread relative to average measurement error, transformed to
  separation reliability and strata. These are not inter-rater agreement and
  not proof of substantive validity by themselves.
- Package QC thresholds in `run_qc_pipeline()` and design-simulation summaries
  are policy overlays on top of the literature-grounded formulas. They belong
  in operational triage and design planning, not in the source-validation layer
  as if they were direct literature cut points.
- Recovery-validation fit/separation summaries are therefore retained as
  diagnostic operating characteristics. They can show how MnSq, ZSTD,
  separation, reliability, and strata behave under a simulated condition, but
  they do not enter the top-line release-recovery status.
- Recovery-assessment and recovery-validation summaries also expose
  `diagnostic_reporting_notes`, a reporter-facing table for zero
  separation/reliability, abs-ZSTD flags, and df-sensitive ZSTD flags. These
  notes are caveat prompts, not psychometric adequacy decisions.

0.2.1 therefore adds `precision_review_report()$fit_separation_basis`, a
compact table that states each topic's source basis, package surface,
interpretation, validation use, and availability in the current run.

## Peer-review simulation boundary

Peer-review or peer-assessment simulations have a stronger design constraint
than ordinary sparse rater-mediated designs: the submission pool and reviewer
pool can be the same people. The 0.2.1 helper therefore uses a fixed
person-by-reviewer-by-criterion skeleton with shared IDs and optional
structural self-review exclusion.

- Peer-assessor severity/leniency is a legitimate MFRM use case, following
  Farrokhi et al. (2012), but the simulation helper does not claim that peer
  ratings are automatically fair, interchangeable, or substantively valid.
- Common-link anchor submissions follow the design logic emphasized by DeMars
  et al. (2023) for sparse or limited-rating settings. They are not universal
  thresholds for fit, separation, rater quality, or recovery.
- Uto and Ueno (2020) remain a source for generalized MFRM and Bayesian
  peer-assessment modeling context. The current package route stays within the
  existing `RSM` / `PCM` / bounded-`GPCM` generator and does not add HMC,
  Bayesian priors, Markov rater drift, or multidimensional peer-assessment
  estimation.
- Peer-review metadata (`SelfReviews`, reviewer load, reciprocal review
  pairs, and common submissions per reviewer pair) is therefore reported as
  design evidence. It can be shown through `build_mfrm_network_review()`, but
  it stays separate from MFRM estimates, fit, separation, and recovery gates.

## Implementation boundary

The `slope_regime` field (`unit_slopes`, `near_flat`, `moderate`,
`high_dispersion`) is not a psychometric adequacy threshold from the cited
literature. It is a package validation label that summarizes centered
log-slope spread so recovery simulations can be read against the intended
generator stress level.

The current operational cut points are intentionally conservative engineering
bins:

- `unit_slopes`: all identified slopes are effectively 1.
- `near_flat`: maximum absolute centered log slope is no larger than
  `log(1.05)`.
- `moderate`: maximum absolute centered log slope is no larger than
  `log(1.50)`.
- `high_dispersion`: centered log-slope spread exceeds the `moderate` bin.

These labels should be reported as generator-condition metadata, not as
model-choice evidence and not as release-level recovery success or failure by
themselves.

`mfrm_results()` is also an implementation-layer UX wrapper, not a new
estimator, diagnostic standard, or validation argument. Its purpose is to make
the existing `fit_mfrm()` -> `diagnose_mfrm()` -> table/review/report/plot
surfaces easier to inspect for users who expect a comprehensive first screen.
Section-level failures are reported as `not_available` so unsupported routes
remain visible without being confused with psychometric evidence. Direct
`data.frame` input is intentionally conservative and should be replaced by an
explicit `fit_mfrm()` call in reproducible analysis scripts whenever column
roles, model settings, anchors, or missing-data rules require documentation.
`mfrm_results_interactive()` is opt-in only; it exists for guided column
selection and code generation in interactive sessions, not for package checks,
batch scripts, or manuscript replay.
`mfrmr_output_guide("entry")` is the companion public-API map for this route;
it documents where to start and where to branch next, but it does not add a
new analysis method. Its lifecycle, user-level, and recommended-entry labels
are navigation metadata, not validation grades. The related
`mfrmr_output_guide("viewer")` maps optional local-viewer workflows back to
the `mfrm_results(include = ...)` object that must be created first. It is a
navigation table for inspection tasks, not a new reporting standard or
interactive analysis record. Likewise,
`mfrm_results(include = "publication"|"validation"|"facets"|"network"|
"gpcm_review")` expands to existing sections; the presets do not add new
estimators, diagnostics, or acceptance thresholds.
`include = "bias"` and `include = "misfit_review"` follow the same pattern:
they expose existing bias-screen, unexpected-response, displacement, and
pathway-map review surfaces without selecting fairness contrasts or turning
case-review prompts into validity decisions.
`include = "linking"` follows the same wrapper rule. It exposes the fitted
object's stored anchor-review evidence and the existing `build_linking_review()`
operational surface when available. It does not infer drift from one fit and
does not replace the explicit multi-fit `detect_anchor_drift()` or
`build_equating_chain()` workflows needed for wave/form comparison. Bounded
`GPCM` linking synthesis remains caveated by the existing capability matrix.
The `summary(mfrm_results(...))$triage` table is likewise a reading-order
surface over existing diagnostics, status, plot, table, precision, reporting,
model-scope, and network-review availability. Its `Severity` labels route
attention; they are not new psychometric cut points or validation decisions.
`launch_mfrmr_viewer()` follows the same boundary. It is an optional Shiny
reader for an already-created `mfrm_results` object, not a fitting wizard,
external-web-app bridge, or new validation layer. The viewer's Replay tab is
intended to keep GUI inspection subordinate to the explicit
`fit_mfrm()` -> `mfrm_results()` workflow.
The QC, report, bias, and pathway/misfit tabs remain display surfaces over
existing components. In particular, bias-interaction review still requires an
explicit facet-pair choice; the viewer does not select contrasts or promote
screening signals into fairness conclusions. Tab-local section-status tables
are navigation aids over `summary(res)$status`, `plot_map`, and collected
tables; they make unavailable and not-requested sections visible, but they are
not additional diagnostics.
`export_mfrm_results()` follows the same object-first boundary. It writes the
already-created `mfrm_results` summary tables, collected tables, HTML, RDS,
replay scaffold, and manifest to disk. When `include = "report"` is requested,
it also writes `mfrm_report()` CSV, Markdown, and HTML artifacts from the same
stored result object. It does not recompute diagnostics, change section
availability, or create a new validation layer. Optional plot export is a
graphics handoff over existing `plot(res, type = ...)` routes.
`mfrm_report()` follows the same boundary for report drafting. It converts an
existing `mfrm_results` object into QC, APA, validation, reviewer, or technical
section plans with claim-readiness, report-gap, evidence-boundary, and
next-action tables. Its `report_index` table is the compact entry point that
lists each major evidence area, evidence status, readiness label, review-signal
count, and the primary/template tables to inspect next. Its fit-specific
`fit_criteria`, `zstd_conventions`, and `fit_decision_policy` tables expose
multiple MnSq threshold profiles and the engine-vs-FACETS-style df/ZSTD
convention boundary before report wording.
Its stored-result evidence tables (`fit_evidence_summary`,
`fit_threshold_sensitivity`, `fit_df_sensitivity_summary`, and
`fit_df_sensitive_rows`) summarize observed fit-status counts and
df-sensitive ZSTD prompts from `res$components$fit_measures`.
`fit_reporting_templates` then converts those stored counts into cautious
APA/QC/validation/reviewer wording scaffolds while preserving separate
sentences for MnSq status, ZSTD standardization, df sensitivity, and
separation/reliability. `precision_evidence_summary`, `precision_basis`, and
`precision_reporting_templates` extend the same boundary to separation,
reliability, and strata, using the Wright/Masters G/R/H convention and the
stored precision review. These rows explicitly prevent Rasch/FACETS-style
separation reliability from being reported as inter-rater agreement, model
fit, or standalone validity evidence. `bias_evidence_summary` and
`bias_reporting_templates` extend the same boundary to bias, DFF, and fairness
wording when the result was explicitly built with `include = "bias"`. They
keep facet-level screens, interaction-bias contrasts, group-by-facet DFF
claims, and fairness conclusions separate. `misfit_evidence_summary` and
`misfit_reporting_templates` extend the same boundary to local misfit review
when the result was explicitly built with `include = "misfit_review"`. They
keep unexpected-response rows, displacement review, pathway-map evidence, and
case-review wording separate. `linking_evidence_summary` and
`linking_reporting_templates` extend the same boundary to anchor readiness,
drift review, and screened equating-chain review when the result was explicitly
built with `include = "linking"`. They keep single-fit anchor-readiness
support separate from multi-fit drift/equating claims. It does not refit the
model, recompute diagnostics, select fairness contrasts, infer drift, or turn
fit, separation, bias-screen, pathway, local misfit, or anchor evidence into
automatic acceptance or exclusion decisions.

Binary person-item data remain inside the same ordered-response contract. When
`person` names the person column, `facets = "Item"` names the single
non-person item facet, and the observed score support has exactly two ordered
integer categories, `model = "RSM"` is the ordinary dichotomous Rasch route.
The package does not need a separate exported `rasch()` front end for this
case; instead, `mfrmr_output_guide("binary")` documents the first-screen route
and checks. `PCM` is still computationally available for two-category data,
but for ordinary person-item binary tests the RSM route is the cleaner default
unless item-specific step structure is part of the intended design.

The bounded-`GPCM` capability matrix is now also the runtime guard contract.
`gpcm_runtime_guard_coverage()` records how each blocked or deferred
capability row is enforced: public helpers either stop with
`mfrmr_gpcm_scope_error` or remain explicitly marked as roadmap-only when no
public runtime surface exists yet. The structured error carries the matrix
area, status, recommended route, and next validation step, so wrappers can
route unsupported paths without parsing prose. The release-readiness protocol
checks that this coverage table stays synchronized with
`gpcm_capability_matrix()` and that every outstanding row is represented.
`mfrmr_output_guide("gpcm")` is the short user-facing route map to both
helpers; it is navigation metadata and does not broaden the bounded-`GPCM`
support contract.

## 0.2.1 checks added

- Recovery assessment now separates unavailable SE/coverage evidence from
  intentionally unassessed coverage thresholds in `uncertainty_review`.
- Bounded-`GPCM` simulation specifications now carry `slope_regime`, and
  `evaluate_mfrm_recovery()` carries it into `settings` and ADEMP method
  metadata.
- Recovery assessment now carries `condition_review`, a generator-condition
  table that reports bounded-`GPCM` slope regime, stress level, slope count,
  centered log-slope spread, and generated score-category support before
  metric and uncertainty evidence are interpreted.
- Recovery assessment and validation summaries now also expose
  `condition_reporting_notes`, a reporter-facing table that separates
  high-dispersion slope stress and sparse generated score support from
  recovery-metric failures.
- Recovery simulation replication overviews now retain generated score-support
  metadata so sparse-category stress can be separated from slope-regime stress
  when reading high-dispersion bounded-`GPCM` runs.
- `evaluate_mfrm_recovery(include_diagnostics = TRUE)` now retains
  replication-by-facet and facet-level fit/separation operating
  characteristics. `assess_mfrm_recovery()` carries these into
  `diagnostic_review` with `ValidationUse =
  "diagnostic_only_not_release_gate"`.
- `assess_mfrm_recovery()` also carries diagnostic operating characteristics
  into `diagnostic_reporting_notes`, a reporter-facing table that should be
  read before the raw diagnostic review when drafting fit, separation, or
  reliability caveats.
- The recovery-validation protocol includes an extended
  `gpcm_high_dispersion_sparse` case. It is not part of the default core tier;
  it is a targeted sensitivity fixture for checking that `condition_review`
  and `condition_summary` distinguish slope-regime stress from sparse
  generated score support.
- Top-line validation output now separates core release status from extended
  sensitivity status. This keeps high-dispersion/sparse-category stress
  evidence visible without treating it as a core release blocker by default.
- Recovery-validation summaries now include condition reporting notes,
  diagnostic operating characteristics, `diagnostic_reporting_notes`, and
  `DiagnosticStatus`, while `PrimaryDecisionBasis` remains recovery metrics,
  convergence, and Monte Carlo precision. `DiagnosticStatus` is an
  availability/status-routing field, not a judgement that fit or separation
  values are psychometrically adequate.
- Summary-table bundles now accept `summary.mfrmr_recovery_validation` objects
  and register recovery-assessment and recovery-validation appendix roles for
  the top-line decision, release decisions, case summary, condition summary,
  condition reporting notes, diagnostic reporting notes, diagnostic summary,
  and domain decisions.
- Regression tests cover all slope-regime classes, including the exact
  unit-slope PCM-reduction case.
- GPCM scope regression tests now verify that every blocked or deferred
  capability row appears in `gpcm_runtime_guard_coverage()` and that each
  runtime-guarded helper returns an `mfrmr_gpcm_scope_error` whose area,
  status, recommended route, and next validation step match the capability
  matrix exactly. `mfrmr_output_guide("gpcm")` also has regression coverage so
  users can find both the capability matrix and guard coverage table from the
  public route guide.
- `mfrmr_interval_guide()` now maps public 95% CI and interval-like
  uncertainty routes across fit-measure tables, Wright maps, fair averages,
  bias screens, displacement, DFF/DIF summaries, anchor drift, rater
  severity profiles, rater trajectories, manuscript Figure 1 composites,
  shrinkage, and ICC review. The guide records the interval basis and
  interpretation boundary so confidence displays remain precision, screening,
  or design-review evidence rather than automatic fit, fairness, or validity
  decisions.
- `build_mfrm_sim_spec()` and `simulate_mfrm_data()` now support
  `assignment = "sparse_linked"` for planned-missing sparse rating designs.
  The generated data retain sparse-design metadata for design density,
  planned missingness, rater coverage, and rater-pair common-person links.
- `evaluate_mfrm_design()` and `evaluate_mfrm_recovery()` carry those
  sparse linked generators into repeated simulation workflows. Sparse columns
  remain design diagnostics; they are deliberately separate from fit,
  separation, and parameter-recovery decisions.
- `build_summary_table_bundle()` now gives sparse linked design diagnostics a
  separate `sparse_design` appendix table for both design-evaluation and
  recovery-simulation summaries. This follows the same boundary: design
  density, planned missingness, and rater-pair linkage are reported as design
  evidence, not as recovery or psychometric adequacy gates.
- `sparse_design$LinkReviewStatus` labels zero common-person rater pairs,
  target shortfalls, or unavailable rater-pair common-person counts as
  design-review issues. The status uses Wind and Jones' linking-set framing
  as a reporting boundary; it is not a literature-derived universal cut point
  for fit, separation, or recovery adequacy.
- Design and recovery summaries now add `sparse_review`, a compact count
  surface for the same review statuses. `plot.mfrm_design_evaluation()` can
  draw planned missingness and rater-link metrics, but these plots remain
  design diagnostics and do not change the recovery-validation gate.
- Diagnostic-screening simulations can optionally retain report-layer signals
  with `evaluate_mfrm_diagnostic_screening(include_report = TRUE)`. The added
  `report_signal_summary` rows summarize `mfrm_report()` `report_index`
  availability, readiness, and review-signal counts across scenarios. They are
  report-routing operating characteristics and remain separate from
  diagnostic-screening Type I proxies, sensitivity proxies, and
  recovery-validation decisions.
- `plot.mfrm_diagnostic_screening()` adds an integrated review surface over
  the same simulation output. Overview plots combine legacy ZSTD, strict
  marginal, strict pairwise, strict combined, and optional report-review rates;
  report, contrast, and runtime views are navigation aids for interpreting
  operating characteristics, not additional validation gates. Draw-free calls
  return `mfrm_plot_data`, so `plot_data(..., component = "plot_long")` can be
  used for external visualization or export without changing the evidence
  boundary.
- `build_mfrm_network_review()` now packages the existing co-observation
  graph analysis into a review surface for report handoff. Its overview labels
  disconnected components as warnings and articulation points / bridge edges
  as review issues, while explicitly marking `ReviewUse =
  "design_diagnostic_not_measurement_gate"`.
- Network-review table bundles expose graph-level, facet-level, central-node,
  articulation-node, bridge-edge, optional sparse-review, caveat, and
  reporting-map tables. These tables follow Wind and Jones' sparse
  rater-mediated assessment network framing and the DeMars et al. common-link
  design recommendation: they support design-link interpretation, not
  replacement of MFRM estimates or fit/separation diagnostics.
- `build_peer_review_sim_spec()` now creates fixed peer-review design skeletons
  for peer-assessment or peer-review scoring studies. The generated data carry
  `mfrm_peer_review_design` metadata with self-review counts, reviewer load,
  reciprocal review pairs, design density, and reviewer-pair common-submission
  counts.
- `build_peer_review_design_review()` now converts the same metadata into a
  reportable design-review object with front-door status, load summaries,
  low-common reviewer-pair tables, reciprocal review pairs, caveats, and a
  reporting map. Its statuses route assignment follow-up; they are not peer
  quality, fairness, fit, separation, or recovery decisions.
- `build_mfrm_network_review(peer_review_design = ...)` can carry those
  peer-review design diagnostics into the same appendix-ready design-review
  surface used for graph connectedness and sparse-design metadata. The table
  role is `peer_review_design_diagnostics`, and its `ReviewUse` remains
  `design_diagnostic_not_measurement_gate`.
- `mfrm_results()` now provides a comprehensive first-screen object with
  automatic diagnostics, section-status reporting, collected data frames,
  plot-route metadata, next-action routing, replay-code scaffolding, optional
  temporary HTML rendering, and `build_summary_table_bundle()` support. This
  consolidates UX without replacing the lower-level table/report/analysis/review
  helpers.
- `mfrmr_output_guide("entry")` now exposes the same UX boundary as a compact
  public-API map, separating explicit-fit, comprehensive-results,
  purpose-specific-guide, and opt-in interactive routes.
- `mfrmr_output_guide("viewer")` now records which `mfrm_results(include = ...)`
  object should be built before launching the optional local viewer for
  publication, validation, bias-screen, pathway/misfit, or combined review.
- `mfrmr_output_guide()` now labels lifecycle, user level, and recommended
  entry status, and adds simulation/network route rows as navigation metadata.
- `mfrm_results()` purpose presets now expand common first-screen workflows
  without changing the underlying delegated methods.
- `mfrm_results(include = "bias")` and
  `mfrm_results(include = "misfit_review")` now provide purpose-specific
  first-screen objects for bias-screen and pathway/misfit review. They are
  display and routing presets over existing diagnostics and table helpers, not
  new inference rules.
- `mfrm_results(include = "linking")` now provides the corresponding
  first-screen object for anchor readiness and operational linking review,
  including an anchor plot route. Drift and screened-chain review remain
  explicit multi-fit workflows.
- `summary.mfrm_results` now includes `triage`, and the summary-table bundle
  route exports it as a review-status table.
- `launch_mfrmr_viewer()` now gives users an optional local Shiny reader for
  existing `mfrm_results` objects, including QC, report, bias-screen, and
  pathway/misfit tabs. Regression tests cover its `mfrm_results`-first input
  contract, payload construction, tab-local section-status metadata, and
  optional-dependency guard.
- `export_mfrm_results()` now provides the corresponding lightweight file
  handoff for `mfrm_results` objects, with summary CSVs, collected result
  tables, HTML, RDS, replay code, manifest files, optional plot PNGs, and
  optional zip creation.
- `mfrm_report()` now provides the corresponding report-drafting surface for
  `mfrm_results` objects, with QC, APA, validation, reviewer, and technical
  styles plus claim-readiness, report-gap, evidence-boundary, and next-action
  tables. It also includes fit-specific tables for threshold-profile and
  ZSTD-convention reporting, plus stored-result summaries and wording
  scaffolds for observed fit status, threshold-profile sensitivity, and
  df/ZSTD sensitivity. Precision-specific summaries and wording scaffolds now
  cover separation, reliability, strata, and precision-tier caveats without
  treating them as inter-rater agreement or standalone validity gates.
  Bias-specific summaries and wording scaffolds now cover facet-level bias
  screens and explicit interaction/DFF follow-up without treating screening
  rows as final fairness or invariance conclusions.
- `mfrmr_output_guide("binary")` now records the two-category person-item
  route so users coming from binary Rasch examples can find the correct
  `fit_mfrm()` call, score-support check, and first-screen `mfrm_results()`
  workflow without adding another exported estimation function.
- `release-evidence-checklist-0.2.1.csv` now records the 0.2.1 release
  blocker, caveat, and roadmap checks separately from the historical 0.2.0
  checklist, including the CRAN-time lightweight-test boundary and the
  non-CRAN regression-evidence route.
- `release-readiness.R` now resolves the target release from `DESCRIPTION`,
  selects versioned evidence-map and checklist files when available, and treats
  a stale `R CMD check` log whose package version does not match the target
  release as a package-check review item.
- `release-readiness.R` also checks bounded-`GPCM` runtime guard coverage:
  every blocked/deferred capability row must be covered by a runtime guard or
  roadmap-only entry, runtime guards must use `mfrmr_gpcm_scope_error`, and
  the guard metadata must match the capability matrix.

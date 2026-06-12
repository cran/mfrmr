# Post-0.2.1 bounded GPCM roadmap

This note tracks bounded-`GPCM` work that remains caveated, blocked, or
deferred after the 0.2.1 release boundary. It is a maintenance roadmap, not a
public support promise. The current public contract remains
`gpcm_capability_matrix()`.

## Current release boundary

0.2.1 supports bounded `GPCM` fitting, core summaries, fixed-calibration
posterior scoring, information, curve/category views, direct simulation-spec
generation, direct parameter-recovery checks, summary-table appendix routing,
fair-average review, residual-bias screening, package-native scorefile export,
report/QC bundles, linking synthesis, role-based design forecasting, and
diagnostic/signal-detection design screening within the caveats documented in
the help pages.

The bounded route still requires an explicit step facet and the current
`slope_facet == step_facet` contract. Direct recovery evidence is not design
operating-characteristic evidence, and exploratory diagnostic screens are not
standalone fairness or validity decisions.

Every blocked or deferred capability row is tracked in
`gpcm_runtime_guard_coverage()`. Rows with public helper surfaces must stop
with `mfrmr_gpcm_scope_error`; rows without a public runtime surface are
marked `roadmap_only`. Caveated rows are not guard rows, but this roadmap
still records the evidence needed before their wording can become stronger.
This roadmap should stay aligned with that coverage table and
`gpcm_capability_matrix()`.

## Roadmap work packages

### GPCM score-side export contract

Capability rows:

- `FACETS output-contract score-side review` (`blocked`)
- `Score-side scorefile export under bounded GPCM`
  (`supported_with_caveat`; keep FACETS-equivalence wording out of scope)
- `APA writer and fit-based export bundles` (`supported_with_caveat`; keep
  operational wording constrained by this score-side contract)

Surface to keep blocked until the full contract is validated:

- `facets_output_contract_review()`

Caveated surfaces that must continue to carry `gpcm_boundary` until the
score-side contract is validated:

- `facets_output_file_bundle(include = "score")`
- `build_apa_outputs()`
- `build_visual_summaries()`
- `run_qc_pipeline()`
- `build_mfrm_manifest()`
- `build_mfrm_replay_script()`
- `export_mfrm_bundle()`

Required evidence before unblocking:

- Keep `gpcm_score_side_contract()` synchronized with this roadmap and
  `gpcm_runtime_guard_coverage()`.
- Define the bounded-`GPCM` score-side estimand separately from Rasch-family
  measure-to-score semantics.
- Keep native structural expected-score SEs and selectable score-side delta
  SEs in the scorefile route synchronized with the MML diagnostics contract,
  and separately define the FACETS-compatible uncertainty contract needed for
  full output-contract review.
- Preserve unit-slope `GPCM` reduction tests against the `PCM` route.
- Add negative tests that fail if unsupported score-side rows are silently
  emitted.
- Add release-note wording that keeps sensitivity-model output separate from
  operational scoring claims.

Exit criteria:

- `gpcm_capability_matrix()` can move exactly scoped score-side rows from
  `blocked` to `supported_with_caveat` or `supported`.
- Export bundles identify unsupported sections explicitly when partial support
  remains.

### GPCM design operating characteristics and forecasting

Capability rows:

- `Design evaluation and population forecasting under bounded GPCM`
  (`supported_with_caveat`)
- `Diagnostic and signal-detection design screening under bounded GPCM`
  (`supported_with_caveat`)
- `Differential facet functioning screening under bounded GPCM`
  (`supported_with_caveat`)

Surfaces available as bounded-`GPCM` sensitivity evidence:

- `evaluate_mfrm_design()`
- `predict_mfrm_population()`
- `evaluate_mfrm_diagnostic_screening()`
- `evaluate_mfrm_signal_detection()`
- `analyze_dff()` / `analyze_dif()`
- `dif_interaction_table()`
- `dif_report()`
- `plot_dif_heatmap()` / `plot_dif_summary()`

These helpers are available only for the current role-based design layer and
only when the requested design matches the bounded slope structure carried by
the simulation specification. They report design-level sensitivity evidence
and slope-aware operating-characteristic readouts, not direct posterior
scoring for observed units, calibrated inferential tests, operational scoring
or screening adequacy, or full arbitrary-facet planning.

DFF/DIF helpers are available as direct slope-aware screening and reporting
surfaces. Their `gpcm_boundary` rows must remain visible, and DFF wording must
not imply fairness, invariance, or operational subgroup-decision evidence
without external study-design support.

Required evidence before unblocking:

- Preserve the separation between direct parameter-recovery checks and design
  operating characteristics in ADEMP terms.
- Validate bounded-`GPCM` data-generating conditions across slope regimes,
  sparse linkage patterns, sample sizes, and score-support stress.
- Define which performance measures are release gates and which are
  diagnostic-only summaries.
- Extend diagnostic-screening and signal-detection evidence across larger
  slope regimes, sparse linkage patterns, sample sizes, and score-support
  stress before using the readouts as stronger screening recommendations.
- Add subgroup DFF fixtures and simulation operating-characteristic evidence
  before using DFF rows as fairness, invariance, or bias claims.
- Keep CRAN-time tests lightweight while storing longer design evidence under
  `inst/validation`.

Exit criteria:

- Design summaries report bounded-`GPCM` slope-regime and score-support
  conditions before performance metrics.
- Forecasting output does not imply direct posterior scoring for observed
  units or operational adequacy of a sensitivity model.
- Diagnostic-screening and signal-detection helpers remain
  `supported_with_caveat` until larger slope-aware operating-characteristic
  evidence is stored and linked from the capability matrix.
- DFF helpers retain `gpcm_boundary` in summaries, reports, and plot payloads
  until larger subgroup and external-fixture evidence is available.

### GPCM linking synthesis

Capability row:

- `Operational linking synthesis` (`supported_with_caveat`)

Surface available as bounded-`GPCM` exploratory synthesis:

- `build_linking_review()`

Required evidence before stronger operational wording:

- Define how anchor, drift, and chain evidence should behave when
  discrimination is free.
- Keep direct anchor/drift helpers available as exploratory inputs and keep
  `gpcm_boundary` visible on the combined review.
- Add examples that distinguish sparse-link design problems from fitted-model
  recovery failures.

Exit criteria:

- A bounded-`GPCM` linking review reports its assumptions and caveats before
  any combined decision or route recommendation.

### GPCM posterior predictive checks

Capability row:

- `MCMC and heavy-backend extensions` (`deferred`) for posterior predictive
  computation in the current matrix wording.

Required evidence before unblocking:

- Define bounded-`GPCM` posterior predictive discrepancy measures for marginal,
  pairwise, residual, and category-support checks.
- Document the replication mechanism and the conditioning set used for each
  check.
- Run false-positive and sensitivity reviews outside CRAN-time tests.
- Keep current strict marginal diagnostics labelled as exploratory screens
  until posterior predictive computation exists.

Exit criteria:

- Posterior predictive tables and plots are computed rather than merely named,
  and their interpretation remains separate from automatic pass/fail QC.

### GPCM engine and model-structure extensions

Capability row:

- `MCMC and heavy-backend extensions` (`deferred`)

Potential future scope:

- `slope_facet != step_facet`
- latent-regression bounded `GPCM`
- multidimensional population models
- MCMC or HMC engines
- compiled backend promotion where it changes runtime feasibility rather than
  the statistical contract

Required evidence before unblocking:

- Identification and covariance-basis tests for every additional model
  structure.
- Reduction tests for the unit-slope and constrained-step cases.
- User-facing documentation that separates the core package route from optional
  heavy-backend routes.

Exit criteria:

- New engine or model structures have explicit capability-matrix rows instead
  of silently broadening the current bounded route.

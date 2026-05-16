# External parameter-recovery simulation evidence for mfrmr 0.2.0

Review date: 2026-05-17.

This file records the package-maintainer review of the separate
`Parameter_Recovery_Simulation` workflow. The workflow itself is not bundled
with the package because it contains large generated datasets, engine outputs,
HTML reports, and optional FACETS batch files. This summary preserves the
release-relevant evidence and the limits of that evidence. The companion script
`external-recovery-audit.R` can be sourced to re-audit a local copy of the
external workflow outputs, check the expected CSV schemas, record file
fingerprints, and regenerate the compact summary tables.

## Scope reviewed

The external workflow is a common-data simulation framework for comparing
many-facet Rasch-family recovery across R, Python, Julia, and optionally FACETS.
The reviewed first-phase outputs focus on `RSM` / `PCM` fits with `JMLE`.

The reviewed analysis outputs came from:

- `analysis/dataset_manifest.csv`
- `analysis/smoke_check_summary.csv`
- `analysis/engine_status_summary.csv`
- `analysis/engine_parity_overview.csv`
- `analysis/key_findings.csv`
- `analysis/key_findings_counts.csv`
- `sample_size_dstudy/analysis/sample_size_decision_summary.csv`
- `sample_size_dstudy/analysis/sample_size_classification_summary.csv`

To refresh the review from a local external-workflow directory:

```r
source(system.file("validation", "external-recovery-audit.R", package = "mfrmr"))
review <- mfrmr_review_external_recovery_simulation(
  "/path/to/Parameter_Recovery_Simulation"
)
summary(review)
```

The main analysis manifest contains five classroom-writing datasets:
`baseline`, `nonrandom_missing`, `rater_drift`, `central_tendency`, and
`weak_bridge`, each with one rater per person. Dataset sizes range from 434 to
473 observed ratings, with observed density from 0.904 to 0.985. The generator
includes criterion-level discrimination differences, while the reviewed `RSM`
and `PCM` fits intentionally do not estimate discrimination. This is useful
misspecification evidence, not a claim that the equal-discrimination models are
adequate under all stress patterns.

## Structural checks

The smoke summary reported 52 of 52 expected output checks passing. The engine
status summary covered 30 engine/model/design groups:

- engines: R, Python, Julia
- models: `RSM`, `PCM`
- design patterns: five classroom-writing stress patterns
- error runs: 0
- minimum convergence rate: 1.00

Observed median runtime ranges in the reviewed output were:

- R: 0.205 to 0.398 seconds
- Python: 0.038 to 0.069 seconds
- Julia: 0.009 to 5.559 seconds, with first-run environment/startup cost visible

The runtime evidence supports the smoke/agreement workflow as a fast independent
check. It should not be read as CRAN-time package evidence, because the
workflow is external and uses separate generated data.

## Engine agreement evidence

R/Python/Julia agreement is strong for centered estimates, steps, standard errors,
and most separation summaries. In the reviewed `engine_parity_overview.csv`:

| Source | Maximum RMSE across agreement groups | Review groups |
|---|---:|---:|
| level recovery | 0.0052 | 0 / 480 |
| step recovery | 0.0034 | 0 / 360 |
| fit statistics | 0.0162 | 32 / 1200 |
| separation | 0.0127 | 4 / 1200 |

The larger review counts occur mainly in fit-standardization and separation
summaries, especially standardized fit quantities where degrees-of-freedom
conventions are known to matter. This supports the package's documentation
choice to treat FACETS-style fit comparison as a convention-aware review rather
than a simple numerical identity claim.

## Stress-pattern findings

The automatically ranked key findings identify the expected weak points of the
one-rater sparse classroom-writing design:

- rater drift plus one rating per person can produce person recovery RMSE near
  0.557 logits and true-band bias near 0.516 logits.
- extreme persons show strong centerward shrinkage under rater-drift stress,
  with directional extremity bias near 0.560 logits.
- nonrandom missingness can create large rater coverage-bias contrasts, with
  absolute contrast to high-coverage rows near 0.683 logits.
- weak bridge designs can produce low-discrimination criterion role-bias
  contrasts near 0.256 logits.
- task separation and reliability can collapse to zero in some one-rater
  classroom-writing stress conditions.

These findings are not package failures. They are design and misspecification
warnings that should remain visible in user guidance: successful convergence
and cross-engine agreement do not imply that a sparse or stressed design gives
adequate recovery for every facet.

## Sample-size D-study evidence

The reviewed sample-size D-study outputs are a compact smoke run over person
count and rater-pool count, again for classroom-writing `RSM` fits across
R/Python/Julia. It is useful for checking that the reporting pipeline produces
decision-oriented summaries, but it is not a final operating-characteristic
study because it uses one replication per cell.

In the reviewed outputs:

- increasing persons from 30 to 60 reduced person RMSE from about 0.366 to
  0.307 and increased cut-score accuracy from about 0.83 to 0.93-0.95.
- increasing the rater pool from 4 to 8 did not uniformly improve every metric;
  some rater and task summaries were noisier at the larger pool size in the
  one-replication smoke run.
- classification risk was split between `low` and `moderate` across person and
  rater sample-size conditions.

The package should therefore describe sample-size output as a decision support
surface that needs replication, uncertainty bands, and design-specific
thresholds before being used for operational planning.

## Release implications

This external simulation evidence supports the 0.2.0 release boundary in three
ways.

First, it strengthens the distinction between direct recovery and design
endorsement. `evaluate_mfrm_recovery()` and `assess_mfrm_recovery()` should stay
framed as parameter-recovery checks under a stated data-generating design.

Second, it supports the warning that model fit, convergence, and agreement are not
enough. The reviewed stress cases converge and agree across engines while still
showing recovery, coverage, precision, and role-bias risks.

Third, it supports keeping GPCM planning and broader score-side workflows
deferred until their estimands and uncertainty behavior are validated. The
reviewed first-phase evidence is strongest for `RSM` / `PCM` JMLE agreement and
stress-pattern sensitivity, not for unrestricted GPCM operational reporting.

## Review limits

- The main reviewed analysis manifest contains five datasets, not the full
  generated-data inventory.
- The sample-size D-study evidence is a smoke run with one replication per
  cell; uncertainty intervals are therefore not stable.
- The reviewed first-phase fits use `RSM` / `PCM` with `JMLE`. They do not
  validate bounded-`GPCM` score-side exports, posterior predictive checks,
  planning/forecasting helpers, or APA/QC pipelines.
- FACETS files are present in the external workflow, but this summary does not
  make a FACETS numerical-reproduction claim. FACETS comparisons still require
  imported external outputs and explicit comparison contracts.

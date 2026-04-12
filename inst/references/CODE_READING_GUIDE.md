# Developer Code Reading Guide

This guide is for developers who want the shortest reliable route through the
current `mfrmr` codebase.

It reflects the package structure as of 2026-04-10. Earlier notes that point
to `/R/api.R` are obsolete.

## 1) Start with the public entry points

Read these files in order:

1. `/R/api-estimation.R`
2. `/R/api-prediction.R`
3. `/R/api-simulation.R`
4. `/R/api-reports.R`
5. `/R/api-advanced.R`

Use them to answer:

- what the public arguments are;
- which workflows are intended to be called directly;
- which boundaries are documented versus merely internal.

If you only need the ordinary estimation route, stop after
`api-estimation.R`, `mfrm_core.R`, and the workflow vignette/help pages.

## 2) Core estimation flow

Core numeric logic lives in `/R/mfrm_core.R`.

The practical public flow is:

1. `fit_mfrm()` in `/R/api-estimation.R`
2. `prepare_mfrm_population_scaffold()` in `/R/api-estimation.R`
3. `mfrm_estimate()` in `/R/mfrm_core.R`
4. `resolve_person_quadrature_basis()` in `/R/mfrm_core.R`
5. `finalize_mfrm_population_fit()` in `/R/api-estimation.R`

`mfrm_estimate()` is still the central orchestration point. Read it as:

1. resolve model/configuration choices;
2. prepare design objects and constraints;
3. choose the estimation engine;
4. optimize;
5. convert optimized parameters into public tables and summaries.

## 3) Latent-regression route

The current latent-regression branch is implemented, but only for a narrow
first-version `MML` path.

Read these pieces together:

- `/R/api-estimation.R`
  - `fit_mfrm()`
  - `prepare_mfrm_population_scaffold()`
  - `finalize_mfrm_population_fit()`
- `/R/mfrm_core.R`
  - `resolve_person_quadrature_basis()`
  - population-aware `MML` helpers around `mfrm_estimate()`
- `/R/api-prediction.R`
  - `predict_mfrm_units()`
  - `sample_mfrm_plausible_values()`

Current boundary:

- ordered-response `RSM` / `PCM`;
- unidimensional `MML`;
- conditional-normal person population model;
- person-level covariates from explicit one-row-per-person `person_data`.

Do not read the current branch as ConQuest numerical equivalence. Use
`build_conquest_overlap_bundle()`, `normalize_conquest_overlap_*()`, and
`audit_conquest_overlap()` only for the documented scoped comparison workflow.

## 4) Bounded GPCM route

The bounded `GPCM` branch is real, but deliberately limited.

Read these together:

- `/R/help_gpcm_scope.R`
- `/R/api-estimation.R`
- `/R/api-prediction.R`
- `/R/api-simulation.R`
- `/R/api-reports.R`
- `/R/api-export-bundles.R`
- `/R/mfrm_core.R` for `stop_if_gpcm_out_of_scope()`

Use `gpcm_capability_matrix()` as the release-scope truth source.

Important distinction:

- lower-level probability, scoring, and direct simulation paths exist;
- score-side reporting semantics, broader export bundles, and planning /
  forecasting remain blocked or deferred.

## 5) Diagnostics, reporting, and operational review

The package is no longer organized around one reporting file.

Read by workflow instead:

- diagnostics and fit objects:
  - `/R/api-estimation.R`
  - `/R/mfrm_core.R`
- reporting and checklists:
  - `/R/api-reports.R`
  - `/R/api-reporting-checklist.R`
- dashboards and plotting:
  - `/R/api-dashboards.R`
  - `/R/api-plotting.R`
  - `/R/api-plotting-anchor-qc.R`
  - `/R/api-plotting-fit-family.R`
- linking, DFF, misfit, weighting:
  - `/R/api-advanced.R`
- compatibility wrappers:
  - `/R/facets_mode_api.R`
  - `/R/facets_mode_methods.R`

## 6) Simulation and planning

Simulation support is broad, but the planner is narrower than the estimation
core.

Read:

- `/R/api-simulation-spec.R`
- `/R/api-simulation.R`
- `/tests/testthat/test-simulation-design.R`

Current boundary:

- direct simulation-spec generation is broad and includes bounded `GPCM`;
- design planning / forecasting is still validated for the role-based
  two-non-person-facet `RSM` / `PCM` route;
- arbitrary-facet planning is still scaffold/spec work, not finished public
  behavior.

## 7) Where to edit safely

If you need to change behavior, start here:

- public signatures and help wording:
  - `/R/api-estimation.R`
  - `/R/api-prediction.R`
  - `/R/help_gpcm_scope.R`
  - `/R/mfrmr-package.R`
- core estimation math:
  - `/R/mfrm_core.R`
- simulation/planning behavior:
  - `/R/api-simulation-spec.R`
  - `/R/api-simulation.R`
- reporting / export contracts:
  - `/R/api-reports.R`
  - `/R/api-export-bundles.R`
- linking / weighting / case-review behavior:
  - `/R/api-advanced.R`

Before changing a large surface, check for local edits and avoid mixing
unrelated documentation and computation changes.

## 8) Verification route

Prefer targeted package-aware test runs, not raw `testthat::test_dir()` from an
unloaded session.

Typical route:

1. `Rscript -e 'devtools::test(filter = "estimation-core|prediction")'`
2. `Rscript -e 'devtools::test(filter = "gpcm-capability-matrix|mml-cpp11-backend|reporting-checklist")'`
3. `Rscript -e 'devtools::test(filter = "simulation-design|reference-benchmark")'`
4. `R CMD check --as-cran` only after the targeted surfaces are stable

If you touch documentation-only markdown under `/inst/references`, you do not
need to regenerate package docs immediately. If you touch roxygen comments in
`/R`, regenerate `.Rd` files before a release-facing check.

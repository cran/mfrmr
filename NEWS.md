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

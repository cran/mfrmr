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

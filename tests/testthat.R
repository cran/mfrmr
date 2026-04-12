library(testthat)
library(mfrmr)

is_cran_check <- local({
  env <- Sys.getenv("NOT_CRAN")
  if (identical(env, "")) {
    !interactive()
  } else {
    !isTRUE(as.logical(env))
  }
})

cran_filter <- paste(
  c(
    "anchor-equating",
    "api-coverage-gaps",
    "api-final-coverage",
    "bubble-chart",
    "bundle-coverage",
    "core-coverage",
    "core-coverage-gaps",
    "core-workflow",
    "coverage-push-95",
    "data-processing",
    "diagnostic-screening-validation",
    "dif-module",
    "draw-coverage",
    "edge-cases",
    "error-handling",
    "estimation-core",
    "expanded-summary-plot",
    "exception-regression",
    "export-bundles",
    "facet-dashboard",
    "facet-equivalence",
    "facets-column-contract",
    "facets-metric-contract",
    "facets-parity-report",
    "facets-mode-api",
    "final-coverage-boost",
    "identifiability-constraints",
    "marginal-fit-diagnostics",
    "marginal-fit-plots",
    "misfit-casebook",
    "numerical-validation",
    "output-stability",
    "parameter-recovery",
    "plot-customization",
    "prediction",
    "qc-pipeline",
    "reference-benchmark",
    "remaining-coverage",
    "report-functions",
    "reporting-coverage",
    "reporting-checklist",
    "reporting-gaps",
    "simulation-design",
    "summary-table-bundle"
  ),
  collapse = "|"
)

if (is_cran_check) {
  # Keep CRAN checks under the check-farm time budget by skipping
  # long integration and coverage-expansion suites.
  test_check("mfrmr", filter = cran_filter, invert = TRUE)
} else {
  test_check("mfrmr")
}

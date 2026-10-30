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

# Keep live estimation, independent numerical checks and representative public paths.
# Broader report/plot/portable-scoring combinations run in the complete suite.
# See CONTRIBUTING.md for representative cases and full-suite commands.
# Always retain posterior-summary: archive encoding can break its public reader
# even when loaded values are identical.
cran_light_tests <- c(
  "cran-smoke",
  "static-table-isolation",
  "model-export",
  "posterior-summary",
  "adaptive-quadrature-review",
  "calibration-public-api",
  "classification",
  "bias-numerics",
  "marginal-fit-availability",
  "interaction-em",
  "likelihood-numerics",
  "source-posterior-numerics",
  "preparation-invariants",
  "anchor-migration",
  "cml",
  "compiled-header-contract",
  "data-and-citation",
  "results-csv-encoding",
  "feature-clustering",
  "linking-numerics",
  "random-rater",
  "random-rater-profile",
  "testlet",
  "gauss-hermite-weights",
  "generalizability-missing",
  "multivariate-gtheory",
  "gstudy-sampling-vcov",
  "gstudy-example-replay",
  "dstudy-delta",
  "multivariate-nested-gtheory",
  "icc-input",
  "mantel-haenszel",
  "ordinal-logistic",
  "mml-cpp11-backend",
  "bundle-summary-privacy",
  "gpcm-ic-eligibility",
  "gpcm-lrt",
  "gpcm-uncertainty-readiness",
  "gpcm-inference-extensions",
  "gpcm-pattern-cache",
  "gmfrm-public-workflow",
  "gmfrm-em",
  "jml-public-workflow",
  "residual-network",
  "prefit-network",
  "krippendorff-alpha",
  "namespace-contract",
  "optimizer-curvature",
  "vignette-artifacts",
  "facets-results"
)

# These detailed files remain in the unfiltered NOT_CRAN=true suite.
extended_regression_tests <- c(
  "api-migration",
  "conditional-scoring",
  "mml-starts",
  "parameter-map-invariants",
  "anchor-centering",
  "anchored-calibration",
  "anchored-gpcm",
  "cml-comparisons",
  "comparison-review",
  "design-evaluation-denominators",
  "response-imputation",
  "missing-codes-integration",
  "study-report",
  "network-inference",
  "model-export-affine",
  "model-export-summary",
  "decision-probabilities",
  "facets-manual-consistency",
  "extended-example-replay",
  "simulation-example-replay",
  "dstudy-summary",
  "gmfrm-portable",
  "planning-report-connections",
  "report-modules"
)

# Research reference code stays in the repository and runs separately after
# full package checks, using the same installed package.
reference_tests <- c(
  "gmfrm-bounded-search",
  "gmfrm-candidate",
  "gmfrm-coordinates",
  "gmfrm-interior-finish",
  "gmfrm-natural-search",
  "gmfrm-natural",
  "gmfrm-procedure",
  "gmfrm-reference-patterns",
  "gmfrm-result",
  "gmfrm-run",
  "gmfrm-segmented-optimizer",
  "gmfrm-segmented",
  "gmfrm-selection",
  "gmfrm-wald",
  "gmfrm-zero-location",
  "gmfrm-zero-profile-search",
  "gmfrm-zero-profile",
  "jml-gpcm-binary-closure-envelope",
  "jml-gpcm-higher-order-face-lifts",
  "jml-gpcm-response-image-face-chart",
  "jml-gpcm-response-quotient-closure",
  "jml-gpcm-saturated-response-envelope",
  "jml-gpcm-boundary-compactification-reference",
  "jml-gpcm-exponential-balance-reference",
  "jml-gpcm-global-existence-reference",
  "jml-gpcm-lexicographic-limit-reference",
  "jml-gpcm-parameter-path-reachability-reference",
  "jml-gpcm-parameter-sequence-flag-reference",
  "jml-gpcm-rate-hierarchy-reference",
  "jml-gpcm-sequence-remainder-diagnostic-reference",
  "fixed-calibration-anchors-reference"
)

cran_light_filter <- paste0(
  "(^|/)(test-)?(",
  paste(cran_light_tests, collapse = "|"),
  ")$"
)

if (is_cran_check) {
  # Exercise the public data -> declared-design review -> fit -> summary ->
  # Wright/pathway plot -> export route once, plus lightweight compatibility,
  # backend, privacy, and installed-artifact contracts. Portable scoring,
  # dashboard restrictions, consolidated plot dispatch on saved results, and
  # guarded optimizer restarts are included.
  # Detailed plotting, repeated estimation/diagnostics, simulation/recovery,
  # documentation scans and broad regression tests remain in the complete
  # local/CI suite, which must be run with NOT_CRAN=true.
  local({
    # Short files share one session; do not start parallel worker libraries.
    withr::local_options(Ncpus = 1L)
    withr::local_envvar(c(
      TESTTHAT_PARALLEL = "false", TESTTHAT_CPUS = "1",
      OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
      BLIS_NUM_THREADS = "1"
    ))
    test_check("mfrmr", filter = cran_light_filter)
  })
} else {
  test_check("mfrmr")
}

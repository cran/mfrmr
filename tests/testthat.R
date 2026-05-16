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

cran_light_tests <- c(
  "compatibility-aliases",
  "data-and-citation",
  "gpcm-capability-matrix",
  "namespace-contract"
)

cran_light_filter <- paste0(
  "(^|/)(test-)?(",
  paste(cran_light_tests, collapse = "|"),
  ")$"
)

if (is_cran_check) {
  # Keep CRAN checks under the check-farm time budget by running only
  # lightweight metadata/package-contract checks. Run the complete suite
  # locally/CI with NOT_CRAN=true.
  test_check("mfrmr", filter = cran_light_filter)
} else {
  test_check("mfrmr")
}

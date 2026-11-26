library(testthat)
library(mfrmr)

is_cran_check <- local({
  env <- Sys.getenv("NOT_CRAN")
  if (identical(env, "")) !interactive() else !isTRUE(as.logical(env))
})

# Live native numerics, S3 contracts, inference guards and one complete workflow.
# Broader regressions remain in this directory and run with NOT_CRAN=true.
cran_light_tests <- c(
  "engine-032",
  "rating-design-display",
  "tidymfrm",
  "report-statistics",
  "mml-population-boundary",
  "step-constraints",
  "classification",
  "markov-drift",
  "unified-contracts",
  "public-polish"
)

cran_light_filter <- paste0("(^|/)(test-)?(", paste(cran_light_tests, collapse = "|"), ")$")

if (is_cran_check) {
  local({
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

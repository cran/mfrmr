# Regenerate the synthetic results used by the extended-model help examples.
# This script defines a function; sourcing it does not start any calculation.
# Requires mfrmr and optional RTMB >= 2.0; no software is installed.
# Full calibration, shared-rater scoring and bootstrap refits can take minutes.
#
# source(system.file("examples", "extended-models.R", package = "mfrmr"))
# RNGkind("Mersenne-Twister", "Inversion", "Rejection") # recorded bootstrap RNG
# example <- extended_model_examples()
# saveRDS(example, "extended-models.rds", compress = "xz")
extended_model_examples <- function() {
  if (!requireNamespace("RTMB", quietly = TRUE) ||
      utils::packageVersion("RTMB") < "2.0") {
    stop("Install optional RTMB >= 2.0 to regenerate the shared-rater results.")
  }
  ratings <- mfrmr::load_mfrmr_data("example_core")
  testlet_fit <- mfrmr::fit_mfrm_testlet(
    ratings, "Person", "Score", "Rater", c("Rater", "Criterion"), 1:4,
    quad_points = 121
  )
  random_fit <- mfrmr::fit_mfrm_random_rater(
    ratings, "Person", "Rater", "Score", "Criterion", 1:4,
    quad_points = 121
  )
  stopifnot(isTRUE(testlet_fit$checks$NumericalReady),
            isTRUE(testlet_fit$checks$InformationPositive),
            isTRUE(random_fit$checks$NumericalReady),
            isTRUE(random_fit$checks$InformationPositive))
  list(
    testlet = list(
      fit = testlet_fit,
      scores = stats::predict(testlet_fit,
        persons = as.character(unique(ratings$Person)[1:4])),
      diagnostics = mfrmr::mfrm_response_diagnostics(testlet_fit, group_by = "Rater")
    ),
    random_rater = list(
      fit = random_fit,
      scores = mfrmr::score_mfrm_random_rater(random_fit,
        persons = as.character(unique(ratings$Person)[1:2])),
      # Nineteen refits illustrate reuse and failure accounting only.
      # They cannot resolve 2.5% tails accurately; no failures are discarded.
      intervals = mfrmr::mfrm_random_rater_intervals(random_fit,
        nsim = 19, seed = 923701)
    )
  )
}

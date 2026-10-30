# Analytic likelihoods exercise the public search independently of RTMB.
random_rater_profile_fixture <- function(center, nll = 0) structure(list(
  settings = list(quad_points = 31L, maxit = 100L, fixed_person_sd = 1),
  checks = list(NumericalReady = TRUE, InformationPositive = TRUE),
  calibration = list(beta = numeric(), steps = 0, rater_sd = center, person_sd = 1),
  loglik = -nll, input = list()), class = "mfrm_random_rater")

random_rater_quadratic_profile <- function(center, se, nll = 0, level = .95) {
  objective <- function(...) list(
    fn = function(p) nll + p[1]^2/2 + (p[2] - center)^2/(2*se^2),
    gr = function(p) c(p[1], (p[2] - center)/se^2))
  testthat::with_mocked_bindings(confint(random_rater_profile_fixture(center, nll), level = level),
    mfrm_random_rater_objective = objective, .package = "mfrmr")
}

test_that("rater SD profile endpoints attain analytic likelihood-ratio cutoffs", {
  for (se in c(.001, 1e-7, 1e-10)) {
    ci <- random_rater_quadratic_profile(.01, se)
    expect_equal((as.numeric(ci) - .01)/se, c(-1,1)*qnorm(.975), tolerance = 2e-5)
    expect_lt(ci[1,1], .01)
    expect_gt(ci[1,2], .01)
    tab <- attr(ci, "profile")
    ends <- match(as.numeric(ci), tab$SD)
    expect_equal(tab$LR[ends], rep(qchisq(.95,1), 2), tolerance = 1e-4)
    expect_true(all(tab$InformationPositive))
    expect_true(all(tab$NuisanceMinEigenvalue > 0 & tab$ReferenceMinEigenvalue > 0))
  }
  # Small requested mass is not rounded into the source NLL before subtraction.
  ci <- random_rater_quadratic_profile(1, .1, level = 1e-7)
  tab <- attr(ci, "profile")
  expect_equal(tab$LR[match(as.numeric(ci), tab$SD)]/qchisq(1e-7,1), c(1,1),
    tolerance = 1e-4)
})

test_that("unresolved profile precision retains evidence instead of coincident bounds", {
  error <- tryCatch(random_rater_quadratic_profile(1, .1, nll = 1000, level = 1e-12),
    mfrm_profile_error = identity)
  expect_s3_class(error, "mfrm_profile_error")
  expect_match(conditionMessage(error), "likelihood-ratio cutoff.*numerical resolution")
  expect_true(is.data.frame(error$profile))
  expect_gt(nrow(error$profile), 1)
  expect_true(all(error$profile$NumericalReady))
  expect_true(all(is.finite(error$profile$LR)))
  # The failure and its candidate diagnostics survive a fresh read.
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(error, path)
  expect_identical(readRDS(path)$profile, error$profile)
  expect_error(random_rater_quadratic_profile(1, .1, level = .Machine$double.xmin),
    "cutoff.*numerical resolution", class = "mfrm_profile_error")
  # A baseline mismatch small on the usual NLL scale can still dominate a
  # small requested LR cutoff; do not silently widen that interval.
  objective <- function(...) list(fn = function(p) 1000 + sum((p-c(0,1))^2)/2,
    gr = function(p) p-c(0,1))
  testthat::with_mocked_bindings(
    expect_error(confint(random_rater_profile_fixture(1, 1000+1e-8), level = 1e-7),
      "inconsistent with the source likelihood", class = "mfrm_profile_error"),
    mfrm_random_rater_objective = objective, .package = "mfrmr")
})

test_that("rater profiles do not accept a stationary nuisance maximum", {
  # At the source (nuisance=0, SD=2) the joint information is positive.
  # At |SD-2| > pi/2 the nuisance optimum is +/-sqrt(-cos(SD-2)).
  # Starting at zero still reports convergence and a zero gradient.
  objective <- function(...) list(
    fn = function(p) (p[2]-2)^2/2 + cos(p[2]-2)*p[1]^2/2 + p[1]^4/4,
    gr = function(p) c(cos(p[2]-2)*p[1]+p[1]^3,
      p[2]-2-sin(p[2]-2)*p[1]^2/2))
  error <- testthat::with_mocked_bindings(
    tryCatch(confint(random_rater_profile_fixture(2)), mfrm_profile_error = identity),
    mfrm_random_rater_objective = objective, .package = "mfrmr")
  expect_s3_class(error, "mfrm_profile_error")
  bad <- subset(error$profile, !NumericalReady)
  expect_identical(bad$OptimizerCode, 0L)
  expect_equal(bad$MaxGradient, 0)
  # The default central difference of the cubic gradient adds h^2 = 1e-6.
  expect_lte(abs(bad$NuisanceMinEigenvalue - cos(2)), 1.01e-6)
  expect_false(bad$InformationPositive)
  expect_match(bad$Detail, "curvature.*not positive")
})

test_that("known-ability profiles verify curvature at both integration orders", {
  for (case in c("negative", "flat", "unavailable")) {
    objective <- function(input, quad_points, ...) {
      sign <- if (quad_points <= 31) 1 else switch(case, negative = -1, flat = 0, unavailable = 1)
      list(fn = function(p) (p[2]-1)^2/2 + sign*p[1]^2/2,
        gr = function(p) {
          if (quad_points > 31 && case == "unavailable" && p[1] != 0) stop("missing derivative")
          c(sign*p[1], p[2]-1)
        })
    }
    error <- testthat::with_mocked_bindings(
      tryCatch(confint(random_rater_profile_fixture(1)), mfrm_profile_error = identity),
      mfrm_random_rater_objective = objective, .package = "mfrmr")
    expect_s3_class(error, "mfrm_profile_error")
    tab <- error$profile
    expect_equal(tab$LogLikDifference, 0)
    expect_equal(tab$GradientDifference, 0)
    expect_false(tab$InformationPositive)
    expect_false(tab$NumericalReady)
    expect_match(tab$Detail, if (case == "unavailable") "could not be verified: missing derivative" else "not positive definite")
  }
})

test_that("a supported zero lower limit does not have to cross the cutoff", {
  for (center in c(0, .1)) {
    ci <- random_rater_quadratic_profile(center, .2)
    expect_identical(unname(ci[1,1]), 0)
    expect_lte(abs(ci[1,2] - (center + qnorm(.975)*.2)), 6e-6)
    tab <- attr(ci, "profile")
    expect_lte(tab$LR[tab$SD == 0], qchisq(.95,1))
    expect_lte(abs(tab$LR[tab$SD == ci[1,2]] - qchisq(.95,1)), 1e-4)
    expect_true(all(tab$NumericalReady))
  }
})

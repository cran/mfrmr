coefficient_interval_fixture <- function() {
  readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits[["PCM-MML-fixed"]]
}

test_that("default intervals use the same free coefficients and covariance", {
  f <- coefficient_interval_fixture()
  b <- coef(f); V <- vcov(f)
  ci <- confint(f, level = .90)
  expect_identical(rownames(ci), names(b))
  expect_identical(colnames(ci), c("5%", "95%"))
  expect_equal(as.numeric(ci[, 1]), unname(b - qnorm(.95) * sqrt(diag(V))))
  expect_equal(as.numeric(ci[, 2]), unname(b + qnorm(.95) * sqrt(diag(V))))
  expect_true(all(attr(ci, "diagnostics")$CIEligible))
  expect_identical(attr(ci, "checks"), attr(V, "checks"))
  expect_identical(attr(ci, "numerical_checks"), attr(V, "numerical_checks"))
  expect_match(attr(ci, "target"), "Free model coordinates")
  expect_match(attr(ci, "method"), "MML")
  selected <- c(3L, 1L)
  a <- confint(f, selected, .90)
  named <- confint(f, names(b)[selected], .90)
  expect_identical(a, named)
  expect_equal(as.numeric(a), as.numeric(ci[selected, , drop = FALSE]))
  expect_identical(attr(a, "diagnostics")$Parameter, names(b)[selected])
})

test_that("CML coefficient intervals agree with the closed-form paired binary solution", {
  skip_if_not_installed("lpSolve")
  patterns <- rbind(c(1, 0), c(0, 1), c(0, 0), c(1, 1))
  x <- patterns[rep(seq_len(4L), c(3L, 7L, 2L, 3L)), , drop = FALSE]
  d <- data.frame(Person = rep(seq_len(nrow(x)), each = 2L),
    Rater = rep(c("A", "B"), nrow(x)), Score = as.vector(t(x)))
  f <- fit_mfrm(d, "Person", "Rater", "Score", method = "CML", category_policy = "preserve")
  ci <- confint(f, level = .90)
  estimate <- log(7 / 3) / 2
  se <- sqrt(10 / (4 * 3 * 7))
  expect_equal(as.numeric(ci), estimate + c(-1, 1) * qnorm(.95) * se, tolerance = 1e-6)
  expect_identical(rownames(ci), names(coef(f)))
  expect_match(attr(ci, "method"), "Conditional-information")
  expect_identical(attr(ci, "checks"), f$checks)
  expect_identical(attr(ci, "cautions"), f$notes)
  expect_true(all(attr(ci, "diagnostics")$CIEligible))
  for (variance in c(0, -1, NA_real_, Inf)) {
    bad <- f; bad$covariance[1, 1] <- variance
    withheld <- confint(bad)
    expect_true(all(is.na(withheld)))
    expect_false(any(attr(withheld, "diagnostics")$CIEligible))
    expect_match(attr(withheld, "diagnostics")$InferenceReview[1], "positive finite sampling variance")
  }
  expect_error(confint(f, scale = "standardized"), "Slope-specific options require")
  expect_error(confint(f, typo = TRUE), "must be empty")
})

test_that("conditional comparison intervals retain unsupported targets and their reasons", {
  f <- readRDS(test_path("fixtures", "cml-comparison.rds"))
  ci <- confint(f, level = .80)
  existing <- mfrm_facet_intervals(f, level = .80)$table
  expect_identical(rownames(ci), names(coef(f)))
  expect_equal(as.numeric(ci[, 1]), existing$Lower)
  expect_equal(as.numeric(ci[, 2]), existing$Upper)
  expect_true(all(is.na(ci[1, ])))
  expect_false(attr(ci, "diagnostics")$CIEligible[1])
  expect_identical(attr(ci, "diagnostics")$InferenceReview[1], f$comparisons$EstimateReason[1])
  expect_identical(rownames(confint(f, 2L)), names(coef(f))[2L])
  withheld <- f; withheld$comparisons$FormalInferenceEligible[2] <- FALSE
  withheld$comparisons$EstimateReason[2] <- "Stored target did not pass inference checks."
  expect_true(all(is.na(confint(withheld)[2, ])))
  expect_match(attr(confint(withheld), "diagnostics")$InferenceReview[2], "did not pass")
  unknown <- f; unknown$comparisons$FormalInferenceEligible <- NULL
  expect_true(all(is.na(confint(unknown))))
  expect_match(attr(confint(unknown), "diagnostics")$InferenceReview[2], "eligibility is unavailable")
  unknown$comparisons$FormalInferenceEligible <- rep("TRUE", nrow(unknown$comparisons))
  expect_true(all(is.na(confint(unknown))))
  expect_match(attr(confint(unknown), "diagnostics")$InferenceReview[2], "eligibility is unavailable")
  reordered <- f; reordered$covariance <- f$covariance[nrow(f$covariance):1, ncol(f$covariance):1]
  expect_identical(confint(reordered, level = .80), ci)
  bad_labels <- f; rownames(bad_labels$covariance)[1] <- "wrong coordinate"
  expect_error(confint(bad_labels), "dimensions and names must match")
  for (case in readRDS(test_path("fixtures", "partial-comparisons.rds")))
    expect_error(confint(case$fit), "native CML.*descriptive estimates")
})

test_that("GPCM slopes remain an explicitly selected transformed target", {
  f <- readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds"))
  expect_warning(free <- confint(f),
    "pointwise Wald inference assumes.*does not guarantee finite-sample coverage")
  slope <- confint(f, parm = "slopes")
  expect_s3_class(slope, "mfrm_slope_intervals")
  expect_false(inherits(free, "mfrm_slope_intervals"))
  expect_identical(rownames(free), names(coef(f)))
  expect_equal(nrow(slope), nrow(f$slopes))
  select <- grep("^log_slopes\\[", rownames(free))
  expect_equal(exp(as.numeric(free[select, ])),
    as.numeric(slope[seq_along(select), ]), tolerance = 1e-7)
  expect_true(any(grepl("coverage", attr(free, "cautions"))))
  for (args in list(list(scale = "standardized"), list(method = "profile"),
      list(method = "model"), list(contrasts = diag(3)), list(simultaneous = "none"),
      list(adjust = FALSE), list(slope = "A"), list(profile_control = list()),
      list(clusters = NULL), list(contrast_scale = "ratio")))
    expect_error(do.call(confint, c(list(object = f), args)), "Slope-specific options require parm = 'slopes'")
  # A profile request still forwards the same object, target and controls;
  # it does not accidentally take the free-coordinate Wald route.
  local_mocked_bindings(mfrm_gpcm_profile_interval = function(object, slope, level, profile_control)
    list(object = object, slope = slope, level = level, control = profile_control))
  prof <- confint(f, parm = "slopes", method = "profile", slope = "C01", profile_control = list(maxit = 80))
  expect_identical(prof$object, f)
  expect_identical(prof$slope, "C01")
  expect_identical(prof$control, list(maxit = 80))
})

test_that("invalid requests and unsupported inference do not fall through", {
  f <- coefficient_interval_fixture()
  for (parm in list(character(), numeric(), c(1, 1), 0, -1, 1.5, Inf, NA_real_,
      NA_character_, "unknown", "all", TRUE, factor("Rater[1]"), matrix(1L)))
    expect_error(confint(f, parm = parm), "parm")
  for (level in list(0, 1, NA_real_, Inf, c(.8, .9), "95%", matrix(.95), 1i))
    expect_error(confint(f, level = level), "level")
  expect_error(confint(f, typo = TRUE), "must be empty")
  j <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits[["PCM-JML"]]
  expect_error(confint(j), "JML point estimates")
  imported <- f; class(imported) <- c("mfrm_imported_fit", class(f))
  expect_error(confint(imported), "native fit_mfrm")
  local_mocked_bindings(mfrm_facet_intervals = function(...) list(
    parameter_covariance = diag(length(f$opt$par)),
    checks = data.frame(Passed = FALSE, Detail = "Unqualified numerical information.")))
  expect_error(confint(f), "Joint covariance is unavailable.*Unqualified")
})

test_that("coefficient interval labels and availability survive serialization", {
  f <- readRDS(test_path("fixtures", "cml-comparison.rds"))
  before <- serialize(f, NULL)
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(list(fit = f, ci = confint(f)), path)
  saved <- readRDS(path)
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected refit"),
    mfrmr_cml_evaluate = function(...) stop("Unexpected CML reevaluation"))
  expect_identical(confint(saved$fit), saved$ci)
  expect_identical(serialize(f, NULL), before)
})

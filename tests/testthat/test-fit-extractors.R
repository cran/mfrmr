extractor_fits <- function() readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits

test_that("ordinary coefficients and covariance use matching free coordinates", {
  f <- extractor_fits()[["PCM-MML-fixed"]]
  expect_warning(b <- coef(f), NA)
  expect_equal(unname(b), unname(f$opt$par))
  expect_true(all(is.finite(b)))
  V <- vcov(f)
  expect_identical(rownames(V), names(b))
  expect_identical(colnames(V), names(b))
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))
  # A difference of the first two raters is independent of the centering
  # constraint. Its variance must include their joint covariance.
  C <- matrix(c(1, -1, 0), 1L, dimnames = list("first-second", f$prep$levels$Rater))
  ci <- mfrm_facet_intervals(f, "Rater", C)
  expect_equal(unname(b["Rater[1]"] - b["Rater[2]"]), ci$table$Estimate)
  expect_equal(unname(V[1, 1] + V[2, 2] - 2 * V[1, 2]), ci$table$SE^2)
  for (type in c("locations", "steps", "persons")) {
    tab <- switch(type, locations = f$facets$others, steps = f$steps, persons = f$facets$person)
    expect_identical(unname(coef(f, type)), as.numeric(tab$Estimate))
    expect_false(anyDuplicated(names(coef(f, type))) > 0)
  }
  expect_error(coef(f, "slopes"), "no saved slopes")
  expect_error(coef(f, typo = TRUE), "must be empty")
  expect_error(nobs(f, type = "rows"), "arg")
})

test_that("population and slope coefficients retain their fitted scales", {
  f <- readRDS(test_path("fixtures", "mml-comparison-intervals.rds"))$RSM$fit
  b <- coef(f)
  expect_equal(unname(b["log_sigma2[1]"]), log(f$config$population_spec$sigma2))
  expect_equal(unname(b["beta[1]"]), unname(f$config$population_spec$coefficients))
  expect_error(vcov(f), "no interactions")
  g <- readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds"))
  expect_equal(unname(coef(g, "slopes")), g$slopes$Estimate)
  # The centered log-slope block omits the last dependent level.
  free <- coef(g)[grepl("^log_slopes\\[", names(coef(g)))]
  expect_equal(unname(free), head(log(g$slopes$Estimate), -1L))
  expect_warning(V <- vcov(g), "numerical eligibility does not guarantee finite-sample coverage")
  expect_identical(rownames(V), names(coef(g)))
  expect_true(all(attr(V, "checks")$Passed))
  expect_true(is.data.frame(attr(V, "numerical_checks")))
  expect_true(any(grepl("coverage", attr(V, "cautions"))))
})

test_that("likelihood extraction preserves the MML IC sample unit and refusals", {
  fits <- extractor_fits(); f <- fits[["PCM-MML-fixed"]]
  ll <- logLik(f)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -f$opt$value)
  expect_identical(attr(ll, "likelihood_basis"), "marginal")
  expect_identical(nobs(f), length(unique(f$prep$data$Person)))
  expect_identical(attr(ll, "nobs"), nobs(f))
  expect_identical(nobs(f, "ratings"), nrow(f$prep$data))
  expect_equal(AIC(f), f$summary$AIC)
  expect_equal(BIC(f), f$summary$BIC)
  expect_equal(AIC(f, k = 3), -2 * as.numeric(ll) + 3 * length(coef(f)))
  expect_gt(abs(BIC(f) - (-2 * as.numeric(ll) + log(nobs(f, "ratings")) * attr(ll, "df"))), 1)
  expect_error(AIC(f, f), "compare_mfrm")
  expect_error(BIC(f, f), "compare_mfrm")
  expect_error(AIC(f, k = NA_real_), "finite numeric")
  expect_error(logLik(f, typo = TRUE), "must be empty")
  changed <- f; changed$opt$value <- changed$opt$value + 1
  expect_error(logLik(changed), "missing or inconsistent")
  legacy <- f; legacy$summary$ICContractVersion <- NULL
  expect_error(logLik(legacy), "legacy")
  coarse <- fit_mfrm(f$prep$data, "Person", c("Rater", "Criterion"), "Score",
    method = "MML", model = "PCM", step_facet = "Criterion", quad_points = 7,
    maxit = 300, category_policy = "preserve")
  expect_true(is.na(attr(logLik(coarse), "df")))
  expect_true(is.na(AIC(coarse)))
  expect_match(attr(logLik(coarse), "ICReason"), "quadrature")
  j <- fits[["PCM-JML"]]
  expect_identical(nobs(j), nrow(j$prep$data))
  expect_equal(as.numeric(logLik(j)), -j$opt$value)
  expect_identical(attr(logLik(j), "likelihood_basis"), "joint")
  expect_true(is.na(AIC(j)) && is.na(BIC(j)))
  expect_match(attr(logLik(j), "ICReason"), "JML.*unsupported")
  expect_error(vcov(j), "JML point estimates")
})

test_that("qualification failures cannot be hidden by generic extraction", {
  f <- extractor_fits()[["PCM-MML-fixed"]]
  testthat::local_mocked_bindings(
    mfrm_ic_fit_check = function(...) list(eligible = FALSE, review = "Source solution failed numerical checks."),
    mfrm_facet_intervals = function(...) list(parameter_covariance = diag(length(f$opt$par)),
      checks = data.frame(Passed = FALSE, Detail = "Source solution failed numerical checks.")),
    .package = "mfrmr")
  expect_false(attr(logLik(f), "ICEligible"))
  expect_true(is.na(AIC(f)) && is.na(BIC(f)))
  expect_match(attr(logLik(f), "ICReason"), "failed numerical")
  expect_error(vcov(f), "Joint covariance is unavailable.*Source solution failed numerical checks")
  # Check the weight contract separately from numerical readiness. These
  # synthetic metadata exercise the adapter, not weighted-fit accuracy.
  weighted <- f; weighted$config$weight_col <- "Weight"
  weighted$prep$data$Weight[] <- 2
  contract <- build_mfrm_ic_contract(-weighted$opt$value, length(weighted$opt$par),
    weighted$prep, weighted$config, "MML")
  for (key in names(contract)) weighted$summary[[key]] <- contract[[key]]
  expect_identical(attr(logLik(weighted), "likelihood_basis"), "weighted_marginal_objective")
  expect_true(is.na(AIC(weighted)))
  expect_match(attr(logLik(weighted), "ICReason"), "Nonunit rating weights")
})

test_that("primary boundary estimates and special-object restrictions are preserved", {
  f <- extractor_fits()[["PCM-JML"]]
  extreme <- f; extreme$facets$person$Estimate[1] <- -Inf
  expect_error(coef(extreme), "Nonfinite primary")
  expect_identical(unname(coef(extreme, "persons"))[1], -Inf)
  expect_true(all(is.finite(extreme$opt$par)))
  adjusted <- f; adjusted$config$jml_adjustment <- TRUE
  expect_error(coef(adjusted), "corrected JML")
  expect_error(logLik(adjusted), "corrected JML")
  expect_equal(coef(adjusted, "locations"), coef(f, "locations"))
  penalized <- f; penalized$config$facet_shrinkage <- "ridge"
  expect_error(coef(penalized), "Penalized")
  expect_error(logLik(penalized), "penalized objective")
  imported <- f; class(imported) <- c("mfrm_imported_fit", class(f))
  expect_error(coef(imported), "native fit_mfrm")
  expect_error(nobs(imported), "native fit_mfrm")
})

test_that("CML extractors agree with a closed-form paired binary model", {
  skip_if_not_installed("lpSolve")
  patterns <- rbind(c(1, 0), c(0, 1), c(0, 0), c(1, 1))
  x <- patterns[rep(seq_len(4L), c(3L, 7L, 2L, 3L)), , drop = FALSE]
  d <- data.frame(Person = rep(seq_len(nrow(x)), each = 2L),
    Rater = rep(c("A", "B"), nrow(x)), Score = as.vector(t(x)))
  f <- fit_mfrm(d, "Person", "Rater", "Score", method = "CML", category_policy = "preserve")
  expect_equal(unname(coef(f)), log(7 / 3) / 2, tolerance = 1e-6)
  expect_identical(names(coef(f)), rownames(vcov(f)))
  expect_equal(as.numeric(vcov(f)), 10 / (4 * 3 * 7), tolerance = 1e-6)
  ll <- logLik(f)
  expect_equal(as.numeric(ll), 3 * log(.3) + 7 * log(.7), tolerance = 1e-10)
  expect_identical(attr(ll, "likelihood_basis"), "conditional")
  expect_identical(nobs(f), 15L)
  expect_identical(nobs(f, "ratings"), 30L)
  expect_equal(attr(ll, "informative_persons"), 10)
  expect_equal(AIC(f), -2 * as.numeric(ll) + 2)
  expect_equal(BIC(f), -2 * as.numeric(ll) + log(15))
  expect_error(coef(f, "persons"), "no saved persons")
  m <- extractor_fits()[["PCM-MML-fixed"]]
  for (generic in list(AIC, BIC)) {
    expect_error(generic(f, m), "conditional and marginal")
    expect_error(generic(m, f), "conditional and marginal")
  }
  expect_error(vcov(f, typo = TRUE), "must be empty")
})

test_that("comparison-only extraction retains unidentified targets as missing", {
  cml <- readRDS(test_path("fixtures", "cml-comparison.rds"))
  b <- coef(cml); V <- vcov(cml)
  expect_identical(names(b), rownames(V))
  expect_identical(names(b), colnames(V))
  expect_true(is.na(b[1]) && all(is.na(V[1, ])) && all(is.na(V[, 1])))
  expect_equal(unname(V[2, 2]), cml$comparisons$SE[2]^2)
  expect_error(logLik(cml), "Comparison-only")
  expect_equal(nobs(cml), length(unique(cml$prep$data$Person)))
  for (case in readRDS(test_path("fixtures", "partial-comparisons.rds"))) {
    f <- case$fit
    expect_identical(unname(coef(f)), f$comparisons$Estimate)
    expect_error(vcov(f), "descriptive estimates")
    expect_error(logLik(f), "Comparison-only")
    expect_equal(nobs(f, "ratings"), nrow(f$prep$data))
  }
})

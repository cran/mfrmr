test_that("normal shrinkage matches independently integrated Gaussian posteriors", {
  estimate <- c(-1.2, .4, 1.6)
  se <- c(.35, .8, 1.1)
  tau <- .7
  out <- mfrmr:::.compute_facet_shrinkage(estimate, se, prior_sd = tau)
  for (j in seq_along(estimate)) {
    density <- function(theta) dnorm(theta, 0, tau) * dnorm(estimate[j], theta, se[j])
    norm <- integrate(density, -Inf, Inf, rel.tol = 1e-11)$value
    mu <- integrate(function(t) t * density(t), -Inf, Inf, rel.tol = 1e-11)$value / norm
    variance <- integrate(function(t) (t - mu)^2 * density(t), -Inf, Inf,
      rel.tol = 1e-11)$value / norm
    expect_equal(out$shrunk_estimates[j], mu, tolerance = 1e-10)
    expect_equal(out$shrunk_ses[j], sqrt(variance), tolerance = 1e-10)
  }
})

test_that("rounded shrinkage factors do not imply full pooling or zero SD", {
  x <- c(-1, 0, 1)
  out <- mfrmr:::.compute_facet_shrinkage(x, rep(1, 3), prior_sd = 1e-9)
  expect_identical(out$shrinkage_factors, rep(1, 3))
  expect_equal(out$shrunk_estimates / 1e-18, x, tolerance = 1e-14)
  expect_equal(out$shrunk_ses / 1e-9, rep(1, 3), tolerance = 1e-14)
  expect_equal(out$retained_weights / 1e-18, rep(1, 3), tolerance = 1e-14)
  fit <- structure(list(config = list(method = "JML"),
    facets = list(others = data.frame(Facet = "Rater", Level = letters[1:3],
      Estimate = x, ModelSE = 1))), class = "mfrm_fit")
  adjusted <- apply_empirical_bayes_shrinkage(fit, facet_prior_sd = 1e-9)
  expect_equal(adjusted$shrinkage_report$EffectiveDF / 1e-18, 3, tolerance = 1e-14)
  expect_identical(adjusted$shrinkage_report$CalculationVersion, 2L)
  expect_identical(adjusted$facets$others$Estimate, fit$facets$others$Estimate)
})

test_that("shrinkage is stable when intermediate squares or sums overflow", {
  x <- c(-1, 0, 1)
  fn <- mfrmr:::.compute_facet_shrinkage
  small <- fn(x, rep(.1, 3))
  large <- fn(x * 1e154, rep(1e153, 3))
  expect_true(is.finite(large$tau2))
  expect_equal(large$tau2 / 1e308, small$tau2, tolerance = 1e-14)
  expect_equal(large$shrunk_estimates / 1e154, small$shrunk_estimates, tolerance = 1e-14)
  expect_equal(large$shrunk_ses / 1e154, small$shrunk_ses, tolerance = 1e-14)
  equal <- fn(x, rep(1e154, 3), prior_sd = 1e154)
  expect_equal(equal$shrinkage_factors, rep(.5, 3))
  expect_equal(equal$shrunk_estimates, x / 2)
  expect_equal(equal$shrunk_ses / 1e154, rep(1 / sqrt(2), 3))
  # The retained variance weight underflows, but the adjusted estimate and
  # posterior SD are still representable when evaluated without that square.
  extreme <- fn(x * 1e300, rep(1e150, 3), prior_sd = 1e-150)
  expect_equal(extreme$shrunk_estimates / 1e-300, x, tolerance = 1e-14)
  expect_equal(extreme$shrunk_ses / 1e-150, rep(1, 3), tolerance = 1e-14)
})

test_that("unrepresentable moment variances preserve points without invented factors", {
  x <- c(-1, 0, 1)
  fn <- mfrmr:::.compute_facet_shrinkage
  for (scale in c(1e200, 1e-200)) {
    out <- fn(x * scale, rep(scale / 10, 3))
    expect_identical(out$shrunk_estimates, x * scale)
    expect_identical(out$shrunk_ses, rep(scale / 10, 3))
    expect_true(all(is.na(out$shrinkage_factors)))
    expect_true(is.na(out$tau2))
    expect_match(out$note, "shrinkage not applied", fixed = TRUE)
  }
  for (prior in c(1e-200, 1e200)) {
    expect_error(fn(x, rep(1, 3), prior_sd = prior), "square is representable", fixed = TRUE)
  }
})

test_that("Person shrinkage does not treat EAP posterior SDs as sampling SEs", {
  fit <- structure(list(config = list(method = "MML"),
    facets = list(
      person = data.frame(Person = letters[1:3], Estimate = c(-1, 0, 1),
        SE = .5, PosteriorSD = .5),
      others = data.frame(Facet = "Rater", Level = letters[1:3],
        Estimate = c(-1, 0, 1), ModelSE = .3))), class = "mfrm_fit")
  expect_error(apply_empirical_bayes_shrinkage(fit, shrink_person = TRUE),
    "EAP posterior SDs are not sampling SEs", fixed = TRUE)
  adjusted <- apply_empirical_bayes_shrinkage(fit, shrink_person = FALSE)
  expect_identical(adjusted$facets$person, fit$facets$person)
  fit$config$method <- "JML"
  fit$config$jml_adjustment <- TRUE
  expect_error(apply_empirical_bayes_shrinkage(fit, shrink_person = TRUE), "corrected JML")

  local_mocked_bindings(mfrm_estimate = function(...) stop("Unexpected fitting"),
    .package = "mfrmr")
  data <- load_mfrmr_data("example_core")
  expect_error(fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
    method = "MML", facet_shrinkage = "empirical_bayes", shrink_person = TRUE),
    "EAP posterior SDs are not sampling SEs", fixed = TRUE)
  expect_error(fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", facet_shrinkage = "empirical_bayes", shrink_person = TRUE,
    jml_correction_order = 1), "ordinary JML")
  expect_error(fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
    facet_prior_sd = 1e-200), "square is representable", fixed = TRUE)
})

test_that("saved shrinkage reports require refreshed calculation metadata", {
  fit <- structure(list(config = list(method = "JML"),
    facets = list(others = data.frame(Facet = "Rater", Level = letters[1:3],
      Estimate = c(-1, 0, 1), ModelSE = .2))), class = "mfrm_fit")
  adjusted <- apply_empirical_bayes_shrinkage(fit)
  expect_false(shrinkage_report(adjusted)$SupportsFormalInference)
  adjusted$shrinkage_report$CalculationVersion <- NULL
  expect_error(shrinkage_report(adjusted), "Reapply apply_empirical_bayes_shrinkage", fixed = TRUE)
  expect_error(mfrmr:::.build_shrinkage_plot_data(adjusted), "no MFRM refit", fixed = TRUE)
  expect_identical(apply_empirical_bayes_shrinkage(adjusted)$shrinkage_report,
    apply_empirical_bayes_shrinkage(fit)$shrinkage_report)
})

test_that("ggplot shrinkage conversion retains levels, both points and saved bands", {
  skip_if_not_installed("ggplot2")
  fit <- structure(list(config = list(method = "JML"),
    facets = list(others = data.frame(Facet = "Rater", Level = letters[1:3],
      Estimate = c(-1, 0, 1), ModelSE = c(.2, .4, .6)))), class = "mfrm_fit")
  adjusted <- apply_empirical_bayes_shrinkage(fit, facet_prior_sd = .7)
  table <- mfrmr:::.add_shrinkage_ci_columns(
    mfrmr:::.build_shrinkage_plot_data(adjusted)$table)
  # Use stored bounds deliberately different from a recomputed normal band.
  table$ShrunkCI_Lower[1] <- -.8
  table$ShrunkCI_Upper[1] <- -.4
  table$CI_Lower[2] <- table$CI_Upper[2] <- NA_real_
  p <- mfrmr:::new_mfrm_plot_data("shrinkage",
    list(data = table, shrinkage_report = adjusted$shrinkage_report,
      show_ci = TRUE, interpretation = mfrmr:::shrinkage_interpretation()))
  g <- as_ggplot(p)
  layers <- ggplot2::ggplot_build(g)$data
  # Reference, arrows, saved bands, and both sets of points.
  expect_identical(vapply(layers, nrow, 1L), c(1L, 3L, 5L, 6L))
  expect_equal(sort(layers[[4]]$x), sort(c(table$Estimate, table$ShrunkEstimate)))
  expect_true(-.8 %in% layers[[3]]$x)
  expect_true(-.4 %in% layers[[3]]$xend)
  expect_identical(attr(g, "mfrmr_plot_data"), p)
  expect_identical(attr(as_ggplot(p, component = "data"), "mfrmr_plot_data"), p)
  expect_null(g$labels$title)
  expect_error(as_ggplot(p, component = "shrinkage_report"), "component = 'data'", fixed = TRUE)
  p$data$show_ci <- FALSE
  expect_identical(vapply(ggplot2::ggplot_build(as_ggplot(p))$data, nrow, 1L), c(1L, 3L, 6L))
  p$data$data <- table[FALSE, ]
  expect_s3_class(as_ggplot(p), "ggplot")
  p$data$shrinkage_report$CalculationVersion <- NULL
  expect_error(as_ggplot(p), "Reapply apply_empirical_bayes_shrinkage", fixed = TRUE)
})

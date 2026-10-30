# GENOVA 3.1 (genova36.exe), directly executed 2026-10-24 on deterministic
# scalar orthogonal-score inputs. References are literal printed output, not
# calculations from the implementation. The executable is not a test dependency.
genova_policy_fixture <- function(name) {
  rows <- read.csv(testthat::test_path("fixtures", "genova-component-policy-input.csv"))
  rows <- rows[rows$Case == name, ]
  mfrm_g_study(rows, "Score", rater = if (name == "one_negative") NULL else "Rater",
    nesting = if (name == "nested_negative") c(Rater = "Task") else NULL)
}

test_that("GENOVA treatments reproduce directly executed component and D-study output", {
  reference <- read.csv(testthat::test_path("fixtures", "genova-component-policy-reference.csv"))
  components <- read.csv(testthat::test_path("fixtures", "genova-component-policy-components.csv"))
  for (name in unique(reference$Case)) {
    g <- genova_policy_fixture(name)
    before <- serialize(g, NULL)
    for (policy in c("genova_algorithm", "genova_ems")) {
      expected <- reference[reference$Case == name & reference$Policy == policy, ]
      grid <- expected[if (name == "one_negative") "Tasks" else c("Raters", "Tasks")]
      d <- mfrm_d_study(g, grid, negative_components = policy)
      metrics <- c("UniverseVariance", "RelativeErrorVariance", "AbsoluteErrorVariance",
        "G", "Phi", "RelativeSEM", "AbsoluteSEM")
      # Printed GENOVA quantities have five decimal places.
      expect_lt(max(abs(as.matrix(d$coefficients[metrics]) - as.matrix(expected[metrics]))), 5.1e-6)
      effective <- components[components$Case == name & components$Policy == policy, ]
      expect_lt(max(abs(unname(unlist(d$effective_components)) - effective$Estimate)), 5.1e-6)
      expect_identical(serialize(g, NULL), before)
      expect_identical(d$gstudy, g)
      expect_true(all(d$coefficients$ComponentTreatment == policy))
      expect_true(all(d$effective_component_diagnostics$PositiveSemidefinite))
    }
  }
})

test_that("raw retention and source-count EMS remain distinct from zero substitution", {
  g <- genova_policy_fixture("nested_negative")
  raw <- mfrm_d_study(g, data.frame(Raters = 3, Tasks = 6))
  algorithm <- mfrm_d_study(g, raw$design_grid, negative_components = "genova_algorithm")
  ems <- mfrm_d_study(g, raw$design_grid, negative_components = "genova_ems")
  expect_identical(raw$effective_components, g$components)
  expect_equal(unname(unlist(g$components)), c(1.3, .075, -.1, -.6, 1.4), tolerance = 1e-10)
  expect_equal(unname(unlist(algorithm$effective_components)), c(1.3, .075, 0, 0, 1.4), tolerance = 1e-10)
  expect_equal(unname(unlist(ems$effective_components)), c(1.1, 0, 0, 0, 1.4), tolerance = 1e-10)
  expect_false(ems$coefficients$ComponentPSD)
  positive <- genova_policy_fixture("crossed_positive")
  for (policy in c("genova_algorithm", "genova_ems"))
    expect_equal(mfrm_d_study(positive, negative_components = policy)$effective_components,
      positive$components, tolerance = 1e-12)
})

test_that("GENOVA policy survives saved consumers and uses effective plot contributions", {
  g <- genova_policy_fixture("crossed_negative")
  d <- mfrm_d_study(g, data.frame(Raters = c(2, 4), Tasks = c(6, 2)),
    weights = c(Score = -2), negative_components = "genova_ems")
  restored <- unserialize(serialize(d, NULL))
  expect_identical(summary(restored), d$coefficients)
  expect_output(print(restored), "GENOVA EMS")
  p <- plot(restored, draw = FALSE)
  expect_identical(p$data$negative_components, "genova_ems")
  expect_true(any(grepl("GENOVA EMS", plot_data(p, "notes")$Text)))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(restored)
    expect_true(any(grepl("GENOVA EMS", plot_data(gg, "notes")$Text)))
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
  e <- plot(restored, type = "error_components", draw = FALSE)
  z <- e$data$table
  sums <- aggregate(z$Contribution, z[c("Scenario", "Score", "Metric")], sum)
  for (i in seq_len(nrow(sums))) {
    expected <- subset(d$coefficients, Scenario == sums$Scenario[i] & Score == sums$Score[i])
    expect_equal(sums$x[i], expected[[sums$Metric[i]]], tolerance = 1e-12)
  }
  for (field in c("effective_components", "negative_components")) {
    altered <- restored; altered[[field]] <- NULL
    expect_error(summary(altered), "treatment")
  }
  altered <- restored; altered$negative_components <- "genova_algorithm"
  expect_error(plot(altered, draw = FALSE), "do not match")
  altered <- restored; altered$negative_components <- altered$effective_components <- NULL
  expect_error(confint(altered, assumption = "normal"), "treatment is missing")
  altered <- restored; altered$effective_components[[1]][] <- 99
  expect_error(print(altered), "do not match")
  old <- mfrm_d_study(g)
  old$negative_components <- old$effective_components <- old$effective_component_diagnostics <- NULL
  old$coefficients$ComponentTreatment <- NULL
  expect_identical(summary(old), old$coefficients)
  expect_s3_class(plot(old, draw = FALSE), "mfrm_plot_data")
})

test_that("unsupported sources and mixed interval treatments fail explicitly", {
  g <- genova_policy_fixture("one_negative")
  grid <- data.frame(Tasks = c(2, 6))
  raw <- mfrm_d_study(g, grid)
  ci <- confint(raw, assumption = "normal")
  for (policy in c("genova_algorithm", "genova_ems")) {
    d <- mfrm_d_study(g, grid, negative_components = policy)
    expect_error(summary(d, intervals = ci), "treatment does not match")
    expect_error(plot(d, metric = "G", intervals = ci, draw = FALSE), "treatment does not match")
    expect_error(mfrm_multivariate_d_compare(d, assumption = "normal"), "Adjusted-estimator delta")
    data <- g$data; data$Score2 <- data$Score * 2
    multi <- mfrm_g_study(data, c("Score", "Score2"), rater = NULL)
    expect_error(mfrm_d_study(multi, negative_components = policy), "univariate ANOVA")
    incomplete <- mfrm_g_study(data[-1L, ], "Score", rater = NULL, method = "minque0")
    expect_error(mfrm_d_study(incomplete, grid, negative_components = policy), "univariate ANOVA")
    altered <- g; altered$components[[1L]][] <- 99
    expect_error(mfrm_d_study(altered, negative_components = policy), "disagree")
  }
})

test_that("adjusted points retain independent population pivots and explicit interval bases", {
  g <- genova_policy_fixture("nested_negative")
  # Source children per parent = 2: the one-facet pivot uses MS(P)/MS(P:T).
  grid <- data.frame(Raters = c(2, 2), Tasks = c(3, 6))
  raw <- mfrm_d_study(g, grid)
  d <- mfrm_d_study(g, grid, negative_components = "genova_ems")
  before <- serialize(d, NULL)
  ci <- confint(d, assumption = "normal")
  ratio <- (8/.2) / qf(c(.975, .025), 3, 6)
  independent <- t(vapply(grid$Tasks, function(k)
    (ratio - 1) / (ratio - 1 + 3/k), numeric(2)))
  expect_equal(ci$Lower, independent[, 1], tolerance = 1e-12)
  expect_equal(ci$Upper, independent[, 2], tolerance = 1e-12)
  expect_identical(ci$Estimate, d$coefficients$G)
  expect_identical(ci$IntervalEstimate, raw$coefficients$G)
  expect_true(any(ci$Estimate != ci$IntervalEstimate))
  expect_true(all(ci$ComponentTreatment == "genova_ems"))
  expect_true(all(grepl("population parameter", ci$IntervalBasis)))
  expect_identical(serialize(d, NULL), before)
  saved <- unserialize(serialize(ci, NULL))
  expect_identical(summary(d, intervals = saved)$intervals[[1]], saved)
  p <- plot(d, metric = "G", intervals = saved, draw = FALSE)
  expect_true(any(grepl("not recentered", plot_data(p, "notes")$Text)))
  expect_s3_class(plot(saved, draw = FALSE), "mfrm_plot_data")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- as_ggplot(d, metric = "G", intervals = saved)
    expect_true(any(grepl("not recentered", plot_data(gg, "notes")$Text)))
    expect_no_warning(ggplot2::ggplotGrob(gg))
  }
  for (field in c("IntervalEstimate", "IntervalBasis", "ComponentTreatment")) {
    bad <- saved; bad[[field]] <- NULL
    expect_error(summary(d, intervals = bad), "interval basis")
  }
  bad <- saved; bad$IntervalEstimate[] <- .1
  expect_error(summary(d, intervals = bad), "point estimates/statuses")
  bad <- saved; attr(bad, "negative_components") <- NULL
  expect_error(summary(d, intervals = bad), "treatment does not match")
  expect_error(summary(raw, intervals = saved), "treatment does not match")
})

test_that("GENOVA points never recenter intervals or silently alter interval availability", {
  for (name in c("one_negative", "nested_negative", "crossed_negative")) {
    g <- genova_policy_fixture(name)
    grid <- if (name == "one_negative") data.frame(Tasks = c(2, 6)) else
      data.frame(Raters = c(2, 2), Tasks = c(3, 6))
    raw <- mfrm_d_study(g, grid, weights = c(Score = -2))
    d <- mfrm_d_study(g, grid, weights = raw$weights, negative_components = "genova_ems")
    for (metric in c("G", "Phi", "RelativeSEM", "AbsoluteSEM")) {
      a <- confint(raw, parm = metric, assumption = "normal", draws = 5000L, seed = 41)
      b <- confint(d, parm = metric, assumption = "normal", draws = 5000L, seed = 41)
      for (column in c("Lower", "Upper", "Status", "Available", "Validation", "MethodCode"))
        if (is.numeric(a[[column]]))
          expect_equal(b[[column]], a[[column]], tolerance = 16 * .Machine$double.eps) else
            expect_identical(b[[column]], a[[column]])
      expect_identical(b$IntervalEstimate, a$Estimate)
      expect_identical(b$IntervalPointStatus, a$PointStatus)
      expect_identical(summary(d, intervals = b)$intervals[[1]], b)
    }
    delta <- confint(d, assumption = "normal", method = "delta")
    original <- confint(raw, assumption = "normal", method = "delta")
    expect_identical(delta$IntervalSE, original$SE)
    expect_null(delta$SE)
    expect_identical(delta$Available, original$Available)
    expect_identical(summary(d, intervals = delta)$intervals[[1]], delta)
    # Joint-region comparisons can vary both facets; the original interval
    # remains a population-difference interval, even when raw points fail.
    method <- if (name == "crossed_negative") "mean_square_region" else "exact_f"
    a <- confint(raw, assumption = "normal", method = method, reference = 1)
    b <- confint(d, assumption = "normal", method = method, reference = 1)
    expect_identical(b$Lower, a$Lower)
    expect_identical(b$Upper, a$Upper)
    expect_identical(b$IntervalEstimate, a$Estimate)
    tab <- subset(d$coefficients, Kind == "Composite")
    expect_identical(b$Estimate, tab$G[2] - tab$G[1])
    expect_identical(summary(d, intervals = b)$intervals[[1]], b)
  }
  # A separate residual projection keeps its own center and diagnostics.
  g <- genova_policy_fixture("one_negative")
  d <- mfrm_d_study(g, negative_components = "genova_algorithm")
  s <- summary(d, intervals = TRUE, assumption = "normal", draws = 5000L)
  expect_true(all(s$table$ComponentTreatment == "genova_algorithm"))
  expect_output(print(s), "IntervalEstimate")
  ci <- confint(d, parm = "RelativeSEM", assumption = "normal", method = "residual_contrast")
  expect_equal(ci$IntervalCenter, ci$IntervalEstimate, tolerance = 1e-12)
  expect_identical(summary(d, intervals = ci)$intervals[[1]], ci)
})

test_that("native GENOVA SEs describe the original estimator despite adjusted points", {
  g <- genova_policy_fixture("nested_negative")
  C <- vcov(g, assumption = "normal")
  # Literal GENOVA 3.1 output for BOTH ALGORITHM and EMS, future T=6,R:T=3:
  # universe SE=.84344, relative error-variance SE=.01856, absolute SE=.01574.
  # Zeroing components does not zero their original sampling covariances.
  relative <- c(0, 0, 0, 1/6, 1/18)
  absolute <- c(0, 1/6, 1/18, 1/6, 1/18)
  se <- c(sqrt(C[1, 1]), sqrt(drop(crossprod(relative, C %*% relative))),
    sqrt(drop(crossprod(absolute, C %*% absolute))))
  expect_lt(max(abs(se - c(.84344, .01856, .01574))), 5.1e-6)
  # Deleting truncated terms/covariances gives a different, unsupported SE.
  expect_gt(abs(sqrt(C[5, 5])/18 - se[2]), .01)
  for (policy in c("genova_algorithm", "genova_ems")) {
    d <- mfrm_d_study(g, data.frame(Raters = 3, Tasks = 6), negative_components = policy)
    rel <- confint(d, parm = "RelativeSEM", assumption = "normal")
    ci_abs <- confint(d, parm = "AbsoluteSEM", assumption = "normal")
    expect_lt(abs(rel$ErrorVarianceSE - .01856), 5.1e-6)
    expect_lt(abs(ci_abs$ErrorVarianceSE - .01574), 5.1e-6)
    # A finite adjusted point must not convert an unavailable raw delta
    # approximation into a successful interval or zero sampling uncertainty.
    delta <- confint(d, assumption = "normal", method = "delta")
    expect_true(is.finite(delta$Estimate))
    expect_true(is.na(delta$IntervalEstimate))
    expect_false(delta$Available)
    expect_true(is.na(delta$IntervalSE))
    expect_identical(summary(d, intervals = delta)$intervals[[1]], delta)
  }
})

max_points_product_arguments <- function() {
  list(data = readRDS(test_path("fixtures", "gmfrm-joint-information.rds"))$data,
    person = "Person", facets = c("Task", "Rater"), score = "Score",
    model = "GPCM", method = "MML", slope_facet = c("Task", "Rater"),
    step_facet = "Rater", noncenter_facet = "Rater", mml_engine = "em",
    mml_integration = "auto", gpcm_mml_identification = "fixed_standard_normal",
    maxit = 500L, rating_min = 0, rating_max = 2, category_policy = "preserve",
    mml_accuracy_control = list(max_points = 241L))
}

test_that("maximum controls include reference orders and retain legacy controls", {
  expect_identical(mfrmr_mml_accuracy_orders(121, offset = -1L), c(31L, 61L))
  expect_identical(mfrmr_mml_accuracy_orders(100, offset = -1L), c(31L, 50L))
  expect_identical(mfrmr_mml_accuracy_orders(121), c(31L, 60L))
  expect_identical(mfrmr_mml_accuracy_orders(7), 3L)
  for (cap in c(7L, 20L, 61L, 100L, 121L, 241L)) {
    for (offset in c(-1L, 1L)) {
      q <- mfrmr_mml_accuracy_orders(cap, offset = offset)
      expect_true(all(2*q + offset <= cap))
      expect_true(all(diff(q) > 0))
    }
  }
  old <- mfrmr_mml_accuracy_control()
  expect_identical(old$adaptive_points, c(31L, 61L, 121L))
  expect_identical(mfrmr_mml_accuracy_control(old), old)
  capped <- mfrmr_mml_accuracy_control(list(max_points = 121))
  expect_null(capped$adaptive_points)
  expect_identical(mfrmr_mml_accuracy_control(capped), capped)
  expect_error(mfrmr_mml_accuracy_control(list(max_points = 121, adaptive_points = 31)), "not both")
  for (bad in list(NULL, 0, 3.5, NA_real_, Inf, 10+1i, TRUE, c(61, 121))) {
    expect_error(mfrmr_mml_accuracy_control(list(max_points = bad)), "max_points")
    expect_error(mfrmr_gmfrm_accuracy_control(list(max_points = bad)), "max_points")
  }
  expect_error(mfrmr_mml_accuracy_orders(121, initial = 1), "quad_points")
  expect_error(mfrmr_gmfrm_accuracy_control(list(gradient_tol = 1)), "stationarity")
})

test_that("ordinary auto never evaluates a rule beyond the maximum", {
  actual <- integer()
  original <- gauss_hermite_normal
  local_mocked_bindings(gauss_hermite_normal = function(n) {
    actual <<- c(actual, as.integer(n))
    if (n > 63) stop("Point ceiling exceeded")
    original(n)
  }, .package = "mfrmr")
  args <- list(data = data.frame(Person = "P1", Rater = rep(c("A", "B"), each = 4),
    Score = c(1, 1, 1, 0, 1, 0, 0, 0), Weight = 50),
    person = "Person", facets = "Rater", score = "Score", weight = "Weight",
    repeated_ratings = "independent", anchor_policy = "silent", reltol = 1e-12,
    mml_integration = "auto", mml_accuracy_control = list(max_points = 63L))
  fit <- suppressWarnings(do.call(fit_mfrm, args))
  expect_identical(fit$config$mml_accuracy$status, "criteria_met")
  expect_identical(fit$config$mml_accuracy$history$Integration, c("fixed", "adaptive"))
  expect_true(length(actual) > 0L && max(actual) <= 63L)
  expect_true(all(fit$config$mml_accuracy$history$ReferencePoints <= 63L))
  expect_identical(fit$config$replay_inputs$mml_accuracy_control$max_points, 63L)
  args$quad_points <- 81L
  expect_error(do.call(fit_mfrm, args), "exceeds.*max_points")
})

test_that("a low ordinary ceiling reduces the default start and reports exhaustion", {
  orders <- integer()
  original <- gauss_hermite_normal
  local_mocked_bindings(gauss_hermite_normal = function(n) {
    orders <<- c(orders, as.integer(n))
    if (n > 7) stop("Point ceiling exceeded")
    original(n)
  }, .package = "mfrmr")
  error <- suppressWarnings(tryCatch(fit_mfrm(load_mfrmr_data("example_operational"),
    "Person", c("Rater", "Criterion"), "Score", mml_integration = "auto",
    mml_accuracy_control = list(max_points = 7)), error = identity))
  expect_s3_class(error, "mfrmr_mml_accuracy_error")
  expect_match(conditionMessage(error), "max_points was exhausted")
  expect_identical(error$history$QuadPoints, c(7L, 3L))
  expect_true(length(orders) > 0L && max(orders) <= 7L)
})

max_points_product_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      args <- max_points_product_arguments()
      args$quad_points <- 15L
      rules <- integer(); fitted <- integer()
      original_rule <- gauss_hermite_normal
      original_fit <- mfrm_fit_product_slopes
      local_mocked_bindings(gauss_hermite_normal = function(n) {
        rules <<- c(rules, as.integer(n))
        if (n > args$mml_accuracy_control$max_points) stop("Point ceiling exceeded")
        original_rule(n)
      }, mfrm_fit_product_slopes = function(args, supplied) {
        fitted <<- c(fitted, as.integer(args$quad_points))
        original_fit(args, supplied)
      }, .package = "mfrmr")
      fit <- suppressWarnings(do.call(fit_mfrm, args))
      cached <<- list(fit = fit, rules = rules, fitted = fitted)
    }
    cached
  }
})

test_that("two-family EM refines while preserving the engine and replay contract", {
  x <- max_points_product_fixture(); fit <- x$fit
  h <- fit$config$mml_accuracy$history
  expect_identical(fit$config$mml_accuracy$status, "criteria_met")
  expect_gt(nrow(h), 1L)
  expect_identical(tail(h$Status, 1), "criteria_met")
  expect_true(all(head(h$Status, -1) == "refinement_required"))
  expect_true(all(h$Integration == "fixed"))
  expect_true(all(h$ReferencePoints <= 241))
  expect_identical(x$fitted, h$QuadPoints)
  expect_equal(max(x$rules), max(h$ReferencePoints))
  expect_identical(fit$config$estimation_control$mml_engine_used, "em")
  expect_identical(fit$config$replay_inputs$mml_integration, "auto")
  expect_identical(fit$config$replay_inputs$quad_points, 15L)
  expect_identical(fit$summary$MMLAccuracyStatus, "criteria_met")
  expect_output(print(summary(fit)), "Automatic integration checks: criteria met")
  expect_true(all(fit$config$mml_accuracy$checks$Passed))
  expect_equal(nrow(fit$config$mml_accuracy$per_person), fit$prep$n_person)
  refit <- mfrmr_gqs_refit_arguments(fit, max_points_product_arguments()$data, 31L)
  expect_identical(refit$mml_integration, "fixed")
  expect_null(refit$mml_accuracy_control)
  expect_identical(refit$em_score_tol, 1e-6)
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(restored, fit)
  expect_identical(attr(suppressWarnings(confint(restored, parm = "slopes")),
    "checks"), fit$config$mml_accuracy$checks)
})

test_that("a point cap and a failed check retain distinct candidate histories", {
  args <- max_points_product_arguments()
  args$quad_points <- 15L
  args$mml_accuracy_control$max_points <- 29L
  cap <- suppressWarnings(tryCatch(do.call(fit_mfrm, args), error = identity))
  expect_s3_class(cap, "mfrmr_mml_accuracy_error")
  expect_match(conditionMessage(cap), "max_points was exhausted")
  expect_identical(cap$history$Status, "refinement_required")
  expect_identical(cap$history$ReferencePoints, 29L)
  expect_identical(cap$last_quad_points, 15L)
  expect_identical(cap$last_candidate$parameters, cap$last_parameters)
  expect_identical(cap$last_candidate$config$estimation_control$quad_points, 15L)
  expect_true(is.data.frame(cap$last_candidate$accuracy$checks))
  expect_identical(unserialize(serialize(cap, NULL)), cap)
  args$mml_accuracy_control$max_points <- 121L
  args$maxit <- 1L
  refused <- suppressWarnings(tryCatch(do.call(fit_mfrm, args), error = identity))
  expect_s3_class(refused, "mfrmr_mml_accuracy_error")
  expect_identical(refused$history$Status, "check_refused")
  expect_equal(nrow(refused$history), 1L)
  expect_match(refused$history$Detail, "EM must meet")
  expect_false(inherits(refused, "mfrm_fit"))
})

test_that("adaptive two-family auto preserves explicit initialization policy", {
  args <- max_points_product_arguments()
  args$mml_engine <- "direct"
  args$gpcm_mml_start <- "neutral"
  actual <- integer()
  original <- gauss_hermite_normal
  local_mocked_bindings(gauss_hermite_normal = function(n) {
    actual <<- c(actual, as.integer(n))
    if (n > 241) stop("Point ceiling exceeded")
    original(n)
  }, .package = "mfrmr")
  fit <- suppressWarnings(do.call(fit_mfrm, args))
  expect_identical(fit$config$estimation_control$mml_integration, "adaptive")
  expect_true(all(fit$config$mml_accuracy$history$Integration == "adaptive"))
  expect_identical(fit$opt$mml_initialization$policy, "neutral")
  expect_identical(fit$config$replay_inputs$gpcm_mml_start, "neutral")
  refit <- mfrmr_gqs_refit_arguments(fit, args$data, 31L)
  expect_identical(refit$mml_integration, "adaptive")
  expect_null(refit$mml_accuracy_control)
  expect_identical(refit$gpcm_mml_start, "neutral")
  expect_equal(max(actual), max(fit$config$mml_accuracy$history$ReferencePoints))
  expect_equal(nrow(fit$config$mml_accuracy$history), 1L)
})

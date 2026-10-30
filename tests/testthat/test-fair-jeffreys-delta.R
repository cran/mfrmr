jeffreys_delta_case <- function(model = "RSM", scope = "facet") {
  d <- load_mfrmr_data("example_core")
  d <- d[d$Person %in% unique(d$Person)[1:6], ]
  args <- list(person = "Person", facets = c("Rater", "Criterion"), score = "Score",
    method = "JML", jml_penalty = "jeffreys", rating_min = 1, rating_max = 4,
    category_policy = "preserve", model = model,
    step_facet = if (model == "PCM") "Criterion" else NULL, maxit = 2000,
    centering_scope = scope)
  cal <- do.call(fit_mfrm, c(list(data = d), args))
  anchors <- data.frame(Facet = "Rater", Level = "R01",
    Anchor = subset(coef(cal)$facets, Facet == "Rater" & Level == "R01")$Estimate)
  scored <- do.call(fit_mfrm, c(list(data = d, anchors = anchors), args))
  plan <- unique(d[c("Rater", "Criterion")]); plan$Weight <- seq_len(nrow(plan))
  list(data = d, args = args, calibration = cal, fit = scored, anchors = anchors, plan = plan)
}

test_that("penalized delta recovers the exact binomial score variance", {
  n <- 6L
  intervals <- matrix(NA_real_, n + 1L, 2L)
  for (s in 0:n) {
    d <- data.frame(Person = "P1", Rater = "R1", Score = c(rep(1, s), rep(0, n - s)))
    f <- fit_mfrm(d, "Person", "Rater", "Score", method = "JML", jml_penalty = "jeffreys",
      rating_min = 0, rating_max = 1, category_policy = "preserve", repeated_ratings = "independent")
    fair <- fair_average_table(f, reference_plan = data.frame(Rater = "R1"))
    set.seed(47); rng <- .Random.seed
    ci <- confint(fair, method = "delta")
    expect_identical(.Random.seed, rng)
    p <- (s + .5) / (n + 1)
    # S* ~ Binomial(n,p): exact variance of (S*+.5)/(n+1), not inverse Hessian.
    exact <- n * p * (1 - p) / (n + 1)^2
    expect_true(ci$Available)
    expect_false(ci$CIEligible)
    expect_equal(ci$SE^2, exact, tolerance = 1e-8)
    expect_gt(abs(ci$SE^2 - p * (1 - p) / (n + 1)), 1e-4)
    expect_identical(attr(ci, "settings")$estimator, "jeffreys_canonical_joint")
    expect_match(ci$InferenceReview, "not inverse penalized curvature")
    intervals[s + 1L, ] <- c(ci$CI_Lower, ci$CI_Upper)
  }
  # Exhaustive coverage under the true law, not the fitted generating law:
  # an exact variance identity alone does not qualify normal-reference limits.
  coverage <- function(p) sum(dbinom(0:n, n, p)[intervals[, 1] <= p & p <= intervals[, 2]])
  expect_equal(coverage(.5), .78125, tolerance = 1e-14)
  expect_equal(coverage(.25), .784423828125, tolerance = 1e-14)
})

test_that("analytic penalized curvature and anchor derivatives match numerical differentiation", {
  for (model in c("RSM", "PCM")) for (scope in c("facet", "remainder", "anchor_origin")) {
    x <- jeffreys_delta_case(model, scope); f <- x$fit
    d <- mfrm_fair_jeffreys_design(f); design <- mfrm_fair_delta_design(f)
    A <- mfrm_fair_delta_anchor_map(x$calibration, x$anchors)
    p <- fitted(f)
    J <- mfrm_fair_jeffreys_jacobian(f, design, p, x$anchors, A)
    objective <- mfrm_fair_jeffreys_objective(d$category, d$offset, d$y)
    H <- if (scope == "facet") numDeriv::hessian(objective$fn, f$opt$par) else
      numDeriv::jacobian(objective$gr, f$opt$par)
    expect_equal(unname(as.matrix(J$information)), unname(H), tolerance = 2e-6)
    # Rebuild the actual constraints/offsets for each anchor perturbation.
    B <- numDeriv::jacobian(function(a) {
      changed <- f; anch <- x$anchors; anch$Anchor <- a
      spec <- prepare_constraint_specs(f$prep, anchor_df = anch,
        noncenter_facet = f$config$noncenter_facet, dummy_facets = f$config$dummy_facets,
        centering_scope = scope)
      changed$config$theta_spec <- spec$theta_spec; changed$config$facet_specs <- spec$facet_specs
      dz <- mfrm_fair_jeffreys_design(changed)
      mfrm_fair_jeffreys_objective(dz$category, dz$offset, dz$y)$gr(f$opt$par)
    }, x$anchors$Anchor)
    expect_equal(unname(J$anchors), unname(B), tolerance = 2e-6)
    expect_lt(J$diagnostics$symmetry_residual, 1e-10)
    # The penalty's anchor derivative is material in this fixture.
    Z <- mfrm_fair_delta_center(design$category, p)
    O <- mfrm_fair_delta_anchor_design(f, design$idx, x$anchors, A)
    ordinary_B <- mfrm_fair_delta_information(Z, p, mfrm_fair_delta_center(O, p))
    expect_gt(max(abs(J$anchors - ordinary_B)), 1e-4)
  }
})

test_that("penalized two-stage covariance retains shared events and paired contrasts", {
  for (model in c("RSM", "PCM")) {
    x <- jeffreys_delta_case(model); f <- x$fit; cfit <- x$calibration
    persons <- c("P003", "P004"); n <- nrow(x$data)
    maps <- list(calibration = seq_len(n), scoring = seq_len(n))
    fair <- fair_average_table(f, reference_plan = x$plan)
    ci <- confint(fair, parm = persons, method = "delta", anchor_source = cfit, response_map = maps)
    expect_true(all(ci$Available)); expect_false(any(ci$CIEligible))
    parts <- attr(ci, "covariance_components"); V <- attr(ci, "covariance")
    expect_equal(Reduce(`+`, parts), V, tolerance = 1e-12)
    expect_gt(max(abs(parts$shared_events)), 1e-6)
    L <- matrix(c(1, -1), 1, dimnames = list("difference", persons))
    contrast <- confint(fair, parm = persons, contrasts = L, method = "delta",
      anchor_source = cfit, response_map = maps)
    expect_equal(attr(contrast, "covariance"), L %*% V %*% t(L), tolerance = 1e-12)
    expect_equal(contrast$Estimate, ci$Estimate[1] - ci$Estimate[2], tolerance = 1e-12)
    id <- diag(2); dimnames(id) <- list(persons, persons)
    independent <- mfrm_fair_delta_covariance(f, cfit, fair$reference_plan, persons, id,
      x$anchors, list(calibration = seq_len(n), scoring = n + seq_len(n)),
      rbind(fitted(cfit), fitted(cfit)))
    expect_equal(unname(independent$components$shared_events), matrix(0, 2, 2))
    expect_equal(independent$covariance, parts$direct_scoring + parts$through_calibration, tolerance = 1e-12)
    # Perturb the unrounded anchor, then refit through the public fitting API.
    h <- .01
    targets <- lapply(c(-h, h), function(delta) {
      anchors <- x$anchors; anchors$Anchor <- anchors$Anchor + delta
      fit <- do.call(fit_mfrm, c(list(data = x$data, anchors = anchors), x$args))
      t <- fair_average_table(fit, reference_plan = x$plan)$reference_scores
      t$Estimate[match(persons, t$Person)]
    })
    empirical <- (targets[[2]] - targets[[1]]) / (2 * h)
    expect_lt(max(abs(unname(attr(ci, "linearization")$anchor_coupling[, 1]) - empirical)), 2e-4)
    path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
    saveRDS(ci, path); expect_identical(readRDS(path), ci)
    expect_equal(plot(ci, draw = FALSE)$data$data$Estimate, ci$Estimate)
    expect_match(plot(ci, draw = FALSE)$data$notes$Text, "Analytic Jeffreys")
    expect_output(print(ci), "Analytic Jeffreys two-stage")
  }
})

test_that("penalized derivative failures retain points and unavailable uncertainty", {
  x <- jeffreys_delta_case()
  fair <- fair_average_table(x$calibration, reference_plan = x$plan)
  local_mocked_bindings(mfrm_fair_jeffreys_jacobian = function(...) stop("Controlled penalized derivative failure"),
    .package = "mfrmr")
  ci <- confint(fair, parm = "P001", method = "delta")
  expect_false(ci$Available); expect_false(ci$CIEligible)
  expect_true(is.finite(ci$Estimate)); expect_true(is.na(ci$SE))
  expect_true(all(is.na(attr(ci, "covariance"))))
  expect_match(ci$Reason, "Controlled penalized derivative failure")
})

test_that("fully anchored scoring propagates the exact enumerated calibration variance", {
  n <- 6L; points <- numeric(n + 1L)
  args <- list(person = "Person", facets = "Rater", score = "Score", method = "JML",
    jml_penalty = "jeffreys", rating_min = 0, rating_max = 1,
    category_policy = "preserve", repeated_ratings = "independent")
  for (s in 0:n) {
    d <- data.frame(Person = "P1", Rater = "R1", Score = c(rep(1, s), rep(0, n - s)))
    cal <- do.call(fit_mfrm, c(list(data = d), args))
    anchor <- data.frame(Facet = "Person", Level = "P1", Anchor = coef(cal)$persons$Estimate)
    score <- do.call(fit_mfrm, c(list(data = d, anchors = anchor), args))
    fair <- fair_average_table(score, reference_plan = data.frame(Rater = "R1"))
    points[s + 1L] <- fair$reference_scores$Estimate
    ci <- confint(fair, method = "delta", anchor_source = cal,
      response_map = list(calibration = seq_len(n), scoring = seq_len(n)))
    expect_length(score$opt$par, 0)
    expect_true(ci$Available)
    expect_equal(attr(ci, "covariance_components")$direct_scoring[1, 1], 0)
    if (s == 2L) observed <- ci
  }
  p <- points[3L]; probability <- dbinom(0:n, n, p)
  mean <- sum(probability * points)
  exact <- sum(probability * (points - mean)^2)
  expect_equal(observed$SE^2, exact, tolerance = 1e-8)
  # This exact variance identity does not imply unbiasedness or normal coverage.
  expect_equal(mean - p, (.5 - p) / (n + 1), tolerance = 1e-7)
  # Exact conditional refit quantiles from the exhaustive PUBLIC fit points.
  # Eliminating inner Monte Carlo error does not fix the basic interval.
  lower <- upper <- numeric(n + 1L)
  for (s in 0:n) {
    cdf <- cumsum(dbinom(0:n, n, points[s + 1L]))
    lower[s + 1L] <- 2 * points[s + 1L] - points[which(cdf >= .975)[1L]]
    upper[s + 1L] <- 2 * points[s + 1L] - points[which(cdf >= .025)[1L]]
  }
  coverage <- function(p) sum(dbinom(0:n, n, p)[lower <= p & p <= upper])
  expect_equal(coverage(.5), .78125, tolerance = 1e-14)
  expect_equal(coverage(.1), .452709, tolerance = 1e-14)
  expect_equal(upper[1], points[1], tolerance = 1e-12)
})

test_that("full two-stage influence agrees with independent tilted-equation refits", {
  x <- jeffreys_delta_case("PCM"); f <- x$fit; cal <- x$calibration
  n <- nrow(x$data); persons <- c("P003", "P004")
  fair <- fair_average_table(f, reference_plan = x$plan)
  ci <- confint(fair, parm = persons, method = "delta", anchor_source = cal,
    response_map = list(calibration = seq_len(n), scoring = seq_len(n)))
  lin <- attr(ci, "linearization")
  dc <- mfrm_fair_jeffreys_design(cal); ds <- mfrm_fair_jeffreys_design(f)
  row <- which(x$data$Person == persons[1])[1]
  vc <- as.numeric(dc$category[[4]][row, ] - dc$category[[1]][row, ])
  vs <- as.numeric(ds$category[[4]][row, ] - ds$category[[1]][row, ])
  predicted <- drop(vc %*% lin$calibration_loading + vs %*% lin$scoring_loading)
  # Continuous sufficient-statistic perturbation is for derivative verification,
  # not a claim that fractional ratings are accepted by the public fitter.
  solve_tilt <- function(d, start, v, h) {
    obj <- mfrm_fair_jeffreys_objective(d$category, d$offset, d$y)
    result <- nlminb(start, function(b) obj$fn(b) - h * sum(v * b),
      function(b) obj$gr(b) - h * v,
      control = list(iter.max = 1000, eval.max = 2000, rel.tol = 1e-14, x.tol = 1e-12))
    # nlminb's objective-change stop can precede equation stationarity. Finish
    # the reference with numerical derivatives, independently of the new code.
    for (i in 1:3) {
      gradient <- obj$gr(result$par) - h * v
      if (max(abs(gradient)) < 1e-8) break
      H <- numDeriv::jacobian(obj$gr, result$par, method = "simple", method.args = list(eps = 1e-5))
      result$par <- result$par - solve(H, gradient)
    }
    expect_lt(max(abs(obj$gr(result$par) - h * v)), 1e-8)
    result$par
  }
  h <- .003
  values <- lapply(c(-h, h), function(delta) {
    bc <- solve_tilt(dc, cal$opt$par, vc, delta)
    a <- x$anchors
    pa <- expand_params(bc, build_param_sizes(cal$config), cal$config)
    a$Anchor <- pa$facets$Rater[match(a$Level, cal$config$facet_levels$Rater)]
    changed <- f
    spec <- prepare_constraint_specs(f$prep, anchor_df = a,
      noncenter_facet = f$config$noncenter_facet, dummy_facets = f$config$dummy_facets,
      centering_scope = f$config$centering_scope)
    changed$config$theta_spec <- spec$theta_spec; changed$config$facet_specs <- spec$facet_specs
    d <- mfrm_fair_jeffreys_design(changed)
    b <- solve_tilt(d, f$opt$par, vs, delta)
    mfrm_fair_jeffreys_target(changed, b, fair$reference_plan, persons)
  })
  expect_lt(max(abs((values[[2]] - values[[1]]) / (2 * h) - predicted)), 2e-4)
})

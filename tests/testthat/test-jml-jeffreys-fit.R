jeffreys_public_case <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) {
      d <- load_mfrmr_data("example_core")
      d <- d[d$Person %in% unique(d$Person)[1:8], ]
      args <- list(person = "Person", facets = c("Rater", "Criterion"), score = "Score",
        method = "JML", jml_penalty = "jeffreys", rating_min = 1, rating_max = 4,
        category_policy = "preserve", maxit = 1000)
      fit <- do.call(fit_mfrm, c(list(data = d), args))
      saved <<- list(data = d, args = args, fit = fit,
        plan = unique(d[c("Rater", "Criterion")]))
    }
    saved
  }
})

test_that("public Jeffreys fits reproduce all binomial count closed forms", {
  n <- 6L
  for (s in 0:n) {
    d <- data.frame(Person = "P1", Rater = "R1", Score = c(rep(1, s), rep(0, n - s)))
    f <- fit_mfrm(d, "Person", "Rater", "Score", method = "JML", jml_penalty = "jeffreys",
      rating_min = 0, rating_max = 1, category_policy = "preserve", repeated_ratings = "independent")
    p <- (s + .5) / (n + 1)
    expect_s3_class(f, "mfrm_jeffreys_fit")
    expect_false(inherits(f, "mfrm_fit"))
    expect_equal(f$facets$person$Estimate, qlogis(p), tolerance = 1e-6)
    expect_equal(unname(fitted(f)[, "1"]), rep(p, n), tolerance = 1e-7)
    expect_equal(f$checks$ordinary_nll, -s * log(p) - (n - s) * log1p(-p), tolerance = 1e-9)
    expect_lt(f$checks$max_gradient, 1e-4)
    fair <- fair_average_table(f, reference_plan = data.frame(Rater = "R1"))
    expect_equal(fair$reference_scores$Estimate, p, tolerance = 1e-7)
  }
  # Ordinary default retains its documented guard for a constant calibration.
  expect_error(fit_mfrm(d, "Person", "Rater", "Score", method = "JML",
    rating_min = 0, rating_max = 1, category_policy = "preserve", repeated_ratings = "independent"),
    "one score category")
})

test_that("absent PCM categories and extreme Persons are fitted on declared support", {
  x <- jeffreys_public_case(); d <- x$data
  d$Score <- ifelse(d$Score > 2, 4L, 1L)
  d$Score[d$Person == "P001"] <- 1L
  args <- x$args; args$model <- "PCM"; args$step_facet <- "Criterion"
  args$anchors <- data.frame(Facet = "Rater", Level = "R01", Anchor = .2)
  f <- do.call(fit_mfrm, c(list(data = d), args))
  expect_true(all(is.finite(c(f$opt$par, f$steps$Estimate, f$facets$person$Estimate))))
  expect_equal(f$prep$score_map$OriginalScore, 1:4)
  expect_equal(f$prep$data$Score, d$Score)
  expect_equal(unname(subset(f$facets$others, Facet == "Rater" & Level == "R01")$Estimate), .2)
  expect_equal(sum(subset(f$facets$others, Facet == "Rater")$Estimate), 0, tolerance = 1e-12)
  expect_true(all(fitted(f) > 0 & fitted(f) < 1))
  design <- mfrm_fair_jeffreys_design(f)
  objective <- mfrm_fair_jeffreys_objective(design$category, design$offset, design$y)
  numerical <- numDeriv::grad(objective$fn, f$opt$par)
  expect_lt(max(abs(numerical)), 1e-4)
  expect_equal(numerical, objective$gr(f$opt$par), tolerance = 1e-6)
  # Independent curvature check at this stationary point, not the unpenalized I.
  H <- numDeriv::hessian(objective$fn, f$opt$par)
  expect_gt(min(eigen(H, symmetric = TRUE, only.values = TRUE)$values), 0)
  idx <- build_indices(f$prep, f$config$step_facet)
  expect_equal(f$checks$ordinary_nll,
    mfrm_loglik_jml(f$opt$par, idx, f$config, build_param_sizes(f$config)), tolerance = 1e-10)
  ordinary <- args; ordinary$jml_penalty <- "none"
  expect_error(do.call(fit_mfrm, c(list(data = d), ordinary)), "unsupported internal-category")
})

test_that("Jeffreys contracts refuse unsupported uncertainty and unidentified designs", {
  x <- jeffreys_public_case(); f <- x$fit
  expect_output(print(f), "Jeffreys-penalized")
  expect_equal(coef(f)$steps, f$steps)
  expect_equal(summary(f)$persons, f$facets$person)
  expect_error(vcov(f), "covariance.*unavailable")
  expect_error(confint(f), "intervals.*unavailable")
  expect_error(logLik(f), "not a maximized ordinary likelihood")
  expect_error(AIC(f), "not a maximized ordinary likelihood")
  expect_error(fair_average_table(f), "reference_plan")
  expect_error(fair_average_table(f, reference_plan = x$plan, xtreme = .3), "FairM/FairZ")
  fair <- fair_average_table(f, reference_plan = x$plan)
  expect_output(print(summary(fair)), "Common-plan Jeffreys")
  expect_equal(summary(fair, include_person = TRUE)$reference_scores, fair$reference_scores)
  delta <- confint(fair, parm = "P001", method = "delta")
  expect_true(delta$Available)
  expect_false(delta$CIEligible)
  expect_identical(attr(delta, "settings")$estimator, "jeffreys_canonical_joint")
  expect_error(confint(fair, nsim = 2, seed = 1, generator_method = "ordinary"), "requires generator_method")
  bounded <- confint(fair, parm = "P001", method = "bounded")
  expect_true(bounded$CIEligible)
  expect_identical(attr(bounded, "settings")$estimator, "jeffreys_canonical_joint")
  changed <- fair; changed$reference_fit$opt$par[1] <- f$opt$par[1] + .1
  expect_error(confint(changed, nsim = 2, seed = 1), "retained reference fit changed")
  check_error <- function(change, pattern, data = x$data) {
    args <- utils::modifyList(x$args, change, keep.null = TRUE)
    expect_error(do.call(fit_mfrm, c(list(data = data), args)), pattern)
  }
  check_error(list(method = "MML"), "requires additive")
  check_error(list(model = "GPCM"), "requires additive")
  check_error(list(facet_interactions = "Rater:Criterion"), "does not accept")
  check_error(list(jml_correction_order = 1), "does not accept")
  check_error(list(facet_shrinkage = "empirical_bayes"), "does not accept")
  check_error(list(attach_diagnostics = TRUE), "attach_diagnostics")
  check_error(list(rating_min = NULL), "declared rubric")
  check_error(list(category_policy = "collapse"), "category_policy")
  check_error(list(anchors = data.frame(Facet = "Rater", Level = "bad", Anchor = 0)), "valid.*anchors")
  weighted <- x$data; weighted$w <- 2
  check_error(list(weight = "w"), "unit weights", weighted)
  disconnected <- x$data
  disconnected <- disconnected[(disconnected$Person %in% c("P001", "P002", "P003", "P004")) ==
    (disconnected$Rater %in% c("R01", "R02")), ]
  check_error(list(), "identified.*rank-stable", disconnected)
})

test_that("overlapping calibration, anchors, points and refits use one estimator", {
  x <- jeffreys_public_case(); d <- x$data
  maps <- list(calibration = which(d$Person %in% paste0("P00", 1:6)),
    scoring = which(d$Person %in% paste0("P00", 3:8)))
  cal <- do.call(fit_mfrm, c(list(data = d[maps$calibration, ]), x$args))
  anchors <- data.frame(Facet = "Rater", Level = "R01",
    Anchor = subset(coef(cal)$facets, Facet == "Rater" & Level == "R01")$Estimate)
  args <- x$args; args$anchors <- anchors
  score <- do.call(fit_mfrm, c(list(data = d[maps$scoring, ]), args))
  fair <- fair_average_table(score, reference_plan = x$plan)
  persons <- c("P004", "P005"); L <- matrix(c(1, -1), 1, dimnames = list("difference", persons))
  set.seed(81); rng <- .Random.seed
  ci <- confint(fair, parm = persons, contrasts = L, nsim = 2, seed = 82,
    anchor_source = cal, generating_fit = x$fit, response_map = maps)
  expect_identical(.Random.seed, rng)
  expect_identical(attr(ci, "settings")$estimator, "jeffreys_canonical_joint")
  expect_identical(attr(ci, "settings")$generator_method, "jeffreys")
  expect_match(ci$InferenceReview, "exact refit quantiles")
  expect_equal(attr(ci, "settings")$shared_responses, length(intersect(maps$calibration, maps$scoring)))
  point <- fair$reference_scores$Estimate[match(persons, fair$reference_scores$Person)]
  expect_equal(ci$Estimate, point[1] - point[2])
  expect_equal(ci$Resolved, 2)
  delta <- confint(fair, parm = persons, contrasts = L, method = "delta",
    anchor_source = cal, generating_fit = x$fit, response_map = maps)
  expect_true(delta$Available)
  expect_equal(delta$Estimate, ci$Estimate)
  expect_equal(attr(delta, "settings")$shared_responses,
    length(intersect(maps$calibration, maps$scoring)))
  expect_equal(Reduce(`+`, attr(delta, "covariance_components")), attr(delta, "covariance"), tolerance = 1e-12)
  expect_false(ci$CIEligible)
  expect_output(print(ci), "Jeffreys estimation throughout")
  expect_match(plot(ci, draw = FALSE)$data$notes$Text, "Jeffreys estimation throughout")
  # Replay draw 1 via public fit_mfrm and direct table extraction; one uniform
  # per union event, with the SAME score supplied to both overlapping stages.
  p <- fitted(x$fit); set.seed(82)
  y <- pmin(rowSums(runif(nrow(p)) > t(apply(p, 1, cumsum))), ncol(p) - 1L) + 1L
  dc <- d[maps$calibration, ]; dc$Score <- y[maps$calibration]
  ds <- d[maps$scoring, ]; ds$Score <- y[maps$scoring]
  cal1 <- do.call(fit_mfrm, c(list(data = dc), x$args))
  a1 <- anchors; a1$Anchor <- subset(coef(cal1)$facets, Facet == "Rater" & Level == "R01")$Estimate
  args$anchors <- a1
  score1 <- do.call(fit_mfrm, c(list(data = ds), args))
  ref <- fair_average_table(score1, reference_plan = x$plan)$reference_scores
  expect_equal(unname(attr(ci, "person_draws")[1, ]), ref$Estimate[match(persons, ref$Person)], tolerance = 1e-10)
  expect_equal(unname(attr(ci, "anchor_draws")[1, ]), a1$Anchor, tolerance = 1e-10)
  roots <- attr(ci, "draws")[, 1] - attr(ci, "generating_target")
  expect_equal(ci$CI_Lower, unname(ci$Estimate - quantile(roots, .975, type = 1)), tolerance = 1e-12)
  path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
  saveRDS(ci, path); expect_identical(readRDS(path), ci)
  ordinary_args <- x$args; ordinary_args$jml_penalty <- "none"
  ordinary <- do.call(fit_mfrm, c(list(data = d), ordinary_args))
  expect_error(confint(fair, parm = persons, nsim = 2, seed = 82,
    anchor_source = cal, generating_fit = ordinary, response_map = maps), "same estimator")
  actual <- mfrm_fair_refit; calls <- 0L
  local_mocked_bindings(mfrm_fair_refit = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop("Controlled penalized calibration failure")
    actual(...)
  }, .package = "mfrmr")
  failed <- confint(fair, parm = persons, contrasts = L, nsim = 2, seed = 82,
    anchor_source = cal, generating_fit = x$fit, response_map = maps)
  expect_equal(failed$Unresolved, 1L)
  expect_equal(failed$Planned, 2L)
  expect_false(failed$Available)
  expect_match(attr(failed, "trials")$Reason[1], "Controlled penalized")
})

test_that("penalized finishing preserves the objective and unchanged stationarity gate", {
  x <- jeffreys_public_case(); f <- x$fit
  d <- mfrm_fair_jeffreys_design(f)
  objective <- mfrm_fair_jeffreys_objective(d$category, d$offset, d$y)
  unchanged <- mfrm_fair_jeffreys_polish(f, f$opt$par, objective, 1e-4)
  expect_identical(unchanged$par, f$opt$par)
  expect_false(unchanged$diagnostics$attempted)
  par <- f$opt$par; par[1] <- par[1] + .002
  before <- serialize(f, NULL)
  finished <- mfrm_fair_jeffreys_polish(f, par, objective, 1e-4)
  expect_identical(serialize(f, NULL), before)
  expect_true(finished$diagnostics$attempted)
  expect_gte(finished$diagnostics$accepted_steps, 1L)
  expect_lte(finished$diagnostics$accepted_steps, 3L)
  expect_lt(finished$diagnostics$final_gradient, 1e-4)
  expect_lte(objective$fn(finished$par), objective$fn(par))
  expect_lt(max(abs(numDeriv::grad(objective$fn, finished$par))), 1e-4)
  expect_lte(finished$diagnostics$max_step, 1)
  # A failed/indefinite Jacobian must not become an accepted approximate fit.
  local_mocked_bindings(mfrm_fair_jeffreys_jacobian = function(...) stop("Controlled curvature failure"),
    .package = "mfrmr")
  rejected <- mfrm_fair_jeffreys_polish(f, par, objective, 1e-4)
  expect_identical(rejected$par, par)
  expect_equal(rejected$diagnostics$accepted_steps, 0L)
  expect_gt(rejected$diagnostics$final_gradient, 1e-4)
  expect_match(rejected$diagnostics$reason, "Controlled curvature failure")
})

fair_reference_case <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      d <- load_mfrmr_data("example_core")
      fit <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
        method = "JML", optimizer = "BFGS", maxit = 500, reltol = 1e-11)
      plan <- unique(d[c("Rater", "Criterion")])
      plan$Weight <- seq_len(nrow(plan))
      fair <- fair_average_table(fit, reference_plan = plan)
      a <- data.frame(Facet = "Rater", Level = "R01",
        Anchor = subset(fit$facets$others, Facet == "Rater" & Level == "R01")$Estimate)
      args <- fit$config$replay_inputs; args$package_version <- NULL
      args$data <- d; args$anchors <- a
      anchored <- do.call(fit_mfrm, args)
      cached <<- list(data = d, fit = fit, plan = plan, fair = fair, anchors = a,
        anchored = fair_average_table(anchored, reference_plan = plan))
    }
    cached
  }
})

test_that("common-plan scores match independent polynomials and preserve legacy targets", {
  case <- fair_reference_case(); fit <- case$fit
  before <- serialize(fit, NULL)
  legacy <- fair_average_table(fit)
  fair <- case$fair
  expect_equal(fair$raw_by_facet, legacy$raw_by_facet)
  expect_equal(fair$by_facet, legacy$by_facet)
  pa <- mfrmr:::expand_params(fit$opt$par, mfrmr:::build_param_sizes(fit$config), fit$config)
  plan <- case$plan
  independent <- function(theta) {
    eta <- theta - pa$facets$Rater[plan$Rater] - pa$facets$Criterion[plan$Criterion]
    logits <- outer(eta, 0:3) - matrix(c(0, cumsum(pa$steps)), nrow(plan), 4, byrow = TRUE)
    p <- exp(logits - apply(logits, 1, max)); p <- p / rowSums(p)
    ey <- drop(p %*% (0:3))
    c(mean = 1 + weighted.mean(ey, plan$Weight),
      derivative = weighted.mean(drop(p %*% (0:3)^2) - ey^2, plan$Weight))
  }
  ref <- vapply(pa$theta, function(t) independent(t)[1], 0)
  expect_equal(fair$reference_scores$Estimate, unname(ref), tolerance = 1e-12)
  h <- 1e-5; t <- pa$theta[1]
  fd <- (independent(t+h)[1] - independent(t-h)[1])/(2*h)
  expect_equal(unname(fd), unname(independent(t)[2]), tolerance = 1e-9)
  # With a complete equal-weight administered plan, the JML score equation
  # independently equates each expected average to that Person's observed mean.
  equal <- fair_average_table(fit, reference_plan = plan[c("Rater", "Criterion")])
  observed <- tapply(case$data$Score, case$data$Person, mean)
  expect_equal(equal$reference_scores$Estimate, as.numeric(observed[equal$reference_scores$Person]), tolerance = 5e-6)
  expect_equal(summary(fair, include_person = TRUE)$reference_scores, fair$reference_scores)
  expect_identical(serialize(fit, NULL), before)
  expect_error(fair_average_table(fit, reference_plan = plan["Rater"]), "every fitted facet")
  bad <- plan; bad$Weight[1] <- 0
  expect_error(fair_average_table(fit, reference_plan = bad), "positive")
  expect_error(fair_average_table(fit, reference_plan = rbind(plan, plan[1, ])), "duplicate")
  expect_error(confint(legacy, nsim = 2, seed = 1), "reference_plan")
  expect_error(confint(fair, nsim = 2, seed = 1, mc_level = NA_real_), "mc_level")
  expect_error(confint(fair, method = "delta", mc_level = .95), "only")
  expect_error(confint(fair, method = "bounded", mc_level = .95), "only")
  changed <- fair; changed$reference_fit$opt$par[1] <- changed$reference_fit$opt$par[1] + .1
  expect_error(confint(changed, nsim = 2, seed = 1), "retained reference fit changed")
})

test_that("common PCM exposure uses fitted profiles and is origin invariant", {
  d <- read.csv(test_path("fixtures", "tam-jml-0.2.4-response.csv"))
  d <- d[d$Case == "PCM_partial_anchor", ]
  fit <- fit_mfrm(d, "Person", c("Rater", "Task"), "Score", model = "PCM", step_facet = "Task",
    method = "JML", rating_min = 0, rating_max = 3, optimizer = "BFGS", maxit = 500, reltol = 1e-11)
  plan <- unique(d[c("Rater", "Task")]); plan$Weight <- seq_len(nrow(plan))
  fair <- fair_average_table(fit, reference_plan = plan)
  pa <- mfrmr:::expand_params(fit$opt$par, mfrmr:::build_param_sizes(fit$config), fit$config)
  expected <- vapply(pa$theta, function(t) {
    v <- vapply(seq_len(nrow(plan)), function(j) {
      eta <- t - pa$facets$Rater[plan$Rater[j]] - pa$facets$Task[plan$Task[j]]
      z <- (0:3)*eta - c(0, cumsum(pa$steps_mat[plan$Task[j], ]))
      p <- exp(z-max(z)); sum((0:3)*p)/sum(p)
    }, 0)
    weighted.mean(v, plan$Weight)
  }, 0)
  expect_equal(fair$reference_scores$Estimate, unname(expected), tolerance = 1e-12)
  anchored <- fit_mfrm(d, "Person", c("Rater", "Task"), "Score", model = "PCM", step_facet = "Task",
    method = "JML", rating_min = 0, rating_max = 3, optimizer = "BFGS", maxit = 500, reltol = 1e-11,
    anchors = data.frame(Facet = "Rater", Level = names(pa$facets$Rater)[1], Anchor = .7),
    centering_scope = "anchor_origin")
  translated <- fair_average_table(anchored, reference_plan = plan)
  expect_equal(translated$reference_scores$Estimate, fair$reference_scores$Estimate, tolerance = 2e-5)
  expect_gt(max(abs(anchored$facets$person$Estimate - fit$facets$person$Estimate)), .1)
  ci <- confint(fair, parm = c("P03", "P04"), nsim = 2, seed = 99)
  mc <- attr(ci, "monte_carlo")
  expect_equal(ci$LowerMCLower, mc$Lower[mc$Endpoint == "lower"])
  expect_equal(ci$UpperMCUpper, mc$Upper[mc$Endpoint == "upper"])
  expect_false(any(ci$MCFinite))
  expect_equal(ci$MCLevel, rep(.95, 2))
  expect_identical(plot(ci, draw = FALSE)$data$monte_carlo, mc)
  expect_match(plot(ci, draw = FALSE)$data$notes$Text, "Monte Carlo")
  expect_output(print(ci), "Monte Carlo")
  expect_identical(unserialize(serialize(ci, NULL)), ci)
  expect_equal(as.data.frame(ci)$LowerMCLower, ci$LowerMCLower)
  expect_equal(ci$Resolved + ci$Unresolved, c(2, 2))
  expect_true(all(attr(ci, "trials")$Available))
  expect_true(all(attr(ci, "person_draws") >= 0 & attr(ci, "person_draws") <= 3))
})

test_that("same-response anchor refits reproduce the estimator and paired covariance", {
  case <- fair_reference_case(); n <- nrow(case$fit$prep$data)
  L <- matrix(c(1, -1), 1, dimnames = list("difference", c("P001", "P002")))
  set.seed(981); rng <- .Random.seed
  joint <- confint(case$fair, parm = c("P001", "P002"), contrasts = L, nsim = 3, seed = 71)
  expect_identical(.Random.seed, rng)
  two <- confint(case$anchored, parm = c("P001", "P002"), contrasts = L, nsim = 3, seed = 71,
    anchor_source = case$fit, response_map = list(calibration = seq_len(n), scoring = seq_len(n)))
  expect_identical(.Random.seed, rng)
  expect_true(all(attr(two, "trials")$Available))
  expect_equal(attr(two, "draws"), attr(joint, "draws"), tolerance = 2e-5)
  expect_equal(attr(two, "settings")$shared_responses, n)
  expect_gt(stats::sd(attr(two, "anchor_draws")[, 1]), .01)
  p <- attr(two, "person_draws"); v <- stats::cov(p)
  expect_equal(stats::var(attr(two, "draws")[, 1]), v[1, 1]+v[2, 2]-2*v[1, 2], tolerance = 1e-12)
  expect_gt(abs(v[1, 2]), 1e-5)
  truth <- attr(two, "generating_target")
  root <- attr(two, "draws")[, 1] - truth
  expect_equal(two$CI_Lower, unname(two$Estimate - quantile(root, .975, type = 1)), tolerance = 1e-12)
  expect_false(any(two$CIEligible))
  expect_identical(two, summary(two))
  expect_equal(as.data.frame(two)$Status, two$Status)
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(two, path); saved <- readRDS(path)
  plotted <- plot(saved, draw = FALSE)
  expect_equal(plotted$data$data$CI_Lower, two$CI_Lower)
  expect_equal(plotted$data$data$Unresolved, two$Unresolved)
  expect_match(plotted$data$notes$Text, "Coverage is not qualified", fixed = TRUE)
  expect_error(confint(case$anchored, nsim = 2, seed = 1), "without the calibration fit")
  expect_error(confint(case$anchored, nsim = 2, seed = 1, anchor_source = case$fit), "explicit")
  expect_error(confint(case$anchored, nsim = 2, seed = 1, anchor_source = case$fit,
    response_map = list(calibration = seq_len(n), scoring = rev(seq_len(n)))), "does not match")
})

test_that("overlap maps retain event identity rather than assuming independence", {
  case <- fair_reference_case(); fit <- case$fit
  keep <- which(as.character(fit$prep$data$Person) %in% fit$config$person_levels[1:24])
  args <- fit$config$replay_inputs; args$package_version <- NULL
  args$data <- fit$prep$data[keep, ]; args$anchors <- case$anchors
  args$reltol <- 1e-12
  scored <- do.call(fit_mfrm, args)
  fair <- fair_average_table(scored, reference_plan = case$plan)
  maps <- list(calibration = seq_len(nrow(fit$prep$data)), scoring = keep)
  out <- confint(fair, parm = c("P001", "P002"), nsim = 2, seed = 71,
    anchor_source = fit, response_map = maps)
  expect_identical(attr(out, "response_map"), maps)
  expect_equal(attr(out, "settings")$shared_responses, length(keep))
  expect_equal(attr(out, "generating_target"),
    setNames(case$fair$reference_scores$Estimate[1:2], c("P001", "P002")), tolerance = 1e-12)
  expect_gt(max(abs(out$Estimate - attr(out, "generating_target"))), 1e-5)
  # Check re-centering against generating truth, not the distinct scoring point.
  b <- mfrmr:::mfrm_fair_bootstrap_bounds(out$Estimate, attr(out, "generating_target"),
    attr(out, "draws"), .95)
  expect_equal(out$CI_Lower[out$Available], unname(b$lower[out$Available]))
  expect_equal(out$Resolved + out$Unresolved, rep(2, 2))
})

test_that("unresolved errors widen empirical tails without deleting trials", {
  draws <- matrix(c(1:7, NA, 9, 10), ncol = 1)
  out <- mfrmr:::mfrm_fair_bootstrap_bounds(4, 3, draws, .95)
  expect_equal(out$available, 9)
  expect_identical(unname(out$lower), -Inf)
  expect_identical(unname(out$upper), Inf)
  # At a less extreme quantile, outer bounds enclose every possible completion.
  outer <- mfrmr:::mfrm_fair_bootstrap_bounds(4, 3, draws, .5)
  for (replacement in c(-1e5, 5, 1e5)) {
    complete <- draws; complete[8, 1] <- replacement
    b <- mfrmr:::mfrm_fair_bootstrap_bounds(4, 3, complete, .5)
    expect_lte(outer$lower, b$lower); expect_gte(outer$upper, b$upper)
  }
})

test_that("bounded score limits survive infinite replicate Person measures", {
  d <- expand.grid(Person = paste0("P", 1:12), Rater = c("R1", "R2"),
    Task = c("T1", "T2"), stringsAsFactors = FALSE)
  d$Score <- (as.integer(sub("P", "", d$Person)) + match(d$Rater, c("R1", "R2")) +
    match(d$Task, c("T1", "T2"))) %% 2
  fit <- fit_mfrm(d, "Person", c("Rater", "Task"), "Score", method = "JML",
    rating_min = 0, rating_max = 1, optimizer = "BFGS", maxit = 500, reltol = 1e-12)
  fair <- fair_average_table(fit, reference_plan = unique(d[c("Rater", "Task")]))
  ci <- confint(fair, nsim = 2, seed = 827)
  trials <- attr(ci, "trials"); draws <- attr(ci, "person_draws")
  expect_true(all(trials$Available))
  expect_gt(sum(trials$BoundaryPersons), 0)
  expect_equal(sum(draws == 0 | draws == 1), sum(trials$BoundaryPersons))
  expect_true(all(is.finite(draws)))
  expect_false(any(ci$CIEligible))
})

test_that("refitting retains user-specified factor coordinates", {
  d <- load_mfrmr_data("example_core")
  d$Person <- factor(d$Person, levels = rev(unique(d$Person)))
  d$Rater <- factor(d$Rater, levels = rev(unique(d$Rater)))
  fit <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method = "JML",
    optimizer = "BFGS", maxit = 500, reltol = 1e-12)
  refit <- mfrmr:::mfrm_fair_refit(fit, fit$prep$data$Score)
  expect_identical(refit$prep$levels, fit$prep$levels)
  expect_equal(refit$opt$par, fit$opt$par, tolerance = 1e-10)
})

test_that("a separate union generator supports partially overlapping samples", {
  case <- fair_reference_case(); union <- case$fit
  person <- as.character(union$prep$data$Person)
  calibration_rows <- which(person %in% union$config$person_levels[1:32])
  scoring_rows <- which(person %in% union$config$person_levels[17:48])
  args <- union$config$replay_inputs; args$package_version <- NULL; args$reltol <- 1e-12
  args$data <- union$prep$data[calibration_rows, ]
  calibration <- do.call(fit_mfrm, args)
  args$anchors <- case$anchors
  args$anchors$Anchor <- mfrmr:::mfrm_fair_anchor_values(calibration, args$anchors)
  args$data <- union$prep$data[scoring_rows, ]
  scoring <- do.call(fit_mfrm, args)
  fair <- fair_average_table(scoring, reference_plan = case$plan)
  ci <- confint(fair, parm = c("P025", "P040"), nsim = 2, seed = 71,
    anchor_source = calibration, generating_fit = union,
    response_map = list(calibration = calibration_rows, scoring = scoring_rows))
  expect_true(all(attr(ci, "trials")$Available))
  expect_equal(attr(ci, "settings")$shared_responses, length(intersect(calibration_rows, scoring_rows)))
  expect_equal(attr(ci, "generating_target"), setNames(case$fair$reference_scores$Estimate[
    match(c("P025", "P040"), case$fair$reference_scores$Person)], c("P025", "P040")))
})

test_that("a failed refit remains in the requested interval denominator", {
  case <- fair_reference_case()
  real_refit <- mfrmr:::mfrm_fair_refit; calls <- 0L
  testthat::local_mocked_bindings(mfrm_fair_refit = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop("Controlled numerical refit failure")
    real_refit(...)
  }, .package = "mfrmr")
  ci <- confint(case$fair, parm = "P001", nsim = 2, seed = 71)
  expect_equal(ci$Planned, 2)
  expect_equal(ci$Unresolved, 1)
  expect_match(attr(ci, "trials")$Reason[1], "Controlled numerical refit failure")
  expect_identical(ci$Status, "unresolved_tail")
  expect_true(is.na(ci$CI_Lower) && is.na(ci$CI_Upper))
  expect_equal(nrow(attr(ci, "draws")), 2)
})

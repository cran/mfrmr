# Independent enumeration uses the model equation and complete response
# patterns, not the production probability kernel or convolution.
classification_patterns <- function(theta, offsets, cumulative) {
  k <- seq_len(ncol(cumulative)) - 1L
  patterns <- as.matrix(expand.grid(rep(list(k), length(offsets))))
  prob <- apply(patterns, 1, function(y) prod(vapply(seq_along(y), function(j) {
    numerator <- exp(y[j] * (theta + offsets[j]) - cumulative[j, y[j] + 1L])
    numerator / sum(exp(k * (theta + offsets[j]) - cumulative[j, ]))
  }, 0)))
  list(patterns = patterns, probability = prob, sum = rowSums(patterns))
}

classification_fixture <- function() {
  readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits
}

test_that("the public RSM route preserves shared polytomous thresholds", {
  patterns <- expand.grid(A = 0:2, B = 0:2, C = 0:2)
  d <- data.frame(Person = rep(seq_len(nrow(patterns)), each = 3),
    Item = rep(c("A", "B", "C"), nrow(patterns)), Score = as.vector(t(patterns)))
  for (method in c("MML", "JML")) {
    fit <- fit_mfrm(d, "Person", "Item", "Score", method = method)
    # All-low/all-high patterns are excluded by JML calibration. Preserve
    # that source restriction while checking the conditional arithmetic.
    if (method == "JML") expect_error(classification_table(fit,
      data.frame(Item = c("A", "C")), 0, 2, c(-1, 0, 1)), "boundary_not_scoring_ready:has_exclusions")
    x <- classification_table(fit, data.frame(Item = c("A", "C")), 0, 2, c(-1, 0, 1),
      readiness_policy = if (method == "JML") "review" else "error")
    expect_identical(x$source$ready, method == "MML")
    params <- expand_params(fit$opt$par, build_param_sizes(fit$config), fit$config)
    expect_equal(x$calibration$cumulative_steps,
      matrix(c(0, cumsum(params$steps)), 2, 3, byrow = TRUE))
    for (i in 1:3) {
      ref <- classification_patterns(c(-1, 0, 1)[i], x$calibration$offsets, x$calibration$cumulative_steps)
      expect_equal(subset(x$probabilities, Point == i & Method == "lee_sum" & Category == 2)$Probability,
        sum(ref$probability[ref$sum >= 2]), tolerance = 1e-13)
    }
  }
})

test_that("finite recursion agrees with independently enumerated binary and polytomous patterns", {
  for (k in c(2L, 4L)) for (shared in c(TRUE, FALSE)) for (theta in c(-3, 0, 2)) {
    offsets <- c(-.7, .2, 1)
    steps <- if (shared) matrix(seq(-.8, .8, length.out = k - 1), 3, k - 1, byrow = TRUE) else
      matrix(seq(-1, 1.5, length.out = 3 * (k - 1)), 3, k - 1)
    cumulative <- t(apply(steps, 1, function(s) c(0, cumsum(s))))
    cuts <- c(1, 2 * (k - 1))
    x <- mfrm_classification_point(theta, offsets, cumulative, cuts, c(-.5, .5))
    ref <- classification_patterns(theta, offsets, cumulative)
    expected_pmf <- vapply(0:(3 * (k - 1)), function(s) sum(ref$probability[ref$sum == s]), 0)
    band <- 1L + (ref$sum >= cuts[1]) + (ref$sum >= cuts[2])
    expected_bands <- vapply(1:3, function(h) sum(ref$probability[band == h]), 0)
    expect_equal(x$pmf, expected_pmf, tolerance = 2e-14)
    expect_equal(x$lee, expected_bands, tolerance = 2e-14)
    # Explicit paired patterns verify repeated-administration agreement.
    pair_probability <- outer(ref$probability, ref$probability)
    expect_equal(sum(x$lee^2), sum(pair_probability[outer(band, band, "==")]), tolerance = 2e-14)
    total_mean <- sum(ref$sum * ref$probability)
    total_var <- sum((ref$sum - total_mean)^2 * ref$probability)
    expect_equal(x$expected, total_mean, tolerance = 2e-14)
    expect_equal(x$information, total_var, tolerance = 2e-14)
    expect_equal(x$sem, 1 / sqrt(total_var), tolerance = 2e-14)
  }
})

test_that("short and long binary plans agree with the binomial reference", {
  for (n in c(3, 20, 80)) for (theta in c(-8, 0, 8)) {
    cut <- ceiling(n / 2)
    x <- mfrm_classification_point(theta, rep(0, n), matrix(0, n, 2), cut, 0)
    expect_equal(x$pmf, dbinom(0:n, n, plogis(theta)), tolerance = 3e-14)
    expect_equal(x$lee[2], pbinom(cut - 1, n, plogis(theta), lower.tail = FALSE), tolerance = 3e-14)
    expect_equal(x$information, n * plogis(theta) * plogis(-theta), tolerance = 3e-14)
  }
})

test_that("score endpoints and repeated cuts preserve every category", {
  x <- mfrm_classification_point(0, rep(0, 3), matrix(0, 3, 2), c(0, 2, 2, 4), -2:1)
  expect_equal(x$lee, c(0, .5, 0, .5, 0))
  expect_equal(sum(x$pmf), 1)
  y <- mfrm_classification_point(0, rep(0, 3), matrix(0, 3, 2), 3, 0)
  expect_equal(y$lee[2], 1 / 8) # The maximum score belongs to the upper band.
})

test_that("information equals the derivative of the expected sum and preserves tails", {
  offsets <- c(-.5, .8)
  cumulative <- rbind(c(0, -.3, .6), c(0, .2, 1.3))
  t <- .4; h <- 1e-5
  a <- mfrm_classification_response(t, offsets, cumulative)
  derivative <- (mfrm_classification_response(t + h, offsets, cumulative)$expected -
                 mfrm_classification_response(t - h, offsets, cumulative)$expected) / (2 * h)
  expect_equal(a$information, derivative, tolerance = 1e-9)
  tail <- mfrm_classification_response(40, 0, matrix(c(0, 0), 1))
  expect_gt(tail$information, 0)
  expect_equal(tail$information / (plogis(40) * plogis(-40)), 1, tolerance = 1e-14)
  expect_error(mfrm_classification_response(1e308, 1e308, matrix(c(0, 0), 1)), "numerically unresolved")
})

test_that("location changes preserve the scale and severity shifts change decisions in the expected direction", {
  offsets <- c(-.6, .3, .9)
  cumulative <- matrix(c(0, -.4, .3), 3, 3, byrow = TRUE)
  original <- mfrm_classification_point(.2, offsets, cumulative, c(2, 5), c(-.5, .5))
  shifted <- mfrm_classification_point(10.2, offsets - 10, cumulative, c(2, 5), c(9.5, 10.5))
  expect_equal(original, shifted, tolerance = 1e-14)
  severe <- mfrm_classification_point(.2, offsets - 1, cumulative, c(2, 5), c(-.5, .5))
  expect_lt(severe$expected, original$expected)
  expect_lt(severe$lee[3], original$lee[3])
})

test_that("normal comparison uses sampling information and stable upper tails", {
  x <- mfrm_classification_point(0, rep(0, 4), matrix(0, 4, 2), c(2, 3), c(8, 9))
  expect_equal(x$sem, 1)
  expect_equal(x$rudner[2] / (pnorm(8, lower.tail = FALSE) - pnorm(9, lower.tail = FALSE)), 1)
  narrow <- mfrm_classification_point(0, rep(0, 4), matrix(0, 4, 2), c(2, 3), c(-1, 1) * qnorm(.51))
  expect_equal(narrow$rudner, c(.49, .02, .49), tolerance = 1e-14)
  expect_lt(narrow$rudner[2], sum(narrow$rudner^2))
})

test_that("narrow normal bands retain probability before CDF and location cancellation", {
  cases <- list(c(0, 0, 1e-18), c(1, 0, 1e-18),
    c(0, 8, 8 + 1e-10), c(0, -8, -8 + 1e-10), c(0, -1e-18, 1e-18))
  for (value in cases) {
    theta <- value[1L]; cuts <- value[2:3]
    x <- mfrm_classification_point(theta, rep(-theta, 4), matrix(0, 4, 2), c(2, 3), cuts)
    # The midpoint integral has relative error O(width^2) here (width <= 1e-10).
    expected <- dnorm(mean(cuts) - theta) * diff(cuts)
    expect_equal(x$sem, 1)
    expect_gt(x$rudner[2L], 0)
    expect_equal(x$rudner[2L] / expected, 1, tolerance = 1e-12)
    expect_equal(sum(x$rudner), 1, tolerance = 1e-14)
  }
  fit <- classification_fixture()[["PCM-MML-fixed"]]
  plan <- expand.grid(lapply(fit$config$facet_levels, head, 2), stringsAsFactors = FALSE)
  x <- classification_table(fit, plan, cut_scores = c(0, 1e-18), score_cuts = c(1, 2),
    theta = 0, theta_weights = 1, readiness_policy = "review")
  middle <- subset(x$probabilities, Method == "rudner_normal" & Category == 2)$Probability
  expect_gt(middle, 0)
  expect_equal(subset(x$table, Method == "rudner_normal")$Accuracy, middle)
  expect_equal(subset(x$marginal, Method == "rudner_normal")$Accuracy, middle)
  expect_equal(subset(x$classification_matrix,
    Method == "rudner_normal" & TrueCategory == 2 & ObservedCategory == 2)$Probability, middle)
})

test_that("ordinary ML and sum cuts agree when their actual decision rules match", {
  offsets <- c(-1, .2, .7)
  cumulative <- rbind(c(0, -.5, .2), c(0, .3, .9), c(0, -.2, 0))
  ref <- classification_patterns(0, offsets, cumulative)
  # Maximize the raw pattern log likelihood independently of the production
  # moment function; endpoints use the infinite-ML convention in this check.
  estimates <- apply(ref$patterns, 1, function(y) {
    if (sum(y) == 0) return(-Inf)
    if (sum(y) == 6) return(Inf)
    optimize(function(t) sum(vapply(1:3, function(j) {
      y[j] * (t + offsets[j]) - cumulative[j, y[j] + 1L] -
        log(sum(exp((0:2) * (t + offsets[j]) - cumulative[j, ])))
    }, 0)), c(-20, 20), maximum = TRUE, tol = 1e-10)$maximum
  })
  for (cut in c(-.8, .6)) {
    # Independent expected sum from complete patterns at the ability cut.
    at_cut <- classification_patterns(cut, offsets, cumulative)
    raw_cut <- ceiling(sum(at_cut$sum * at_cut$probability))
    expect_identical(estimates >= cut, ref$sum >= raw_cut)
  }
})

test_that("native MML and JML plans retain their calibration and exact probabilities", {
  for (fit in classification_fixture()) {
    plan <- data.frame(Rater = c("1", "2", "3"), Criterion = c("3", "1", "2"))
    x <- classification_table(fit, plan, c(Low = -.5, High = .5), c(3, 7), c(-1, -.5, 0, .5, 1))
    expect_true(x$source$ready)
    expect_equal(x$table$TrueCategory, rep(c(1L, 2L, 2L, 3L, 3L), each = 2))
    expect_equal(nrow(x$marginal), 0L)
    expect_equal(nrow(x$classification_matrix), 0L)
    expect_true(all(is.na(x$population$Weight)))
    p <- expand_params(fit$opt$par, build_param_sizes(fit$config), fit$config)
    direct_offset <- -p$facets$Rater[match(plan$Rater, fit$config$facet_levels$Rater)] -
      p$facets$Criterion[match(plan$Criterion, fit$config$facet_levels$Criterion)]
    expect_equal(unname(x$calibration$offsets), unname(direct_offset))
    for (i in 1:5) {
      ref <- classification_patterns(x$population$Theta[i], direct_offset, x$calibration$cumulative_steps)
      bands <- 1L + (ref$sum >= 3) + (ref$sum >= 7)
      expected <- vapply(1:3, function(b) sum(ref$probability[bands == b]), 0)
      actual <- subset(x$probabilities, Point == i & Method == "lee_sum")$Probability
      expect_equal(actual, expected, tolerance = 1e-13)
      true <- x$table$TrueCategory[x$table$Point == i][1]
      expect_equal(subset(x$table, Point == i & Method == "lee_sum")$Accuracy, expected[true])
    }
    expect_equal(x$table$Accuracy + x$table$Overclassification + x$table$Underclassification,
                 rep(1, nrow(x$table)), tolerance = 1e-14)
    expect_equal(x$cuts$Cut, c("Low", "High"))
    for (i in 1:2) {
      at <- subset(x$table, Theta == x$cuts$TrueCut[i])
      expect_equal(at$SEM, rep(x$cuts$SEMAtTrueCut[i], 2))
    }
    swapped <- classification_table(fit, plan[3:1, ], c(-.5, .5), c(3, 7), c(-1, -.5, 0, .5, 1))
    expect_equal(x$table, swapped$table, tolerance = 1e-14)
  }
})

test_that("population weights give auditable joint cells and preserve unavailable points", {
  fit <- classification_fixture()[[2]]
  plan <- data.frame(Rater = "1", Criterion = c("1", "2"))
  x <- classification_table(fit, plan, c(-.5, .5), c(2, 4), c(-1, 0, 1),
    theta_weights = c(1e308, 1e308, 5e307), operational_cuts = c(-.2, .7))
  expect_equal(x$population$NormalizedWeight, c(.4, .4, .2))
  expect_equal(x$cuts$OperationalAbilityCut, c(-.2, .7))
  for (method in c("lee_sum", "rudner_normal")) {
    j <- subset(x$classification_matrix, Method == method)
    m <- subset(x$marginal, Method == method)
    expect_equal(sum(j$Probability), 1)
    expect_equal(sum(j$Probability[j$TrueCategory == j$ObservedCategory]), m$Accuracy)
    expect_equal(sum(j$Probability[j$TrueCategory < j$ObservedCategory]), m$Overclassification)
    expect_equal(sum(j$Probability[j$TrueCategory > j$ObservedCategory]), m$Underclassification)
    expect_equal(m$Consistency, sum(subset(x$table, Method == method)$Consistency * c(.4, .4, .2)))
  }
  extreme <- classification_table(fit, plan, 0, 3, c(0, 1000), theta_weights = c(1, 0))
  expect_equal(extreme$marginal$Status, c("available", "available"))
  expect_equal(extreme$table$Status, c("available", "available", "available", "zero_information"))
  expect_true(is.na(extreme$table$Underclassification[4]))
  expect_true(is.infinite(extreme$table$SEM[4]))
  positive <- classification_table(fit, plan, 0, 3, c(0, 1000), theta_weights = c(1, 1))
  expect_equal(positive$marginal$Status, c("available", "positive_weight_point_unavailable"))
  expect_true(all(is.na(subset(positive$classification_matrix, Method == "rudner_normal")$Probability)))
  tiny_positive <- classification_table(fit, plan, 0, 3, c(0, 1000), theta_weights = c(1e308, 1e-308))
  expect_equal(tiny_positive$marginal$Status[2], "positive_weight_point_unavailable")
})

test_that("classification tables retain source review separately from joint error probabilities", {
  fit <- classification_fixture()[[2]]
  fit$config$estimability_audit$nonlinear_local_estimability$local_first_order_rank_deficient <- TRUE
  plan <- data.frame(Rater = "1", Criterion = c("1", "2"))
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"), .package = "mfrmr")
  x <- classification_table(fit, plan, 0, 3, c(-1, 1), theta_weights = c(3, 1),
    readiness_policy = "review")
  expect_identical(x$source$status, "review_only")
  for (name in c("table", "probabilities", "scores", "cuts", "marginal", "classification_matrix")) {
    expect_true(all(x[[name]]$SourceStatus == "review_only"), info = name)
    csv <- tempfile(fileext = ".csv")
    write.csv(x[[name]], csv, row.names = FALSE)
    expect_identical(read.csv(csv)$SourceStatus, x[[name]]$SourceStatus)
    unlink(csv)
  }
  expect_true(all(x$probabilities$Status == "available"))
  expect_true(all(x$classification_matrix$Status == "available"))
  expect_match(paste(x$notes, collapse = " "), "joint population probabilities")
  for (method in c("lee_sum", "rudner_normal")) {
    joint <- subset(x$classification_matrix, Method == method)
    expect_equal(sum(joint$Probability[joint$TrueCategory == 1]), .75)
    expect_equal(sum(joint$Probability[joint$TrueCategory == 2]), .25)
    conditional_pass <- subset(x$probabilities, Method == method & Point == 1 & Category == 2)$Probability
    conditional_fail <- subset(x$probabilities, Method == method & Point == 2 & Category == 1)$Probability
    marginal <- subset(x$marginal, Method == method)
    expect_equal(marginal$Overclassification, .75 * conditional_pass)
    expect_equal(marginal$Underclassification, .25 * conditional_fail)
  }
  without_population <- classification_table(fit, plan, 0, 3, 0, readiness_policy = "review")
  for (name in c("marginal", "classification_matrix")) {
    expect_equal(nrow(without_population[[name]]), 0L)
    expect_identical(without_population[[name]]$SourceStatus, character())
  }
})

test_that("plan mapping and score recoding never change model category sums", {
  fit <- classification_fixture()[[2]]
  plan <- data.frame(Rater = c("1", "2"), Criterion = c("1", "2"))
  original <- classification_table(fit, plan, 0, 3, 0)
  names(plan) <- c("Judge", "Task")
  fit$prep$score_map$OriginalScore <- c(2, 5, 10, 20)
  recoded <- classification_table(fit, plan, 0, 3, 0, facets = c(Rater = "Judge", Criterion = "Task"))
  expect_equal(original$table, recoded$table)
  expect_equal(recoded$calibration$score_map$CategoryScore, 0:3)
  expect_equal(recoded$calibration$score_map$OriginalScore, c(2, 5, 10, 20))
})

test_that("public inputs and source eligibility fail explicitly", {
  fit <- classification_fixture()[[2]]
  plan <- data.frame(Rater = "1", Criterion = "1")
  run <- function(f = fit, p = plan, cuts = 0, scores = 2, theta = 0, ...) classification_table(f, p, cuts, scores, theta, ...)
  expect_error(run(cuts = NULL), "at least one")
  expect_error(run(cuts = c(1, 0), scores = c(1, 2)), "strictly increasing")
  expect_error(run(operational_cuts = c(0, 1)), "one.*per")
  for (bad in list(NA_real_, numeric(), Inf, "0", matrix(0), factor("0")))
    expect_error(run(theta = bad), "theta")
  for (bad in list(-1, .5, 5, NA_real_, "2", c(1, 2)))
    expect_error(run(scores = bad), "score_cuts")
  for (bad in list(0, -1, NA_real_, c(1, 2), Inf))
    expect_error(run(theta_weights = bad), "theta_weights")
  expect_error(run(p = rbind(plan, plan)), "Duplicate")
  expect_error(run(p = transform(plan, Weight = 1)), "exactly")
  expect_error(run(p = transform(plan, Rater = "unknown")), "Unknown")
  expect_error(run(p = transform(plan, Rater = NA_character_)), "nonmissing")
  expect_error(run(p = transform(plan, Rater = " ")), "nonblank")
  expect_error(run(facets = c(NA_character_, "Criterion")), "missing")
  bad <- fit; bad$config$model <- "GPCM"
  expect_error(run(f = bad), "ordinary RSM/PCM")
  bad <- fit; bad$config$interaction_specs <- list(list())
  expect_error(run(f = bad), "without interactions")
  bad <- fit; class(bad) <- c("mfrm_random_rater", class(bad))
  expect_error(run(f = bad), "random effects")
  bad <- fit; bad$opt$par[1] <- NA_real_
  expect_error(run(f = bad, readiness_policy = "review"), "invalid")
  review <- fit
  review$config$estimability_audit$nonlinear_local_estimability$local_first_order_rank_deficient <- TRUE
  expect_error(run(f = review), "not ready")
  x <- run(f = review, readiness_policy = "review")
  expect_false(x$source$ready)
  expect_equal(x$table$SourceStatus, rep(x$source$status, nrow(x$table)))
  expect_true(length(x$source$reason_codes) > 0)
  expect_equal(x$settings$readiness_policy, "review")
})

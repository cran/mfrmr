shared_rater_example <- function(persons = 120, raters = 8, sd = .7, seed = 92351) {
  withr::with_seed(seed, {
    p <- rep(seq_len(persons), each = 6)
    r <- (p - 1 + rep(rep(0:1, each = 3), persons)) %% raters + 1
    criterion <- rep(1:3, 2 * persons)
    theta <- rnorm(persons); severity <- rnorm(raters, sd = sd)
    eta <- theta[p] - severity[r] - c(-.3, 0, .3)[criterion]
    w <- cbind(0, eta + .6, 2 * eta); w <- exp(w - apply(w, 1, max)); w <- w / rowSums(w)
    y <- rowSums(runif(length(p)) > t(apply(w, 1, cumsum)))
    data.frame(Person = p, Rater = r, Criterion = criterion, Score = y)
  })
}

shared_rater_fit <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) saved <<- fit_mfrm_random_rater(shared_rater_example(), "Person", "Rater", "Score", "Criterion", 0:2, person_sd = 1)
    saved
  }
})

test_that("random-rater input preserves roles, categories and omitted assigned scores", {
  d <- shared_rater_example()
  x <- mfrm_random_rater_data(d, "Person", "Rater", "Criterion", "Score", 0:2)
  expect_equal(colSums(x$basis$Criterion), c(C1 = 0, C2 = 0), tolerance = 1e-12)
  expect_equal(crossprod(x$basis$Criterion), diag(2), ignore_attr = TRUE)
  d$Score[1] <- NA
  expect_error(mfrm_random_rater_data(d, "Person", "Rater", "Criterion", "Score", 0:2), "Assigned scores")
  x <- mfrm_random_rater_data(d, "Person", "Rater", "Criterion", "Score", 0:2, "omit")
  expect_equal(x$omitted_rows, 1)
  expect_equal(nrow(x$data), nrow(d) - 1)
  bad <- d; bad$Person[1] <- NA
  expect_error(mfrm_random_rater_data(bad, "Person", "Rater", "Criterion", "Score", 0:2, "omit"), "identifiers")
  bad <- d; bad$Other <- bad$Criterion
  expect_error(mfrm_random_rater_data(bad, "Person", "Rater", c("Criterion", "Other"), "Score", 0:2, "omit"), "aliased")
  bad <- d; bad$Rater <- (bad$Person - 1) %% 8
  expect_error(mfrm_random_rater_data(bad, "Person", "Rater", "Criterion", "Score", 0:2, "omit"), "connected")
  expect_error(fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:2, rater_sd = -1), "nonnegative")
  expect_error(fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:2, quad_points = 3), "quad_points")
  expect_error(mfrm_random_rater_data(d, "Person", "Rater", "Criterion", "Score", c(0, 2)), "consecutive")
})

test_that("shared likelihood, derivatives and variance-zero score have independent references", {
  skip_if_not_installed("RTMB", "2.0")
  d <- expand.grid(Person = 1:2, Rater = 1:2, Criterion = 1:2)
  d$Score <- c(0, 1, 1, 2, 0, 0, 1, 1)
  input <- mfrm_random_rater_data(d, "Person", "Rater", "Criterion", "Score", 0:2)
  beta <- as.vector(crossprod(input$basis$Criterion, c(-.3, .3)))
  fixed <- c(beta, -.6, .6)
  joint <- mfrm_random_rater_objective(input, 41, fixed_sd = .7, random = FALSE)
  # Direct probabilities and integration, without RTMB or the backend's logsum operations.
  rule <- gauss_hermite_normal(41)
  direct <- function(par) {
    total <- sum(stats::dnorm(par[4:5], log = TRUE))
    for (p in 1:2) {
      likelihood <- rep(1, length(rule$nodes))
      for (i in which(input$person == p)) {
        eta <- rule$nodes - input$X[i, 1] * par[1] - .7 * par[3 + input$rater[i]]
        w <- cbind(1, exp(eta - par[2]), exp(2 * eta - sum(par[2:3])))
        likelihood <- likelihood * (w / rowSums(w))[, input$y[i] + 1]
      }
      total <- total + log(sum(rule$weights * likelihood))
    }
    -total
  }
  par <- c(fixed, .3, -.4)
  expect_equal(as.numeric(joint$fn(par)), direct(par), tolerance = 1e-10)
  gradient <- vapply(seq_along(par), function(i) {
    h <- rep(0, length(par)); h[i] <- 1e-5
    (direct(par + h) - direct(par - h)) / 2e-5
  }, numeric(1))
  expect_equal(as.vector(joint$gr(par)), gradient, tolerance = 1e-7)
  inner <- stats::optim(c(0, 0), function(z) direct(c(fixed, z)), method = "BFGS",
    control = list(reltol = 1e-12))
  H <- stats::optimHess(inner$par, function(z) direct(c(fixed, z)))
  laplace <- inner$value + as.numeric(determinant(H, logarithm = TRUE)$modulus) / 2 - log(2 * pi)
  marginal <- mfrm_random_rater_objective(input, 41, fixed_sd = .7)
  expect_equal(as.numeric(marginal$fn(fixed)), laplace, tolerance = 2e-6)
  # Archived independent tensor integral is -7.61845115123505. Approximation
  # error is real and is not confused with the backend implementation error.
  expect_gt(abs(as.numeric(marginal$fn(fixed)) - 7.61845115123505), .001)
  zero <- mfrm_random_rater_objective(input, 41, fixed_sd = 0)
  h <- 1e-4
  q1 <- mfrm_random_rater_objective(input, 41, fixed_sd = sqrt(h))
  q2 <- mfrm_random_rater_objective(input, 41, fixed_sd = sqrt(h / 2))
  slope <- 2 * (zero$fn(fixed) - q2$fn(fixed)) / (h / 2) - (zero$fn(fixed) - q1$fn(fixed)) / h
  expect_equal(mfrm_random_rater_zero_score(input, beta, c(-.6, .6), 41), as.numeric(slope), tolerance = 1e-5)
  relabeled <- input; relabeled$rater <- 3L - input$rater
  swapped <- mfrm_random_rater_objective(relabeled, 41, fixed_sd = .7)
  expect_equal(as.numeric(swapped$fn(fixed)), as.numeric(marginal$fn(fixed)), tolerance = 1e-9)
})

test_that("shared predictor cells retain unequal counts, focal abilities and replicates", {
  skip_if_not_installed("RTMB", "2.0")
  d <- expand.grid(Person = 1:4, Rater = 1:2,
    Task = c("a:1", "a"), Criterion = c("b", "1:b"))
  d$Score <- rep(1:4, length.out = nrow(d))
  d <- rbind(d, d[c(1, 1, 5, 9, 13, 17), ])
  d$Score[c(2, 11)] <- NA
  input <- mfrm_random_rater_data(d, "Person", "Rater", c("Task", "Criterion"),
    "Score", 1:4, "omit")
  rule <- gauss_hermite_normal(21L)
  # Direct row-wise probabilities deliberately do not aggregate design cells.
  for (case in list(c(.7, 1.2, NA, NA), c(.7, 1.2, 1, 3), c(0, 0, 1, 1))) {
    focal <- if (is.na(case[3])) NULL else case[3]
    replicate <- if (is.na(case[4])) NULL else case[4]
    direct <- function(par) {
      beta <- par[1:2]; steps <- par[3:5]; z <- par[6:7]
      sd <- par[8]; ps <- par[9]
      focal_ability <- if (is.null(focal)) NULL else par[10]
      extra <- if (is.null(replicate)) NULL else par[11:14]
      logprob <- function(theta, i, score, step_weights) {
        eta <- theta - sum(input$X[i, ] * beta) - sd * z[input$rater[i]]
        lw <- (0:3) * eta - c(0, cumsum(steps))
        peak <- max(lw)
        score * eta - sum(step_weights * steps) - peak - log(sum(exp(lw - peak)))
      }
      likelihood <- vapply(1:4, function(p) {
        at <- which(input$person == p)
        prob <- function(theta) exp(sum(vapply(at, function(i)
          logprob(theta, i, input$y[i], as.numeric(1:3 <= input$y[i])), numeric(1))) +
          if (!is.null(replicate) && input$person[replicate] == p)
            logprob(theta, replicate, extra[1], extra[-1]) else 0)
        if (identical(as.numeric(p), focal)) prob(focal_ability) else
          sum(rule$weights * vapply(ps * rule$nodes, prob, numeric(1)))
      }, numeric(1))
      -sum(log(likelihood)) - sum(dnorm(z, log = TRUE))
    }
    obj <- mfrm_random_rater_objective(input, 21L, fixed_sd = NULL,
      fixed_person_sd = NULL, random = FALSE, focal_person = focal,
      replicate_row = replicate)
    par <- c(.2, -.3, -.8, .1, .7, .25, -.4, case[1:2],
      if (!is.null(focal)) .35, if (!is.null(replicate)) c(2, 1, 1, 0))
    expect_equal(as.numeric(obj$fn(par)), direct(par), tolerance = 1e-10)
    gradient <- vapply(seq_along(par), function(i) {
      h <- numeric(length(par)); h[i] <- 1e-5
      (direct(par + h) - direct(par - h)) / 2e-5
    }, numeric(1))
    expect_equal(as.vector(obj$gr(par)), gradient, tolerance = 1e-7)
  }
})

test_that("fitting, covariance and replacement-rater prediction keep their targets", {
  skip_if_not_installed("RTMB", "2.0")
  fit <- shared_rater_fit()
  expect_s3_class(fit, "mfrm_random_rater")
  expect_true(fit$checks$NumericalReady)
  expect_true(fit$checks$InformationPositive)
  expect_gt(fit$calibration$rater_sd, 0)
  expect_true(all(is.finite(fit$raters$PredictionSE)))
  expect_true(all(is.na(fit$raters$Lower) & is.na(fit$raters$Upper)))
  expect_true(all(fit$raters$PredictionSE >= fit$raters$ConditionalSD - 1e-8))
  expect_equal(fit$rater_covariance, t(fit$rater_covariance), tolerance = 1e-10)
  expect_identical(dimnames(fit$rater_covariance), list(fit$raters$Rater, fit$raters$Rater))
  expect_identical(dimnames(fit$conditional_rater_covariance), dimnames(fit$rater_covariance))
  expect_gt(max(abs(fit$rater_covariance[row(fit$rater_covariance) != col(fit$rater_covariance)])), 1e-4)
  # Derive mode sensitivity by implicit differentiation of the joint Hessian;
  # this checks the covariance meaning independently of sdreport's assembly.
  joint <- mfrm_random_rater_objective(fit$input, 31, random = FALSE)
  sd <- fit$calibration$rater_sd; z <- fit$rater_mode / sd
  par <- c(fit$calibration$beta, fit$calibration$steps, z, sd)
  H <- stats::optimHess(par, joint$fn, joint$gr)
  iz <- length(fit$calibration$beta) + length(fit$calibration$steps) + seq_along(z)
  ip <- setdiff(seq_along(par), iz)
  conditional <- sd^2 * solve(H[iz, iz])
  J <- -sd * solve(H[iz, iz], H[iz, ip])
  J[, ncol(J)] <- J[, ncol(J)] + z
  expect_equal(conditional, fit$conditional_rater_covariance, tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(conditional + J %*% fit$covariance %*% t(J), fit$rater_covariance,
    tolerance = 1e-5, ignore_attr = TRUE)
  expect_error(confint(fit, parm = "steps"), "rater_sd")
  expect_error(confint(fit, level = 1), "level")
  f <- tempfile(); saveRDS(fit, f)
  expect_identical(readRDS(f), fit)
  new <- data.frame(Rater = "New", Criterion = "2", Ability = c(-1, 0, 1))
  pred <- predict(fit, new, "Ability", rater = "new")
  expect_equal(rowSums(pred$probabilities), rep(1, 3), tolerance = 1e-12)
  expect_true(all(diff(pred$expected_scores$ExpectedScore) > 0))
  expect_identical(predict(readRDS(f), new, "Ability", rater = "new"), pred)
  expect_error(predict(fit, new, "Ability"), "known rater")
  known <- transform(new, Rater = "1")
  expect_error(predict(fit, known, "Ability", rater = "new"), "new rater")
  observed <- predict(fit, known, "Ability")
  expect_gt(max(abs(pred$probabilities - observed$probabilities)), .01)
  point <- fit; point$calibration$rater_sd <- 0
  expect_gt(max(abs(pred$probabilities - predict(point, new, "Ability", rater = "new")$probabilities)), .001)
  expect_error(predict(fit, transform(known, Criterion = "Unknown"), 0), "unknown fixed")
  before <- grDevices::dev.cur()
  plot <- plot(fit, draw = FALSE)
  expect_identical(grDevices::dev.cur(), before)
  expect_equal(plot_data(plot)$table[names(fit$raters)], fit$raters)
  expect_true(all(plot_data(plot)$table$IntervalStatus == "not_requested"))
  expect_false(any(plot_data(plot)$table$CIEligible))
  grDevices::pdf(tempfile(fileext = ".pdf")); before <- graphics::par("mar")
  plot(fit); expect_equal(graphics::par("mar"), before); grDevices::dev.off()
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(as_ggplot(plot), "ggplot")
})

test_that("positive random-rater variance profiles attain the likelihood cutoff", {
  # Keep the expensive nuisance refits in the full suite; ordinary checks retain
  # a live boundary profile and independent analytic endpoint/curvature tests.
  skip_on_cran()
  skip_if_not_installed("RTMB", "2.0")
  fit <- shared_rater_fit()
  interval <- confint(fit)
  expect_lt(interval[1, "Lower"], fit$calibration$rater_sd)
  expect_gt(interval[1, "Upper"], fit$calibration$rater_sd)
  profile <- attr(interval, "profile")
  expect_true(all(profile$NumericalReady))
  ends <- vapply(interval[1, ], function(sd) profile$NLL[which.min(abs(profile$SD - sd))], numeric(1))
  expect_equal(unname(2 * (ends + fit$loglik)), rep(qchisq(.95, 1), 2), tolerance = 2e-4)
})

test_that("an estimated zero variance does not manufacture certainty", {
  # Extended numerical/reference case; shorter live routes remain on CRAN.
  skip_on_cran()
  skip_if_not_installed("RTMB", "2.0")
  d <- expand.grid(Person = 1:12, Rater = 1:4, Criterion = 1:3)
  d$Score <- (d$Person + d$Criterion) %% 3
  expect_warning(coarse <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:2, person_sd = 1), "Person quadrature is not stable.*Refit")
  expect_false(coarse$checks$NumericalReady)
  expect_false(coarse$checks$PersonQuadratureStable)
  expect_error(confint(coarse), "numerical")
  expect_error(predict(coarse, data.frame(Rater = "1", Criterion = "1"), 0), "numerical")
  fit <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:2, quad_points = 121, person_sd = 1)
  expect_true(fit$checks$NumericalReady)
  expect_true(fit$checks$PersonQuadratureStable)
  expect_true(fit$checks$EstimatedVarianceBoundary)
  expect_equal(fit$calibration$rater_sd, 0)
  expect_lte(fit$checks$HigherOrderZeroVarianceScore, 1e-5)
  expect_true(all(is.na(fit$raters$PredictionSE)))
  expect_true(all(is.na(fit$calibration_table$SE)))
  expect_equal(fit$raters$Estimate, rep(0, 4))
  expect_true(any(grepl("withheld", summary(fit)$notes)))
  interval <- confint(fit)
  expect_equal(interval[1, "Lower"], 0)
  expect_gt(interval[1, "Upper"], 0)
  expect_true(all(attr(interval, "profile")$NumericalReady))
  known <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:2,
    rater_sd = 0, quad_points = 121, person_sd = 1)
  expect_false(known$checks$EstimatedVarianceBoundary)
  expect_equal(known$loglik, fit$loglik, tolerance = 1e-8)
  expect_true(all(is.finite(known$calibration_table$SE[known$calibration_table$Parameter != "Population SD"])))
  expect_equal(known$raters$PredictionSE, rep(0, 4))
  expect_error(confint(known), "fixed")
})

test_that("a false-convergence zero submodel receives a verified curvature restart", {
  # Extended numerical/reference case; shorter live routes remain on CRAN.
  skip_on_cran()
  skip_if_not_installed("RTMB", "2.0")
  d <- shared_rater_example(240, 6, 0, 9236206)
  fit <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:2, quad_points = 61, person_sd = 1)
  expect_true(fit$checks$NumericalReady)
  expect_true(fit$checks$EstimatedVarianceBoundary)
  expect_identical(fit$zero_optimization$convergence, 0L)
  expect_lt(max(abs(fit$zero_optimization$gradient)), 1e-4)
  expect_true(all(c("original_code", "polished") %in% names(fit$zero_optimization)))
})

test_that("a separated fixed facet cannot pass on a small gradient alone", {
  skip_if_not_installed("RTMB", "2.0")
  # Reproduce the Windows failure mode: roundoff can leave a tiny positive
  # information matrix even though the likelihood has no finite maximum.
  local_mocked_bindings(optimHess = function(par, ...) diag(1e-11, length(par)),
    .package = "mfrmr")
  d <- expand.grid(Person = 1:24, Rater = 1:3, Criterion = 1:2)
  d$Score <- d$Criterion - 1
  expect_warning(fit <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion",
    0:1, rater_sd = .4, person_sd = 1), "information")
  expect_true(all(eigen(fit$information, symmetric = TRUE)$values > 0))
  expect_true(fit$checks$FixedFacetSeparation)
  expect_false(fit$checks$NumericalReady)
  expect_false(fit$checks$InformationPositive)
  expect_true(all(is.na(fit$calibration_table$SE)))
  expect_error(predict(fit, data.frame(Rater = "New", Criterion = "1"), 0, rater = "new"), "information")
  # A free-SD fit with unresolved information is also ineligible for profiling.
  fit$settings$fixed_rater_sd <- NULL
  expect_error(confint(fit), "information")
})

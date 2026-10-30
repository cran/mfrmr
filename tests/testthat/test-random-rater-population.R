# Independent population checks use small crossed designs, not archived fits.
random_population_data <- function() {
  withr::with_seed(923811, {
    d <- expand.grid(Person = 1:80, Rater = 1:4, Criterion = 1:2)
    theta <- rnorm(80, sd = 1.25)
    eta <- theta[d$Person] - c(-.8, -.3, .2, .9)[d$Rater] - c(-.2, .2)[d$Criterion]
    d$Score <- rbinom(nrow(d), 1, plogis(eta - .15))
    d
  })
}

test_that("nonunit ability likelihood, derivatives and Laplace boundary have independent references", {
  skip_if_not_installed("RTMB", "2.0")
  d <- expand.grid(Person = 1:3, Rater = 1:2)
  d$Score <- c(0, 1, 0, 1, 1, 0)
  input <- mfrm_random_rater_data(d, "Person", "Rater", character(), "Score", 0:1)
  # Adaptive normal integration and a separate optimizer/Hessian, no RTMB
  # quadrature or backend probability helper in this reference.
  direct <- function(par) {
    step <- par[1]; z <- par[2:3]; sd <- par[4]; ps <- par[5]
    ll <- sum(dnorm(z, log = TRUE))
    for (p in 1:3) {
      rows <- which(input$person == p)
      probability <- function(a) vapply(a, function(theta) {
        pr <- plogis(ps * theta - sd * z[input$rater[rows]] - step)
        prod(ifelse(input$y[rows] == 1, pr, 1 - pr)) * dnorm(theta)
      }, numeric(1))
      ll <- ll + log(integrate(probability, -Inf, Inf, rel.tol = 1e-11)$value)
    }
    -ll
  }
  joint <- mfrm_random_rater_objective(input, 61, random = FALSE, fixed_person_sd = NULL)
  par <- c(.2, .3, -.4, .7, 1.4)
  expect_equal(as.numeric(joint$fn(par)), direct(par), tolerance = 1e-9)
  grad <- vapply(seq_along(par), function(i) {
    h <- rep(0, length(par)); h[i] <- 1e-5
    (direct(par + h) - direct(par - h)) / 2e-5
  }, numeric(1))
  expect_equal(as.vector(joint$gr(par)), grad, tolerance = 1e-7)
  laplace <- function(v, sd = .7) {
    fn <- function(z) direct(c(.2, z, sd, sqrt(v)))
    opt <- optim(c(0, 0), fn, method = "BFGS", control = list(reltol = 1e-13))
    H <- optimHess(opt$par, fn)
    opt$value + as.numeric(determinant(H, logarithm = TRUE)$modulus) / 2 - log(2 * pi)
  }
  obj <- mfrm_random_rater_objective(input, 61, fixed_sd = .7, fixed_person_sd = 1.4)
  expect_equal(as.numeric(obj$fn(.2)), laplace(1.4^2), tolerance = 1e-6)
  # Boundary derivative of the whole Laplace likelihood, including its determinant.
  h <- .002
  reference <- 2 * (laplace(0) - laplace(h / 2)) / (h / 2) - (laplace(0) - laplace(h)) / h
  score <- mfrm_random_rater_person_zero_score(input, numeric(), .2, .7, 31)
  expect_equal(unname(score["score"]), reference, tolerance = 2e-4)
  expect_lt(score["difference"], 1e-4)
  refined <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", score_levels = 0:1,
    rater_sd = .7, person_sd = 1.4, quad_points = 181)
  expect_true(refined$checks$NumericalReady)
  expect_equal(refined$checks$CheckPoints, 363)
  expect_error(fit_mfrm_random_rater(d, "Person", "Rater", "Score", score_levels = 0:1,
    quad_points = 242), "7 to 241")
  # The other variance face retains the nonunit Person distribution.
  zero <- mfrm_random_rater_objective(input, 61, fixed_sd = 0, fixed_person_sd = 1.4)
  a <- mfrm_random_rater_objective(input, 61, fixed_sd = sqrt(1e-4), fixed_person_sd = 1.4)
  b <- mfrm_random_rater_objective(input, 61, fixed_sd = sqrt(5e-5), fixed_person_sd = 1.4)
  slope <- 2 * (zero$fn(.2) - b$fn(.2)) / 5e-5 - (zero$fn(.2) - a$fn(.2)) / 1e-4
  expect_equal(mfrm_random_rater_zero_score(input, numeric(), .2, 61, 1.4), as.numeric(slope), tolerance = 1e-5)
})

test_that("estimated population enters generalized rater uncertainty, profiling and generation", {
  skip_if_not_installed("RTMB", "2.0")
  fit <- fit_mfrm_random_rater(random_population_data(), "Person", "Rater", "Score", "Criterion", 0:1, quad_points = 61)
  expect_true(fit$checks$NumericalReady)
  expect_true(fit$checks$InformationPositive)
  expect_null(fit$settings$fixed_person_sd)
  expect_gt(fit$calibration$person_sd, 1)
  expect_equal(fit$calibration$person_variance, fit$calibration$person_sd^2)
  # Joint-Hessian implicit differentiation includes the ability-population coordinate.
  joint <- mfrm_random_rater_objective(fit$input, 61, random = FALSE, fixed_person_sd = NULL)
  sd <- fit$calibration$rater_sd; z <- fit$rater_mode / sd
  par <- c(fit$calibration$beta, fit$calibration$steps, z, sd, fit$calibration$person_sd)
  H <- optimHess(par, joint$fn, joint$gr)
  iz <- 2 + seq_along(z); ip <- setdiff(seq_along(par), iz)
  conditional <- sd^2 * solve(H[iz, iz])
  J <- -sd * solve(H[iz, iz], H[iz, ip]); J[, 3] <- J[, 3] + z
  expect_equal(conditional + J %*% fit$covariance %*% t(J), fit$rater_covariance,
    tolerance = 1e-5, ignore_attr = TRUE)
  profile <- confint(fit, level = .8)
  tab <- attr(profile, "profile")
  expect_true(all(tab$NumericalReady))
  expect_gt(diff(range(tab$PersonSD)), .001)
  expect_equal(vapply(profile[1, ], function(sd) 2 * (tab$NLL[which.min(abs(tab$SD - sd))] + fit$loglik), numeric(1)),
    rep(qchisq(.8, 1), 2), tolerance = 2e-4, ignore_attr = TRUE)
  generated <- mfrm_random_rater_generate(fit, 92815)
  expected <- withr::with_seed(92815, rnorm(80, sd = fit$calibration$person_sd))
  expect_identical(generated$theta, expected)
  # A real refit must re-estimate, rather than silently fix, that generating SD.
  bs <- mfrm_random_rater_intervals(fit, nsim = 2, seed = 923835)
  expect_true(all(bs$trials$FitReady))
  expect_true(any(abs(bs$trials$PersonSD - fit$calibration$person_sd) > .01))
  expect_null(bs$settings$fixed_person_sd)
  expect_equal(bs$settings$person_sd, fit$calibration$person_sd)
  expect_equal(bs$error, bs$truth - bs$estimates)
  rows <- data.frame(Rater = "new", Criterion = "1")
  pred <- predict(fit, rows, 1, rater = "new")
  # The supplied value is one logit, not one SD.
  beta <- as.vector(fit$input$basis$Criterion %*% fit$calibration$beta)[1]
  reference <- integrate(function(u) plogis(1 - beta - u - fit$calibration$steps) *
    dnorm(u, sd = sd), -Inf, Inf, rel.tol = 1e-10)$value
  expect_equal(unname(pred$probabilities[1, 2]), reference, tolerance = 1e-8)
  res <- mfrm_results(fit, predictions = pred, intervals = bs)
  expect_equal(subset(res$tables$variance, Effect == "Person ability")$Variance, fit$calibration$person_variance)
  expect_true(subset(res$tables$variance, Effect == "Person ability")$Estimated)
  f <- tempfile(); saveRDS(fit, f)
  expect_identical(predict(readRDS(f), rows, 1, rater = "new"), pred)
  changed <- fit; changed$calibration$person_sd <- 1
  expect_error(mfrm_results(changed, predictions = pred), "matching source")
})

test_that("zero-rater free-ability fit matches a separately integrated binary model", {
  skip_if_not_installed("RTMB", "2.0")
  d <- random_population_data()
  fit <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:1,
    rater_sd = 0, quad_points = 61)
  input <- fit$input
  direct <- function(par) {
    -sum(vapply(split(seq_len(nrow(d)), input$person), function(rows) {
      log(integrate(function(a) vapply(a, function(t) {
        p <- plogis(par[3] * t - input$X[rows, 1] * par[1] - par[2])
        prod(ifelse(input$y[rows] == 1, p, 1 - p)) * dnorm(t)
      }, numeric(1)), -Inf, Inf, rel.tol = 1e-9)$value)
    }, numeric(1)))
  }
  opt <- nlminb(c(0, 0, 1), direct, lower = c(-Inf, -Inf, 0),
    control = list(rel.tol = 1e-11))
  expect_true(fit$checks$NumericalReady)
  expect_equal(unname(fit$coefficients), unname(opt$par), tolerance = 2e-5)
  expect_equal(fit$loglik, -opt$objective, tolerance = 1e-7)
  expect_equal(unname(fit$covariance), unname(solve(optimHess(opt$par, direct))), tolerance = 2e-5)
})

test_that("ability boundaries, known populations and older saved populations remain explicit", {
  skip_if_not_installed("RTMB", "2.0")
  d <- expand.grid(Person = 1:12, Rater = 1:4, Criterion = 1:3)
  d$Score <- (d$Person + d$Criterion) %% 3
  fit <- fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:2)
  expect_true(fit$checks$NumericalReady)
  expect_true(fit$checks$EstimatedPersonVarianceBoundary)
  expect_true(fit$checks$EstimatedVarianceBoundary)
  expect_equal(fit$calibration$person_sd, 0)
  expect_lte(fit$checks$PersonZeroVarianceScore, 0)
  expect_lt(fit$checks$PersonZeroScoreDifference, 1e-4)
  expect_true(all(is.na(fit$raters$PredictionSE)))
  expect_true(all(is.na(fit$calibration_table$SE)))
  expect_error(confint(fit), "Person variance boundary")
  expect_error(mfrm_random_rater_intervals(fit, nsim = 2, seed = 1), "boundary")
  known <- fit_mfrm_random_rater(random_population_data(), "Person", "Rater", "Score", "Criterion", 0:1,
    person_sd = .7, quad_points = 61)
  expect_true(known$checks$NumericalReady)
  expect_identical(known$calibration$person_sd, .7)
  expect_equal(subset(mfrm_results(known)$tables$variance, Effect == "Person ability")$Estimated, FALSE)
  old <- known; old$calibration$person_sd <- old$calibration$person_variance <- NULL
  old$settings$fixed_person_sd <- NULL
  expect_identical(mfrm_random_rater_fixed_person_sd(old), 1)
  expect_equal(subset(mfrm_results(old)$tables$variance, Effect == "Person ability")$Variance, 1)
  expect_identical(mfrm_random_rater_generate(old, 7)$theta, withr::with_seed(7, rnorm(80)))
  for (sd in list(0, -1, Inf, NA, 1e-300, 1e200, "1"))
    expect_error(fit_mfrm_random_rater(d, "Person", "Rater", "Score", "Criterion", 0:2, person_sd = sd), "person_sd")
  expect_warning(capped <- fit_mfrm_random_rater(random_population_data(), "Person", "Rater", "Score", "Criterion",
    0:1, rater_sd = .4, person_variance_max = .01), "checks")
  expect_true(capped$checks$PersonVarianceUpperBoundary)
  expect_false(capped$checks$NumericalReady)
  expect_error(predict(capped, data.frame(Rater = "1", Criterion = "1"), 0), "numerical")
})

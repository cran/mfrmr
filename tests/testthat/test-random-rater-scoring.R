# A fixed calibration reference, not evidence that six observations estimate
# a rater population well. Actual fitted tutorial data exercise that workflow.
random_scoring_fixture <- function(sd = .7, ps = 1.2) {
  d <- expand.grid(Person = c("A", "B", "C"), Rater = c("R1", "R2"))
  d$Score <- c(0, 1, 2, 1, 0, 2)
  input <- mfrm_random_rater_data(d, "Person", "Rater", character(), "Score", 0:2)
  structure(list(input = input, calibration = list(beta = numeric(), steps = c(-.6, .6),
    rater_sd = sd, person_sd = ps, person_variance = ps^2),
    calibration_table = data.frame(Parameter = "Step", Facet = "Score", Level = c("1", "2"),
      Estimate = c(-.6, .6), SE = NA_real_, Lower = NA_real_, Upper = NA_real_),
    raters = data.frame(Rater = c("R1", "R2"), Estimate = 0, PredictionSE = NA_real_,
      Lower = NA_real_, Upper = NA_real_),
    checks = data.frame(NumericalReady = TRUE, InformationPositive = TRUE,
      EstimatedVarianceBoundary = FALSE),
    settings = list(method = "Fixed reference calibration", quad_points = 41L,
      fixed_rater_sd = sd, fixed_person_sd = ps, missing = "fail")), class = "mfrm_random_rater")
}

# Direct normalized adjacent-category probabilities, outside the AD backend.
random_scoring_probability <- function(theta, u, y) {
  eta <- theta - u
  lw <- cbind(0, eta + .6, 2 * eta)
  w <- exp(lw - apply(lw, 1L, max))
  (w / rowSums(w))[, y + 1]
}

test_that("focal ability likelihood keeps the joint shared-rater integral", {
  skip_if_not_installed("RTMB", "2.0")
  f <- random_scoring_fixture(); input <- f$input
  direct <- function(z, theta) {
    likelihood <- vapply(1:3, function(p) {
      rows <- which(input$person == p)
      conditional <- function(a) vapply(a, function(t) {
        prod(vapply(rows, function(i) random_scoring_probability(t,
          .7 * z[input$rater[i]], input$y[i]), numeric(1)))
      }, numeric(1))
      if (p == 1) conditional(theta) else integrate(function(a)
        conditional(1.2 * a) * dnorm(a), -Inf, Inf, rel.tol = 1e-10)$value
    }, numeric(1))
    -sum(log(likelihood)) - sum(dnorm(z, log = TRUE))
  }
  joint <- mfrm_random_rater_objective(input, 61, fixed_sd = .7,
    fixed_person_sd = 1.2, focal_person = 1, random = FALSE)
  expect_equal(as.numeric(joint$fn(c(-.6, .6, .2, -.3, .5))), direct(c(.2, -.3), .5), tolerance = 1e-8)
  obj <- mfrm_random_rater_objective(input, 61, fixed_sd = .7,
    fixed_person_sd = 1.2, focal_person = 1)
  for (theta in c(-1, .5, 2)) {
    opt <- optim(c(0, 0), direct, theta = theta, method = "BFGS", control = list(reltol = 1e-12))
    H <- optimHess(opt$par, function(z) direct(z, theta))
    laplace <- opt$value + as.numeric(determinant(H, logarithm = TRUE)$modulus) / 2 - log(2 * pi)
    expect_equal(as.numeric(obj$fn(c(-.6, .6, theta))), laplace, tolerance = 2e-6)
  }
  # Interpolation must agree with integration of the actual conditional Laplace
  # density, rather than validating itself with the same interpolated values.
  raw <- function(z) vapply(z, function(a)
    -as.numeric(obj$fn(c(-.6, .6, 1.2 * a))), numeric(1))
  reference <- mfrm_person_posterior_interval(raw, 1.2, .95)
  scored <- mfrm_random_rater_person_interval(input, f$calibration, 1, 61, .95)
  expect_equal(scored$value, reference, tolerance = 1e-6)
})

test_that("known zero raters give independently integrated normal-population scores", {
  skip_if_not_installed("RTMB", "2.0")
  f <- random_scoring_fixture(sd = 0, ps = .7)
  s <- score_mfrm_random_rater(f, persons = "A")
  expect_identical(s$table$Status, "available_conditional")
  density <- function(t) random_scoring_probability(t, 0, 0) *
    random_scoring_probability(t, 0, 1) * dnorm(t, sd = .7)
  mass <- integrate(density, -Inf, Inf, rel.tol = 1e-10)$value
  mean <- integrate(function(t) t * density(t), -10, 10, rel.tol = 1e-10)$value / mass
  variance <- integrate(function(t) (t - mean)^2 * density(t), -10, 10, rel.tol = 1e-10)$value / mass
  expect_equal(s$table$Estimate, mean, tolerance = 1e-7)
  expect_equal(s$table$ConditionalSD, sqrt(variance), tolerance = 1e-7)
  expect_equal(integrate(density, -Inf, s$table$Lower, rel.tol = 1e-10)$value / mass, .025, tolerance = 1e-7)
  expect_equal(integrate(density, s$table$Upper, Inf, rel.tol = 1e-10)$value / mass, .025, tolerance = 1e-7)
})

test_that("output selection retains other Persons and shared evidence is counted once", {
  skip_if_not_installed("RTMB", "2.0")
  f <- random_scoring_fixture()
  a <- score_mfrm_random_rater(f, persons = "A")
  all <- score_mfrm_random_rater(f)
  expect_equal(a$table, all$table[all$table$Person == "A", ], ignore_attr = TRUE)
  expect_identical(a$scoring_data, f$input$assigned_data)
  d <- f$input$assigned_data
  changed <- d; changed$Score[changed$Person != "A"] <- 2
  other <- score_mfrm_random_rater(f, changed, persons = "A")
  expect_identical(other$table$Status, "available_conditional")
  expect_gt(abs(other$table$Estimate - a$table$Estimate), .03)
  # Subsetting the roster is deliberately different from selecting output IDs.
  alone <- score_mfrm_random_rater(f, d[d$Person == "A", ])
  expect_gt(abs(alone$table$ConditionalSD - a$table$ConditionalSD), .001)
  # All identical rater IDs are shared; labels carry no numeric meaning.
  relabeled <- d; relabeled$Rater <- ifelse(d$Rater == "R1", "z", "a")
  renamed <- score_mfrm_random_rater(f, relabeled, persons = "A")
  expect_equal(renamed$table[1, 4:7], a$table[1, 4:7], tolerance = 1e-7)
  # Newdata replaces rather than appends the complete roster.
  expect_equal(score_mfrm_random_rater(f, d, persons = "A")$table, a$table)
  old <- f; old$input$assigned_data <- NULL
  old$calibration$person_sd <- old$calibration$person_variance <- NULL
  old$settings$fixed_person_sd <- NULL
  fixed <- f; fixed$calibration$person_sd <- 1; fixed$calibration$person_variance <- 1
  expect_equal(score_mfrm_random_rater(old, persons = "A")$table,
    score_mfrm_random_rater(fixed, persons = "A")$table)
})

test_that("prior-only, zero-population and failed scoring rows stay explicit", {
  skip_if_not_installed("RTMB", "2.0")
  f <- random_scoring_fixture()
  d <- f$input$assigned_data; extra <- d[1, ]; extra$Person <- "unobserved"; extra$Score <- NA
  d <- rbind(d, extra)
  expect_error(score_mfrm_random_rater(f, d), "missing")
  s <- score_mfrm_random_rater(f, d, persons = c("unobserved", "A"), missing = "omit")
  expect_identical(s$table$Status, c("prior_only", "available_conditional"))
  expect_equal(s$table$ConditionalSD[1], 1.2)
  expect_equal(s$table$Lower[1], 1.2 * qnorm(.025))
  expect_identical(s$omitted_rows, 7L)
  expect_equal(unname(s$data_usage[c("Input", "Observed", "RosterPersons", "RequestedPersons")]), c(7, 6, 4, 2))
  zero <- f; zero$calibration$person_sd <- zero$calibration$person_variance <- 0
  z <- score_mfrm_random_rater(zero, persons = "A")
  expect_identical(z$table$Status, "unavailable"); expect_true(is.na(z$table$Estimate))
  expect_match(z$table$Reason, "ability variance is zero")
  expect_error(score_mfrm_random_rater(f, persons = "absent"), "Person IDs")
  expect_error(score_mfrm_random_rater(f, persons = c("A", "A")), "distinct")
  expect_error(score_mfrm_random_rater(f, quad_points = 242), "quad_points")
  old <- f; old$input$assigned_data <- NULL; old$input$omitted_rows <- 7L
  expect_error(score_mfrm_random_rater(old), "missing-score roster")
  local_mocked_bindings(mfrm_random_rater_person_interval = function(...) stop("Unresolved inner integral"))
  bad <- score_mfrm_random_rater(f, persons = c("C", "A"))
  expect_identical(bad$table$Person, c("C", "A"))
  expect_true(all(bad$table$Status == "unavailable"))
  expect_true(all(grepl("Unresolved inner integral", bad$table$Reason)))
})

test_that("saved Person scoring connects to plots, reports and source guards", {
  skip_if_not_installed("RTMB", "2.0")
  f <- random_scoring_fixture()
  d <- f$input$assigned_data; d$Score[d$Person == "C"] <- NA
  s <- score_mfrm_random_rater(f, d, persons = c("C", "A"), missing = "omit")
  p <- predict(f, data.frame(Rater = "new"), ability = 0, rater = "new")
  res <- mfrm_results(f, predictions = p, scores = s)
  expect_equal(res$tables$person_scores, s$table)
  expect_equal(res$tables$scoring_roster, s$scoring_data)
  expect_equal(res$tables$category_probabilities$Probability, as.vector(p$probabilities))
  expect_true(res$plot_map$Available[res$plot_map$Type == "scores"])
  expect_equal(plot_data(plot(res, type = "scores", draw = FALSE))$table, s$table)
  expect_equal(summary(s)$status_counts[c("available_conditional", "prior_only")], c(available_conditional = 1L, prior_only = 1L), ignore_attr = TRUE)
  expect_match(mfrm_report(res)$markdown, "joint conditional rater Laplace")
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(as_ggplot(plot(s, draw = FALSE)), "ggplot")
  tmp <- tempfile(); saveRDS(res, tmp); expect_identical(readRDS(tmp), res)
  different <- f; different$calibration$person_sd <- .8
  expect_error(mfrm_results(different, scores = s), "matching calibration")
  expect_error(mfrm_results(s), "saved Person scores")
  bad <- s; bad$table$Estimate <- NA_real_; bad$table$Lower <- bad$table$Upper <- NA_real_
  bad$table$Status <- "unavailable"
  before <- grDevices::dev.cur(); plot(bad, draw = FALSE); expect_identical(grDevices::dev.cur(), before)
  grDevices::pdf(tempfile(fileext = ".pdf")); expect_silent(plot(bad)); grDevices::dev.off()
})

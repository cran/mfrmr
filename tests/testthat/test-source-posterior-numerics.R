test_that("continuous source scoring preserves a normal prior under constant likelihoods", {
  cuts <- c(-1, 0, 1)
  expected <- c(Estimate = 0, ConditionalSD = 1.7,
    Lower = 1.7 * qnorm(.025), Upper = 1.7 * qnorm(.975))
  for (offset in c(0, -1e12, -1e16)) {
    out <- mfrm_person_posterior_interval(function(z) rep(offset, length(z)),
      1.7, .95, cuts)
    expect_equal(out[1:4], expected, tolerance = 1e-9)
    expect_equal(unname(out["LogNormalizer"]), offset, tolerance = 1e-12)
    expect_equal(unname(out[paste0("Above", 1:3)]), pnorm(cuts, lower.tail = FALSE), tolerance = 1e-10)
    expect_equal(unname(out[paste0("Band", 1:4)]), diff(pnorm(c(-Inf, cuts, Inf))), tolerance = 1e-10)
  }
  expect_error(mfrm_person_posterior_interval(function(z) rep(NA_real_, length(z)), 1, .95), "nonfinite")
  expect_error(mfrm_person_posterior_interval(function(z) 30 * z, 1, .95), "outside")
})

test_that("testlet and shared-rater consumers center likelihoods before adding the prior", {
  local_mocked_bindings(
    mfrm_testlet_kernel = function(input, par, theta, rule, ...) list(
      loglik = matrix(-1e16, length(theta), 1L),
      score = array(0, c(length(theta), 1L, 2L)),
      theta_score = matrix(0, length(theta), 1L), theta_curvature = matrix(0, length(theta), 1L)),
    mfrm_random_rater_objective = function(...) list(fn = function(...) 1e16), .package = "mfrmr")
  input <- list(X = matrix(numeric(), 1L, 0L), score_levels = 0:1)
  expected <- c(Estimate = 0, ConditionalSD = 1, Lower = qnorm(.025), Upper = qnorm(.975))
  rule <- gauss_hermite_normal(7L)
  testlet <- mfrm_testlet_person_interval(input, c(0, 0, 1), 1L, rule, .95)
  expect_equal(testlet[1:4], expected, tolerance = 1e-9)
  moments <- mfrm_testlet_evaluate(input, c(0, 0, 1), rule)
  expect_equal(unname(moments$moments), matrix(c(0, 1), 1L), tolerance = 1e-12)
  expect_equal(moments$gradient, rep(0, 3L))
  # A weighted SD remains nonzero even when E(theta^2)-E(theta)^2 cancels.
  narrow <- list(nodes = c(1, 1 + 1e-8), weights = c(.5, .5))
  moments <- mfrm_testlet_evaluate(input, c(0, 0, 1), narrow)$moments
  expect_equal(unname(moments[1, "ConditionalSD"]), diff(narrow$nodes) / 2, tolerance = 1e-12)
  random <- mfrm_random_rater_person_interval(input,
    list(beta = numeric(), steps = 0, rater_sd = 0, person_sd = 1), 1L, 7L, .95)
  expect_equal(random$value[1:4], expected, tolerance = 1e-9)
  expect_lt(random$difference, 1e-8)
})

test_that("testlet finite category ties have the correct likelihood and scalar-node derivatives", {
  d <- expand.grid(Repeat = 1:2, Block = c("X", "Y"), Person = c("A", "B", "C"))
  d$Score <- ifelse(d$Person == "C", 0L, d$Repeat)
  input <- mfrm_testlet_data(d, "Person", "Score", "Block", character(), 0:2, "fail")
  for (big in c(1e12, 1e16)) for (n in c(1L, 2L)) {
    out <- mfrm_testlet_kernel(input, c(-big/2, big/2, 0), rep(big/2, n),
      gauss_hermite_normal(7L), gradient = TRUE, persons = 1:2)
    expect_equal(out$loglik, matrix(-4 * log(2), n, 2L), tolerance = 1e-12)
    expect_equal(out$theta_score, matrix(0, n, 2L))
    expect_equal(out$theta_curvature, matrix(-1, n, 2L))
    expect_equal(as.numeric(out$score[, , 1:2]), rep(0, 4L * n))
    expect_equal(as.numeric(out$score[, , 3L]), rep(-.5, 2L * n))
  }
  # A single latent point occurs at zero variance; it must work away from ties too.
  out <- mfrm_testlet_kernel(input, c(-.5, .5, 0), 0, gauss_hermite_normal(7L), gradient = TRUE)
  direct <- function(par) sum(mfrm_testlet_kernel(input, par, 0, gauss_hermite_normal(7L))$loglik)
  for (k in 1:2) {
    delta <- numeric(3); delta[k] <- 1e-5
    expect_equal(sum(out$score[, , k]), (direct(c(-.5,.5,0)+delta)-direct(c(-.5,.5,0)-delta))/2e-5,
      tolerance = 1e-8)
  }
  diag <- mfrm_testlet_response_probabilities(input, c(-1e16, 0, 0, 0), 1:4, 7L)
  expect_equal(diag$probabilities, matrix(c(0, .5, .5), 4L, 3L, byrow = TRUE), tolerance = 1e-12)
  ordinary <- mfrm_ordinary_response_probabilities(list(mean = 0, sd = 0,
    person = c("A", "A"), y = c(1L, 2L), offset = c(0, 0), steps = c(-1e16, 0), score_levels = 0:2), 1:2, 7L)
  expect_equal(ordinary$probabilities, matrix(c(0, .5, .5), 2L, 3L, byrow = TRUE), tolerance = 1e-12)
})

test_that("shared-rater AD probabilities preserve ties and adjacent-step differences", {
  skip_if_not_installed("RTMB", "2.0")
  d <- expand.grid(Person = c("A", "B", "C"), Rater = c("r1", "r2"))
  d$Score <- c(0, 1, 2, 1, 2, 0)
  reference <- mfrm_random_rater_data(d, "Person", "Rater", character(), "Score", 0:2)
  d$Score <- pmax(1L, d$Score)
  input <- mfrm_random_rater_data(d, "Person", "Rater", character(), "Score", 0:2, reference = reference)
  joint <- mfrm_random_rater_objective(input, 7L, fixed_sd = 0, fixed_person_sd = 0, random = FALSE)
  for (step in c(0, .1, 1)) {
    par <- c(-1e16, step, 0, 0)
    probability <- plogis(-step)
    expected <- -4 * log1p(-probability) - 2 * log(probability) + log(2*pi)
    expect_equal(as.numeric(joint$fn(par)), expected, tolerance = 1e-12)
    expect_equal(as.vector(joint$gr(par)), c(0, 2 - 6 * probability, 0, 0), tolerance = 1e-12)
  }
})

test_that("shared-rater replicate normalization does not hide unresolved likelihood totals", {
  local_mocked_bindings(mfrm_random_rater_objective = function(...) list(fn = function(...) 1e16),
    .package = "mfrmr")
  out <- mfrm_random_rater_response_probabilities(list(),
    list(beta = numeric(), steps = c(-.5, .5), rater_sd = 0), 1L, 7L)
  expect_equal(out$probabilities, matrix(1/3, 1L, 3L), tolerance = 1e-12)
  expect_equal(out$normalization_error, 2, tolerance = 1e-12)
})

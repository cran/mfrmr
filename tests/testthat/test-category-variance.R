test_that("category variance preserves tail information and score reflection", {
  eta <- c(-40, -20, -1, 0, 1, 20, 40)
  p <- mfrmr:::category_prob_rsm(eta, c(0, 0))
  q <- exp(-abs(eta))
  reference <- q / (1 + q)^2
  variance <- mfrmr:::mfrm_category_variance(p)
  # Relative errors matter: an absolute tolerance would accept zero tail variance.
  expect_equal(variance / reference, rep(1, length(eta)), tolerance = 1e-12)
  expect_equal(mfrmr:::mfrm_category_variance(p[, 2:1]) / variance,
    rep(1, length(eta)), tolerance = 1e-12)
  expect_equal(mfrmr:::mfrm_category_variance(p, c(5, 6)) / variance,
    rep(1, length(eta)), tolerance = 1e-12)
  expect_equal(mfrmr:::mfrm_category_variance(matrix(c(0, 1), 1)), 0)
  expect_length(mfrmr:::mfrm_category_variance(matrix(numeric(), 0, 3)), 0)

  p <- mfrmr:::category_prob_rsm(eta, c(0, -.4, .7, 0))
  # Independent pairwise identity: sum_{i<j} p_i p_j (i-j)^2.
  reference <- Reduce(`+`, lapply(combn(1:4, 2, simplify = FALSE), function(ij)
    p[, ij[1]] * p[, ij[2]] * diff(ij)^2))
  expect_equal(mfrmr:::mfrm_category_variance(p) / reference,
    rep(1, length(eta)), tolerance = 1e-12)
})

test_that("R and compiled posterior category moments preserve tail variance", {
  p <- mfrmr:::category_prob_rsm(c(-40, 40), c(0, 0))
  logprob <- list(prob_list = list(p, p))
  posterior <- list(obs_posterior = matrix(c(.3, .7), 2, 2, byrow = TRUE))
  reference <- exp(-40) / (1 + exp(-40))^2
  r <- mfrmr:::mfrm_mml_expected_category_bundle_r(logprob, posterior)
  expect_equal(r$var_k / reference, c(1, 1), tolerance = 1e-12)
  expect_true(mfrmr:::mfrm_cpp11_backend_available())
  cpp <- mfrmr:::mfrm_cpp_expected_category_bundle(logprob$prob_list,
    posterior$obs_posterior, FALSE)
  expect_equal(cpp$var_k / reference, c(1, 1), tolerance = 1e-12)
  expect_equal(cpp$expected_k, r$expected_k)
})

test_that("public GPCM information routes agree in both tails", {
  fit <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit
  grid <- expand.grid(Theta = c(-40, 40),
    Rater = fit$prep$levels$Rater, Criterion = fit$prep$levels$Criterion)
  curves <- mfrm_curve_intervals(fit, grid, type = "information")
  probability <- mfrm_curve_intervals(fit, grid, type = "probability")
  reference <- vapply(split(probability$table$Estimate, probability$table$InputRow),
    function(p) sum(vapply(combn(seq_along(p), 2, simplify = FALSE),
      function(ij) p[ij[1]] * p[ij[2]] * diff(ij)^2, numeric(1))), numeric(1))
  reference <- reference[as.character(seq_len(nrow(grid)))]
  slope <- fit$slopes$Estimate[match(grid[[fit$config$slope_facet]], fit$slopes$SlopeFacet)]
  expect_equal(curves$table$Estimate / unname(reference * slope^2),
    rep(1, nrow(grid)), tolerance = 1e-10)

  data <- fit$prep$data
  exposure <- vapply(seq_len(nrow(grid)), function(i)
    sum(data$Weight[as.character(data$Rater) == grid$Rater[i] &
      as.character(data$Criterion) == grid$Criterion[i]]), numeric(1))
  total <- compute_information(fit, theta_range = c(-40, 40), theta_points = 2)$tif
  expected <- vapply(total$Theta, function(theta)
    sum(curves$table$Estimate[grid$Theta == theta] * exposure[grid$Theta == theta]), numeric(1))
  expect_true(all(total$Information > 0))
  expect_equal(as.numeric(total$Information / expected), c(1, 1), tolerance = 1e-10)

  for (model in c("RSM", "PCM", "GPCM")) {
    spec <- list(model = model, categories = 0:2,
      groups = list(A = list(name = "A", step_cum = c(0, 0, 0), slope = 1)))
    tables <- mfrmr:::build_curve_tables(spec, c(-40, 40))
    parts <- tapply(tables$probabilities$CategoryInformation,
      tables$probabilities$Theta, sum)
    expect_true(all(tables$expected$Information > 0))
    expect_equal(tables$expected$Information / as.numeric(parts), c(1, 1), tolerance = 1e-12)
  }
})

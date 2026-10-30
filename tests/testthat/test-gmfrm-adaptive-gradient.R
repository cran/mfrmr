test_that("two-family adaptive gradients include moving nodes in free coordinates", {
  for (maximum in c(1L, 3L)) {
    data <- expand.grid(Candidate = paste0("p", 1:4),
      A = c("same", "a2", "a3"), B = c("same", "b2", "b3"))
    data <- data[!(data$A == "a3" & data$B == "b2"), ]
    data$Rating <- rep(0:maximum, length.out = nrow(data))
    names(data)[2:3] <- c("First owner", "Second owner")
    problem <- mfrm_gmfrm_problem(data, maximum, gauss_hermite_normal(31L),
      c("First owner", "Second owner"), "Candidate", "Rating")
    setup <- problem$common
    par <- unname(problem$start + .3 * sin(seq_along(problem$start)))
    for (order in c(1L, 3L, 15L)) {
      evaluate <- mfrmr_make_adaptive_mml_evaluator(
        setup$idx, setup$config, setup$sizes, order)
      # Independently assemble the finite sum through the diagnostic path;
      # it does not call the gradient evaluator. Low orders expose missing
      # derivatives of the posterior mode, width or coordinate Jacobian.
      objective <- function(x) {
        review <- mfrmr_adaptive_quadrature_review(setup$idx, setup$config,
          expand_params(x, setup$sizes, setup$config),
          problem$specification$quadrature, setup$prep$levels$Person, order)
        stopifnot(all(review$Status == "computed"))
        -sum(review$AdaptiveLogMarginal)
      }
      difference <- function(h) vapply(seq_along(par), function(j) {
        plus <- minus <- par
        plus[j] <- plus[j] + h
        minus[j] <- minus[j] - h
        (objective(plus) - objective(minus)) / (2 * h)
      }, 0)
      result <- evaluate(par)
      numeric <- (4 * difference(5e-5) - difference(1e-4)) / 3
      expect_lt(abs(result$value - objective(par)), 1e-10)
      expect_lt(max(abs(result$gradient - numeric)), 1e-6)
    }
    # At high order the moving-node information must also approach the
    # independent complete-data-minus-missing-information calculation.
    cfg <- setup$config; cfg$estimation_control <- list(mml_integration = "adaptive")
    params <- expand_params(par, setup$sizes, cfg)
    basis <- mfrmr_adaptive_quadrature_basis(setup$idx, cfg, params,
      gauss_hermite_normal(61L), compute_base_eta(setup$idx, params, cfg))
    reference <- gmfrm_louis_reference(problem$specification, par, basis$nodes, basis$log_weights)
    evaluate <- mfrmr_make_adaptive_mml_evaluator(setup$idx, cfg, setup$sizes, 61L)
    hessian <- optimHess(par, function(x) evaluate(x)$value, function(x) evaluate(x)$gradient,
      control = list(ndeps = rep(1e-4, length(par))))
    expect_equal(hessian, reference$information, tolerance = 1e-6)
    expect_equal(evaluate(par)$gradient, -colSums(reference$scores), tolerance = 1e-7)
  }
})

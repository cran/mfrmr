test_that("high-order normal quadrature retains relatively accurate tail weights", {
  # 100-digit mpmath Hermite roots and derivative weights; see the repository
  # gauss-hermite-weight-repair-0.2.4 validation record for full-grid comparison.
  orders <- c(31L, 61L, 121L, 181L, 301L)
  tail_weights <- c(2.6059738548930083e-22, 9.371228767883570e-47,
                    7.899361179421854e-97, 1.646586394927613e-147,
                    1.242297899579240e-249)
  for (i in seq_along(orders)) {
    q <- mfrmr:::gauss_hermite_normal(orders[i])
    expect_true(all(is.finite(q$weights) & q$weights > 0))
    expect_equal(sum(q$weights), 1, tolerance = 1e-14)
    expect_equal(q$nodes, -rev(q$nodes), tolerance = 1e-14)
    expect_equal(q$weights / rev(q$weights), rep(1, orders[i]), tolerance = 1e-12)
    expect_equal(q$weights[1] / tail_weights[i], 1, tolerance = 2e-10)
    expect_equal(vapply(1:4, function(k) sum(q$weights * q$nodes^(2*k)), 0),
                 c(1, 3, 15, 105), tolerance = 1e-12)
  }

  # A tilted normal puts material posterior mass on formerly discarded nodes.
  # Ordinary low-order normal moments alone cannot detect that loss.
  q <- mfrmr:::gauss_hermite_normal(181L)
  joint <- log(q$weights) + 14 * q$nodes - 14^2 / 2
  normalizer <- mfrmr:::logsumexp(joint)
  posterior <- exp(joint - normalizer)
  expect_equal(normalizer, 0, tolerance = 1e-10)
  expect_equal(sum(posterior * q$nodes), 14, tolerance = 1e-10)
  expect_equal(sum(posterior * (q$nodes - 14)^2), 1, tolerance = 1e-10)
})

test_that("quadrature refuses invalid orders and unrepresentable weights", {
  for (n in list(0, -1, 1.5, NA_real_, Inf, numeric(), c(1, 2), "31")) {
    expect_error(mfrmr:::gauss_hermite_normal(n), "integer n >= 1", fixed = TRUE)
  }
  expect_identical(mfrmr:::gauss_hermite_normal(1L), list(nodes = 0, weights = 1))
  expect_error(mfrmr:::gauss_hermite_normal(400L),
               "cannot all be represented", fixed = TRUE)
})

test_that("endpoint rank errors obey the exact binomial allocation", {
  for (n in c(2, 3, 49, 173, 174, 199, 999)) {
    for (level in c(.5, .8, .95, .99)) {
      m <- mfrm_fair_bootstrap_mc(0, matrix(seq_len(n), ncol = 1), level, .95)
      # Independent probability check, including ranks outside observed data.
      left <- pbinom(m$RankLower - 1, n, m$Probability)
      right <- pbinom(m$RankUpper - 1, n, m$Probability, lower.tail = FALSE)
      expect_true(all(left <= .0125 + 1e-14))
      expect_true(all(right <= .0125 + 1e-14))
      expect_lte(sum(left + right), .05 + 1e-14)
    }
  }
  small <- mfrm_fair_bootstrap_mc(0, matrix(1:173, ncol = 1), .95, .95)
  enough <- mfrm_fair_bootstrap_mc(0, matrix(1:174, ncol = 1), .95, .95)
  expect_true(all(small$RankLimited))
  expect_false(any(enough$RankLimited))
})

test_that("discrete quantile endpoint enclosures pass exact multinomial enumeration", {
  n <- 12L; level <- .5; mc_level <- .8; estimate <- 2
  counts <- expand.grid(a = 0:n, b = 0:n)
  counts <- counts[rowSums(counts) <= n, ]; counts$c <- n - rowSums(counts)
  for (prob in list(c(.25, .5, .25), c(.2, .5, .3))) {
    q <- vapply(c(.75, .25), function(p) (0:2)[which(cumsum(prob) >= p)[1]], 0)
    endpoints <- estimate - q
    covered <- weights <- numeric(nrow(counts))
    for (i in seq_len(nrow(counts))) {
      k <- as.integer(counts[i, ])
      roots <- matrix(rep(0:2, k), ncol = 1)
      m <- mfrm_fair_bootstrap_mc(estimate, roots, level, mc_level)
      covered[i] <- all(m$Lower <= endpoints & endpoints <= m$Upper)
      weights[i] <- dmultinom(k, prob = prob)
    }
    expect_equal(sum(weights), 1, tolerance = 1e-13)
    expect_gte(sum(weights * covered), mc_level)
    expect_lt(sum(weights * covered), 1)
  }
})

test_that("unknown refits enclose completions and are never removed from ranks", {
  roots <- cbind(first = seq(-1, 1, length.out = 199), second = rep(0:2, length.out = 199))
  roots[c(1, 70, 199), 1] <- NA_real_
  m <- mfrm_fair_bootstrap_mc(c(2, 3), roots, .95, .95)
  expect_equal(m$Planned, rep(199, 4))
  expect_equal(m$Unresolved, c(3, 3, 0, 0))
  expect_false(any(m$RankLimited))
  expect_identical(m$Lower[1], -Inf)
  expect_identical(m$Upper[2], Inf)
  for (replacement in c(-1e6, -.2, 1e6)) {
    complete <- roots; complete[!is.finite(complete)] <- replacement
    cmp <- mfrm_fair_bootstrap_mc(c(2, 3), complete, .95, .95)
    expect_true(all(m$Lower <= cmp$Lower))
    expect_true(all(m$Upper >= cmp$Upper))
    # Changing only target 1 cannot change target 2's uncertainty.
    expect_identical(m[3:4, ], cmp[3:4, ])
  }
  broad <- mfrm_fair_bootstrap_mc(c(2, 3), roots, .95, .99)
  expect_true(all(broad$Lower <= m$Lower))
  expect_true(all(broad$Upper >= m$Upper))
})

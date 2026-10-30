test_that("documented bootstrap tail limits retain unresolved roots as sample size grows", {
  # These deliberately constructed roots check the completion rule, not coverage.
  for (b in c(19L, 99L, 499L)) {
    cap <- switch(as.character(b), `19` = 0L, `99` = 2L, `499` = 12L)
    x <- structure(list(source = list(raters = data.frame(Rater = "R",
      Estimate = .2, PredictionSE = .5)),
      studentized = matrix(seq(-2, 2, length.out = b), ncol = 1),
      settings = list(nsim = b, level = .95)), class = "mfrm_random_rater_intervals")
    if (cap > 0) x$studentized[seq_len(cap), 1] <- NA_real_
    finite <- confint(x, level = .95)
    expect_true(all(is.finite(finite)))
    expect_equal(attr(finite, "availability")$Planned, b)
    x$studentized[cap + 1L, 1] <- NA_real_
    infinite <- confint(x, level = .95)
    expect_identical(unname(infinite[1, ]), c(-Inf, Inf))
    expect_equal(attr(infinite, "availability")$Unresolved, cap + 1L)
    # A persistent 3% unresolved fraction does not become finite at a larger B.
    x$studentized[seq_len(ceiling(.03 * b)), 1] <- NA_real_
    expect_identical(unname(confint(x, level = .95)[1, ]), c(-Inf, Inf))
  }
})

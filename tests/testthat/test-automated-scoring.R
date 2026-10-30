test_that("QWK agrees with a hand contingency calculation and independent MSE", {
  # Marginals (13,27,60), (13,24,63); observed squared disagreement .11,
  # independent-marginal disagreement 1.02, hence QWK = 91/102.
  cells <- expand.grid(a = 0:2, h = 0:2)
  freq <- c(13, 0, 0, 0, 20, 4, 0, 7, 56)
  a <- rep(cells$a, freq); h <- rep(cells$h, freq)
  x <- automated_scoring_table(a, h, h, categories = 0:2)
  expect_equal(x$table$QWK, 91 / 102, tolerance = 1e-14)
  expect_equal(x$table$HumanQWK, 1)
  expect_equal(x$table$QWKDegradation, 11 / 102)
  expect_equal(x$table$SMDPooled, (mean(a) - mean(h)) /
    sqrt((var(a) + var(h)) / 2))
  expect_equal(x$table$SMDHuman, (mean(a) - mean(h)) / sd(h))
  expect_identical(dim(x$counts$automated_human), c(3L, 3L))

  # An independent response-pair sum implements note 1 directly. Uneven
  # categories distinguish actual score distances from category-rank distances.
  points <- c(0, 1, 5)
  a <- points[c(1, 1, 2, 2, 3, 3, 3)]
  h <- points[c(1, 2, 1, 3, 2, 3, 3)]
  expected <- 1 - mean((a - h)^2) / mean(outer(a, h, "-")^2)
  x <- automated_scoring_table(a, h, categories = points)
  expect_equal(x$table$QWK, expected)
  expect_equal(automated_scoring_table(2 + 3 * a, 2 + 3 * h,
    categories = 2 + 3 * points)$table$QWK, expected)
  expect_equal(automated_scoring_table(a, h, categories = c(0, 1, 3, 5))$table$QWK, expected)
  expect_equal(automated_scoring_table(a, h, categories = c(0, 1, 5, 1e200))$table$QWK, expected)
  nearby <- c(1, 1 + .Machine$double.eps, 1 + 2 * .Machine$double.eps)
  z <- automated_scoring_table(nearby, nearby, categories = nearby)
  expect_equal(z$table$QWK, 1)
  expect_equal(unname(diag(z$counts$automated_human)), rep(1, 3))
})

test_that("rounding, SMD denominator and degradation have explicit meanings", {
  a <- c(.2, .9, 1.3, 2.1, 2.7, 2.8)
  ac <- c(0, 1, 1, 2, 3, 3)
  h <- c(0, 1, 1, 2, 3, 3); hr <- c(0, 1, 2, 2, 2, 3)
  x <- automated_scoring_table(a, h, hr, categories = 0:3, automated_categories = ac)
  expect_equal(x$table$QWK, 1)
  expect_equal(x$table$Correlation, cor(a, h))
  expect_equal(x$table$HumanCorrelation, cor(hr, h))
  expect_lt(x$table$QWKDegradation, 0)
  expect_equal(x$table$CorrelationDegradation, cor(hr, h) - cor(a, h))
  expect_equal(x$table$SMD, (mean(a) - mean(h)) / sqrt((var(a) + var(h)) / 2))
  y <- automated_scoring_table(a, h, hr, categories = 0:3,
    automated_categories = ac, smd = "human")
  expect_equal(y$table$SMD, (mean(a) - mean(h)) / sd(h))
  expect_identical(y$table$SMDBasis, "human")
  expect_identical(tail(y$criteria$Basis, 1), "unrounded_scores_smd_human")
  expect_error(automated_scoring_table(a, h, categories = 0:3), "outside")
  # Exactly at the configured limits is included; display rounding must not
  # rescue a score just below the limit.
  z <- automated_scoring_table(a, h, hr, categories = 0:3, automated_categories = ac,
    qwk_min = x$table$QWK, correlation_min = x$table$Correlation,
    max_abs_smd = abs(x$table$SMD))
  expect_true(all(z$criteria$Status[c(1, 2, 5)] == "within_threshold"))
  z <- automated_scoring_table(a, h, hr, categories = 0:3, automated_categories = ac,
    correlation_min = x$table$Correlation + .00001)
  expect_identical(z$criteria$Status[2], "review")
  expect_identical(z$settings$threshold_rounding, "none")
})

test_that("a single complete-case roster governs both baselines and comparisons", {
  a <- c(0, 1, 2, 3, NA, 1)
  h <- c(0, 2, 2, 3, 1, NA)
  hr <- c(0, 1, NA, 3, 2, 1)
  x <- automated_scoring_table(a, h, hr, categories = 0:3)
  y <- automated_scoring_table(a[c(1, 2, 4)], h[c(1, 2, 4)],
    hr[c(1, 2, 4)], categories = 0:3)
  expect_identical(x$cases$Included, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(x$table$NInput, 6)
  expect_equal(x$table$N, 3)
  expect_equal(x$table$NExcluded, 3)
  expect_identical(x$counts, y$counts)
  expect_equal(x$criteria, y$criteria)
  expect_equal(sum(x$counts$automated_human), sum(x$counts$human_human))
  ac <- a; ac[2] <- NA
  z <- automated_scoring_table(a, h, hr, categories = 0:3, automated_categories = ac)
  expect_equal(z$table$N, 2)
  expect_true(z$cases$MissingAutomatedCategory[2])
})

test_that("constant, missing and insufficient ratings remain unavailable", {
  x <- automated_scoring_table(rep(1, 3), rep(1, 3), categories = 0:2)
  expect_true(all(is.na(x$criteria$Value)))
  expect_true(all(x$criteria$Status == "unavailable"))
  expect_match(x$criteria$Reason[1], "zero")
  expect_match(x$criteria$Reason[3], "not supplied")
  expect_false(x$table$HumanReferenceSupplied)
  y <- automated_scoring_table(rep(2, 3), rep(1, 3), categories = 0:2)
  expect_equal(y$table$QWK, 0)
  expect_true(is.na(y$table$SMD))
  for (a in list(numeric(), NA_real_, 1, c(NA_real_, NA_real_))) {
    z <- automated_scoring_table(a, a, a, categories = 0:2)
    expect_true(all(z$criteria$Status == "unavailable"))
    expect_true(all(nzchar(z$criteria$Reason)))
  }
  z <- automated_scoring_table(0:3, 3:0, 0:3, categories = 0:3)
  expect_equal(z$table$QWK, -1)
  expect_equal(z$table$Correlation, -1)
  expect_identical(z$criteria$Status[1], "review")
  expect_error(automated_scoring_table(0:2, 0:1, categories = 0:2), "same length")
  expect_error(automated_scoring_table(c(0, Inf), 0:1, categories = 0:2), "finite")
  expect_error(automated_scoring_table(factor(0:2), 0:2, categories = 0:2), "numeric")
  expect_error(automated_scoring_table(0:2, 0:2, categories = c(0, 2, 1)), "increasing")
  expect_error(automated_scoring_table(0:2, 0:2, categories = 0:1), "outside")
  expect_error(automated_scoring_table(0:2, 0:2, categories = 0:2, max_abs_smd = -1), "threshold")
  expect_error(automated_scoring_table(0:2, 0:2, categories = 0:2, qwk_min = NA), "threshold")
})

test_that("optional MFRM values reuse matching saved diagnostics and retain scope", {
  fit <- make_toy_fit()
  diag <- make_toy_diagnostics(fit)
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fit"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"))
  build <- function(...) automated_scoring_table(0:3, 0:3, 0:3, categories = 0:3, ...)
  x <- build(fit = fit, diagnostics = diag, rater_facet = "Rater", automated_rater = "R01")
  r <- which(fit$facets$others$Facet == "Rater" & fit$facets$others$Level == "R01")
  d <- which(diag$measures$Facet == "Rater" & diag$measures$Level == "R01")
  expect_equal(x$table$MFRMSeverity,
    unname(-fit$config$facet_signs["Rater"] * fit$facets$others$Estimate[r]))
  expect_equal(x$table$MFRMInfit, diag$measures$Infit[d])
  expect_equal(x$table$MFRMOutfit, diag$measures$Outfit[d])
  expect_equal(x$table$MFRMRatingRows, sum(fit$prep$data$Rater == "R01"))
  expect_identical(x$table$MFRMScope, "all_retained_source_ratings_for_selected_rater")
  expect_identical(x$mfrm$config, fit$config)
  expect_identical(x$table$MFRMRater, "R01")
  expect_equal(x$table$N, 4)
  y <- build(fit = fit, rater_facet = "Rater", automated_rater = "R01")
  expect_true(is.na(y$table$MFRMInfit))
  expect_match(y$mfrm$notes[["fit"]], "not supplied")
  expect_identical(y$criteria, x$criteria)
  bad <- diag; bad$measures$Estimate[d] <- bad$measures$Estimate[d] + .1
  expect_error(build(fit = fit, diagnostics = bad, rater_facet = "Rater",
    automated_rater = "R01"), "mismatch")
  expect_error(build(fit = fit, rater_facet = "Rater", automated_rater = "absent"), "exactly one")
  expect_error(build(diagnostics = diag), "require `fit`", fixed = TRUE)
})

test_that("saved tables replay and export without invoking model work", {
  x <- automated_scoring_table(c(0, 1, 1, 2), c(0, 0, 1, 2),
    c(0, 1, 2, 2), categories = 0:2)
  rds <- tempfile(fileext = ".rds"); csv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(rds, csv)))
  saveRDS(x, rds)
  expect_identical(readRDS(rds), x)
  utils::write.csv(x$criteria, csv, row.names = FALSE)
  expect_equal(utils::read.csv(csv, colClasses = vapply(x$criteria, class, "")), x$criteria)
})

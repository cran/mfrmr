ol_data <- function(binary = FALSE, crossing = FALSE) {
  cells <- expand.grid(Match = c(-1, 0, 1), Group = c("R", "F"))
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    eta <- cells$Match[i] + (cells$Group[i] == "F") *
      (if (crossing) 1.2 * cells$Match[i] else .4 + .6 * cells$Match[i])
    prob <- diff(c(0, plogis((if (binary) 0 else c(-.6, .9)) - eta), 1))
    y <- rep(seq_along(prob) - 1L, pmax(1, round(80 * prob)))
    data.frame(Person = paste(i, seq_along(y), sep = "_"), Item = "A",
      Match = cells$Match[i], Group = as.character(cells$Group[i]), Score = y)
  }))
}

ol_call <- function(d, categories = 0:2, groups = c("R", "F"), ...) {
  ordinal_logistic_table(d, "Person", "Item", "Score", "Group", "Match", groups, categories, ...)
}

test_that("ordinal likelihoods match an independent ordered-threshold optimizer", {
  skip_if_not_installed("ordinal")
  d <- ol_data(); a <- ol_call(d)
  matrices <- list(cbind(d$Match), cbind(d$Match, d$Group == "F"),
    cbind(d$Match, d$Group == "F", d$Match * (d$Group == "F")))
  ll <- vapply(matrices, function(X) {
    objective <- function(par) {
      cuts <- c(par[1], par[1] + exp(par[2]))
      eta <- drop(X %*% par[-(1:2)])
      cdf <- plogis(outer(eta, cuts, function(e, c) c - e))
      prob <- cbind(cdf[, 1], cdf[, 2] - cdf[, 1], 1 - cdf[, 2])
      -sum(log(prob[cbind(seq_len(nrow(d)), d$Score + 1L)]))
    }
    f <- optim(c(-.5, log(1.5), rep(0, ncol(X))), objective, method = "BFGS",
      control = list(reltol = 1e-12, maxit = 2000))
    expect_equal(f$convergence, 0L)
    -f$value
  }, 0)
  expect_true(all(a$table$Status == "available"))
  expect_equal(a$models$LogLik[1:3], ll, tolerance = 1e-9)
  expect_equal(a$table$ChiSquared, 2 * c(ll[2] - ll[1], ll[3] - ll[2], ll[3] - ll[1]), tolerance = 1e-9)
  expect_equal(a$table$DF, c(1, 1, 2))
  n <- nrow(d); f <- as.numeric(table(d$Score)); null <- sum(f * log(f / n))
  r2 <- (1 - exp(2 * (null - ll) / n)) / (1 - exp(2 * null / n))
  expect_equal(a$models$NagelkerkeR2[1:3], r2, tolerance = 1e-9)
  expect_equal(a$table$DeltaNagelkerkeR2, c(r2[2] - r2[1], r2[3] - r2[2], r2[3] - r2[1]), tolerance = 1e-9)
  expect_equal(a$table$PAdjusted, p.adjust(a$table$PValue, "holm"))
})

test_that("binary limiting case matches binomial GLM and probability direction", {
  skip_if_not_installed("ordinal")
  d <- ol_data(TRUE); d$g <- as.integer(d$Group == "F")
  a <- ol_call(d, 0:1)
  f <- lapply(list(Score ~ Match, Score ~ Match + g, Score ~ Match * g), function(f)
    glm(f, data = d, family = binomial()))
  expect_equal(a$models$LogLik, vapply(f, function(x) as.numeric(logLik(x)), 0), tolerance = 1e-10)
  p <- predict(f[[3]], newdata = expand.grid(Match = c(-1, 0, 1), g = 0:1), type = "response")
  expect_equal(a$contrasts$OddsRatio, unname((p[1:3] / (1 - p[1:3])) / (p[4:6] / (1 - p[4:6]))), tolerance = 1e-7)
  expect_identical(a$proportional_odds$Status, "not_applicable_binary")
  expect_true(is.na(a$proportional_odds$PValue))
})

test_that("global proportional odds diagnostic agrees with a saturated likelihood", {
  skip_if_not_installed("ordinal")
  # Two matching levels in each group: cut-specific intercepts/slopes saturate four cells.
  d <- subset(ol_data(), Match != 0)
  a <- ol_call(d)
  counts <- xtabs(~ Match + Group + Score, d)
  saturated <- sum(apply(counts, c(1, 2), function(v) sum(v * log(v / sum(v)))))
  expect_equal(a$models$LogLik[a$models$Model == "NonPO"], saturated, tolerance = 1e-10)
  expect_identical(a$proportional_odds$Status, "available")
  expect_equal(a$proportional_odds$DF, 3)
  expect_equal(a$proportional_odds$ChiSquared, 2 * (saturated - a$models$LogLik[3]), tolerance = 1e-9)
})

test_that("five categories preserve order, ignore spacing and use the correct diagnostic df", {
  skip_if_not_installed("ordinal")
  d <- expand.grid(Match = c(-1, 1), Group = c("R", "F"), Score = 0:4, Replicate = 1:10)
  d$Person <- seq_len(nrow(d)); d$Item <- "A"
  a <- ol_call(d, 0:4)
  expect_true(all(a$models$Status == "available"))
  expect_equal(a$models$Parameters, c(5, 6, 7, 16))
  expect_equal(a$proportional_odds$DF, 9)
  expect_equal(a$table$ChiSquared, rep(0, 3), tolerance = 1e-10)
  labels <- c(0, 1, 2, 10, 100)
  d$Score <- labels[d$Score + 1]
  expect_equal(ol_call(d, labels)$table, a$table)
})

test_that("group/category reversal and affine matching changes preserve comparisons", {
  skip_if_not_installed("ordinal")
  d <- ol_data(); a <- ol_call(d)
  reverse <- ol_call(d, groups = c("F", "R"))
  expect_equal(reverse$table$PValue, a$table$PValue, tolerance = 1e-8)
  expect_equal(reverse$contrasts$LogOddsRatio, -a$contrasts$LogOddsRatio, tolerance = 1e-7)
  reverse_score <- ol_call(d, categories = 2:0)
  expect_equal(reverse_score$table$PValue, a$table$PValue, tolerance = 1e-8)
  expect_equal(reverse_score$contrasts$LogOddsRatio, -a$contrasts$LogOddsRatio, tolerance = 1e-7)
  d$Match <- 1e100 + d$Match * 2e99
  scaled <- ol_call(d)
  expect_equal(scaled$table$PValue, a$table$PValue, tolerance = 1e-8)
  expect_equal(scaled$contrasts$LogOddsRatio, a$contrasts$LogOddsRatio, tolerance = 1e-7)
  negative <- ol_data(); negative$Match <- 10 - 2 * negative$Match
  negative <- ol_call(negative)
  expect_equal(negative$table$PValue, a$table$PValue, tolerance = 1e-8)
  expect_equal(negative$contrasts$LogOddsRatio, rev(a$contrasts$LogOddsRatio), tolerance = 1e-7)
  d$Score <- factor(c("low", "medium", "high")[d$Score + 1], levels = c("high", "low", "medium"))
  labels <- ol_call(d, c("low", "medium", "high"))
  expect_equal(labels$table, scaled$table)
  expect_equal(ol_call(dplyr::group_by(d, Group), c("low", "medium", "high"))$table, scaled$table)
})

test_that("missing rows, target failures and multiplicity families remain visible", {
  skip_if_not_installed("ordinal")
  d <- ol_data(); a <- ol_call(d)
  extra <- d[1:3, ]; extra$Person <- c("missing_y", "missing_x", "missing_both")
  extra$Score[c(1, 3)] <- NA; extra$Match[c(2, 3)] <- NA
  bad <- d; bad$Item <- "unused"; bad$Score <- pmin(1, bad$Score)
  x <- ol_call(rbind(d, extra, bad))
  expect_equal(x$table$PValue[1:3], a$table$PValue)
  expect_equal(x$table$PAdjusted, p.adjust(x$table$PValue, "holm", n = 6))
  expect_equal(x$settings$adjustment_family_size, 6L)
  expect_equal(x$settings$diagnostic_family_size, 2L)
  expect_identical(tail(x$matching$Status, 1), "unused_category")
  expect_identical(x$cases$Reason[nrow(d) + 1:3], c("missing_score", "missing_matching", "missing_score_and_matching"))
  expect_false(any(x$cases$Included[x$cases$Target == "unused"]))
  expect_true(all(is.na(x$table$PValue[4:6])))
  d$Score <- NA_integer_
  empty <- ol_call(d)
  expect_equal(nrow(empty$coefficients), 0L)
  expect_identical(empty$matching$Status, "no_complete_responses")
  expect_true(all(is.na(empty$table$PValue)))
})

test_that("overlap and within-group matching variation are required", {
  skip_if_not_installed("ordinal")
  d <- ol_data(); d$Match[d$Group == "F"] <- d$Match[d$Group == "F"] + 2
  expect_identical(ol_call(d)$matching$Status, "no_matching_range_overlap")
  d$Match[d$Group == "F"] <- 0
  expect_identical(ol_call(d)$matching$Status, "constant_matching_within_group")
  expect_identical(ol_call(subset(d, Group == "R"))$matching$Status, "missing_group")
})

test_that("separation is unavailable and numerical failure is not a null finding", {
  skip_if_not_installed("ordinal")
  d <- ol_data(); d$Score <- as.integer(d$Match > 0)
  x <- ol_call(d, 0:1)
  expect_false(any(x$models$Status == "available"))
  expect_true(all(is.na(x$table$PValue)))
  expect_true(all(is.na(x$contrasts$OddsRatio)))
  expect_true(all(nzchar(x$models$Message)))
  expect_true(all(x$cases$Included))
})

test_that("crossing group effects remain visible in conditional contrasts", {
  skip_if_not_installed("ordinal")
  x <- ol_call(ol_data(crossing = TRUE))
  expect_lt(x$table$PValue[x$table$Comparison == "Nonuniform"], .001)
  expect_gt(x$contrasts$LogOddsRatio[1], 0)
  expect_lt(x$contrasts$LogOddsRatio[3], 0)
})

test_that("input guards precede exclusions and avoid silent category recoding", {
  d <- ol_data()
  dup <- d[1, ]; dup$Score <- NA_integer_
  expect_error(ol_call(rbind(d, dup)), "Duplicate Person-target")
  changed <- d[1, ]; changed$Item <- "other"; changed$Group <- "F"; changed$Score <- NA_integer_
  expect_error(ol_call(rbind(d, changed)), "exactly one group")
  expect_error(ol_call(d, 0:1), "outside declared")
  expect_error(ol_call(d, c(0, 1, 1)), "distinct nonmissing")
  expect_error(ol_call(d, c("0", "1", "2")), "matching numeric or character")
  expect_error(ol_call(d, groups = c("F", "F")), "two distinct")
  d$Match[1] <- Inf
  expect_error(ol_call(d), "finite numeric")
  d$Match <- factor(d$Match)
  expect_error(ol_call(d), "finite numeric")
})

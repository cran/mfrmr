mh_rows <- function(tables, name = "Item") {
  do.call(rbind, lapply(seq_along(tables), function(k) {
    z <- tables[[k]]
    cells <- data.frame(Group = c("R", "F", "R", "F"), Score = c(1, 1, 0, 0))
    d <- cells[rep(seq_len(4), as.vector(z)), , drop = FALSE]
    d$Person <- paste(name, k, seq_len(nrow(d)), sep = "_")
    d$Target <- name; d$Match <- k
    rownames(d) <- NULL
    d
  }))
}

mh_call <- function(d, ...) mantel_haenszel_table(d, "Person", "Target", "Score",
  "Group", "Match", c("R", "F"), ...)

# Independently enumerate the null law of the sum of reference successes.
mh_null <- function(tables) {
  supports <- lapply(tables, function(z) seq(max(0, sum(z[1, ]) - sum(z[, 2])),
                                             min(sum(z[1, ]), sum(z[, 1]))))
  patterns <- as.matrix(expand.grid(supports))
  probs <- apply(patterns, 1, function(a) prod(vapply(seq_along(tables), function(k) {
    z <- tables[[k]]
    dhyper(a[k], sum(z[, 1]), sum(z[, 2]), sum(z[1, ]))
  }, 0)))
  support <- sort(unique(rowSums(patterns)))
  list(s = support, p = vapply(support, function(s) sum(probs[rowSums(patterns) == s]), 0))
}

test_that("exact inference matches independently enumerated conditional laws", {
  tables <- list(matrix(c(3, 1, 2, 4), 2), matrix(c(4, 2, 1, 3), 2))
  x <- mh_call(mh_rows(tables))
  law <- mh_null(tables)
  observed <- sum(vapply(tables, function(z) z[1, 1], 0))
  p <- sum(law$p[law$p <= law$p[law$s == observed] * (1 + 1e-7)])
  expect_equal(x$table$PValue, p, tolerance = 1e-13)
  expect_equal(x$table$MHOddsRatio, 12 / 2, tolerance = 1e-13)
  expect_identical(x$table$EstimateMethod, "conditional_mle")
  tilted <- function(or) { p <- law$p * or^law$s; p / sum(p) }
  expect_equal(sum(law$s * tilted(x$table$OddsRatio)), observed, tolerance = 1e-4)
  # R's uniroot default uses ~1e-4 coordinate precision; check absolute tail error.
  expect_lt(abs(sum(tilted(x$table$Lower)[law$s >= observed]) - .025), 1e-4)
  expect_lt(abs(sum(tilted(x$table$Upper)[law$s <= observed]) - .025), 1e-4)
  expect_true(all(x$cases$Included))
  expect_true(is.na(x$table$ChiSquared))
  expect_false(x$settings$continuity_correction)
})

test_that("asymptotic statistics and log-odds intervals agree with formulas", {
  tables <- list(matrix(c(3, 1, 2, 4), 2), matrix(c(4, 2, 1, 3), 2))
  d <- mh_rows(tables)
  x <- mh_call(d, inference = "asymptotic", correct = FALSE)$table
  abcd <- do.call(rbind, lapply(tables, as.vector))
  a <- abcd[, 1]; c <- abcd[, 2]; b <- abcd[, 3]; e <- abcd[, 4]
  n <- a + b + c + e
  delta <- sum(a - (a + b) * (a + c) / n)
  variance <- sum((a + b) * (c + e) * (a + c) * (b + e) / (n^2 * (n - 1)))
  or <- sum(a * e / n) / sum(b * c / n)
  log_var <- sum((a * e + or * b * c) * (a + e + or * (b + c)) / n^2) /
    (2 * sum(a * e / n)^2)
  expect_equal(x$ChiSquared, delta^2 / variance, tolerance = 1e-13)
  expect_equal(x$PValue, pchisq(delta^2 / variance, 1, lower.tail = FALSE))
  expect_equal(c(x$Lower, x$Upper), exp(log(or) + c(-1, 1) * qnorm(.975) * sqrt(log_var)))
  corrected <- mh_call(d, inference = "asymptotic", correct = TRUE)$table
  expect_equal(corrected$ChiSquared, (abs(delta) - .5)^2 / variance)
  tiny <- mh_rows(list(matrix(c(1, 1, 1, 2), 2)))
  expect_equal(mh_call(tiny, inference = "asymptotic", correct = TRUE)$table$ChiSquared,
               mh_call(tiny, inference = "asymptotic", correct = FALSE)$table$ChiSquared)
})

test_that("one stratum, reversed orientation and boundary estimates are supported", {
  tables <- list(matrix(c(3, 1, 2, 4), 2))
  d <- mh_rows(tables)
  x <- mh_call(d)
  reference <- fisher.test(tables[[1]])
  expect_equal(x$table$PValue, reference$p.value, tolerance = 1e-12)
  expect_equal(x$table$OddsRatio, unname(reference$estimate), tolerance = 1e-10)
  reverse <- mantel_haenszel_table(d, "Person", "Target", "Score", "Group", "Match", c("F", "R"))
  expect_equal(reverse$table$MHOddsRatio, 1 / x$table$MHOddsRatio)
  expect_equal(reverse$table$PValue, x$table$PValue)
  expect_equal(c(reverse$table$Lower, reverse$table$Upper), 1 / c(x$table$Upper, x$table$Lower), tolerance = 1e-4)
  edge <- mh_rows(list(matrix(c(0, 5, 5, 0), 2)))
  ex <- mh_call(edge)$table
  expect_equal(ex$MHOddsRatio, 0); expect_equal(ex$OddsRatio, 0)
  expect_equal(ex$Lower, 0); expect_true(is.finite(ex$Upper))
  asy <- mh_call(edge, inference = "asymptotic")$table
  expect_true(is.finite(asy$PValue))
  expect_true(is.na(asy$Lower) && is.na(asy$Upper))
  expect_identical(asy$IntervalStatus, "unavailable_boundary_or_numerical")
  edge$Score <- 1 - edge$Score
  ex <- mh_call(edge)$table
  expect_equal(ex$OddsRatio, Inf); expect_equal(ex$Upper, Inf)
})

test_that("missingness and uninformative strata are audited without padding", {
  d <- mh_rows(list(matrix(c(3, 1, 2, 4), 2)))
  extras <- data.frame(Group = c("R", "F", "R", "R", "F", "R"),
    Score = c(0, 0, 1, NA, 1, NA), Person = paste0("X", 1:6), Target = "Item",
    Match = c(2, 2, 3, 4, NA, NA))
  x <- mh_call(rbind(d, extras))
  expect_equal(x$table$PValue, mh_call(d)$table$PValue)
  expect_equal(x$table$ExcludedRows, 6)
  expect_identical(tail(x$cases$Reason, 6), c("constant_response", "constant_response",
    "no_group_overlap", "missing_score", "missing_matching", "missing_score_and_matching"))
  expect_identical(x$strata$Reason, c("informative", "constant_response", "no_group_overlap", "no_complete_responses"))
  extras$Target <- "Unavailable"
  x <- mh_call(rbind(d, extras))
  expect_equal(x$table$PAdjusted[1], min(1, 2 * x$table$PValue[1]))
  expect_true(is.na(x$table$PValue[2]))
  expect_equal(x$settings$adjustment_family_size, 2)
  extras$Match <- NA_real_
  empty <- mh_call(extras)
  expect_equal(nrow(empty$strata), 0)
  expect_identical(empty$table$Status, "no_informative_strata")
})

test_that("joint strata retain plan support and arbitrary column names", {
  d <- mh_rows(list(matrix(c(3, 1, 2, 4), 2), matrix(c(4, 2, 1, 3), 2)))
  d$Plan <- ifelse(d$Match == 1, "A:B", "A")
  d$Match <- ifelse(d$Match == 1, "C", "B:C")
  x <- mantel_haenszel_table(d, "Person", "Target", "Score", "Group", c("Plan", "Match"), c("R", "F"))
  expect_equal(nrow(x$matching), 2)
  expect_equal(x$table$InformativeStrata, 2)
  shuffled <- d[rev(seq_len(nrow(d))), ]
  y <- mantel_haenszel_table(shuffled, "Person", "Target", "Score", "Group", c("Plan", "Match"), c("R", "F"))
  expect_equal(x$table, y$table)
  names(d) <- paste0("Column ", seq_along(d))
  z <- mantel_haenszel_table(d, "Column 3", "Column 4", "Column 2", "Column 1",
    c("Column 6", "Column 5"), c("R", "F"))
  expect_equal(z$table, x$table)
  d2 <- rbind(shuffled, transform(shuffled, Target = "Other"))
  expect_equal(nrow(mantel_haenszel_table(d2, "Person", "Target", "Score", "Group",
    c("Plan", "Match"), c("R", "F"))$table), 2)
  grouped <- dplyr::group_by(shuffled, Group)
  expect_equal(mantel_haenszel_table(grouped, "Person", "Target", "Score", "Group",
    c("Plan", "Match"), c("R", "F"))$table, x$table)
})

test_that("input validation refuses hidden dependence and recoding", {
  d <- mh_rows(list(matrix(c(3, 1, 2, 4), 2)))
  expect_error(mh_call(rbind(d, d[1, ])), "Duplicate Person-target")
  other <- d; other$Target <- "Other"; other$Group <- ifelse(other$Group == "R", "F", "R")
  expect_error(mh_call(rbind(d, other)), "exactly one group")
  bad <- d; bad$Score[1] <- 2
  expect_error(mh_call(bad), "numeric 0, 1")
  bad$Score <- factor(bad$Score)
  expect_error(mh_call(bad), "numeric 0, 1")
  bad <- d; bad$Group[1] <- "Other"
  expect_error(mh_call(bad), "outside")
  bad <- d; bad$Person[1] <- NA
  expect_error(mh_call(bad), "Identifiers")
  bad <- d; bad$Match[1] <- Inf
  expect_error(mh_call(bad), "Identifiers")
  expect_error(mh_call(d, conf_level = 1), "strictly")
  expect_error(mh_call(d, correct = NA), "TRUE or FALSE")
  expect_error(mantel_haenszel_table(d, "Person", "Target", "Score", "Group", "Score", c("R", "F")), "distinct")
})

test_that("opposing stratum effects remain visible when MH cancels", {
  d <- mh_rows(list(matrix(c(10, 40, 40, 10), 2), matrix(c(40, 10, 10, 40), 2)))
  x <- mh_call(d)
  expect_equal(x$table$MHOddsRatio, 1)
  expect_equal(x$table$PValue, 1)
  expect_equal(x$strata$OddsRatio, c(1 / 16, 16))
})

test_that("large integer strata retain finite MH estimates and intervals", {
  for (n in c(2000L, 50000L)) {
    d <- mh_rows(list(matrix(rep(n, 4), 2)))
    expect_no_warning(x <- mh_call(d, inference = "asymptotic"))
    expect_equal(x$strata$OddsRatio, 1)
    expect_equal(x$table$MHOddsRatio, 1)
    expect_equal(x$table$OddsRatio, 1)
    expect_equal(x$table$PValue, 1)
    # A single stratum reduces to the four-cell log-odds variance.
    expect_equal(c(x$table$Lower, x$table$Upper),
      exp(c(-1, 1) * qnorm(.975) * sqrt(4 / n)), tolerance = 1e-12)
    expect_identical(x$table$Status, "available")
    expect_identical(x$table$IntervalStatus, "available")
    expect_identical(x$table$Message, "")
  }
})

test_that("double-precision MH agrees with the standard numeric-table backend", {
  # The array interface preserves double storage, unlike the factor interface.
  tab <- array(as.double(c(50000, 49980, 50020, 50000,
    2000, 1800, 1900, 2100, 0, 60, 50, 10)), c(2, 2, 3))
  counts <- data.frame(ReferenceSuccess = tab[1, 1, ],
    ReferenceFailure = tab[1, 2, ], FocalSuccess = tab[2, 1, ],
    FocalFailure = tab[2, 2, ])
  for (correct in c(FALSE, TRUE)) {
    a <- mfrmr:::.mfrm_mh_asymptotic(counts, correct, .99)
    b <- stats::mantelhaen.test(tab, correct = correct, conf.level = .99)
    for (name in c("estimate", "statistic", "parameter", "p.value", "conf.int")) {
      expect_equal(unname(a[[name]]), as.numeric(b[[name]]), tolerance = 1e-12)
    }
  }
})

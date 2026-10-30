icc_pair_data <- function(a, b) {
  data.frame(Unit = rep(seq_along(a), 2), Judge = rep(c("A", "B"), each = length(a)),
    Score = c(a, b))
}

test_that("two-way pair ICCs agree with an independent ANOVA decomposition", {
  # ANOVA through a fitted linear model, independently of the API's matrix sums.
  d <- icc_pair_data(c(2, 4, 5, 3, 8, 7, 10), c(3, 6, 4, 5, 9, 7, 8))
  a <- stats::anova(stats::lm(Score ~ factor(Unit) + factor(Judge), data = d))
  r <- a[1, "Mean Sq"]; c <- a[2, "Mean Sq"]; e <- a[3, "Mean Sq"]; n <- 7
  denominators <- c(agreement_single = r + e + 2 * (c - e) / n,
    consistency_single = r + e, agreement_average = r + (c - e) / n,
    consistency_average = r)
  for (definition in c("agreement", "consistency")) for (type in c("single", "average")) {
    x <- pairwise_icc_table(d, "Judge", "Unit", "Score", definition = definition, type = type)
    expect_equal(x$table$ICC, unname((r - e) / denominators[paste(definition, type, sep = "_")]),
      tolerance = 1e-12)
    z <- x$anova
    expect_equal(unname(unlist(z[c("ScaledMSUnit", "ScaledMSRater", "ScaledMSError")])) *
      z$ScoreScale^2, c(r, c, e), tolerance = 1e-12)
    expect_equal(unname(unlist(z[c("DFUnit", "DFRater", "DFError")])), c(6, 1, 6))
  }
  mixed <- pairwise_icc_table(d, "Judge", "Unit", "Score")
  random <- pairwise_icc_table(d, "Judge", "Unit", "Score", model = "random")
  expect_identical(mixed$table$ICC, random$table$ICC)
  expect_identical(mixed$table$Model, "two_way_mixed")
  expect_identical(random$table$Model, "two_way_random")
  expect_identical(random$settings$model, "random")
  expect_match(random$notes[4], "sampling design")
})

test_that("absolute agreement retains offsets and average refers to exactly two ratings", {
  d <- icc_pair_data(1:6, 2:7)
  absolute <- pairwise_icc_table(d, "Judge", "Unit", "Score")
  consistent <- pairwise_icc_table(d, "Judge", "Unit", "Score", definition = "consistency")
  average <- pairwise_icc_table(d, "Judge", "Unit", "Score", type = "average")
  expect_equal(absolute$table$ICC, 7 / 8)
  expect_equal(consistent$table$ICC, 1)
  expect_equal(average$table$ICC, 14 / 15)
  expect_identical(c(absolute$table$Form, consistent$table$Form, average$table$Form),
    c("ICC(A,1)", "ICC(C,1)", "ICC(A,2)"))
  expect_identical(average$settings$raters_per_pair, 2L)
  expect_match(average$notes[5], "mean of these two ratings")
  identical_scores <- pairwise_icc_table(icc_pair_data(1:6, 1:6), "Judge", "Unit", "Score")
  expect_equal(identical_scores$table$ICC, 1)
})

test_that("negative estimates and degenerate denominators are not made reassuring", {
  d <- icc_pair_data(1:6, 6:1)
  x <- pairwise_icc_table(d, "Judge", "Unit", "Score", definition = "consistency")
  expect_equal(x$table$ICC, -1)
  expect_identical(x$table$Status, "available")
  x <- pairwise_icc_table(d, "Judge", "Unit", "Score", type = "average")
  expect_true(is.na(x$table$ICC))
  expect_true(is.finite(x$table$RawICC) && x$table$RawICC > 1)
  expect_identical(x$table$Status, "nonpositive_denominator")
  constant <- icc_pair_data(rep(2, 6), rep(2, 6))
  x <- pairwise_icc_table(constant, "Judge", "Unit", "Score")
  expect_true(is.na(x$table$ICC) && is.na(x$table$RawICC))
  expect_identical(x$table$Status, "nonpositive_denominator")
  different <- icc_pair_data(rep(2, 6), rep(3, 6))
  x <- pairwise_icc_table(different, "Judge", "Unit", "Score")
  expect_equal(x$table$ICC, 0)
  x <- pairwise_icc_table(different, "Judge", "Unit", "Score", definition = "consistency")
  expect_true(is.na(x$table$ICC))
})

test_that("pair-specific matching preserves missing scores and absent assignments", {
  d <- data.frame(Judge = c(rep("A", 4), rep("B", 4), "C", "D"),
    Unit = c(1:4, 1:3, 5, 1, 9), Score = c(1, 2, NA, 4, 1, 3, 3, 5, 1, NA))
  x <- pairwise_icc_table(d, "Judge", "Unit", "Score")
  expect_equal(nrow(x$table), 6)
  expect_equal(x$table$N, c(2, 1, 0, 1, 0, 0))
  expect_equal(x$table$CandidateUnits, c(5, 4, 5, 4, 5, 2))
  expect_equal(x$table$ExcludedUnits, x$table$CandidateUnits - x$table$N)
  expect_identical(x$table$Status[-1], rep("fewer_than_two_complete_units", 5))
  ab <- x$cases[x$cases$Pair == 1, ]
  expect_equal(x$units$Unit[ab$Unit], 1:5)
  expect_identical(ab$Missing1, c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(ab$Assigned1, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(ab$Assigned2, c(TRUE, TRUE, TRUE, FALSE, TRUE))
  expect_identical(ab$Included, c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(ab$Score1, d$Score[ab$Row1])
  expect_equal(ab$Score2, d$Score[ab$Row2])
  alone <- pairwise_icc_table(d[d$Judge %in% c("A", "B"), ], "Judge", "Unit", "Score")
  expect_equal(x$table$ICC[1], alone$table$ICC)
})

test_that("composite unit keys and row permutations never silently merge ratings", {
  d <- data.frame(Judge = rep(c("A", "B"), each = 3),
    Person = rep(c("a|b", "a", "z"), 2), Task = rep(c("c", "b|c", "x"), 2),
    Score = c(1, 4, 5, 2, 3, 5))
  x <- pairwise_icc_table(d, "Judge", c("Person", "Task"), "Score")
  expect_equal(x$table$N, 3)
  expect_equal(nrow(x$units), 3)
  expect_identical(names(x$units), c("Person", "Task"))
  expect_equal(x$units[x$cases$Unit, ], d[x$cases$Row1, c("Person", "Task")],
    ignore_attr = TRUE)
  reversed <- pairwise_icc_table(d[6:1, ], "Judge", c("Person", "Task"), "Score")
  expect_equal(x$table$ICC, reversed$table$ICC, tolerance = 1e-12)
  expect_identical(reversed$table$Rater1, "B")
  d$Judge <- factor(d$Judge, levels = c("unused", "B", "A"))
  expect_equal(pairwise_icc_table(d, "Judge", c("Person", "Task"), "Score")$table, x$table)
})

test_that("ambiguous cells, invalid scores and invalid identities fail explicitly", {
  d <- icc_pair_data(1:6, 2:7)
  bad <- rbind(d, d[1, ]); bad$Score[nrow(bad)] <- NA
  expect_error(pairwise_icc_table(bad, "Judge", "Unit", "Score"), "Duplicate")
  for (value in list(factor(d$Score), as.character(d$Score), rep(Inf, nrow(d)))) {
    bad <- d; bad$Score <- value
    expect_error(pairwise_icc_table(bad, "Judge", "Unit", "Score"), "finite numeric")
  }
  for (value in c(NA_character_, " ")) {
    bad <- d; bad$Judge[1] <- value
    expect_error(pairwise_icc_table(bad, "Judge", "Unit", "Score"), "Identifiers must")
  }
  bad <- d; bad$Unit[1] <- Inf
  expect_error(pairwise_icc_table(bad, "Judge", "Unit", "Score"), "Identifiers must")
  # Distinct doubles can print as the same text label: never coalesce those raters.
  bad <- d; bad$Judge <- rep(c(1, 1 + .Machine$double.eps), each = 6)
  expect_error(pairwise_icc_table(bad, "Judge", "Unit", "Score"), "identical text labels")
  expect_error(pairwise_icc_table(d, "Judge", "Score", "Score"), "distinct existing")
  expect_error(pairwise_icc_table(d, "Judge", character(), "Score"), "must name")
  expect_error(pairwise_icc_table(d[1:6, ], "Judge", "Unit", "Score"), "two observed raters")
  expect_error(pairwise_icc_table(d[FALSE, ], "Judge", "Unit", "Score"), "nonempty")
})

test_that("stable rescaling avoids overflow and underflow without changing the target", {
  d <- icc_pair_data(c(-3, -1, 2, 3, 4), c(-2, 0, 1, 4, 3))
  baseline <- pairwise_icc_table(d, "Judge", "Unit", "Score")
  for (multiplier in c(1e200, 1e-200, -1e200)) {
    z <- d; z$Score <- z$Score * multiplier
    x <- pairwise_icc_table(z, "Judge", "Unit", "Score")
    expect_equal(x$table$ICC, baseline$table$ICC, tolerance = 1e-12)
    expect_identical(x$table$Status, "available")
  }
  z <- d; z$Score <- z$Score + 1e12
  expect_equal(pairwise_icc_table(z, "Judge", "Unit", "Score")$table$ICC,
    baseline$table$ICC, tolerance = 1e-12)
  z <- icc_pair_data(c(-1e308, 0, 1e308), c(-1e308, 0, 1e308))
  expect_equal(pairwise_icc_table(z, "Judge", "Unit", "Score")$table$ICC, 1)
})

test_that("observed-score results need no fitted model and save as ordinary tables", {
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fitting"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"))
  d <- icc_pair_data(1:6, 2:7)
  x <- pairwise_icc_table(d, "Judge", "Unit", "Score")
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(x, path)
  expect_identical(readRDS(path), x)
  expect_identical(x$settings$uncertainty, "not_computed")
  expect_true(all(vapply(x[c("table", "anova", "cases", "units")], is.data.frame, logical(1))))
})

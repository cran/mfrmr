alpha_paper_data <- function() {
  d <- expand.grid(Unit = 1:12, Judge = LETTERS[1:4])
  # Krippendorff (2011/2013), physical pp. 4 and 9; 41 scores, 40 pairable.
  d$Score <- c(1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
               1, 2, 3, 3, 2, 2, 4, 1, 2, 5, NA, 3,
               NA, 3, 3, 3, 2, 3, 4, 2, 2, 5, 1, NA,
               1, 2, 3, 3, 2, 4, 4, 1, 2, 5, 1, NA)
  d
}

alpha_pair_reference <- function(d, distance) {
  # Enumerate ordered observation pairs, independently of the API's moments/counts.
  groups <- split(d$Score, d$Unit)
  groups <- lapply(groups, function(z) z[!is.na(z)])
  groups <- groups[lengths(groups) >= 2]
  numerator <- 0
  for (z in groups) {
    for (i in seq_along(z)) for (j in seq_along(z)) if (i != j)
      numerator <- numerator + distance(z[i], z[j]) / (length(z) - 1)
  }
  pooled <- unlist(groups, use.names = FALSE); n <- length(pooled)
  expected <- 0
  for (i in seq_along(pooled)) for (j in seq_along(pooled)) if (i != j)
    expected <- expected + distance(pooled[i], pooled[j])
  c(Alpha = 1 - numerator * (n - 1) / expected,
    ObservedDisagreement = numerator / n, ExpectedDisagreement = expected / (n * (n - 1)))
}

test_that("all four alpha levels reproduce the primary source's incomplete-panel example", {
  d <- alpha_paper_data()
  count <- c(9, 13, 10, 5, 3)
  ordinal_position <- (cumsum(count) - count / 2) / 40
  distances <- list(nominal = function(a, b) as.numeric(a != b),
    ordinal = function(a, b) (ordinal_position[a] - ordinal_position[b])^2,
    interval = function(a, b) ((a - b) / 4)^2,
    ratio = function(a, b) ((a - b) / (a + b))^2)
  printed <- c(nominal = .743, ordinal = .815, interval = .849, ratio = .797)
  for (level in names(distances)) {
    x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = level,
      categories = if (level == "ordinal") 1:5 else NULL)
    ref <- alpha_pair_reference(d, distances[[level]])
    expect_equal(unname(unlist(x$table[names(ref)])), unname(ref), tolerance = 1e-12)
    expect_lt(abs(x$table$Alpha - printed[level]), .0005)
    z <- x$counts[match(1:5, x$counts$Category), ]
    expect_equal(z$PairableRatings, count)
    expect_equal(x$table$PairableRatings, 40)
    expect_equal(x$table$NonmissingRatings, 41)
    expect_equal(x$table$MissingScores, 7)
    expect_equal(x$table$SingletonUnits, 1)
    expect_identical(x$table$Status, "available")
  }
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score")
  expect_equal(x$table$Alpha, 113 / 152, tolerance = 1e-12)
})

test_that("binary and nominal source examples use finite-sample pooled chance disagreement", {
  d <- data.frame(Unit = rep(1:10, 2), Judge = rep(c("Meg", "Owen"), each = 10),
    Score = c(0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0))
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score")
  expect_equal(x$table$Alpha, 2 / 21)
  expect_equal(x$table$ObservedDisagreement, 8 / 20)
  expect_equal(x$table$ExpectedDisagreement, 2 * 14 * 6 / (20 * 19))
  d <- data.frame(Unit = rep(1:12, 2), Judge = rep(c("Ben", "Gerry"), each = 12),
    Score = c("a", "a", "b", "b", "d", "c", "c", "c", "e", "d", "d", "a",
              "b", "a", "b", "b", "b", "c", "c", "c", "e", "d", "d", "d"))
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score")
  expect_equal(x$table$Alpha, 155 / 224)
  expect_equal(x$table$ObservedDisagreement, .25)
})

test_that("unpaired scores never enter margins and assignment audit distinguishes missingness", {
  d <- alpha_paper_data()
  base <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval")
  d$Score[d$Unit == 12 & d$Judge == "B"] <- 1e308
  changed <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval")
  expect_identical(changed$table, base$table)
  expect_identical(changed$settings$score_half_range, base$settings$score_half_range)
  expect_equal(changed$counts$PairableRatings[changed$counts$Category == 1e308], 0)
  expect_identical(unique(changed$cases$Reason[changed$cases$Unit == 12 & !changed$cases$MissingScore]),
    "fewer_than_two_ratings")
  d <- d[!is.na(d$Score), ]
  sparse <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval")
  expect_equal(sparse$table$Alpha, base$table$Alpha)
  expect_equal(sum(sparse$unit_summary$AbsentAssignments), 7)
  expect_equal(sum(base$unit_summary$MissingScores), 7)
  expect_equal(sum(base$unit_summary$AbsentAssignments), 0)
  expect_equal(sparse$cases$Score, d$Score[sparse$cases$Row])
  expect_equal(sparse$units$Unit[sparse$cases$Unit], d$Unit)
  d <- rbind(d, data.frame(Unit = 13, Judge = c("A", "B"), Score = NA_real_))
  z <- krippendorff_alpha_table(d, "Judge", "Unit", "Score")
  expect_equal(z$table$UnratedUnits, 1)
  expect_false(z$unit_summary$Included[z$units$Unit == 13])
})

test_that("ordinal alpha uses declared order and frequencies rather than equal rank spacing", {
  d <- alpha_paper_data()
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "ordinal", categories = 1:5)
  unused <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "ordinal",
    categories = c(-1, 1:4, 4.5, 5, 6))
  expect_equal(x$table, unused$table)
  expect_equal(x$counts$OrdinalPosition, (cumsum(c(9, 13, 10, 5, 3)) - c(9, 13, 10, 5, 3)/2)/40)
  expect_equal(unused$counts$PairableRatings[c(1, 6, 8)], c(0, 0, 0))
  reverse <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "ordinal", categories = 5:1)
  expect_equal(x$table$Alpha, reverse$table$Alpha)
  interval <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval")
  expect_gt(abs(x$table$Alpha - interval$table$Alpha), .02)
  labels <- c("Low", "Developing", "Middle", "High", "Advanced")
  d$Score <- factor(labels[d$Score], levels = rev(labels))
  labeled <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "ordinal", categories = labels)
  expect_equal(x$table, labeled$table)
  expect_identical(labeled$counts$Category, labels)
  expect_error(krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "ordinal"), "substantive ordinal order")
})

test_that("negative, perfect, single-unit and undefined alpha estimates remain distinct", {
  d <- data.frame(Unit = c(1, 2, 1, 2), Judge = c("A", "A", "B", "B"), Score = c(0, 1, 1, 0))
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score")
  expect_equal(x$table$Alpha, -.5)
  expect_identical(x$table$Status, "available")
  x <- krippendorff_alpha_table(d[d$Unit == 1, ], "Judge", "Unit", "Score")
  expect_equal(x$table$Alpha, 0)
  d$Score <- c(0, 1, 0, 1)
  expect_equal(krippendorff_alpha_table(d, "Judge", "Unit", "Score")$table$Alpha, 1)
  d$Score <- rep(1, 4)
  for (level in c("nominal", "ordinal", "interval", "ratio")) {
    x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = level,
      categories = if (level == "ordinal") 1:3 else NULL)
    expect_true(is.na(x$table$Alpha))
    expect_identical(x$table$Status, "zero_expected_disagreement")
  }
  d <- d[c(1, 4), ]
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score")
  expect_identical(x$table$Status, "no_pairable_ratings")
  expect_equal(x$table$PairableRatings, 0)
  d$Score <- NA_real_
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score")
  expect_equal(nrow(x$counts), 0)
  expect_true(is.na(x$table$Alpha) && is.na(x$table$ExpectedDisagreement))
})

test_that("numeric distances remain stable at extreme scales and preserve ratio zero", {
  d <- alpha_paper_data()
  for (level in c("interval", "ratio")) {
    base <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = level)
    for (scale in c(1e-200, 1e200)) {
      z <- d; z$Score <- z$Score * scale
      x <- krippendorff_alpha_table(z, "Judge", "Unit", "Score", level = level)
      expect_equal(x$table$Alpha, base$table$Alpha, tolerance = 1e-12)
      expect_identical(x$table$Status, "available")
    }
  }
  base <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval")
  d$Score <- d$Score + 1e12
  expect_equal(krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval")$table$Alpha,
    base$table$Alpha, tolerance = 1e-12)
  d <- data.frame(Unit = rep(1:3, 2), Judge = rep(c("A", "B"), each = 3),
    Score = c(0, 1, 2, 0, 2, 2))
  ref <- alpha_pair_reference(d, function(a, b) if (a + b == 0) 0 else ((a-b)/(a+b))^2)
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "ratio")
  expect_equal(x$table$Alpha, unname(ref["Alpha"]), tolerance = 1e-12)
  d$Score <- d$Score * 5e307
  expect_equal(krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "ratio")$table$Alpha,
    x$table$Alpha, tolerance = 1e-12)
  d$Score <- rep(c(-1e308, 0, 1e308), 2)
  expect_equal(krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval")$table$Alpha, 1)
  d$Score <- rep(c(1, 1 + .Machine$double.eps, 1), 2)
  expect_equal(krippendorff_alpha_table(d, "Judge", "Unit", "Score")$counts$PairableRatings, c(4, 2))
})

test_that("matching rejects repetitions and invalid scales before calculating alpha", {
  d <- alpha_paper_data()
  bad <- rbind(d, d[1, ]); bad$Score[nrow(bad)] <- NA
  expect_error(krippendorff_alpha_table(bad, "Judge", "Unit", "Score"), "Duplicate")
  for (value in c(Inf, -Inf)) {
    bad <- d; bad$Score[1] <- value
    expect_error(krippendorff_alpha_table(bad, "Judge", "Unit", "Score"), "finite numbers")
  }
  bad <- d; bad$Score <- as.character(bad$Score); bad$Score[1] <- " "
  expect_error(krippendorff_alpha_table(bad, "Judge", "Unit", "Score"), "nonblank")
  bad <- d; bad$Score <- factor(bad$Score)
  expect_error(krippendorff_alpha_table(bad, "Judge", "Unit", "Score", level = "interval"), "numeric")
  bad <- d; bad$Score[1] <- -1
  expect_error(krippendorff_alpha_table(bad, "Judge", "Unit", "Score", level = "ratio"), "nonnegative")
  expect_error(krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval", categories = 1:5), "Leave")
  expect_error(krippendorff_alpha_table(d, "Judge", "Unit", "Score", categories = 1:4), "outside")
  expect_error(krippendorff_alpha_table(d, "Judge", "Unit", "Score", categories = c(1:5, 5)), "distinct")
  expect_error(krippendorff_alpha_table(d, "Judge", "Unit", "Score", categories = as.character(1:5)), "representation")
  expect_error(krippendorff_alpha_table(d, "Judge", "Unit", "Score", categories = c(1:5, NA)), "nonmissing")
  bad <- d; bad$Unit[1] <- NA
  expect_error(krippendorff_alpha_table(bad, "Judge", "Unit", "Score"), "Identifiers")
})

test_that("composite units preserve identity and saved results need no fitted model", {
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fitting"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"))
  d <- data.frame(Judge = rep(c("A", "B"), each = 3),
    Person = rep(c("a|b", "a", "z"), 2), Task = rep(c("c", "b|c", "x"), 2),
    Score = c(1, 4, 5, 2, 3, 5))
  x <- krippendorff_alpha_table(d, "Judge", c("Person", "Task"), "Score", level = "interval")
  expect_equal(x$table$PairableUnits, 3)
  expect_equal(x$units[x$cases$Unit, ], d[x$cases$Row, c("Person", "Task")], ignore_attr = TRUE)
  reversed <- krippendorff_alpha_table(d[6:1, ], "Judge", c("Person", "Task"), "Score", level = "interval")
  expect_equal(x$table, reversed$table, tolerance = 1e-12)
  expect_identical(x$settings$uncertainty, "not_computed")
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(x, path)
  expect_identical(x, readRDS(path))
})

test_that("many distinct values require no category-square matrix", {
  # Also protects against a category-length count allocation for every tiny unit.
  d <- data.frame(Unit = rep(seq_len(5000), 2), Judge = rep(c("A", "B"), each = 5000),
    Score = seq_len(10000))
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score")
  expect_equal(x$table$Alpha, 0)
  expect_equal(x$table$PairableRatings, 10000)
  expect_equal(nrow(x$counts), 10000)
  x <- krippendorff_alpha_table(d, "Judge", "Unit", "Score", level = "interval")
  expect_true(is.finite(x$table$Alpha))
})

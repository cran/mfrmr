test_that("fractional nesting weights follow probabilities and are scale invariant", {
  data <- data.frame(A = c("a", "a", "b", "b"), B = c("x", "y", "x", "y"),
                     w = c(.9, .1, .1, .9))
  expected <- round(1 - (-.9 * log(.9) - .1 * log(.1)) / log(2), 4)
  for (scale in c(1, 10, 1e-200, 1e300)) {
    data$w <- c(.9, .1, .1, .9) * scale
    result <- detect_facet_nesting(data, c("A", "B"), weight_col = "w")$pairwise_table
    expect_equal(result$NestingIndex_AinB, expected)
    expect_equal(result$NestingIndex_BinA, expected)
    expect_equal(result$UsedRows, 4)
    expect_equal(result$ExcludedRows, 0)
  }
})

test_that("nesting indices agree with an independent joint mutual-information calculation", {
  counts <- matrix(c(7, 2, 1, 0, 3, 5, 4, 1, 2), nrow = 3)
  data <- expand.grid(A = letters[1:3], B = LETTERS[1:3])
  data$w <- as.vector(counts)
  p <- counts / sum(counts)
  product <- outer(rowSums(p), colSums(p))
  mutual_information <- sum(p[p > 0] * log(p[p > 0] / product[p > 0]))
  ha <- -sum(rowSums(p) * log(rowSums(p)))
  hb <- -sum(colSums(p) * log(colSums(p)))
  result <- detect_facet_nesting(data, c("A", "B"), weight_col = "w")$pairwise_table
  expect_equal(result$NestingIndex_AinB, round(mutual_information / hb, 4))
  expect_equal(result$NestingIndex_BinA, round(mutual_information / ha, 4))
  expect_equal(result$UsedRows, 8)
  expect_equal(result$ExcludedRows, 1)
  expanded <- data[rep(seq_len(nrow(data)), data$w), c("A", "B")]
  unweighted <- detect_facet_nesting(expanded, c("A", "B"))$pairwise_table
  expect_equal(result[c("NestingIndex_AinB", "NestingIndex_BinA")],
               unweighted[c("NestingIndex_AinB", "NestingIndex_BinA")])
})

test_that("marginal and conditional nesting entropies use the same observed pair", {
  data <- data.frame(A = c("a", "a", "b", "b"), B = c("x", "y", "x", "y"))
  incomplete <- rbind(data, data.frame(A = rep(NA_character_, 100), B = "x"))
  result <- detect_facet_nesting(incomplete, c("A", "B"))$pairwise_table
  expect_equal(result$NestingIndex_AinB, 0)
  expect_equal(result$NestingIndex_BinA, 0)
  expect_equal(result$UsedRows, 4)
  expect_equal(result$ExcludedRows, 100)
  expect_equal(result[c("LevelsA", "LevelsB")], data.frame(LevelsA = 2L, LevelsB = 2L))
  incomplete$B <- NA_character_
  absent <- detect_facet_nesting(incomplete, c("A", "B"))$pairwise_table
  expect_true(is.na(absent$NestingIndex_AinB))
  expect_true(is.na(absent$NestingIndex_BinA))
  expect_identical(absent$Direction, "insufficient_data")
  expect_equal(absent$UsedRows, 0)
  data$B <- "x"
  constant <- detect_facet_nesting(data, c("A", "B"))$pairwise_table
  expect_true(is.na(constant$NestingIndex_AinB))
  expect_equal(constant$NestingIndex_BinA, 0)
})

test_that("weight input is explicit and repeated person roles do not create self pairs", {
  data <- data.frame(Person = c("p", "p", "q", "q"),
    A = c("a", "a", "b", "b"), B = c("x", "y", "x", "y"), w = 1)
  expect_error(detect_facet_nesting(data, c("A", "B"), weight_col = "missing"), "must name a column")
  for (weights in list(c(1, 1, 1, -1), c(1, 1, 1, NA), c(1, 1, 1, Inf),
                       rep(0, 4), rep("1", 4), factor(rep(1, 4)))) {
    data$w <- weights
    expect_error(detect_facet_nesting(data, c("A", "B"), weight_col = "w"),
                 "finite, nonnegative numeric")
  }
  result <- detect_facet_nesting(data, c("Person", "A", "B"), person = "Person")
  expect_equal(result$summary$NPairs, 3)
  expect_false(any(result$pairwise_table$FacetA == result$pairwise_table$FacetB))
})

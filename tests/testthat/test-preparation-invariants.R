.preparation_input <- function() {
  d <- expand.grid(Person = paste0("p", 1:6), Rater = c("r1", "r2"),
                   Item = c("i1", "i2"), stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  d$Score <- rep(0:3, length.out = nrow(d))
  d
}

.prepare_input <- function(data, ...) prepare_mfrm_data(data, "Person",
  c("Rater", "Item"), "Score", keep_original = TRUE, ...)

test_that("near-integer scores are rounded to the intended category on both sides", {
  for (origin in c(-2L, 0L, 1L)) {
    raw <- .preparation_input()
    raw$Score <- raw$Score + origin
    expected <- .prepare_input(raw)
    for (coerce in list(identity, as.character, factor)) {
      input <- raw
      input$Score <- coerce(input$Score + rep(c(-1e-7, 1e-7), length.out = nrow(input)))
      actual <- .prepare_input(input, rating_min = origin, rating_max = origin + 3L)
      expect_identical(actual$data, expected$data)
      expect_identical(actual$score_map, expected$score_map)
      expect_identical(build_indices(actual)$score_k, as.integer(raw$Score - origin))
      note <- subset(actual$preparation_notes, Condition == "score_roundoff_normalized")
      expect_equal(note$Count, nrow(input))
    }
  }
})

test_that("prepared category indices select the correct likelihood contributions", {
  data <- .preparation_input()
  expected_score <- data$Score
  data$Score <- data$Score - 1e-7
  data$weight <- rep(c(.3, pi), length.out = nrow(data))
  prep <- .prepare_input(data, weight_col = "weight")
  idx <- build_indices(prep, step_facet = "Item", slope_facet = "Rater")
  eta <- seq(-1, 1, length.out = nrow(data))
  steps <- rbind(c(0, -.4, -.2, .5), c(0, -.1, .1, .8))
  slopes <- c(.6, 1.5)
  reference <- function(model) sum(vapply(seq_len(nrow(data)), function(i) {
    step <- if (model == "RSM") steps[1, ] else steps[match(data$Item[i], c("i1","i2")), ]
    a <- if (model == "GPCM") slopes[match(data$Rater[i], c("r1","r2"))] else 1
    kernel <- exp(a * ((0:3) * eta[i] - step))
    data$weight[i] * log(kernel[expected_score[i] + 1L] / sum(kernel))
  }, numeric(1)))
  expect_equal(loglik_rsm(eta, idx$score_k, steps[1, ], idx$weight), reference("RSM"), tolerance = 1e-12)
  expect_equal(loglik_pcm(eta, idx$score_k, steps, idx$step_idx, idx$weight), reference("PCM"), tolerance = 1e-12)
  expect_equal(loglik_gpcm(eta, idx$score_k, steps, idx$step_idx, slopes,
    idx$slope_idx, idx$weight), reference("GPCM"), tolerance = 1e-12)
  expect_identical(idx$weight, data$weight)
})

test_that("dplyr grouping cannot add unselected columns or alter row retention", {
  data <- tibble::as_tibble(.preparation_input())
  data$Unused <- c(NA_real_, rep(1, nrow(data) - 1L))
  baseline <- .prepare_input(data)
  for (input in list(dplyr::group_by(data, Unused), dplyr::group_by(data, Item),
                    dplyr::rowwise(data))) {
    actual <- .prepare_input(input)
    expect_identical(actual$data, baseline$data)
    expect_identical(actual$row_retention, baseline$row_retention)
    expect_identical(actual$omitted_input_rows, integer(0))
    expect_false("Unused" %in% names(actual$data))
  }
  grouped <- dplyr::group_by(data, Unused)
  plain_review <- describe_mfrm_data(data, "Person", c("Rater", "Item"), "Score", include_agreement = FALSE)
  grouped_review <- describe_mfrm_data(grouped, "Person", c("Rater", "Item"), "Score", include_agreement = FALSE)
  expect_identical(grouped_review, plain_review)
})

test_that("internal names cannot overwrite facets and numeric weights retain precision", {
  data <- .preparation_input()
  for (name in c("Person", "Score", "Weight", "score_k")) {
    changed <- data
    names(changed)[names(changed) == "Person"] <- "Examinee"
    names(changed)[names(changed) == "Score"] <- "Rating"
    names(changed)[names(changed) == "Item"] <- name
    expect_error(prepare_mfrm_data(changed, "Examinee", c("Rater", name), "Rating"),
                 "reserved internal names")
    expect_error(describe_mfrm_data(changed, "Examinee", c("Rater", name), "Rating"),
                 "reserved internal names")
  }
  data$w <- rep(pi, nrow(data))
  expect_identical(.prepare_input(data, weight_col = "w")$data$Weight, data$w)
  data$w <- factor(rep(c("0.25", "1.75"), length.out = nrow(data)), levels = c("1.75", "0.25"))
  expect_identical(.prepare_input(data, weight_col = "w")$data$Weight, as.numeric(as.character(data$w)))
})

test_that("category and weight overflow fail before coercion or range allocation", {
  data <- .preparation_input()
  for (value in c(.Machine$integer.max + 1, -.Machine$integer.max - 1)) {
    changed <- data; changed$Score[1] <- value
    expect_error(.prepare_input(changed), "category codes exceed R's integer range")
    expect_error(.prepare_input(data, rating_min = value), "rating_min.*single finite integer")
    expect_error(.prepare_input(data, rating_max = value), "rating_max.*single finite integer")
  }
  expect_error(.prepare_input(data, rating_min = -.Machine$integer.max,
    rating_max = .Machine$integer.max), "integer indexing capacity")
  expect_error(.prepare_input(data, rating_min = 0+1i), "rating_min.*single finite integer")
  data$Score <- as.complex(data$Score)
  expect_error(.prepare_input(data), "vector of real ordered integer")
  data <- .preparation_input()
  data$w <- rep(.Machine$double.xmax, nrow(data))
  expect_error(.prepare_input(data, weight_col = "w"), "sum of retained weights exceeds")
  data$w <- 1+1i
  expect_error(.prepare_input(data, weight_col = "w"), "vector of real numeric")
  data$Score[1] <- 1.001
  expect_error(.prepare_input(data), "Fractional value")
})

test_that("public estimators receive the same ratings after normalization", {
  skip_if_not_installed("lpSolve")
  patterns <- as.matrix(expand.grid(rep(list(0:1), 4L)))
  patterns <- patterns[rep(seq_len(nrow(patterns)), 1L + patterns[, 1L]), , drop = FALSE]
  data <- data.frame(Person = rep(seq_len(nrow(patterns)), each = 4L),
    Rater = rep(c("A", "A", "B", "B"), nrow(patterns)),
    Item = rep(c("X", "Y", "X", "Y"), nrow(patterns)),
    Score = as.vector(t(patterns)), Unused = NA_real_)
  data <- tibble::as_tibble(data)
  changed <- data
  changed$Score <- changed$Score + rep(c(-1e-7, 1e-7), length.out = nrow(data))
  changed <- dplyr::group_by(changed, Unused)
  for (method in c("CML", "MML", "JML")) {
    fit <- function(d) fit_mfrm(d, "Person", c("Rater", "Item"), "Score",
      method = method, category_policy = "preserve", rating_min = 0, rating_max = 1)
    expected <- fit(data)
    actual <- fit(changed)
    expect_identical(actual$prep$data, expected$prep$data)
    expect_gt(length(actual$opt$par), 0L)
    expect_true(all(is.finite(c(actual$opt$par, actual$opt$value))))
    expect_equal(actual$opt$par, expected$opt$par, tolerance = 1e-12)
    expect_equal(actual$opt$value, expected$opt$value, tolerance = 1e-12)
    likelihood <- function(x) if (method == "CML") as.numeric(logLik(x)) else x$summary$LogLik
    expect_equal(likelihood(actual), likelihood(expected), tolerance = 1e-12)
    expect_true("score_roundoff_normalized" %in% actual$prep$preparation_notes$Condition)
  }
})

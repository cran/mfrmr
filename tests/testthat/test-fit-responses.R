response_extractor_cases <- function() {
  out <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits
  out$GPCM <- readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds"))
  out
}

test_that("plug-in responses agree with independently assembled category probabilities", {
  for (f in response_extractor_cases()) {
    data <- f$prep$data; cfg <- f$config
    theta <- f$facets$person$Estimate[match(data$Person, f$facets$person$Person)]
    eta <- theta
    for (facet in cfg$facet_names) {
      tab <- f$facets$others[f$facets$others$Facet == facet, ]
      eta <- eta + cfg$facet_signs[[facet]] * tab$Estimate[match(data[[facet]], tab$Level)]
    }
    independent <- t(vapply(seq_len(nrow(data)), function(i) {
      steps <- f$steps$Estimate[f$steps$StepFacet == as.character(data[[cfg$step_facet]][i])]
      slope <- if (cfg$model == "GPCM")
        f$slopes$Estimate[match(data[[cfg$slope_facet]][i], f$slopes$SlopeFacet)] else 1
      w <- slope * (seq.int(0, length(steps)) * eta[i] - c(0, cumsum(steps)))
      p <- exp(w - max(w)); p / sum(p)
    }, numeric(f$config$n_cat)))
    probabilities <- fitted(f, "prob")
    expect_equal(matrix(as.numeric(probabilities), nrow(data)), independent, tolerance = 1e-12)
    scores <- f$prep$rating_min + seq_len(ncol(independent)) - 1L
    mean <- drop(independent %*% scores)
    variance <- rowSums(independent * (matrix(scores, nrow(data), length(scores), byrow = TRUE) - mean)^2)
    expect_equal(as.numeric(fitted(f)), mean, tolerance = 1e-12)
    expect_equal(as.numeric(residuals(f)), data$Score - mean, tolerance = 1e-12)
    expect_equal(as.numeric(residuals(f, "pearson")), (data$Score - mean) / sqrt(variance), tolerance = 1e-12)
    expect_equal(attr(probabilities, "diagnostics")$Variance, variance, tolerance = 1e-12)
    expect_identical(predict(f), fitted(f))
    expect_identical(predict(f, type = "prob"), probabilities)
    expect_identical(attr(probabilities, "identifiers"), as.data.frame(data[c("Person", cfg$facet_names)]))
    expect_identical(names(fitted(f)), as.character(seq_len(nrow(data))))
    expect_true(all(attr(probabilities, "diagnostics")$ProbabilityAvailable))
  }
})

test_that("primary Person identities and unavailable estimates are respected", {
  f <- response_extractor_cases()[[1]]
  expected <- fitted(f)
  reordered <- f; reordered$facets$person <- f$facets$person[nrow(f$facets$person):1, ]
  expect_identical(fitted(reordered), expected)
  missing <- f; missing$facets$person$Estimate[1] <- NA_real_
  x <- fitted(missing); rows <- f$prep$data$Person == f$facets$person$Person[1]
  expect_true(all(is.na(x[rows])))
  expect_equal(as.numeric(x[!rows]), as.numeric(expected[!rows]))
  expect_true(all(!attr(x, "diagnostics")$ProbabilityAvailable[rows]))
  expect_match(attr(x, "diagnostics")$Reason[which(rows)[1]], "no optimizer or display proxy")
  duplicate <- f; duplicate$facets$person$Person[2] <- duplicate$facets$person$Person[1]
  expect_error(fitted(duplicate), "distinct, nonmissing identifiers")
  # Conditional per-rating probabilities do not multiply by rating weights.
  weighted <- f; weighted$prep$data$Weight <- seq_len(nrow(weighted$prep$data)) / 3
  expect_equal(as.numeric(residuals(weighted)), as.numeric(residuals(f)))
  expect_identical(attr(expected, "source_readiness"), mfrmr_get_readiness_record(f)$fit)
})

test_that("extreme JML uses its true boundary limit and leaves Pearson undefined", {
  pattern <- rbind(c(0, 0), c(1, 1), c(0, 1), c(1, 0))
  x <- pattern[rep(seq_len(4L), c(2, 2, 4, 6)), ]
  data <- data.frame(Person = rep(seq_len(nrow(x)), each = 2),
    Rater = rep(c("A", "B"), nrow(x)), Score = as.vector(t(x)))
  f <- fit_mfrm(data, "Person", "Rater", "Score", method = "JML",
    maxit = 300, category_policy = "preserve")
  expect_true(any(is.infinite(f$facets$person$Estimate)))
  expect_true(all(is.finite(f$opt$par)))
  theta <- f$facets$person$Estimate[match(f$prep$data$Person, f$facets$person$Person)]
  low <- theta == -Inf; high <- theta == Inf
  p <- fitted(f, "prob"); expected <- fitted(f); raw <- residuals(f); pearson <- residuals(f, "pearson")
  severity <- f$facets$others$Estimate[match(f$prep$data$Rater, f$facets$others$Level)]
  expect_equal(as.numeric(p[, 2]), unname(plogis(theta - severity)), tolerance = 1e-12)
  expect_true(all(p[low, 1] == 1) && all(p[low, 2] == 0))
  expect_true(all(p[high, 1] == 0) && all(p[high, 2] == 1))
  expect_true(all(expected[low] == 0) && all(expected[high] == 1))
  expect_true(all(raw[low | high] == 0))
  expect_true(all(is.na(pearson[low | high])))
  expect_true(all(attr(pearson, "diagnostics")$Status[low | high] == "zero_variance"))
  expect_true(all(grepl("infinite JML Person", attr(pearson, "diagnostics")$Reason[low | high])))
  # Tiny positive variance must not be rounded up or replaced with NA.
  near <- f; near$facets$person$Estimate[] <- -30
  near_p <- fitted(near, "prob")
  v <- attr(near_p, "diagnostics")$Variance
  expect_true(all(v > 0 & v < 1e-10))
  expect_true(all(is.finite(residuals(near, "pearson"))))
})

test_that("omitted input rows, repeated events and collapsed scale remain explicit", {
  # Retain repeated events and two original category labels with a gap.
  data <- data.frame(Person = rep(1:4, each = 4), Rater = rep(c("A", "B"), 8),
    Score = rep(c(2, 4, 4, 2), 4))
  data$Score[c(3, 14)] <- NA_real_
  f <- suppressWarnings(fit_mfrm(data, "Person", "Rater", "Score", method = "JML",
    repeated_ratings = "independent", category_policy = "collapse", maxit = 200))
  out <- fitted(f)
  expect_length(out, 14L)
  expect_equal(attr(out, "rows")$InputRow, setdiff(1:16, c(3, 14)))
  expect_identical(attr(out, "score_map"), f$prep$score_map)
  expect_equal(colnames(fitted(f, "prob")), as.character(sort(unique(f$prep$data$Score))))
  expect_equal(as.numeric(residuals(f)) + as.numeric(out), f$prep$data$Score)
  expect_false(identical(f$prep$score_map$OriginalScore, f$prep$score_map$InternalScore))
  legacy <- f; legacy$prep$omitted_input_rows <- NULL
  expect_true(all(is.na(attr(fitted(legacy), "rows")$InputRow)))
})

test_that("unsupported targets and extra arguments do not silently fall through", {
  f <- response_extractor_cases()[[1]]
  expect_error(fitted(f, typo = TRUE), "must be empty")
  expect_error(residuals(f, "deviance"), "arg")
  expect_error(predict(f, newdata = f$prep$data), "retained source ratings only")
  expect_error(predict(f, se.fit = TRUE), "must be empty")
  for (class in c("mfrm_cml_fit", "mfrm_comparison_fit")) {
    x <- structure(list(), class = class)
    expect_error(fitted(x), "do not supply primary Person estimates")
    expect_error(residuals(x), "do not supply primary Person estimates")
    expect_error(predict(x), "do not supply primary Person estimates")
  }
  imported <- f; class(imported) <- c("mfrm_imported_fit", class(f))
  expect_error(fitted(imported), "native fit_mfrm")
  adjusted <- f; adjusted$config$jml_adjustment <- TRUE
  expect_error(fitted(adjusted), "mfrm_response_diagnostics")
  penalized <- f; penalized$config$facet_shrinkage <- "ridge"
  expect_error(residuals(penalized), "unpenalized")
  nonfinite <- f; nonfinite$facets$others$Estimate[1] <- Inf
  expect_error(fitted(nonfinite), "Finite structural calibration")
  zero <- response_extractor_cases()$GPCM; zero$slopes$Estimate[1] <- 0
  expect_error(fitted(zero), "strictly positive primary slopes")
})

test_that("saved extraction does not invoke fitting, scoring or covariance", {
  f <- response_extractor_cases()[[1]]
  before <- serialize(f, NULL)
  expected <- list(fitted(f), fitted(f, "prob"), residuals(f), residuals(f, "pearson"))
  blocked <- function(...) stop("Unexpected numerical estimation.")
  testthat::local_mocked_bindings(fit_mfrm = blocked, predict_mfrm_units = blocked,
    score_mfrm_persons = blocked, compute_mml_parameter_covariance = blocked,
    diagnose_mfrm = blocked, .package = "mfrmr")
  saved <- unserialize(serialize(f, NULL))
  expect_identical(list(predict(saved), predict(saved, type = "prob"),
    residuals(saved), residuals(saved, "pearson")), expected)
  expect_identical(serialize(f, NULL), before)
})

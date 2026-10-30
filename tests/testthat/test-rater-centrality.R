# Reuse the same native fixture across contracts; these are arithmetic and
# availability checks, not a simulation of centrality detection performance.
centrality_pcm_fixture <- function(binary = FALSE) {
  key <- paste0("centrality_pcm_", binary)
  if (!exists(key, envir = .mfrmr_test_cache, inherits = FALSE)) {
    d <- load_mfrmr_data("example_core")
    if (binary) d$Score <- as.integer(d$Score > 2)
    x <- suppressMessages(suppressWarnings(fit_mfrm(d, "Person",
      c("Rater", "Criterion"), "Score", method = "JML", model = "PCM",
      step_facet = "Rater", maxit = 100)))
    assign(key, x, envir = .mfrmr_test_cache)
  }
  get(key, envir = .mfrmr_test_cache, inherits = FALSE)
}

test_that("RSM indices agree with a separate adjacent-logit calculation", {
  f <- make_toy_fit()
  d <- f$prep$data
  # Eq. 1: P(k) proportional to exp(k * eta - sum(tau[1:k])).
  # Build this from the public coefficient tables, without response helpers.
  eta <- f$facets$person$Estimate[match(d$Person, f$facets$person$Person)]
  for (facet in f$config$facet_names) {
    co <- subset(f$facets$others, Facet == facet)
    eta <- eta + unname(f$config$facet_signs[facet]) *
      co$Estimate[match(d[[facet]], co$Level)]
  }
  k <- seq.int(0, f$config$n_cat - 1L)
  logw <- outer(eta, k) - matrix(c(0, cumsum(f$steps$Estimate)),
    nrow(d), length(k), byrow = TRUE)
  p <- exp(logw - apply(logw, 1, max))
  p <- p / rowSums(p)
  e <- as.vector(p %*% k) + f$config$rating_min
  v <- as.vector(p %*% (k^2)) - (e - f$config$rating_min)^2
  residual <- d$Score - e
  expected <- lapply(levels(d$Rater), function(r) {
    at <- which(d$Rater == r)
    er <- e[at] - mean(e[at]); rr <- residual[at] - mean(residual[at])
    c(sum(er * rr) / sqrt(sum(er^2) * sum(rr^2)),
      sum(residual[at]^2) / sum(v[at]), mean(residual[at]^2 / v[at]))
  })
  x <- rater_centrality_table(f)
  expect_equal(as.matrix(x$table[c("ResidualExpectedCorrelation", "Infit", "Outfit")]),
    do.call(rbind, expected), ignore_attr = TRUE, tolerance = 1e-11)
  expect_true(all(x$table$ResidualReason == "available"))
  expect_true(all(x$table$FitReason == "available"))
  expect_true(all(is.na(x$table$ThresholdSD)))
  expect_true(all(grepl("not supplied", x$table$ThresholdReason)))
  expect_equal(sum(x$categories$Ratings), nrow(d))
  expect_equal(as.numeric(tapply(x$categories$Proportion, x$categories$Rater, sum)),
    rep(1, nlevels(d$Rater)))
  expect_identical(x$settings$residual_basis, "RSM_point_score_plugin_raw_residuals")
  # The short toy optimizer is deliberately provisional, even with finite indices.
  expect_true(all(x$table$ResidualSourceStatus == "source_needs_review"))
  expect_identical(x$sources$rsm$config, f$config)
})

test_that("PCM spreads use all steps with a declared sample divisor", {
  f <- centrality_pcm_fixture()
  x <- rater_centrality_table(threshold_fit = f)
  by_rater <- split(f$steps$Estimate, f$steps$StepFacet)
  # Pairwise squared distances provide an independent variance identity.
  expected <- vapply(by_rater, function(s) {
    sum(outer(s, s, "-")^2) / (2 * length(s) * (length(s) - 1))
  }, numeric(1))
  expect_equal(x$table$ThresholdVariance, unname(expected))
  expect_equal(x$table$ThresholdSD, sqrt(unname(expected)))
  expect_equal(x$table$ThresholdCount, rep(3L, 4))
  expect_true(all(x$table$ThresholdsOrdered))
  expect_true(all(x$thresholds$SDDivisor == 2))
  expect_equal(x$thresholds$Estimate, f$steps$Estimate)
  expect_true(all(is.na(x$table$ResidualExpectedCorrelation)))
  expect_true(all(is.na(x$table$Infit)))
  expect_true(all(grepl("not supplied", x$table$ResidualReason)))
  expect_identical(x$settings$sd_divisor, "count_minus_one")
  expect_equal(x$table$ObservedSD, unname(vapply(split(f$prep$data$Score,
    f$prep$data$Rater), sd, numeric(1))))
})

test_that("two sources require the same roster including repeated observations", {
  rsm <- make_toy_fit(); pcm <- centrality_pcm_fixture()
  x <- rater_centrality_table(rsm, pcm)
  expect_true(x$settings$matched_roster)
  expect_identical(x$table$Infit, rater_centrality_table(rsm)$table$Infit)
  expect_identical(x$table$ThresholdSD,
    rater_centrality_table(threshold_fit = pcm)$table$ThresholdSD)
  permuted <- pcm
  permuted$prep$data <- permuted$prep$data[rev(seq_len(nrow(pcm$prep$data))), ]
  expect_identical(rater_centrality_table(rsm, permuted)$table, x$table)
  bad <- pcm; bad$prep$data$Score[1] <- bad$prep$data$Score[1] + 1
  expect_error(rater_centrality_table(rsm, bad), "same retained response roster")
  bad <- pcm; bad$prep$data <- rbind(bad$prep$data, bad$prep$data[1, ])
  expect_error(rater_centrality_table(rsm, bad), "same retained response roster")
  bad <- pcm; bad$config$score_map$OriginalScore[1] <- 99
  expect_error(rater_centrality_table(rsm, bad), "category mapping")
  bad <- pcm; bad$config$facet_signs["Rater"] <- 1
  expect_error(rater_centrality_table(rsm, bad), "facet constraints")
  expect_error(rater_centrality_table(rsm, make_toy_fit(model = "PCM")),
    "steps owned by")
})

test_that("unsupported sources and altered step tables cannot masquerade as indices", {
  rsm <- make_toy_fit(); pcm <- centrality_pcm_fixture()
  expect_error(rater_centrality_table(), "Supply an RSM")
  expect_error(rater_centrality_table(rsm, rater_facet = NA_character_), "facet name")
  expect_error(rater_centrality_table(rsm, rater_facet = "Person"), "non-Person")
  expect_error(rater_centrality_table(pcm), "native additive RSM")
  expect_error(rater_centrality_table(threshold_fit = rsm), "native additive PCM")
  bad <- rsm; bad$prep$data$Weight[1] <- 2
  expect_error(rater_centrality_table(bad), "unit weights")
  bad <- rsm; class(bad) <- c("mfrm_imported_fit", class(bad))
  expect_error(rater_centrality_table(bad), "native additive")
  bad <- rsm; bad$config$interaction_specs <- list(example = 1)
  expect_error(rater_centrality_table(bad), "without shrinkage or interactions")
  bad <- rsm; bad$config$facet_shrinkage <- "ridge"
  expect_error(rater_centrality_table(bad), "without shrinkage or interactions")
  bad <- pcm; bad$steps$Estimate[1] <- bad$steps$Estimate[1] + .1
  expect_error(rater_centrality_table(threshold_fit = bad), "match the fitted")
  bad <- pcm; bad$steps <- bad$steps[-1, ]
  expect_error(rater_centrality_table(threshold_fit = bad), "match the fitted")
})

test_that("binary and disordered steps retain their distinct meanings", {
  binary <- rater_centrality_table(threshold_fit = centrality_pcm_fixture(TRUE))
  expect_true(all(binary$table$ThresholdCount == 1L))
  expect_true(all(is.na(binary$table$ThresholdSD)))
  expect_true(all(grepl("binary", binary$table$ThresholdReason)))
  f <- centrality_pcm_fixture()
  params <- expand_params(f$opt$par, build_param_sizes(f$config), f$config)
  # A coherent parameter vector with reversed category steps, not a sorted
  # display-table edit. Spreads must stay unchanged and disorder stay visible.
  params$steps_mat <- params$steps_mat[, 3:1]
  f$opt$par <- collapse_expanded_params(params, f$config)
  f$steps <- build_step_table(f$config, f$prep, params)
  x <- rater_centrality_table(threshold_fit = f)
  expect_true(all(!x$table$ThresholdsOrdered))
  expect_equal(x$table$ThresholdSD,
    rater_centrality_table(threshold_fit = centrality_pcm_fixture())$table$ThresholdSD)
  expect_identical(x$thresholds$Estimate, f$steps$Estimate)
})

test_that("boundary and constant rows are not silently dropped or labeled zero", {
  f <- make_toy_fit()
  # All raters scored this Person. An infinite public JML estimate must not
  # leak a finite internal optimizer proxy into the descriptive calculation.
  f$facets$person$Estimate[1] <- Inf
  x <- rater_centrality_table(f)
  expect_true(all(x$table$ResidualUnavailableRows > 0))
  expect_true(all(is.na(x$table$ResidualExpectedCorrelation)))
  expect_true(all(is.na(x$table$Infit)))
  expect_true(all(grepl("nonfinite primary", x$table$ResidualReason)))
  expect_equal(x$table$ResidualRows + x$table$ResidualUnavailableRows, x$table$RatingRows)
  f <- make_toy_fit()
  # With all free locations and steps zero, every expected score is identical.
  f$opt$par[] <- 0
  f$facets$person$Estimate[] <- 0; f$facets$others$Estimate[] <- 0
  f$steps$Estimate[] <- 0
  x <- rater_centrality_table(f)
  expect_true(all(is.na(x$table$ResidualExpectedCorrelation)))
  expect_true(all(grepl("constant", x$table$ResidualReason)))
  expect_true(all(is.finite(x$table$Infit)))
  expect_true(all(x$table$FitReason == "available"))
})

test_that("MML uses saved EAP values and all outputs replay without model work", {
  f <- make_toy_fit(method = "MML", maxit = 500)
  pcm <- centrality_pcm_fixture()
  # Diagnostics are only a comparator prepared before forbidding their use.
  obs <- compute_obs_table(f)
  local_mocked_bindings(fit_mfrm = function(...) stop("Unexpected fit"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"))
  x <- rater_centrality_table(f, pcm)
  expected <- vapply(split(seq_len(nrow(obs)), obs$Rater), function(i)
    cor(obs$Residual[i], obs$Expected[i]), numeric(1))
  expect_equal(x$table$ResidualExpectedCorrelation, unname(expected))
  expect_true(all(x$table$ResidualMethod == "MML"))
  expect_true(all(x$table$ThresholdMethod == "JML"))
  expect_true(any(grepl("saved Person EAP", x$notes)))
  f$facets$person$Estimate <- f$facets$person$Estimate + .5
  changed <- rater_centrality_table(f)
  expect_false(isTRUE(all.equal(changed$table$Infit, x$table$Infit)))
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(x, file)
  expect_identical(readRDS(file), x)
  for (name in c("table", "thresholds", "categories")) {
    expect_s3_class(x[[name]], "data.frame")
    expect_false(any(vapply(x[[name]], is.list, logical(1))))
  }
})

marginal_pair_fixture <- function() {
  list(prep = list(data = data.frame(Person = c("P1", "P1", "P2", "P2"), Rater = c("A", "B", "A", "B"))),
       config = list(facet_names = "Rater"),
       logprob_bundle = list(prob_list = list(matrix(0.5, 4, 2))),
       posterior_bundle = list(obs_posterior = matrix(1, 4, 1)),
       observed_cat = c(0, 1, 0, 1), weights = rep(1, 4), categories = 0:1)
}

test_that("pairwise integration agrees with literal weighted category-pair probabilities", {
  d <- expand.grid(Person = c("P1", "P2"), Rater = c("A", "B", "C"),
    Criterion = c("X", "Y"), stringsAsFactors = FALSE)
  d <- d[c(1:11, 1), ] # Unequal contexts and a duplicate level within a context.
  n <- nrow(d)
  probabilities <- lapply(1:3, function(q) {
    p <- matrix(seq_len(n * 4) %% (5 + q) + 1, n, 4)
    p / rowSums(p)
  })
  posterior <- t(vapply(d$Person, function(p)
    if (p == "P1") c(.2, .3, .5) else c(.6, .3, .1), numeric(3)))
  core <- list(prep = list(data = d), config = list(facet_names = c("Rater", "Criterion")),
    logprob_bundle = list(prob_list = probabilities),
    posterior_bundle = list(obs_posterior = posterior),
    observed_cat = rep(0:3, length.out = n), weights = rep(c(0, .5, 1, 2), length.out = n),
    categories = 0:3)
  reference <- list()
  for (facet in core$config$facet_names) {
    context <- setdiff(names(d), facet)
    for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n)) {
      if (d[[facet]][i] == d[[facet]][j] ||
          !all(vapply(context, function(name) d[[name]][i] == d[[name]][j], TRUE))) next
      w <- core$weights[i] * core$weights[j]
      if (w == 0) next
      exact <- adjacent <- 0
      # Literal joint probability, independent of the matrix implementation.
      for (q in 1:3) for (a in 0:3) for (b in 0:3) {
        mass <- posterior[i, q] * probabilities[[q]][i, a + 1L] * probabilities[[q]][j, b + 1L]
        if (a == b) exact <- exact + mass
        if (abs(a - b) <= 1) adjacent <- adjacent + mass
      }
      levels <- sort(c(d[[facet]][i], d[[facet]][j]))
      reference[[length(reference) + 1L]] <- data.frame(Facet = facet,
        Level1 = levels[1], Level2 = levels[2], LevelPairCount = 1,
        OpportunityWeight = w,
        ObservedExactCount = w * (core$observed_cat[i] == core$observed_cat[j]),
        ExpectedExactCount = w * exact, ExactVarianceCount = w^2 * exact * (1 - exact),
        ObservedAdjacentCount = w * (abs(core$observed_cat[i] - core$observed_cat[j]) <= 1),
        ExpectedAdjacentCount = w * adjacent, AdjacentVarianceCount = w^2 * adjacent * (1 - adjacent))
    }
  }
  expected <- do.call(rbind, reference)
  keys <- c("Facet", "Level1", "Level2")
  expected <- aggregate(expected[setdiff(names(expected), keys)], expected[keys], sum)
  actual <- calc_marginal_pairwise_bundle(core)$pair_stats
  key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
  expect_equal(as.data.frame(actual[match(key(expected), key(actual)), names(expected)]),
    expected, tolerance = 1e-13, ignore_attr = TRUE)
  expect_true(all(actual$UnavailableContextPairs == 0))

  # Node and observation order must not alter the diagnostic target.
  rows <- rev(seq_len(n)); nodes <- c(3L, 1L, 2L)
  core$prep$data <- d[rows, ]
  core$logprob_bundle$prob_list <- lapply(probabilities[nodes], function(p) p[rows, ])
  core$posterior_bundle$obs_posterior <- posterior[rows, nodes]
  core$observed_cat <- core$observed_cat[rows]
  core$weights <- core$weights[rows]
  expect_equal(calc_marginal_pairwise_bundle(core)$pair_stats, actual, tolerance = 1e-13)
})

test_that("invalid node probabilities and inconsistent posteriors remain unavailable", {
  core <- marginal_pair_fixture()
  core$logprob_bundle$prob_list <- rep(core$logprob_bundle$prob_list, 2)
  core$posterior_bundle$obs_posterior <- matrix(rep(c(1, 0), each = 4), 4, 2)
  for (bad in list(NA_real_, Inf, -.1, 1.1, .6)) {
    x <- core
    # An invalid zero-posterior node is still invalid, not silently omitted.
    x$logprob_bundle$prob_list[[2]][3, 1] <- bad
    result <- calc_marginal_pairwise_bundle(x)$pair_stats
    expect_equal(result$UnavailableContextPairs, 1)
    expect_true(is.na(result$ExpectedExactCount) && is.na(result$ExpectedAdjacentCount))
    expect_true(is.na(result$Flagged))
  }
  core$posterior_bundle$obs_posterior[3, ] <- c(.3, .7)
  expect_equal(calc_marginal_pairwise_bundle(core)$pair_stats$UnavailableContextPairs, 1)
})

test_that("approximate posterior agreement is not shared transitively across a context", {
  delta <- .3 * sqrt(.Machine$double.eps)
  posterior <- rbind(c(.5 + delta, .5 - delta), c(.5, .5), c(.5 + 2 * delta, .5 - 2 * delta))
  core <- marginal_pair_fixture()
  core$prep$data <- data.frame(Person = "P1", Rater = c("B", "A", "C"))
  core$posterior_bundle$obs_posterior <- posterior
  core$logprob_bundle$prob_list <- rep(list(matrix(.5, 3, 2)), 2)
  core$observed_cat <- c(0, 1, 1); core$weights <- rep(1, 3)
  x <- calc_marginal_pairwise_bundle(core)$pair_stats
  expect_identical(x$UnavailableContextPairs, c(0L, 1L, 0L))
  expect_identical(paste(x$Level1, x$Level2), c("A B", "A C", "B C"))
  expect_true(is.na(x$ExpectedExactCount[2]))
  expect_equal(x$ExpectedExactCount[c(1, 3)], c(.5, .5))
})

test_that("precomputed level order retains literal radix pairs and missing levels", {
  core <- marginal_pair_fixture()
  levels <- c("z", "\u00e9", "A", NA_character_)
  core$prep$data <- data.frame(Person = "P1", Rater = factor(levels))
  x <- calc_marginal_pairwise_bundle(core)$pair_stats
  pairs <- utils::combn(seq_along(levels), 2)
  expected <- apply(pairs, 2, function(pair) {
    sorted <- sort(levels[pair], method = "radix")
    paste(sorted[1], sorted[2], sep = "/")
  })
  expect_setequal(paste(x$Level1, x$Level2, sep = "/"), expected)
  expect_equal(nrow(x), 6)
  expect_equal(x$ExpectedExactCount, rep(.5, 6))
  expect_equal(x$UnavailableContextPairs, rep(0, 6))
})

test_that("marginal cells do not omit missing probabilities or fabricate classifications", {
  cells <- function(p, w = c(1, 1), obs = c(0, 1)) summarize_marginal_fit_grid(NULL, obs, p, w, 0:1)
  good <- cells(matrix(0.5, 2, 2))
  expect_equal(good$cell_stats$ExpectedCount, c(1, 1))
  expect_equal(good$cell_stats$VarianceCount, c(0.5, 0.5))
  expect_false(good$summary_stats$Flagged)
  for (p in list(matrix(NA_real_, 2, 2), matrix(c(NA, 0.5, 0.5, 0.5), 2),
                 matrix(0.6, 2, 2), matrix(c(-0.1, 0.5, 1.1, 0.5), 2))) {
    bad <- cells(p)
    expect_true(all(is.na(bad$cell_stats$ExpectedCount)))
    expect_true(all(is.na(bad$cell_stats$FlaggedAbsZ)))
    expect_true(is.na(bad$summary_stats$RMSD))
    expect_true(is.na(bad$summary_stats$FlaggedCellCount))
    expect_true(is.na(bad$summary_stats$Flagged))
    expect_equal(bad$summary_stats$UnclassifiedCells, 2)
  }
  for (w in list(c(NA, 1), c(-1, 1))) {
    bad <- cells(matrix(0.5, 2, 2), w)
    expect_true(all(is.na(bad$cell_stats$ObservedCount)))
    expect_true(all(is.na(bad$cell_stats$FlaggedAbsZ)))
  }
  zero <- cells(rbind(c(NA, NA), c(0.5, 0.5)), c(0, 2))
  expect_equal(zero$cell_stats$ExpectedCount, c(1, 1))
  expect_equal(zero$summary_stats$UnavailableObservations, 0)
  expect_true(is.na(cells(matrix(0.5, 2, 2), c(0, 0))$summary_stats$Flagged))
  missing_score <- cells(matrix(0.5, 2, 2), obs = c(NA, 1))
  expect_true(all(is.na(missing_score$cell_stats$ObservedCount)))
  expect_equal(missing_score$cell_stats$ExpectedCount, c(1, 1))
  expect_error(cells(matrix(0.5, 1, 2)), "must match")
})

test_that("pairwise marginal screens keep unavailable context opportunities", {
  core <- marginal_pair_fixture()
  complete <- calc_marginal_pairwise_bundle(core)
  expect_equal(complete$pair_stats$LevelPairCount, 2)
  expect_equal(complete$pair_stats$ExpectedExactCount, 1)
  # Adjacent agreement is deterministic on a binary scale, so its z score is unavailable.
  expect_true(is.na(complete$pair_stats$AdjacentStdResidual))
  expect_true(is.na(complete$pair_stats$FlaggedAdjacent))
  expect_true(complete$pair_stats$FlaggedExact)
  expect_true(complete$pair_stats$Flagged)
  core$posterior_bundle$obs_posterior[3:4, ] <- NA_real_
  partial <- calc_marginal_pairwise_bundle(core)
  expect_equal(partial$pair_stats$LevelPairCount, 2)
  expect_equal(partial$pair_stats$AvailableContextPairs, 1)
  expect_equal(partial$pair_stats$UnavailableContextPairs, 1)
  expect_true(is.na(partial$pair_stats$ExpectedExactCount))
  expect_true(is.na(partial$pair_stats$Flagged))
  expect_true(is.na(partial$facet_summary$FlaggedLevelPairs))
  expect_equal(partial$coverage$Unclassified, 1)
  core <- marginal_pair_fixture()
  core$weights[1] <- NA_real_
  expect_equal(calc_marginal_pairwise_bundle(core)$pair_stats$UnavailableContextPairs, 1)
  core$weights[1:2] <- 0
  expect_equal(calc_marginal_pairwise_bundle(core)$pair_stats$LevelPairCount, 1)
  core <- marginal_pair_fixture()
  core$logprob_bundle$prob_list[[1]][3, 1] <- NA_real_
  expect_equal(calc_marginal_pairwise_bundle(core)$pair_stats$UnavailableContextPairs, 1)
})

test_that("marginal coverage reaches summaries, reports, exports and saved-result guards", {
  # Detailed propagation remains in the extended local/CI suite.
  skip_on_cran()
  fit <- make_toy_fit(method = "MML")
  diag <- make_toy_diagnostics(fit, diagnostic_mode = "both")
  core <- compute_mml_expected_category_diagnostics(fit)
  core$expected_bundle$posterior_prob[,] <- NA_real_
  core$posterior_bundle$obs_posterior[,] <- NA_real_
  testthat::local_mocked_bindings(compute_mml_expected_category_diagnostics = function(...) core)
  diag$marginal_fit <- calc_marginal_fit_bundle(fit)
  expect_true(all(is.na(diag$marginal_fit$coverage$Flagged)))
  s <- summary(diag)
  expect_equal(s$marginal_coverage, diag$marginal_fit$coverage)
  expect_equal(build_summary_table_bundle(s)$tables$marginal_coverage, as.data.frame(s$marginal_coverage))
  expect_true(all(is.na(s$flags$Count[grepl("Marginal", s$flags$Metric)])))
  printed <- gsub("[[:space:]]+", " ", paste(capture.output(print(summary(diag, detail = "full"))), collapse = " "))
  expect_match(printed, "no classifications available")
  expect_match(printed, "some classifications unavailable")
  expect_false(grepl("partly_classified", printed))
  expect_false(grepl("limited_information_inspired|latent_integrated_first_order_counts|screening_only|generalized_residual_logic", printed))
  contract <- build_apa_reporting_contract(fit, diag)
  text <- gsub("[[:space:]]+", " ", paste(contract$report_text, collapse = " "))
  expect_match(text, "no classifications available")
  expect_false(grepl("largest strict marginal cell|largest strict pairwise signal", text))
  expect_true(is.na(rating_scale_table(fit, diag)$summary$MarginalFlaggedCategories))
  expect_equal(rating_scale_table(fit, diag)$summary$MarginalUnclassifiedCategories, fit$config$n_cat)
  expect_true(is.na(summary(category_structure_report(fit, diag))$summary$MarginalFlaggedCategories))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- plot_marginal_fit(diag, draw = TRUE)
  expect_equal(p$data$retention$Available, 0)
  expect_equal(p$data$retention$Unavailable, p$data$retention$Candidates)
  expect_true(all(is.na(p$data$table$StdResidual)))
  expect_match(p$data$subtitle, "unavailable")
  expect_no_error(plot_marginal_pairwise(diag, draw = TRUE))
  mismatch <- fit
  mismatch$facets$others$Estimate[1] <- mismatch$facets$others$Estimate[1] + 1
  expect_error(plot_marginal_fit(mismatch, diagnostics = diag, draw = FALSE), "fit/diagnostics mismatch", fixed = TRUE)
  diag$marginal_fit$coverage_version <- NULL
  expect_error(summary(diag), "Recreate diagnostics")
  expect_error(plot_marginal_fit(diag, draw = FALSE), "Recreate diagnostics")
  expect_error(rating_scale_table(fit, diag), "Recreate diagnostics")
  expect_error(build_apa_reporting_contract(fit, diag), "Recreate diagnostics")
  s$marginal_coverage <- NULL
  expect_error(print(s), "Recreate")
  expect_error(build_summary_table_bundle(s), "marginal_coverage")
})

test_that("marginal plots rank the full requested facet and metric", {
  fit <- make_toy_fit(method = "MML")
  diag <- make_toy_diagnostics(fit, diagnostic_mode = "both")
  diag$marginal_fit$top_cells <- diag$marginal_fit$top_cells[0, ]
  diag$marginal_fit$pairwise$top_pairs <- diag$marginal_fit$pairwise$top_pairs[0, ]
  p <- plot_marginal_fit(diag, plot_type = "prop_diff", facet = "Rater", top_n = 100, draw = FALSE)
  expect_true(all(p$data$table$Facet == "Rater"))
  expect_equal(nrow(p$data$table), sum(diag$marginal_fit$facet_level$cell_stats$Facet == "Rater"))
  pair <- plot_marginal_pairwise(diag, metric = "adjacent", facet = "Rater", top_n = 100, draw = FALSE)
  expect_equal(nrow(pair$data$table), sum(diag$marginal_fit$pairwise$pair_stats$Facet == "Rater"))
})

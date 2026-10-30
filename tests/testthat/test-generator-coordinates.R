test_that("normal latent draws are assigned to IDs without severity ordering", {
  seed <- 542
  reference <- mfrmr:::with_preserved_rng_seed(seed, list(
    person = rnorm(7, 0, .8), rater = rnorm(4, 0, .4), criterion = rnorm(4, 0, .6)))
  for (model in c("RSM", "PCM", "GPCM")) {
    sim <- simulate_mfrm_data(n_person = 7, n_rater = 4, n_criterion = 4,
      theta_sd = .8, rater_sd = .4, criterion_sd = .6, model = model, seed = seed)
    truth <- attr(sim, "mfrm_truth")
    expect_equal(unname(truth$person), reference$person, tolerance = 0)
    expect_equal(unname(truth$facets$Rater), reference$rater, tolerance = 0)
    expect_equal(unname(truth$facets$Criterion), reference$criterion, tolerance = 0)
    expect_true(is.unsorted(truth$facets$Rater))
    expect_true(is.unsorted(truth$facets$Criterion))
  }
})

test_that("generator probabilities follow independently evaluated adjacent logits", {
  captured <- NULL
  implementations <- lapply(c("category_prob_rsm", "category_prob_pcm", "category_prob_gpcm"),
    function(n)get(n, envir = asNamespace("mfrmr")))
  local_mocked_bindings(
    category_prob_rsm = function(...) { captured <<- implementations[[1]](...); captured },
    category_prob_pcm = function(...) { captured <<- implementations[[2]](...); captured },
    category_prob_gpcm = function(...) { captured <<- implementations[[3]](...); captured },
    .package = "mfrmr")
  for (model in c("RSM", "PCM", "GPCM")) {
    args <- list(n_person = 6, n_rater = 3, n_criterion = 3, score_levels = 3,
      model = model, seed = 73)
    args$thresholds <- if (model == "RSM") c(-1, 2) else
      list(C01 = c(-1, 3), C02 = c(-3, 0), C03 = c(-.5, 1))
    if (model == "GPCM") args$slopes <- c(.5, 1, 2)
    sim <- do.call(simulate_mfrm_data, args)
    truth <- attr(sim, "mfrm_truth")
    expected <- t(vapply(seq_len(nrow(sim)), function(i) {
      row <- sim[i, ]
      eta <- truth$person[row$Person] - truth$facets$Rater[row$Rater] -
        truth$facets$Criterion[row$Criterion]
      key <- if (model == "RSM") "Common" else row$Criterion
      steps <- truth$step_table$Estimate[truth$step_table$StepFacet == key]
      slope <- if (model == "GPCM")
        truth$slope_table$Estimate[truth$slope_table$SlopeFacet == row$Criterion] else 1
      weights <- c(1, cumprod(exp(slope * (eta - steps))))
      weights / sum(weights)
    }, numeric(3)))
    expect_equal(unname(captured), unname(expected), tolerance = 1e-13)
    expect_equal(unname(rowSums(captured)), rep(1, nrow(sim)), tolerance = 1e-14)
    expect_true(all(sim$Score %in% 1:3))
  }
})

# Structural fixtures isolate nesting from numerical solution qualification.
rsm_pcm_nesting_fixture <- function(rsm_interactions = NULL, pcm_interactions = NULL) {
  data <- expand.grid(Person = paste0("P", 1:12), Rater = paste0("R", 1:3),
    Criterion = paste0("C", 1:4), Task = paste0("T", 1:2))
  data$Score <- (seq_len(nrow(data)) - 1L) %% 5L
  prep <- prepare_mfrm_data(data, "Person", c("Rater", "Criterion", "Task"), "Score")
  Map(function(model, interactions) {
    specs <- build_facet_interaction_specs(prep, interactions, interaction_policy = "silent")
    config <- build_estimation_config(prep, model, "MML",
      step_facet = if (model == "PCM") "Criterion" else NULL, slope_facet = NULL,
      interaction_specs = specs, weight_col = NULL,
      facet_signs = build_facet_signs(prep$facet_names)$signs,
      positive_facets = character(), noncenter_facet = "Person",
      dummy_facets = character(), anchor_df = NULL, group_anchor_df = NULL)$config
    config$estimation_control <- list(quad_points = 31L, mml_integration = "fixed")
    par <- rep(0, sum(unlist(build_param_sizes(config))))
    opt <- list(par = par, value = if (model == "RSM") 100 else 97)
    contract <- build_mfrm_ic_contract(-opt$value, length(par), prep, config, "MML")
    structure(list(config = config, prep = prep, opt = opt, summary = contract),
      class = "mfrm_fit")
  }, c("RSM", "PCM"), list(rsm_interactions, pcm_interactions))
}

test_that("RSM/PCM nesting requires inclusion of the RSM interaction terms", {
  cases <- list(
    list(NULL, NULL, TRUE),
    list(NULL, "Rater:Criterion", TRUE),
    list("Rater:Criterion", "Criterion:Rater", TRUE),
    list("Rater:Criterion", c("Rater:Criterion", "Rater:Task"), TRUE),
    list("Rater:Criterion", NULL, FALSE),
    list("Rater:Criterion", "Rater:Task", FALSE),
    list(c("Rater:Criterion", "Rater:Task"), "Rater:Criterion", FALSE))
  for (case in cases) {
    pair <- rsm_pcm_nesting_fixture(case[[1]], case[[2]])
    for (reverse in c(FALSE, TRUE)) {
      fits <- if (reverse) rev(pair) else pair
      review <- audit_compare_mfrm_nesting(fits, c("first", "second"))
      expect_identical(review$eligible, case[[3]])
      if (case[[3]]) {
        expect_identical(review$relation, "RSM_in_PCM")
      } else {
        expect_identical(review$relation, "RSM_PCM_interactions_not_nested")
        expect_match(review$reason, "RSM interaction terms")
      }
    }
  }
})

test_that("the rejected RSM interaction direction lies outside the PCM response map", {
  design <- function(fit) {
    config <- fit$config
    idx <- build_indices(fit$prep, step_facet = config$step_facet,
      interaction_specs = config$interaction_specs)
    as.matrix(mfrmr_estimability_adjacent_design(fit$prep, idx, config,
      build_param_sizes(config))$design)
  }
  bad <- rsm_pcm_nesting_fixture("Rater:Criterion", NULL)
  x <- design(bad[[1]]); y <- design(bad[[2]])
  expect_lt(ncol(x), ncol(y))
  expect_gt(qr(cbind(y, x), tol = 1e-10)$rank, qr(y, tol = 1e-10)$rank)
  good <- rsm_pcm_nesting_fixture("Rater:Criterion", c("Rater:Criterion", "Rater:Task"))
  x <- design(good[[1]]); y <- design(good[[2]])
  expect_equal(qr(cbind(y, x), tol = 1e-10)$rank, qr(y, tol = 1e-10)$rank)
})

test_that("a nonnested RSM/PCM pair retains IC comparison but has no chi-square test", {
  local_mocked_bindings(
    mfrm_convergence_state = function(...) list(code_converged = TRUE, inference_ready = TRUE, severity = "pass"),
    mfrm_ic_fit_check = function(...) list(eligible = TRUE, basis = "test", review = "test"),
    mfrmr_get_readiness_record = function(...) list(fit = data.frame(ReasonCodes = "test")),
    .package = "mfrmr")
  pair <- rsm_pcm_nesting_fixture("Rater:Criterion", NULL)
  for (reverse in c(FALSE, TRUE)) {
    fits <- if (reverse) rev(pair) else pair
    expect_warning(result <- do.call(compare_mfrm, c(fits, list(nested = TRUE))),
      "RSM interaction terms")
    expect_true(all(result$table$ICComparable))
    expect_true(all(is.finite(result$table$Delta_AIC)))
    expect_null(result$lrt)
    expect_identical(result$comparison_basis$lrt_status, "not_computed")
    expect_match(result$comparison_basis$lrt_reason, "RSM interaction terms")
    expect_match(paste(capture.output(print(result)), collapse = " "), "RSM interaction terms")
    expect_identical(summary(result)$comparison_basis, result$comparison_basis)
  }
  good <- rsm_pcm_nesting_fixture("Rater:Criterion", c("Rater:Criterion", "Rater:Task"))
  result <- do.call(compare_mfrm, c(good, list(nested = TRUE)))
  expect_identical(result$comparison_basis$lrt_status, "computed")
  expect_equal(result$lrt$ChiSq, 6)
  expect_equal(result$lrt$df, length(good[[2]]$opt$par) - length(good[[1]]$opt$par))
})

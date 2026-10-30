fitted_comparison_cases <- function() readRDS(test_path("fixtures", "fitted-comparisons.rds"))

test_that("joint-fit contrasts agree with independently weighted interaction cells", {
  cases <- fitted_comparison_cases()
  local_mocked_bindings(fit_mfrm = function(...) stop("No refit"), .package = "mfrmr")
  for (name in setdiff(names(cases), "coarse")) {
    x <- cases[[name]]
    out <- interaction_effect_table(x$fit, comparison_review = x$review)
    cells <- interaction_effect_table(x$fit)
    question <- x$review$comparison_review
    weights <- question$weights
    # Cell labels and weights supply an independent path from free coordinates.
    rater <- if (cells$FacetA[1] == "Rater") cells$FacetA_Level else cells$FacetB_Level
    group <- if (cells$FacetA[1] == "Group") cells$FacetA_Level else cells$FacetB_Level
    expected <- vapply(seq_len(nrow(weights)), function(i)
      sum(weights[i, rater] * ifelse(group == "A", 1, -1) * cells$Estimate), numeric(1))
    expect_equal(out$Estimate, expected, tolerance = 1e-12)
    expect_true(all(out$EstimateStatus == "descriptive_only"))
    expect_false(any(out$FormalInferenceEligible))
    expect_true(all(out$IdentificationBasis == "free_person_response_model"))
    expect_true(all(is.na(out$SE) & is.na(out$Lower) & is.na(out$Upper)))
    expect_equal(attr(out, "comparison_review")$weights, weights)
    expect_true(all(is.na(attr(out, "comparison_review")$table$Estimate)))
    if (x$fit$config$method == "MML") {
      expect_true(all(out$State == "not_evaluated"))
      expect_true(all(out$QuadratureState == "stable_at_retained_point"))
      expect_equal(attr(out, "fit_checks")$quadrature, c(63, 127))
    }
    reversed <- question; reversed$specification$groups <- c("B", "A")
    expect_equal(interaction_effect_table(x$fit, reversed)$Estimate, -out$Estimate, tolerance = 1e-12)
  }
})

test_that("fitted comparisons preserve output and source identity through summary and storage", {
  x <- fitted_comparison_cases()[["RSM-JML"]]
  table <- interaction_effect_table(x$fit, x$review)
  result <- summary(x$fit, comparison_review = x$review, compute = "never")
  expect_identical(result$interaction_comparisons, table)
  expect_output(print(result), "Requested joint-fit interaction comparisons")
  expect_output(print(summary(x$fit, comparison_review = x$review, detail = "brief")), "descriptive_only")
  expect_null(summary(x$fit)$interaction_comparisons)
  restored <- unserialize(serialize(result, NULL))
  expect_identical(restored$interaction_comparisons, table)
  shuffled <- x$fit
  shuffled$prep$data <- shuffled$prep$data[rev(seq_len(nrow(shuffled$prep$data))), ]
  expect_equal(interaction_effect_table(shuffled, x$review)$Estimate, table$Estimate, tolerance = 1e-12)
  review <- x$review$comparison_review
  review$specification$contrasts <- review$specification$contrasts[, 4:1]
  expect_equal(interaction_effect_table(x$fit, review)$Estimate, table$Estimate, tolerance = 1e-12)
})

test_that("model or retained-data mismatches cannot borrow a prospective review", {
  x <- fitted_comparison_cases()[["PCM-MML"]]
  review <- x$review$comparison_review
  for (field in c("model", "method", "step_facet", "population_formula", "positive_facets")) {
    bad <- review
    replacement <- list(model="RSM", method="JML", step_facet="Rater",
      population_formula=~1, positive_facets="Rater")
    bad$specification[[field]] <- replacement[[field]]
    expect_error(interaction_effect_table(x$fit, bad))
  }
  bad <- review; bad$source$rows$Score[1] <- bad$source$rows$Score[1] + 1
  expect_error(interaction_effect_table(x$fit, bad), "retained responses")
  bad <- review; bad$source$rows$Group[1] <- "changed"
  expect_error(interaction_effect_table(x$fit, bad), "retained responses")
  bad <- x$fit; bad$interactions$effects$Estimate[1] <- bad$interactions$effects$Estimate[1] + .1
  expect_error(interaction_effect_table(bad, review), "Saved interaction")
  bad <- x$fit; bad$config$facet_shrinkage <- "empirical_bayes"
  expect_error(interaction_effect_table(bad, review), "without shrinkage")
  expect_error(interaction_effect_table(x$fit, list()), "Supply the comparison_review")
})

test_that("unready and coarse-grid fits retain unavailable rows and optimizer traces", {
  cases <- fitted_comparison_cases()
  x <- cases[["RSM-JML"]]
  for (state in c("blocked", "legacy_unknown", "ready_with_exclusions")) {
    bad <- x$fit; bad$readiness$fit$FitReadiness <- state
    out <- interaction_effect_table(bad, x$review)
    expect_true(all(is.na(out$Estimate)))
    expect_true(all(is.finite(out$OptimizerEstimate)))
  }
  bad <- x$fit; bad$opt$convergence <- 1L
  expect_true(all(is.na(interaction_effect_table(bad, x$review)$Estimate)))
  bad <- x$fit; bad$opt$value <- bad$opt$value + 1
  expect_true(all(is.na(interaction_effect_table(bad, x$review)$Estimate)))
  # Alter a Person parameter, leaving the interaction cells unchanged.
  bad <- x$fit; bad$opt$par[1] <- bad$opt$par[1] + .1
  out <- interaction_effect_table(bad, x$review)
  expect_true(all(is.na(out$Estimate)))
  expect_match(out$EstimateReason[1], "objective/gradient")
  coarse <- cases$coarse
  out <- interaction_effect_table(coarse$fit, coarse$review)
  expect_true(all(is.na(out$Estimate)))
  expect_true(all(out$QuadratureState == "review"))
  expect_true(all(is.finite(out$OptimizerEstimate)))
  expect_false(any(out$FormalInferenceEligible))
})

test_that("a partially identified target does not relax the full-fit guard", {
  x <- fitted_comparison_cases()[["RSM-JML"]]
  d <- subset(x$data, !(Rater == "R04" & Group == "B"))
  review <- describe_mfrm_data(d, "Person", c("Rater", "Criterion", "Group"), "Score",
    rating_min=1, rating_max=4, category_policy="preserve", include_agreement=FALSE,
    comparisons=x$review$comparison_review$specification)
  expect_identical(review$comparison_review$table$State,
    c("target_unidentified", "target_identified_only"))
  expect_error(do.call(fit_mfrm, c(list(data=d, interaction_policy="silent"),
    review$comparison_review$fit_arguments)), class="mfrmr_estimability_error")
  expect_error(interaction_effect_table(x$fit, review), "retained responses")
})

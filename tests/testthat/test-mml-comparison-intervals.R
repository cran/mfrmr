mml_comparison_intervals_fixture <- function() readRDS(test_path("fixtures", "mml-comparison-intervals.rds"))

mml_comparison_saved_information <- function(ci, eligible = TRUE) list(
  covariance = if (eligible) ci$parameter_covariance else NULL,
  information = list(status = if (eligible) "ok" else "unavailable"),
  checks = ci$checks, numerical_checks = ci$numerical_checks,
  check = list(eligible = eligible, review = if (eligible) "Saved joint information qualification" else "Finer-grid information failed"))

test_that("interaction contrasts use the complete population-SD derivative and joint covariance", {
  local_mocked_bindings(fit_mfrm = function(...) stop("No refit"), .package = "mfrmr")
  for (x in mml_comparison_intervals_fixture()) {
    local_mocked_bindings(mfrm_native_location_inference = function(...) mml_comparison_saved_information(x$intervals), .package = "mfrmr")
    regression <- length(x$fit$config$population_spec$design_columns) > 1L
    scale <- if (regression) "native" else "standardized"
    ci <- suppressWarnings(mfrm_facet_intervals(x$fit, comparison_review=x$review, scale=scale, level=.9))
    sizes <- build_param_sizes(x$fit$config)
    spec <- x$review$comparison_review$specification
    W <- spec$contrasts
    # Differentiate weighted expanded cells, not the free-coordinate target map.
    value <- function(z) {
      p <- expand_params(z, sizes, x$fit$config)
      cells <- build_interaction_effect_table(x$fit$config, x$fit$prep, p)
      rater <- if (cells$FacetA[1] == spec$facet) cells$FacetA_Level else cells$FacetB_Level
      group <- if (cells$FacetA[1] == spec$group) cells$FacetA_Level else cells$FacetB_Level
      d <- vapply(seq_len(nrow(W)), function(i) sum(W[i, rater] *
        ((group == spec$groups[1]) - (group == spec$groups[2])) * cells$Estimate), numeric(1))
      if (scale == "standardized") d / sqrt(p$population$sigma2 %||% 1) else d
    }
    h <- 1e-5
    J <- vapply(seq_along(x$fit$opt$par), function(j) {
      step <- numeric(length(x$fit$opt$par)); step[j] <- h
      (value(x$fit$opt$par + step) - value(x$fit$opt$par - step)) / (2*h)
    }, numeric(nrow(W)))
    expect_equal(unname(ci$target_jacobian), J, ignore_attr=TRUE, tolerance=1e-8)
    expect_equal(ci$table$Estimate, value(x$fit$opt$par), tolerance=1e-12)
    expect_equal(ci$covariance, J %*% x$intervals$parameter_covariance %*% t(J),
      ignore_attr=TRUE, tolerance=1e-8)
    expect_equal(ci$table$Upper-ci$table$Estimate, qnorm(.95)*ci$table$SE)
    expect_true(all(ci$table$CIEligible))
    expect_false(any(ci$table$CoverageQualified))
    expect_true(all(ci$table$MarginalIdentificationState == "not_evaluated"))
    expect_true(all(ci$table$IdentificationBasis == "free_person_response_model"))
    expect_identical(ci$table$IdentificationDetail, ci$comparison_review$table$IdentificationDetail)
    expect_identical(ci$fit, x$fit)
    expect_identical(ci$settings$procedure, "mml_interaction_comparison_model_v1")
    if (regression) expect_error(mfrm_facet_intervals(x$fit, comparison_review=x$review,
      scale="standardized"), "residual SD is not marginal")
    if (!regression && isTRUE(x$fit$population$active)) {
      native <- suppressWarnings(mfrm_facet_intervals(x$fit, comparison_review=x$review))
      sigma <- sqrt(x$fit$population$sigma2)
      lambda <- build_param_slices(sizes)$log_sigma2
      V <- ci$parameter_covariance; d <- native$table$Estimate
      covariance <- drop(native$target_jacobian %*% V[,lambda])
      expected <- (diag(native$covariance) + d^2*V[lambda,lambda]/4 - d*covariance)/sigma^2
      expect_equal(ci$table$SE^2, unname(expected), tolerance=1e-10)
      expect_gt(max(abs(ci$table$SE - native$table$SE/sigma)), 1e-4)
      expect_equal(ci$settings$population_sd, sigma)
      expect_equal(ci$table$NativeEstimate, native$table$Estimate)
    }
  }
})

test_that("review direction, ordering and correlated target rows remain explicit", {
  x <- mml_comparison_intervals_fixture()$RSM
  local_mocked_bindings(mfrm_native_location_inference = function(...) mml_comparison_saved_information(x$intervals), .package = "mfrmr")
  review <- x$review$comparison_review
  review$specification$groups <- rev(review$specification$groups)
  W <- review$specification$contrasts
  review$specification$contrasts <- rbind(W[,4:1], double=2*W[2,4:1])
  review$target_map[,] <- 0
  ci <- suppressWarnings(mfrm_facet_intervals(x$fit, comparison_review=review, scale="standardized"))
  expect_equal(ci$table$Estimate[1:2], -x$intervals$table$Estimate, tolerance=1e-12)
  expect_equal(ci$table$Lower[1:2], -x$intervals$table$Upper, tolerance=1e-12)
  expect_equal(ci$table$SE[3], 2*ci$table$SE[2], tolerance=1e-12)
  expect_equal(ci$covariance[3,], 2*ci$covariance[2,], ignore_attr=TRUE, tolerance=1e-12)
  expect_true(all(ci$table$Group1 == "B" & ci$table$Group2 == "A"))
})

test_that("invalid sources and failed qualification cannot yield interaction intervals", {
  x <- mml_comparison_intervals_fixture()$RSM
  local_mocked_bindings(mfrm_native_location_inference = function(...) stop("No information for rejected sources"), .package = "mfrmr")
  bad <- x$fit; bad$interactions$effects$Estimate[1] <- bad$interactions$effects$Estimate[1]+.1
  expect_error(mfrm_facet_intervals(bad, comparison_review=x$review), "Saved interactions")
  bad <- x$fit; bad$population$sigma2 <- bad$population$sigma2*2
  expect_error(mfrm_facet_intervals(bad, comparison_review=x$review), "population estimates")
  bad <- x$fit; bad$prep$data$Score[1] <- 99
  expect_error(mfrm_facet_intervals(bad, comparison_review=x$review), "retained responses")
  bad <- x$fit; bad$config$estimation_control$mml_engine_used <- "em"
  expect_error(mfrm_facet_intervals(bad, comparison_review=x$review), "direct MML")
  expect_error(mfrm_facet_intervals(x$fit, "Criterion", comparison_review=x$review), "review's contrast")
  expect_error(mfrm_facet_intervals(x$fit, contrasts=x$review$comparison_review$weights,
    comparison_review=x$review), "review's contrast")
  expect_error(mfrm_facet_intervals(x$fit, comparison_review=x$review, method="sandwich"), "method = 'model'")
  expect_error(mfrm_facet_intervals(x$fit, comparison_review=x$review, level=1), "0 < level < 1")
  expect_error(mfrm_facet_intervals(x$fit, comparison_review=x$review, clusters=data.frame()), "clusters")
  partial <- readRDS(test_path("fixtures", "partial-comparisons.rds"))[["RSM-MML"]]
  expect_error(mfrm_facet_intervals(partial$fit, comparison_review=partial$review), "full native")
  local_mocked_bindings(mfrm_native_location_inference = function(...) mml_comparison_saved_information(x$intervals, FALSE), .package = "mfrmr")
  ci <- suppressWarnings(mfrm_facet_intervals(x$fit, comparison_review=x$review, scale="standardized"))
  expect_equal(ci$table$Estimate, x$intervals$table$Estimate)
  expect_true(all(is.na(ci$table$SE) & is.na(ci$table$Lower) & is.na(ci$table$Upper)))
  expect_true(all(ci$table$Status == "numerical_review_failed"))
  expect_false(any(ci$table$CIEligible))
})

test_that("saved interaction intervals reach plots, reports and exports without recomputation", {
  cases <- mml_comparison_intervals_fixture()
  local_mocked_bindings(fit_mfrm=function(...)stop("No refit"),
    compute_mml_parameter_covariance=function(...)stop("No new information"), .package="mfrmr")
  for (name in c("RSM", "RSM-regression")) {
    x <- cases[[name]]; ci <- x$intervals
    p <- plot(ci, draw=FALSE)
    expect_identical(plot_data(p, "table"), ci$table)
    expect_match(plot_data(p,"subtitle"), "interaction differences")
    expect_false(grepl("Weak information", plot_data(p,"subtitle")))
    if (name == "RSM") expect_identical(plot_data(p,"xlab"), "Contrast (population SD units)")
    else expect_identical(plot_data(p,"xlab"), "Interaction contrast (native ability units)")
    if (requireNamespace("ggplot2",quietly=TRUE)) expect_s3_class(as_ggplot(p), "ggplot")
    out <- mfrm_results(x$fit, intervals=list(school=ci), include=c("fit","plots"), compute="never")
    expect_identical(out$tables$facet_school_comparison_review, ci$comparison_review$table)
    expect_identical(out$tables$facet_school_comparison_cells, ci$comparison_review$cells)
    expect_identical(plot_data(out,type="facet_school",component="table"),ci$table)
    expect_identical(mfrm_report(out)$tables$facet_school_intervals$ScaleReference, ci$table$ScaleReference)
    expect_identical(unserialize(serialize(out,NULL))$facet_intervals$school,ci)
    bad <- x$fit; bad$interactions$effects$Estimate[1] <- bad$interactions$effects$Estimate[1]+.1
    expect_error(mfrm_results(bad, intervals=ci, compute="never"), "source coordinates")
  }
})

test_that("new identification explanations reach saved MML interval reports", {
  cases <- mml_comparison_intervals_fixture()
  local_mocked_bindings(fit_mfrm=function(...)stop("No refit"), .package="mfrmr")
  for (name in c("RSM", "RSM-regression")) {
    x <- cases[[name]]
    local_mocked_bindings(mfrm_native_location_inference=function(...)
      mml_comparison_saved_information(x$intervals), .package="mfrmr")
    ci <- suppressWarnings(mfrm_facet_intervals(x$fit, comparison_review=x$review,
      scale=if(name=="RSM")"standardized" else "native"))
    expect_identical(ci$table$Estimate,x$intervals$table$Estimate)
    # Recomputed matrix products can differ in their last bits across platforms.
    expect_equal(ci$covariance,x$intervals$covariance,
      tolerance=64*.Machine$double.eps)
    expect_identical(ci$table$IdentificationBasis,rep("free_person_response_model",2))
    out <- mfrm_results(x$fit,intervals=list(school=ci),include=c("fit","plots"),compute="never")
    expect_identical(out$tables$facet_school_intervals$IdentificationBasis,ci$table$IdentificationBasis)
    expect_identical(mfrm_report(out)$tables$facet_school_intervals$IdentificationDetail,ci$table$IdentificationDetail)
    expect_identical(plot_data(out,type="facet_school",component="table")$IdentificationBasis,ci$table$IdentificationBasis)
    expect_identical(unserialize(serialize(out,NULL))$facet_intervals$school,ci)
  }
})

context_curve_fits <- function() {
  list(readRDS(test_path("fixtures", "native-location-fixed-gpcm.rds")),
    readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))$fit,
    readRDS(test_path("fixtures", "gmfrm-explicit-start.rds"))$fit)
}

test_that("explicit GPCM contexts preserve the full predictor and cumulative target", {
  local_mocked_bindings(fit_mfrm = function(...) stop("No refitting"),
    compute_mml_parameter_covariance = function(...) stop("No covariance"),
    diagnose_mfrm = function(...) stop("No diagnostics"))
  for (fit in context_curve_fits()) {
    nd <- unique(fit$prep$data[fit$config$facet_names])
    nd <- head(nd, 2)
    out <- category_curves_report(fit, newdata = nd, theta_points = 51, digits = 12)
    expect_identical(out$newdata, as.data.frame(nd))
    expect_true(all(out$curve_profiles$ObservedContext))
    expect_identical(out$settings$curve_basis, "fitted_rating_context")
    expect_true(isTRUE(out$settings$point_only))
    bounds <- out$cumulative_boundaries
    expect_true(all(bounds$BoundaryStatus == "in_range"))
    expect_true(all(abs(bounds$RootResidual) < 1e-8))
    expect_true(all(is.na(bounds$SE) & is.na(bounds$CI_Lower) & is.na(bounds$CI_Upper)))
    expect_false(any(bounds$CIEligible))
    for (i in seq_len(nrow(out$curve_profiles))) {
      group <- out$curve_profiles$CurveGroup[i]
      # Independent softmax using the publicly stored facet, step and slope tables.
      locations <- fit$facets$others
      offset <- sum(vapply(fit$config$facet_names, function(facet) {
        fit$config$facet_signs[[facet]] * locations$Estimate[
          locations$Facet == facet & locations$Level == as.character(nd[[facet]][i])]
      }, numeric(1)))
      slopes <- fit$slopes
      a <- prod(vapply(fit$config$slope_facet, function(owner) {
        keep <- slopes$SlopeFacet == as.character(nd[[owner]][i])
        if ("SlopeOwner" %in% names(slopes)) keep <- keep & slopes$SlopeOwner == owner
        slopes$Estimate[keep]
      }, numeric(1)))
      steps <- fit$steps[fit$steps$StepFacet == as.character(nd[[fit$config$step_facet]][i]), ]
      steps <- steps[order(as.integer(sub("Step_", "", steps$Step))), ]
      tau <- c(0, cumsum(steps$Estimate))
      probability <- function(theta) {
        eta <- a * ((seq_along(tau) - 1) * (theta + offset) - tau)
        exp(eta - max(eta)) / sum(exp(eta - max(eta)))
      }
      tab <- out$probabilities[out$probabilities$CurveGroup == group, ]
      at_zero <- tab[tab$Theta == 0, ]
      expect_equal(at_zero$Probability, probability(0), tolerance = 1e-11)
      expect_equal(out$curve_profiles$EffectiveSlope[i], a)
      expect_equal(out$curve_profiles$PredictorOffset[i], offset)
      b <- bounds[bounds$CurveGroup == group, ]
      for (j in seq_len(nrow(b))) {
        expect_equal(sum(probability(b$ThurstonianThreshold[j])[seq_len(j)]), .5, tolerance = 1e-8)
        adjacent <- probability(b$AdjacentThreshold[j])
        expect_equal(adjacent[j], adjacent[j + 1], tolerance = 1e-10)
      }
      expect_true(all(diff(b$ThurstonianThreshold) > 0))
      info <- out$category_information[out$category_information$CurveGroup == group, ]
      expect_equal(aggregate(CategoryInformation ~ Theta, info, sum)$CategoryInformation,
        out$expected_ogive$Information[out$expected_ogive$CurveGroup == group], tolerance = 1e-10)
    }
    refined <- category_curves_report(fit, newdata = nd, theta_points = 71, digits = 2)
    expect_identical(refined$cumulative_boundaries, bounds)
    expect_match(paste(summary(out)$notes, collapse = " "), "Point calculations only")
    p <- plot(out, type = "cumulative", draw = FALSE)
    expect_identical(p$data$fit_readiness, out$fit_readiness)
    expect_identical(p$data$newdata, out$newdata)
    expect_identical(p$data$notes, out$notes)
    path <- tempfile(); on.exit(unlink(path), add = TRUE)
    saveRDS(out, path)
    expect_identical(readRDS(path), out)
    expect_identical(plot_data(readRDS(path), "cumulative_boundaries"), bounds)
  }
})

test_that("context boundaries handle endpoints, outside ranges, labels and weak source states", {
  fit <- context_curve_fits()[[3]]
  nd <- head(fit$prep$data[fit$config$facet_names], 1)
  spec <- mfrmr:::category_context_curve_spec(fit, nd)
  boundary <- unname(-spec$groups[[1]]$predictor_offset) # Binary step is centered at zero.
  endpoint <- category_curves_report(fit, newdata = nd, theta_range = c(boundary, boundary + 2), theta_points = 51)
  expect_equal(endpoint$cumulative_boundaries$ThurstonianThreshold, boundary, tolerance = 1e-10)
  outside <- category_curves_report(fit, newdata = nd, theta_range = c(boundary + 1, boundary + 2), theta_points = 51)
  expect_true(all(is.na(outside$cumulative_boundaries$ThurstonianThreshold)))
  expect_identical(outside$cumulative_boundaries$BoundaryStatus, "outside_theta_range")
  expect_identical(outside$cumulative_boundaries$RootMethod, "analytic_binary")
  expect_true(all(outside$cumulative_boundaries$LowerCumulativeProbability < .5))
  fit$summary$Converged <- FALSE
  fit$summary$InferenceReady <- FALSE
  fit$readiness$fit$NumericalState <- "failed"
  fit$readiness$fit$InferenceReady <- FALSE
  fit$prep$score_map$OriginalScore <- c(2, 8)
  out <- category_curves_report(fit, newdata = nd, theta_points = 51, include_fixed = TRUE)
  expect_match(out$fixed, "Point calculations only")
  expect_match(out$fixed, "Source convergence: FALSE")
  expect_true(is.finite(out$cumulative_boundaries$ThurstonianThreshold))
  expect_false(out$settings$source_converged)
  expect_false(out$settings$source_inference_ready)
  expect_identical(out$fit_readiness$NumericalState, "failed")
  expect_setequal(out$probabilities$CategoryLabel, c("2", "8"))
  expect_identical(out$cumulative_boundaries$UpperCategoryLabel, "8")
  expect_setequal(na.omit(plot(out, draw = FALSE)$data$plot_long$Category), c("2", "8"))
  expect_match(paste(summary(out)$notes, collapse = " "), "Source convergence: FALSE")
  expect_match(paste(summary(out)$notes, collapse = " "), "Two-family GPCM values are provisional")
})

test_that("context validation distinguishes unknown levels and unobserved combinations", {
  fit <- context_curve_fits()[[3]]
  nd <- head(fit$prep$data[fit$config$facet_names], 1)
  bad <- nd; bad[[1]] <- "unknown"
  expect_error(category_curves_report(fit, newdata = bad), "Unknown or missing")
  bad[[1]] <- NA_character_
  expect_error(category_curves_report(fit, newdata = bad), "Unknown or missing")
  expect_error(category_curves_report(fit, newdata = nd[-1]), "every non-Person facet")
  expect_error(category_curves_report(fit, newdata = transform(nd, Theta = 0)), "no Theta")
  expect_error(category_curves_report(fit, newdata = nd[FALSE, ]), "every non-Person facet")
  expect_error(category_curves_report(fit), "two|Two")
  pcm <- readRDS(test_path("fixtures", "native-location-pcm-regression.rds"))$fit
  expect_error(category_curves_report(pcm, newdata = nd), "native GPCM MML")
  rows <- Reduce(`&`, lapply(names(nd), function(n) fit$prep$data[[n]] == nd[[n]][1]))
  fit$prep$data <- fit$prep$data[!rows, ]
  out <- category_curves_report(fit, newdata = rbind(nd, nd), theta_points = 51)
  expect_false(any(out$curve_profiles$ObservedContext))
  expect_equal(length(unique(out$curve_profiles$CurveGroup)), 2)
  expect_true(all(is.finite(out$cumulative_boundaries$ThurstonianThreshold)))
  expect_match(paste(summary(out)$notes, collapse = " "), "not jointly observed")
})

test_that("tiny positive slopes retain correct boundaries without rounded probability roots", {
  fit <- context_curve_fits()[[3]]
  nd <- head(fit$prep$data[fit$config$facet_names], 1)
  out <- category_curves_report(fit, newdata = nd, theta_points = 51)
  spec <- mfrmr:::category_context_curve_spec(fit, nd)
  expected <- unname(-spec$groups[[1]]$predictor_offset)
  for (a in c(1e-12, 1e-16, 1e-20, 1e-300, 1e-320)) {
    spec$groups[[1]]$slope <- a
    binary <- mfrmr:::category_context_boundaries(out$cumulative_boundaries, spec, c(-6, 6), fit)
    expect_identical(binary$BoundaryStatus, "in_range")
    expect_equal(binary$ThurstonianThreshold, expected, tolerance = 1e-12)
    expect_identical(binary$RootMethod, "analytic_binary")
    if (a < 1e-300) next
    # Four-category central boundary: 1 + exp(.9a)*x = exp(.7a)*x^2 + x^3,
    # x = exp(a*theta). Independent high-precision roots tend to .05 + .03*a.
    poly_spec <- spec
    poly_spec$groups[[1]]$step_cum <- c(0, -.9, -.7, 0)
    poly_spec$groups[[1]]$predictor_offset <- 0
    bounds <- out$cumulative_boundaries
    bounds$BoundaryOrder <- 2L
    poly <- mfrmr:::category_context_boundaries(bounds, poly_spec, c(-6, 6), fit)
    expect_identical(poly$BoundaryStatus, "in_range")
    expect_equal(poly$ThurstonianThreshold, .05 + .03 * a, tolerance = 1e-12)
    expect_identical(poly$RootMethod, "paired_mass_uniroot")
  }
})

test_that("unrepresentable ability precision is reported instead of a false root", {
  fit <- context_curve_fits()[[3]]
  nd <- head(fit$prep$data[fit$config$facet_names], 1)
  out <- category_curves_report(fit, newdata = nd, theta_points = 51)
  spec <- mfrmr:::category_context_curve_spec(fit, nd)
  spec$groups[[1]]$predictor_offset <- -1e8
  failed <- mfrmr:::category_context_boundaries(out$cumulative_boundaries, spec, c(1e8 - 1, 1e8 + 1), fit)
  expect_identical(failed$BoundaryStatus, "numerical_failure")
  expect_true(is.na(failed$ThurstonianThreshold))
  expect_false(failed$InThetaRange)
  expect_match(failed$BoundaryDetail, "ability tolerance")
})

test_that("both cumulative directions display the same complementary 50 percent boundaries", {
  fit <- context_curve_fits()[[1]]
  nd <- head(fit$prep$data[fit$config$facet_names], 1)
  fit$prep$score_map$OriginalScore <- seq_len(fit$config$n_cat) * 10
  for (contexts in list(NULL, nd)) {
    out <- category_curves_report(fit, newdata = contexts, theta_points = 51)
    below <- plot(out, type = "cumulative", cumulative_direction = "at_or_below", draw = FALSE)
    above <- plot(out, type = "cumulative", cumulative_direction = "at_or_above", draw = FALSE)
    expect_gt(nrow(below$data$boundary_lines), 0)
    expect_equal(above$data$boundary_lines$ThurstonianThreshold,
      below$data$boundary_lines$ThurstonianThreshold)
    expect_true(all(above$data$boundary_lines$CumulativeDirection == "at_or_above"))
    label <- if (is.null(contexts)) out$cumulative_boundaries$AboveCategory else
      out$cumulative_boundaries$UpperCategoryLabel
    expect_identical(above$data$boundary_lines$BoundaryLabel, paste0("P(X >= ", label, ") = 0.5"))
    expect_identical(above$data$cumulative_boundaries, out$cumulative_boundaries)
    empty <- plot(out, type = "cumulative", cumulative_direction = "at_or_above",
      boundary_status = "none", draw = FALSE)
    expect_equal(nrow(empty$data$boundary_lines), 0)
  }
})

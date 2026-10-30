# Known perturbations test display arithmetic, not estimator recovery.
comparison_plot <- function(...) .mfrmr_muffle_expected_warnings(
  plot_compare_mfrm(..., draw = FALSE), "^Review-only display:")

test_that("paired Wright plots match identities and retain unavailable differences", {
  a <- make_toy_fit(maxit = 10)
  b <- a
  b$facets$person <- b$facets$person[nrow(b$facets$person):1, ]
  b$facets$person$Estimate <- b$facets$person$Estimate + 0.25
  removed <- b$facets$person$Person[1]
  b$facets$person <- b$facets$person[-1, ]
  missing <- b$facets$person$Person[1]
  b$facets$person$Estimate[1] <- NA_real_
  p <- comparison_plot(a, b)
  d <- p$data$differences
  expect_s3_class(p, "mfrm_plot_data")
  expect_equal(d$Difference[d$Kind == "Person" & is.finite(d$Difference)],
    rep(0.25, sum(d$Kind == "Person" & is.finite(d$Difference))))
  expect_equal(d$Status_Comparison[d$Kind == "Person" & d$Level == removed], "unmatched")
  expect_equal(d$Status_Comparison[d$Kind == "Person" & d$Level == missing], "nonfinite")
  expect_true(all(is.na(d$Difference[d$Kind == "Person" & d$Level %in% c(removed, missing)])))
  expect_true(all(d$Difference[d$Kind != "Person"] == 0))
  expect_equal(p$data$settings$alignment, "none")
  expect_equal(length(p$data$source_plots), 2L)
  expect_identical(a, make_toy_fit(maxit = 10))
  b <- a; b$facets$person$ParameterStatus <- "unbounded_high"
  p <- comparison_plot(a, b)
  expect_true(all(is.na(p$data$differences$Difference[p$data$differences$Kind == "Person"])))
})

test_that("CCC differences use shared engine coordinates with a known direction", {
  a <- make_toy_fit(maxit = 10)
  b <- a; b$steps$Estimate <- b$steps$Estimate + 0.4
  p <- comparison_plot(a, b, type = "ccc", theta_range = c(-2, 2), theta_points = 41)
  d <- p$data$differences
  lowest <- min(as.numeric(d$Category)); highest <- max(as.numeric(d$Category))
  expect_true(all(d$Difference[as.numeric(d$Category) == lowest] > 0))
  expect_true(all(d$Difference[as.numeric(d$Category) == highest] < 0))
  expect_lt(max(abs(tapply(d$Difference, d$Theta, sum))), 1e-12)
  expect_equal(p$data$summary$MaxAbsProbabilityDifference, max(abs(d$Difference)))
  expect_equal(p$data$summary$Difference, d$Difference[which.max(abs(d$Difference))])
  expect_true(all(d$ExpectedScoreDifference < 0))
  same <- comparison_plot(a, a, type = "ccc")
  expect_true(all(same$data$differences$Difference == 0))
  expect_true(all(same$data$summary$MaxAbsProbabilityDifference == 0))
})

test_that("common RSM curves are paired explicitly with selected PCM groups", {
  a <- make_toy_fit(maxit = 10)
  b <- make_toy_fit(model = "PCM", maxit = 10)
  groups <- unique(b$steps$StepFacet)[1:2]
  p <- comparison_plot(a, b, type = "ccc", curve_groups = groups)
  expect_setequal(unique(p$data$probabilities$CurveGroup), groups)
  expect_true(all(p$data$probabilities$SourceCurveGroup[p$data$probabilities$Fit == "Reference"] == "Common"))
  expect_equal(sum(p$data$group_selection$Selected), 2)
  expect_true(any(grepl("Selected 2", p$data$notes$Text)))
  w <- comparison_plot(a, b, curve_groups = groups)
  steps <- w$data$differences[w$data$differences$Kind == "Step", ]
  expect_equal(nrow(steps), 2 * nrow(a$steps))
  expect_setequal(steps$Facet, groups)
  # Hand-check one paired transition against the actual source tables.
  expect_equal(steps$Difference[1], b$steps$Estimate[b$steps$StepFacet == groups[1]][1] - a$steps$Estimate[1])
  b$steps$StepFacet[b$steps$StepFacet == groups[1]] <- "unmatched group"
  b$prep$levels[[b$config$step_facet]] <- unique(b$steps$StepFacet)
  expect_error(comparison_plot(make_toy_fit(model = "PCM", maxit = 10), b, type = "ccc"), "missing in one fit")
  expect_no_error(comparison_plot(make_toy_fit(model = "PCM", maxit = 10), b, type = "ccc", curve_groups = groups[2]))
})

test_that("comparison validates scale, categories, selection and identity", {
  a <- make_toy_fit(maxit = 10)
  for (field in c("noncenter_facet", "positive_facets", "method")) {
    b <- a; b$config[[field]] <- "different"
    expect_error(comparison_plot(a, b), "Comparison basis differs")
  }
  b <- a; b$prep$rating_max <- b$prep$rating_max + 1
  expect_error(comparison_plot(a, b), "Comparison basis differs")
  b <- a; b$steps$Estimate[1] <- NA_real_
  expect_error(comparison_plot(a, b, type = "ccc"), "complete finite step")
  b <- a; b$facets$others <- rbind(b$facets$others, b$facets$others[1, ])
  expect_error(suppressWarnings(comparison_plot(a, b)), "unique non-missing")
  expect_error(comparison_plot(a, a, labels = c("same", "same")), "distinct")
  expect_error(comparison_plot(a, a, curve_groups = "absent"), "known")
  expect_error(comparison_plot(a, a, theta_points = 1), "theta_points")
  expect_error(comparison_plot(a, a, show_notes = NA), "TRUE or FALSE")
  expect_error(comparison_plot(list(), a), "mfrm_fit")
})

test_that("original scores, eleven categories and GPCM slopes survive comparison", {
  a <- make_toy_fit(maxit = 10)
  # Controlled curve fixture: 0:10 internally, 5:15 on the original scale.
  a$config$n_cat <- 11L; a$prep$rating_min <- 0; a$prep$rating_max <- 10
  a$prep$score_map <- data.frame(InternalScore = 0:10, OriginalScore = 5:15)
  a$steps <- data.frame(Step = paste0("Step_", 1:10), Estimate = seq(-2, 2, length.out = 10))
  b <- a
  b$config$model <- "GPCM"; b$config$step_facet <- "Criterion"
  b$config$slope_facet <- "Criterion"
  b$prep$levels$Criterion <- "C1"
  b$steps$StepFacet <- "C1"
  b$slopes <- data.frame(SlopeFacet = "C1", Estimate = 1.4)
  p <- comparison_plot(a, b, type = "ccc", theta_points = 21)
  expect_equal(p$data$category_labels$OriginalScore, as.character(5:15))
  expect_setequal(p$data$probabilities$OriginalCategory, as.character(5:15))
  expect_equal(as.numeric(p$data$summary$OriginalCategory), as.numeric(p$data$summary$Category) + 5)
  t <- p$data$probabilities[p$data$probabilities$Fit == "Comparison", ]
  expect_true(all(t$Slope == 1.4))
  theta <- t$Theta[t$Category == "0"]
  logits <- 1.4 * (outer(theta, 0:10) - matrix(c(0, cumsum(b$steps$Estimate)), length(theta), 11, byrow = TRUE))
  probs <- exp(logits - apply(logits, 1, max)); probs <- probs / rowSums(probs)
  expect_equal(t$Probability, as.vector(probs), tolerance = 1e-12)
  w <- comparison_plot(a, b)
  expect_true("5 -> 6" %in% w$data$locations$Level)
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  expect_warning(g <- as_ggplot(p), "Many comparison panels")
  built <- ggplot2::ggplot_build(g)
  expect_equal(nrow(built$layout$layout), 11)
  expect_setequal(as.character(built$layout$layout$Category), as.character(5:15))
})

test_that("comparison renderers retain notes, layout and monochrome encodings", {
  skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  a <- make_toy_fit(maxit = 10)
  b <- a; b$facets$person <- b$facets$person[1, ]
  grDevices::pdf(NULL, width = 10, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2))
  old <- graphics::par(c("mfrow", "mar", "bg"))
  for (type in c("wright", "ccc")) for (view in c("comparison", "difference")) {
    p <- comparison_plot(a, b, type = type, view = view, preset = "monochrome",
      show_title = FALSE, show_notes = FALSE)
    g <- as_ggplot(p)
    expect_null(g$labels$title)
    expect_null(g$labels$subtitle)
    expect_identical(utils::head(plot_data(g, "notes"), nrow(p$data$notes)), p$data$notes)
    expect_true(nrow(p$data$fit_readiness) > 0)
    expect_no_warning(ggplot2::ggplot_build(g))
    expect_no_warning(print(g))
    actual <- .mfrmr_muffle_expected_warnings(plot_compare_mfrm(a, b, type = type, view = view,
      preset = "monochrome", show_title = FALSE, show_notes = FALSE), "^Review-only display:")
    expect_identical(actual, p)
    expect_equal(graphics::par(c("mfrow", "mar", "bg")), old)
    if (type == "ccc" && view == "comparison") {
      expect_equal(nrow(ggplot2::ggplot_build(g)$layout$layout), a$config$n_cat)
      expect_equal(unname(ggplot2::ggplot_build(g)$plot$scales$get_scales("linetype")$map(c("Reference", "Comparison"))), c("solid", "dashed"))
    }
  }
  p <- comparison_plot(a, b, type = "ccc", preset = "monochrome", panel = "group")
  expect_warning(as_ggplot(p), "indistinguishable")
  a$facets$person$Estimate <- a$facets$others$Estimate <- a$steps$Estimate <- 0
  p <- comparison_plot(a, a)
  built <- ggplot2::ggplot_build(as_ggplot(p))
  text_layer <- which(vapply(built$plot$layers, function(layer) inherits(layer$geom, "GeomLabel"), logical(1)))
  expect_equal(length(unique(built$data[[text_layer]]$y)), nrow(a$steps))
})

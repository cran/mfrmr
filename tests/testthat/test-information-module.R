test_that("compute_information returns stable precision-curve structures", {
  toy <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = c("Content", "Organization", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
    2 * as.integer(factor(toy$Rater)) +
    as.integer(factor(toy$Criterion))
  ) %% 3

  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", model = "RSM", maxit = 25)
  info <- compute_information(fit, theta_points = 11)

  expect_s3_class(info, "mfrm_information")
  expect_equal(nrow(info$tif), 11)
  expect_true(all(c("Theta", "Information", "SE") %in% names(info$tif)))
  expect_true(all(c("Theta", "ConditionalSEM", "Information") %in% names(info$conditional_sem)))
  expect_true(all(c("PlotType", "Metric", "Theta", "Value", "ValueName") %in%
                    names(info$information_long)))
  expect_true(all(c("ThetaAtMaxInformation", "MaxInformation", "MinConditionalSEM") %in%
                    names(info$summary)))
  expect_equal(nrow(info$iif), nrow(fit$facets$others) * 11)
  expect_true(all(c("Theta", "Facet", "Level", "Information", "Exposure") %in% names(info$iif)))
  expect_true(all(is.finite(info$tif$Information)))
  expect_true(all(info$tif$Information >= 0))
})

test_that("compute_information reflects realized observation exposure", {
  toy <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = c("Content", "Organization", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
      2 * as.integer(factor(toy$Rater)) +
      as.integer(factor(toy$Criterion))
  ) %% 3

  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", model = "RSM", maxit = 25)
  fit_dup <- fit
  fit_dup$prep$data <- rbind(fit$prep$data, fit$prep$data)

  info <- compute_information(fit, theta_points = 11)
  info_dup <- compute_information(fit_dup, theta_points = 11)

  expect_equal(info_dup$tif$Information, 2 * info$tif$Information)
})

test_that("plot_information returns plot data for tif and iif views", {
  toy <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = c("Content", "Organization", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
    2 * as.integer(factor(toy$Rater)) +
    as.integer(factor(toy$Criterion))
  ) %% 3

  fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", model = "RSM", maxit = 25)
  info <- compute_information(fit, theta_points = 11)

  tif_plot <- plot_information(info, type = "tif", draw = FALSE)
  sem_plot <- plot_information(info, type = "sem", draw = FALSE)
  both_plot <- plot_information(info, type = "both", draw = FALSE)
  iif_plot <- plot_information(info, type = "iif", facet = "Rater",
                               draw = FALSE)

  expect_s3_class(tif_plot, "mfrm_plot_data")
  expect_s3_class(sem_plot, "mfrm_plot_data")
  expect_s3_class(both_plot, "mfrm_plot_data")
  expect_s3_class(iif_plot, "mfrm_plot_data")
  expect_true(is.data.frame(tif_plot$data$plot) || tibble::is_tibble(tif_plot$data$plot))
  expect_equal(nrow(tif_plot$data$plot), 11)
  expect_true(all(c("plot_long", "conditional_sem", "summary", "settings") %in%
                    names(sem_plot$data)))
  expect_true(any(sem_plot$data$plot_long$ValueName == "ConditionalSEM" &
                    sem_plot$data$plot_long$DisplayedByDefault))
  expect_true(all(c("Information", "ConditionalSEM") %in% both_plot$data$series))
  expect_true(is.data.frame(iif_plot$data$plot) || tibble::is_tibble(iif_plot$data$plot))
  expect_true(all(iif_plot$data$plot$Facet == "Rater"))
  expect_equal(length(unique(iif_plot$data$plot$Theta)), 11)
  expect_true(all(iif_plot$data$plot_long$ValueName == "InformationContribution"))
})

test_that("compute_information supports PCM fits with custom step facet names", {
  toy <- load_mfrmr_data("example_core")
  names(toy)[names(toy) == "Rater"] <- "Judge"
  names(toy)[names(toy) == "Criterion"] <- "Task"
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Judge", "Task"), "Score",
             method = "JML", model = "PCM", step_facet = "Task", maxit = 25)
  )

  info <- compute_information(fit, theta_points = 11)
  task_plot <- plot_information(info, type = "iif", facet = "Task", draw = FALSE)
  surface_plot <- plot(fit, type = "ccc_surface", draw = FALSE, theta_points = 55)

  expect_s3_class(info, "mfrm_information")
  expect_equal(nrow(info$tif), 11)
  expect_true(all(c("Judge", "Task") %in% unique(info$iif$Facet)))
  expect_true(all(is.finite(info$tif$Information)))
  expect_true(all(info$tif$Information >= 0))
  expect_s3_class(task_plot, "mfrm_plot_data")
  expect_true(all(task_plot$data$plot$Facet == "Task"))
  expect_equal(length(unique(task_plot$data$plot$Theta)), 11)
  expect_s3_class(surface_plot, "mfrm_plot_data")
  expect_true(all(c("SurfaceX", "SurfaceY", "SurfaceZ") %in% names(surface_plot$data$surface)))
  expect_true(all(sort(unique(surface_plot$data$surface$Category)) == sort(as.character(fit$prep$rating_min:fit$prep$rating_max))))
})

test_that("compute_information supports bounded GPCM fits", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(
    fit_mfrm(
      toy, "Person", c("Rater", "Criterion"), "Score",
      method = "JML", model = "GPCM", step_facet = "Criterion", maxit = 20
    )
  )

  info <- compute_information(fit, theta_points = 11)
  criterion_plot <- plot_information(info, type = "iif", facet = "Criterion", draw = FALSE)
  surface_plot <- plot(fit, type = "ccc_surface", draw = FALSE, theta_points = 55)
  curves <- category_curves_report(fit, theta_points = 55, digits = 8)
  curve_cumulative_plot <- plot(curves, type = "cumulative", draw = FALSE)
  curve_info_plot <- plot(curves, type = "information", draw = FALSE)
  curve_cat_info_plot <- plot(curves, type = "category_information", draw = FALSE)

  expect_s3_class(info, "mfrm_information")
  expect_equal(nrow(info$tif), 11)
  expect_true(all(is.finite(info$tif$Information)))
  expect_true(all(info$tif$Information >= 0))
  expect_true(all(c("Rater", "Criterion") %in% unique(info$iif$Facet)))
  expect_s3_class(criterion_plot, "mfrm_plot_data")
  expect_true(all(criterion_plot$data$plot$Facet == "Criterion"))
  expect_equal(length(unique(criterion_plot$data$plot$Theta)), 11)
  expect_s3_class(surface_plot, "mfrm_plot_data")
  expect_true(all(c("SurfaceX", "SurfaceY", "SurfaceZ") %in% names(surface_plot$data$surface)))
  expect_true(all(sort(unique(surface_plot$data$surface$Category)) == sort(as.character(fit$prep$rating_min:fit$prep$rating_max))))
  expect_true(all(c("ScoreVariance", "Information", "Slope") %in%
                    names(curves$expected_ogive)))
  expect_true(all(c("CumulativeProbability", "Direction", "BoundaryCategory") %in%
                    names(curves$cumulative_probabilities)))
  expect_true(all(c("ThurstonianThreshold", "InThetaRange", "CrossingCount", "BoundaryStatus") %in%
                    names(curves$cumulative_boundaries)))
  cumulative_max <- curves$cumulative_probabilities[
    curves$cumulative_probabilities$Direction == "at_or_below" &
      curves$cumulative_probabilities$BoundaryOrder == max(curves$cumulative_probabilities$BoundaryOrder, na.rm = TRUE),
    "CumulativeProbability",
    drop = TRUE
  ]
  expect_equal(cumulative_max, rep(1, length(cumulative_max)), tolerance = 1e-6)
  cumulative_min <- curves$cumulative_probabilities[
    curves$cumulative_probabilities$Direction == "at_or_above" &
      curves$cumulative_probabilities$BoundaryOrder == min(curves$cumulative_probabilities$BoundaryOrder, na.rm = TRUE),
    "CumulativeProbability",
    drop = TRUE
  ]
  expect_equal(cumulative_min, rep(1, length(cumulative_min)), tolerance = 1e-6)
  below <- curves$cumulative_probabilities[
    curves$cumulative_probabilities$Direction == "at_or_below",
    c("CurveGroup", "Theta", "BoundaryOrder", "CumulativeProbability"),
    drop = FALSE
  ]
  below <- below[order(below$CurveGroup, below$Theta, below$BoundaryOrder), , drop = FALSE]
  below_groups <- split(seq_len(nrow(below)), list(below$CurveGroup, below$Theta), drop = TRUE)
  below_diff <- unlist(lapply(below_groups, function(i) diff(below$CumulativeProbability[i])), use.names = FALSE)
  expect_true(all(below_diff >= -1e-6))
  above <- curves$cumulative_probabilities[
    curves$cumulative_probabilities$Direction == "at_or_above",
    c("CurveGroup", "Theta", "BoundaryOrder", "CumulativeProbability"),
    drop = FALSE
  ]
  above <- above[order(above$CurveGroup, above$Theta, above$BoundaryOrder), , drop = FALSE]
  above_groups <- split(seq_len(nrow(above)), list(above$CurveGroup, above$Theta), drop = TRUE)
  above_diff <- unlist(lapply(above_groups, function(i) diff(above$CumulativeProbability[i])), use.names = FALSE)
  expect_true(all(above_diff <= 1e-6))
  below_next <- below
  below_next$NextBoundaryOrder <- below_next$BoundaryOrder + 1L
  names(below_next)[names(below_next) == "CumulativeProbability"] <- "BelowProbability"
  names(above)[names(above) == "CumulativeProbability"] <- "AboveProbability"
  cumulative_identity <- merge(
    below_next,
    above,
    by.x = c("CurveGroup", "Theta", "NextBoundaryOrder"),
    by.y = c("CurveGroup", "Theta", "BoundaryOrder"),
    all = FALSE
  )
  expect_equal(
    cumulative_identity$BelowProbability + cumulative_identity$AboveProbability,
    rep(1, nrow(cumulative_identity)),
    tolerance = 1e-6
  )
  expect_true(all(curves$expected_ogive$Slope > 0))
  expect_equal(
    curves$expected_ogive$Information,
    curves$expected_ogive$ScoreVariance * curves$expected_ogive$Slope^2,
    tolerance = 1e-6
  )
  category_info_sum <- stats::aggregate(
    CategoryInformation ~ CurveGroup + Theta,
    data = curves$category_information,
    FUN = sum
  )
  category_info_sum <- merge(
    category_info_sum,
    curves$expected_ogive[, c("CurveGroup", "Theta", "Information")],
    by = c("CurveGroup", "Theta"),
    all.x = TRUE
  )
  expect_equal(
    category_info_sum$CategoryInformation,
    category_info_sum$Information,
    tolerance = 1e-5
  )
  expect_s3_class(curve_info_plot, "mfrm_plot_data")
  expect_s3_class(curve_cumulative_plot, "mfrm_plot_data")
  expect_identical(curve_cumulative_plot$data$plot, "cumulative")
  expect_true(all(c("cumulative_probabilities", "cumulative_boundaries") %in%
                    names(curve_cumulative_plot$data)))
  expect_identical(curve_info_plot$data$plot, "information")
  expect_true(all(c("Information", "Slope") %in%
                    names(curve_info_plot$data$expected_ogive)))
  expect_s3_class(curve_cat_info_plot, "mfrm_plot_data")
  expect_identical(curve_cat_info_plot$data$plot, "category_information")
  expect_true(all(c("CategoryInformation", "CategoryInformationShare") %in%
                    names(curve_cat_info_plot$data$category_information)))
})

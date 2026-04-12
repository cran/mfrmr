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
  iif_plot <- plot_information(info, type = "iif", facet = "Rater",
                               draw = FALSE)

  expect_s3_class(tif_plot, "mfrm_plot_data")
  expect_s3_class(iif_plot, "mfrm_plot_data")
  expect_true(is.data.frame(tif_plot$data$plot) || tibble::is_tibble(tif_plot$data$plot))
  expect_equal(nrow(tif_plot$data$plot), 11)
  expect_true(is.data.frame(iif_plot$data$plot) || tibble::is_tibble(iif_plot$data$plot))
  expect_true(all(iif_plot$data$plot$Facet == "Rater"))
  expect_equal(length(unique(iif_plot$data$plot$Theta)), 11)
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
})

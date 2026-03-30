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

  expect_true(is.data.frame(tif_plot) || tibble::is_tibble(tif_plot))
  expect_equal(nrow(tif_plot), 11)
  expect_true(is.data.frame(iif_plot) || tibble::is_tibble(iif_plot))
  expect_true(all(iif_plot$Facet == "Rater"))
  expect_equal(length(unique(iif_plot$Theta)), 11)
})

test_that("compute_information rejects PCM fits until step-facet information is implemented", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", model = "PCM", maxit = 25)
  )

  expect_error(
    compute_information(fit),
    "currently supports only `model = \"RSM\"` fits",
    fixed = TRUE
  )
})

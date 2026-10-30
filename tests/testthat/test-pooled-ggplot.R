pooled_plot_fixture <- function() {
  estimate <- c(-.4, 0, .3)
  se <- c(.2, 0, .4)
  df <- c(3, NA, Inf)
  margin <- stats::qt(.95, df) * se
  structure(list(table = data.frame(Target = c("R02 minus R01", "Constrained mean", "R03 minus R01"),
      Estimate = estimate, SE = se, Lower = estimate - margin, Upper = estimate + margin,
      DF = df, Status = c("pooled", "fixed", "pooled")),
    contrasts = matrix(c(-1, 1, 0, 1/3, 1/3, 1/3, -1, 0, 1), 3, byrow = TRUE),
    settings = list(facet = "Rater", imputations = 2L, ci_level = .9, df_complete = Inf,
      intervals = "pointwise multiple-imputation t intervals"),
    cautions = "Imputation(s) 2: weak complete-data information.",
    information_review = list(list(status = "ready"), list(status = "weak"))), class = "mfrm_pooled")
}

test_that("pooled conversion preserves saved t intervals, fixed targets and interpretation", {
  skip_if_not_installed("ggplot2")
  x <- pooled_plot_fixture()
  p <- plot(x, preset = "monochrome", draw = FALSE)
  g <- as_ggplot(p, component = "table")
  b <- ggplot2::ggplot_build(g)
  expect_equal(b$data[[1]]$xintercept, 0)
  expect_equal(b$data[[2]]$x, x$table$Lower[c(1, 3)])
  expect_equal(b$data[[2]]$xend, x$table$Upper[c(1, 3)])
  expect_equal(b$data[[3]]$x, x$table$Estimate)
  expect_equal(b$data[[3]]$y, 3:1)
  expect_equal(b$data[[3]]$shape, c(16, 5, 16))
  expect_true(abs(b$data[[2]]$x[1] - (x$table$Estimate[1] - stats::qnorm(.95) * x$table$SE[1])) > .1)
  expect_identical(b$layout$panel_params[[1]]$y$get_labels()[
    match(3:1, b$layout$panel_params[[1]]$y$get_breaks())], x$table$Target)
  expect_identical(plot_data(g), p$data)
  expect_identical(plot_data(g)$table$DF, c(3, NA, Inf))
  expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "90% pointwise MI", fixed = TRUE)
  expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "Weak information", fixed = TRUE)
  expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "fixed targets", fixed = TRUE)
  expect_silent(ggplot2::ggplotGrob(g))
})

test_that("unbounded and missing intervals stay distinct and do not become finite bounds", {
  skip_if_not_installed("ggplot2")
  x <- pooled_plot_fixture()
  x$table$Lower[1] <- -Inf; x$table$Upper[1] <- Inf
  x$table$Lower[3] <- NA_real_
  x$table$InferenceCaution <- c("Extreme degrees of freedom", "Fixed by constraint", "Unavailable endpoint")
  p <- plot(x, draw = FALSE)
  g <- as_ggplot(p)
  b <- ggplot2::ggplot_build(g)
  expect_equal(nrow(b$data[[2]]), 1)
  expect_equal(b$data[[3]]$shape, c(1, 5, 4))
  expect_length(g$layers, 5)
  expect_true(b$data[[4]]$xend < b$data[[4]]$x)
  expect_true(b$data[[5]]$xend > b$data[[5]]$x)
  expect_false(is.null(g$layers[[4]]$geom_params$arrow))
  expect_identical(plot_data(g)$table, x$table)
  expect_identical(plot_data(g)$table$Lower[1], -Inf)
  expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "infinite endpoints", fixed = TRUE)
  expect_silent(ggplot2::ggplotGrob(g))
  fixed <- x; fixed$table <- x$table[2, , drop = FALSE]
  expect_silent(ggplot2::ggplotGrob(as_ggplot(fixed)))
  unavailable <- x; unavailable$table <- x$table[3, , drop = FALSE]
  expect_silent(ggplot2::ggplotGrob(as_ggplot(unavailable)))
})

test_that("pooled replay and display controls do not refit, repool or recalculate intervals", {
  skip_if_not_installed("ggplot2")
  x <- pooled_plot_fixture()
  p <- plot(x, title = NULL, subtitle = NULL, preset = "monochrome", draw = FALSE)
  direct <- as_ggplot(x, title = NULL, subtitle = NULL, preset = "monochrome")
  local_mocked_bindings(pool_mfrm_imputed = function(...) stop("Unexpected pooling"),
    fit_mfrm_imputed = function(...) stop("Unexpected fitting"),
    compute_mml_parameter_covariance = function(...) stop("Unexpected covariance"))
  local_mocked_bindings(qt = function(...) stop("Unexpected interval calculation"), .package = "mfrmr")
  withr::local_options(mfrmr.plot_preset = "invalid")
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(p, file)
  device <- grDevices::dev.cur()
  theme <- ggplot2::theme_get()
  g <- as_ggplot(readRDS(file))
  expect_equal(ggplot2::ggplot_build(g)$data, ggplot2::ggplot_build(direct)$data)
  expect_null(g$labels$title); expect_null(g$labels$subtitle)
  plain <- g + ggplot2::labs(caption = NULL)
  expect_identical(plot_data(plain), p$data)
  expect_identical(plot_data(plain)$cautions, x$cautions)
  expect_identical(plot_data(plain)$information_review, x$information_review)
  expect_identical(grDevices::dev.cur(), device)
  expect_identical(ggplot2::theme_get(), theme)
  expect_error(as_ggplot(p, component = "contrasts"), "plot_data")
  expect_error(as_ggplot(p, level = .95), "empty")
  expect_error(as_ggplot(p, type = "forest"), "already selects a view")
  bad <- p; bad$data$table$DF <- NULL
  expect_error(as_ggplot(bad), "saved interval fields")
  bad <- p; bad$data$table$Lower[1] <- 4
  expect_error(as_ggplot(bad), "ordered")
  bad <- p; bad$data$table$Estimate[1] <- NA_real_
  expect_error(as_ggplot(bad), "finite estimates")
})

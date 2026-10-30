co_membership_plot_fixture <- function() {
  ids <- c("001", "1", "NA", "excluded")
  memberships <- list(c(1, 1, 2, NA), c(1, 1, 2, NA),
    c(1, 2, 2, NA), c(1, 2, 3, NA))
  mat <- Reduce(`+`, lapply(memberships, function(g) outer(g, g, "=="))) / 4
  dimnames(mat) <- list(ids, ids)
  structure(list(co_membership = mat, settings = list(imputations = 4L)),
    class = "mfrm_imputed_clusters")
}

test_that("co-membership conversion preserves fractions, order and unavailable cells", {
  skip_if_not_installed("ggplot2")
  x <- co_membership_plot_fixture()
  ids <- c("excluded", "NA", "001", "1")
  p <- plot(x, ids = ids, preset = "monochrome", draw = FALSE)
  g <- as_ggplot(p, component = "matrix")
  b <- ggplot2::ggplot_build(g)
  expect_equal(nrow(b$data[[1]]), 16)
  # Interpret plotted coordinates independently in terms of requested IDs.
  for (i in seq_len(nrow(g$data))) {
    row <- ids[length(ids) + 1 - g$data$Y[i]]
    col <- ids[g$data$Column[i]]
    expect_identical(g$data$Fraction[i], x$co_membership[row, col])
  }
  expect_identical(b$layout$panel_params[[1]]$x$get_labels(), ids)
  expect_identical(b$layout$panel_params[[1]]$y$get_labels(), rev(ids))
  expect_equal(b$plot$scales$get_scales("fill")$limits, c(0, 1))
  expect_equal(b$plot$scales$get_scales("fill")$breaks, seq(0, 1, .25))
  expect_equal(nrow(b$data[[2]]), sum(is.na(p$data$matrix)))
  expect_true(all(b$data[[2]]$shape == 4))
  missing_xy <- g$data[is.na(g$data$Fraction), c("Column", "Y")]
  expect_equal(b$data[[2]]$x, missing_xy$Column)
  expect_equal(b$data[[2]]$y, missing_xy$Y)
  key <- p$data$legend
  fill <- b$plot$scales$get_scales("fill")
  # Both grey and a non-colour symbol distinguish unavailable cells from zero.
  expect_false(identical(fill$map(0), fill$map(NA_real_)))
  expect_identical(fill$map(NA_real_), key$value[key$role == "missing"])
  expect_identical(plot_data(g), p$data)
  expect_identical(plot_data(g)$imputations, 4L)
  expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "all supplied imputations", fixed = TRUE)
  expect_silent(ggplot2::ggplotGrob(g))
  expect_equal(ggplot2::ggplot_build(as_ggplot(p))$data, b$data)
})

test_that("selection and hidden labels do not renormalize or remove cells", {
  skip_if_not_installed("ggplot2")
  x <- co_membership_plot_fixture()
  g <- as_ggplot(x, ids = c("1", "001"), labels = FALSE)
  expect_equal(sort(g$data$Fraction), c(.5, .5, 1, 1))
  expect_null(g$scales$get_scales("x")$breaks)
  expect_null(g$scales$get_scales("y")$breaks)
  expect_length(g$layers, 1)
  expect_identical(plot_data(g)$excluded_ids, "excluded")
  for (id in c("excluded", "001")) {
    g <- as_ggplot(x, ids = id)
    expect_equal(nrow(g$data), 1)
    expect_identical(g$data$Fraction, x$co_membership[id, id])
    expect_equal(g$scales$get_scales("fill")$limits, c(0, 1))
    expect_silent(ggplot2::ggplotGrob(g))
  }
  ids <- paste0("R", 1:51)
  x$co_membership <- matrix(1, 51, 51, dimnames = list(ids, ids))
  p <- plot(x, draw = FALSE)
  expect_false(p$data$labels)
  g <- as_ggplot(p)
  expect_equal(nrow(g$data), 51^2)
  expect_identical(plot_data(g)$ids, ids)
})

test_that("saved co-membership conversion uses saved colours and avoids new analyses", {
  skip_if_not_installed("ggplot2")
  x <- co_membership_plot_fixture()
  p <- plot(x, preset = "monochrome", draw = FALSE)
  direct <- as_ggplot(x, preset = "monochrome")
  local_mocked_bindings(mfrm_cluster_imputed = function(...) stop("Unexpected imputation"),
    mfrm_cluster = function(...) stop("Unexpected clustering"))
  withr::local_options(mfrmr.plot_preset = "invalid")
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(p, file)
  device <- grDevices::dev.cur()
  theme <- ggplot2::theme_get()
  g <- as_ggplot(readRDS(file))
  expect_equal(ggplot2::ggplot_build(g)$data, ggplot2::ggplot_build(direct)$data)
  expect_identical(ggplot2::theme_get(), theme)
  expect_identical(grDevices::dev.cur(), device)
  plain <- g + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
  expect_null(plain$labels$caption)
  expect_identical(plot_data(plain), p$data)
  expect_error(as_ggplot(p, component = "legend"), "plot_data")
  expect_error(as_ggplot(p, ids = "001"), "empty")
  expect_error(as_ggplot(p, type = "heatmap"), "already selects a view")
  bad <- p; bad$data$matrix[1, 2] <- 2
  expect_error(as_ggplot(bad), "zero-to-one matrix")
  bad <- p; bad$data$ids <- rev(bad$data$ids)
  expect_error(as_ggplot(bad), "matching saved IDs")
  bad <- p; bad$data$imputations <- NULL
  expect_error(as_ggplot(bad), "number of imputations")
  bad <- p; bad$data$legend <- NULL
  expect_error(as_ggplot(bad), "saved proportion")
})

cluster_summary_plot_fixture <- function() {
  data <- data.frame(ID = paste0("R", 1:8), Years = c(1, 2, 3, 12, 13, 18, 20, NA),
    Training = ordered(c(rep("Basic", 3), rep("Advanced", 5)),
      levels = c("Basic", "Intermediate", "Advanced")))
  mfrm_cluster(mfrm_features(data, "ID", c("Years", "Training")), 2, missing = "omit")
}

test_that("silhouette conversion preserves negative widths, order, references and IDs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  p <- plot(cluster_summary_plot_fixture(), preset = "monochrome", draw = FALSE)
  # Negative widths are valid, including in a saved partition not sorted again.
  p$data$table$Silhouette[2] <- -.4
  p$data$reference_lines$value <- .123
  g <- as_ggplot(p)
  b <- ggplot2::ggplot_build(g)
  tab <- p$data$table
  expect_equal(b$data[[1]]$xmin, pmin(0, tab$Silhouette))
  expect_equal(b$data[[1]]$xmax, pmax(0, tab$Silhouette))
  expect_equal((b$data[[1]]$ymin + b$data[[1]]$ymax)/2, nrow(tab) + 1 - seq_len(nrow(tab)))
  expect_equal(b$data[[2]]$xintercept, 0)
  expect_equal(b$data[[3]]$xintercept, .123)
  expect_equal(b$plot$scales$get_scales("x")$limits, c(-1, 1))
  expect_identical(b$layout$panel_params[[1]]$y.sec$get_labels(), rev(tab$ID))
  expect_identical(plot_data(g), p$data)
  p$data$labels <- FALSE
  no_labels <- as_ggplot(p)
  expect_equal(nrow(ggplot2::ggplot_build(no_labels)$data[[1]]), nrow(tab))
  expect_null(no_labels$scales$get_scales("y")$secondary.axis$breaks)
  expect_silent(ggplot2::ggplotGrob(no_labels))
})

test_that("numeric profiles preserve original values, counts and separate mean/median symbols", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  x <- cluster_summary_plot_fixture()
  p <- plot(x, type = "profile", feature = "Years", preset = "monochrome", draw = FALSE)
  g <- as_ggplot(p, component = "table")
  b <- ggplot2::ggplot_build(g)
  expect_equal(b$data[[1]]$x, c(p$data$table$Mean, p$data$table$Median))
  expect_equal(b$data[[1]]$shape, c(16, 16, 17, 17))
  expect_equal(diff(b$data[[1]]$y[c(3, 1)]), .12)
  expect_identical(g$labels$x, "Years (original units)")
  expected <- rev(sprintf("Group %d (n = %d)", p$data$table$Cluster, p$data$table$N))
  expect_identical(b$layout$panel_params[[1]]$y$get_labels(), expected)
  expect_identical(plot_data(g), p$data)
  expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "without confidence intervals", fixed = TRUE)
  expect_silent(ggplot2::ggplotGrob(g))
})

test_that("categorical profiles preserve unused levels, original order and fixed proportions", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  p <- plot(cluster_summary_plot_fixture(), type = "profile", feature = "Training", draw = FALSE)
  g <- as_ggplot(p)
  b <- ggplot2::ggplot_build(g)
  expect_identical(g$data$Proportion, as.vector(p$data$matrix))
  expect_identical(colnames(p$data$matrix), c("Basic", "Intermediate", "Advanced"))
  expect_true(all(g$data$Proportion[g$data$Column == 2] == 0))
  expect_identical(b$layout$panel_params[[1]]$x$get_labels(), colnames(p$data$matrix))
  expect_identical(b$layout$panel_params[[1]]$y$get_labels(), rev(rownames(p$data$matrix)))
  expect_equal(g$scales$get_scales("fill")$limits, c(0, 1))
  expect_identical(plot_data(g), p$data)
  for (component in c("table", "matrix")) {
    expect_equal(ggplot2::ggplot_build(as_ggplot(p, component = component))$data, b$data)
  }
  expect_silent(ggplot2::ggplotGrob(g))
})

test_that("partition renderers replay saved data and reject generic fallbacks", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  x <- cluster_summary_plot_fixture()
  cases <- list(list(), list(type = "profile", feature = "Years"), list(type = "profile", feature = "Training"))
  plots <- lapply(cases, function(args) do.call(plot, c(list(x, draw = FALSE, preset = "monochrome"), args)))
  direct <- lapply(cases, function(args) do.call(as_ggplot, c(list(x, preset = "monochrome"), args)))
  local_mocked_bindings(mfrm_cluster = function(...) stop("Unexpected clustering"),
    cluster_external_features = function(...) stop("Unexpected clustering"))
  withr::local_options(mfrmr.plot_preset = "invalid")
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(plots, file)
  saved <- readRDS(file)
  theme <- ggplot2::theme_get()
  device <- grDevices::dev.cur()
  for (i in seq_along(saved)) {
    p <- saved[[i]]
    g <- as_ggplot(p)
    expect_equal(ggplot2::ggplot_build(g)$data, ggplot2::ggplot_build(direct[[i]])$data)
    expect_identical(plot_data(g)$excluded_ids, "R8")
    plain <- g + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
    expect_null(plain$labels$title)
    expect_identical(plot_data(plain), p$data)
    expect_error(as_ggplot(p, component = "legend"), "plot_data")
    expect_error(as_ggplot(p, labels = FALSE), "empty")
    expect_error(as_ggplot(p, type = "profile"), "already selects a view")
  }
  expect_identical(ggplot2::theme_get(), theme)
  expect_identical(grDevices::dev.cur(), device)
  bad <- saved[[1]]; bad$data$table$Silhouette[1] <- NA_real_
  expect_error(as_ggplot(bad), "finite")
  bad <- saved[[1]]; bad$data$reference_lines <- NULL
  expect_error(as_ggplot(bad), "overall-mean reference")
  bad <- saved[[2]]; bad$data$table$Mean <- NULL
  expect_error(as_ggplot(bad), "saved summary table")
  bad <- saved[[3]]; bad$data$matrix[1, 1] <- 2
  expect_error(as_ggplot(bad), "zero-to-one")
})

test_that("k-means and hierarchy share profile conversion without forcing silhouettes", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  f <- mfrm_features(data.frame(ID = letters[1:6], Years = c(1, 2, 4, 10, 14, 16)), "ID", "Years")
  k <- mfrm_cluster_kmeans(f, 2, seed = 5, silhouette = FALSE)
  expect_error(as_ggplot(k), "Silhouettes are unavailable")
  expect_silent(ggplot2::ggplotGrob(as_ggplot(k, type = "profile", feature = "Years")))
  h <- mfrm_cluster_hierarchical(f, 2)
  expect_silent(ggplot2::ggplotGrob(as_ggplot(h, type = "silhouette")))
  expect_silent(ggplot2::ggplotGrob(as_ggplot(h, type = "profile", feature = "Years")))
})

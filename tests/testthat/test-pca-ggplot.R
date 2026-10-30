pca_plot_fixture <- function() {
  features <- mfrm_features(data.frame(ID = c(paste0("R", 1:8), "excluded"),
    Years = c(1, 2, 4, 6, 8, 10, 12, 15, NA),
    Hours = c(8, 16, 12, 24, 16, 32, 28, 40, 20),
    Ratings = c(30, 80, 60, 150, 120, 200, 170, 250, 100)),
    "ID", c("Years", "Hours", "Ratings"))
  pca <- mfrm_pca(features, components = 2, missing = "omit")
  list(pca = pca, groups = mfrm_cluster_kmeans(pca, 3, seed = 42, silhouette = FALSE))
}

test_that("PCA conversions keep selected axes, IDs, group encodings and saved settings", {
  skip_if_not_installed("ggplot2")
  f <- pca_plot_fixture()
  # Reorder memberships to expose any positional rather than ID-based matching.
  f$groups$membership <- f$groups$membership[9:1, ]
  p <- plot(f$pca, type = "scores", components = c(2, 1), groups = f$groups,
    preset = "monochrome", draw = FALSE)
  withr::local_options(mfrmr.plot_preset = "invalid")
  local_mocked_bindings(mfrm_pca = function(...) stop("Unexpected PCA"),
    mfrm_cluster_kmeans = function(...) stop("Unexpected clustering"))
  device <- grDevices::dev.cur()
  theme <- ggplot2::theme_get()
  g <- as_ggplot(p, component = "table")
  built <- ggplot2::ggplot_build(g)
  tab <- p$data$table
  enc <- p$data$encoding
  expect_equal(built$data[[1]]$x, tab$PC2)
  expect_equal(built$data[[1]]$y, tab$PC1)
  expect_equal(built$data[[1]]$colour, enc$Colour[match(tab$Cluster, enc$Cluster)])
  expect_equal(built$data[[1]]$shape, enc$Shape[match(tab$Cluster, enc$Cluster)])
  expect_identical(built$data[[2]]$label, tab$ID)
  expect_identical(tab$Cluster, f$groups$membership$Cluster[match(tab$ID, f$groups$membership$ID)])
  expect_equal(g$coordinates$ratio, 1)
  expect_equal(diff(g$coordinates$limits$x), diff(g$coordinates$limits$y))
  expect_identical(g$labels$x, p$data$xlab)
  expect_identical(g$labels$y, p$data$ylab)
  expect_identical(plot_data(g), p$data)
  expect_identical(plot_data(g)$excluded_ids, "excluded")
  expect_match(g$labels$alt, "not latent ability", fixed = TRUE)
  expect_equal(ggplot2::ggplot_build(as_ggplot(p))$data, built$data)
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(p, file)
  expect_equal(ggplot2::ggplot_build(as_ggplot(readRDS(file)))$data, built$data)
  expect_identical(ggplot2::theme_get(), theme)
  expect_identical(grDevices::dev.cur(), device)
  plain <- g + ggplot2::labs(title = NULL, subtitle = NULL)
  expect_null(plain$labels$title)
  expect_null(plain$labels$subtitle)
  expect_identical(plot_data(plain), p$data)
  expect_silent(ggplot2::ggplotGrob(plain))
})

test_that("PCA scree and loadings conversions retain units, order and component status", {
  skip_if_not_installed("ggplot2")
  f <- pca_plot_fixture()
  for (preset in c("standard", "monochrome")) {
    p <- plot(f$pca, draw = FALSE, preset = preset)
    g <- as_ggplot(p)
    b <- ggplot2::ggplot_build(g)
    expect_equal(b$data[[1]]$x, 1:3)
    expect_equal(b$data[[1]]$y, 100 * f$pca$variance$Proportion)
    expect_equal(b$data[[2]]$shape, c(16, 16, 1))
    expect_identical(plot_data(g), p$data)
    expect_identical(b$plot$scales$get_scales("shape")$get_labels(), c("Retained", "Not retained"))
    p <- plot(f$pca, type = "loadings", components = 2, preset = preset, draw = FALSE)
    b <- ggplot2::ggplot_build(as_ggplot(p))
    expect_equal(b$data[[1]]$xintercept, 0)
    expect_equal(b$data[[2]]$x, unname(f$pca$loadings[, "PC2"]))
    expect_equal(b$data[[2]]$y, seq_len(nrow(f$pca$loadings)))
    expect_identical(b$layout$panel_params[[1]]$y$get_labels(), rownames(f$pca$loadings))
    expect_equal(b$plot$scales$get_scales("x")$limits, c(-1, 1))
  }
  p <- plot(f$pca, type = "scores", labels = FALSE, draw = FALSE)
  g <- as_ggplot(p)
  expect_length(g$layers, 1)
  expect_equal(nrow(ggplot2::ggplot_build(g)$data[[1]]), 8)
  expect_true(all(ggplot2::ggplot_build(g)$data[[1]]$shape == 16))
  expect_identical(g$labels$colour, NULL)
  expect_equal(ggplot2::ggplot_build(as_ggplot(f$pca, type = "scores", labels = FALSE))$data,
    ggplot2::ggplot_build(g)$data)
})

test_that("PCA payloads cannot fall back to misleading generic graphics", {
  skip_if_not_installed("ggplot2")
  f <- pca_plot_fixture()
  p <- plot(f$pca, type = "scores", groups = f$groups, draw = FALSE)
  expect_error(as_ggplot(p, component = "encoding"), "plot_data")
  expect_error(as_ggplot(p, type = "scree"), "already selects a view")
  expect_error(as_ggplot(p, labels = FALSE), "empty")
  bad <- p; bad$data$components <- NULL
  expect_error(as_ggplot(bad, component = "table"), "selected axes")
  bad <- p; bad$data$encoding <- NULL
  expect_error(as_ggplot(bad), "group colours and shapes")
  bad <- p; bad$data$encoding <- bad$data$encoding[-1, ]
  expect_error(as_ggplot(bad), "group colours and shapes")
  bad <- p; bad$data$table$PC2 <- NULL
  expect_error(as_ggplot(bad), "selected axes")
})

test_that("a one-feature PCA has a drawable scree plot without a spurious line warning", {
  skip_if_not_installed("ggplot2")
  x <- mfrm_pca(mfrm_features(data.frame(ID = letters[1:4], Years = c(1, 2, 5, 9)), "ID", "Years"))
  g <- as_ggplot(x)
  expect_silent(ggplot2::ggplotGrob(g))
  expect_equal(ggplot2::ggplot_build(g)$data[[1]]$y, 100)
  expect_identical(g$guides$guides$colour, "none")
})

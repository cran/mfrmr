residual_plot_fixture <- function() {
  obs <- expand.grid(Person = 1:20, Criterion = letters[1:4])
  obs$StdResidual <- c(sin(1:20), cos(1:20), sin(2 * (1:20)), c(-10:-1, 10:1) / 10)
  analyze_residual_pca(list(obs = obs, facet_names = "Criterion"), mode = "both",
    parallel = TRUE, parallel_reps = 5, seed = 9028)
}

test_that("residual ggplots preserve eigenvalues, cutoffs and saved display encodings", {
  skip_if_not_installed("ggplot2")
  pca <- residual_plot_fixture()
  for (preset in c("standard", "monochrome")) {
    saved <- plot_residual_pca(pca, plot_type = "parallel_scree", preset = preset, draw = FALSE)
    g <- as_ggplot(saved)
    b <- ggplot2::ggplot_build(g)
    tab <- saved$data$data
    expect_equal(b$data[[1]]$x, rep(tab$Component, 3))
    expect_equal(b$data[[1]]$y, c(tab$Eigenvalue, tab$ParallelCutoff, tab$ParallelMean))
    expect_identical(unique(b$data[[1]]$linetype), c("solid", "dashed", "dotted"))
    expect_identical(unique(b$data[[2]]$shape), c(16, 17, NA_real_))
    expect_equal(b$data[[3]]$yintercept, 1)
    expect_identical(g$labels$x, "Component")
    expect_equal(b$layout$panel_params[[1]]$x$get_breaks(), tab$Component)
    expect_identical(g$labels$y, "Eigenvalue")
    expect_identical(plot_data(g), saved$data)
    labels <- b$plot$scales$get_scales("colour")$get_labels()
    expect_identical(labels, saved$data$legend$label[1:3])
    expect_false(any(grepl("residual_permutation|model_bootstrap", labels)))
    expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "no adjustment for scanning")
    expect_silent(grob <- ggplot2::ggplotGrob(g))
    guides <- Filter(function(z) inherits(z, "gtable"),
      grob$grobs[grepl("^guide-box", grob$layout$name)])
    expect_length(guides, 1L)
    expect_lte(grid::convertWidth(sum(guides[[1L]]$widths), "in", valueOnly = TRUE), 4.8)
    file <- tempfile(fileext = ".rds")
    saveRDS(saved, file)
    expect_equal(ggplot2::ggplot_build(as_ggplot(readRDS(file)))$data, b$data)
    unlink(file)
    expect_null((g + ggplot2::labs(title = NULL, subtitle = NULL))$labels$subtitle)
  }
})

test_that("residual scree, excess and loadings keep their own axes and units", {
  skip_if_not_installed("ggplot2")
  pca <- residual_plot_fixture()
  for (type in c("scree", "parallel_excess", "loadings")) {
    saved <- plot_residual_pca(pca, mode = "facet", facet = "Criterion",
      plot_type = type, draw = FALSE)
    b <- ggplot2::ggplot_build(as_ggplot(saved))
    tab <- saved$data$data
    if (type == "loadings") {
      expect_equal(b$data[[1]]$x, tab$Loading)
      expect_equal(as.numeric(b$data[[1]]$y), seq_len(nrow(tab)))
      expect_identical(b$layout$panel_params[[1]]$y$get_labels(), tab$Variable)
    } else {
      expect_equal(as.numeric(b$data[[1]]$x), tab$Component)
      expect_equal(b$data[[1]]$y,
        if (type == "scree") tab$Eigenvalue else tab$ExcessOverParallelCutoff)
    }
  }
})

test_that("incomplete residual payloads cannot fall through to generic bars", {
  skip_if_not_installed("ggplot2")
  pca <- residual_plot_fixture()
  saved <- plot_residual_pca(pca, plot_type = "parallel_scree", draw = FALSE)
  expect_error(as_ggplot(saved, component = "legend"), "plot_data")
  bad <- saved; bad$data$data$Eigenvalue <- NULL
  expect_error(as_ggplot(bad), "saved values")
  bad <- saved; bad$data$legend <- NULL
  expect_error(as_ggplot(bad), "display encoding")
  expect_error(as_ggplot(saved, unknown = TRUE), "empty")
})

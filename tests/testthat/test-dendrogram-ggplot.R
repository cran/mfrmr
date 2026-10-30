hierarchy_plot_fixture <- function(linkage = "average", tied = FALSE) {
  d <- if (tied) data.frame(ID = letters[1:4], X = c(0, 0, 1, 1), Y = c(0, 1, 0, 1)) else
    data.frame(ID = c("z", "a", "R|3", "d", "e", "excluded"),
      X = c(0, 1, 3, 9, 10, NA), Y = c(1, 4, 2, 6, 8, 5))
  mfrm_cluster_hierarchical(mfrm_features(d, "ID", c("X", "Y")),
    k = 3, linkage = linkage, missing = "omit")
}

test_that("dendrogram coordinates match stats dendrogram midpoints, heights and order", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  for (linkage in c("average", "complete")) for (tied in c(FALSE, TRUE)) {
    x <- hierarchy_plot_fixture(linkage, tied)
    p <- plot(x, draw = FALSE, preset = "monochrome")
    g <- as_ggplot(p)
    built <- ggplot2::ggplot_build(g)
    seg <- built$data[[1]]
    tree <- stats::as.dendrogram(x$tree)
    # stats supplies independent node midpoints and member counts.
    expected <- list()
    visit <- function(node, offset = 0) {
      if (is.leaf(node)) return(invisible(NULL))
      offsets <- offset + c(0, attr(node[[1]], "members"))
      positions <- vapply(1:2, function(j)
        offsets[j] + 1 + (attr(node[[j]], "midpoint") %||% 0), numeric(1))
      heights <- vapply(node, attr, numeric(1), which = "height")
      parent <- attr(node, "height")
      expected[[length(expected) + 1L]] <<- rbind(
        c(positions[1], heights[1], positions[1], parent),
        c(positions[2], heights[2], positions[2], parent),
        c(positions[1], parent, positions[2], parent))
      visit(node[[1]], offsets[1]); visit(node[[2]], offsets[2])
    }
    visit(tree)
    canonical <- function(m) m[do.call(order, as.data.frame(m)), , drop = FALSE]
    expect_equal(canonical(unname(as.matrix(seg[c("x", "y", "xend", "yend")]))),
      canonical(do.call(rbind, expected)), ignore_attr = TRUE)
    expect_identical(built$layout$panel_params[[1]]$x$get_labels(), p$data$leaf_order)
    expect_identical(plot_data(g), p$data)
    expect_equal(nrow(seg), 3 * (length(x$tree$order) - 1))
    boxes <- built$data[[2]]
    for (j in seq_len(nrow(boxes))) {
      positions <- which(seq_along(p$data$leaf_order) > boxes$xmin[j] &
        seq_along(p$data$leaf_order) < boxes$xmax[j])
      expect_length(unique(p$data$table$Cluster[positions]), 1)
    }
    expect_equal(nrow(boxes), p$data$k)
    expect_true(all(boxes$linetype == "dashed"))
    expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "Heights do not measure branch support", fixed = TRUE)
    if (tied) expect_match(paste(plot_data(g, "notes")$Text, collapse = "\n"), "Tied heights", fixed = TRUE)
    expect_silent(ggplot2::ggplotGrob(g))
  }
})

test_that("saved trees convert without clustering, recutting or changing session styles", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  x <- hierarchy_plot_fixture()
  p <- plot(x, labels = FALSE, preset = "monochrome", draw = FALSE)
  direct <- as_ggplot(x, labels = FALSE, preset = "monochrome")
  local_mocked_bindings(mfrm_cluster_hierarchical = function(...) stop("Unexpected clustering"),
    cluster_external_features = function(...) stop("Unexpected clustering"))
  local_mocked_bindings(hclust = function(...) stop("Unexpected tree"),
    cutree = function(...) stop("Unexpected partition"), .package = "mfrmr")
  withr::local_options(mfrmr.plot_preset = "invalid")
  device <- grDevices::dev.cur()
  theme <- ggplot2::theme_get()
  file <- tempfile(fileext = ".rds"); on.exit(unlink(file), add = TRUE)
  saveRDS(p, file)
  g <- as_ggplot(readRDS(file), component = "tree")
  expect_equal(ggplot2::ggplot_build(g)$data, ggplot2::ggplot_build(direct)$data)
  expect_identical(plot_data(g), p$data)
  expect_identical(plot_data(g)$excluded_ids, "excluded")
  expect_null(g$scales$get_scales("x")$breaks)
  expect_equal(nrow(ggplot2::ggplot_build(g)$data[[1]]), 12)
  plain <- g + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
  expect_null(plain$labels$caption)
  expect_identical(plot_data(plain), p$data)
  expect_identical(grDevices::dev.cur(), device)
  expect_identical(ggplot2::theme_get(), theme)
  expect_error(as_ggplot(p, component = "table"), "plot_data")
  expect_error(as_ggplot(p, type = "silhouette"), "already selects a view")
  expect_error(as_ggplot(p, labels = TRUE), "empty")
  bad <- p; bad$data$leaf_order <- rev(bad$data$leaf_order)
  expect_error(as_ggplot(bad), "matching saved tree")
  bad <- p; bad$data$tree$merge[1, 1] <- 1
  expect_error(as_ggplot(bad), "invalid merge indices")
  bad <- p; bad$data$k <- 4
  expect_error(as_ggplot(bad), "group count")
})

test_that("zero-height trees retain all leaves and visible group boxes", {
  skip_if_not_installed("ggplot2")
  # A valid saved tree from identical distances; no feature preprocessing needed.
  tree <- stats::hclust(stats::dist(matrix(0, 4, 1)))
  tree$labels <- letters[1:4]
  ids <- tree$labels[tree$order]
  p <- new_mfrm_plot_data("cluster_dendrogram", list(tree = tree, leaf_order = ids,
    table = data.frame(ID = ids, Cluster = unname(stats::cutree(tree, 3)[tree$order])),
    k = 3, labels = TRUE, preset = "monochrome"))
  b <- ggplot2::ggplot_build(as_ggplot(p))
  expect_equal(nrow(b$data[[1]]), 9)
  expect_true(all(b$data[[1]]$y == 0 & b$data[[1]]$yend == 0))
  expect_true(all(b$data[[2]]$ymin < b$data[[2]]$ymax))
  expect_silent(ggplot2::ggplotGrob(as_ggplot(p)))
})


test_that("height ties away from the selected cut do not trigger the cut caution", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cluster")
  x <- hierarchy_plot_fixture(tied = TRUE)
  p <- plot(x, draw = FALSE)
  p$data$k <- 2L
  ids <- p$data$leaf_order
  membership <- stats::cutree(p$data$tree, 2)
  p$data$table$Cluster <- unname(membership[match(ids, p$data$tree$labels)])
  expect_true(anyDuplicated(p$data$tree$height) > 0L)
  expect_false(grepl("Tied heights", paste(plot_data(as_ggplot(p), "notes")$Text, collapse = "\n"), fixed = TRUE))
})

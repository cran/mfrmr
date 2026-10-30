test_that("cluster plots use stored silhouettes and feature summaries without refitting", {
  skip_if_not_installed("cluster")
  input <- data.frame(ID = paste0("R", 1:7), Years = c(1, 2, 3, 12, 13, 14, NA),
    Training = ordered(c(rep("Basic", 3), rep("Advanced", 4)),
      levels = c("Basic", "Intermediate", "Advanced")))
  x <- mfrm_cluster(mfrm_features(input, "ID", c("Years", "Training")), 2, missing = "omit")
  local_mocked_bindings(mfrm_cluster = function(...) stop("Unexpected refit"))
  before <- grDevices::dev.cur()
  s <- plot_data(plot(x, draw = FALSE))
  expect_identical(grDevices::dev.cur(), before)
  expect_identical(s$excluded_ids, "R7")
  expect_equal(s$table$Silhouette, x$membership$Silhouette[match(s$table$ID, x$membership$ID)])
  expect_equal(s$reference_lines$value, mean(x$membership$Silhouette, na.rm = TRUE))
  p <- plot_data(plot(x, type = "profile", feature = "Years", draw = FALSE))
  expect_equal(p$table$Mean, c(2, 13))
  expect_equal(p$table$Median, c(2, 13))
  expect_equal(p$table$N, c(3, 3))
  c <- plot_data(plot(x, type = "profile", feature = "Training", draw = FALSE))
  expect_identical(colnames(c$matrix), levels(input$Training))
  expect_equal(unname(c$matrix), rbind(c(1, 0, 0), c(0, 0, 1)))
  expect_identical(c$excluded_ids, "R7")
  expect_error(plot(x, type = "profile", draw = FALSE), "one selected feature")
  expect_error(plot(x, type = "profile", feature = "ID", draw = FALSE), "one selected feature")
  expect_error(plot(x, feature = "Years", draw = FALSE), "only used")
  expect_error(plot(x, draw = NA), "TRUE or FALSE")
  expect_error(plot(x, labels = 1, draw = FALSE), "TRUE or FALSE")
  expect_error(plot(x, unexpected = TRUE, draw = FALSE), "must be empty")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    for (view in list(list(), list(type = "profile", feature = "Years"),
                     list(type = "profile", feature = "Training"))) {
      expect_s3_class(do.call(as_ggplot, c(list(x), view)), "ggplot")
      payload <- do.call(plot, c(list(x, draw = FALSE), view))
      expect_s3_class(as_ggplot(payload), "ggplot")
    }
    hierarchy <- mfrm_cluster_hierarchical(
      mfrm_features(input, "ID", c("Years", "Training")), 2, missing = "omit")
    expect_s3_class(as_ggplot(hierarchy), "ggplot")
    expect_s3_class(as_ggplot(plot(hierarchy, draw = FALSE)), "ggplot")
  }

  grDevices::pdf(NULL, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(c("mar", "bg", "fg", "cex.axis", "mfrow"))
  expect_no_warning(plot(x, preset = "publication"))
  expect_no_warning(plot(x, type = "profile", feature = "Years"))
  expect_no_warning(plot(x, type = "profile", feature = "Training", preset = "monochrome"))
  expect_equal(graphics::par(names(old)), old)
})

test_that("imputation heatmaps preserve pair fractions, unavailable cells, and requested order", {
  # Fixed formal-result payload: two partitions, one excluded entity.
  ids <- c("b", "a", "c", "excluded")
  mat <- matrix(c(1, .5, 0, NA, .5, 1, .5, NA, 0, .5, 1, NA,
    rep(NA_real_, 4)), 4, dimnames = list(ids, ids))
  x <- structure(list(co_membership = mat, settings = list(imputations = 2L)),
    class = "mfrm_imputed_clusters")
  local_mocked_bindings(mfrm_cluster = function(...) stop("Unexpected refit"))
  all <- plot_data(plot(x, draw = FALSE))
  expect_identical(all$matrix, mat)
  expect_identical(all$excluded_ids, "excluded")
  expect_identical(all$imputations, 2L)
  selection <- c("excluded", "c", "b")
  p <- plot_data(plot(x, ids = selection, draw = FALSE))
  expect_identical(p$ids, selection)
  expect_identical(p$matrix, mat[selection, selection])
  expect_true(all(is.na(p$matrix[1, ])))
  expect_identical(p$matrix[2, 3], 0)
  expect_error(plot(x, ids = c("b", "b")), "distinct entity IDs")
  expect_error(plot(x, ids = "unknown"), "distinct entity IDs")
  expect_error(plot(x, ids = character()), "distinct entity IDs")
  expect_error(plot(x, labels = NA), "TRUE or FALSE")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(as_ggplot(x), "ggplot")
    expect_s3_class(as_ggplot(plot(x, ids = selection, draw = FALSE)), "ggplot")
  }
  grDevices::pdf(NULL, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(c("mar", "bg", "fg", "cex.axis", "mfrow"))
  expect_no_warning(plot(x))
  expect_no_warning(plot(x, ids = "excluded", preset = "monochrome"))
  expect_no_warning(plot(x, ids = "b", labels = FALSE))
  expect_equal(graphics::par(names(old)), old)
  big_ids <- paste0("R", 1:51)
  x$co_membership <- matrix(1, 51, 51, dimnames = list(big_ids, big_ids))
  big <- plot_data(plot(x, draw = FALSE))
  expect_false(big$labels)
  expect_identical(big$ids, big_ids)
  expect_equal(dim(big$matrix), c(51, 51))
})

test_that("hierarchies match weighted mixed-feature distances and retain omitted IDs", {
  skip_if_not_installed("cluster")
  data <- data.frame(ID = c("z", "a", "R|3", "d", "e", "excluded"),
    Years = c(0, 1, 3, 9, 10, NA),
    Specialty = factor(c("A", "A", "B", "B", "B", "A")),
    Training = ordered(c("L", "H", "H", "L", "H", "L"), levels = c("L", "H")),
    Certified = c(FALSE, FALSE, TRUE, TRUE, FALSE, TRUE))
  x <- mfrm_features(data, "ID", names(data)[-1])
  weights <- c(Years = 3, Specialty = 2, Training = 1, Certified = 1)
  # Independent mixed-distance calculation with declared ordinal/binary coding.
  d <- (3 * abs(outer(data$Years[1:5], data$Years[1:5], "-")) / 10 +
    2 * outer(data$Specialty[1:5], data$Specialty[1:5], "!=") +
    outer(data$Training[1:5], data$Training[1:5], "!=") +
    outer(data$Certified[1:5], data$Certified[1:5], "!=")) / 7
  for (linkage in c("average", "complete")) {
    out <- mfrm_cluster_hierarchical(x, 3, weights, missing = "omit", linkage = linkage)
    reference <- stats::hclust(stats::as.dist(d), method = linkage)
    expected <- stats::cutree(reference, k = 3)
    expect_s3_class(out, "mfrm_hierarchical_clusters")
    expect_s3_class(out, "mfrm_clusters")
    expect_identical(out$tree$merge, reference$merge)
    expect_equal(out$tree$height, reference$height, tolerance = 1e-12)
    expect_identical(out$tree$labels, data$ID[1:5])
    expect_identical(out$membership$ID, data$ID)
    expect_equal(out$membership$Cluster[1:5], unname(expected))
    expect_equal(out$membership$Silhouette[1:5],
      unname(cluster::silhouette(expected, stats::as.dist(d))[, "sil_width"]))
    expect_true(all(is.na(out$membership[6, c("Cluster", "Silhouette")])))
    expect_null(out$medoids)
    expect_false("Medoid" %in% names(out$membership))
    expect_equal(out$settings$excluded, 1)
    expect_identical(out$settings$linkage, linkage)
    expect_equal(sum(summary(out)$N), 5)
    expect_equal(out$profiles$numeric$Mean,
      as.numeric(tapply(data$Years[1:5], expected, mean)))
    expect_output(print(out), paste0(linkage, " linkage"))
  }
  expect_error(mfrm_cluster_hierarchical(x, 2), "missing values")
  expect_error(mfrm_cluster_hierarchical(x, 2, linkage = "ward.D2"), "average.*complete")
  expect_error(mfrm_cluster_hierarchical(x, 2, weights = weights * 0, missing = "omit"), "positive finite")
  comparison <- mfrm_cluster_compare(list(
    PAM = mfrm_cluster(x, 3, weights, missing = "omit"), Hierarchical = out))
  expect_identical(comparison$analysis_summary$Method, c("PAM", "Hierarchical"))
  expect_identical(comparison$analysis_summary$Linkage, c(NA_character_, "complete"))
  expect_equal(comparison$comparisons$Pairs, 10)
})

test_that("dendrograms preserve tied cuts and use the saved tree without refitting", {
  skip_if_not_installed("cluster")
  square <- data.frame(ID = letters[1:4], X = c(0, 0, 1, 1), Y = c(0, 1, 0, 1))
  x <- mfrm_cluster_hierarchical(mfrm_features(square, "ID", c("X", "Y")), 3)
  expect_true(anyDuplicated(x$tree$height) > 0L)
  expect_equal(x$membership$Cluster, unname(stats::cutree(x$tree, 3)))
  expect_equal(sum(x$cluster_summary$N == 1L), 2L)
  expect_true(all(x$membership$Silhouette[x$membership$Cluster %in%
    x$cluster_summary$Cluster[x$cluster_summary$N == 1L]] == 0))
  local_mocked_bindings(cluster_external_features = function(...) stop("Unexpected refit"))
  device <- grDevices::dev.cur()
  view <- plot_data(plot(x, draw = FALSE))
  expect_identical(grDevices::dev.cur(), device)
  expect_identical(view$tree, x$tree)
  expect_identical(view$leaf_order, x$tree$labels[x$tree$order])
  expect_identical(view$table$ID, view$leaf_order)
  expect_equal(view$table$Cluster, x$membership$Cluster[match(view$table$ID, square$ID)])
  expect_equal(view$k, 3)
  expect_match(plot_data(plot(x, type = "silhouette", draw = FALSE))$subtitle, "average linkage")
  expect_match(plot_data(plot(x, type = "profile", feature = "X", draw = FALSE))$subtitle, "average linkage")
  expect_error(plot(x, labels = NA), "TRUE or FALSE")
  expect_error(plot(x, feature = "X"), "only used")
  expect_error(plot(x, ids = "a"), "must be empty")
  grDevices::pdf(NULL, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(c("mar", "bg", "fg", "cex.axis", "mfrow", "xpd"))
  expect_no_warning(plot(x))
  expect_no_warning(plot(x, labels = FALSE, preset = "monochrome"))
  expect_no_warning(plot(x, type = "silhouette"))
  expect_no_warning(plot(x, type = "profile", feature = "X"))
  expect_equal(graphics::par(names(old)), old)
})

test_that("hierarchical imputations retain each tree and all-imputation pair denominators", {
  skip_if_not_installed("cluster")
  skip_if_not_installed("mice")
  original <- data.frame(ID = 1:6, X = c(0, 1, NA, 9, 10, NA))
  first <- second <- original
  first$X[3] <- 2
  second$X[3] <- 8
  where <- is.na(original)
  where[6, "X"] <- FALSE
  model <- mice::as.mids(rbind(cbind(.imp = 0L, original),
    cbind(.imp = 1L, first), cbind(.imp = 2L, second)), where = where)
  x <- mfrm_features(original, "ID", "X")
  cells <- subset(x$missing, ID == "3")
  out <- mfrm_cluster_imputed(x, model, cells, 2, missing = "omit", method = "hierarchical")
  complete <- mfrm_cluster_imputed(x, model, cells, 2, missing = "omit",
    method = "hierarchical", linkage = "complete")
  pam <- mfrm_cluster_imputed(x, model, cells, 2, missing = "omit")
  expect_identical(out$settings$method, "Hierarchical")
  expect_identical(out$settings$linkage, "average")
  expect_output(print(out), "Hierarchical.*average linkage")
  expect_true(all(vapply(out$analyses, inherits, logical(1), "mfrm_hierarchical_clusters")))
  expect_true(all(vapply(complete$analyses, function(a) a$tree$method == "complete", logical(1))))
  expect_identical(out$imputation_model, model)
  expect_equal(out$co_membership["1", "3"], 0.5)
  expect_equal(out$co_membership["3", "4"], 0.5)
  expect_equal(out$co_membership["1", "5"], 0)
  expect_true(all(is.na(out$co_membership["6", ])))
  for (i in 1:2) {
    expect_equal(out$analyses[[i]]$membership$Cluster[1:5],
      unname(stats::cutree(out$analyses[[i]]$tree, k = 2)))
  }
  comparison <- mfrm_cluster_compare(list(PAM = pam, Average = out, Complete = complete))
  expect_equal(nrow(comparison$comparisons), 6)
  expect_equal(comparison$analysis_summary$Linkage, rep(c(NA_character_, "average", "complete"), each = 2))
  expect_equal(plot_data(plot(out, draw = FALSE))$matrix, out$co_membership)
  expect_match(plot_data(plot(out, draw = FALSE))$title, "average linkage")
  expect_identical(plot_data(plot(out, draw = FALSE))$method, "Hierarchical")
  expect_identical(plot_data(plot(out, draw = FALSE))$linkage, "average")
  expect_identical(plot_data(plot(out$analyses[[1]], draw = FALSE))$excluded_ids, "6")
  expect_error(mfrm_cluster_imputed(x, model, cells, 2, linkage = "average"), "PAM has no linkage")
  expect_error(mfrm_cluster_imputed(x, model, cells, 2, method = "hierarchical", linkage = "ward.D2"), "average.*complete")
  calls <- 0L
  fit <- mfrm_cluster_hierarchical
  local_mocked_bindings(mfrm_cluster_hierarchical = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("completion failed")
    fit(...)
  })
  expect_error(mfrm_cluster_imputed(x, model, cells, 2, missing = "omit", method = "hierarchical"),
    "Imputation 2: completion failed")
})

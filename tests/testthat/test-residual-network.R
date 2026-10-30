test_that("observation-only consumers do not compute unrelated diagnostics", {
  x <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
  fit <- x$fit
  # Compare the full and observation-only paths on the same platform. Saved
  # diagnostics were generated elsewhere and can differ in floating-point bits.
  diagnostics <- diagnose_mfrm(fit)
  routes <- list(
    pca = function(d) analyze_residual_pca(if (is.null(d)) fit else d, mode = "both"),
    q3 = function(d) q3_statistic(fit, d),
    network = function(d) mfrm_residual_network(fit, d, facet = "Criterion"),
    matrix = function(d) plot_residual_matrix(fit, d, draw = FALSE),
    qq = function(d) plot_residual_qq(fit, d, draw = FALSE))
  if (requireNamespace("igraph", quietly = TRUE)) {
    routes$halo <- function(d) rater_halo_network_analysis(fit, d)
  }
  reference <- lapply(routes, function(f) f(diagnostics))
  local_mocked_bindings(
    diagnose_mfrm = function(...) stop("Unrelated full diagnostics requested"),
    compute_mml_parameter_covariance = function(...) stop("Unused covariance requested"))
  for (name in names(routes)) {
    expect_identical(routes[[name]](NULL), reference[[name]], info = name)
  }
  # Normalization of the historical JMLE label must not mutate the caller's fit.
  fit$config$method <- "JMLE"
  before <- serialize(fit, NULL)
  expect_identical(routes$pca(NULL), reference$pca)
  expect_identical(serialize(fit, NULL), before)
  class(fit) <- c("mfrm_imported_fit", class(fit))
  for (name in names(routes)) expect_error(routes[[name]](NULL), "imported fit statistics")
})

test_that("agreement consumers omit unused marginal diagnostics", {
  x <- readRDS(system.file("examples", "operational-mml.rds", package = "mfrmr"))
  diagnostics <- diagnose_mfrm(x$fit)
  expected <- interrater_agreement_table(x$fit, diagnostics)
  if (requireNamespace("igraph", quietly = TRUE))
    expected_network <- rater_network_analysis(x$fit, diagnostics)
  local_mocked_bindings(calc_marginal_fit_bundle = function(...) stop("Unused marginal diagnostics"))
  expect_identical(interrater_agreement_table(x$fit), expected)
  if (requireNamespace("igraph", quietly = TRUE))
    expect_identical(rater_network_analysis(x$fit), expected_network)
})

residual_network_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    dat <- load_mfrmr_data("example_core")
    dat <- dat[dat$Person %in% unique(dat$Person)[1:15], ]
    fit <- suppressWarnings(fit_mfrm(dat, "Person", c("Rater", "Criterion"), "Score",
      method = "JML", maxit = 100))
    diagnostics <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "legacy"))
    cached <<- list(fit = fit, diagnostics = diagnostics)
    cached
  }
})

test_that("residual network agrees with independent aggregation and common-person geometry", {
  z <- residual_network_fixture()
  n <- mfrm_residual_network(z$fit, z$diagnostics, facet = "Criterion",
    missing = "complete", clusters = 2)
  o <- z$diagnostics$obs
  wide <- tapply(o$StdResidual, list(as.character(o$Person), as.character(o$Criterion)), mean)
  reference <- stats::cor(wide[complete.cases(wide), ])
  expect_equal(n$matrix, reference, tolerance = 1e-12)
  normalized <- scale(wide[complete.cases(wide), ], center = TRUE, scale = FALSE)
  normalized <- sweep(normalized, 2, sqrt(colSums(normalized^2)), "/")
  distance <- dist(t(normalized))
  expect_equal(as.matrix(n$clustering$distance), as.matrix(distance), tolerance = 1e-12)
  expect_equal(n$clustering$membership, cutree(hclust(distance, method = "average"), k = 2))
  expect_equal(unname(n$cell_counts), unname(table(as.character(o$Person), as.character(o$Criterion))), ignore_attr = TRUE)
  p <- plot(z$fit, type = "local_dependence", diagnostics = z$diagnostics, facet = "Criterion", draw = FALSE)
  expect_equal(p$data$matrix, n$matrix, tolerance = 1e-12)
  threshold <- mfrm_residual_network(z$fit, z$diagnostics, facet = "Criterion",
    missing = "complete", clusters = 2, min_abs = 1)
  expect_equal(n$clustering, threshold$clustering)
  expect_identical(n$matrix, threshold$matrix)
  expect_equal(nrow(threshold$nodes), ncol(reference))
  expect_equal(nrow(threshold$edges), 0)
  expect_equal(threshold$nodes$DisplayedDegree, rep(0L, ncol(reference)))
})

test_that("missing associations and constant levels are not replaced by zero", {
  z <- residual_network_fixture()
  d <- z$diagnostics
  levels <- sort(unique(as.character(d$obs$Criterion)))
  persons <- sort(unique(as.character(d$obs$Person)))
  d$obs$StdResidual[d$obs$Criterion == levels[1] & !d$obs$Person %in% persons[1:2]] <- NA_real_
  a <- mfrm_residual_network(z$fit, d, facet = "Criterion")
  bad <- a$pairs$Level1 == levels[1] | a$pairs$Level2 == levels[1]
  expect_true(all(is.na(a$pairs$ResidualCor[bad])))
  expect_true(all(a$pairs$N[bad] == 2L))
  expect_true(all(a$pairs$Reason[bad] == "Insufficient shared persons"))
  b <- mfrm_residual_network(z$fit, d, facet = "Criterion", missing = "complete", clusters = 2)
  expect_equal(sum(b$person_usage$Used), 2L)
  expect_identical(b$clustering$status, "unavailable")
  expect_null(b$clustering$tree)
  expect_error(plot(b, type = "dendrogram", draw = FALSE), "No clustering tree")
  d <- z$diagnostics; d$obs$StdResidual[d$obs$Criterion == levels[1]] <- 0
  b <- mfrm_residual_network(z$fit, d, facet = "Criterion", missing = "complete", clusters = 2)
  expect_identical(b$clustering$status, "unavailable")
  expect_true(any(grepl("No residual variation", b$pairs$Reason)))
  expect_error(mfrm_residual_network(z$fit, d, clusters = 2), "explicit")
  expect_error(mfrm_residual_network(z$fit, d, min_pairs = 2), "at least 3")
  expect_error(mfrm_residual_network(z$fit, d, min_abs = 2), "min_abs")
})

test_that("saved plots and report exports retain residual evidence", {
  z <- residual_network_fixture()
  n <- mfrm_residual_network(z$fit, z$diagnostics, facet = "Criterion", missing = "complete", clusters = 2)
  path <- tempfile(); saveRDS(n, path); saved <- readRDS(path)
  for (type in c("network", "matrix", "dendrogram")) {
    fig <- plot(saved, type = type, draw = FALSE)
    expect_s3_class(fig, "mfrm_plot_data")
    expect_identical(fig$data$matrix, n$matrix)
    expect_identical(fig$data$clustering, n$clustering)
  }
  pdf <- tempfile(fileext = ".pdf"); grDevices::pdf(pdf)
  tryCatch(for (type in c("network", "matrix", "dendrogram")) plot(saved, type = type),
    finally = grDevices::dev.off())
  expect_true(file.info(pdf)$size > 0)
  if (requireNamespace("ggplot2", quietly = TRUE))
    expect_s3_class(as_ggplot(plot(saved, type = "matrix", draw = FALSE)), "ggplot")
  expect_error(as_ggplot(plot(saved, type = "network", draw = FALSE)), "matrix.*only")
  expect_error(as_ggplot(plot(saved, type = "dendrogram", draw = FALSE)), "matrix.*only")
  r <- mfrm_report(saved)
  expect_identical(r$tables$pairs, n$pairs)
  expect_identical(build_summary_table_bundle(r)$tables$nodes, n$nodes)
  expect_identical(build_summary_table_bundle(saved)$tables$pairs, n$pairs)
  expect_s3_class(apa_table(r, which = "pairs"), "apa_table")
  out <- tempfile(); dir.create(out)
  expect_silent(export_summary_appendix(r, output_dir = out))
})

test_that("opposite residual profiles remain negative edges and distant nodes", {
  z <- residual_network_fixture(); d <- z$diagnostics
  persons <- sort(unique(as.character(d$obs$Person)))
  levels <- sort(unique(as.character(d$obs$Criterion)))
  signal <- match(as.character(d$obs$Person), persons)
  d$obs$StdResidual <- signal * ifelse(d$obs$Criterion == levels[1], -1, 1)
  n <- mfrm_residual_network(z$fit, d, facet = "Criterion", missing = "complete", clusters = 2)
  expect_equal(n$matrix[1, 2], -1, tolerance = 1e-12)
  expect_equal(as.matrix(n$clustering$distance)[1, 2], 2, tolerance = 1e-12)
  expect_true(any(n$edges$ResidualCor < 0))
  expect_false(n$nodes$Cluster[1] == n$nodes$Cluster[2])
})

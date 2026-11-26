test_that("Fair Average prints the selected reference plan and native intervals", {
  d <- load_mfrmr_data("example_core")
  fits <- list(read_example_bundle("operational-mml")$fit,
    read_example_bundle()$fit,
    fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method = "JML_jeffreys"))
  for (f in fits) {
    d <- f$encoding$retained_data
    fair <- fair_average_table(f, reference_plan = unique(d[c("Rater", "Criterion")]), fair_se = TRUE)
    output <- paste(capture.output(print(fair)), collapse = "\n")
    expect_match(output, "Common reference plan")
    expect_match(output, "Lower")
    expect_equal(confint(fair), fair$reference_scores)
    expect_identical(fitted(f, type = "prob"), fitted(f, type = "probs"))
    a <- data.frame(theta = 0, Rater = "R01", Criterion = "Content")
    b <- a; names(b)[1L] <- "Theta"
    expect_identical(predict(f, a), predict(f, b))
    b$theta <- 1
    expect_error(predict(f, b), "different abilities")
    expect_equal(confint(f, "facets", .9), confint(f, component = "facets", level = .9))
  }
})

test_that("cut inversion tolerates unavailable facet rows and excluded persons", {
  d <- load_mfrmr_data("example_core")
  d$Weight <- ifelse(d$Person == d$Person[1L], 0, 1)
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", weight = "Weight")
  out <- .tidymfrm_cut_decisions(f, 2.5)
  expect_true(is.finite(out$CutLogit))
  expect_equal(out$CutLogit, out$CutLogitMin)
  # Controlled display-state check: one unavailable level must not erase all thresholds.
  f$facets$Estimate[1L] <- NA_real_
  out <- .tidymfrm_cut_decisions(f, 2.5)
  base <- sum(vapply(split(f$facets, f$facets$Facet), function(z)
    z$Sign[1L] * mean(z$Estimate[is.finite(z$Estimate)]), 0))
  expected <- function(theta) {
    logits <- (0:3) * (theta + base) - f$step_cumulative[1L, ]
    prob <- exp(logits - max(logits)); sum((1:4) * prob) / sum(prob)
  }
  expect_equal(expected(out$CutLogit), 2.5, tolerance = 1e-8)
})

test_that("all packaged datasets share the public loader without losing their table roles", {
  catalogue <- load_mfrmr_data("list", details = TRUE)
  expect_equal(nrow(catalogue), 13L)
  expect_identical(catalogue$Key, load_mfrmr_data("list"))
  for (key in catalogue$Key) {
    d <- load_mfrmr_data(key)
    expect_identical(d, getExportedValue("mfrmr", .mfrm_dataset_names[[key]]))
    expect_equal(nrow(d), catalogue$Rows[catalogue$Key == key])
  }
  expect_false("Score" %in% names(load_mfrmr_data("example_operational_design")))
  expect_error(load_mfrmr_data("example_core", details = TRUE), "only available")
  expect_false("read_mfrm_posterior" %in% getNamespaceExports("mfrmr"))
})

test_that("base descriptive statistics and residual principal components retain their definitions", {
  x <- .mfrm_describe_numeric(c(1, 2, 3, NA))
  expect_equal(unname(unlist(x[c("n", "mean", "sd", "median", "min", "max", "range", "skew", "kurtosis", "se")])),
    c(3, 2, 1, 2, 1, 3, 2, 0, -7/3, 1/sqrt(3)))
  expect_true(is.na(.mfrm_describe_numeric(rep(1, 3))$skew))
  C <- matrix(c(1, .4, .4, 1), 2L, dimnames = list(c("A", "B"), c("A", "B")))
  pc <- compute_principal_components(C, 2L)$pca
  expect_equal(pc$values, c(1.4, .6))
  L <- unclass(pc$loadings)
  expect_equal(tcrossprod(L), C)
})

test_that("reports reuse matching reliability and include current bias outputs", {
  d <- load_mfrmr_data("example_core")
  wf <- step_mfrm_diagnose(mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score",
    fit = read_example_bundle()$fit))
  expected <- wf$diagnostics$report_statistics$reliability
  testthat::local_mocked_bindings(diagnose_mfrm = function(...) stop("unexpected recomputation"))
  expect_equal(.mfrm_diagnostic_report_statistics(wf)$reliability, expected)
  out <- paste(capture.output(print(estimate_bias(wf$fit))), collapse = "\n")
  expect_match(out, "Pairwise contrasts")
  expect_match(out, "RajuSA-Welch")
})

test_that("the explicit check-core limit caps workers without changing user options", {
  before <- options()
  withr::local_envvar(c(`_R_CHECK_LIMIT_CORES_` = "true"))
  expect_identical(.mfrm_threads(8L), 2L)
  expect_identical(.mfrm_threads(1L), 1L)
  expect_error(.mfrm_threads(1.5), "positive integer")
  expect_error(.mfrm_threads(NA_real_), "positive integer")
  expect_identical(options(), before)
})


test_that("ordinary saved calibrations agree with their generation recipes", {
  skip_on_cran()
  recipes <- new.env(parent = baseenv())
  sys.source(system.file("examples", "ordinary-fits.R", package = "mfrmr",
    mustWork = TRUE), envir = recipes)
  regenerated <- recipes$ordinary_fit_examples()
  for (name in names(regenerated)) {
    saved <- read_example_bundle(name)
    current <- regenerated[[name]]
    expect_identical(saved$dataset, current$dataset)
    expect_identical(saved$fit$schema_version, current$fit$schema_version)
    expect_equal(coef(saved$fit), coef(current$fit), tolerance = 1e-6)
    expect_equal(vcov(saved$fit), vcov(current$fit), tolerance = 1e-6)
    expect_equal(predict(saved$fit, type = "probs"),
      predict(current$fit, type = "probs"), tolerance = 1e-6)
    expect_equal(saved$diagnostics$measures, current$diagnostics$measures,
      tolerance = 1e-6)
  }
})

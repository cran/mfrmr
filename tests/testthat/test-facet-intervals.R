local({
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score")
  intervals <- mfrm_facet_intervals(f, "Rater", method = "sandwich")

  test_that("known unit population SD preserves ordinary fixed-population contrasts", {
    local_mocked_bindings(compute_mml_parameter_covariance = function(...)
      list(status = "ok", cov = intervals$model_parameter_covariance))
    C <- matrix(c(1, -1, 0, 0), 1, dimnames = list("Difference", rownames(intervals$contrasts)))
    native <- mfrm_facet_intervals(f, "Rater", C)
    std <- suppressWarnings(mfrm_facet_intervals(f, "Rater", C, scale = "standardized"))
    common <- setdiff(names(native$table), "InferenceReview")
    expect_identical(std$table[common], native$table[common])
    expect_match(std$table$InferenceReview, "Population-SD-standardized")
    expect_equal(std$covariance, native$covariance, ignore_attr = TRUE)
    expect_false(std$settings$population_sd_estimated)
    expect_match(plot_data(plot(std, draw = FALSE))$subtitle, "Fixed population SD = 1")
  })

  test_that("Person scores and sandwich reproduce direct full-covariance calculations", {
    h <- compute_mml_parameter_covariance(f)
    expect_equal(intervals$parameter_covariance,
      h$cov %*% crossprod(intervals$person_scores) %*% h$cov, ignore_attr = TRUE)
    idx <- build_indices(f$prep)
    gradient <- mfrm_grad_mml(f$opt$par, idx, f$config, h$sizes, gauss_hermite_normal(31))
    expect_lt(max(abs(colSums(intervals$person_scores) + gradient)), 1e-9)
    ordinary <- mfrm_facet_intervals(f, "Rater")
    expect_equal(ordinary$table$SE, intervals$table$ModelSE)
    expect_null(ordinary$person_scores)
    j <- matrix(c(1, -1, 0, 0), 1, dimnames = list("R01 minus R02", rownames(intervals$contrasts)))
    pair <- mfrm_facet_intervals(f, "Rater", j, method = "sandwich")
    expect_equal(pair$table$SE^2, drop(j %*% intervals$covariance %*% t(j)))
    expect_equal(pair$table$Estimate, sum(j * intervals$table$Estimate))
    expect_equal(mfrm_facet_intervals(f, "Rater", j[, 4:1, drop = FALSE],
      method = "sandwich")$table, pair$table)
    expect_equal(pair$table$Lower, pair$table$Estimate - qnorm(.975) * pair$table$SE)
    expect_equal(summary(intervals), intervals$table)
    saved <- tempfile(fileext = ".rds"); saveRDS(intervals, saved)
    expect_identical(readRDS(saved), intervals)
  })

  test_that("larger clusters aggregate scores before their outer products", {
    ids <- rownames(intervals$person_scores)
    groups <- data.frame(Person = rev(ids), Cluster = rep(seq_len(16), each = 3))
    clustered <- mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups)
    sums <- rowsum(intervals$person_scores, groups$Cluster[match(ids, groups$Person)], reorder = FALSE)
    expect_equal(clustered$cluster_scores, sums)
    v <- intervals$model_parameter_covariance
    expect_equal(clustered$parameter_covariance, v %*% crossprod(sums) %*% v, ignore_attr = TRUE)
    adjusted <- mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups, adjust = TRUE)
    expect_equal(adjusted$table$SE, sqrt(16 / 15) * clustered$table$SE)
    expect_equal(adjusted$settings$adjustment_factor, 16 / 15)
    groups$Cluster <- rep(1:2, each = 24)
    too_few <- mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups)
    expect_true(all(too_few$table$Status == "insufficient_cluster_rank"))
    expect_false(any(too_few$table$CIEligible))
    expect_true(all(grepl("do not split dependent clusters", too_few$table$InferenceReview)))
    expect_true(all(is.na(too_few$table$Lower)))
    expect_equal(too_few$table$Estimate, intervals$table$Estimate)
    expect_true(all(is.finite(too_few$table$ModelLower)))
    expect_error(mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups[-1, ]), "exactly once")
    groups$Cluster <- 1
    expect_error(mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups), "at least two")
    groups$Cluster[1] <- NA
    expect_error(mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups), "nonmissing")
  })

  test_that("unsupported inference cannot acquire a robust label", {
    expect_error(mfrm_facet_intervals(f, "Person"), "non-person")
    expect_error(mfrm_facet_intervals(f, "Rater", adjust = TRUE), "only to method")
    expect_error(mfrm_facet_intervals(f, "Rater", level = 1), "level")
    bad <- f; bad$config$method <- "JML"
    expect_error(mfrm_facet_intervals(bad, "Rater"), "inference-ready")
    bad <- f; bad$config$population_spec$active <- TRUE
    expect_error(mfrm_facet_intervals(bad, "Rater"), "Population coding")
    bad <- f; bad$config$model <- "GPCM"
    # Changing only the model label cannot create a valid GPCM calibration.
    expect_error(mfrm_facet_intervals(bad, "Rater"))
    bad <- f; bad$prep$data$Weight[1] <- 2
    expect_error(mfrm_facet_intervals(bad, "Rater"), "unit weights")
    bad <- f; bad$readiness$fit$InferenceReady <- FALSE
    expect_error(mfrm_facet_intervals(bad, "Rater"), "inference-ready")
    covariance <- compute_mml_parameter_covariance(f); covariance$status <- "regularized"
    local_mocked_bindings(compute_mml_parameter_covariance = function(...) covariance)
    expect_warning(unavailable <- mfrm_facet_intervals(f, "Rater"), "Unregularized")
    expect_true(all(unavailable$table$Status == "information_unavailable"))
    expect_true(all(is.na(unavailable$table$ModelSE)))
  })

  test_that("information-workspace refusal preserves points and reasons across accessors", {
    p <- length(coef(f))
    withr::local_options(mfrmr.max_information_bytes = 16 * p^2)
    expect_warning(ci <- mfrm_facet_intervals(f, "Rater"), "workspace estimate")
    expect_identical(ci$table$Estimate, intervals$table$Estimate)
    expect_true(all(ci$table$Status == "information_unavailable"))
    expect_false(any(ci$table$CIEligible))
    expect_true(all(is.na(ci$table[c("SE", "Lower", "Upper", "ModelSE", "ModelLower", "ModelUpper")])))
    expect_true(all(is.na(ci$covariance)))
    expect_identical(rownames(ci$covariance), ci$table$Target)
    expect_null(ci$parameter_covariance)
    expect_match(ci$unavailable_reason, "workspace estimate")
    expect_false(any(ci$checks$Passed))
    expect_warning(V <- vcov(f), "workspace estimate")
    expect_true(all(is.na(V)))
    expect_identical(rownames(V), names(coef(f)))
    expect_identical(attr(V, "unavailable_reason"), ci$unavailable_reason)
    expect_warning(bounds <- confint(f), "workspace estimate")
    expect_true(all(is.na(bounds)))
    diagnostics <- attr(bounds, "diagnostics")
    expect_identical(diagnostics$Estimate, unname(coef(f)))
    expect_true(all(grepl("workspace estimate", diagnostics$InferenceReview)))
    expect_false(any(diagnostics$CIEligible))
    expect_identical(plot_data(ci)$table, ci$table)
    expect_identical(apa_table(ci)$table$InferenceReview, ci$table$InferenceReview)
    res <- mfrm_results(f, intervals = ci, include = c("fit", "plots"), compute = "never")
    expect_match(mfrm_report(res)$markdown, "workspace estimate")
    withr::local_options(mfrmr.max_information_bytes = 1)
    expect_error(suppressWarnings(vcov(f)), "requested NA covariance matrix exceeds")
    withr::local_options(mfrmr.max_information_bytes = -1)
    expect_error(mfrm_facet_intervals(f, "Rater"), "finite positive byte count")
  })

  test_that("failed information cannot be bypassed by sandwich or bad cluster inputs", {
    withr::local_options(mfrmr.max_information_bytes = 1)
    local_mocked_bindings(mfrm_person_likelihood_scores = function(...) stop("unexpected score calculation"))
    ids <- f$prep$levels$Person
    groups <- data.frame(Person = rev(ids), Cluster = rep(seq_len(16), each = 3))
    expect_warning(ci <- mfrm_facet_intervals(f, "Rater", method = "sandwich",
      clusters = groups, adjust = TRUE), "workspace estimate")
    expect_identical(ci$table$Estimate, intervals$table$Estimate)
    expect_true(all(ci$table$Status == "information_unavailable"))
    expect_true(all(is.na(ci$table$Lower)))
    expect_true(all(is.na(ci$table$ModelLower)))
    expect_null(ci$person_scores)
    expect_null(ci$cluster_scores)
    expect_true(is.na(ci$settings$cluster_score_rank))
    expect_equal(ci$settings$clusters, 16)
    expect_equal(ci$settings$adjustment_factor, 16 / 15)
    expect_identical(ci$clusters$Person, ids)
    expect_identical(ci$clusters$Cluster, as.character(groups$Cluster[match(ids, groups$Person)]))
    expect_error(mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups[-1, ]), "exactly once")
    groups$Cluster <- 1
    expect_error(mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups), "at least two")
    groups$Cluster[1] <- NA
    expect_error(mfrm_facet_intervals(f, "Rater", method = "sandwich", clusters = groups), "nonmissing")
  })

  test_that("rejected covariance never supplies model or standardized intervals", {
    p <- length(f$opt$par)
    bad_covariances <- list(diag(p - 1L), matrix(NA_real_, p, p), -diag(p),
      matrix(0, p, p), NULL)
    for (bad in bad_covariances) local({
      local_mocked_bindings(compute_mml_parameter_covariance = function(...)
        list(status = "ok", cov = bad, detail = "Controlled covariance failure."))
      C <- matrix(c(1, -1, 0, 0), 1, dimnames = list("Difference", rownames(intervals$contrasts)))
      ci <- suppressWarnings(mfrm_facet_intervals(f, "Rater", C, scale = "standardized"))
      expect_identical(ci$table$Status, "information_unavailable")
      expect_false(ci$table$CIEligible)
      expect_true(is.na(ci$table$ModelLower))
      expect_match(ci$table$InferenceReview, "Controlled covariance failure")
      expect_null(ci$model_parameter_covariance)
    })
    local_mocked_bindings(compute_mml_parameter_covariance = function(...)
      list(status = "ok", regularized = TRUE, cov = diag(p)))
    expect_warning(ci <- mfrm_facet_intervals(f, "Rater"), "Unregularized")
    expect_true(all(is.na(ci$table$ModelSE)))
  })

  test_that("PCM, anchors and interaction gradients preserve their fitted constraints", {
    for (args in list(list(model = "PCM", step_facet = "Criterion"),
                     list(anchors = data.frame(Facet = "Rater", Level = "R01", Anchor = .1)),
                     list(facet_interactions = "Rater:Criterion"))) {
      fit <- suppressWarnings(do.call(fit_mfrm, c(list(data = d, person = "Person",
        facets = c("Rater", "Criterion"), score = "Score"), args)))
      result <- mfrm_facet_intervals(fit, "Rater", method = "sandwich")
      expect_true(all(result$table$Status %in% c("available", "fixed")))
      if (!is.null(args$anchors)) {
        fixed <- result$table$Target == "R01"
        expect_equal(result$table$Estimate[fixed], .1)
        expect_equal(result$table$SE[fixed], 0)
        expect_true(is.na(result$table$Lower[fixed]))
        expect_false(result$table$CIEligible[fixed])
        expect_match(result$table$InferenceReview[fixed], "Constraint-fixed")
      }
    }
  })

  test_that("plots and summaries retain methods without refitting", {
    local_mocked_bindings(fit_mfrm = function(...) stop("unexpected refit"),
      compute_mml_parameter_covariance = function(...) stop("unexpected recomputation"))
    device <- grDevices::dev.cur()
    result <- plot(intervals, draw = FALSE)
    expect_identical(grDevices::dev.cur(), device)
    expect_equal(plot_data(result)$table, intervals$table)
    expect_true(plot_data(result)$comparison)
    expect_false(plot_data(plot(intervals, comparison = FALSE, draw = FALSE))$comparison)
    expect_error(plot(intervals, draw = NA), "TRUE or FALSE")
    if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(as_ggplot(result), "ggplot")
    file <- tempfile(fileext = ".pdf"); grDevices::pdf(file)
    before <- graphics::par("mar"); plot(intervals)
    expect_equal(graphics::par("mar"), before)
    grDevices::dev.off()
    expect_gt(file.info(file)$size, 0)
  })
})

test_that("individual marginal scores match independent continuous integration", {
  with_preserved_rng_seed(9812, {
    d <- expand.grid(Person = paste0("P", 1:40), Rater = c("A", "B"),
      Task = c("X", "Y"), stringsAsFactors = FALSE)
    theta <- rnorm(40)
    r <- ifelse(d$Rater == "A", 1, -1); t <- ifelse(d$Task == "X", 1, -1)
    eta <- theta[match(d$Person, unique(d$Person))] - .3 * r - .4 * t
    logits <- cbind(0, eta + .6, 2 * eta)
    prob <- exp(logits - apply(logits, 1, max)); prob <- prob / rowSums(prob)
    d$Score <- vapply(seq_len(nrow(d)), function(i) sample(0:2, 1, prob = prob[i, ]), integer(1))
  })
  for (model in c("RSM", "PCM")) {
    f <- fit_mfrm(d, "Person", c("Rater", "Task"), "Score", model = model,
      step_facet = if (model == "PCM") "Task" else NULL, quad_points = 61)
    scores <- mfrm_person_likelihood_scores(f)
    par <- f$opt$par
    for (person in c("P1", "P2", "P3")) {
      rows <- which(d$Person == person)
      loglik <- function(p) {
        density <- function(z) {
          eta <- outer(z, r[rows] * p[1] + t[rows] * p[2], "-")
          steps <- if (model == "RSM") rep(p[3], length(rows)) else p[2 + match(d$Task[rows], c("X", "Y"))]
          middle <- sweep(eta, 2, steps, "-")
          hi <- pmax(0, middle, 2 * eta)
          normalizer <- hi + log(exp(-hi) + exp(middle - hi) + exp(2 * eta - hi))
          numerator <- sweep(eta, 2, d$Score[rows], "*") -
            matrix(steps * (d$Score[rows] == 1), nrow(eta), length(rows), byrow = TRUE)
          exp(rowSums(numerator - normalizer)) * dnorm(z)
        }
        log(integrate(density, -Inf, Inf, rel.tol = 1e-10)$value)
      }
      reference <- vapply(seq_along(par), function(j) {
        delta <- numeric(length(par)); delta[j] <- 1e-4
        (loglik(par + delta) - loglik(par - delta)) / 2e-4
      }, numeric(1))
      expect_equal(unname(scores[person, ]), reference, tolerance = 1e-6)
    }
  }
})

test_that("both plotted covariance methods mark constraint-fixed targets", {
  x <- readRDS(test_path("fixtures", "anchored-facet-information.rds"))$RSM_group
  local_mocked_bindings(compute_mml_parameter_covariance = function(...)
    list(status = "ok", cov = x$covariance),
    mfrm_person_likelihood_scores = function(...) x$scores, .package = "mfrmr")
  ci <- mfrm_facet_intervals(x$fit, "Rater", method = "sandwich")
  shapes <- list()
  draw_points <- graphics::points
  local_mocked_bindings(points = function(...) {
    args <- list(...); shapes[[length(shapes) + 1L]] <<- args$pch
    draw_points(...)
  }, .package = "mfrmr")
  grDevices::pdf(NULL)
  tryCatch(plot(ci), finally = grDevices::dev.off())
  expect_gte(length(shapes), 2L)
  expect_identical(shapes[[1]], c(5, 1, 1, 5))
  expect_identical(shapes[[2]], c(5, 16, 16, 5))
})

test_that("information failure distinguishes anchored values from missing uncertainty", {
  x <- readRDS(test_path("fixtures", "anchored-facet-information.rds"))$RSM_group
  withr::local_options(mfrmr.max_information_bytes = 1)
  expect_warning(ci <- mfrm_facet_intervals(x$fit, "Rater"), "workspace estimate")
  fixed <- ci$table$Status == "fixed"
  expect_identical(fixed, c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(ci$table$SE[fixed], c(0, 0))
  expect_equal(ci$table$ModelSE[fixed], c(0, 0))
  expect_true(all(is.na(ci$table$Lower)))
  expect_true(all(ci$covariance[fixed, ] == 0))
  expect_true(all(is.na(ci$covariance[!fixed, !fixed])))
  expect_true(all(grepl("Constraint-fixed", ci$table$InferenceReview[fixed])))
  expect_true(all(grepl("workspace estimate", ci$table$InferenceReview[!fixed])))
  ci <- suppressWarnings(mfrm_facet_intervals(x$fit, "Rater", method = "sandwich"))
  shapes <- list()
  draw_points <- graphics::points
  local_mocked_bindings(points = function(...) {
    args <- list(...); shapes[[length(shapes) + 1L]] <<- args$pch
    draw_points(...)
  }, .package = "mfrmr")
  grDevices::pdf(NULL)
  tryCatch(plot(ci, show_legend = FALSE), finally = grDevices::dev.off())
  expect_identical(shapes, rep(list(c(5, 4, 4, 5)), 2L))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    built <- ggplot2::ggplot_build(as_ggplot(ci))
    expect_equal(nrow(built$data[[1]]), 0L)
    expect_identical(built$data[[2]]$shape, rep(c(5, 4, 4, 5), 2L))
  }
})

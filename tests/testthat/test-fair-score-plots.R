test_that("fair-score derivatives match the actual reference transformation", {
  for (model in c("RSM", "PCM")) {
    fit <- make_toy_fit(model = model, maxit = 20)
    # Nonzero person mean distinguishes FairM from FairZ on non-person rows.
    fit$facets$person$Estimate <- fit$facets$person$Estimate + 0.7
    dx <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "legacy")
    for (metric in c("FairM", "FairZ")) {
      p <- plot_fair_average(fit, diagnostics = dx, metric = metric, show_ci = TRUE, draw = FALSE)
      d <- p$data$data
      derivative <- mfrmr:::.fair_average_delta_variance(fit, d, metric)
      for (i in seq_len(nrow(d))) {
        plus <- minus <- dx
        at <- which(dx$measures$Facet == d$Facet[i] & dx$measures$Level == d$Level[i])
        if (!length(at)) next
        plus$measures$Estimate[at] <- plus$measures$Estimate[at] + 1e-5
        minus$measures$Estimate[at] <- minus$measures$Estimate[at] - 1e-5
        # Perturb only the internal transformation: public reports reject
        # diagnostic measures that no longer match the fitted analysis.
        a <- mfrmr:::calc_fair_average_bundle(fit, plus)$raw_by_facet[[d$Facet[i]]]
        b <- mfrmr:::calc_fair_average_bundle(fit, minus)$raw_by_facet[[d$Facet[i]]]
        numerical <- abs((a[[metric]][match(d$Level[i], a$Level)] - b[[metric]][match(d$Level[i], b$Level)]) / 2e-5)
        expect_equal(derivative[i], unname(numerical), tolerance = 1e-7)
      }
      scaled <- plot_fair_average(fit, diagnostics = dx, metric = metric, show_ci = TRUE,
        umean = 50, uscale = -10, draw = FALSE)$data$data
      at <- match(paste(d$Facet, d$Level), paste(scaled$Facet, scaled$Level))
      expect_equal(d$CI_SE, scaled$CI_SE[at], tolerance = 1e-8)
      expect_equal(d$CI_Lower, scaled$CI_Lower[at], tolerance = 1e-8)
    }
  }
})

test_that("PCM conditional intervals follow administered scales, weights and score selection", {
  fit <- readRDS(test_path("fixtures", "anchored-facet-information.rds"))$PCM$fit
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  original_dx <- dx
  pars <- mfrmr:::expand_params(fit$opt$par, mfrmr:::build_param_sizes(fit$config), fit$config)
  facet_means <- vapply(fit$config$facet_names, function(f) mean(pars$facets[[f]]), numeric(1))
  signs <- fit$config$facet_signs
  theta <- fit$facets$person$Estimate
  theta_mean <- mean(theta[is.finite(theta)])
  # Keep the saved calibration fixed while exercising unequal reporting weights
  # and the explicit totalscore exclusion rule with one extreme Person.
  dx$obs$Weight <- 1 + seq_len(nrow(dx$obs)) %% 7
  dx$obs$Observed[dx$obs$Person == dx$obs$Person[1]] <- fit$prep$rating_min
  facets <- c("Person", fit$config$facet_names)
  extreme <- vapply(facets, function(f) {
    vapply(as.character(dx$obs[[f]]), function(level) {
      y <- dx$obs$Observed[as.character(dx$obs[[f]]) == level]
      length(unique(y)) == 1L && y[1] %in% c(fit$prep$rating_min, fit$prep$rating_max)
    }, logical(1))
  }, logical(nrow(dx$obs)))
  k <- 0:ncol(pars$steps_mat)
  score <- function(eta, owner) {
    log_p <- k * eta - c(0, cumsum(pars$steps_mat[owner, ]))
    p <- exp(log_p - max(log_p)); p <- p / sum(p)
    sum(k * p) + fit$prep$rating_min
  }
  results <- list()
  for (totalscore in c(TRUE, FALSE)) {
    raw <- mfrmr:::stack_fair_raw_tables(mfrmr:::calc_facets_report_tbls(fit, dx,
      totalscore = totalscore))
    for (metric in c("FairM", "FairZ")) {
      actual <- mfrmr:::.fair_average_delta_variance(fit, raw, metric,
        totalscore = totalscore, obs = dx$obs)
      reference <- vapply(seq_len(nrow(raw)), function(i) {
        facet <- raw$Facet[i]
        use <- as.character(dx$obs[[facet]]) == raw$Level[i]
        if (!totalscore) use <- use & rowSums(extreme[, setdiff(facets, facet), drop = FALSE]) == 0
        owner <- match(as.character(dx$obs[[fit$config$step_facet]][use]),
          fit$prep$levels[[fit$config$step_facet]])
        w <- dx$obs$Weight[use]
        if (!length(w)) return(NA_real_)
        focal <- if (facet == "Person") raw$PrimaryMeasure[i] else signs[[facet]] * raw$PrimaryMeasure[i]
        eta <- focal + if (metric == "FairZ") 0 else {
          other <- setdiff(names(facet_means), facet)
          sum(signs[other] * facet_means[other]) + if (facet == "Person") 0 else theta_mean
        }
        average <- function(e) sum(w * vapply(owner, function(j) score(e, j), numeric(1))) / sum(w)
        expect_equal(raw[[metric]][i], average(eta), tolerance = 1e-12)
        (average(eta + 1e-5) - average(eta - 1e-5)) / 2e-5
      }, numeric(1))
      expect_equal(actual, reference, tolerance = 1e-7)
      results[[paste(totalscore, metric)]] <- actual
    }
  }
  expect_gt(max(abs(results[["TRUE FairM"]] - results[["FALSE FairM"]]), na.rm = TRUE), 1e-6)
  # Public plot dispatch must forward the same observation scope and reverse
  # reporting-scale changes once, while retaining diagnostic-only intervals.
  # Restore the actual saved analysis: controlled internal inputs above are
  # deliberately not acceptable to the public diagnostics-identity guard.
  dx <- original_dx
  for (metric in c("FairM", "FairZ")) {
    p <- plot(fit, type = "fair_average", diagnostics = dx, metric = metric,
      totalscore = FALSE, show_ci = TRUE, draw = FALSE)$data$data
    scaled <- plot(fit, type = "fair_average", diagnostics = dx, metric = metric,
      totalscore = FALSE, show_ci = TRUE, umean = 50, uscale = -10, draw = FALSE)$data$data
    at <- match(paste(p$Facet, p$Level), paste(scaled$Facet, scaled$Level))
    reference <- mfrmr:::.fair_average_delta_variance(fit, p, metric,
      totalscore = FALSE, obs = dx$obs)
    expect_equal(p$CI_SE, reference * p$ModelSE, tolerance = 1e-7)
    expect_equal(p$CI_SE, scaled$CI_SE[at], tolerance = 1e-8)
    expect_equal(p$CI_Lower, scaled$CI_Lower[at], tolerance = 1e-8)
    expect_true(all(!p$CI_Eligible))
  }
})

test_that("fair-score relationship plots preserve data, notes and grayscale encodings", {
  fit <- make_toy_fit(maxit = 10)
  bundle <- fair_average_table(fit)
  bundle$raw_by_facet$Person$ObservedAverage[1] <- NA_real_
  p <- plot_fair_average(bundle, plot_type = "measure", facet = "Person", draw = FALSE,
    metric = "FairZ", preset = "monochrome", show_title = FALSE, show_notes = FALSE)
  expect_equal(nrow(p$data$data), nrow(bundle$raw_by_facet$Person))
  expect_equal(p$data$plot_data$X, p$data$data$Measure)
  expect_equal(p$data$plot_data$Y, p$data$data$FairZ)
  expect_equal(nrow(p$data$reference_lines), 0L)
  expect_true(any(grepl("does not mean z-score", p$data$notes$Text)))
  all <- plot_fair_average(fit, plot_type = "scatter", draw = FALSE, preset = "monochrome")
  rgb <- grDevices::col2rgb(all$data$encoding$Color)
  expect_equal(unname(rgb[1, ]), unname(rgb[2, ]))
  expect_equal(unname(rgb[2, ]), unname(rgb[3, ]))
  expect_equal(length(unique(all$data$encoding$Shape)), nrow(all$data$encoding))
  ranked <- plot_fair_average(fit, top_n = 2, draw = FALSE)
  expect_equal(nrow(ranked$data$plot_data), 2L)
  expect_true(all(diff(abs(ranked$data$plot_data$Gap)) <= 0))
  expect_equal(nrow(ranked$data$data), nrow(all$data$data))
  unavailable <- plot_fair_average(fit, xtreme = 0.25, show_ci = TRUE, draw = FALSE)
  expect_true(all(is.na(unavailable$data$data$CI_Lower)))
  expect_true(all(unavailable$data$data$CI_Status == "unavailable"))
  expect_error(plot_fair_average(fit, top_n = NA, draw = FALSE), "top_n")
  skip_if_not_installed("ggplot2")
  g <- as_ggplot(p)
  expect_null(g$labels$title)
  expect_null(g$labels$subtitle)
  expect_identical(utils::head(plot_data(g, "notes"), nrow(p$data$notes)), p$data$notes)
  expect_no_warning(ggplot2::ggplot_build(g))
  built <- ggplot2::ggplot_build(g)$data[[1]]
  expect_equal(built$x, unname(p$data$plot_data$X))
  expect_equal(built$y, unname(p$data$plot_data$Y))
  expect_equal(nrow(built), nrow(p$data$plot_data))
})

test_that("stored GPCM fair-score intervals honor the requested confidence level", {
  # A precomputed bundle tests interval conversion without another costly fit.
  fit <- make_toy_fit(maxit = 10)
  b <- fair_average_table(fit)
  b$settings$model <- "GPCM"
  for (facet in names(b$raw_by_facet)) {
    d <- b$raw_by_facet[[facet]]
    for (metric in c("FairM", "FairZ")) {
      d[[paste0(metric, "SE")]] <- 0.1
      d[[paste0(metric, "_SE_Status")]] <- "ok"
      d[[paste0(metric, "_CI_Level")]] <- 0.95
      d[[paste0(metric, "_CI_Lower")]] <- d[[metric]] - qnorm(.975) * .1
      d[[paste0(metric, "_CI_Upper")]] <- d[[metric]] + qnorm(.975) * .1
    }
    b$raw_by_facet[[facet]] <- d
  }
  p <- plot_fair_average(b, show_ci = TRUE, ci_level = .9, draw = FALSE)
  expect_equal(p$data$data$CI_Lower, pmax(b$settings$rating_min, p$data$data$FairM - qnorm(.95)*.1))
  expect_true(all(p$data$data$CI_Level == .9))
  b$settings$rating_min <- b$settings$rating_max <- NULL
  p <- plot_fair_average(b, show_ci = TRUE, ci_level = .9, draw = FALSE)
  expect_true(all(is.na(p$data$data$CI_Lower)))
})

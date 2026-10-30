posterior_example <- function() read_mfrm_posterior(system.file(
  "extdata", "posterior-example", package = "mfrmr"))

test_that("saved posterior summaries connect to methods without a backend", {
  x <- posterior_example()
  expect_s3_class(x, "mfrm_posterior_summary")
  expect_identical(x$source$Estimation, "External Bayesian sampling")
  expect_identical(x$source$ExportSourceMethod, "MML")
  expect_identical(summary(x), x)
  before <- serialize(x, NULL)
  expect_output(print(x), "Sampling checks do not establish model validity", fixed = TRUE)
  ci <- confint(x)
  p <- x$parameters[x$parameters$role == "reported_parameter", ]
  expect_identical(ci[names(p)], p)
  fixed <- ci$value_role %in% c("fixed_anchor", "fixed_by_constraint")
  expect_true(any(fixed))
  expect_true(all(is.na(ci$lower[fixed]) & is.na(ci$upper[fixed])))
  expect_true(all(ci$interval_status[!fixed] == "computed_sampling_review_required"))
  expect_match(unique(ci$interval_type), "posterior credible", fixed = TRUE)
  expect_identical(confint(x, parm = "locations[1]")$variable, "locations[1]")
  expect_identical(confint(x, parm = 1L)$variable, x$parameters$variable[1L])
  expect_error(confint(x, level = .9), "backend's draws", fixed = TRUE)
  for (parm in list("unknown", c(1L, 1L), 0, 1.5, NA, logical()))
    expect_error(confint(x, parm = parm), "parm")
  expect_identical(serialize(x, NULL), before)
  expect_error(mfrm_results(x), "directly to mfrm_report", fixed = TRUE)
})

test_that("checksum, fixed-value and diagnostic inconsistencies fail before reporting", {
  dir <- tempfile(); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  src <- system.file("extdata", "posterior-example", package = "mfrmr")
  file.copy(file.path(src, c("model-data.rds", "posterior-summary.rds")), dir)
  expect_s3_class(read_mfrm_posterior(file.path(dir, "posterior-summary.rds")), "mfrm_posterior_summary")
  b <- readRDS(file.path(dir, "model-data.rds"))
  # Recompression preserves R values but breaks the saved file identity.
  saveRDS(b, file.path(dir, "model-data.rds"), compress = FALSE)
  expect_identical(readRDS(file.path(dir, "model-data.rds")), b)
  expect_error(read_mfrm_posterior(dir), "checksum")
  b$priors$location_marginal_sd <- 2
  saveRDS(b, file.path(dir, "model-data.rds"))
  expect_error(read_mfrm_posterior(dir), "checksum")
  expect_error(read_mfrm_posterior(tempfile()), "together")
  x <- posterior_example()
  mutations <- list(
    function(z) {z$parameters$level[1] <- "another person"; z},
    function(z) {z$parameters$mean[z$parameters$value_role == "fixed_anchor"][1] <- 5; z},
    function(z) {z$parameters$lower[z$parameters$value_role == "fixed_anchor"][1] <- 0; z},
    function(z) {z$decision$Status <- "no_sampling_flags_detected"; z},
    function(z) {z$checks$Status[] <- "no_flag"; z},
    function(z) {z$sampler <- z$sampler[-1, ]; z},
    function(z) {z$parameters$variable[2] <- z$parameters$variable[1]; z},
    function(z) {z$parameters$lower[z$parameters$quantity == "slope"][1] <- -.1; z},
    function(z) {z$parameters$interval_status[] <- "computed_model_review_required"; z})
  for (change in mutations) {
    bad <- change(x)
    expect_error(summary(bad), "Inconsistent")
    expect_error(confint(bad), "Inconsistent")
    expect_error(plot(bad, draw = FALSE), "Inconsistent")
    expect_error(mfrm_report(bad), "Inconsistent")
  }
})

test_that("reports, APA tables and appendix exports preserve posterior evidence", {
  x <- posterior_example()
  context <- c(score_direction = "Higher means harder to understand")
  r <- mfrm_report(x, style = "apa", study_context = context)
  expect_identical(r$tables$parameters, x$parameters)
  expect_identical(r$tables$sampling_review, x$checks)
  expect_identical(r$tables$source_model$ModelDataMD5, x$model_data_md5)
  expect_match(r$first_screen$NextAction, "sampling_review", fixed = TRUE)
  expect_identical(r$tables$study_context$Value, unname(context))
  for (obj in list(x, summary(x), r, summary(r))) {
    b <- build_summary_table_bundle(obj, top_n = 1)
    expect_identical(b$tables$parameters, x$parameters)
  }
  expect_s3_class(apa_table(x, which = "parameters"), "apa_table")
  h <- mfrm_report(x, output = "html", study_context = context)
  on.exit(unlink(h$path), add = TRUE)
  text <- paste(readLines(h$path, warn = FALSE), collapse = "\n")
  expect_match(text, "posterior credible", fixed = TRUE)
  expect_match(text, "fixed_anchor", fixed = TRUE)
  expect_match(text, unname(context), fixed = TRUE)
  dir <- tempfile(); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  export_summary_appendix(r, output_dir = dir, prefix = "posterior")
  file <- list.files(dir, "parameters.*csv$", full.names = TRUE)
  expect_length(file, 1L)
  p <- read.csv(file)
  expect_equal(p$lower, x$parameters$lower, tolerance = 1e-12)
  expect_identical(p$interval_status, x$parameters$interval_status)
  saveRDS(x, file.path(dir, "portable.rds"))
  expect_identical(summary(readRDS(file.path(dir, "portable.rds"))), x)
  expect_error(mfrm_report(x, style = "rater"), "rater")
  expect_error(mfrm_report(x, intervals = confint(x)), "only to G/D")
})

test_that("posterior plots preserve units, fixed values and exact interval endpoints", {
  x <- posterior_example()
  for (quantity in c("location", "ability", "step_offset", "slope", "log_slope")) {
    p <- plot(x, quantity = quantity, draw = FALSE)
    tab <- plot_data(p, "table")
    expect_identical(tab[names(x$parameters)], x$parameters[x$parameters$quantity == quantity, ])
    expect_null(p$data$title)
    expect_identical(p$data$notes, x$notes)
    expect_true(all(tab$Display[tab$value_role == "fixed_anchor"] == "Fixed value"))
  }
  expect_error(plot(x, facet = "unknown", draw = FALSE), "facet")
  expect_error(plot(x, reference = Inf, draw = FALSE), "reference")
  expect_error(plot(x, draw = NA), "draw")
  expect_error(plot(x, title = NA_character_, draw = FALSE), "text")
  p <- plot(x, draw = FALSE)
  skip_if_not_installed("ggplot2")
  gg <- as_ggplot(p)
  expect_identical(attr(gg, "mfrmr_plot_data"), p)
  built <- ggplot2::ggplot_build(gg)
  finite <- p$data$table[is.finite(p$data$table$lower), ]
  expect_equal(built$data[[1]]$x, finite$lower)
  expect_equal(built$data[[1]]$xend, finite$upper)
  expect_true(all(built$data[[1]]$linetype == "dashed"))
  expect_true(5 %in% built$data[[2]]$shape)
  expect_error(as_ggplot(p, component = "table"), "plot_data")
})

test_that("RSM, PCM and product-GPCM summary formats preserve model roles", {
  fixtures <- readRDS(test_path("fixtures", "anchored-facet-information.rds"))
  fits <- list(fixtures$RSM$fit, fixtures$PCM$fit,
    readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))$fit)
  helpers <- mfrmr:::.mfrm_posterior_helpers()
  for (fit in fits) {
    bundle <- mfrmr:::mfrm_model_export_data(fit, 1, 1, 1)
    p <- helpers$mfrm_parameter_map(bundle)
    fixed <- !is.na(p$fixed_value)
    # Fabricated summaries exercise the reader contract, not posterior accuracy.
    p$mean <- p$median <- ifelse(fixed, p$fixed_value, .2)
    p$sd <- ifelse(fixed, NA_real_, .3)
    p$lower <- ifelse(fixed, NA_real_, -.4)
    p$upper <- ifelse(fixed, NA_real_, .8)
    p$lower[p$quantity == "slope" & !fixed] <- .1
    p$rhat <- ifelse(fixed, NA_real_, 1)
    p$ess_bulk <- p$ess_tail <- ifelse(fixed, NA_real_, 800)
    p$mcse_lower <- p$mcse_upper <- ifelse(fixed, NA_real_, .01)
    p$role <- "reported_parameter"; p$interval_computed <- !fixed
    p$sampling_review <- "no_sampling_flags_detected"
    p$interval_status <- ifelse(fixed, "fixed_value_no_interval", "computed_model_review_required")
    sampler <- data.frame(chain = 1:4, num_divergent = 0, num_max_treedepth = 0, ebfmi = .8)
    dir <- tempfile(); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    saveRDS(bundle, file.path(dir, "model-data.rds"))
    x <- list(schema_version = 1L, chains = 4L,
      model_data_md5 = unname(tools::md5sum(file.path(dir, "model-data.rds"))),
      parameters = p, sampler = sampler, checks = helpers$mfrm_sampling_review(p, sampler, 4L),
      decision = data.frame(Status = "no_sampling_flags_detected", Interpretation = "Model review remains",
        Why = "Fabricated numerical summaries", NextAction = "Do not use this test as sampling evidence"),
      backend = "Stan", interval = "equal-tailed posterior credible interval for varying parameters",
      level = .95, notes = "Fabricated format test; no sampling evidence")
    saveRDS(x, file.path(dir, "posterior-summary.rds"))
    result <- read_mfrm_posterior(dir)
    expect_identical(result$source$Model, fit$config$model)
    expect_identical(confint(result)$lower, p$lower)
    report <- mfrm_report(result)
    expect_identical(report$first_screen$Status, "review")
    expect_match(report$first_screen$NextAction, "priors", fixed = TRUE)
    if (fit$config$model == "GPCM")
      expect_equal(length(unique(plot_data(plot(result, quantity = "slope", draw = FALSE), "table")$facet)), 2L)
    else expect_error(plot(result, quantity = "slope", draw = FALSE), "No saved parameters")
  }
})

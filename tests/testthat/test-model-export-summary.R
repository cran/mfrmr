export_summary_helpers <- function() {
  env <- new.env(parent = globalenv())
  sys.source(system.file("model-export", "summarize.R", package = "mfrmr"), env)
  env
}

test_that("posterior labels and fixed status come from exported constraints", {
  env <- export_summary_helpers()
  fit <- readRDS(test_path("fixtures", "anchored-gpcm.rds"))$fixed_population$fit
  bundle <- mfrmr:::mfrm_model_export_data(fit, 1, 1, 1)
  map <- env$mfrm_parameter_map(bundle)
  locations <- map[match(bundle$levels$Parameter, map$variable), ]
  expect_identical(locations$level, bundle$levels$Level)
  expect_equal(!is.na(locations$fixed_value), bundle$levels$Fixed)
  expect_true(any(locations$value_role == "fixed_anchor"))
  expect_true(all(locations$value_role[locations$quantity == "ability"] == "estimated"))
  expect_match(locations$direction[locations$quantity == "ability"][1], "higher coded scores", fixed = TRUE)
  expect_true(all(grepl("lower coded scores", locations$direction[locations$quantity == "location"], fixed = TRUE)))
  expect_false(anyDuplicated(map$variable) > 0L)

  # A singleton slope and a binary centered step are identified constants, not anchors.
  bundle$stan_data$slope_basis <- matrix(numeric(), 1L, 0L)
  bundle$slope_levels <- bundle$slope_levels[1L, ]
  bundle$stan_data$C <- 2L
  bundle$stan_data$step_basis <- matrix(numeric(), bundle$stan_data$G, 0L)
  bundle$stan_data$step_offset <- rep(0, bundle$stan_data$G)
  bundle$constraints$steps <- NULL
  singleton <- env$mfrm_parameter_map(bundle)
  expect_true(all(singleton$value_role[singleton$quantity %in% c("slope", "step_offset")] == "fixed_by_constraint"))
  expect_equal(singleton$fixed_value[singleton$quantity == "slope"], 1)
  expect_equal(singleton$fixed_value[singleton$quantity == "log_slope"], 0)
})

test_that("missing or poor diagnostics cannot become a clean sampling screen", {
  env <- export_summary_helpers()
  p <- data.frame(value_role = c("estimated", "fixed_anchor", "sampling_coordinate"),
    rhat = c(1, NA, 1), ess_bulk = c(800, NA, 700), ess_tail = c(700, NA, 600))
  sampler <- data.frame(chain = 1:4, num_divergent = 0, num_max_treedepth = 0, ebfmi = .8)
  clean <- env$mfrm_sampling_review(p, sampler, 4L)
  expect_true(all(clean$Status == "no_flag"))
  p$rhat[3] <- NA_real_
  missing <- env$mfrm_sampling_review(p, sampler, 4L)
  expect_identical(missing$Status[missing$Check == "rhat"], "not_available")
  p$rhat[3] <- 1.02
  sampler$num_divergent[2] <- 1
  sampler$ebfmi[3] <- .2
  bad <- env$mfrm_sampling_review(p, sampler, 4L)
  expect_true(all(bad$Status[bad$Check %in% c("rhat", "num_divergent", "ebfmi")] == "review"))
  incomplete <- env$mfrm_sampling_review(p, sampler[-1L, ], 4L)
  expect_true(all(incomplete$Status[incomplete$Check %in% c("num_divergent", "num_max_treedepth", "ebfmi")] == "not_available"))
  sampler$chain[2] <- 1
  duplicate <- env$mfrm_sampling_review(p, sampler, 4L)
  expect_true(all(duplicate$Status[duplicate$Check %in% c("num_divergent", "num_max_treedepth", "ebfmi")] == "not_available"))
})

test_that("legacy unanchored maps retain binary steps and both slope owners", {
  env <- export_summary_helpers()
  fit <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))$fit
  bundle <- mfrmr:::mfrm_model_export_data(fit, 1, 1, 1)
  expect_false(isTRUE(bundle$affine))
  map <- env$mfrm_parameter_map(bundle)
  slopes <- map[map$quantity == "slope", ]
  expect_identical(slopes$level, bundle$slope_levels$Level)
  expect_identical(slopes$facet, bundle$slope_levels$Owner)
  expect_true(all(slopes$value_role == "estimated"))
  expect_true(any(grepl("geometric mean slope = 1", slopes$reference, fixed = TRUE)))
  expect_true(any(grepl("Component scale identified", slopes$reference, fixed = TRUE)))
  # RSM/PCM without affine anchors use the same shared step basis.
  bundle$slope_levels <- NULL
  bundle$stan_data$C <- 2L
  bundle$stan_data$step_basis <- matrix(numeric(), 1L, 0L)
  binary <- env$mfrm_parameter_map(bundle)
  expect_true(all(binary$value_role[binary$quantity == "step_offset"] == "fixed_by_constraint"))
  expect_true(all(binary$fixed_value[binary$quantity == "step_offset"] == 0))
})

test_that("both bundle runners can be sourced without executing estimation", {
  for (format in c("stan", "brms")) {
    env <- new.env(parent = globalenv())
    expect_silent(sys.source(system.file("model-export", paste0("run-", format, ".R"), package = "mfrmr"), env))
    expect_true(is.function(env[[paste0("summarize_mfrm_", format)]]))
    expect_error(env[[paste0("run_mfrm_", format)]]("--unknown"), "Usage:")
  }
})

test_that("estimated-population summaries retain units and reach all reporting consumers", {
  skip_if_not_installed("posterior")
  withr::local_seed(27)
  fit <- readRDS(test_path("fixtures", "anchored-gpcm.rds"))$shared$fit
  dir <- tempfile(); dir.create(dir)
  withr::defer(unlink(dir, recursive = TRUE))
  export_mfrm_model(fit, dir, acknowledge_sensitive = TRUE)
  bundle <- readRDS(file.path(dir, "model-data.rds"))
  env <- export_summary_helpers()
  map <- env$mfrm_parameter_map(bundle)
  expect_equal(map$quantity[map$facet == "Population"], c("population_mean", "population_sd"))
  # Synthetic draws exercise reporting, not posterior validity or MCMC performance.
  vars <- c(map$variable, "population_log_sd", "theta_raw[1]")
  draws <- array(rnorm(200L * 4L * length(vars)), c(200L, 4L, length(vars)),
    dimnames = list(NULL, NULL, vars))
  for (i in which(map$quantity %in% c("slope", "population_sd"))) draws[, , i] <- exp(draws[, , i] * .2)
  for (i in which(!is.na(map$fixed_value))) draws[, , i] <- map$fixed_value[i]
  sampler <- data.frame(chain = 1:4, num_divergent = 0, num_max_treedepth = 0, ebfmi = .8)
  withr::local_dir(dir)
  env$summarize_mfrm_draws(draws, bundle, sampler, "Stan")
  x <- read_mfrm_posterior(dir)
  expect_identical(x$schema_version, 2L)
  expect_identical(summary(x), x)
  expect_identical(x$source$Population, "Normal population with estimated mean and SD")
  p <- x$parameters
  expect_identical(p$role[p$variable == "theta_raw[1]"], "free_coordinate")
  expect_identical(p$value_role[p$quantity == "population_sd"], "estimated")
  for (quantity in c("population_mean", "population_sd", "ability", "location")) {
    plot <- plot(x, quantity = quantity, draw = FALSE)
    expect_identical(plot_data(plot, "table")[names(p)], p[p$quantity == quantity, ])
    expect_false(grepl("fixed N(0,1)", plot$data$xlab, fixed = TRUE))
  }
  expect_equal(confint(x, parm = "population_sd")$lower, p$lower[p$variable == "population_sd"])
  report <- mfrm_report(x)
  expect_identical(report$tables$parameters, p)
  expect_identical(report$tables$source_model$Population, x$source$Population)
  export_summary_appendix(report, output_dir = file.path(dir, "appendix"))
  files <- list.files(file.path(dir, "appendix"), "parameters.*csv$", full.names = TRUE)
  expect_length(files, 1L)
  csv <- read.csv(files)
  expect_equal(csv$lower[csv$variable == "population_sd"], p$lower[p$variable == "population_sd"])
  bad <- x; bad$parameters$lower[bad$parameters$quantity == "population_sd"] <- -1
  expect_error(summary(bad), "Inconsistent")
  # A badly mixed internal population coordinate must reach every consumer,
  # even when the reported SD draws themselves have no obvious chain shift.
  draws[, 4L, "population_log_sd"] <- draws[, 4L, "population_log_sd"] + 3
  env$summarize_mfrm_draws(draws, bundle, sampler, "Stan")
  flagged <- read_mfrm_posterior(dir)
  expect_identical(flagged$checks$Status[flagged$checks$Check == "rhat"], "review")
  expect_identical(flagged$decision$Status, "sampling_review_required")
  expect_true(all(confint(flagged, parm = "population_sd")$interval_status ==
    "computed_sampling_review_required"))
  expect_identical(mfrm_report(flagged)$tables$sampling_review, flagged$checks)
  expect_identical(plot_data(plot(flagged, quantity = "population_sd", draw = FALSE),
    "table")$sampling_review, "sampling_review_required")
  saved <- readRDS(file.path(dir, "posterior-summary.rds"))
  saved$schema_version <- 1L
  saveRDS(saved, file.path(dir, "posterior-summary.rds"))
  expect_error(read_mfrm_posterior(dir), "Incompatible population")
})

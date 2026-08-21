gpcm_identification_fixture <- function() {
  data <- expand.grid(
    Person = sprintf("P%02d", seq_len(36L)),
    Item = paste0("I", seq_len(4L)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  set.seed(20260811)
  theta <- stats::rnorm(36L, mean = 0.25, sd = 0.8)
  location <- c(-0.8, -0.2, 0.3, 0.7)
  slope <- c(0.65, 0.9, 1.2, 1.55)
  step <- c(-0.9, 0, 0.9)
  person_index <- match(data$Person, sprintf("P%02d", seq_len(36L)))
  item_index <- match(data$Item, paste0("I", seq_len(4L)))
  probability <- t(vapply(seq_len(nrow(data)), function(row) {
    log_kernel <- c(0, cumsum(
      slope[item_index[row]] *
        (theta[person_index[row]] - location[item_index[row]] - step)
    ))
    value <- exp(log_kernel - max(log_kernel))
    value / sum(value)
  }, numeric(4L)))
  data$Score <- vapply(seq_len(nrow(data)), function(row) {
    sample.int(4L, size = 1L, prob = probability[row, ]) - 1L
  }, integer(1L))
  data
}

fit_gpcm_identification_fixture <- function(data, ...) {
  suppressMessages(suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = "Item",
    score = "Score",
    model = "GPCM",
    method = "MML",
    step_facet = "Item",
    slope_facet = "Item",
    quad_points = 7L,
    maxit = 30L,
    ...
  )))
}

test_that("default GPCM MML estimates the common population scale", {
  fit <- fit_gpcm_identification_fixture(gpcm_identification_fixture())

  expect_true(isTRUE(fit$population$active))
  expect_identical(fit$population$source, "gpcm_mml_default_identification")
  expect_identical(
    fit$population$identification_role,
    "free_population_scale_with_geometric_mean_one_relative_slopes"
  )
  expect_identical(fit$config$gpcm_mml_identification, "free_population")
  expect_identical(
    fit$config$gpcm_estimator_family,
    "marginal_maximum_likelihood"
  )
  expect_identical(fit$config$gpcm_statistical_penalty, "none")
  expect_false(fit$config$gpcm_finite_parameter_box)
  expect_identical(
    fit$config$gpcm_model_family,
    "aligned_single_owner_relative_slope_gpcm"
  )
  expect_identical(
    fit$config$gpcm_slope_action,
    "complete_adjacent_predictor"
  )
  expect_identical(
    fit$config$gpcm_slope_composition,
    "single_owner_relative_gm1"
  )
  expect_identical(fit$config$gpcm_latent_dimension_count, 1L)
  expect_identical(
    fit$config$gpcm_extreme_person_policy,
    "posterior_eap_under_population_model"
  )
  expect_identical(
    fit$config$gpcm_common_discrimination,
    "estimated_via_population_sd"
  )
  expect_identical(paste(deparse(fit$population$formula), collapse = " "), "~1")
  expect_true(is.finite(fit$population$sigma2))
  expect_gt(fit$population$sigma2, 0)

  relative <- as.numeric(fit$slopes$OptimizerEstimate)
  fixed_sd <- as.numeric(fit$slopes$FixedLatentSDOptimizerEstimate)
  population_sd <- sqrt(fit$population$sigma2)
  expect_equal(exp(mean(log(relative))), 1, tolerance = 1e-8)
  expect_equal(fixed_sd, population_sd * relative, tolerance = 1e-12)
  expect_equal(exp(mean(log(fixed_sd))), population_sd, tolerance = 1e-8)

  summary_fit <- summary(fit)
  expect_identical(
    summary_fit$settings_overview$GpcmMmlIdentification[1],
    "free_population"
  )
  expect_identical(
    summary_fit$settings_overview$GpcmEstimatorFamily[1],
    "marginal_maximum_likelihood"
  )
  expect_identical(
    summary_fit$settings_overview$GpcmStatisticalPenalty[1],
    "none"
  )
  expect_identical(
    summary_fit$settings_overview$GpcmModelFamily[1],
    "aligned_single_owner_relative_slope_gpcm"
  )
  expect_identical(
    summary_fit$settings_overview$GpcmSlopeAction[1],
    "complete_adjacent_predictor"
  )
  expect_identical(
    summary_fit$settings_overview$GpcmSlopeComposition[1],
    "single_owner_relative_gm1"
  )
  expect_identical(
    summary_fit$settings_overview$GpcmLatentDimensionCount[1],
    1L
  )
  expect_false(summary_fit$settings_overview$GpcmFiniteParameterBox[1])
  expect_identical(
    summary_fit$settings_overview$GpcmExtremePersonPolicy[1],
    "posterior_eap_under_population_model"
  )
  expect_equal(
    summary_fit$slope_overview$FixedLatentSDOptimizerGeometricMean[1],
    population_sd,
    tolerance = 1e-8
  )
  fit_console <- capture.output(print(fit))
  fit_console_text <- gsub(
    "[[:space:]]+",
    " ",
    paste(fit_console, collapse = " ")
  )
  expect_true(any(grepl(
    "Scale: estimated_population_scale",
    fit_console_text,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Discrimination: geometric_mean_one_relative_discrimination",
    fit_console_text,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "GPCM estimator: marginal_maximum_likelihood",
    fit_console_text,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Statistical penalty: none | Finite parameter box: no",
    fit_console_text,
    fixed = TRUE
  )))
  kernel_console_line <- paste0(
    "GPCM kernel: aligned_single_owner_relative_slope_gpcm | ",
    "Slope action: complete_adjacent_predictor"
  )
  expect_true(any(grepl(kernel_console_line, fit_console_text, fixed = TRUE)))
  summary_console <- capture.output(print(summary_fit))
  summary_console_text <- gsub(
    "[[:space:]]+",
    " ",
    paste(summary_console, collapse = " ")
  )
  expect_true(any(grepl(
    "GPCM estimator: marginal_maximum_likelihood",
    summary_console_text,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    kernel_console_line,
    summary_console_text,
    fixed = TRUE
  )))
  for (plot_type in c("wright", "pathway", "ccc")) {
    plot_payload <- suppressWarnings(plot(fit, type = plot_type, draw = FALSE))
    expect_s3_class(plot_payload, "mfrm_plot_data")
    expect_identical(
      plot_payload$data$scale_contract$GpcmMmlIdentification[1],
      "free_population"
    )
    expect_equal(
      plot_payload$data$scale_contract$PopulationSD[1],
      population_sd,
      tolerance = 1e-8
    )
    expect_identical(
      plot_payload$data$scale_contract$SlopeBasis[1],
      "geometric_mean_one_relative_discrimination"
    )
    expect_identical(
      plot_payload$data$scale_contract$GpcmEstimatorFamily[1],
      "marginal_maximum_likelihood"
    )
    expect_identical(
      plot_payload$data$scale_contract$GpcmStatisticalPenalty[1],
      "none"
    )
    expect_identical(
      plot_payload$data$scale_contract$GpcmModelFamily[1],
      "aligned_single_owner_relative_slope_gpcm"
    )
    expect_identical(
      plot_payload$data$scale_contract$GpcmSlopeAction[1],
      "complete_adjacent_predictor"
    )
    expect_identical(
      plot_payload$data$scale_contract$GpcmSlopeComposition[1],
      "single_owner_relative_gm1"
    )
    expect_identical(
      plot_payload$data$scale_contract$GpcmLatentDimensionCount[1],
      1L
    )
    expect_false(
      plot_payload$data$scale_contract$GpcmFiniteParameterBox[1]
    )
    expect_identical(
      plot_payload$data$scale_contract$GpcmExtremePersonPolicy[1],
      "posterior_eap_under_population_model"
    )
  }
  legacy_fit <- fit
  legacy_fit$config$gpcm_estimator_family <- NULL
  legacy_fit$config$gpcm_statistical_penalty <- NULL
  legacy_fit$config$gpcm_finite_parameter_box <- NULL
  legacy_fit$config$gpcm_extreme_person_policy <- NULL
  legacy_fit$config$gpcm_model_family <- NULL
  legacy_fit$config$gpcm_slope_action <- NULL
  legacy_fit$config$gpcm_slope_composition <- NULL
  legacy_fit$config$gpcm_latent_dimension_count <- NULL
  legacy_contract <- mfrmr:::mfrm_fit_scale_contract(legacy_fit)
  expect_identical(
    legacy_contract$GpcmEstimatorFamily[1],
    "marginal_maximum_likelihood"
  )
  expect_identical(legacy_contract$GpcmStatisticalPenalty[1], "none")
  expect_false(legacy_contract$GpcmFiniteParameterBox[1])
  expect_identical(
    legacy_contract$GpcmExtremePersonPolicy[1],
    "posterior_eap_under_population_model"
  )
  expect_identical(
    legacy_contract$GpcmModelFamily[1],
    "aligned_single_owner_relative_slope_gpcm"
  )
  expect_identical(
    legacy_contract$GpcmSlopeAction[1],
    "complete_adjacent_predictor"
  )
  expect_identical(
    mfrmr:::resolve_dff_refit_controls(fit)$gpcm_mml_identification,
    "free_population"
  )

  replay <- build_mfrm_replay_script(fit)
  expect_match(
    replay$script,
    'gpcm_mml_identification = "free_population"',
    fixed = TRUE
  )
  manifest <- build_mfrm_manifest(fit)
  manifest_value <- function(setting) {
    as.character(manifest$model_settings$Value[
      match(setting, manifest$model_settings$Setting)
    ])
  }
  expect_identical(
    manifest_value("gpcm_mml_identification"),
    "free_population"
  )
  expect_identical(
    manifest_value("gpcm_estimator_family"),
    "marginal_maximum_likelihood"
  )
  expect_identical(manifest_value("gpcm_statistical_penalty"), "none")
  expect_identical(
    manifest_value("gpcm_model_family"),
    "aligned_single_owner_relative_slope_gpcm"
  )
  expect_identical(
    manifest_value("gpcm_slope_action"),
    "complete_adjacent_predictor"
  )
  expect_identical(
    manifest_value("gpcm_slope_composition"),
    "single_owner_relative_gm1"
  )
  expect_identical(manifest_value("gpcm_latent_dimension_count"), "1")
  expect_identical(manifest_value("gpcm_finite_parameter_box"), "FALSE")
  expect_identical(
    manifest_value("gpcm_extreme_person_policy"),
    "posterior_eap_under_population_model"
  )
  expect_identical(
    manifest_value("population_source"),
    "gpcm_mml_default_identification"
  )

  sizes <- mfrmr:::build_param_sizes(fit$config)
  slices <- mfrmr:::build_param_slices(sizes)
  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = fit$config$step_facet,
    slope_facet = fit$config$slope_facet,
    interaction_specs = fit$config$interaction_specs
  )
  quad <- mfrmr:::gauss_hermite_normal(7L)
  objective <- function(par) {
    mfrmr:::mfrm_loglik_mml(par, idx, fit$config, sizes, quad)
  }
  probe <- as.numeric(fit$opt$par) +
    0.025 * sin(seq_along(fit$opt$par)) * pmax(1, abs(fit$opt$par))
  analytic <- mfrmr:::mfrm_grad_mml(probe, idx, fit$config, sizes, quad)
  step <- 1e-5 * pmax(1, abs(probe))
  numeric <- vapply(seq_along(probe), function(index) {
    high <- low <- probe
    high[index] <- high[index] + step[index]
    low[index] <- low[index] - step[index]
    (objective(high) - objective(low)) / (2 * step[index])
  }, numeric(1L))
  expect_lt(max(abs(analytic - numeric)), 1e-5)
  expect_lt(max(abs(
    analytic[c(slices$beta, slices$log_sigma2)] -
      numeric[c(slices$beta, slices$log_sigma2)]
  )), 1e-6)
})

test_that("automatic intercept-only population wiring matches an explicit model", {
  data <- gpcm_identification_fixture()
  person_data <- unique(data["Person"])
  automatic <- fit_gpcm_identification_fixture(data)
  explicit <- fit_gpcm_identification_fixture(
    data,
    population_formula = ~ 1,
    person_data = person_data,
    person_id = "Person"
  )

  expect_equal(automatic$opt$value, explicit$opt$value, tolerance = 1e-10)
  expect_equal(automatic$opt$par, explicit$opt$par, tolerance = 1e-10)
  expect_identical(explicit$population$source, "user_supplied")
  expect_identical(
    explicit$config$gpcm_mml_identification,
    "free_population"
  )
})

test_that("legacy fixed-standard-normal GPCM remains an explicit restriction", {
  data <- gpcm_identification_fixture()
  conventional <- fit_gpcm_identification_fixture(data)
  legacy <- fit_gpcm_identification_fixture(
    data,
    gpcm_mml_identification = "fixed_standard_normal"
  )

  expect_false(isTRUE(legacy$population$active))
  expect_identical(legacy$config$posterior_basis, "legacy_mml")
  expect_identical(
    legacy$config$gpcm_common_discrimination,
    "fixed_to_one_legacy_restriction"
  )
  expect_equal(unique(legacy$slopes$PopulationSD), 1)
  expect_equal(
    legacy$slopes$FixedLatentSDOptimizerEstimate,
    legacy$slopes$OptimizerEstimate,
    tolerance = 1e-12
  )
  expect_equal(length(conventional$opt$par), length(legacy$opt$par) + 2L)
  expect_identical(
    mfrmr:::resolve_dff_refit_controls(legacy)$gpcm_mml_identification,
    "fixed_standard_normal"
  )
  expect_match(
    build_mfrm_replay_script(legacy)$script,
    'gpcm_mml_identification = "fixed_standard_normal"',
    fixed = TRUE
  )

  person_data <- unique(data["Person"])
  expect_error(
    fit_gpcm_identification_fixture(
      data,
      gpcm_mml_identification = "fixed_standard_normal",
      population_formula = ~ 1,
      person_data = person_data
    ),
    "requires `population_formula = NULL`",
    fixed = TRUE
  )
})

test_that("automatic GPCM population IDs follow ordinary preparation", {
  data <- gpcm_identification_fixture()
  data$Person[1] <- NA_character_
  data$Person[2] <- paste0(" ", data$Person[2], " ")

  fit <- fit_gpcm_identification_fixture(data)

  expect_equal(fit$prep$row_retention$Rows[2], nrow(data) - 1L)
  expect_true(any(
    fit$prep$preparation_notes$Condition == "trimmed_person_ids"
  ))
  expect_false(anyNA(fit$config$population_spec$person_lookup))
  expect_equal(nrow(fit$population$person_table), 36L)
  expect_equal(fit$population$response_rows_retained, nrow(data) - 1L)
  expect_equal(fit$population$response_rows_omitted, 1L)
})

test_that("default GPCM location identification fails closed when facets are uncentered", {
  expect_error(
    fit_gpcm_identification_fixture(
      gpcm_identification_fixture(),
      noncenter_facet = "Item"
    ),
    "requires `noncenter_facet = \"Person\"`",
    fixed = TRUE
  )
})

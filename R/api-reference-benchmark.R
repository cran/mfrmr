# Package-native reference benchmark helpers

reference_benchmark_dataset_specs <- function() {
  tibble::tibble(
    Dataset = c(
      "example_core", "example_bias",
      "study1", "study2", "combined",
      "study1_itercal", "study2_itercal", "combined_itercal",
      "synthetic_truth", "synthetic_latent_regression", "synthetic_latent_regression_omit", "synthetic_gpcm"
    ),
    Persons = c(48L, 48L, 307L, 206L, 307L, 307L, 206L, 307L, 36L, 60L, 60L, 600L),
    Raters = c(4L, 4L, 18L, 12L, 18L, 18L, 12L, 18L, 3L, NA_integer_, NA_integer_, NA_integer_),
    Criteria = c(4L, 4L, 3L, 9L, 12L, 3L, 9L, 12L, 3L, 6L, 6L, 4L),
    Tasks = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, 4L, NA_integer_, NA_integer_, NA_integer_),
    Rows = c(768L, 384L, 1842L, 3287L, 5129L, 1842L, 3341L, 5183L, 1296L, 360L, 360L, 2400L),
    ScoreMin = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L),
    ScoreMax = c(4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 1L, 1L, 3L)
  )
}

reference_benchmark_case_specs <- function() {
  tibble::tibble(
    Case = c(
      "synthetic_truth",
      "synthetic_latent_regression",
      "synthetic_latent_regression_omit",
      "synthetic_conquest_overlap_dry_run",
      "synthetic_gpcm",
      "synthetic_bias_contract",
      "study1_itercal_pair",
      "study2_itercal_pair",
      "combined_itercal_pair"
    ),
    CaseType = c(
      "truth_recovery",
      "latent_recovery",
      "latent_omission_contract",
      "conquest_overlap_dry_run",
      "gpcm_recovery",
      "bias_contract",
      "pair_stability",
      "pair_stability",
      "pair_stability"
    ),
    PrimaryDataset = c(
      "synthetic_truth",
      "synthetic_latent_regression",
      "synthetic_latent_regression_omit",
      "synthetic_latent_regression",
      "synthetic_gpcm",
      "synthetic_truth",
      "study1",
      "study2",
      "combined"
    ),
    ReferenceDataset = c(
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_,
      "study1_itercal",
      "study2_itercal",
      "combined_itercal"
    )
  )
}

synthetic_truth_targets <- function() {
  list(
    Rater = c(-0.4, 0, 0.4),
    Task = seq(-0.5, 0.5, length.out = 4),
    Criterion = c(-0.3, 0, 0.3)
  )
}

sample_mfrm_latent_benchmark_data <- function(seed = 20240403) {
  with_preserved_rng_seed(seed, {
    persons <- paste0("P", sprintf("%03d", seq_len(60L)))
    criteria <- paste0("C", seq_len(6L))
    x <- seq(-1.5, 1.5, length.out = length(persons))
    beta <- c(`(Intercept)` = 0.25, X = 0.80)
    sigma2 <- 0.36
    theta <- beta[["(Intercept)"]] + beta[["X"]] * x + stats::rnorm(length(persons), sd = sqrt(sigma2))
    criterion_eff <- seq(-1.0, 1.0, length.out = length(criteria))

    dat <- expand.grid(Person = persons, Criterion = criteria, stringsAsFactors = FALSE)
    eta <- theta[match(dat$Person, persons)] - criterion_eff[match(dat$Criterion, criteria)]
    dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))

    list(
      data = dat,
      person = "Person",
      facets = "Criterion",
      score = "Score",
      person_data = data.frame(Person = persons, X = x, stringsAsFactors = FALSE),
      person_id = "Person",
      population_formula = ~ X,
      truth = list(
        coefficients = beta,
        sigma2 = sigma2,
        criterion = stats::setNames(criterion_eff, criteria)
      )
    )
  })
}

sample_mfrm_latent_omit_benchmark_data <- function(seed = 20240403) {
  cfg <- sample_mfrm_latent_benchmark_data(seed = seed)
  omitted_person <- as.character(cfg$person_data$Person[1])
  cfg$person_data$X[match(omitted_person, cfg$person_data$Person)] <- NA_real_
  rows_omitted <- sum(as.character(cfg$data$Person) == omitted_person)
  cfg$population_policy <- "omit"
  cfg$truth$omitted_persons <- omitted_person
  cfg$truth$response_rows_omitted <- as.integer(rows_omitted)
  cfg$truth$response_rows_retained <- as.integer(nrow(cfg$data) - rows_omitted)
  cfg
}

sample_mfrm_gpcm_benchmark_data <- function(seed = 20260404) {
  with_preserved_rng_seed(seed, {
    persons <- paste0("P", sprintf("%03d", seq_len(600L)))
    criteria <- paste0("C", seq_len(4L))
    theta <- stats::rnorm(length(persons), mean = 0, sd = 1)
    criterion_eff <- c(C1 = -0.8, C2 = -0.2, C3 = 0.3, C4 = 0.7)
    slope_raw <- c(C1 = 0.55, C2 = 0.85, C3 = 1.20, C4 = 1.75)
    slope_truth <- slope_raw / exp(mean(log(slope_raw)))
    step_mat <- rbind(
      C1 = c(-1.30, -0.20, 1.50),
      C2 = c(-1.00, 0.00, 1.00),
      C3 = c(-1.45, -0.15, 1.60),
      C4 = c(-0.70, 0.10, 0.60)
    )
    step_mat <- t(apply(step_mat, 1, center_sum_zero))
    step_cum_mat <- t(apply(step_mat, 1, function(x) c(0, cumsum(x))))

    dat <- expand.grid(Person = persons, Criterion = criteria, stringsAsFactors = FALSE)
    criterion_idx <- match(dat$Criterion, criteria)
    eta <- theta[match(dat$Person, persons)] - criterion_eff[dat$Criterion]
    prob_mat <- category_prob_gpcm(
      eta = eta,
      step_cum_mat = step_cum_mat,
      criterion_idx = criterion_idx,
      slopes = unname(slope_truth),
      slope_idx = criterion_idx
    )
    dat$Score <- apply(prob_mat, 1, function(p) sample.int(4L, size = 1L, prob = p) - 1L)

    step_truth <- expand.grid(
      StepFacet = criteria,
      Step = paste0("Step_", seq_len(ncol(step_mat))),
      stringsAsFactors = FALSE
    )
    step_truth$Estimate <- as.vector(step_mat)

    list(
      data = dat,
      person = "Person",
      facets = "Criterion",
      score = "Score",
      step_facet = "Criterion",
      slope_facet = "Criterion",
      truth = list(
        slopes = slope_truth,
        steps = step_truth,
        criterion = criterion_eff
      )
    )
  })
}

reference_benchmark_source_profile <- function() {
  tibble::tibble(
    RuleID = c(
      "bias_obs_exp_average",
      "bias_local_measure",
      "bias_pairwise_welch",
      "linking_common_elements",
      "latent_population_omission_check",
      "conquest_overlap_bundle_check",
      "conquest_overlap_audit_check"
    ),
    Domain = c("bias", "bias", "bias", "linking", "latent_regression", "conquest_overlap", "conquest_overlap"),
    SourceLabel = c(
      "Facets Tutorial 3",
      "FACETS Table 14",
      "FACETS Table 14 / change log",
      "FACETS equating guidance",
      "mfrmr population-policy check",
      "ACER ConQuest Command Reference / mfrmr overlap boundary",
      "ACER ConQuest show-cases workflow / mfrmr audit workflow"
    ),
    SourceURL = c(
      "https://www.winsteps.com/a/ftutorial3.pdf",
      "https://www.winsteps.com/facetman/table14.htm",
      "https://www.winsteps.com/facetman/table14.htm",
      "https://www.winsteps.com/facetman/equating.htm",
      NA_character_,
      "https://conquestmanual.acer.org/s4-00.html",
      "https://conquestmanual.acer.org/s4-00.html"
    ),
    Detail = c(
      "Observed-minus-expected average should equal (Observed Score - Expected Score) / Observed Count.",
      "Local target measure in a context should equal the target's overall measure plus the context-specific bias term.",
      "Pairwise local contrasts are reported with a Rasch-Welch t statistic and approximate degrees of freedom.",
      "A practical linking audit should confirm at least 5 common elements per linking facet when equating forms or datasets.",
      "Population-model complete-case omission should be explicit in fitted metadata, response-row counts, and active person estimates.",
      "The ConQuest-overlap package-side fixture checks bundle, normalized-table, and audit preparation without claiming an executed ConQuest comparison.",
      "Case-level EAP and item/population tables must be normalized before numerical audit against the exact-overlap bundle."
    )
  )
}

resolve_reference_benchmark_data <- function(dataset) {
  if (identical(dataset, "synthetic_truth")) {
    return(list(
      data = sample_mfrm_data(seed = 20240131),
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score"
    ))
  }
  if (identical(dataset, "synthetic_latent_regression")) {
    return(sample_mfrm_latent_benchmark_data(seed = 20240403))
  }
  if (identical(dataset, "synthetic_latent_regression_omit")) {
    return(sample_mfrm_latent_omit_benchmark_data(seed = 20240403))
  }
  if (identical(dataset, "synthetic_gpcm")) {
    return(sample_mfrm_gpcm_benchmark_data(seed = 20260404))
  }
  list(
    data = load_mfrmr_data(dataset),
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score"
  )
}

score_reference_metric <- function(value, pass_max = 1e-8, warn_max = 1e-6) {
  if (!is.finite(value)) {
    return("Fail")
  }
  if (value <= pass_max) {
    return("Pass")
  }
  if (value <= warn_max) {
    return("Warn")
  }
  "Fail"
}

collect_reference_dataset_design <- function(dataset, data) {
  out <- tibble::tibble(
    Dataset = dataset,
    Rows = nrow(data),
    Persons = length(unique(as.character(data$Person))),
    Raters = if ("Rater" %in% names(data)) length(unique(as.character(data$Rater))) else NA_integer_,
    Criteria = if ("Criterion" %in% names(data)) length(unique(as.character(data$Criterion))) else NA_integer_,
    Tasks = if ("Task" %in% names(data)) length(unique(as.character(data$Task))) else NA_integer_,
    ScoreMin = suppressWarnings(min(as.numeric(data$Score), na.rm = TRUE)),
    ScoreMax = suppressWarnings(max(as.numeric(data$Score), na.rm = TRUE))
  )
  out
}

build_reference_design_checks <- function(actual_row, spec_row, case_id) {
  metrics <- c("Rows", "Persons", "Raters", "Criteria", "Tasks", "ScoreMin", "ScoreMax")
  rows <- lapply(metrics, function(metric) {
    expected <- spec_row[[metric]][1]
    actual <- actual_row[[metric]][1]
    if (is.na(expected)) {
      status <- "Skip"
      detail <- "No fixed reference target for this metric."
    } else if (identical(as.integer(actual), as.integer(expected))) {
      status <- "Pass"
      detail <- "Observed design matched the reference target."
    } else {
      status <- "Fail"
      detail <- "Observed design did not match the reference target."
    }
    tibble::tibble(
      Case = case_id,
      Dataset = as.character(actual_row$Dataset[1]),
      Metric = metric,
      Actual = as.character(actual),
      Expected = if (is.na(expected)) NA_character_ else as.character(expected),
      Status = status,
      Detail = detail
    )
  })
  dplyr::bind_rows(rows)
}

fit_reference_benchmark_dataset <- function(dataset,
                                            method = "MML",
                                            model = "RSM",
                                            quad_points = 7,
                                            maxit = 40,
                                            reltol = 1e-6,
                                            mml_engine = "direct") {
  cfg <- resolve_reference_benchmark_data(dataset)
  fit_args <- list(
    data = cfg$data,
    person = cfg$person,
    facets = cfg$facets,
    score = cfg$score,
    method = method,
    model = model,
    maxit = maxit,
    reltol = reltol
  )
  if (identical(toupper(method), "MML")) {
    fit_args$quad_points <- quad_points
    fit_args$mml_engine <- mml_engine
  }
  if (!is.null(cfg$population_formula)) {
    fit_args$population_formula <- cfg$population_formula
    fit_args$person_data <- cfg$person_data
    fit_args$person_id <- cfg$person_id %||% cfg$person
    fit_args$population_policy <- cfg$population_policy %||% "error"
  }
  if (!is.null(cfg$step_facet)) {
    fit_args$step_facet <- cfg$step_facet
  }
  if (!is.null(cfg$slope_facet)) {
    fit_args$slope_facet <- cfg$slope_facet
  }
  fit <- suppressWarnings(do.call(fit_mfrm, fit_args))
  diag <- if (identical(toupper(as.character(model)), "GPCM")) {
    list(
      overall_fit = tibble::tibble(Infit = NA_real_, Outfit = NA_real_),
      precision_profile = tibble::tibble(
        PrecisionTier = if (identical(toupper(as.character(method)), "MML")) "model_based" else "exploratory",
        SupportsFormalInference = identical(toupper(as.character(method)), "MML")
      )
    )
  } else {
    suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  }
  list(
    dataset = dataset,
    data = cfg$data,
    truth = cfg$truth %||% NULL,
    fit = fit,
    diagnostics = diag,
    design = collect_reference_dataset_design(dataset, cfg$data)
  )
}

collect_reference_fit_run <- function(case_id, fit_obj) {
  fit <- fit_obj$fit
  diag <- fit_obj$diagnostics
  design <- fit_obj$design
  population <- fit$population %||% list()
  population_formula <- population$formula %||% fit$config$population_formula %||% NULL
  population_formula_label <- if (is.null(population_formula)) {
    NA_character_
  } else {
    paste(deparse(population_formula), collapse = " ")
  }
  population_design_columns <- as.character(population$design_columns %||% character(0))
  population_coefficients <- population$coefficients %||% numeric(0)
  population_xlevel_variables <- compact_population_coding_variables(population$xlevels)
  population_contrast_variables <- compact_population_coding_variables(population$contrasts)
  tibble::tibble(
    Case = case_id,
    Dataset = fit_obj$dataset,
    Method = as.character(fit$config$method %||% NA_character_),
    Model = as.character(fit$config$model %||% NA_character_),
    Rows = as.integer(design$Rows[1]),
    Persons = as.integer(design$Persons[1]),
    Raters = as.integer(design$Raters[1]),
    Criteria = as.integer(design$Criteria[1]),
    Tasks = if ("Tasks" %in% names(design)) as.integer(design$Tasks[1]) else NA_integer_,
    Converged = isTRUE(fit$summary$Converged),
    LogLik = suppressWarnings(as.numeric(fit$summary$LogLik %||% NA_real_)),
    MMLEngineRequested = as.character(fit$summary$MMLEngineRequested[1] %||% NA_character_),
    MMLEngineUsed = as.character(fit$summary$MMLEngineUsed[1] %||% NA_character_),
    EMIterations = suppressWarnings(as.integer(fit$summary$EMIterations[1] %||% NA_integer_)),
    PosteriorBasis = as.character(
      fit$config$posterior_basis %||%
        population$posterior_basis %||%
        if (isTRUE(population$active)) "population_model" else "legacy_mml"
    ),
    PopulationModelActive = isTRUE(population$active),
    PopulationFormula = population_formula_label,
    PopulationPolicy = as.character(population$policy %||% NA_character_),
    PopulationDesignColumns = paste(population_design_columns, collapse = ", "),
    PopulationXlevelVariables = population_xlevel_variables,
    PopulationContrastVariables = population_contrast_variables,
    PopulationCoefficientCount = as.integer(length(population_coefficients)),
    PopulationResidualVariance = suppressWarnings(as.numeric(population$sigma2 %||% NA_real_)),
    PopulationIncludedPersons = as.integer(length(population$included_persons %||% character(0))),
    PopulationOmittedPersons = as.integer(length(population$omitted_persons %||% character(0))),
    PopulationResponseRowsRetained = suppressWarnings(as.integer(population$response_rows_retained %||% NA_integer_)),
    PopulationResponseRowsOmitted = suppressWarnings(as.integer(population$response_rows_omitted %||% NA_integer_)),
    Infit = suppressWarnings(as.numeric(diag$overall_fit$Infit[1] %||% NA_real_)),
    Outfit = suppressWarnings(as.numeric(diag$overall_fit$Outfit[1] %||% NA_real_)),
    PrecisionTier = as.character(diag$precision_profile$PrecisionTier[1] %||% NA_character_),
    SupportsFormalInference = isTRUE(diag$precision_profile$SupportsFormalInference[1] %||% FALSE)
  )
}

build_truth_recovery_checks <- function(case_id, fit_obj) {
  fit_tbl <- fit_obj$fit$facets$others
  truth_targets <- synthetic_truth_targets()

  rows <- lapply(names(truth_targets), function(facet_name) {
    truth <- truth_targets[[facet_name]]
    est_tbl <- fit_tbl[fit_tbl$Facet == facet_name, c("Level", "Estimate"), drop = FALSE]
    est_tbl <- est_tbl[order(est_tbl$Level), , drop = FALSE]
    est <- suppressWarnings(as.numeric(est_tbl$Estimate))
    est_centered <- est - mean(est, na.rm = TRUE)
    truth_centered <- truth - mean(truth)
    corr <- suppressWarnings(stats::cor(est_centered, truth_centered))
    mae <- mean(abs(est_centered - truth_centered), na.rm = TRUE)

    status <- if (is.finite(corr) && corr >= 0.95 && is.finite(mae) && mae <= 0.30) {
      "Pass"
    } else if (is.finite(corr) && corr >= 0.90 && is.finite(mae) && mae <= 0.45) {
      "Warn"
    } else {
      "Fail"
    }

    tibble::tibble(
      Case = case_id,
      Facet = facet_name,
      Correlation = corr,
      MeanAbsoluteDeviation = mae,
      Status = status,
      Detail = if (identical(status, "Pass")) {
        "Recovered facet ordering and spacing were close to the known generating values."
      } else if (identical(status, "Warn")) {
        "Recovered facet ordering was acceptable, but spacing deviated from the reference profile."
      } else {
        "Recovered facet ordering or spacing missed the reference profile."
      }
    )
  })

  dplyr::bind_rows(rows)
}

build_latent_recovery_checks <- function(case_id, fit_obj) {
  truth <- fit_obj$truth %||% NULL
  if (is.null(truth)) {
    return(tibble::tibble())
  }

  coeff_est <- as.numeric(fit_obj$fit$population$coefficients %||% numeric(0))
  names(coeff_est) <- names(fit_obj$fit$population$coefficients %||% NULL)
  coeff_truth <- truth$coefficients %||% numeric(0)
  coeff_names <- intersect(names(coeff_truth), names(coeff_est))

  coeff_rows <- lapply(coeff_names, function(term) {
    err <- abs(as.numeric(coeff_est[term]) - as.numeric(coeff_truth[term]))
    status <- if (err <= 0.20) {
      "Pass"
    } else if (err <= 0.35) {
      "Warn"
    } else {
      "Fail"
    }
    tibble::tibble(
      Case = case_id,
      Facet = paste0("Population:", term),
      Correlation = NA_real_,
      MeanAbsoluteDeviation = err,
      Status = status,
      Detail = "Population-model coefficient recovery is monitored by absolute error against the known generating value."
    )
  })

  sigma_err <- abs(as.numeric(fit_obj$fit$population$sigma2) - as.numeric(truth$sigma2))
  sigma_row <- tibble::tibble(
    Case = case_id,
    Facet = "Population:sigma2",
    Correlation = NA_real_,
    MeanAbsoluteDeviation = sigma_err,
    Status = if (sigma_err <= 0.20) {
      "Pass"
    } else if (sigma_err <= 0.35) {
      "Warn"
    } else {
      "Fail"
    },
    Detail = "Residual latent-variance recovery is monitored by absolute error against the known generating value."
  )

  criterion_tbl <- as.data.frame(fit_obj$fit$facets$others, stringsAsFactors = FALSE)
  criterion_tbl <- criterion_tbl[criterion_tbl$Facet == "Criterion", c("Level", "Estimate"), drop = FALSE]
  criterion_tbl <- criterion_tbl[order(criterion_tbl$Level), , drop = FALSE]
  truth_criterion <- truth$criterion[criterion_tbl$Level]
  est_centered <- suppressWarnings(as.numeric(criterion_tbl$Estimate))
  est_centered <- est_centered - mean(est_centered, na.rm = TRUE)
  truth_centered <- as.numeric(truth_criterion) - mean(as.numeric(truth_criterion), na.rm = TRUE)
  criterion_corr <- suppressWarnings(stats::cor(est_centered, truth_centered))
  criterion_mae <- mean(abs(est_centered - truth_centered), na.rm = TRUE)
  criterion_row <- tibble::tibble(
    Case = case_id,
    Facet = "Criterion",
    Correlation = criterion_corr,
    MeanAbsoluteDeviation = criterion_mae,
    Status = if (is.finite(criterion_corr) && criterion_corr >= 0.95 && criterion_mae <= 0.20) {
      "Pass"
    } else if (is.finite(criterion_corr) && criterion_corr >= 0.90 && criterion_mae <= 0.30) {
      "Warn"
    } else {
      "Fail"
    },
    Detail = "Criterion recovery is reviewed after centering because the Rasch location origin remains unidentified."
  )

  pred_new <- dplyr::bind_rows(lapply(c("LOW", "HIGH"), function(id) {
    data.frame(
      Person = id,
      Criterion = names(truth$criterion),
      Score = c(0L, 0L, 0L, 1L, 1L, 1L),
      stringsAsFactors = FALSE
    )
  }))
  pred_person <- data.frame(Person = c("LOW", "HIGH"), X = c(-1.5, 1.5), stringsAsFactors = FALSE)
  pred_obj <- predict_mfrm_units(
    fit = fit_obj$fit,
    new_data = pred_new,
    person_data = pred_person,
    n_draws = 0
  )
  pred_tbl <- as.data.frame(pred_obj$estimates, stringsAsFactors = FALSE)
  shift <- with(pred_tbl, Estimate[Person == "HIGH"][1] - Estimate[Person == "LOW"][1])
  shift_row <- tibble::tibble(
    Case = case_id,
    Facet = "Population:posterior_shift",
    Correlation = NA_real_,
    MeanAbsoluteDeviation = NA_real_,
    Status = if (is.finite(shift) && shift > 0) "Pass" else "Fail",
    Detail = "Posterior scoring for matched response patterns should shift upward when the scored covariate value is higher."
  )

  dplyr::bind_rows(dplyr::bind_rows(coeff_rows), sigma_row, criterion_row, shift_row)
}

build_gpcm_recovery_checks <- function(case_id, fit_obj) {
  truth <- fit_obj$truth %||% NULL
  if (is.null(truth)) {
    return(tibble::tibble())
  }

  slope_truth <- truth$slopes %||% numeric(0)
  slope_est_tbl <- as.data.frame(fit_obj$fit$slopes %||% data.frame(), stringsAsFactors = FALSE)
  slope_est_tbl <- slope_est_tbl[, c("SlopeFacet", "Estimate"), drop = FALSE]
  slope_tbl <- merge(
    data.frame(SlopeFacet = names(slope_truth), Truth = as.numeric(slope_truth), stringsAsFactors = FALSE),
    slope_est_tbl,
    by = "SlopeFacet",
    sort = FALSE
  )
  slope_truth_log <- log(slope_tbl$Truth) - mean(log(slope_tbl$Truth), na.rm = TRUE)
  slope_est_log <- log(slope_tbl$Estimate) - mean(log(slope_tbl$Estimate), na.rm = TRUE)
  slope_corr <- suppressWarnings(stats::cor(slope_truth_log, slope_est_log))
  slope_mae <- mean(abs(slope_truth_log - slope_est_log), na.rm = TRUE)
  slope_row <- tibble::tibble(
    Case = case_id,
    Facet = "GPCM:slopes",
    Correlation = slope_corr,
    MeanAbsoluteDeviation = slope_mae,
    Status = if (is.finite(slope_corr) && slope_corr >= 0.95 && slope_mae <= 0.15) {
      "Pass"
    } else if (is.finite(slope_corr) && slope_corr >= 0.85 && slope_mae <= 0.25) {
      "Warn"
    } else {
      "Fail"
    },
    Detail = "Slope recovery is compared on the centered log-discrimination scale implied by the package's geometric-mean-one identification."
  )

  step_truth_tbl <- as.data.frame(truth$steps %||% data.frame(), stringsAsFactors = FALSE)
  step_est_tbl <- as.data.frame(fit_obj$fit$steps %||% data.frame(), stringsAsFactors = FALSE)
  step_tbl <- merge(
    step_truth_tbl[, c("StepFacet", "Step", "Estimate"), drop = FALSE],
    step_est_tbl[, c("StepFacet", "Step", "Estimate"), drop = FALSE],
    by = c("StepFacet", "Step"),
    sort = FALSE,
    suffixes = c(".Truth", ".Estimate")
  )
  step_corr <- suppressWarnings(stats::cor(step_tbl$Estimate.Truth, step_tbl$Estimate.Estimate))
  step_mae <- mean(abs(step_tbl$Estimate.Truth - step_tbl$Estimate.Estimate), na.rm = TRUE)
  step_row <- tibble::tibble(
    Case = case_id,
    Facet = "GPCM:steps",
    Correlation = step_corr,
    MeanAbsoluteDeviation = step_mae,
    Status = if (is.finite(step_corr) && step_corr >= 0.98 && step_mae <= 0.20) {
      "Pass"
    } else if (is.finite(step_corr) && step_corr >= 0.93 && step_mae <= 0.30) {
      "Warn"
    } else {
      "Fail"
    },
    Detail = "Step recovery is reviewed on the row-centered step scale stored by the current bounded GPCM branch."
  )

  criterion_tbl <- as.data.frame(fit_obj$fit$facets$others, stringsAsFactors = FALSE)
  criterion_tbl <- criterion_tbl[criterion_tbl$Facet == "Criterion", c("Level", "Estimate"), drop = FALSE]
  criterion_tbl <- criterion_tbl[order(criterion_tbl$Level), , drop = FALSE]
  truth_criterion <- truth$criterion[criterion_tbl$Level]
  est_centered <- suppressWarnings(as.numeric(criterion_tbl$Estimate))
  est_centered <- est_centered - mean(est_centered, na.rm = TRUE)
  truth_centered <- as.numeric(truth_criterion) - mean(as.numeric(truth_criterion), na.rm = TRUE)
  criterion_corr <- suppressWarnings(stats::cor(est_centered, truth_centered))
  criterion_mae <- mean(abs(est_centered - truth_centered), na.rm = TRUE)
  criterion_row <- tibble::tibble(
    Case = case_id,
    Facet = "Criterion",
    Correlation = criterion_corr,
    MeanAbsoluteDeviation = criterion_mae,
    Status = if (is.finite(criterion_corr) && criterion_corr >= 0.98 && criterion_mae <= 0.15) {
      "Pass"
    } else if (is.finite(criterion_corr) && criterion_corr >= 0.93 && criterion_mae <= 0.25) {
      "Warn"
    } else {
      "Fail"
    },
    Detail = "Criterion recovery is reviewed after centering because the location origin remains unidentified."
  )

  dplyr::bind_rows(slope_row, step_row, criterion_row)
}

build_pair_stability_checks <- function(case_id, primary_fit_obj, reference_fit_obj) {
  primary_fit <- primary_fit_obj$fit
  reference_fit <- reference_fit_obj$fit
  primary_diag <- primary_fit_obj$diagnostics
  reference_diag <- reference_fit_obj$diagnostics

  primary_est <- as.data.frame(primary_fit$facets$others, stringsAsFactors = FALSE)
  reference_est <- as.data.frame(reference_fit$facets$others, stringsAsFactors = FALSE)
  shared <- merge(
    primary_est[, c("Facet", "Level", "Estimate"), drop = FALSE],
    reference_est[, c("Facet", "Level", "Estimate"), drop = FALSE],
    by = c("Facet", "Level"),
    suffixes = c("_primary", "_reference")
  )

  primary_rel <- as.data.frame(primary_diag$reliability, stringsAsFactors = FALSE)
  primary_rel <- primary_rel[primary_rel$Facet != "Person", c("Facet", "Reliability", "Separation"), drop = FALSE]
  reference_rel <- as.data.frame(reference_diag$reliability, stringsAsFactors = FALSE)
  reference_rel <- reference_rel[reference_rel$Facet != "Person", c("Facet", "Reliability", "Separation"), drop = FALSE]
  rel_shared <- merge(primary_rel, reference_rel, by = "Facet", suffixes = c("_primary", "_reference"))

  facet_rows <- lapply(unique(shared$Facet), function(facet_name) {
    facet_tbl <- shared[shared$Facet == facet_name, , drop = FALSE]
    est_primary <- suppressWarnings(as.numeric(facet_tbl$Estimate_primary))
    est_reference <- suppressWarnings(as.numeric(facet_tbl$Estimate_reference))
    pearson <- suppressWarnings(stats::cor(est_primary, est_reference))
    spearman <- suppressWarnings(stats::cor(est_primary, est_reference, method = "spearman"))
    mae <- mean(abs(est_primary - est_reference), na.rm = TRUE)

    rel_row <- rel_shared[rel_shared$Facet == facet_name, , drop = FALSE]
    reliability_gap <- if (nrow(rel_row) > 0) {
      abs(suppressWarnings(as.numeric(rel_row$Reliability_primary[1])) -
            suppressWarnings(as.numeric(rel_row$Reliability_reference[1])))
    } else {
      NA_real_
    }
    separation_gap <- if (nrow(rel_row) > 0) {
      abs(suppressWarnings(as.numeric(rel_row$Separation_primary[1])) -
            suppressWarnings(as.numeric(rel_row$Separation_reference[1])))
    } else {
      NA_real_
    }

    if (identical(facet_name, "Criterion")) {
      status <- if (is.finite(pearson) && pearson >= 0.95 && is.finite(mae) && mae <= 0.10) {
        "Pass"
      } else if (is.finite(pearson) && pearson >= 0.90 && is.finite(mae) && mae <= 0.15) {
        "Warn"
      } else {
        "Fail"
      }
      detail <- "Criterion measures should remain highly aligned across the paired calibration datasets."
    } else {
      status <- if (is.finite(spearman) && spearman >= 0.50 && is.finite(mae) && mae <= 0.35) {
        "Pass"
      } else if (is.finite(spearman) && spearman >= 0.40 && is.finite(mae) && mae <= 0.45) {
        "Warn"
      } else {
        "Fail"
      }
      detail <- "Rater measures may shift more under iterative recalibration, so rank stability and average deviation are tracked together."
    }

    tibble::tibble(
      Case = case_id,
      Facet = facet_name,
      Pearson = pearson,
      Spearman = spearman,
      MeanAbsoluteDifference = mae,
      ReliabilityGap = reliability_gap,
      SeparationGap = separation_gap,
      Status = status,
      Detail = detail
    )
  })

  overall_row <- tibble::tibble(
    Case = case_id,
    Facet = "OverallFit",
    Pearson = NA_real_,
    Spearman = NA_real_,
    MeanAbsoluteDifference = NA_real_,
    ReliabilityGap = NA_real_,
    SeparationGap = NA_real_,
    InfitDelta = abs(suppressWarnings(as.numeric(primary_diag$overall_fit$Infit[1])) -
                       suppressWarnings(as.numeric(reference_diag$overall_fit$Infit[1]))),
    OutfitDelta = abs(suppressWarnings(as.numeric(primary_diag$overall_fit$Outfit[1])) -
                        suppressWarnings(as.numeric(reference_diag$overall_fit$Outfit[1]))),
    Status = NA_character_,
    Detail = "Global fit deltas summarize whether the paired reference datasets stay in the same fit regime."
  )
  overall_row$Status <- if (is.finite(overall_row$InfitDelta) &&
                              is.finite(overall_row$OutfitDelta) &&
                              overall_row$InfitDelta <= 0.10 &&
                              overall_row$OutfitDelta <= 0.10) {
    "Pass"
  } else if (is.finite(overall_row$InfitDelta) &&
             is.finite(overall_row$OutfitDelta) &&
             overall_row$InfitDelta <= 0.20 &&
             overall_row$OutfitDelta <= 0.20) {
    "Warn"
  } else {
    "Fail"
  }

  dplyr::bind_rows(dplyr::bind_rows(facet_rows), overall_row)
}

build_linking_guideline_checks <- function(case_id, primary_fit_obj, reference_fit_obj, guideline_min = 5L) {
  primary_data <- primary_fit_obj$data
  reference_data <- reference_fit_obj$data
  shared_facets <- intersect(
    intersect(names(primary_data), names(reference_data)),
    c("Rater", "Criterion", "Task")
  )
  if (length(shared_facets) == 0) {
    return(tibble::tibble())
  }

  dplyr::bind_rows(lapply(shared_facets, function(facet_name) {
    common_count <- length(intersect(
      unique(as.character(primary_data[[facet_name]])),
      unique(as.character(reference_data[[facet_name]]))
    ))
    status <- if (common_count >= guideline_min) {
      "Pass"
    } else if (common_count >= 1) {
      "Warn"
    } else {
      "Fail"
    }
    tibble::tibble(
      Case = case_id,
      Facet = facet_name,
      CommonElements = as.integer(common_count),
      GuidelineMinimum = as.integer(guideline_min),
      Status = status,
      Detail = if (identical(status, "Pass")) {
        "Common-element coverage satisfied the package's linking audit rule."
      } else if (identical(status, "Warn")) {
        "Common-element coverage fell below the preferred linking rule-of-thumb."
      } else {
        "No common elements were available for this linking facet."
      }
    )
  }))
}

build_bias_contract_checks <- function(case_id,
                                       fit_obj,
                                       interaction_facets = c("Rater", "Task"),
                                       target_facet = "Rater",
                                       context_facet = "Task") {
  data_names <- names(fit_obj$data)
  if (!all(interaction_facets %in% data_names)) {
    return(tibble::tibble())
  }

  bias_obj <- suppressWarnings(estimate_bias(
    fit_obj$fit,
    fit_obj$diagnostics,
    interaction_facets = interaction_facets,
    max_iter = 2
  ))
  if (!inherits(bias_obj, "mfrm_bias") || is.null(bias_obj$table) || nrow(bias_obj$table) == 0) {
    return(tibble::tibble())
  }

  bias_tbl <- as.data.frame(bias_obj$table, stringsAsFactors = FALSE)
  pair_tbl <- as.data.frame(
    bias_pairwise_report(
      bias_obj,
      target_facet = target_facet,
      context_facet = context_facet,
      top_n = 200
    )$table,
    stringsAsFactors = FALSE
  )

  obs_exp_err <- with(
    bias_tbl,
    abs(`Obs-Exp Average` - ((`Observd Score` - `Expctd Score`) / `Observd Count`))
  )
  df_err <- with(bias_tbl, abs(`d.f.` - (`Observd Count` - 1)))
  finite_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    max(x)
  }
  obs_exp_max <- finite_max(obs_exp_err)
  df_max <- finite_max(df_err)

  rows <- list(
    tibble::tibble(
      Case = case_id,
      Metric = "ObsExpAverageIdentity",
      MaxError = obs_exp_max,
      Status = score_reference_metric(obs_exp_max),
      Detail = "Observed-minus-expected averages matched the score/count identity used in package-native bias tables."
    ),
    tibble::tibble(
      Case = case_id,
      Metric = "BiasDFIdentity",
      MaxError = df_max,
      Status = score_reference_metric(df_max),
      Detail = "Cell-level bias degrees of freedom matched the observed-count minus 1 approximation."
    )
  )

  if (nrow(pair_tbl) > 0) {
    use_a <- isTRUE(bias_tbl$FacetA[1] == target_facet)
    target_level_col <- if (use_a) "FacetA_Level" else "FacetB_Level"
    context_level_col <- if (use_a) "FacetB_Level" else "FacetA_Level"
    bias_lookup <- bias_tbl[, c(target_level_col, context_level_col, "Bias Size", "S.E.", "d.f."), drop = FALSE]
    names(bias_lookup) <- c("Target", "Context", "BiasSize", "BiasSE", "BiasDF")

    lookup_key <- paste(bias_lookup$Target, bias_lookup$Context, sep = "\r")
    pair1_idx <- match(paste(pair_tbl$Target, pair_tbl$Context1, sep = "\r"), lookup_key)
    pair2_idx <- match(paste(pair_tbl$Target, pair_tbl$Context2, sep = "\r"), lookup_key)
    pair1_lookup <- bias_lookup[pair1_idx, , drop = FALSE]
    pair2_lookup <- bias_lookup[pair2_idx, , drop = FALSE]

    local1_err <- abs((suppressWarnings(as.numeric(pair_tbl$`Local Measure1`)) -
                         suppressWarnings(as.numeric(pair_tbl$`Target Measure`))) -
                        suppressWarnings(as.numeric(pair1_lookup$BiasSize)))
    local2_err <- abs((suppressWarnings(as.numeric(pair_tbl$`Local Measure2`)) -
                         suppressWarnings(as.numeric(pair_tbl$`Target Measure`))) -
                        suppressWarnings(as.numeric(pair2_lookup$BiasSize)))
    contrast_err <- with(pair_tbl, abs(Contrast - (`Local Measure1` - `Local Measure2`)))
    se_expected <- sqrt(suppressWarnings(as.numeric(pair1_lookup$BiasSE))^2 +
                          suppressWarnings(as.numeric(pair2_lookup$BiasSE))^2)
    se_err <- abs(suppressWarnings(as.numeric(pair_tbl$SE)) - se_expected)
    df_expected <- mapply(
      function(se1, se2, df1, df2) {
        welch_satterthwaite_df(c(se1^2, se2^2), c(df1, df2))
      },
      se1 = suppressWarnings(as.numeric(pair1_lookup$BiasSE)),
      se2 = suppressWarnings(as.numeric(pair2_lookup$BiasSE)),
      df1 = suppressWarnings(as.numeric(pair1_lookup$BiasDF)),
      df2 = suppressWarnings(as.numeric(pair2_lookup$BiasDF))
    )
    df_pair_err <- abs(suppressWarnings(as.numeric(pair_tbl$`d.f.`)) - df_expected)
    local_max <- finite_max(c(local1_err, local2_err))
    contrast_max <- finite_max(contrast_err)
    se_max <- finite_max(se_err)
    df_pair_max <- finite_max(df_pair_err)

    rows <- c(rows, list(
      tibble::tibble(
        Case = case_id,
        Metric = "LocalMeasureIdentity",
        MaxError = local_max,
        Status = score_reference_metric(local_max),
        Detail = "Local target measures matched overall target measures plus the corresponding context-specific bias terms."
      ),
      tibble::tibble(
        Case = case_id,
        Metric = "PairContrastIdentity",
        MaxError = contrast_max,
        Status = score_reference_metric(contrast_max),
        Detail = "Pairwise contrasts matched the difference between the two local target measures."
      ),
      tibble::tibble(
        Case = case_id,
        Metric = "PairSEIdentity",
        MaxError = se_max,
        Status = score_reference_metric(se_max),
        Detail = "Pairwise standard errors matched the joint local-measure standard error identity."
      ),
      tibble::tibble(
        Case = case_id,
        Metric = "PairDFIdentity",
        MaxError = df_pair_max,
        Status = score_reference_metric(df_pair_max),
        Detail = "Pairwise degrees of freedom matched the Rasch-Welch Satterthwaite approximation."
      )
    ))
  }

  dplyr::bind_rows(rows)
}

build_conquest_overlap_dry_run_checks <- function(case_id, fit_obj) {
  bundle <- build_conquest_overlap_bundle(fit = fit_obj$fit)
  normalized <- normalize_conquest_overlap_tables(
    conquest_population = data.frame(
      Term = bundle$mfrmr_population$Parameter,
      Est = bundle$mfrmr_population$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_item_estimates = data.frame(
      ItemID = bundle$mfrmr_item_estimates$ResponseVar,
      Est = bundle$mfrmr_item_estimates$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_case_eap = data.frame(
      PID = bundle$mfrmr_case_eap$Person,
      EAP = bundle$mfrmr_case_eap$Estimate,
      stringsAsFactors = FALSE
    ),
    conquest_population_term = "Term",
    conquest_population_estimate = "Est",
    conquest_item_id = "ItemID",
    conquest_item_estimate = "Est",
    conquest_case_person = "PID",
    conquest_case_estimate = "EAP"
  )
  audit <- audit_conquest_overlap(bundle, normalized)
  overall <- as.data.frame(audit$overall, stringsAsFactors = FALSE)

  row_metric <- function(metric, value, pass_max, warn_max, detail) {
    tibble::tibble(
      Case = case_id,
      Metric = metric,
      Actual = suppressWarnings(as.numeric(value[1] %||% NA_real_)),
      PassMax = pass_max,
      WarnMax = warn_max,
      Status = score_reference_metric(
        abs(suppressWarnings(as.numeric(value[1] %||% NA_real_))),
        pass_max = pass_max,
        warn_max = warn_max
      ),
      Detail = detail
    )
  }

  dplyr::bind_rows(
    row_metric(
      "AttentionItems",
      overall$AttentionItems[1],
      pass_max = 0,
      warn_max = 0,
      detail = "The package-side check should produce no missing, duplicate, or non-numeric audit attention items."
    ),
    row_metric(
      "PopulationMaxAbsDifference",
      overall$PopulationMaxAbsDifference[1],
      pass_max = 1e-10,
      warn_max = 1e-8,
      detail = "Copied population parameters should round-trip through the normalized-table audit workflow."
    ),
    row_metric(
      "ItemCenteredMaxAbsDifference",
      overall$ItemCenteredMaxAbsDifference[1],
      pass_max = 1e-10,
      warn_max = 1e-8,
      detail = "Copied item estimates should agree after centering under the overlap audit workflow."
    ),
    row_metric(
      "CaseMaxAbsDifference",
      overall$CaseMaxAbsDifference[1],
      pass_max = 1e-10,
      warn_max = 1e-8,
      detail = "Copied case EAP estimates should round-trip through the overlap audit workflow."
    )
  )
}

build_latent_omission_contract_checks <- function(case_id, fit_obj) {
  fit <- fit_obj$fit
  pop <- fit$population %||% list()
  truth <- fit_obj$truth %||% list()
  expected_omitted <- as.character(truth$omitted_persons %||% character(0))
  actual_omitted <- as.character(pop$omitted_persons %||% character(0))

  person_tbl <- as.data.frame(fit$facets$person %||% data.frame(), stringsAsFactors = FALSE)
  fitted_persons <- as.character(person_tbl$Person %||% character(0))
  replay_tbl <- as.data.frame(pop$person_table_replay %||% data.frame(), stringsAsFactors = FALSE)
  replay_id <- as.character(pop$person_id %||% "Person")
  replay_persons <- if (replay_id %in% names(replay_tbl)) as.character(replay_tbl[[replay_id]]) else character(0)

  check_row <- function(metric, actual, expected, pass, detail) {
    tibble::tibble(
      Case = case_id,
      Metric = metric,
      Actual = paste(as.character(actual), collapse = ", "),
      Expected = paste(as.character(expected), collapse = ", "),
      Status = if (isTRUE(pass)) "Pass" else "Fail",
      Detail = detail
    )
  }

  dplyr::bind_rows(
    check_row(
      "PopulationPolicy",
      pop$policy %||% NA_character_,
      "omit",
      identical(as.character(pop$policy %||% NA_character_), "omit"),
      "The omission fixture should run under the documented complete-case omission policy."
    ),
    check_row(
      "PopulationOmittedPersons",
      length(actual_omitted),
      length(expected_omitted),
      setequal(actual_omitted, expected_omitted),
      "The fitted population scaffold should omit exactly the background-incomplete person(s)."
    ),
    check_row(
      "PopulationResponseRowsOmitted",
      as.integer(pop$response_rows_omitted %||% NA_integer_),
      as.integer(truth$response_rows_omitted %||% NA_integer_),
      identical(
        as.integer(pop$response_rows_omitted %||% NA_integer_),
        as.integer(truth$response_rows_omitted %||% NA_integer_)
      ),
      "All response rows for omitted persons should be excluded from the active latent-regression fit."
    ),
    check_row(
      "PopulationResponseRowsRetained",
      as.integer(pop$response_rows_retained %||% NA_integer_),
      as.integer(truth$response_rows_retained %||% NA_integer_),
      identical(
        as.integer(pop$response_rows_retained %||% NA_integer_),
        as.integer(truth$response_rows_retained %||% NA_integer_)
      ),
      "The retained-response count should match the complete-case scaffold after omission."
    ),
    check_row(
      "OmittedPersonExcludedFromEstimates",
      !any(expected_omitted %in% fitted_persons),
      TRUE,
      length(expected_omitted) > 0L && !any(expected_omitted %in% fitted_persons),
      "Omitted persons should not receive active person estimates from the fitted calibration."
    ),
    check_row(
      "OmittedPersonPreservedForReplay",
      all(expected_omitted %in% replay_persons),
      TRUE,
      length(expected_omitted) > 0L && all(expected_omitted %in% replay_persons),
      "The replay provenance table should preserve omitted person IDs for auditability."
    )
  )
}

summarize_reference_benchmark_case <- function(case_id, case_type, fit_runs, design_checks, recovery_checks, pair_checks, bias_checks, linking_checks, conquest_overlap_checks, population_policy_checks) {
  subset_reference_case <- function(tbl, case_id, drop_skip = FALSE) {
    if (!is.data.frame(tbl) || !("Case" %in% names(tbl)) || nrow(tbl) == 0L) {
      return(tibble::tibble())
    }
    out <- tbl[tbl$Case == case_id, , drop = FALSE]
    if (drop_skip && "Status" %in% names(out)) {
      out <- out[out$Status != "Skip", , drop = FALSE]
    }
    out
  }

  case_statuses <- function(tbl) {
    if (!is.data.frame(tbl) || !("Status" %in% names(tbl)) || nrow(tbl) == 0L) {
      return(character(0))
    }
    as.character(tbl$Status)
  }

  case_fit <- subset_reference_case(fit_runs, case_id)
  case_design <- subset_reference_case(design_checks, case_id, drop_skip = TRUE)
  case_recovery <- subset_reference_case(recovery_checks, case_id)
  case_pairs <- subset_reference_case(pair_checks, case_id)
  case_bias <- subset_reference_case(bias_checks, case_id)
  case_link <- subset_reference_case(linking_checks, case_id)
  case_conquest <- subset_reference_case(conquest_overlap_checks, case_id)
  case_population_policy <- subset_reference_case(population_policy_checks, case_id)

  statuses <- c(
    case_fit$Converged == TRUE,
    case_statuses(case_design),
    case_statuses(case_recovery),
    case_statuses(case_pairs),
    case_statuses(case_bias),
    case_statuses(case_link),
    case_statuses(case_conquest),
    case_statuses(case_population_policy)
  )
  normalized <- ifelse(statuses %in% c(TRUE, "Pass"), "Pass",
                       ifelse(statuses %in% c("Warn"), "Warn",
                              ifelse(statuses %in% c(FALSE, "Fail"), "Fail", NA_character_)))
  missing_expected_checks <- switch(
    case_type,
    truth_recovery = nrow(case_recovery) == 0L,
    latent_recovery = nrow(case_recovery) == 0L,
    gpcm_recovery = nrow(case_recovery) == 0L,
    bias_contract = nrow(case_bias) == 0L,
    conquest_overlap_dry_run = nrow(case_conquest) == 0L,
    latent_omission_contract = nrow(case_population_policy) == 0L,
    pair_stability = nrow(case_pairs) == 0L && nrow(case_link) == 0L,
    FALSE
  )

  overall_status <- if (any(normalized == "Fail", na.rm = TRUE)) {
    "Fail"
  } else if (any(normalized == "Warn", na.rm = TRUE)) {
    "Warn"
  } else if (isTRUE(missing_expected_checks)) {
    "Warn"
  } else {
    "Pass"
  }

  key_signal <- if (identical(case_type, "truth_recovery")) {
    if (nrow(case_recovery) > 0) {
      paste0(
        "Min recovery correlation = ",
        formatC(min(case_recovery$Correlation, na.rm = TRUE), format = "f", digits = 3)
      )
    } else {
      "No recovery checks were produced."
    }
  } else if (identical(case_type, "latent_recovery")) {
    slope_row <- case_recovery[case_recovery$Facet == "Population:X", , drop = FALSE]
    if (nrow(slope_row) > 0) {
      paste0(
        "Population slope absolute error = ",
        formatC(slope_row$MeanAbsoluteDeviation[1], format = "f", digits = 3)
      )
    } else {
      "No latent-recovery checks were produced."
    }
  } else if (identical(case_type, "gpcm_recovery")) {
    slope_row <- case_recovery[case_recovery$Facet == "GPCM:slopes", , drop = FALSE]
    if (nrow(slope_row) > 0) {
      paste0(
        "GPCM slope log-scale correlation = ",
        formatC(slope_row$Correlation[1], format = "f", digits = 3)
      )
    } else {
      "No GPCM-recovery checks were produced."
    }
  } else if (identical(case_type, "bias_contract")) {
    if (nrow(case_bias) > 0) {
      paste0(
        "Max bias-identity error = ",
        formatC(max(case_bias$MaxError, na.rm = TRUE), format = "f", digits = 6)
      )
    } else {
      "No bias checks were produced."
    }
  } else if (identical(case_type, "conquest_overlap_dry_run")) {
    attention_row <- case_conquest[case_conquest$Metric == "AttentionItems", , drop = FALSE]
    if (nrow(attention_row) > 0) {
      paste0(
        "ConQuest-overlap package-side attention items = ",
        formatC(attention_row$Actual[1], format = "f", digits = 0)
      )
    } else {
      "No ConQuest-overlap package-side checks were produced."
    }
  } else if (identical(case_type, "latent_omission_contract")) {
    omitted_row <- case_population_policy[case_population_policy$Metric == "PopulationResponseRowsOmitted", , drop = FALSE]
    if (nrow(omitted_row) > 0) {
      paste0(
        "Population response rows omitted = ",
        as.character(omitted_row$Actual[1])
      )
    } else {
      "No latent-regression omission checks were produced."
    }
  } else {
    facet_rows <- case_pairs[case_pairs$Facet != "OverallFit", , drop = FALSE]
    if (nrow(facet_rows) > 0) {
      paste0(
        "Min paired rank correlation = ",
        formatC(min(facet_rows$Spearman, na.rm = TRUE), format = "f", digits = 3)
      )
    } else {
      "No pair-stability checks were produced."
    }
  }

  tibble::tibble(
    Case = case_id,
    CaseType = case_type,
    Status = overall_status,
    Fits = nrow(case_fit),
    DesignChecks = nrow(case_design),
    RecoveryChecks = nrow(case_recovery),
    BiasChecks = nrow(case_bias),
    LinkingChecks = nrow(case_link),
    ConQuestOverlapChecks = nrow(case_conquest),
    PopulationPolicyChecks = nrow(case_population_policy),
    StabilityChecks = nrow(case_pairs),
    KeySignal = key_signal
  )
}

#' Benchmark packaged reference cases
#'
#' @param cases Reference cases to run. Defaults to the standard
#'   `RSM`-compatible reference suite. Specialized `GPCM` and
#'   ConQuest-overlap package-side cases can be requested explicitly.
#' @param method Estimation method passed to [fit_mfrm()]. Defaults to `"MML"`.
#' @param model Model family passed to [fit_mfrm()]. Defaults to `"RSM"`.
#' @param quad_points Quadrature points for `method = "MML"`.
#' @param maxit Maximum optimizer iterations passed to [fit_mfrm()].
#' @param reltol Convergence tolerance passed to [fit_mfrm()].
#' @param mml_engine MML optimization engine passed to [fit_mfrm()]. Applies
#'   only when `method = "MML"`.
#'
#' @details
#' This function checks `mfrmr` against the package's curated reference case
#' families:
#' - `synthetic_truth`: checks whether recovered facet measures align with the
#'   known generating values from the package's synthetic design.
#' - `synthetic_latent_regression`: checks whether the first-version
#'   latent-regression `MML` branch recovers known population coefficients,
#'   residual latent variance, criterion ordering, and posterior-shift
#'   direction from a synthetic overlap case.
#' - `synthetic_latent_regression_omit`: checks whether the population-model
#'   complete-case omission policy is reflected in the fitted metadata,
#'   response-row audit, active person estimates, and replay provenance.
#' - `synthetic_conquest_overlap_dry_run`: builds the narrow ConQuest-overlap
#'   bundle for the latent-regression synthetic case, round-trips package tables
#'   through the normalization/audit helpers, and confirms the package-side
#'   workflow without claiming that ConQuest itself was executed.
#' - `synthetic_gpcm`: checks whether the bounded `GPCM` branch recovers
#'   known criterion-specific slopes, row-centered step parameters, and
#'   criterion ordering from a synthetic overlap case. This case
#'   currently requires `model = "GPCM"` and is intended for `method = "MML"`.
#' - `synthetic_bias_contract`: checks whether package bias tables and
#'   pairwise local comparisons satisfy the identities documented in the bias
#'   help workflow.
#' - `*_itercal_pair`: compares a baseline packaged dataset with its iterative
#'   recalibration counterpart to review fit stability, facet-measure
#'   alignment, and linking coverage together.
#'
#' The resulting object is intended as a reference-case check for package
#' behavior. It does not by itself establish
#' external validity against FACETS, ConQuest, or published calibration
#' studies, and it does not assume any familiarity with external table
#' numbering or printer layouts.
#' When specialized latent-regression omission or ConQuest-overlap package-side
#' cases are requested, `summary(bench)` prints preview rows from
#' `population_policy_checks` and `conquest_overlap_checks` alongside the
#' reference notes so the package-versus-external validation boundary remains
#' visible.
#'
#' @section Interpreting output:
#' - `overview`: one-row reference-case summary.
#' - `case_summary`: pass/warn/fail triage by reference case.
#' - `fit_runs`: fitted-run metadata (fit, precision tier, convergence, and
#'   latent-regression population-model/posterior-basis fields, including
#'   categorical-coding details when present).
#' - `design_checks`: exact design recovery checks for each dataset.
#' - `recovery_checks`: known-truth recovery metrics for the synthetic cases,
#'   including the latent-regression reference case.
#' - `bias_checks`: source-backed bias/local-measure identity checks.
#' - `pair_checks`: paired-dataset stability screens for the iterated cases.
#' - `linking_checks`: common-element audits for paired calibration datasets.
#' - `conquest_overlap_checks`: package-side checks for the
#'   ConQuest-overlap bundle/normalization/audit workflow; this remains a
#'   package-side check until actual ConQuest output tables are supplied.
#' - `population_policy_checks`: complete-case omission checks for population
#'   model benchmark fixtures.
#' - `source_profile`: source-backed rules used by the reference checks.
#'
#' @return An object of class `mfrm_reference_benchmark`.
#' @examples
#' \donttest{
#' bench <- reference_case_benchmark(
#'   cases = "synthetic_truth",
#'   method = "JML",
#'   maxit = 30
#' )
#' summary(bench)
#' }
#' @export
reference_case_benchmark <- function(cases = c(
                                       "synthetic_truth",
                                       "synthetic_latent_regression",
                                       "synthetic_bias_contract",
                                       "study1_itercal_pair",
                                       "study2_itercal_pair",
                                       "combined_itercal_pair"
                                     ),
                                     method = "MML",
                                     model = "RSM",
                                     quad_points = 7,
                                     maxit = 40,
                                     reltol = 1e-6,
                                     mml_engine = c("direct", "em", "hybrid")) {
  case_specs <- reference_benchmark_case_specs()
  selected_cases <- match.arg(as.character(cases), choices = case_specs$Case, several.ok = TRUE)
  case_specs <- case_specs[match(selected_cases, case_specs$Case), , drop = FALSE]
  dataset_specs <- reference_benchmark_dataset_specs()
  mml_engine <- tolower(match.arg(mml_engine))
  if (any(case_specs$CaseType == "gpcm_recovery")) {
    if (!identical(toupper(as.character(model)), "GPCM")) {
      stop(
        "The `synthetic_gpcm` benchmark case currently requires `model = \"GPCM\"`.",
        call. = FALSE
      )
    }
    if (!identical(toupper(as.character(method)), "MML")) {
      stop(
        "The `synthetic_gpcm` benchmark case is currently validated only for `method = \"MML\"`.",
        call. = FALSE
      )
    }
  }
  if (any(case_specs$CaseType == "conquest_overlap_dry_run")) {
    if (!identical(toupper(as.character(method)), "MML")) {
      stop(
        "The `synthetic_conquest_overlap_dry_run` benchmark case requires `method = \"MML\"`.",
        call. = FALSE
      )
    }
    if (!(toupper(as.character(model)) %in% c("RSM", "PCM"))) {
      stop(
        "The `synthetic_conquest_overlap_dry_run` benchmark case requires `model = \"RSM\"` or `model = \"PCM\"`.",
        call. = FALSE
      )
    }
  }

  fit_cache <- list()
  get_fit_obj <- function(dataset) {
    if (!dataset %in% names(fit_cache)) {
      fit_cache[[dataset]] <<- fit_reference_benchmark_dataset(
        dataset = dataset,
        method = method,
        model = model,
        quad_points = quad_points,
        maxit = maxit,
        reltol = reltol,
        mml_engine = mml_engine
      )
    }
    fit_cache[[dataset]]
  }

  fit_runs <- list()
  design_checks <- list()
  recovery_checks <- list()
  bias_checks <- list()
  pair_checks <- list()
  linking_checks <- list()
  conquest_overlap_checks <- list()
  population_policy_checks <- list()

  for (i in seq_len(nrow(case_specs))) {
    case_id <- as.character(case_specs$Case[i])
    case_type <- as.character(case_specs$CaseType[i])
    primary_dataset <- as.character(case_specs$PrimaryDataset[i])
    reference_dataset <- as.character(case_specs$ReferenceDataset[i] %||% NA_character_)

    primary_fit_obj <- get_fit_obj(primary_dataset)
    fit_runs[[length(fit_runs) + 1L]] <- collect_reference_fit_run(case_id, primary_fit_obj)
    spec_primary <- dataset_specs[dataset_specs$Dataset == primary_dataset, , drop = FALSE]
    design_checks[[length(design_checks) + 1L]] <- build_reference_design_checks(primary_fit_obj$design, spec_primary, case_id)

    if (identical(case_type, "truth_recovery")) {
      recovery_checks[[length(recovery_checks) + 1L]] <- build_truth_recovery_checks(case_id, primary_fit_obj)
    } else if (identical(case_type, "latent_recovery")) {
      recovery_checks[[length(recovery_checks) + 1L]] <- build_latent_recovery_checks(case_id, primary_fit_obj)
    } else if (identical(case_type, "gpcm_recovery")) {
      recovery_checks[[length(recovery_checks) + 1L]] <- build_gpcm_recovery_checks(case_id, primary_fit_obj)
    } else if (identical(case_type, "bias_contract")) {
      bias_checks[[length(bias_checks) + 1L]] <- build_bias_contract_checks(case_id, primary_fit_obj)
    } else if (identical(case_type, "conquest_overlap_dry_run")) {
      conquest_overlap_checks[[length(conquest_overlap_checks) + 1L]] <- build_conquest_overlap_dry_run_checks(case_id, primary_fit_obj)
    } else if (identical(case_type, "latent_omission_contract")) {
      population_policy_checks[[length(population_policy_checks) + 1L]] <- build_latent_omission_contract_checks(case_id, primary_fit_obj)
    } else {
      reference_fit_obj <- get_fit_obj(reference_dataset)
      fit_runs[[length(fit_runs) + 1L]] <- collect_reference_fit_run(case_id, reference_fit_obj)
      spec_reference <- dataset_specs[dataset_specs$Dataset == reference_dataset, , drop = FALSE]
      design_checks[[length(design_checks) + 1L]] <- build_reference_design_checks(reference_fit_obj$design, spec_reference, case_id)
      pair_checks[[length(pair_checks) + 1L]] <- build_pair_stability_checks(case_id, primary_fit_obj, reference_fit_obj)
      linking_checks[[length(linking_checks) + 1L]] <- build_linking_guideline_checks(case_id, primary_fit_obj, reference_fit_obj)
    }
  }

  fit_runs_tbl <- dplyr::bind_rows(fit_runs)
  design_checks_tbl <- dplyr::bind_rows(design_checks)
  recovery_checks_tbl <- dplyr::bind_rows(recovery_checks)
  bias_checks_tbl <- dplyr::bind_rows(bias_checks)
  pair_checks_tbl <- dplyr::bind_rows(pair_checks)
  linking_checks_tbl <- dplyr::bind_rows(linking_checks)
  conquest_overlap_checks_tbl <- dplyr::bind_rows(conquest_overlap_checks)
  population_policy_checks_tbl <- dplyr::bind_rows(population_policy_checks)
  source_profile_tbl <- reference_benchmark_source_profile()

  case_summary_tbl <- dplyr::bind_rows(lapply(seq_len(nrow(case_specs)), function(i) {
    summarize_reference_benchmark_case(
      case_id = as.character(case_specs$Case[i]),
      case_type = as.character(case_specs$CaseType[i]),
      fit_runs = fit_runs_tbl,
      design_checks = design_checks_tbl,
      recovery_checks = recovery_checks_tbl,
      pair_checks = pair_checks_tbl,
      bias_checks = bias_checks_tbl,
      linking_checks = linking_checks_tbl,
      conquest_overlap_checks = conquest_overlap_checks_tbl,
      population_policy_checks = population_policy_checks_tbl
    )
  }))

  overall <- tibble::tibble(
    Cases = nrow(case_summary_tbl),
    Fits = nrow(fit_runs_tbl),
    Pass = sum(case_summary_tbl$Status == "Pass", na.rm = TRUE),
    Warn = sum(case_summary_tbl$Status == "Warn", na.rm = TRUE),
    Fail = sum(case_summary_tbl$Status == "Fail", na.rm = TRUE),
    Method = as.character(method),
    Model = as.character(model),
    PassRate = if (nrow(case_summary_tbl) > 0) mean(case_summary_tbl$Status == "Pass", na.rm = TRUE) else NA_real_
  )

  notes <- c(
    "Synthetic truth checks compare recovered facet measures against known generating values from the package simulation design.",
    "ConQuest-overlap package-side checks cover only export/normalization/audit preparation; actual external ConQuest output is still required for external validation.",
    "Bias checks review package identities for observed-minus-expected averages, local measures, and pairwise Rasch-Welch contrasts.",
    "Pair stability checks review baseline and iterative-calibration packaged datasets using facet-measure alignment, fit deltas, reliability deltas, and common-element linking coverage.",
    "Use this reference check as package evidence, not as a substitute for external validation against commercial software or published studies."
  )
  if (!identical(toupper(method), "MML")) {
    notes <- c(
      notes,
      "Non-MML benchmark runs remain useful for stability auditing, but formal-inference expectations should be interpreted more conservatively."
    )
  }

  out <- list(
    overview = overall,
    summary = case_summary_tbl,
    table = fit_runs_tbl,
    fit_runs = fit_runs_tbl,
    case_summary = case_summary_tbl,
    design_checks = design_checks_tbl,
    recovery_checks = recovery_checks_tbl,
    bias_checks = bias_checks_tbl,
    pair_checks = pair_checks_tbl,
    linking_checks = linking_checks_tbl,
    conquest_overlap_checks = conquest_overlap_checks_tbl,
    population_policy_checks = population_policy_checks_tbl,
    source_profile = source_profile_tbl,
    settings = list(
      cases = selected_cases,
      method = method,
      model = model,
      mml_engine = if (identical(toupper(method), "MML")) mml_engine else NA_character_,
      intended_use = "internal_benchmark",
      external_validation = FALSE,
      quad_points = if (identical(toupper(method), "MML")) as.integer(quad_points) else NA_integer_,
      maxit = as.integer(maxit),
      reltol = reltol
    ),
    notes = notes
  )
  as_mfrm_bundle(out, "mfrm_reference_benchmark")
}

# Package-native reference benchmark helpers

reference_benchmark_dataset_specs <- function() {
  tibble::tibble(
    Dataset = c(
      "example_core", "example_bias",
      "study1", "study2", "combined",
      "study1_itercal", "study2_itercal", "combined_itercal",
      "synthetic_truth"
    ),
    Persons = c(48L, 48L, 307L, 206L, 307L, 307L, 206L, 307L, 36L),
    Raters = c(4L, 4L, 18L, 12L, 18L, 18L, 12L, 18L, 3L),
    Criteria = c(4L, 4L, 3L, 9L, 12L, 3L, 9L, 12L, 3L),
    Tasks = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, 4L),
    Rows = c(768L, 384L, 1842L, 3287L, 5129L, 1842L, 3341L, 5183L, 1296L),
    ScoreMin = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    ScoreMax = c(4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L)
  )
}

reference_benchmark_case_specs <- function() {
  tibble::tibble(
    Case = c("synthetic_truth", "synthetic_bias_contract", "study1_itercal_pair", "study2_itercal_pair", "combined_itercal_pair"),
    CaseType = c("truth_recovery", "bias_contract", "pair_stability", "pair_stability", "pair_stability"),
    PrimaryDataset = c("synthetic_truth", "synthetic_truth", "study1", "study2", "combined"),
    ReferenceDataset = c(NA_character_, NA_character_, "study1_itercal", "study2_itercal", "combined_itercal")
  )
}

synthetic_truth_targets <- function() {
  list(
    Rater = c(-0.4, 0, 0.4),
    Task = seq(-0.5, 0.5, length.out = 4),
    Criterion = c(-0.3, 0, 0.3)
  )
}

reference_benchmark_source_profile <- function() {
  tibble::tibble(
    RuleID = c(
      "bias_obs_exp_average",
      "bias_local_measure",
      "bias_pairwise_welch",
      "linking_common_elements"
    ),
    Domain = c("bias", "bias", "bias", "linking"),
    SourceLabel = c(
      "Facets Tutorial 3",
      "FACETS Table 14",
      "FACETS Table 14 / change log",
      "FACETS equating guidance"
    ),
    SourceURL = c(
      "https://www.winsteps.com/a/ftutorial3.pdf",
      "https://www.winsteps.com/facetman/table14.htm",
      "https://www.winsteps.com/facetman/table14.htm",
      "https://www.winsteps.com/facetman/equating.htm"
    ),
    Detail = c(
      "Observed-minus-expected average should equal (Observed Score - Expected Score) / Observed Count.",
      "Local target measure in a context should equal the target's overall measure plus the context-specific bias term.",
      "Pairwise local contrasts are reported with a Rasch-Welch t statistic and approximate degrees of freedom.",
      "A practical linking audit should confirm at least 5 common elements per linking facet when equating forms or datasets."
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
                                            reltol = 1e-6) {
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
  }
  fit <- suppressWarnings(do.call(fit_mfrm, fit_args))
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  list(
    dataset = dataset,
    data = cfg$data,
    fit = fit,
    diagnostics = diag,
    design = collect_reference_dataset_design(dataset, cfg$data)
  )
}

collect_reference_fit_run <- function(case_id, fit_obj) {
  fit <- fit_obj$fit
  diag <- fit_obj$diagnostics
  design <- fit_obj$design
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

summarize_reference_benchmark_case <- function(case_id, case_type, fit_runs, design_checks, recovery_checks, pair_checks, bias_checks, linking_checks) {
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

  statuses <- c(
    case_fit$Converged == TRUE,
    case_statuses(case_design),
    case_statuses(case_recovery),
    case_statuses(case_pairs),
    case_statuses(case_bias),
    case_statuses(case_link)
  )
  normalized <- ifelse(statuses %in% c(TRUE, "Pass"), "Pass",
                       ifelse(statuses %in% c("Warn"), "Warn",
                              ifelse(statuses %in% c(FALSE, "Fail"), "Fail", NA_character_)))
  missing_expected_checks <- switch(
    case_type,
    truth_recovery = nrow(case_recovery) == 0L,
    bias_contract = nrow(case_bias) == 0L,
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
  } else if (identical(case_type, "bias_contract")) {
    if (nrow(case_bias) > 0) {
      paste0(
        "Max bias-identity error = ",
        formatC(max(case_bias$MaxError, na.rm = TRUE), format = "f", digits = 6)
      )
    } else {
      "No bias-contract checks were produced."
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
    StabilityChecks = nrow(case_pairs),
    KeySignal = key_signal
  )
}

#' Benchmark packaged reference cases
#'
#' @param cases Reference cases to run. Defaults to all package-native
#'   benchmark cases.
#' @param method Estimation method passed to [fit_mfrm()]. Defaults to `"MML"`.
#' @param model Model family passed to [fit_mfrm()]. Defaults to `"RSM"`.
#' @param quad_points Quadrature points for `method = "MML"`.
#' @param maxit Maximum optimizer iterations passed to [fit_mfrm()].
#' @param reltol Convergence tolerance passed to [fit_mfrm()].
#'
#' @details
#' This function audits `mfrmr` against the package's curated internal
#' benchmark cases in three ways:
#' - `synthetic_truth`: checks whether recovered facet measures align with the
#'   known generating values from the package's internal synthetic design.
#' - `synthetic_bias_contract`: checks whether package-native bias tables and
#'   pairwise local comparisons satisfy the identities documented in the bias
#'   help workflow.
#' - `*_itercal_pair`: compares a baseline packaged dataset with its iterative
#'   recalibration counterpart to review fit stability, facet-measure
#'   alignment, and linking coverage together.
#'
#' The resulting object is intended as an internal benchmark harness for
#' package QA and regression auditing. It does not by itself establish
#' external validity against FACETS, ConQuest, or published calibration
#' studies, and it does not assume any familiarity with external table
#' numbering or printer layouts.
#'
#' @section Interpreting output:
#' - `overview`: one-row internal-benchmark summary.
#' - `case_summary`: pass/warn/fail triage by reference case.
#' - `fit_runs`: fitted-run metadata (fit, precision tier, convergence).
#' - `design_checks`: exact design recovery checks for each dataset.
#' - `recovery_checks`: known-truth recovery metrics for the internal synthetic
#'   case.
#' - `bias_checks`: source-backed bias/local-measure identity checks.
#' - `pair_checks`: paired-dataset stability screens for the iterated cases.
#' - `linking_checks`: common-element audits for paired calibration datasets.
#' - `source_profile`: source-backed rules that define the internal benchmark
#'   contract.
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
                                       "synthetic_bias_contract",
                                       "study1_itercal_pair",
                                       "study2_itercal_pair",
                                       "combined_itercal_pair"
                                     ),
                                     method = "MML",
                                     model = "RSM",
                                     quad_points = 7,
                                     maxit = 40,
                                     reltol = 1e-6) {
  case_specs <- reference_benchmark_case_specs()
  selected_cases <- match.arg(as.character(cases), choices = case_specs$Case, several.ok = TRUE)
  case_specs <- case_specs[match(selected_cases, case_specs$Case), , drop = FALSE]
  dataset_specs <- reference_benchmark_dataset_specs()

  fit_cache <- list()
  get_fit_obj <- function(dataset) {
    if (!dataset %in% names(fit_cache)) {
      fit_cache[[dataset]] <<- fit_reference_benchmark_dataset(
        dataset = dataset,
        method = method,
        model = model,
        quad_points = quad_points,
        maxit = maxit,
        reltol = reltol
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
    } else if (identical(case_type, "bias_contract")) {
      bias_checks[[length(bias_checks) + 1L]] <- build_bias_contract_checks(case_id, primary_fit_obj)
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
      linking_checks = linking_checks_tbl
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
    "Synthetic truth checks compare recovered facet measures against known generating values from the package's internal simulation design.",
    "Bias-contract checks audit package-native identities for observed-minus-expected averages, local measures, and pairwise Rasch-Welch contrasts.",
    "Pair stability checks review baseline and iterative-calibration packaged datasets using facet-measure alignment, fit deltas, reliability deltas, and common-element linking coverage.",
    "Use this benchmark as an internal regression and contract audit, not as a substitute for external validation against commercial software or published studies."
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
    source_profile = source_profile_tbl,
    settings = list(
      cases = selected_cases,
      method = method,
      model = model,
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

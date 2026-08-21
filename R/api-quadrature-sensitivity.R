# Explicit same-data quadrature sensitivity for the bounded GPCM-MML route.

mfrmr_gqs_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

mfrmr_gqs_capture_conditions <- function(expression) {
  warnings <- character()
  messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
    },
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
    }
  )
  list(
    value = value,
    warnings = unique(warnings),
    messages = unique(messages)
  )
}

mfrmr_gqs_validate <- function(fit, data, quad_points,
                               theta_range, theta_points) {
  mfrmr_gqs_assert(
    inherits(fit, "mfrm_fit"),
    "`fit` must be an `mfrm_fit` object returned by `fit_mfrm()`."
  )
  model <- toupper(as.character(fit$config$model %||% "")[1L])
  method <- public_mfrm_method_label(fit$config$method %||% "")
  mfrmr_gqs_assert(
    identical(model, "GPCM") && identical(method, "MML"),
    "`gpcm_mml_quadrature_sensitivity()` requires a GPCM MML fit."
  )
  mfrmr_gqs_assert(
    is.data.frame(data) && nrow(data) > 0L,
    "`data` must be the non-empty data.frame used to create `fit`."
  )
  mfrmr_gqs_assert(
    length(fit$config$interaction_specs %||% list()) == 0L,
    paste(
      "Quadrature sensitivity for GPCM facet interactions is not yet",
      "implemented. Refit the additive GPCM or compare interaction fits",
      "explicitly outside this helper."
    )
  )
  population <- fit$population %||% list()
  mfrmr_gqs_assert(
    !isTRUE(population$active) || identical(
      as.character(population$source %||% ""),
      "gpcm_mml_default_identification"
    ),
    paste(
      "Quadrature sensitivity for user-supplied latent-regression",
      "population models is not yet implemented. Compare those fits",
      "explicitly with the same person-data contract."
    )
  )

  integer_tolerance <- sqrt(.Machine$double.eps)
  mfrmr_gqs_assert(
    is.numeric(quad_points) && length(quad_points) >= 2L &&
      all(is.finite(quad_points)) && all(quad_points >= 1) &&
      all(abs(quad_points - round(quad_points)) <= integer_tolerance),
    "`quad_points` must contain at least two finite positive integers."
  )
  quad_points <- sort(unique(as.integer(quad_points)))
  mfrmr_gqs_assert(
    length(quad_points) >= 2L,
    "`quad_points` must contain at least two distinct values."
  )
  reference_nodes <- as.integer(
    fit$config$estimation_control$quad_points %||%
      fit$config$quad_points %||% NA_integer_
  )[1L]
  mfrmr_gqs_assert(
    is.finite(reference_nodes) && reference_nodes %in% quad_points,
    paste0(
      "`quad_points` must include the reference fit's quadrature count (q=",
      reference_nodes, ")."
    )
  )

  theta_range <- as.numeric(theta_range)
  mfrmr_gqs_assert(
    length(theta_range) == 2L && all(is.finite(theta_range)) &&
      theta_range[1L] < theta_range[2L],
    "`theta_range` must contain two increasing finite values."
  )
  mfrmr_gqs_assert(
    is.numeric(theta_points) && length(theta_points) == 1L &&
      is.finite(theta_points) && theta_points >= 21L &&
      abs(theta_points - round(theta_points)) <= integer_tolerance,
    "`theta_points` must be one integer of at least 21."
  )

  replay <- fit$config$replay_inputs %||% list()
  mfrmr_gqs_assert(
    is.list(replay) &&
      all(c("person", "facets", "score", "model", "method") %in%
            names(replay)),
    paste(
      "The fit does not retain the current replay-input contract needed for",
      "same-model refits. Refit it with the current mfrmr version first."
    )
  )

  list(
    quad_points = quad_points,
    reference_nodes = reference_nodes,
    theta_range = theta_range,
    theta_points = as.integer(theta_points),
    replay = replay
  )
}

mfrmr_gqs_refit_arguments <- function(fit, data, nodes) {
  replay <- fit$config$replay_inputs
  config <- fit$config %||% list()
  control <- config$estimation_control %||% list()

  arguments <- list(
    data = data,
    person = as.character(replay$person),
    facets = as.character(replay$facets),
    score = as.character(replay$score),
    rating_min = replay$rating_min %||% config$rating_min,
    rating_max = replay$rating_max %||% config$rating_max,
    weight = replay$weight,
    keep_original = isTRUE(replay$keep_original),
    missing_codes = replay$missing_codes,
    model = "GPCM",
    method = "MML",
    step_facet = replay$step_facet %||% config$step_facet,
    slope_facet = replay$slope_facet %||% config$slope_facet,
    facet_interactions = replay$facet_interactions,
    min_obs_per_interaction = as.numeric(
      replay$min_obs_per_interaction %||% 10
    ),
    interaction_policy = as.character(
      replay$interaction_policy %||% "warn"
    ),
    anchors = replay$anchors,
    group_anchors = replay$group_anchors,
    noncenter_facet = as.character(
      replay$noncenter_facet %||% "Person"
    ),
    dummy_facets = replay$dummy_facets,
    positive_facets = replay$positive_facets,
    anchor_policy = as.character(replay$anchor_policy %||% "warn"),
    min_common_anchors = as.integer(replay$min_common_anchors %||% 5L),
    min_obs_per_element = as.numeric(replay$min_obs_per_element %||% 30),
    min_obs_per_category = as.numeric(
      replay$min_obs_per_category %||% 10
    ),
    quad_points = as.integer(nodes),
    maxit = as.integer(replay$maxit %||% control$maxit %||% 400L),
    reltol = as.numeric(replay$reltol %||% control$reltol %||% 1e-9),
    optimizer = as.character(
      replay$optimizer %||% control$optimizer_requested %||% "auto"
    ),
    mml_engine = as.character(
      replay$mml_engine %||% control$mml_engine_requested %||% "direct"
    ),
    population_policy = as.character(
      replay$population_policy %||% fit$population$policy %||% "error"
    ),
    facet_shrinkage = as.character(replay$facet_shrinkage %||% "none"),
    facet_prior_sd = replay$facet_prior_sd,
    shrink_person = isTRUE(replay$shrink_person),
    attach_diagnostics = FALSE,
    gpcm_mml_identification = as.character(
      replay$gpcm_mml_identification %||%
        config$gpcm_mml_identification %||% "free_population"
    )
  )
  arguments
}

mfrmr_gqs_canonical_prepared_data <- function(fit) {
  prepared <- as.data.frame(fit$prep$data %||% data.frame(),
                            stringsAsFactors = FALSE)
  facet_names <- as.character(fit$config$facet_names %||% character())
  columns <- c("Person", facet_names, "score_k", "Weight")
  mfrmr_gqs_assert(
    nrow(prepared) > 0L && all(columns %in% names(prepared)),
    "A fit lacked the prepared response columns required for same-data review."
  )
  prepared <- prepared[, columns, drop = FALSE]
  for (name in names(prepared)) {
    if (is.factor(prepared[[name]])) {
      prepared[[name]] <- as.character(prepared[[name]])
    }
  }
  ordering_values <- lapply(prepared, function(value) {
    if (is.character(value)) ifelse(is.na(value), "<NA>", value) else value
  })
  ordering <- do.call(
    order,
    c(ordering_values, list(na.last = TRUE, method = "radix"))
  )
  prepared <- prepared[ordering, , drop = FALSE]
  rownames(prepared) <- NULL
  prepared
}

mfrmr_gqs_same_prepared_data <- function(reference, candidate) {
  reference_data <- mfrmr_gqs_canonical_prepared_data(reference)
  candidate_data <- mfrmr_gqs_canonical_prepared_data(candidate)
  identical(names(reference_data), names(candidate_data)) &&
    isTRUE(all.equal(
      reference_data,
      candidate_data,
      check.attributes = FALSE,
      tolerance = sqrt(.Machine$double.eps)
    ))
}

mfrmr_gqs_design_probabilities <- function(fit, theta_range, theta_points) {
  theta <- seq(theta_range[1L], theta_range[2L], length.out = theta_points)
  structure <- information_build_step_structure(fit, "GPCM")
  prepared <- as.data.frame(fit$prep$data, stringsAsFactors = FALSE)
  facet_names <- as.character(fit$config$facet_names %||% character())
  mfrmr_gqs_assert(
    length(facet_names) > 0L && all(facet_names %in% names(prepared)),
    "The fitted facet columns were unavailable for probability comparison."
  )

  cells <- unique(prepared[, facet_names, drop = FALSE])
  ordering_values <- lapply(cells, function(value) as.character(value))
  ordering <- do.call(
    order,
    c(ordering_values, list(na.last = TRUE, method = "radix"))
  )
  cells <- cells[ordering, , drop = FALSE]
  rownames(cells) <- NULL

  facets <- as.data.frame(fit$facets$others %||% data.frame(),
                          stringsAsFactors = FALSE)
  mfrmr_gqs_assert(
    nrow(facets) > 0L &&
      all(c("Facet", "Level", "Estimate") %in% names(facets)),
    "Finite non-Person facet estimates are required for probability review."
  )
  estimate_key <- paste(facets$Facet, facets$Level, sep = "\r")
  mfrmr_gqs_assert(
    !anyDuplicated(estimate_key),
    "The fitted facet table contained duplicate Facet/Level rows."
  )
  estimates <- stats::setNames(as.numeric(facets$Estimate), estimate_key)
  signs <- fit$config$facet_signs %||%
    stats::setNames(rep(-1, length(facet_names)), facet_names)
  signs <- as.numeric(signs[facet_names])
  signs[!is.finite(signs)] <- -1

  offsets <- numeric(nrow(cells))
  for (facet_index in seq_along(facet_names)) {
    facet <- facet_names[facet_index]
    keys <- paste(facet, as.character(cells[[facet]]), sep = "\r")
    values <- as.numeric(estimates[keys])
    mfrmr_gqs_assert(
      all(is.finite(values)),
      "The fitted facet table did not cover every observed design cell."
    )
    offsets <- offsets + signs[facet_index] * values
  }

  step_index <- match(
    as.character(cells[[structure$step_facet]]),
    structure$step_levels
  )
  slope_index <- match(
    as.character(cells[[structure$slope_facet]]),
    structure$slope_levels
  )
  mfrmr_gqs_assert(
    all(is.finite(step_index)) && all(is.finite(slope_index)),
    "Step or slope levels were not aligned with the observed design cells."
  )

  probability_blocks <- lapply(seq_len(nrow(cells)), function(cell_index) {
    probabilities <- structure$compute(
      theta + offsets[cell_index],
      rep(step_index[cell_index], length(theta)),
      rep(slope_index[cell_index], length(theta))
    )
    mfrmr_gqs_assert(
      is.matrix(probabilities) && all(is.finite(probabilities)) &&
        max(abs(rowSums(probabilities) - 1)) < 1e-10,
      "A fitted design cell produced invalid category probabilities."
    )
    as.numeric(t(probabilities))
  })

  list(
    cells = cells,
    theta = theta,
    categories = structure$categories,
    values = unlist(probability_blocks, use.names = FALSE)
  )
}

mfrmr_gqs_probability_difference <- function(reference, candidate) {
  mfrmr_gqs_assert(
    identical(names(reference$cells), names(candidate$cells)) &&
      isTRUE(all.equal(reference$cells, candidate$cells,
                       check.attributes = FALSE)) &&
      isTRUE(all.equal(reference$theta, candidate$theta,
                       tolerance = 1e-14)) &&
      identical(reference$categories, candidate$categories) &&
      length(reference$values) == length(candidate$values),
    "Probability grids could not be aligned across quadrature fits."
  )
  max(abs(reference$values - candidate$values))
}

mfrmr_gqs_raw_information <- function(fit) {
  covariance <- compute_mml_parameter_covariance(fit)
  hessian <- covariance$hessian
  eigenvalues <- if (is.matrix(hessian) && nrow(hessian) > 0L) {
    eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
  } else {
    numeric()
  }
  minimum <- if (length(eigenvalues) > 0L) min(eigenvalues) else NA_real_
  maximum <- if (length(eigenvalues) > 0L) max(eigenvalues) else NA_real_
  condition_number <- if (is.finite(minimum) && minimum > 0 &&
                          is.finite(maximum)) {
    maximum / minimum
  } else {
    NA_real_
  }

  slopes <- as.data.frame(fit$slopes %||% data.frame(),
                          stringsAsFactors = FALSE)
  slope_values <- as.numeric(
    slopes$OptimizerEstimate %||% slopes$Estimate %||% numeric()
  )
  slope_se <- rep(NA_real_, length(slope_values))
  slope_slice <- covariance$param_slices$log_slopes %||% integer()
  if (covariance$status %in% c("ok", "regularized") &&
      is.matrix(covariance$cov) && length(slope_values) > 0L &&
      length(slope_slice) == length(slope_values) - 1L) {
    jacobian <- sum_zero_jacobian(length(slope_values))
    log_covariance <- jacobian %*%
      covariance$cov[slope_slice, slope_slice, drop = FALSE] %*%
      t(jacobian)
    slope_se <- slope_values * sqrt(pmax(diag(log_covariance), 0))
  }

  population_sd <- if (isTRUE(fit$population$active)) {
    sqrt(as.numeric(fit$population$sigma2 %||% NA_real_)[1L])
  } else {
    1
  }
  population_sd_se <- NA_real_
  sigma_slice <- covariance$param_slices$log_sigma2 %||% integer()
  if (covariance$status %in% c("ok", "regularized") &&
      is.matrix(covariance$cov) && length(sigma_slice) == 1L &&
      is.finite(population_sd)) {
    population_sd_se <- 0.5 * population_sd *
      sqrt(pmax(covariance$cov[sigma_slice, sigma_slice], 0))
  }

  list(
    covariance_status = as.character(covariance$status),
    rank = as.integer(covariance$rank),
    dimension = as.integer(length(fit$opt$par %||% numeric())),
    minimum_eigenvalue = minimum,
    maximum_eigenvalue = maximum,
    condition_number = condition_number,
    slope_values = slope_values,
    slope_se = slope_se,
    population_sd = population_sd,
    population_sd_se = population_sd_se
  )
}

mfrmr_gqs_max_named_difference <- function(reference, candidate, label) {
  mfrmr_gqs_assert(
    setequal(names(reference), names(candidate)),
    paste0("The ", label, " parameter identities changed across refits.")
  )
  candidate <- candidate[names(reference)]
  max(abs(reference - candidate))
}

mfrmr_gqs_extract_run <- function(fit, nodes, theta_range, theta_points) {
  information <- mfrmr_gqs_raw_information(fit)
  readiness <- as.data.frame(fit$readiness$fit %||% data.frame(),
                             stringsAsFactors = FALSE)
  mfrmr_gqs_assert(
    nrow(readiness) == 1L,
    "A quadrature fit did not retain exactly one fit-readiness row."
  )
  slopes <- as.data.frame(fit$slopes, stringsAsFactors = FALSE)
  slope_labels <- as.character(slopes$SlopeFacet)
  mfrmr_gqs_assert(
    length(slope_labels) == length(information$slope_values) &&
      !anyDuplicated(slope_labels),
    "Slope labels were not aligned with the observed-information coordinates."
  )

  list(
    run = data.frame(
      Nodes = as.integer(nodes),
      Persons = as.integer(fit$summary$Persons[1L] %||%
                             fit$prep$n_person),
      ResponseRows = as.integer(fit$summary$N[1L] %||% fit$prep$n_obs),
      Npar = as.integer(fit$summary$Npar[1L]),
      NegativeLogLikelihood = -as.numeric(fit$summary$LogLik[1L]),
      NLLPerPerson = -as.numeric(fit$summary$LogLik[1L]) /
        as.numeric(fit$summary$Persons[1L] %||% fit$prep$n_person),
      TerminalGradientSupNorm = as.numeric(
        fit$summary$TerminalGradientSupNorm[1L]
      ),
      OptimizerConvergenceCode = as.integer(fit$opt$convergence %||%
                                              NA_integer_),
      EstimationConverged = isTRUE(
        fit$population$estimation_converged %||%
          identical(as.integer(fit$opt$convergence), 0L)
      ),
      CovarianceStatus = information$covariance_status,
      HessianRank = information$rank,
      HessianDimension = information$dimension,
      MinimumEigenvalue = information$minimum_eigenvalue,
      MaximumEigenvalue = information$maximum_eigenvalue,
      HessianConditionNumber = information$condition_number,
      PopulationSD = information$population_sd,
      RawPopulationSDSE = information$population_sd_se,
      FitReadiness = as.character(readiness$FitReadiness[1L]),
      InferenceReady = isTRUE(readiness$InferenceReady[1L]),
      ReadinessReasons = as.character(readiness$ReasonCodes[1L]),
      stringsAsFactors = FALSE
    ),
    slopes = data.frame(
      Nodes = as.integer(nodes),
      SlopeFacet = slope_labels,
      OptimizerEstimate = information$slope_values,
      RawObservedInformationSE = information$slope_se,
      SEStatus = information$covariance_status,
      PublicSEEligible = as.logical(slopes$SEEligible %||%
                                      rep(FALSE, nrow(slopes))),
      stringsAsFactors = FALSE
    ),
    slope_parameters = stats::setNames(
      information$slope_values,
      slope_labels
    ),
    slope_se = stats::setNames(information$slope_se, slope_labels),
    population_sd = information$population_sd,
    population_sd_se = information$population_sd_se,
    probabilities = mfrmr_gqs_design_probabilities(
      fit, theta_range, theta_points
    )
  )
}

mfrmr_gqs_condition_rows <- function(nodes, capture) {
  rows <- list()
  if (length(capture$warnings) > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      Nodes = as.integer(nodes),
      Type = "warning",
      Message = capture$warnings,
      stringsAsFactors = FALSE
    )
  }
  if (length(capture$messages) > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      Nodes = as.integer(nodes),
      Type = "message",
      Message = capture$messages,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(data.frame(
      Nodes = integer(), Type = character(), Message = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

#' Review GPCM-MML sensitivity to the quadrature grid
#'
#' Refit one bounded GPCM-MML model to the same supplied response data at two
#' or more Gauss--Hermite quadrature counts. The original fit is reused at its
#' own quadrature count; all other fits reuse its stored model, identification,
#' anchor, optimizer, and population settings.
#'
#' @param fit A GPCM MML `mfrm_fit` returned by [fit_mfrm()].
#' @param data The original response data.frame used to create `fit`. Prepared
#'   response rows are compared semantically after refitting; row order may
#'   differ, but changed observations fail closed.
#' @param quad_points At least two distinct positive integers, including the
#'   quadrature count stored on `fit`. A common choice is `c(31, 41)`.
#' @param theta_range Two finite values defining the common ability grid used
#'   for fitted category-probability comparison.
#' @param theta_points Number of common-grid ability points; at least 21.
#'
#' @details
#' This is an explicit refit diagnostic: neither [summary.mfrm_fit()] nor
#' [diagnose_mfrm()] invokes it automatically. It reports continuous changes
#' rather than classifying a fit as quadrature-stable or unstable. In
#' particular, it does not set a practical cutoff, promote slope standard
#' errors, or override the fit-readiness record.
#'
#' The probability comparison evaluates every observed combination of
#' non-Person facet levels on the same theta grid. It therefore exercises the
#' complete-predictor GPCM slope action, including additive Rater and Criterion
#' locations. The first version fails closed for fitted facet interactions
#' rather than approximating their contribution.
#'
#' Raw slope and population-SD standard errors are computed from each local
#' observed-information Hessian for diagnostic comparison only. The public
#' parameter-level `SEEligible` state remains unchanged.
#'
#' @return An object of class `mfrm_quadrature_sensitivity` containing:
#' - `summary`: one comparison row per quadrature grid relative to the original
#'   fit;
#' - `runs`: likelihood, gradient, curvature, population-scale, and readiness
#'   details for each fit;
#' - `slopes`: relative-slope estimates and raw diagnostic SEs;
#' - `conditions`: warnings and messages emitted by the explicit refits;
#' - `fits`: the reference and refitted `mfrm_fit` objects;
#' - `settings` and `notes`: the fixed comparison contract and interpretation
#'   boundary.
#'
#' The object supports `print()`, `summary()`, `as.data.frame()`, and the
#' existing `apa_table()` list route (for example `apa_table(out)`).
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [apa_table()]
#' @examples
#' \donttest{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", model = "GPCM",
#'   step_facet = "Criterion", slope_facet = "Criterion",
#'   quad_points = 31
#' )
#' sensitivity <- gpcm_mml_quadrature_sensitivity(
#'   fit, toy, quad_points = c(31, 41)
#' )
#' summary(sensitivity)
#' # Preserve small numerical differences when preparing the sensitivity table.
#' apa_table(sensitivity, digits = 5)
#' }
#' @export
gpcm_mml_quadrature_sensitivity <- function(
    fit,
    data,
    quad_points = c(31L, 41L),
    theta_range = c(-4, 4),
    theta_points = 161L) {
  contract <- mfrmr_gqs_validate(
    fit, data, quad_points, theta_range, theta_points
  )

  fits <- list()
  extracted <- list()
  condition_rows <- list()
  reference_name <- paste0("q", contract$reference_nodes)
  fits[[reference_name]] <- fit

  for (nodes in contract$quad_points) {
    name <- paste0("q", nodes)
    if (nodes == contract$reference_nodes) {
      candidate <- fit
      capture <- list(warnings = character(), messages = character())
    } else {
      capture <- mfrmr_gqs_capture_conditions(
        do.call(
          fit_mfrm,
          mfrmr_gqs_refit_arguments(fit, data, nodes)
        )
      )
      candidate <- capture$value
      mfrmr_gqs_assert(
        mfrmr_gqs_same_prepared_data(fit, candidate),
        paste0(
          "The q=", nodes,
          " refit did not reproduce the reference fit's prepared response rows."
        )
      )
      fits[[name]] <- candidate
    }
    extracted[[name]] <- mfrmr_gqs_extract_run(
      candidate,
      nodes,
      contract$theta_range,
      contract$theta_points
    )
    condition_rows[[name]] <- mfrmr_gqs_condition_rows(nodes, capture)
  }

  reference <- extracted[[reference_name]]
  comparison_rows <- lapply(contract$quad_points, function(nodes) {
    candidate <- extracted[[paste0("q", nodes)]]
    slope_difference <- mfrmr_gqs_max_named_difference(
      reference$slope_parameters,
      candidate$slope_parameters,
      "slope"
    )
    slope_se_difference <- if (
      all(is.finite(reference$slope_se)) &&
        all(is.finite(candidate$slope_se))
    ) {
      mfrmr_gqs_max_named_difference(
        reference$slope_se,
        candidate$slope_se,
        "slope-SE"
      )
    } else {
      NA_real_
    }
    reference_run <- reference$run
    candidate_run <- candidate$run
    nll_change_per_person <- candidate_run$NLLPerPerson -
      reference_run$NLLPerPerson
    data.frame(
      ReferenceNodes = contract$reference_nodes,
      Nodes = as.integer(nodes),
      IsReference = nodes == contract$reference_nodes,
      NLLChangePerPerson = nll_change_per_person,
      NLLAbsChangePerPerson = abs(nll_change_per_person),
      SlopeMaxAbsChange = slope_difference,
      RawSlopeSEMaxAbsChange = slope_se_difference,
      PopulationSDAbsChange = abs(
        candidate$population_sd - reference$population_sd
      ),
      RawPopulationSDSEAbsChange = abs(
        candidate$population_sd_se - reference$population_sd_se
      ),
      ProbabilityMaxAbsChange = mfrmr_gqs_probability_difference(
        reference$probabilities,
        candidate$probabilities
      ),
      stringsAsFactors = FALSE
    )
  })

  comparison <- do.call(rbind, comparison_rows)
  runs <- do.call(rbind, lapply(extracted, `[[`, "run"))
  slopes <- do.call(rbind, lapply(extracted, `[[`, "slopes"))
  conditions <- do.call(rbind, condition_rows)
  rownames(comparison) <- NULL
  rownames(runs) <- NULL
  rownames(slopes) <- NULL
  rownames(conditions) <- NULL

  out <- list(
    summary = comparison,
    runs = runs,
    slopes = slopes,
    conditions = conditions,
    fits = fits,
    settings = list(
      reference_nodes = contract$reference_nodes,
      quad_points = contract$quad_points,
      theta_range = contract$theta_range,
      theta_points = contract$theta_points,
      response_contract = "same_semantically_canonical_prepared_rows",
      probability_contract =
        "all_observed_additive_facet_cells_on_common_theta_grid",
      gpcm_slope_action = as.character(
        fit$config$gpcm_slope_action %||% "complete_adjacent_category_predictor"
      ),
      stability_threshold_frozen = FALSE,
      readiness_effect = "none_diagnostic_only"
    ),
    notes = c(
      paste(
        "Continuous q-grid differences are reported without an automatic",
        "stable/unstable classification."
      ),
      paste(
        "Raw observed-information SEs are diagnostic and do not change",
        "parameter-level SE eligibility."
      ),
      paste(
        "Local positive curvature cannot establish absence of a remote",
        "likelihood-recession path or promote inference readiness."
      )
    )
  )
  class(out) <- c("mfrm_quadrature_sensitivity", "list")
  out
}

#' @export
as.data.frame.mfrm_quadrature_sensitivity <- function(x, row.names = NULL,
                                                       optional = FALSE, ...) {
  as.data.frame(x$summary, row.names = row.names, optional = optional,
                stringsAsFactors = FALSE)
}

#' @export
summary.mfrm_quadrature_sensitivity <- function(object, ...) {
  comparison <- as.data.frame(object$summary, stringsAsFactors = FALSE)
  nonreference <- comparison[!comparison$IsReference, , drop = FALSE]
  max_or_na <- function(value) {
    value <- as.numeric(value)
    value <- value[is.finite(value)]
    if (length(value) == 0L) NA_real_ else max(value)
  }
  runs <- as.data.frame(object$runs, stringsAsFactors = FALSE)
  overview <- data.frame(
    ReferenceNodes = as.integer(object$settings$reference_nodes),
    ComparedGrids = nrow(runs),
    AllEstimationConverged = all(runs$EstimationConverged),
    AllHessiansFullRank = all(
      runs$HessianRank == runs$HessianDimension &
        is.finite(runs$MinimumEigenvalue) & runs$MinimumEigenvalue > 0
    ),
    MaxNLLAbsChangePerPerson = max_or_na(
      nonreference$NLLAbsChangePerPerson
    ),
    MaxSlopeAbsChange = max_or_na(nonreference$SlopeMaxAbsChange),
    MaxPopulationSDAbsChange = max_or_na(
      nonreference$PopulationSDAbsChange
    ),
    MaxProbabilityAbsChange = max_or_na(
      nonreference$ProbabilityMaxAbsChange
    ),
    StabilityClassification = "not_assigned_continuous_evidence_only",
    ReadinessEffect = as.character(object$settings$readiness_effect),
    stringsAsFactors = FALSE
  )
  out <- list(
    overview = overview,
    comparison = comparison,
    runs = runs,
    slopes = as.data.frame(object$slopes, stringsAsFactors = FALSE),
    conditions = as.data.frame(object$conditions, stringsAsFactors = FALSE),
    notes = object$notes
  )
  class(out) <- c("summary.mfrm_quadrature_sensitivity", "list")
  out
}

#' @export
print.mfrm_quadrature_sensitivity <- function(x, digits = 5L, ...) {
  cat("GPCM-MML quadrature sensitivity\n")
  cat("Reference grid: q=", x$settings$reference_nodes, "\n", sep = "")
  display <- as.data.frame(x$summary, stringsAsFactors = FALSE)
  numeric <- vapply(display, is.numeric, logical(1))
  display[numeric] <- lapply(display[numeric], round, digits = digits)
  print(display, row.names = FALSE)
  cat(
    "No automatic stability classification or readiness change is applied.\n"
  )
  invisible(x)
}

#' @export
print.summary.mfrm_quadrature_sensitivity <- function(x, digits = 5L, ...) {
  cat("GPCM-MML quadrature sensitivity summary\n")
  overview <- x$overview
  numeric <- vapply(overview, is.numeric, logical(1))
  overview[numeric] <- lapply(overview[numeric], round, digits = digits)
  print(overview, row.names = FALSE)
  cat("\nGrid comparisons\n")
  comparison <- x$comparison
  numeric <- vapply(comparison, is.numeric, logical(1))
  comparison[numeric] <- lapply(comparison[numeric], round, digits = digits)
  print(comparison, row.names = FALSE)
  invisible(x)
}

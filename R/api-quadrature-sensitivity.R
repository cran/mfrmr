# Explicit same-data quadrature sensitivity for ordered MML models.

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
                               theta_range, theta_points,
                               allowed_models = "GPCM") {
  mfrmr_gqs_assert(
    inherits(fit, "mfrm_fit"),
    "`fit` must be an `mfrm_fit` object returned by `fit_mfrm()`."
  )
  model <- toupper(as.character(fit$config$model %||% "")[1L])
  method <- public_mfrm_method_label(fit$config$method %||% "")
  allowed_models <- unique(toupper(as.character(allowed_models)))
  scope_message <- if (identical(allowed_models, "GPCM")) {
    "`gpcm_mml_quadrature_sensitivity()` requires a GPCM MML fit."
  } else {
    paste0(
      "`mml_quadrature_sensitivity()` requires an RSM, PCM, or ",
      "GPCM MML fit."
    )
  }
  mfrmr_gqs_assert(model %in% allowed_models && identical(method, "MML"),
                   scope_message)
  mfrmr_gqs_assert(
    is.data.frame(data) && nrow(data) > 0L,
    "`data` must be the non-empty data.frame used to create `fit`."
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
    model = model,
    replay = replay
  )
}

mfrmr_gqs_refit_arguments <- function(fit, data, nodes) {
  replay <- fit$config$replay_inputs
  if (mfrm_has_product_slopes(fit)) {
    replay$mml_integration <- fit$config$estimation_control$mml_integration
    replay$mml_accuracy_control <- NULL
    if (mfrmr_adaptive_integration(fit$config)) {
      policy <- fit$config$estimation_control$gpcm_mml_start %||% "neutral"
      if (identical(policy, "explicit")) {
        record <- fit$opt$mml_initialization
        starts <- lapply(record$attempts, `[[`, "start")
        mfrmr_gqs_assert(identical(record$policy, "explicit") && length(starts) >= 2L &&
          identical(names(starts), record$table$Start) && identical(replay$gpcm_mml_start, starts),
          "The saved explicit initialization vectors and replay settings must agree.")
        policy <- starts
      }
      mfrmr_gqs_assert(is.null(replay$gpcm_mml_start) || identical(replay$gpcm_mml_start, policy),
        "The saved adaptive initialization and replay settings must agree.")
      replay$gpcm_mml_start <- policy
    }
    replay$package_version <- NULL
    replay$data <- data
    replay$quad_points <- as.integer(nodes)
    return(replay)
  }
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
    model = as.character(replay$model %||% config$model),
    method = "MML",
    mml_starts = replay$mml_starts,
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
    centering_scope = config$centering_scope %||% "remainder",
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
    mml_integration = as.character(control$mml_integration %||% "fixed"),
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
  c(arguments, mfrmr_gqs_population_arguments(fit))
}

mfrmr_gqs_population_arguments <- function(fit) {
  pop <- fit$population %||% list()
  if (!isTRUE(pop$active)) return(list())
  persons <- pop$person_table_replay %||% pop$person_table
  if (!is.data.frame(persons) || is.null(pop$formula) ||
      !pop$person_id %in% names(persons)) {
    stop("The fitted population formula and person data are required for a same-model refit.", call. = FALSE)
  }
  for (name in names(pop$xlevels)) {
    persons[[name]] <- factor(persons[[name]], levels = pop$xlevels[[name]],
      ordered = is.ordered(persons[[name]]))
    if (!is.null(pop$contrasts[[name]])) stats::contrasts(persons[[name]]) <- pop$contrasts[[name]]
  }
  arguments <- list(population_formula = pop$formula, person_data = persons, person_id = pop$person_id)
  fixed_sd <- pop$fixed_sd %||% fit$config$population_spec$fixed_sd
  if (!is.null(fixed_sd)) arguments$population_sd <- fixed_sd
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
  model <- toupper(as.character(fit$config$model %||% "")[1L])
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

  if (mfrm_has_product_slopes(fit)) {
    grid <- cells[rep(seq_len(nrow(cells)), each = length(theta)), , drop = FALSE]
    grid$Theta <- rep(theta, nrow(cells))
    response <- mfrm_gpcm_response_evaluator(fit$prep, fit$config, grid)
    probabilities <- response$evaluate(fit$opt$par)$probabilities
    return(list(cells = cells, theta = theta,
      categories = seq_len(fit$config$n_cat) - 1L,
      values = as.numeric(t(probabilities))))
  }
  structure <- information_build_step_structure(fit, model)

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

  interactions <- as.data.frame(
    fit$interactions$effects %||% data.frame(), stringsAsFactors = FALSE
  )
  if (nrow(interactions) > 0L) {
    required <- c(
      "Interaction", "FacetA", "FacetA_Level", "FacetB", "FacetB_Level",
      "Estimate"
    )
    mfrmr_gqs_assert(
      all(required %in% names(interactions)),
      "The fitted interaction table is incomplete."
    )
    for (interaction in unique(as.character(interactions$Interaction))) {
      table <- interactions[interactions$Interaction == interaction,
                            , drop = FALSE]
      facet_a <- unique(as.character(table$FacetA))
      facet_b <- unique(as.character(table$FacetB))
      mfrmr_gqs_assert(
        length(facet_a) == 1L && length(facet_b) == 1L &&
          all(c(facet_a, facet_b) %in% names(cells)),
        "The fitted interaction could not be aligned to the design cells."
      )
      lookup <- stats::setNames(
        as.numeric(table$Estimate),
        paste(table$FacetA_Level, table$FacetB_Level, sep = "\r")
      )
      values <- as.numeric(lookup[paste(
        cells[[facet_a]], cells[[facet_b]], sep = "\r"
      )])
      mfrmr_gqs_assert(
        all(is.finite(values)),
        "The fitted interaction did not cover every observed design cell."
      )
      offsets <- offsets + values
    }
  }

  step_index <- slope_index <- NULL
  if (!identical(structure$kind, "common")) {
    step_index <- match(
      as.character(cells[[structure$step_facet]]), structure$step_levels
    )
    mfrmr_gqs_assert(
      all(is.finite(step_index)),
      "Step levels were not aligned with the observed design cells."
    )
  }
  if (identical(structure$kind, "step_and_slope_specific")) {
    slope_index <- match(
      as.character(cells[[structure$slope_facet]]), structure$slope_levels
    )
    mfrmr_gqs_assert(
      all(is.finite(slope_index)),
      "Slope levels were not aligned with the observed design cells."
    )
  }

  probability_blocks <- lapply(seq_len(nrow(cells)), function(cell_index) {
    eta <- theta + offsets[cell_index]
    probabilities <- if (identical(structure$kind, "common")) {
      structure$compute(eta)
    } else if (identical(structure$kind, "step_facet_specific")) {
      structure$compute(eta, rep(step_index[cell_index], length(theta)))
    } else {
      structure$compute(
        eta,
        rep(step_index[cell_index], length(theta)),
        rep(slope_index[cell_index], length(theta))
      )
    }
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

mfrmr_gqs_slope_keys <- function(fit) {
  slopes <- fit$slopes
  if (is.null(slopes) || nrow(slopes) == 0L) return(character())
  if (!mfrm_has_product_slopes(fit)) return(as.character(slopes$SlopeFacet))
  metadata <- mfrm_gpcm_product_metadata(fit$config$gpcm_spec)
  index <- mfrm_match_slope_table(metadata$Owner, metadata$Level, slopes)
  mfrmr_gqs_assert(length(index) == nrow(slopes) && !anyNA(index),
    "Slope owners and levels must match the fitted component specification.")
  keys <- character(nrow(slopes))
  keys[index] <- as.character(seq_along(index))
  keys
}

mfrmr_gqs_measurement_parameters <- function(fit) {
  facets <- as.data.frame(fit$facets$others %||% data.frame(),
                          stringsAsFactors = FALSE)
  steps <- as.data.frame(fit$steps %||% data.frame(),
                         stringsAsFactors = FALSE)
  interactions <- as.data.frame(
    fit$interactions$effects %||% data.frame(), stringsAsFactors = FALSE
  )
  slopes <- as.data.frame(fit$slopes %||% data.frame(),
                          stringsAsFactors = FALSE)
  out <- numeric()
  if (nrow(facets) > 0L) {
    value <- as.numeric(facets$Estimate)
    names(value) <- paste("facet", facets$Facet, facets$Level, sep = "::")
    out <- c(out, value)
  }
  if (nrow(steps) > 0L) {
    owner <- if ("StepFacet" %in% names(steps)) {
      as.character(steps$StepFacet)
    } else {
      rep("Common", nrow(steps))
    }
    value <- as.numeric(steps$Estimate)
    names(value) <- paste("step", owner, steps$Step, sep = "::")
    out <- c(out, value)
  }
  if (nrow(interactions) > 0L) {
    value <- as.numeric(interactions$Estimate)
    names(value) <- paste(
      "interaction", interactions$Interaction,
      interactions$FacetA_Level, interactions$FacetB_Level, sep = "::"
    )
    out <- c(out, value)
  }
  if (nrow(slopes) > 0L) {
    value <- as.numeric(slopes$OptimizerEstimate %||% slopes$Estimate)
    # Owner and level identities remain explicit in the output slope table.
    # Numbered coordinates avoid collisions when both families use the same IDs.
    names(value) <- paste("slope", mfrmr_gqs_slope_keys(fit), sep = "::")
    out <- c(out, value)
  }
  mfrmr_gqs_assert(
    length(out) > 0L && all(is.finite(out)) && !anyDuplicated(names(out)),
    "Fitted measurement coordinates were incomplete or duplicated."
  )
  out
}

mfrmr_gqs_person_scores <- function(fit, data, nodes) {
  if (nodes < 2L || mfrm_has_product_slopes(fit)) {
    return(data.frame(
      Person = character(), EAP = numeric(), PosteriorSD = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  estimates <- predict_mfrm_units(
    fit,
    data,
    scoring_quad_points = as.integer(nodes),
    readiness_policy = "review",
    person_data = fit$population$person_table_replay %||% fit$population$person_table,
    person_id = fit$population$person_id
  )$estimates
  estimates <- as.data.frame(estimates, stringsAsFactors = FALSE)
  mfrmr_gqs_assert(
    all(c("Person", "Estimate", "SD") %in% names(estimates)) &&
      nrow(estimates) > 0L && !anyDuplicated(as.character(estimates$Person)) &&
      all(is.finite(estimates$Estimate)) && all(is.finite(estimates$SD)),
    "Same-Person posterior summaries were incomplete."
  )
  data.frame(
    Person = as.character(estimates$Person),
    EAP = as.numeric(estimates$Estimate),
    PosteriorSD = as.numeric(estimates$SD),
    stringsAsFactors = FALSE
  )
}

mfrmr_gqs_score_difference <- function(reference, candidate, column) {
  if (nrow(reference) == 0L || nrow(candidate) == 0L) return(NA_real_)
  index <- match(reference$Person, candidate$Person)
  mfrmr_gqs_assert(
    !anyNA(index) && nrow(reference) == nrow(candidate),
    "Person identities changed across quadrature refits."
  )
  max(abs(candidate[[column]][index] - reference[[column]]))
}

mfrmr_gqs_raw_information <- function(fit) {
  product <- mfrm_has_product_slopes(fit)
  inference <- if (product) mfrm_gpcm_product_inference(fit) else NULL
  # Reuse only within this call. An early interval refusal may have no matrix;
  # raw diagnostic information remains separately available in that case.
  covariance <- inference$information %||% compute_mml_parameter_covariance(fit)
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
  slope_slice <- covariance$param_slices[["log_slopes"]] %||% integer()
  intervals <- if (product) mfrm_gpcm_intervals_from_inference(fit, inference) else NULL
  if (product && identical(covariance$status, "ok") && is.matrix(covariance$cov)) {
    target <- mfrm_gpcm_product_slope_targets(fit$opt$par, fit$config)$components
    index <- mfrm_match_slope_table(target$metadata$Owner, target$metadata$Level, slopes)
    log_covariance <- target$jacobian %*% covariance$cov %*% t(target$jacobian)
    slope_se[index] <- target$estimate * sqrt(pmax(diag(log_covariance), 0))
  } else if (!product && covariance$status %in% c("ok", "regularized") &&
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
  sigma_slice <- covariance$param_slices[["log_sigma2"]] %||% integer()
  if (covariance$status %in% c("ok", "regularized") &&
      is.matrix(covariance$cov) && length(sigma_slice) == 1L &&
      is.finite(population_sd)) {
    population_sd_se <- 0.5 * population_sd *
      sqrt(pmax(covariance$cov[sigma_slice, sigma_slice], 0))
  }

  list(
    intervals = intervals,
    slope_intervals = if (product) attr(intervals, "diagnostics") else if (identical(fit$config$model, "GPCM"))
      compute_mml_structural_parameter_se(fit, covariance = covariance)$slopes else NULL,
    covariance_status = as.character(covariance$status),
    covariance_detail = covariance$detail,
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
  if (length(reference) == 0L) return(0)
  candidate <- candidate[names(reference)]
  max(abs(reference - candidate))
}

mfrmr_gqs_extract_run <- function(fit, data, nodes,
                                  theta_range, theta_points) {
  information <- mfrmr_gqs_raw_information(fit)
  readiness <- as.data.frame(fit$readiness$fit %||% data.frame(),
                             stringsAsFactors = FALSE)
  mfrmr_gqs_assert(
    nrow(readiness) == 1L,
    "A quadrature fit did not retain exactly one fit-readiness row."
  )
  slopes <- as.data.frame(fit$slopes %||% data.frame(),
                          stringsAsFactors = FALSE)
  slope_labels <- if (nrow(slopes) > 0L) {
    as.character(slopes$SlopeFacet)
  } else {
    character()
  }
  product <- mfrm_has_product_slopes(fit)
  slope_keys <- mfrmr_gqs_slope_keys(fit)
  mfrmr_gqs_assert(
    length(slope_labels) == length(information$slope_values) &&
      !anyDuplicated(slope_keys),
    "Slope labels were not aligned with the observed-information coordinates."
  )

  slope_table <- if (length(slope_labels) == 0L) {
    data.frame(
      Nodes = integer(), SlopeFacet = character(),
      OptimizerEstimate = numeric(), RawObservedInformationSE = numeric(),
      SEStatus = character(), PublicSEEligible = logical(),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Nodes = rep(as.integer(nodes), length(slope_labels)),
      SlopeFacet = slope_labels,
      OptimizerEstimate = information$slope_values,
      RawObservedInformationSE = information$slope_se,
      SEStatus = rep(information$covariance_status, length(slope_labels)),
      PublicSEEligible = as.logical(
        slopes$SEEligible %||% rep(FALSE, nrow(slopes))
      ),
      stringsAsFactors = FALSE
    )
  }

  if (product) {
    slope_table$SlopeOwner <- slopes$SlopeOwner
    slope_table$ScaleReference <- slopes$ScaleReference
  }

  ci <- information$slope_intervals
  if (!is.null(ci) && nrow(ci)) {
    if (product) ci$SlopeFacet <- ci$SlopeLevel
    aligned <- if (product) mfrm_match_slope_table(slopes$SlopeOwner, slope_labels,
      ci) else match(slope_labels, ci$SlopeFacet)
    slope_table$PublicSEEligible <- ci$SEEligible[aligned]
    slope_table$CIEligible <- ci$CIEligible[aligned]
    slope_table$CI_Lower <- ci$CI_Lower[aligned]
    slope_table$CI_Upper <- ci$CI_Upper[aligned]
    slope_table$InferenceReview <- ci$InferenceReview[aligned]
  }
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
      EstimationConverged = identical(
        as.integer(fit$opt$convergence), 0L
      ),
      CovarianceStatus = information$covariance_status,
      CovarianceDetail = information$covariance_detail,
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
    slopes = slope_table,
    intervals = information$intervals,
    measurement_parameters = mfrmr_gqs_measurement_parameters(fit),
    slope_parameters = stats::setNames(
      information$slope_values,
      slope_keys
    ),
    slope_se = stats::setNames(information$slope_se, slope_keys),
    population_sd = information$population_sd,
    population_sd_se = information$population_sd_se,
    probabilities = mfrmr_gqs_design_probabilities(
      fit, theta_range, theta_points
    ),
    scores = mfrmr_gqs_person_scores(fit, data, nodes)
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

#' Review MML sensitivity to the quadrature grid
#'
#' Refit one ordered-response MML model to the same supplied response data at
#' two or more Gauss--Hermite quadrature counts. The original fit is reused at
#' its own quadrature count; all other fits reuse its stored model,
#' identification, anchor, interaction, optimizer, and population settings.
#'
#' @param fit An RSM, PCM, or GPCM MML `mfrm_fit` returned by
#'   [fit_mfrm()]. [gpcm_mml_quadrature_sensitivity()] accepts only GPCM fits.
#' @param data The original response data.frame used to create `fit`. Prepared
#'   response rows are compared semantically after refitting; row order may
#'   differ, but changed observations fail closed.
#' @param quad_points At least two distinct positive integers, including the
#'   quadrature count stored on `fit`. The default `c(31, 41)` is a comparison
#'   starting point, not a claim that either grid is adequate for the data.
#'   Refits retain the source fit's fixed or adaptive integration mode.
#' @param theta_range Two finite values defining the common ability grid used
#'   for fitted category-probability comparison.
#' @param theta_points Number of common-grid ability points; at least 21.
#' @param adaptive_quad_points Optional vector of at least two distinct integer
#'   orders >= 3, for example `c(31, 61)`. Adds a per-Person
#'   `quadrature_review` at each fit's fixed parameters, comparing a fixed-prior
#'   grid with mode/curvature-adapted grids. This separates integration error
#'   from parameter changes during refitting. Inspect changes between adaptive
#'   orders too; neither grid is certified exact.
#'   Also available for two slope families as a numerical diagnostic; it does
#'   not add Person scoring or change the fitted integration method.
#'
#' @details
#' This is an explicit refit diagnostic: neither [summary.mfrm_fit()] nor
#' [diagnose_mfrm()] invokes it automatically. It reports continuous changes
#' rather than classifying a fit as quadrature-stable or unstable. In
#' particular, it does not select a cutoff for acceptable changes, make
#' unavailable standard errors valid, or resolve other warnings about the fit.
#'
#' The probability comparison evaluates every observed combination of
#' non-Person facet levels on the same theta grid. It includes fitted two-way
#' facet interactions and, for GPCM, the complete-predictor slope action.
#' Same-Person EAP and posterior-SD changes use each fit's corresponding
#' quadrature count. A one-point grid has no public scoring route, so those two
#' changes are `NA` when it is the reference.
#'
#' Explicit intercept-only and covariate population models reuse the retained
#' person table, factor coding and formula, with design equality checked by Person.
#' Older fits without that data cannot be replayed.
#'
#' Raw slope and population-SD standard errors are computed from each local
#' observed-information Hessian for diagnostic comparison only. The public
#' parameter-level `SEEligible` state remains unchanged.
#'
#' @section Two slope families:
#' For a two-family GPCM, refits preserve the ordered slope owners,
#' fixed N(0,1) population, score ladder, engine, integration method and stopping
#' controls (`em_score_tol` for fixed-grid EM or `reltol` for adaptive direct MML).
#' Each refit reconstructs the source initialization policy from the observed
#' data. Fixed-grid EM uses the neutral vector; adaptive fits retain either
#' neutral-only or neutral-plus-EM initialization. Older adaptive fits without
#' this setting retain the earlier neutral-only procedure. The reference
#' estimate is not used as a new starting vector. Component slopes retain
#' `SlopeOwner` and `ScaleReference`: the
#' first family has geometric mean one; the second is free on the fixed ability
#' scale. Category-probability comparisons use the product of both slopes.
#'
#' For both two-family engines, `intervals` contains the separately checked experimental `confint()` result
#' for each grid, including failed checks, cautions and missing bounds. For
#' example, after comparing `quad_points = c(31, 61)`, use
#' `plot(out$intervals$q61)` and
#' `apa_table(out$intervals$q61, which = "numerical_checks")`. These checks also
#' evaluate q versus 2q-1 at that fit's saved parameters without another refit.
#' Passing these checks does not establish sampling coverage or resolve
#' other warnings about the fitted model. Adaptive fits use their moving-node
#' objective for both the information and the higher-order comparison.
#' Person-score comparisons in `summary` remain unavailable for two families:
#' EAP and posterior-SD changes there are `NA`, not zero. With
#' `adaptive_quad_points`, the separate `quadrature_review` table evaluates
#' posterior moments and log marginal likelihoods at each unchanged calibration.
#' These diagnostic integrals do not enable [predict_mfrm_units()] or supply
#' Person intervals.
#' The effective slope is the product of both facet slopes; moving the
#' integration grid retains the original N(0,1) density and its Jacobian.
#' This does not provide adaptive EM, update estimates or repair their intervals.
#'
#' More quadrature points improve the numerical approximation, not the amount
#' of observed information. In a sparse design, a small group of candidates
#' rated many times may need a finer grid than candidates with few ratings.
#' Inspect interval availability and measured changes; no fixed point count
#' is guaranteed to be adequate for every design. See the GPCM scope vignette
#' for a same-data example and its sampling limitations.
#'
#' @return An object of class `mfrm_quadrature_sensitivity` containing:
#' - `summary`: one comparison row per quadrature grid relative to the original
#'   fit, including likelihood, measurement-coordinate, probability, EAP, and
#'   posterior-SD changes;
#' - `quadrature_review`: when requested, unrounded per-Person fixed/adaptive
#'   log-marginal, EAP and posterior-SD differences, adaptive-order changes,
#'   mode/curvature and computation status. Unavailable rows retain their reason.
#'   `summary()` also supplies a compact `quadrature_overview`;
#' - `runs`: likelihood, gradient, curvature, population-scale, and readiness
#'   details for each fit;
#' - `slopes`: component-slope estimates, raw diagnostic SEs and freshly checked
#'   95% model intervals. The summary reports endpoint changes among jointly
#'   eligible levels and counts changes in interval availability. One-family
#'   slopes have geometric mean one; two-family scale references differ;
#' - `intervals`: for two families, saved `mfrm_slope_intervals` objects named
#'   by quadrature count, with their numerical checks and cautions;
#' - `conditions`: warnings and messages from refitting and result extraction;
#' - `fits`: the reference and refitted `mfrm_fit` objects;
#' - `settings` and `notes`: the fixed comparison contract and interpretation
#'   boundary.
#'
#' The object supports `print()`, `summary()`, `as.data.frame()`, and the
#' existing `apa_table()` list route (for example `apa_table(out)`).
#'
#' @seealso [fit_mfrm()], [diagnose_mfrm()], [apa_table()]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("mml_quadrature_sensitivity", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(
#'   toy, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", model = "RSM",
#'   quad_points = 31
#' )
#' sensitivity <- mml_quadrature_sensitivity(
#'   fit, toy, quad_points = c(31, 41), adaptive_quad_points = c(31, 61)
#' )
#' summary(sensitivity)
#' # Preserve small numerical differences when preparing the sensitivity table.
#' apa_table(sensitivity, digits = 5)
#' }
#' @export
mml_quadrature_sensitivity <- function(
    fit,
    data,
    quad_points = c(31L, 41L),
    theta_range = c(-4, 4),
    theta_points = 161L,
    adaptive_quad_points = NULL) {
  mfrmr_gqs_run(
    fit, data, quad_points, theta_range, theta_points,
    allowed_models = c("RSM", "PCM", "GPCM"),
    adaptive_quad_points = adaptive_quad_points
  )
}

mfrmr_gqs_run <- function(
    fit, data, quad_points, theta_range, theta_points, allowed_models,
    adaptive_quad_points = NULL) {
  adaptive_quad_points <- mfrmr_validate_adaptive_quad_points(adaptive_quad_points)
  product <- mfrm_has_product_slopes(fit)
  contract <- mfrmr_gqs_validate(
    fit, data, quad_points, theta_range, theta_points, allowed_models
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
      mfrmr_gqs_assert(!is.null(mfrm_population_design(fit)) &&
        identical(mfrm_population_design(fit), mfrm_population_design(candidate)),
        "The refit did not preserve the fitted population design by Person.")
      fits[[name]] <- candidate
    }
    details <- mfrmr_gqs_capture_conditions(mfrmr_gqs_extract_run(
      candidate,
      data,
      nodes,
      contract$theta_range,
      contract$theta_points
    ))
    extracted[[name]] <- details$value
    capture$warnings <- unique(c(capture$warnings, details$warnings))
    capture$messages <- unique(c(capture$messages, details$messages))
    condition_rows[[name]] <- mfrmr_gqs_condition_rows(nodes, capture)
  }

  reference <- extracted[[reference_name]]
  comparison_rows <- lapply(contract$quad_points, function(nodes) {
    candidate <- extracted[[paste0("q", nodes)]]
    slope_difference <- if (length(reference$slope_parameters) > 0L) {
      mfrmr_gqs_max_named_difference(
        reference$slope_parameters,
        candidate$slope_parameters,
        "slope"
      )
    } else {
      NA_real_
    }
    slope_se_difference <- if (
      length(reference$slope_se) > 0L &&
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
    interval_index <- if (product) mfrm_match_slope_table(reference$slopes$SlopeOwner,
      reference$slopes$SlopeFacet, candidate$slopes) else
      match(reference$slopes$SlopeFacet, candidate$slopes$SlopeFacet)
    nll_change_per_person <- candidate_run$NLLPerPerson -
      reference_run$NLLPerPerson
    data.frame(
      Model = contract$model,
      ReferenceNodes = contract$reference_nodes,
      Nodes = as.integer(nodes),
      IsReference = nodes == contract$reference_nodes,
      NLLChangePerPerson = nll_change_per_person,
      NLLAbsChangePerPerson = abs(nll_change_per_person),
      MeasurementParameterMaxAbsChange = mfrmr_gqs_max_named_difference(
        reference$measurement_parameters,
        candidate$measurement_parameters,
        "measurement"
      ),
      SlopeMaxAbsChange = slope_difference,
      RawSlopeSEMaxAbsChange = slope_se_difference,
      SlopeIntervalMaxAbsChange = {
        old <- reference$slopes; new <- candidate$slopes
        if (!"CIEligible" %in% names(old)) NA_real_ else {
          new <- new[interval_index, , drop = FALSE]
          usable <- old$CIEligible & new$CIEligible
          if (any(usable)) max(abs(c(old$CI_Lower[usable] - new$CI_Lower[usable],
            old$CI_Upper[usable] - new$CI_Upper[usable]))) else NA_real_
        }
      },
      SlopeIntervalEligibilityChanged = {
        old <- reference$slopes; new <- candidate$slopes
        if (!"CIEligible" %in% names(old)) NA_integer_ else
          sum(old$CIEligible != new$CIEligible[interval_index])
      },
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
      EAPMaxAbsChange = mfrmr_gqs_score_difference(
        reference$scores, candidate$scores, "EAP"
      ),
      PosteriorSDMaxAbsChange = mfrmr_gqs_score_difference(
        reference$scores, candidate$scores, "PosteriorSD"
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
      model = contract$model,
      reference_nodes = contract$reference_nodes,
      quad_points = contract$quad_points,
      mml_integration = fit$config$estimation_control$mml_integration %||% "fixed",
      theta_range = contract$theta_range,
      theta_points = contract$theta_points,
      response_contract = "same_semantically_canonical_prepared_rows",
      probability_contract =
        "all_observed_facet_cells_and_fitted_interactions_on_common_theta_grid",
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
  if (product) {
    out$intervals <- lapply(extracted, `[[`, "intervals")
    out$settings$person_score_comparison <- "unavailable_two_family"
    out$notes <- c(out$notes,
      paste("The summary does not compare two-family Person scores; EAP and posterior-SD changes are NA.",
        "Optional quadrature_review moments diagnose integration at fixed calibration; Person scoring remains unavailable."),
      paste("Component intervals are experimental and keep separate family scale references.",
        "Inspect $intervals for numerical checks, cautions and unavailable bounds; coverage is not established."))
  }
  if (!is.null(adaptive_quad_points)) {
    out$quadrature_review <- do.call(rbind, lapply(fits, function(candidate) {
      config <- candidate$config
      idx <- build_indices(candidate$prep, config$step_facet, config$slope_facet,
                           config$interaction_specs, gpcm_spec = config$gpcm_spec)
      params <- expand_params(candidate$opt$par, build_param_sizes(config), config)
      nodes <- as.integer(config$estimation_control$quad_points)
      mfrmr_adaptive_quadrature_review(
        idx, config, params, gauss_hermite_normal(nodes),
        candidate$prep$levels$Person, adaptive_quad_points
      )
    }))
    rownames(out$quadrature_review) <- NULL
    out$settings$adaptive_quad_points <- adaptive_quad_points
    out$notes <- c(out$notes, mfrmr_adaptive_quadrature_note())
  }
  class(out) <- c("mfrm_quadrature_sensitivity", "list")
  out
}

#' @rdname mml_quadrature_sensitivity
#' @export
gpcm_mml_quadrature_sensitivity <- function(
    fit,
    data,
    quad_points = c(31L, 41L),
    theta_range = c(-4, 4),
    theta_points = 161L,
    adaptive_quad_points = NULL) {
  mfrmr_gqs_run(
    fit, data, quad_points, theta_range, theta_points,
    allowed_models = "GPCM", adaptive_quad_points = adaptive_quad_points
  )
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
    MaxMeasurementParameterAbsChange = max_or_na(
      nonreference$MeasurementParameterMaxAbsChange
    ),
    MaxSlopeAbsChange = max_or_na(nonreference$SlopeMaxAbsChange),
    MaxPopulationSDAbsChange = max_or_na(
      nonreference$PopulationSDAbsChange
    ),
    MaxProbabilityAbsChange = max_or_na(
      nonreference$ProbabilityMaxAbsChange
    ),
    MaxEAPAbsChange = max_or_na(nonreference$EAPMaxAbsChange),
    MaxPosteriorSDAbsChange = max_or_na(
      nonreference$PosteriorSDMaxAbsChange
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
  if (!is.null(object$intervals)) out$intervals <- object$intervals
  if (!is.null(object$quadrature_review)) {
    out$quadrature_review <- object$quadrature_review
    out$quadrature_overview <- mfrmr_adaptive_quadrature_overview(object$quadrature_review)
  }
  class(out) <- c("summary.mfrm_quadrature_sensitivity", "list")
  out
}

#' @export
print.mfrm_quadrature_sensitivity <- function(x, digits = 5L, ...) {
  cat(x$settings$model, "-MML quadrature sensitivity\n", sep = "")
  cat("Reference grid: q=", x$settings$reference_nodes, "\n", sep = "")
  display <- as.data.frame(x$summary, stringsAsFactors = FALSE)
  numeric <- vapply(display, is.numeric, logical(1))
  display[numeric] <- lapply(display[numeric], round, digits = digits)
  print(display, row.names = FALSE)
  if (!is.null(x$quadrature_review)) {
    cat("\nFixed-parameter integration review (adaptive minus fixed)\n")
    print(mfrmr_adaptive_quadrature_overview(x$quadrature_review), row.names = FALSE)
    cat("Inspect $quadrature_review for adaptive-order changes and unavailable rows.\n")
  }
  cat(
    "No automatic stability classification or readiness change is applied.\n"
  )
  if (!is.null(x$intervals)) cat(
    "The summary does not compare two-family Person scores.\n",
    "Any adaptive-review moments are numerical diagnostics; inspect $intervals for experimental interval checks.\n", sep = "")
  invisible(x)
}

#' @export
print.summary.mfrm_quadrature_sensitivity <- function(x, digits = 5L, ...) {
  model <- as.character(x$comparison$Model[1L] %||% "MML")
  cat(model, "-MML quadrature sensitivity summary\n", sep = "")
  overview <- x$overview[setdiff(
    names(x$overview), c("StabilityClassification", "ReadinessEffect")
  )]
  numeric <- vapply(overview, is.numeric, logical(1))
  overview[numeric] <- lapply(overview[numeric], round, digits = digits)
  print(overview, row.names = FALSE)
  cat("\nGrid comparisons\n")
  comparison <- x$comparison
  numeric <- vapply(comparison, is.numeric, logical(1))
  comparison[numeric] <- lapply(comparison[numeric], round, digits = digits)
  print(comparison, row.names = FALSE)
  if (!is.null(x$quadrature_overview)) {
    cat("\nFixed-parameter integration review (adaptive minus fixed)\n")
    print(x$quadrature_overview, row.names = FALSE)
    cat("Inspect $quadrature_review for adaptive-order changes and unavailable rows.\n")
  }
  cat("No automatic stability classification or readiness change is applied.\n")
  if (!is.null(x$intervals)) cat(
    "The summary does not compare two-family Person scores.\n",
    "Any adaptive-review moments are numerical diagnostics; inspect $intervals for experimental interval checks.\n", sep = "")
  invisible(x)
}

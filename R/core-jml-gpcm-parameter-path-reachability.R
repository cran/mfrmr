# JML GPCM parameter-path reachability for declared lexicographic limits
# ==============================================================================
#
# P2b accepts cumulative category-utility directions directly. An inverse
# column-space test would make reachability tolerance-dependent and could alter
# an exact category tie. This contract therefore takes the safer forward route:
# the caller declares directions in the retained constrained free-additive
# coordinates, the exact retained adjacent-utility operator maps them to
# adjacent directions, and cumulative utilities are built with category zero
# fixed at zero. The resulting utility path is reachable by construction in
# the implementation parameter space.
#
# No direction, scale, or asymptotic path is inferred from an optimizer trace.

mfrmr_jml_gpcm_parameter_path_contract_version <- function() {
  "mfrmr-jml-gpcm-parameter-path-reachability-0.2.3-v1"
}

mfrmr_jml_gpcm_parameter_path_empty_stage_table <- function() {
  data.frame(
    Stage = integer(0),
    FreeCoordinates = integer(0),
    CoordinateL1 = numeric(0),
    CoordinateMaxAbs = numeric(0),
    AdjacentDirectionMaxAbs = numeric(0),
    CategoryDirectionMaxAbs = numeric(0),
    NonzeroCategoryContrasts = integer(0),
    ReachabilityState = character(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_parameter_path_empty_loading_table <- function() {
  data.frame(
    Stage = integer(0),
    OptimizerIndex = integer(0),
    Block = character(0),
    Coordinate = character(0),
    Facet = character(0),
    Level = character(0),
    Loading = numeric(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_parameter_path_scalar_numeric <- function(value) {
  converted <- tryCatch(
    suppressWarnings(as.numeric(value)[1]),
    error = function(e) NA_real_
  )
  if (length(converted) == 0L) NA_real_ else converted
}

mfrmr_jml_gpcm_parameter_path_directions <- function(directions) {
  if (is.null(directions) || length(directions) == 0L) return(list())
  if (is.matrix(directions) || is.data.frame(directions)) {
    converted <- tryCatch(
      matrix(
        as.numeric(directions),
        nrow = nrow(directions), ncol = ncol(directions),
        dimnames = dimnames(directions)
      ),
      error = function(e) NULL
    )
    if (is.null(converted)) return(NULL)
    return(lapply(seq_len(nrow(converted)), function(stage) {
      as.numeric(converted[stage, ])
    }))
  }
  if (is.atomic(directions) && is.null(dim(directions))) {
    converted <- tryCatch(
      suppressWarnings(as.numeric(directions)),
      error = function(e) NULL
    )
    if (is.null(converted)) return(NULL)
    return(list(converted))
  }
  if (!is.list(directions)) return(NULL)
  out <- lapply(directions, function(direction) {
    tryCatch(
      suppressWarnings(as.numeric(direction)),
      error = function(e) NULL
    )
  })
  if (any(vapply(out, is.null, logical(1)))) return(NULL)
  out
}

mfrmr_jml_gpcm_additive_utility_path <- function(
    adjacent_design,
    coordinate_stage_directions = list(),
    n_observations,
    n_categories,
    parameter_map = NULL,
    max_stages = 2L,
    max_design_nonzeros = 5e6) {
  empty_stages <- mfrmr_jml_gpcm_parameter_path_empty_stage_table()
  empty_loadings <- mfrmr_jml_gpcm_parameter_path_empty_loading_table()
  n_observations_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(n_observations)
  n_categories_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(n_categories)
  max_stages_value <- mfrmr_jml_gpcm_parameter_path_scalar_numeric(max_stages)
  max_nonzeros_value <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(max_design_nonzeros)
  design <- tryCatch(
    methods::as(Matrix::Matrix(adjacent_design, sparse = TRUE), "dgCMatrix"),
    error = function(e) NULL
  )
  directions <- mfrmr_jml_gpcm_parameter_path_directions(
    coordinate_stage_directions
  )
  design_rows <- if (is.null(design)) 0L else nrow(design)
  free_coordinates <- if (is.null(design)) 0L else ncol(design)
  design_nonzeros <- if (is.null(design)) 0 else length(design@x)

  finish <- function(state, evaluated, certified, detail, reason_codes,
                     stage_table = empty_stages,
                     loading_table = empty_loadings,
                     coordinate_directions = list(),
                     adjacent_directions = list(),
                     category_directions = list(),
                     map = data.frame()) {
    list(
      contract_version =
        mfrmr_jml_gpcm_parameter_path_contract_version(),
      state = state,
      evaluated = isTRUE(evaluated),
      parameter_space_reachability_checked = isTRUE(evaluated),
      parameter_space_reachability_certified = isTRUE(certified),
      reachability_basis =
        "forward_map_from_retained_constrained_free_additive_coordinates",
      inverse_projection_used = FALSE,
      inverse_tolerance_reachability_claim = FALSE,
      exact_tie_contract_preserved_by_forward_construction = isTRUE(certified),
      retained_observations = as.integer(n_observations_value),
      categories = as.integer(n_categories_value),
      adjacent_design_rows = as.integer(design_rows),
      additive_free_coordinates = as.integer(free_coordinates),
      adjacent_design_nonzeros = as.double(design_nonzeros),
      stages = as.integer(length(coordinate_directions)),
      maximum_stages = suppressWarnings(as.integer(max_stages_value)),
      category_zero_gauge = "first_category_utility_fixed_to_zero",
      probability_equivalent_row_shifts_required = FALSE,
      stage_table = stage_table,
      loading_table = loading_table,
      parameter_map = map,
      coordinate_stage_directions = coordinate_directions,
      adjacent_stage_directions = adjacent_directions,
      category_stage_directions = category_directions,
      path_inferred_from_optimizer_trace = FALSE,
      scale_exponents_inferred = FALSE,
      path_search_performed = FALSE,
      arbitrary_utility_direction_reachability_classified = FALSE,
      readiness_effect = "none_diagnostic_only",
      external_comparison_eligible = FALSE,
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "Reachability is certified only for directions declared in the",
        "retained constrained free-additive coordinate basis and mapped",
        "forward on retained response rows. The contract does not solve an",
        "arbitrary inverse utility problem, infer directions or scales from",
        "the optimizer, search paths, establish remainder-stable ties, or",
        "classify a global boundary, MML, uncertainty, readiness, recovery,",
        "or external comparison."
      )
    )
  }

  controls_valid <- length(n_observations_value) == 1L &&
    length(n_categories_value) == 1L && length(max_stages_value) == 1L &&
    length(max_nonzeros_value) == 1L &&
    all(is.finite(c(
      n_observations_value, n_categories_value,
      max_stages_value, max_nonzeros_value
    ))) &&
    n_observations_value >= 1 &&
    n_observations_value == floor(n_observations_value) &&
    n_categories_value >= 2 &&
    n_categories_value == floor(n_categories_value) &&
    max_stages_value >= 0 && max_stages_value == floor(max_stages_value) &&
    max_nonzeros_value >= 0 && max_nonzeros_value == floor(max_nonzeros_value)
  if (!controls_valid) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "Dimensions and workload controls must be finite nonnegative integers.",
      "parameter_path_control_invalid"
    ))
  }
  n_observations_value <- as.integer(n_observations_value)
  n_categories_value <- as.integer(n_categories_value)
  max_stages_value <- as.integer(max_stages_value)

  design_valid <- !is.null(design) &&
    nrow(design) == n_observations_value * (n_categories_value - 1L) &&
    all(is.finite(design@x))
  if (!design_valid) {
    return(finish(
      "not_evaluated_design", FALSE, FALSE,
      "The adjacent-utility design was malformed or dimensionally inconsistent.",
      "parameter_path_adjacent_design_invalid"
    ))
  }
  if (design_nonzeros > max_nonzeros_value) {
    return(finish(
      "not_evaluated_size_limit", FALSE, FALSE,
      "The adjacent-utility operator exceeded the declared sparse workload.",
      "parameter_path_design_nonzero_limit"
    ))
  }
  if (is.null(directions)) {
    return(finish(
      "not_evaluated_coordinate_directions", FALSE, FALSE,
      "Coordinate directions could not be converted to finite numeric vectors.",
      "parameter_path_coordinate_directions_invalid"
    ))
  }
  if (length(directions) > max_stages_value) {
    return(finish(
      "not_evaluated_stage_limit", FALSE, FALSE,
      "The declared additive path exceeded the forward-map stage limit.",
      "parameter_path_stage_limit"
    ))
  }
  direction_valid <- all(vapply(directions, function(direction) {
    length(direction) == free_coordinates && all(is.finite(direction))
  }, logical(1)))
  if (!direction_valid) {
    return(finish(
      "not_evaluated_coordinate_directions", FALSE, FALSE,
      "Every stage must match the free-additive dimension and be finite.",
      "parameter_path_coordinate_directions_invalid"
    ))
  }

  if (is.null(parameter_map)) {
    map <- data.frame(
      OptimizerIndex = seq_len(free_coordinates),
      Block = rep("additive", free_coordinates),
      Coordinate = paste0("additive:", seq_len(free_coordinates)),
      Facet = rep("", free_coordinates),
      Level = rep("", free_coordinates),
      stringsAsFactors = FALSE
    )
  } else {
    map <- tryCatch(
      as.data.frame(parameter_map, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    required <- c("OptimizerIndex", "Block", "Coordinate", "Facet", "Level")
    optimizer_index <- tryCatch(
      suppressWarnings(as.numeric(map$OptimizerIndex)),
      error = function(e) NULL
    )
    map_valid <- !is.null(map) && nrow(map) == free_coordinates &&
      all(required %in% names(map)) &&
      !anyNA(map$OptimizerIndex) && !anyDuplicated(map$OptimizerIndex) &&
      !is.null(optimizer_index) && all(is.finite(optimizer_index)) &&
      all(optimizer_index == floor(optimizer_index))
    if (!map_valid) {
      return(finish(
        "not_evaluated_parameter_map", FALSE, FALSE,
        "The free-additive parameter map was incomplete or duplicated.",
        "parameter_path_parameter_map_invalid"
      ))
    }
    map <- map[, required, drop = FALSE]
  }

  if (length(directions) == 0L) {
    return(finish(
      "no_additive_stage_reachability_vacuous", TRUE, TRUE,
      "No additive direction was declared; reachability is vacuous.",
      "parameter_path_no_additive_stage",
      map = map
    ))
  }

  adjacent_paths <- vector("list", length(directions))
  category_paths <- vector("list", length(directions))
  stage_rows <- vector("list", length(directions))
  loading_rows <- vector("list", length(directions))
  for (stage in seq_along(directions)) {
    direction <- directions[[stage]]
    adjacent <- matrix(
      as.numeric(design %*% direction),
      nrow = n_observations_value,
      ncol = n_categories_value - 1L
    )
    category <- matrix(
      0, nrow = n_observations_value, ncol = n_categories_value
    )
    for (transition in seq_len(n_categories_value - 1L)) {
      category[, transition + 1L] <-
        category[, transition] + adjacent[, transition]
    }
    nonzero_contrasts <- sum(adjacent != 0)
    if (!all(is.finite(adjacent)) || !all(is.finite(category)) ||
        nonzero_contrasts == 0L) {
      return(finish(
        "not_evaluated_null_utility_stage", FALSE, FALSE,
        "A declared free-coordinate stage produced no category contrast.",
        "parameter_path_null_utility_stage",
        coordinate_directions = directions,
        map = map
      ))
    }
    adjacent_paths[[stage]] <- adjacent
    category_paths[[stage]] <- category
    stage_rows[[stage]] <- data.frame(
      Stage = as.integer(stage),
      FreeCoordinates = as.integer(free_coordinates),
      CoordinateL1 = sum(abs(direction)),
      CoordinateMaxAbs = max(abs(direction)),
      AdjacentDirectionMaxAbs = max(abs(adjacent)),
      CategoryDirectionMaxAbs = max(abs(category)),
      NonzeroCategoryContrasts = as.integer(nonzero_contrasts),
      ReachabilityState = "reachable_by_forward_construction",
      stringsAsFactors = FALSE
    )
    loading_rows[[stage]] <- data.frame(
      Stage = rep(as.integer(stage), free_coordinates),
      OptimizerIndex = as.integer(map$OptimizerIndex),
      Block = as.character(map$Block),
      Coordinate = as.character(map$Coordinate),
      Facet = as.character(map$Facet),
      Level = as.character(map$Level),
      Loading = as.numeric(direction),
      stringsAsFactors = FALSE
    )
  }

  finish(
    "reachable_additive_utility_path", TRUE, TRUE,
    paste(
      "Every declared additive stage was mapped forward from the retained",
      "constrained free-coordinate basis to cumulative category utilities."
    ),
    "parameter_path_reachable_by_forward_construction",
    stage_table = do.call(rbind, stage_rows),
    loading_table = do.call(rbind, loading_rows),
    coordinate_directions = directions,
    adjacent_directions = adjacent_paths,
    category_directions = category_paths,
    map = map
  )
}

mfrmr_jml_gpcm_fit_parameter_hierarchy_limit <- function(
    fit,
    additive_coordinate_directions = list(),
    additive_scale_exponents = numeric(0),
    slope_stage_rates = NULL,
    slope_scale_exponents = numeric(0),
    max_slope_stages = 2L,
    max_additive_stages = 2L,
    max_observations = 20000L,
    max_categories = 20L,
    max_design_nonzeros = 5e6) {
  finish <- function(state, evaluated, detail, reason_codes,
                     reachability = NULL, limit = NULL,
                     certified = FALSE) {
    list(
      contract_version =
        mfrmr_jml_gpcm_parameter_path_contract_version(),
      state = state,
      evaluated = isTRUE(evaluated),
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      retained_fit_base_reconstructed = !is.null(limit),
      declared_parameter_path_reachability_checked =
        isTRUE(certified) || isTRUE(reachability$evaluated),
      declared_parameter_path_reachability_certified = isTRUE(
        reachability$parameter_space_reachability_certified
      ),
      parameter_reachable_limit_contract_completed = isTRUE(certified),
      declared_path_limit_classified =
        isTRUE(limit$declared_path_limit_classified),
      log_likelihood_limit = as.numeric(
        limit$log_likelihood_limit %||% NA_real_
      ),
      additive_path = reachability,
      limit = limit,
      parameter_directions_declared_by_caller = TRUE,
      path_inferred_from_optimizer_trace = FALSE,
      scale_exponents_inferred = FALSE,
      path_search_performed = FALSE,
      monotone_tail_certified = FALSE,
      competitive_boundary_certified = FALSE,
      common_subsequence_limit_certified = FALSE,
      arbitrary_path_limit_classified = FALSE,
      global_boundary_classified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "This wrapper classifies only a caller-declared parameter path on the",
        "retained response design. It does not infer a path, prove that an",
        "optimizer sequence follows it, search for a competitive tail, handle",
        "unreported remainders, establish common subsequence limits, or",
        "classify a global boundary, MML, uncertainty, readiness, recovery,",
        "or external comparison."
      )
    )
  }

  if (!inherits(fit, "mfrm_fit") || is.null(fit$config) ||
      is.null(fit$prep) || is.null(fit$opt$par)) {
    return(finish(
      "not_evaluated_fit", FALSE,
      "A current mfrm_fit with retained preparation and optimizer state is required.",
      "parameter_path_fit_invalid"
    ))
  }
  config <- fit$config
  if (!identical(config$model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The parameter-path GPCM limit wrapper applies only to GPCM.",
      "parameter_path_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(config$method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional parameter-path wrapper is not reused for MML.",
      "parameter_path_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction uses the existing PCM geometry.",
      "no_free_log_slope_coordinate"
    ))
  }
  source_contract <- as.character(
    config$boundary_audit$gpcm_lexicographic_limit$contract_version %||% ""
  )[1]
  if (!identical(
    source_contract,
    mfrmr_jml_gpcm_lexicographic_limit_contract_version()
  )) {
    return(finish(
      "not_evaluated_source_contract", FALSE,
      "The fit did not retain the current P2b lexicographic-limit contract.",
      "parameter_path_lexicographic_source_mismatch"
    ))
  }

  built <- tryCatch({
    idx <- build_indices(
      fit$prep,
      step_facet = config$step_facet,
      slope_facet = config$slope_facet,
      interaction_specs = config$interaction_specs
    )
    sizes <- build_param_sizes(config)
    adjacent <- mfrmr_estimability_adjacent_design(
      fit$prep, idx, config, sizes,
      include_person = TRUE, include_population_beta = FALSE
    )
    kernel <- mfrmr_gpcm_response_kernel_design(
      fit$prep, idx, config, sizes, fit$opt$par
    )
    params <- expand_params(fit$opt$par, sizes, config)
    list(
      idx = idx, sizes = sizes, adjacent = adjacent,
      kernel = kernel, params = params
    )
  }, error = function(e) e)
  if (inherits(built, "error")) {
    return(finish(
      "not_evaluated_design", FALSE,
      paste0("The retained additive operator could not be reconstructed: ",
             conditionMessage(built)),
      "parameter_path_retained_design_invalid"
    ))
  }

  n_observations <- length(built$idx$score_k)
  n_categories <- as.integer(config$n_cat)
  reachability <- mfrmr_jml_gpcm_additive_utility_path(
    adjacent_design = built$adjacent$design,
    coordinate_stage_directions = additive_coordinate_directions,
    n_observations = n_observations,
    n_categories = n_categories,
    parameter_map = built$adjacent$map,
    max_stages = max_additive_stages,
    max_design_nonzeros = max_design_nonzeros
  )
  if (!isTRUE(reachability$parameter_space_reachability_certified)) {
    return(finish(
      "not_evaluated_parameter_path", isTRUE(reachability$evaluated),
      "The declared free-coordinate additive path was not certified.",
      reachability$reason_codes,
      reachability = reachability
    ))
  }

  adjacent_base <- matrix(
    as.numeric(built$kernel$eta_minus_step),
    nrow = n_observations,
    ncol = n_categories - 1L
  )
  base_utilities <- matrix(
    0, nrow = n_observations, ncol = n_categories
  )
  for (transition in seq_len(n_categories - 1L)) {
    base_utilities[, transition + 1L] <-
      base_utilities[, transition] + adjacent_base[, transition]
  }
  weight <- as.numeric(
    built$idx$weight %||% rep(1, n_observations)
  )
  limit <- mfrmr_jml_gpcm_declared_hierarchy_limit(
    base_utilities = base_utilities,
    score_k = built$idx$score_k,
    slope_index = built$idx$slope_idx,
    base_log_slopes = built$params$log_slopes,
    slope_stage_rates = slope_stage_rates,
    slope_scale_exponents = slope_scale_exponents,
    additive_stage_directions = reachability$category_stage_directions,
    additive_scale_exponents = additive_scale_exponents,
    weight = weight,
    slope_levels = config$gpcm_spec$levels,
    max_slope_stages = max_slope_stages,
    max_additive_stages = max_additive_stages,
    max_observations = max_observations,
    max_categories = max_categories
  )
  if (!isTRUE(limit$declared_path_limit_classified)) {
    return(finish(
      "not_evaluated_lexicographic_limit", isTRUE(limit$evaluated),
      "The reachable parameter path did not satisfy the P2b limit contract.",
      limit$reason_codes,
      reachability = reachability,
      limit = limit
    ))
  }
  limit$declared_utility_path_parameter_reachability_checked <- TRUE
  limit$declared_utility_path_parameter_reachability_certified <- TRUE
  limit$parameter_reachability_contract <-
    mfrmr_jml_gpcm_parameter_path_contract_version()
  finish(
    paste0("parameter_reachable_", limit$state), TRUE,
    paste(
      "The caller-declared free-coordinate path was mapped through the",
      "retained constrained additive operator and classified by P2b."
    ),
    paste(
      "parameter_path_reachable_by_forward_construction",
      limit$reason_codes,
      sep = ";"
    ),
    reachability = reachability,
    limit = limit,
    certified = TRUE
  )
}

mfrmr_jml_gpcm_parameter_hierarchy_loglik_at <- function(result, distance) {
  valid <- is.list(result) && identical(
    as.character(result$contract_version %||% "")[1],
    mfrmr_jml_gpcm_parameter_path_contract_version()
  ) && isTRUE(result$declared_parameter_path_reachability_certified) &&
    isTRUE(result$limit$declared_path_limit_classified)
  if (!valid) {
    return(list(
      valid = FALSE,
      distance = mfrmr_jml_gpcm_parameter_path_scalar_numeric(distance),
      log_likelihood = NA_real_,
      row_log_probability = numeric(0),
      reason_codes = "parameter_path_result_invalid"
    ))
  }
  path <- result$limit$path_definition
  do.call(
    mfrmr_jml_gpcm_declared_hierarchy_loglik_at,
    c(list(distance = distance), path)
  )
}

audit_mfrm_jml_gpcm_parameter_path_scope <- function(
    prep,
    idx,
    config,
    sizes,
    lexicographic_limit =
      config$boundary_audit$gpcm_lexicographic_limit %||% list(),
    max_design_nonzeros = 5e6) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  source_matches <- identical(
    as.character(lexicographic_limit$contract_version %||% "")[1],
    mfrmr_jml_gpcm_lexicographic_limit_contract_version()
  )
  empty_map <- data.frame()
  finish <- function(state, evaluated, detail, reason_codes,
                     available = FALSE, design = NULL, map = empty_map) {
    list(
      contract_version =
        mfrmr_jml_gpcm_parameter_path_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      state = state,
      evaluated = isTRUE(evaluated),
      declared_free_coordinate_path_wrapper_available = isTRUE(available),
      forward_parameter_reachability_operator_available = isTRUE(available),
      additive_design_rows = if (is.null(design)) 0L else nrow(design),
      additive_free_coordinates = if (is.null(design)) 0L else ncol(design),
      additive_design_nonzeros = if (is.null(design)) 0 else length(design@x),
      coordinate_map = map,
      declared_paths_evaluated = 0L,
      production_parameter_path_reachability_checked = FALSE,
      production_path_limit_classified = FALSE,
      source_lexicographic_limit_contract = as.character(
        lexicographic_limit$contract_version %||% "missing"
      )[1],
      source_lexicographic_limit_contract_matches = isTRUE(source_matches),
      inverse_projection_used = FALSE,
      path_inferred_from_optimizer_trace = FALSE,
      scale_exponents_inferred = FALSE,
      path_search_performed = FALSE,
      arbitrary_path_limit_classified = FALSE,
      common_subsequence_limit_certified = FALSE,
      monotone_tail_certified = FALSE,
      competitive_boundary_certified = FALSE,
      global_boundary_classified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "The fit-level object records a forward operator and wrapper scope but",
        "does not declare, infer, or evaluate any path. Only caller-declared",
        "free-coordinate directions can receive a retained-design reachability",
        "certificate. Remainders, path search, competitiveness, common",
        "subsequence limits, global status, MML, uncertainty, readiness,",
        "recovery, and external comparison remain open."
      )
    )
  }

  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The parameter-path reachability operator applies only to GPCM.",
      "parameter_path_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional free-coordinate operator is not reused for MML.",
      "parameter_path_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction uses the existing PCM geometry.",
      "no_free_log_slope_coordinate"
    ))
  }
  if (!isTRUE(source_matches)) {
    return(finish(
      "not_evaluated_source_contract", FALSE,
      "The current P2b lexicographic-limit contract was unavailable.",
      "parameter_path_lexicographic_source_mismatch"
    ))
  }
  max_nonzeros <-
    mfrmr_jml_gpcm_parameter_path_scalar_numeric(max_design_nonzeros)
  if (!is.finite(max_nonzeros) || max_nonzeros < 0 ||
      max_nonzeros != floor(max_nonzeros)) {
    return(finish(
      "not_evaluated_control", FALSE,
      "The sparse operator workload must be a finite nonnegative integer.",
      "parameter_path_control_invalid"
    ))
  }
  adjacent <- tryCatch(
    mfrmr_estimability_adjacent_design(
      prep, idx, config, sizes,
      include_person = TRUE, include_population_beta = FALSE
    ),
    error = function(e) e
  )
  if (inherits(adjacent, "error")) {
    return(finish(
      "not_evaluated_design", FALSE,
      paste0("The retained additive operator could not be constructed: ",
             conditionMessage(adjacent)),
      "parameter_path_retained_design_invalid"
    ))
  }
  design <- methods::as(adjacent$design, "dgCMatrix")
  if (length(design@x) > max_nonzeros) {
    return(finish(
      "not_evaluated_size_limit", FALSE,
      "The retained additive operator exceeded the sparse workload limit.",
      "parameter_path_design_nonzero_limit",
      design = design,
      map = adjacent$map
    ))
  }
  finish(
    "retained_additive_parameter_operator_available_no_path_declared", TRUE,
    paste(
      "The constrained retained additive operator and coordinate map are",
      "available for caller-declared forward path construction."
    ),
    paste(
      "parameter_path_forward_operator_available",
      "production_path_not_declared",
      "global_boundary_not_classified",
      sep = ";"
    ),
    available = TRUE,
    design = design,
    map = adjacent$map
  )
}

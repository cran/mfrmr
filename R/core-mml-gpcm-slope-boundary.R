# MML GPCM fixed-quadrature log-slope boundary-path audit
# ==============================================================================
#
# For Person p, the implemented fixed-quadrature marginal likelihood is
#
#   L_p(t) = sum_q w_pq exp{ell_pq(t)}.
#
# Along a constant sum-zero expanded log-slope ray, its log derivative is the
# posterior-weighted average of the conditional-node derivatives.  Therefore,
# if every positive-loading response maximizes its unscaled category utility
# at every retained quadrature node, and every negative-loading response
# minimizes it at every node, the fixed-quadrature marginal log likelihood is
# nondecreasing along that ray.  A strict utility span supplies a strict finite
# derivative.  This is a sufficient boundary certificate for the declared
# numerical quadrature objective only.  A negative result is not a proof of a
# finite maximum, and no continuous-integral or joint additive/slope result is
# inferred.

mfrmr_mml_gpcm_slope_empty_group_table <- function() {
  data.frame(
    SlopeFacet = character(0),
    Level = character(0),
    SlopeIndex = integer(0),
    LogSlopeEstimate = numeric(0),
    SlopeEstimate = numeric(0),
    EffectiveObservations = integer(0),
    EffectiveWeight = numeric(0),
    QuadratureNodes = integer(0),
    MaxCompatible = logical(0),
    MinCompatible = logical(0),
    MaxSupportMargin = numeric(0),
    MinSupportMargin = numeric(0),
    StrictMaxRows = integer(0),
    StrictMinRows = integer(0),
    stringsAsFactors = FALSE
  )
}

audit_mfrm_mml_gpcm_slope_boundary <- function(
    prep,
    idx,
    config,
    sizes,
    par,
    quad_points,
    max_observations = 50000L,
    max_utility_elements = 5e6,
    max_slope_levels = 100L,
    max_direction_pairs = 10000L,
    max_pair_node_elements = 2e7,
    support_tolerance = sqrt(.Machine$double.eps),
    likelihood_tolerance = 1e-8) {
  method <- as.character(config$method %||% NA_character_)
  model <- as.character(config$model %||% NA_character_)
  slope_facet <- as.character(config$slope_facet %||% NA_character_)
  empty_groups <- mfrmr_mml_gpcm_slope_empty_group_table()
  empty_targets <- mfrmr_jml_gpcm_slope_empty_target_table()
  empty_pairs <- mfrmr_jml_gpcm_slope_empty_pair_table()
  empty_loadings <- mfrmr_jml_gpcm_slope_empty_loading_table()

  finish <- function(state, complete, scope_complete, detail,
                     target_status = empty_targets,
                     group_support = empty_groups,
                     certificates = empty_pairs,
                     direction_loadings = empty_loadings,
                     dimensions = data.frame(),
                     retained_log_likelihood = NA_real_,
                     optimizer_log_likelihood = NA_real_,
                     likelihood_difference = NA_real_) {
    list(
      contract_version = mfrmr_boundary_contract_version(),
      method = method,
      model = model,
      scope = paste(
        "MML GPCM constant log-slope rays with retained additive",
        "coordinates and the declared finite quadrature rule fixed"
      ),
      state = state,
      complete = isTRUE(complete),
      scope_complete = isTRUE(scope_complete),
      structural_identification_complete = FALSE,
      fixed_quadrature_certificate = TRUE,
      continuous_integral_certificate = FALSE,
      readiness_effect = "none_instrumentation_only",
      target_status = target_status,
      group_support = group_support,
      certificates = certificates,
      direction_loadings = direction_loadings,
      dimensions = dimensions,
      retained_log_likelihood = as.numeric(retained_log_likelihood),
      optimizer_log_likelihood = as.numeric(optimizer_log_likelihood),
      likelihood_difference = as.numeric(likelihood_difference),
      detail = detail,
      limitations = paste(
        "The certificate is sufficient only for the implemented finite-node",
        "quadrature objective while all additive, step, interaction, and",
        "population coordinates remain fixed. It is not a continuous-normal",
        "integral proof. A none-certified result does not establish finite",
        "GPCM slopes, a global finite maximum, structural identification,",
        "standard errors, confidence intervals, or comparison eligibility."
      )
    )
  }

  if (!identical(model, "GPCM")) {
    return(finish(
      "not_applicable_model", TRUE, FALSE,
      "The marginal log-slope boundary audit applies only to GPCM."
    ))
  }
  if (!identical(method, "MML")) {
    return(finish(
      "not_applicable_estimator", TRUE, FALSE,
      "The audit is specific to the fixed-quadrature MML objective."
    ))
  }

  levels <- as.character(config$gpcm_spec$levels %||% character(0))
  n_levels <- length(levels)
  n_obs <- length(idx$score_k %||% integer(0))
  n_categories <- as.integer(config$n_cat %||% 0L)
  n_nodes <- suppressWarnings(as.integer(quad_points)[1])
  n_pairs <- as.double(n_levels) * as.double(max(n_levels - 1L, 0L))
  utility_elements <- as.double(n_obs) * as.double(n_categories) *
    as.double(n_nodes)
  pair_node_elements <- n_pairs * as.double(n_obs) * as.double(n_nodes)
  dimensions <- data.frame(
    Observations = as.integer(n_obs),
    EffectiveObservations = NA_integer_,
    Categories = as.integer(n_categories),
    QuadratureNodes = as.integer(n_nodes),
    SlopeLevels = as.integer(n_levels),
    FreeLogSlopes = max(as.integer(n_levels - 1L), 0L),
    OrderedDirectionPairs = as.double(n_pairs),
    UtilityElements = as.double(utility_elements),
    PairNodeElements = as.double(pair_node_elements),
    CertifiedPairs = NA_integer_,
    stringsAsFactors = FALSE
  )

  control_values <- list(
    max_observations = max_observations,
    max_utility_elements = max_utility_elements,
    max_slope_levels = max_slope_levels,
    max_direction_pairs = max_direction_pairs,
    max_pair_node_elements = max_pair_node_elements,
    support_tolerance = support_tolerance,
    likelihood_tolerance = likelihood_tolerance
  )
  valid_controls <- all(vapply(control_values, function(value) {
    length(value) == 1L && !is.na(value) && is.finite(value) && value >= 0
  }, logical(1))) && length(n_nodes) == 1L && !is.na(n_nodes) && n_nodes >= 1L
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      paste(
        "Quadrature points, marginal slope-path execution limits, and",
        "tolerances must be finite nonnegative scalars."
      ),
      dimensions = dimensions
    ))
  }

  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  cache <- make_param_cache(sizes, config, idx, is_mml = TRUE)
  params <- tryCatch({
    if (is.null(par) || length(par) != free_dimension || any(!is.finite(par))) {
      stop("The retained optimizer vector is unavailable or malformed.")
    }
    cache$ensure(par)
    cache$params()
  }, error = function(e) e)
  if (inherits(params, "error") || length(params$slopes) != n_levels ||
      length(params$log_slopes) != n_levels ||
      any(!is.finite(params$slopes)) || any(params$slopes <= 0)) {
    return(finish(
      "not_evaluated_mapping", FALSE, FALSE,
      if (inherits(params, "error")) {
        paste0("GPCM slope expansion failed: ", conditionMessage(params))
      } else {
        "Expanded GPCM slope coordinates did not match the declared levels."
      },
      dimensions = dimensions
    ))
  }
  targets_not_evaluated <- mfrmr_jml_gpcm_slope_targets(
    levels, slope_facet, params$log_slopes, params$slopes
  )

  if (n_levels <= 1L || as.integer(sizes$log_slopes %||% 0L) == 0L) {
    targets <- mfrmr_jml_gpcm_slope_targets(
      levels, slope_facet, params$log_slopes, params$slopes,
      positive = rep(FALSE, n_levels), negative = rep(FALSE, n_levels),
      evaluation_state = "evaluated",
      reason = "no_free_log_slope_coordinate"
    )
    if (nrow(targets) > 0L) targets$CandidateStatus <- "fixed_unit_slope"
    dimensions$EffectiveObservations <- as.integer(n_obs)
    dimensions$CertifiedPairs <- 0L
    return(finish(
      "no_free_log_slope_coordinates", TRUE, TRUE,
      "The GPCM reduction has no free relative log-slope coordinate.",
      target_status = targets, dimensions = dimensions
    ))
  }

  if (n_obs > as.double(max_observations) ||
      utility_elements > as.double(max_utility_elements) ||
      n_levels > as.double(max_slope_levels) ||
      n_pairs > as.double(max_direction_pairs) ||
      pair_node_elements > as.double(max_pair_node_elements)) {
    targets_not_evaluated$CandidateStatus <- "not_evaluated_size_limit"
    targets_not_evaluated$ReasonCodes <-
      "bounded_marginal_slope_path_size_limit"
    return(finish(
      "not_evaluated_size_limit", FALSE, FALSE,
      paste0(
        "The marginal slope-path audit exceeded its execution limit (",
        n_obs, " observations; ", n_categories, " categories; ", n_nodes,
        " quadrature nodes; ", n_levels, " slope levels; ",
        format(n_pairs, scientific = FALSE), " ordered direction pairs)."
      ),
      target_status = targets_not_evaluated, dimensions = dimensions
    ))
  }

  score_k <- suppressWarnings(as.integer(idx$score_k))
  slope_index <- suppressWarnings(as.integer(idx$slope_idx))
  weight <- as.numeric(idx$weight %||% rep(1, n_obs))
  person_index <- suppressWarnings(as.integer(idx$person))
  valid_observation_map <- length(score_k) == n_obs &&
    length(slope_index) == n_obs && length(weight) == n_obs &&
    length(person_index) == n_obs && !anyNA(score_k) &&
    !anyNA(slope_index) && !anyNA(weight) && !anyNA(person_index) &&
    all(score_k >= 0L & score_k < n_categories) &&
    all(slope_index >= 1L & slope_index <= n_levels) &&
    all(person_index >= 1L) && all(is.finite(weight)) && all(weight >= 0)
  if (!valid_observation_map || n_obs < 1L || n_categories < 2L) {
    return(finish(
      "not_evaluated_design", FALSE, FALSE,
      paste(
        "The retained score, slope, Person, weight, or category map was",
        "malformed."
      ),
      target_status = targets_not_evaluated, dimensions = dimensions
    ))
  }
  effective <- weight > 0
  dimensions$EffectiveObservations <- sum(effective)
  if (!any(effective)) {
    return(finish(
      "not_evaluated_no_positive_weight", FALSE, FALSE,
      "No retained observation had positive likelihood weight.",
      target_status = targets_not_evaluated, dimensions = dimensions
    ))
  }

  quad <- gauss_hermite_normal(n_nodes)
  logprob_bundle <- tryCatch(
    mfrm_mml_logprob_bundle(
      idx = idx,
      config = config,
      quad = quad,
      params = params,
      base_eta = cache$base_eta(),
      step_cum = cache$step_cum(),
      include_probs = FALSE,
      include_linear_part = TRUE
    ),
    error = function(e) e
  )
  if (inherits(logprob_bundle, "error") ||
      length(logprob_bundle$linear_part_list %||% list()) != n_nodes) {
    return(finish(
      "not_evaluated_design", FALSE, FALSE,
      if (inherits(logprob_bundle, "error")) {
        paste0("Marginal GPCM utility construction failed: ",
               conditionMessage(logprob_bundle))
      } else {
        "The marginal GPCM utility bundle did not match the quadrature rule."
      },
      target_status = targets_not_evaluated, dimensions = dimensions
    ))
  }

  person_bundle <- mfrm_mml_person_bundle(
    log_prob_mat = logprob_bundle$log_prob_mat,
    person_int = logprob_bundle$person_int,
    quad_basis = logprob_bundle$quad_basis
  )
  retained_log_likelihood <- sum(person_bundle$log_marginal)
  optimizer_log_likelihood <- tryCatch(
    -mfrm_loglik_mml(par, idx, config, sizes, quad),
    error = function(e) NA_real_
  )
  likelihood_difference <- retained_log_likelihood - optimizer_log_likelihood
  likelihood_scale <- max(
    1, abs(retained_log_likelihood), abs(optimizer_log_likelihood)
  )
  if (!is.finite(retained_log_likelihood) ||
      !is.finite(optimizer_log_likelihood) ||
      !is.finite(likelihood_difference) ||
      abs(likelihood_difference) > likelihood_tolerance * likelihood_scale) {
    return(finish(
      "not_evaluated_likelihood_mismatch", FALSE, FALSE,
      paste(
        "The independently reconstructed fixed-quadrature marginal",
        "likelihood did not match the optimizer objective within tolerance."
      ),
      target_status = targets_not_evaluated, dimensions = dimensions,
      retained_log_likelihood = retained_log_likelihood,
      optimizer_log_likelihood = optimizer_log_likelihood,
      likelihood_difference = likelihood_difference
    ))
  }

  max_margin <- min_margin <- span <- matrix(
    NA_real_, nrow = n_obs, ncol = n_nodes
  )
  max_ties <- matrix(NA_integer_, nrow = n_obs, ncol = n_nodes)
  observation_index <- cbind(seq_len(n_obs), score_k + 1L)
  for (node in seq_len(n_nodes)) {
    utilities <- as.matrix(logprob_bundle$linear_part_list[[node]])
    if (!identical(dim(utilities), c(n_obs, n_categories)) ||
        any(!is.finite(utilities))) {
      return(finish(
        "not_evaluated_numerical", FALSE, FALSE,
        "A retained quadrature-node GPCM utility matrix was non-finite.",
        target_status = targets_not_evaluated, dimensions = dimensions,
        retained_log_likelihood = retained_log_likelihood,
        optimizer_log_likelihood = optimizer_log_likelihood,
        likelihood_difference = likelihood_difference
      ))
    }
    row_maximum <- apply(utilities, 1L, max)
    row_minimum <- apply(utilities, 1L, min)
    observed <- utilities[observation_index]
    max_margin[, node] <- observed - row_maximum
    min_margin[, node] <- row_minimum - observed
    span[, node] <- row_maximum - row_minimum
    max_ties[, node] <- rowSums(
      utilities == matrix(
        row_maximum, nrow = n_obs, ncol = n_categories
      )
    )
  }
  max_compatible_node <- max_margin == 0
  min_compatible_node <- min_margin == 0

  group_rows <- lapply(seq_len(n_levels), function(group) {
    rows <- which(effective & slope_index == group)
    has_rows <- length(rows) > 0L
    max_compatible <- has_rows && all(max_compatible_node[rows, , drop = FALSE])
    min_compatible <- has_rows && all(min_compatible_node[rows, , drop = FALSE])
    data.frame(
      SlopeFacet = slope_facet,
      Level = levels[group],
      SlopeIndex = as.integer(group),
      LogSlopeEstimate = as.numeric(params$log_slopes[group]),
      SlopeEstimate = as.numeric(params$slopes[group]),
      EffectiveObservations = as.integer(length(rows)),
      EffectiveWeight = if (has_rows) sum(weight[rows]) else 0,
      QuadratureNodes = as.integer(n_nodes),
      MaxCompatible = max_compatible,
      MinCompatible = min_compatible,
      MaxSupportMargin = if (has_rows) {
        min(max_margin[rows, , drop = FALSE])
      } else {
        NA_real_
      },
      MinSupportMargin = if (has_rows) {
        min(min_margin[rows, , drop = FALSE])
      } else {
        NA_real_
      },
      StrictMaxRows = if (max_compatible) {
        as.integer(sum(apply(
          span[rows, , drop = FALSE] > support_tolerance, 1L, any
        )))
      } else {
        0L
      },
      StrictMinRows = if (min_compatible) {
        as.integer(sum(apply(
          span[rows, , drop = FALSE] > support_tolerance, 1L, any
        )))
      } else {
        0L
      },
      stringsAsFactors = FALSE
    )
  })
  group_support <- do.call(rbind, group_rows)

  pair_grid <- expand.grid(
    PositiveIndex = seq_len(n_levels),
    NegativeIndex = seq_len(n_levels),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pair_grid <- pair_grid[
    pair_grid$PositiveIndex != pair_grid$NegativeIndex, , drop = FALSE
  ]
  rownames(pair_grid) <- NULL
  pair_rows <- lapply(seq_len(nrow(pair_grid)), function(pair_index) {
    positive_index <- pair_grid$PositiveIndex[pair_index]
    negative_index <- pair_grid$NegativeIndex[pair_index]
    positive_group <- group_support[positive_index, , drop = FALSE]
    negative_group <- group_support[negative_index, , drop = FALSE]
    strict_rows <- positive_group$StrictMaxRows +
      negative_group$StrictMinRows
    compatible <- isTRUE(positive_group$MaxCompatible) &&
      isTRUE(negative_group$MinCompatible) && strict_rows > 0L
    boundary_log_likelihood <- NA_real_
    improvement <- NA_real_
    if (compatible) {
      boundary_log_prob <- logprob_bundle$log_prob_mat
      positive_rows <- which(effective & slope_index == positive_index)
      negative_rows <- which(effective & slope_index == negative_index)
      positive_limit <- -log(max_ties[positive_rows, , drop = FALSE])
      positive_limit <- sweep(
        positive_limit, 1L, weight[positive_rows], `*`
      )
      negative_limit <- matrix(
        -log(n_categories),
        nrow = length(negative_rows), ncol = n_nodes
      )
      negative_limit <- sweep(
        negative_limit, 1L, weight[negative_rows], `*`
      )
      boundary_log_prob[positive_rows, ] <- positive_limit
      boundary_log_prob[negative_rows, ] <- negative_limit
      boundary_person <- mfrm_mml_person_bundle(
        log_prob_mat = boundary_log_prob,
        person_int = logprob_bundle$person_int,
        quad_basis = logprob_bundle$quad_basis
      )
      boundary_log_likelihood <- sum(boundary_person$log_marginal)
      improvement <- boundary_log_likelihood - retained_log_likelihood
    }
    certified <- compatible && is.finite(boundary_log_likelihood) &&
      is.finite(improvement) &&
      improvement >= -likelihood_tolerance * likelihood_scale
    data.frame(
      PairId = sprintf("MSP%05d", pair_index),
      PositiveLevel = levels[positive_index],
      NegativeLevel = levels[negative_index],
      PositiveIndex = as.integer(positive_index),
      NegativeIndex = as.integer(negative_index),
      MaxSupportMargin = positive_group$MaxSupportMargin,
      MinSupportMargin = negative_group$MinSupportMargin,
      StrictRows = as.integer(strict_rows),
      DirectionSum = 0,
      CurrentLogLikelihood = retained_log_likelihood,
      BoundaryLogLikelihood = boundary_log_likelihood,
      BoundaryImprovement = improvement,
      Certified = isTRUE(certified),
      stringsAsFactors = FALSE
    )
  })
  certificates <- if (length(pair_rows) == 0L) {
    empty_pairs
  } else {
    do.call(rbind, pair_rows)
  }
  certified_rows <- which(certificates$Certified)
  dimensions$CertifiedPairs <- as.integer(length(certified_rows))

  positive <- seq_len(n_levels) %in%
    certificates$PositiveIndex[certified_rows]
  negative <- seq_len(n_levels) %in%
    certificates$NegativeIndex[certified_rows]
  targets <- mfrmr_jml_gpcm_slope_targets(
    levels, slope_facet, params$log_slopes, params$slopes,
    positive = positive, negative = negative,
    evaluation_state = "evaluated_fixed_quadrature_marginal",
    reason = "no_marginal_monotone_pair_certified_at_retained_additive_point"
  )
  targets$ReasonCodes[positive | negative] <-
    "certified_fixed_quadrature_marginal_slope_boundary_path"
  targets$ReasonCodes[group_support$EffectiveObservations == 0L] <-
    "no_positive_weight_observations"

  loading_rows <- list()
  loading_cursor <- 0L
  slices <- build_param_slices(sizes)
  optimizer_indices <- as.integer(slices$log_slopes %||% integer(0))
  for (certificate_index in certified_rows) {
    certificate <- certificates[certificate_index, , drop = FALSE]
    expanded_direction <- rep(0, n_levels)
    expanded_direction[certificate$PositiveIndex] <- 1
    expanded_direction[certificate$NegativeIndex] <- -1
    expanded_keep <- which(expanded_direction != 0)
    loading_cursor <- loading_cursor + 1L
    loading_rows[[loading_cursor]] <- data.frame(
      PairId = certificate$PairId,
      CoordinateType = "expanded_log_slope",
      OptimizerIndex = rep(NA_integer_, length(expanded_keep)),
      Coordinate = paste0("expanded_log_slope:", levels[expanded_keep]),
      Level = levels[expanded_keep],
      Loading = expanded_direction[expanded_keep],
      stringsAsFactors = FALSE
    )
    free_direction <- expanded_direction[seq_len(n_levels - 1L)]
    free_keep <- which(free_direction != 0)
    if (length(free_keep) > 0L) {
      loading_cursor <- loading_cursor + 1L
      loading_rows[[loading_cursor]] <- data.frame(
        PairId = certificate$PairId,
        CoordinateType = "optimizer_log_slope",
        OptimizerIndex = optimizer_indices[free_keep],
        Coordinate = paste0("log_slopes:", levels[free_keep]),
        Level = levels[free_keep],
        Loading = free_direction[free_keep],
        stringsAsFactors = FALSE
      )
    }
  }
  direction_loadings <- if (length(loading_rows) == 0L) {
    empty_loadings
  } else {
    do.call(rbind, loading_rows)
  }

  any_certified <- length(certified_rows) > 0L
  finish(
    state = if (any_certified) {
      "certified_fixed_quadrature_marginal_boundary_path"
    } else {
      "none_certified_fixed_quadrature_marginal"
    },
    complete = TRUE,
    scope_complete = TRUE,
    detail = if (any_certified) {
      paste0(
        "Certified ", length(certified_rows), " of ", nrow(certificates),
        " ordered slope-level pairs as nondecreasing constant sum-zero",
        " boundary paths for the declared ", n_nodes,
        "-node marginal objective."
      )
    } else {
      paste0(
        "No constant sum-zero marginal slope boundary path was certified",
        " among ", nrow(certificates), " ordered pairs for the declared ",
        n_nodes, "-node objective."
      )
    },
    target_status = targets,
    group_support = group_support,
    certificates = certificates,
    direction_loadings = direction_loadings,
    dimensions = dimensions,
    retained_log_likelihood = retained_log_likelihood,
    optimizer_log_likelihood = optimizer_log_likelihood,
    likelihood_difference = likelihood_difference
  )
}

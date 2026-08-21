# JML GPCM joint additive/log-slope boundary-path audit
# ==============================================================================
#
# This audit certifies a workload-bounded family of nonlinear paths. Expanded
# log slopes use constant sum-zero rates. Any nonzero rate vector with a
# strictly favorable leading negative tier has a canonical representative:
# positive-rate levels, zero-rate levels, a leading negative group at -1, and
# an unconstrained deeper negative group at -2. Positive rates are equal and
# chosen to make the expanded rates sum to zero. Existing ordered +1/-1 pairs
# are exact special cases. The free additive coordinates move linearly along a
# constrained direction d. If
#
# * every positive-rate observation uniquely favors its observed category in
#   the cumulative additive direction,
# * every zero-rate observation weakly favors its observed category, and
# * the leading-negative observations have a strictly favorable aggregate
#   term as their slopes tend to zero,
#
# then the likelihood derivative is positive on an asymptotic tail.  A path is
# retained as a competitive boundary candidate only when its analytic boundary
# likelihood is no worse than the retained likelihood.  These sufficient
# Deeper negative groups vanish faster and require no directional inequality.
# These conditions detect joint paths missed by additive-only, slope-only, and
# pair-only audits. Failure to certify this family is not evidence of a finite
# MLE, and a positive candidate is not a global identification proof. Curved
# or rate-nonconvergent paths remain outside the audit.

mfrmr_jml_gpcm_joint_empty_pair_table <- function() {
  data.frame(
    PairId = character(0),
    PositiveLevel = character(0),
    NegativeLevel = character(0),
    PositiveIndex = integer(0),
    NegativeIndex = integer(0),
    SolverStatus = integer(0),
    MarginCapacity = numeric(0),
    PositiveMinimumMargin = numeric(0),
    NeutralMinimumMargin = numeric(0),
    NegativeLeadingCoefficient = numeric(0),
    AdditiveDirectionL1 = numeric(0),
    CurrentLogLikelihood = numeric(0),
    BoundaryLogLikelihood = numeric(0),
    BoundaryImprovement = numeric(0),
    AnalyticTailCertified = logical(0),
    CompetitiveBoundary = logical(0),
    Certified = logical(0),
    EvaluationState = character(0),
    ReasonCodes = character(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_joint_empty_target_table <- function() {
  data.frame(
    ParameterId = character(0),
    ParameterClass = character(0),
    Facet = character(0),
    Level = character(0),
    PositiveBoundaryCandidate = logical(0),
    NegativeBoundaryCandidate = logical(0),
    CandidateStatus = character(0),
    EvaluationState = character(0),
    ReasonCodes = character(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_joint_empty_loading_table <- function() {
  data.frame(
    PairId = character(0),
    CoordinateType = character(0),
    OptimizerIndex = integer(0),
    Coordinate = character(0),
    Level = character(0),
    Loading = numeric(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_joint_empty_rate_table <- function() {
  data.frame(
    RateId = character(0),
    PositiveLevels = character(0),
    ZeroLevels = character(0),
    LeadingNegativeLevels = character(0),
    DeeperNegativeLevels = character(0),
    ExpandedRates = character(0),
    SolverStatus = integer(0),
    MarginCapacity = numeric(0),
    PositiveMinimumMargin = numeric(0),
    ZeroMinimumMargin = numeric(0),
    LeadingNegativeCoefficient = numeric(0),
    AdditiveDirectionL1 = numeric(0),
    CurrentLogLikelihood = numeric(0),
    BoundaryLogLikelihood = numeric(0),
    BoundaryImprovement = numeric(0),
    AnalyticTailCertified = logical(0),
    CompetitiveBoundary = logical(0),
    Certified = logical(0),
    EvaluationState = character(0),
    ReasonCodes = character(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_joint_rate_candidate_count <- function(n_levels) {
  n_levels <- as.integer(n_levels)
  if (!is.finite(n_levels) || n_levels < 2L) return(0)
  as.double(4^n_levels - 2 * 3^n_levels + 2^n_levels)
}

mfrmr_jml_gpcm_joint_rate_candidates <- function(n_levels) {
  n_levels <- as.integer(n_levels)
  if (!is.finite(n_levels) || n_levels < 2L) return(list())
  assignments <- expand.grid(
    rep(list(c(1L, 0L, -1L, -2L)), n_levels),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  assignments <- as.matrix(assignments)
  keep <- rowSums(assignments == 1L) > 0L &
    rowSums(assignments == -1L) > 0L
  assignments <- assignments[keep, , drop = FALSE]
  pair <- rowSums(assignments == 1L) == 1L &
    rowSums(assignments == -1L) == 1L &
    rowSums(assignments == -2L) == 0L
  assignments <- assignments[!pair, , drop = FALSE]
  lapply(seq_len(nrow(assignments)), function(index) {
    assignment <- as.integer(assignments[index, ])
    positive <- which(assignment == 1L)
    zero <- which(assignment == 0L)
    leading_negative <- which(assignment == -1L)
    deeper_negative <- which(assignment == -2L)
    rates <- numeric(n_levels)
    rates[positive] <-
      (length(leading_negative) + 2 * length(deeper_negative)) /
      length(positive)
    rates[leading_negative] <- -1
    rates[deeper_negative] <- -2
    list(
      positive = positive,
      zero = zero,
      leading_negative = leading_negative,
      deeper_negative = deeper_negative,
      rates = rates
    )
  })
}

mfrmr_jml_gpcm_joint_targets <- function(levels, slope_facet,
                                          positive = NULL,
                                          negative = NULL,
                                          evaluation_state = "evaluated",
                                          reason = "no_joint_pair_path_certified") {
  levels <- as.character(levels)
  if (length(levels) == 0L) {
    return(mfrmr_jml_gpcm_joint_empty_target_table())
  }
  positive <- as.logical(positive %||% rep(FALSE, length(levels)))
  negative <- as.logical(negative %||% rep(FALSE, length(levels)))
  status <- rep("none_certified_in_joint_pair_paths", length(levels))
  status[positive & !negative] <- "joint_boundary_path_high"
  status[!positive & negative] <- "joint_boundary_path_low"
  status[positive & negative] <- "joint_boundary_path_both"
  data.frame(
    ParameterId = paste0("Slope:", slope_facet, ":", levels),
    ParameterClass = rep("log_slope", length(levels)),
    Facet = rep(as.character(slope_facet), length(levels)),
    Level = levels,
    PositiveBoundaryCandidate = positive,
    NegativeBoundaryCandidate = negative,
    CandidateStatus = status,
    EvaluationState = rep(evaluation_state, length(levels)),
    ReasonCodes = ifelse(positive | negative,
                         "certified_joint_nonlinear_boundary_candidate",
                         reason),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_joint_rate_lp <- function(
    contrast,
    contrast_observation,
    slope_index,
    weight,
    n_categories,
    positive_indices,
    zero_indices,
    leading_negative_indices,
    max_lp_nonzeros,
    lp_timeout,
    certificate_tolerance) {
  contrast <- methods::as(contrast, "dgCMatrix")
  contrast_observation <- as.integer(contrast_observation)
  n_parameters <- ncol(contrast)
  effective_contrast <- weight[contrast_observation] > 0
  positive_rows <- which(
    effective_contrast &
      slope_index[contrast_observation] %in% positive_indices
  )
  negative_rows <- which(
    effective_contrast &
      slope_index[contrast_observation] %in% leading_negative_indices
  )
  neutral_rows <- which(
    effective_contrast & slope_index[contrast_observation] %in% zero_indices
  )
  empty_result <- list(
    evaluated = TRUE,
    certified = FALSE,
    solver_status = NA_integer_,
    margin_capacity = NA_real_,
    positive_minimum_margin = NA_real_,
    neutral_minimum_margin = NA_real_,
    negative_leading_coefficient = NA_real_,
    direction_l1 = NA_real_,
    direction = rep(NA_real_, n_parameters),
    lp_constraints = 0L,
    lp_variables = as.integer(2L * n_parameters + 1L),
    lp_nonzeros = 0,
    reason = "not_evaluated"
  )
  if (length(positive_rows) == 0L || length(negative_rows) == 0L) {
    empty_result$reason <-
      "joint_rate_lacks_positive_and_leading_negative_support"
    return(empty_result)
  }

  positive <- contrast[positive_rows, , drop = FALSE]
  neutral <- contrast[neutral_rows, , drop = FALSE]
  negative <- contrast[negative_rows, , drop = FALSE]
  negative_weight <- weight[contrast_observation[negative_rows]] /
    as.numeric(n_categories)
  negative_vector <- as.numeric(Matrix::colSums(
    Matrix::Diagonal(x = negative_weight) %*% negative
  ))

  n_positive <- nrow(positive)
  n_neutral <- nrow(neutral)
  negative_constraint_row <- n_positive + n_neutral + 1L
  bound_row_offset <- negative_constraint_row
  n_constraints <- n_positive + n_neutral + 1L + n_parameters
  n_variables <- 2L * n_parameters + 1L
  margin_index <- n_variables

  support <- rbind(positive, neutral)
  support_triplet <- Matrix::summary(support)
  triplet <- matrix(numeric(0), nrow = 0L, ncol = 3L)
  if (nrow(support_triplet) > 0L) {
    triplet <- rbind(
      triplet,
      cbind(support_triplet$i, support_triplet$j, support_triplet$x),
      cbind(support_triplet$i,
            n_parameters + support_triplet$j,
            -support_triplet$x)
    )
  }
  if (n_positive > 0L) {
    triplet <- rbind(
      triplet,
      cbind(seq_len(n_positive), margin_index, -1)
    )
  }
  negative_keep <- which(negative_vector != 0)
  if (length(negative_keep) > 0L) {
    triplet <- rbind(
      triplet,
      cbind(negative_constraint_row, negative_keep,
            -negative_vector[negative_keep]),
      cbind(negative_constraint_row,
            n_parameters + negative_keep,
            negative_vector[negative_keep])
    )
  }
  triplet <- rbind(
    triplet,
    c(negative_constraint_row, margin_index, -1),
    cbind(bound_row_offset + seq_len(n_parameters),
          seq_len(n_parameters), 1),
    cbind(bound_row_offset + seq_len(n_parameters),
          n_parameters + seq_len(n_parameters), 1)
  )
  storage.mode(triplet) <- "double"
  empty_result$lp_constraints <- as.integer(n_constraints)
  empty_result$lp_nonzeros <- as.double(nrow(triplet))
  if (nrow(triplet) > max_lp_nonzeros) {
    empty_result$evaluated <- FALSE
    empty_result$reason <- "joint_rate_lp_nonzero_limit"
    return(empty_result)
  }

  fit <- tryCatch(
    lpSolve::lp(
      direction = "max",
      objective.in = c(rep(0, 2L * n_parameters), 1),
      const.dir = c(
        rep(">=", n_positive + n_neutral + 1L),
        rep("<=", n_parameters)
      ),
      const.rhs = c(
        rep(0, n_positive + n_neutral + 1L),
        rep(1, n_parameters)
      ),
      dense.const = triplet,
      timeout = as.integer(lp_timeout)
    ),
    error = function(e) e
  )
  if (inherits(fit, "error") ||
      !identical(as.integer(fit$status %||% -1L), 0L)) {
    empty_result$evaluated <- FALSE
    empty_result$solver_status <- if (inherits(fit, "error")) {
      NA_integer_
    } else {
      as.integer(fit$status)
    }
    empty_result$reason <- "joint_rate_linear_program_failed"
    return(empty_result)
  }

  solution <- as.numeric(fit$solution)
  direction <- solution[seq_len(n_parameters)] -
    solution[n_parameters + seq_len(n_parameters)]
  positive_margin <- as.numeric(positive %*% direction)
  neutral_margin <- if (n_neutral > 0L) {
    as.numeric(neutral %*% direction)
  } else {
    numeric(0)
  }
  negative_leading <- -sum(negative_vector * direction)
  margin_capacity <- as.numeric(fit$objval)
  positive_minimum <- min(positive_margin)
  neutral_minimum <- if (length(neutral_margin) > 0L) {
    min(neutral_margin)
  } else {
    NA_real_
  }
  certified <- is.finite(margin_capacity) &&
    margin_capacity > certificate_tolerance &&
    is.finite(positive_minimum) &&
    positive_minimum > certificate_tolerance &&
    (length(neutral_margin) == 0L ||
       (is.finite(neutral_minimum) &&
        neutral_minimum >= -certificate_tolerance)) &&
    is.finite(negative_leading) &&
    negative_leading > certificate_tolerance

  list(
    evaluated = TRUE,
    certified = isTRUE(certified),
    solver_status = as.integer(fit$status),
    margin_capacity = margin_capacity,
    positive_minimum_margin = positive_minimum,
    neutral_minimum_margin = neutral_minimum,
    negative_leading_coefficient = negative_leading,
    direction_l1 = sum(abs(direction)),
    direction = direction,
    lp_constraints = as.integer(n_constraints),
    lp_variables = as.integer(n_variables),
    lp_nonzeros = as.double(nrow(triplet)),
    reason = if (isTRUE(certified)) {
      "certified_joint_constant_rate_asymptotic_tail"
    } else {
      "joint_rate_failed_postsolve_certificate"
    }
  )
}

mfrmr_jml_gpcm_joint_pair_lp <- function(
    contrast,
    contrast_observation,
    slope_index,
    weight,
    n_categories,
    positive_index,
    negative_index,
    max_lp_nonzeros,
    lp_timeout,
    certificate_tolerance) {
  out <- mfrmr_jml_gpcm_joint_rate_lp(
    contrast = contrast,
    contrast_observation = contrast_observation,
    slope_index = slope_index,
    weight = weight,
    n_categories = n_categories,
    positive_indices = positive_index,
    zero_indices = setdiff(
      seq_len(max(slope_index)), c(positive_index, negative_index)
    ),
    leading_negative_indices = negative_index,
    max_lp_nonzeros = max_lp_nonzeros,
    lp_timeout = lp_timeout,
    certificate_tolerance = certificate_tolerance
  )
  reason_map <- c(
    joint_rate_lacks_positive_and_leading_negative_support =
      "slope_pair_lacks_positive_weight_support",
    joint_rate_lp_nonzero_limit = "joint_pair_lp_nonzero_limit",
    joint_rate_linear_program_failed = "joint_pair_linear_program_failed",
    certified_joint_constant_rate_asymptotic_tail =
      "certified_joint_pair_asymptotic_tail",
    joint_rate_failed_postsolve_certificate =
      "joint_pair_failed_postsolve_certificate"
  )
  if (out$reason %in% names(reason_map)) {
    out$reason <- unname(reason_map[[out$reason]])
  }
  out
}

mfrmr_jml_gpcm_joint_rate_boundary_limit <- function(
    utilities,
    category_direction,
    score_k,
    slope_index,
    slopes,
    weight,
    positive_indices,
    zero_indices,
    negative_indices,
    tolerance) {
  n_obs <- nrow(utilities)
  n_categories <- ncol(utilities)
  row_limit <- numeric(n_obs)
  valid <- TRUE
  for (observation in seq_len(n_obs)) {
    if (weight[observation] <= 0) next
    group <- slope_index[observation]
    observed <- score_k[observation] + 1L
    if (group %in% positive_indices) {
      row_limit[observation] <- 0
    } else if (group %in% negative_indices) {
      row_limit[observation] <- -log(n_categories)
    } else if (group %in% zero_indices) {
      maximum <- max(category_direction[observation, ])
      top <- which(
        category_direction[observation, ] >= maximum - tolerance
      )
      if (!observed %in% top) {
        valid <- FALSE
        row_limit[observation] <- NA_real_
        next
      }
      scaled <- slopes[group] * utilities[observation, top]
      scaled_max <- max(scaled)
      denominator <- scaled_max + log(sum(exp(scaled - scaled_max)))
      row_limit[observation] <-
        slopes[group] * utilities[observation, observed] - denominator
    } else {
      valid <- FALSE
      row_limit[observation] <- NA_real_
    }
  }
  list(
    valid = isTRUE(valid) && all(is.finite(row_limit[weight > 0])),
    row_limit = row_limit,
    log_likelihood = if (isTRUE(valid)) {
      sum(weight * row_limit)
    } else {
      NA_real_
    }
  )
}

mfrmr_jml_gpcm_joint_boundary_limit <- function(
    utilities,
    category_direction,
    score_k,
    slope_index,
    slopes,
    weight,
    positive_index,
    negative_index,
    tolerance) {
  mfrmr_jml_gpcm_joint_rate_boundary_limit(
    utilities = utilities,
    category_direction = category_direction,
    score_k = score_k,
    slope_index = slope_index,
    slopes = slopes,
    weight = weight,
    positive_indices = positive_index,
    zero_indices = setdiff(
      seq_along(slopes), c(positive_index, negative_index)
    ),
    negative_indices = negative_index,
    tolerance = tolerance
  )
}

audit_mfrm_jml_gpcm_joint_boundary <- function(
    prep,
    idx,
    config,
    sizes,
    par,
    max_observations = 20000L,
    max_contrast_rows = 100000L,
    max_additive_coordinates = 250L,
    max_direction_pairs = 200L,
    max_lp_nonzeros = 2e6,
    lp_timeout = 2L,
    certificate_tolerance = 1e-7,
    likelihood_tolerance = 1e-8,
    max_rate_candidates = 1000L) {
  method <- as.character(config$method %||% NA_character_)
  model <- as.character(config$model %||% NA_character_)
  slope_facet <- as.character(config$slope_facet %||% NA_character_)
  empty_pairs <- mfrmr_jml_gpcm_joint_empty_pair_table()
  empty_targets <- mfrmr_jml_gpcm_joint_empty_target_table()
  empty_loadings <- mfrmr_jml_gpcm_joint_empty_loading_table()
  empty_rates <- mfrmr_jml_gpcm_joint_empty_rate_table()

  finish <- function(state, complete, scope_complete, detail,
                     target_status = empty_targets,
                     certificates = empty_pairs,
                     direction_loadings = empty_loadings,
                     rate_certificates = empty_rates,
                     rate_direction_loadings = empty_loadings,
                     rate_scope_state = "not_evaluated",
                     rate_scope_complete = FALSE,
                     positive_certificate_present = FALSE,
                     dimensions = data.frame(),
                     retained_log_likelihood = NA_real_,
                     optimizer_log_likelihood = NA_real_,
                     likelihood_difference = NA_real_) {
    list(
      contract_version = mfrmr_boundary_contract_version(),
      method = method,
      model = model,
      scope = paste(
        "JML GPCM joint linear-additive paths with canonical constant",
        "sum-zero log-slope rate vectors"
      ),
      state = state,
      complete = isTRUE(complete),
      scope_complete = isTRUE(scope_complete),
      structural_identification_complete = FALSE,
      target_status = target_status,
      certificates = certificates,
      direction_loadings = direction_loadings,
      rate_certificates = rate_certificates,
      rate_direction_loadings = rate_direction_loadings,
      rate_scope_state = rate_scope_state,
      rate_scope_complete = isTRUE(rate_scope_complete),
      positive_certificate_present =
        isTRUE(positive_certificate_present),
      dimensions = dimensions,
      retained_log_likelihood = as.numeric(retained_log_likelihood),
      optimizer_log_likelihood = as.numeric(optimizer_log_likelihood),
      likelihood_difference = as.numeric(likelihood_difference),
      detail = detail,
      limitations = paste(
        "The finite candidate family covers constant expanded log-slope",
        "rates through positive, zero, leading-negative, and deeper-negative",
        "groups; the existing +1/-1 pair paths are exact special cases.",
        "It does not cover curved paths or paths without a limiting rate",
        "vector. A positive competitive path is a boundary candidate, not a",
        "proof of global behavior for the non-concave GPCM likelihood. A completed",
        "negative result does not establish a finite maximum or nonlinear",
        "structural identification. MML requires a separate argument."
      )
    )
  }

  if (!identical(model, "GPCM")) {
    return(finish(
      "not_applicable_model", TRUE, FALSE,
      "The joint nonlinear slope-path audit applies only to GPCM."
    ))
  }
  if (!identical(method, "JML")) {
    return(finish(
      "not_applicable_mml", TRUE, FALSE,
      "The conditional joint JML path is not reused for marginal MML."
    ))
  }
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    return(finish(
      "not_evaluated_dependency", FALSE, FALSE,
      "Package 'lpSolve' is required for the bounded joint-path audit."
    ))
  }

  levels <- as.character(config$gpcm_spec$levels %||% character(0))
  n_levels <- length(levels)
  n_observations <- length(idx$score_k %||% integer(0))
  n_categories <- as.integer(config$n_cat %||% 0L)
  n_steps <- max(n_categories - 1L, 0L)
  n_pairs <- as.double(n_levels) * as.double(max(n_levels - 1L, 0L))
  n_rate_candidates <-
    mfrmr_jml_gpcm_joint_rate_candidate_count(n_levels)
  n_additional_rate_candidates <- max(n_rate_candidates - n_pairs, 0)
  controls <- list(
    max_observations, max_contrast_rows, max_additive_coordinates,
    max_direction_pairs, max_lp_nonzeros, lp_timeout,
    certificate_tolerance, likelihood_tolerance, max_rate_candidates
  )
  valid_controls <- all(vapply(controls, function(value) {
    length(value) == 1L && !is.na(value) && is.finite(value) && value >= 0
  }, logical(1)))
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE, FALSE,
      "Joint-path execution limits and tolerances must be finite nonnegative scalars."
    ))
  }

  params <- tryCatch(expand_params(par, sizes, config), error = function(e) e)
  if (inherits(params, "error") || n_levels < 1L ||
      length(params$slopes) != n_levels ||
      any(!is.finite(params$slopes)) || any(params$slopes <= 0)) {
    return(finish(
      "not_evaluated_mapping", FALSE, FALSE,
      if (inherits(params, "error")) {
        paste0("Joint GPCM parameter expansion failed: ",
               conditionMessage(params))
      } else {
        "Expanded slope coordinates did not match the declared levels."
      }
    ))
  }
  if (n_levels <= 1L || as.integer(sizes$log_slopes %||% 0L) == 0L) {
    return(finish(
      "no_free_log_slope_coordinates", TRUE, TRUE,
      "The unit-slope GPCM reduction has no joint log-slope pair path.",
      target_status = mfrmr_jml_gpcm_joint_targets(
        levels, slope_facet, reason = "no_free_log_slope_coordinate"
      )
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
      "not_evaluated_design", FALSE, FALSE,
      paste0("Joint adjacent-design construction failed: ",
             conditionMessage(adjacent))
    ))
  }
  additive_coordinates <- ncol(adjacent$design)
  contrast_rows <- as.double(n_observations) * as.double(n_steps)
  dimensions <- data.frame(
    Observations = as.integer(n_observations),
    EffectiveObservations = NA_integer_,
    Categories = as.integer(n_categories),
    ContrastRows = as.double(contrast_rows),
    AdditiveFreeCoordinates = as.integer(additive_coordinates),
    SlopeLevels = as.integer(n_levels),
    OrderedDirectionPairs = as.double(n_pairs),
    CanonicalRateCandidates = as.double(n_rate_candidates),
    AdditionalRateCandidates = as.double(n_additional_rate_candidates),
    EvaluatedPairs = 0L,
    CertifiedPairs = 0L,
    EvaluatedRateCandidates = 0L,
    CertifiedRateCandidates = 0L,
    MaximumLPVariables = 0L,
    MaximumLPConstraints = 0L,
    MaximumLPNonzeros = 0,
    stringsAsFactors = FALSE
  )
  not_evaluated_targets <- mfrmr_jml_gpcm_joint_targets(
    levels, slope_facet,
    evaluation_state = "not_evaluated",
    reason = "joint_pair_path_not_evaluated"
  )
  not_evaluated_targets$CandidateStatus <- "not_evaluated"
  if (n_observations > max_observations ||
      contrast_rows > max_contrast_rows ||
      additive_coordinates > max_additive_coordinates ||
      n_pairs > max_direction_pairs) {
    return(finish(
      "not_evaluated_size_limit", FALSE, FALSE,
      paste0(
        "The joint GPCM path audit exceeded its bounded workload (",
        n_observations, " observations; ", contrast_rows,
        " contrasts; ", additive_coordinates,
        " additive coordinates; ", n_pairs, " ordered slope pairs)."
      ),
      target_status = not_evaluated_targets,
      dimensions = dimensions
    ))
  }
  if (additive_coordinates == 0L) {
    return(finish(
      "no_free_additive_coordinates", TRUE, TRUE,
      "No free additive coordinate was available for a joint path.",
      target_status = mfrmr_jml_gpcm_joint_targets(levels, slope_facet),
      dimensions = dimensions
    ))
  }

  score_k <- suppressWarnings(as.integer(idx$score_k))
  slope_index <- suppressWarnings(as.integer(idx$slope_idx))
  weight <- as.numeric(idx$weight %||% rep(1, n_observations))
  valid_map <- length(score_k) == n_observations &&
    length(slope_index) == n_observations &&
    length(weight) == n_observations &&
    !anyNA(score_k) && !anyNA(slope_index) && !anyNA(weight) &&
    all(score_k >= 0L & score_k < n_categories) &&
    all(slope_index >= 1L & slope_index <= n_levels) &&
    all(is.finite(weight)) && all(weight >= 0)
  if (!valid_map || n_observations < 1L || n_categories < 2L) {
    return(finish(
      "not_evaluated_design", FALSE, FALSE,
      "The retained score, slope-index, weight, or category map was malformed.",
      target_status = not_evaluated_targets,
      dimensions = dimensions
    ))
  }
  dimensions$EffectiveObservations <- sum(weight > 0)
  if (!any(weight > 0)) {
    return(finish(
      "not_evaluated_no_positive_weight", FALSE, FALSE,
      "No retained observation had positive likelihood weight.",
      target_status = not_evaluated_targets,
      dimensions = dimensions
    ))
  }

  contrast <- tryCatch(
    mfrmr_jml_observed_contrast_design(
      adjacent$design, score_k, n_observations, n_steps
    ),
    error = function(e) e
  )
  built <- tryCatch(
    mfrmr_gpcm_response_kernel_design(prep, idx, config, sizes, par),
    error = function(e) e
  )
  if (inherits(contrast, "error") || inherits(built, "error")) {
    failure <- if (inherits(contrast, "error")) contrast else built
    return(finish(
      "not_evaluated_design", FALSE, FALSE,
      paste0("Joint GPCM response construction failed: ",
             conditionMessage(failure)),
      target_status = not_evaluated_targets,
      dimensions = dimensions
    ))
  }
  optimizer_index <- as.integer(adjacent$map$OptimizerIndex)
  if (length(optimizer_index) != additive_coordinates ||
      anyNA(optimizer_index) || anyDuplicated(optimizer_index)) {
    return(finish(
      "not_evaluated_mapping", FALSE, FALSE,
      "The additive optimizer-coordinate map was incomplete or duplicated.",
      target_status = not_evaluated_targets,
      dimensions = dimensions
    ))
  }

  adjacent_base <- matrix(
    as.numeric(built$eta_minus_step),
    nrow = n_observations, ncol = n_steps
  )
  utilities <- matrix(0, nrow = n_observations, ncol = n_categories)
  for (transition in seq_len(n_steps)) {
    utilities[, transition + 1L] <-
      utilities[, transition] + adjacent_base[, transition]
  }
  scaled <- utilities * matrix(
    params$slopes[slope_index],
    nrow = n_observations, ncol = n_categories
  )
  row_max <- apply(scaled, 1L, max)
  log_denominator <- row_max + log(rowSums(exp(
    scaled - matrix(row_max, nrow = n_observations,
                    ncol = n_categories)
  )))
  observed_scaled <- scaled[cbind(seq_len(n_observations), score_k + 1L)]
  retained_log_likelihood <- sum(weight * (observed_scaled - log_denominator))
  optimizer_log_likelihood <- tryCatch(
    -mfrm_loglik_jml(par, idx, config, sizes),
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
      "The reconstructed retained likelihood did not match the optimizer objective.",
      target_status = not_evaluated_targets,
      dimensions = dimensions,
      retained_log_likelihood = retained_log_likelihood,
      optimizer_log_likelihood = optimizer_log_likelihood,
      likelihood_difference = likelihood_difference
    ))
  }

  contrast_observation <- rep(seq_len(n_observations), each = n_steps)
  pair_grid <- expand.grid(
    PositiveIndex = seq_len(n_levels),
    NegativeIndex = seq_len(n_levels),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pair_grid <- pair_grid[
    pair_grid$PositiveIndex != pair_grid$NegativeIndex, , drop = FALSE
  ]
  pair_rows <- vector("list", nrow(pair_grid))
  pair_directions <- vector("list", nrow(pair_grid))
  for (pair_position in seq_len(nrow(pair_grid))) {
    positive_index <- pair_grid$PositiveIndex[pair_position]
    negative_index <- pair_grid$NegativeIndex[pair_position]
    solved <- mfrmr_jml_gpcm_joint_pair_lp(
      contrast = contrast,
      contrast_observation = contrast_observation,
      slope_index = slope_index,
      weight = weight,
      n_categories = n_categories,
      positive_index = positive_index,
      negative_index = negative_index,
      max_lp_nonzeros = max_lp_nonzeros,
      lp_timeout = lp_timeout,
      certificate_tolerance = certificate_tolerance
    )
    dimensions$MaximumLPVariables <- max(
      dimensions$MaximumLPVariables, solved$lp_variables
    )
    dimensions$MaximumLPConstraints <- max(
      dimensions$MaximumLPConstraints, solved$lp_constraints
    )
    dimensions$MaximumLPNonzeros <- max(
      dimensions$MaximumLPNonzeros, solved$lp_nonzeros
    )
    boundary <- list(valid = FALSE, log_likelihood = NA_real_)
    if (isTRUE(solved$certified)) {
      adjacent_direction <- matrix(
        as.numeric(adjacent$design %*% solved$direction),
        nrow = n_observations, ncol = n_steps
      )
      category_direction <- matrix(
        0, nrow = n_observations, ncol = n_categories
      )
      for (transition in seq_len(n_steps)) {
        category_direction[, transition + 1L] <-
          category_direction[, transition] +
          adjacent_direction[, transition]
      }
      boundary <- mfrmr_jml_gpcm_joint_boundary_limit(
        utilities = utilities,
        category_direction = category_direction,
        score_k = score_k,
        slope_index = slope_index,
        slopes = params$slopes,
        weight = weight,
        positive_index = positive_index,
        negative_index = negative_index,
        tolerance = certificate_tolerance
      )
    }
    improvement <- boundary$log_likelihood - retained_log_likelihood
    competitive <- isTRUE(boundary$valid) && is.finite(improvement) &&
      improvement >= -likelihood_tolerance * likelihood_scale
    certified <- isTRUE(solved$certified) && isTRUE(competitive)
    reason <- if (isTRUE(certified)) {
      "certified_competitive_joint_nonlinear_boundary_candidate"
    } else if (isTRUE(solved$certified)) {
      "joint_tail_boundary_not_competitive_with_retained_fit"
    } else {
      solved$reason
    }
    pair_rows[[pair_position]] <- data.frame(
      PairId = sprintf("JP%05d", pair_position),
      PositiveLevel = levels[positive_index],
      NegativeLevel = levels[negative_index],
      PositiveIndex = as.integer(positive_index),
      NegativeIndex = as.integer(negative_index),
      SolverStatus = as.integer(solved$solver_status),
      MarginCapacity = as.numeric(solved$margin_capacity),
      PositiveMinimumMargin = as.numeric(solved$positive_minimum_margin),
      NeutralMinimumMargin = as.numeric(solved$neutral_minimum_margin),
      NegativeLeadingCoefficient =
        as.numeric(solved$negative_leading_coefficient),
      AdditiveDirectionL1 = as.numeric(solved$direction_l1),
      CurrentLogLikelihood = retained_log_likelihood,
      BoundaryLogLikelihood = as.numeric(boundary$log_likelihood),
      BoundaryImprovement = as.numeric(improvement),
      AnalyticTailCertified = isTRUE(solved$certified),
      CompetitiveBoundary = isTRUE(competitive),
      Certified = isTRUE(certified),
      EvaluationState = if (isTRUE(solved$evaluated)) {
        "evaluated"
      } else {
        "not_evaluated"
      },
      ReasonCodes = reason,
      stringsAsFactors = FALSE
    )
    pair_directions[[pair_position]] <- solved$direction
  }
  certificates <- do.call(rbind, pair_rows)
  evaluated <- certificates$EvaluationState == "evaluated"
  certified_rows <- which(certificates$Certified)
  dimensions$EvaluatedPairs <- sum(evaluated)
  dimensions$CertifiedPairs <- length(certified_rows)
  if (any(!evaluated)) {
    return(finish(
      "not_evaluated_pair_solver", FALSE, FALSE,
      "At least one ordered joint GPCM pair path was not evaluated.",
      target_status = not_evaluated_targets,
      certificates = certificates,
      dimensions = dimensions,
      retained_log_likelihood = retained_log_likelihood,
      optimizer_log_likelihood = optimizer_log_likelihood,
      likelihood_difference = likelihood_difference
    ))
  }

  rate_certificates <- empty_rates
  rate_directions <- list()
  rate_scope_complete <- TRUE
  rate_scope_state <- "no_additional_rate_candidates"
  if (n_additional_rate_candidates > max_rate_candidates) {
    rate_scope_complete <- FALSE
    rate_scope_state <- "not_evaluated_size_limit"
  } else if (n_additional_rate_candidates > 0) {
    rate_candidates <- mfrmr_jml_gpcm_joint_rate_candidates(n_levels)
    rate_rows <- vector("list", length(rate_candidates))
    rate_directions <- vector("list", length(rate_candidates))
    for (rate_position in seq_along(rate_candidates)) {
      candidate <- rate_candidates[[rate_position]]
      solved <- mfrmr_jml_gpcm_joint_rate_lp(
        contrast = contrast,
        contrast_observation = contrast_observation,
        slope_index = slope_index,
        weight = weight,
        n_categories = n_categories,
        positive_indices = candidate$positive,
        zero_indices = candidate$zero,
        leading_negative_indices = candidate$leading_negative,
        max_lp_nonzeros = max_lp_nonzeros,
        lp_timeout = lp_timeout,
        certificate_tolerance = certificate_tolerance
      )
      dimensions$MaximumLPVariables <- max(
        dimensions$MaximumLPVariables, solved$lp_variables
      )
      dimensions$MaximumLPConstraints <- max(
        dimensions$MaximumLPConstraints, solved$lp_constraints
      )
      dimensions$MaximumLPNonzeros <- max(
        dimensions$MaximumLPNonzeros, solved$lp_nonzeros
      )
      boundary <- list(valid = FALSE, log_likelihood = NA_real_)
      if (isTRUE(solved$certified)) {
        adjacent_direction <- matrix(
          as.numeric(adjacent$design %*% solved$direction),
          nrow = n_observations, ncol = n_steps
        )
        category_direction <- matrix(
          0, nrow = n_observations, ncol = n_categories
        )
        for (transition in seq_len(n_steps)) {
          category_direction[, transition + 1L] <-
            category_direction[, transition] +
            adjacent_direction[, transition]
        }
        boundary <- mfrmr_jml_gpcm_joint_rate_boundary_limit(
          utilities = utilities,
          category_direction = category_direction,
          score_k = score_k,
          slope_index = slope_index,
          slopes = params$slopes,
          weight = weight,
          positive_indices = candidate$positive,
          zero_indices = candidate$zero,
          negative_indices = c(
            candidate$leading_negative, candidate$deeper_negative
          ),
          tolerance = certificate_tolerance
        )
      }
      improvement <- boundary$log_likelihood - retained_log_likelihood
      competitive <- isTRUE(boundary$valid) && is.finite(improvement) &&
        improvement >= -likelihood_tolerance * likelihood_scale
      certified <- isTRUE(solved$certified) && isTRUE(competitive)
      reason <- if (isTRUE(certified)) {
        "certified_competitive_joint_constant_rate_boundary_candidate"
      } else if (isTRUE(solved$certified)) {
        "joint_constant_rate_boundary_not_competitive_with_retained_fit"
      } else {
        solved$reason
      }
      collapse_levels <- function(indices) {
        if (length(indices) == 0L) "" else paste(levels[indices], collapse = ";")
      }
      rate_rows[[rate_position]] <- data.frame(
        RateId = sprintf("JR%05d", rate_position),
        PositiveLevels = collapse_levels(candidate$positive),
        ZeroLevels = collapse_levels(candidate$zero),
        LeadingNegativeLevels =
          collapse_levels(candidate$leading_negative),
        DeeperNegativeLevels = collapse_levels(candidate$deeper_negative),
        ExpandedRates = paste(
          formatC(candidate$rates, digits = 16L, format = "fg"),
          collapse = ";"
        ),
        SolverStatus = as.integer(solved$solver_status),
        MarginCapacity = as.numeric(solved$margin_capacity),
        PositiveMinimumMargin =
          as.numeric(solved$positive_minimum_margin),
        ZeroMinimumMargin = as.numeric(solved$neutral_minimum_margin),
        LeadingNegativeCoefficient =
          as.numeric(solved$negative_leading_coefficient),
        AdditiveDirectionL1 = as.numeric(solved$direction_l1),
        CurrentLogLikelihood = retained_log_likelihood,
        BoundaryLogLikelihood = as.numeric(boundary$log_likelihood),
        BoundaryImprovement = as.numeric(improvement),
        AnalyticTailCertified = isTRUE(solved$certified),
        CompetitiveBoundary = isTRUE(competitive),
        Certified = isTRUE(certified),
        EvaluationState = if (isTRUE(solved$evaluated)) {
          "evaluated"
        } else {
          "not_evaluated"
        },
        ReasonCodes = reason,
        stringsAsFactors = FALSE
      )
      rate_directions[[rate_position]] <- solved$direction
    }
    rate_certificates <- do.call(rbind, rate_rows)
    rate_evaluated <- rate_certificates$EvaluationState == "evaluated"
    dimensions$EvaluatedRateCandidates <- sum(rate_evaluated)
    dimensions$CertifiedRateCandidates <-
      sum(rate_certificates$Certified)
    rate_scope_complete <- all(rate_evaluated)
    rate_scope_state <- if (!rate_scope_complete) {
      "not_evaluated_rate_solver"
    } else if (any(rate_certificates$Certified)) {
      "certified_competitive_constant_rate_boundary_path"
    } else {
      "none_certified"
    }
  }

  rate_certified_rows <- which(rate_certificates$Certified)
  split_levels <- function(values) {
    values <- as.character(values)
    unlist(strsplit(values[nzchar(values)], ";", fixed = TRUE),
           use.names = FALSE)
  }
  rate_positive_levels <- split_levels(
    rate_certificates$PositiveLevels[rate_certified_rows]
  )
  rate_negative_levels <- c(
    split_levels(
      rate_certificates$LeadingNegativeLevels[rate_certified_rows]
    ),
    split_levels(
      rate_certificates$DeeperNegativeLevels[rate_certified_rows]
    )
  )
  positive <- seq_len(n_levels) %in%
    certificates$PositiveIndex[certified_rows] |
    levels %in% rate_positive_levels
  negative <- seq_len(n_levels) %in%
    certificates$NegativeIndex[certified_rows] |
    levels %in% rate_negative_levels
  targets <- mfrmr_jml_gpcm_joint_targets(
    levels, slope_facet, positive = positive, negative = negative
  )

  slices <- build_param_slices(sizes)
  optimizer_log_slope_index <- as.integer(slices$log_slopes %||% integer(0))
  loading_rows <- list()
  rate_loading_rows <- list()
  loading_cursor <- 0L
  rate_loading_cursor <- 0L
  for (certificate_row in certified_rows) {
    pair <- certificates[certificate_row, , drop = FALSE]
    direction <- pair_directions[[certificate_row]]
    additive_keep <- which(abs(direction) > certificate_tolerance)
    if (length(additive_keep) > 0L) {
      loading_cursor <- loading_cursor + 1L
      loading_rows[[loading_cursor]] <- data.frame(
        PairId = pair$PairId,
        CoordinateType = "optimizer_additive",
        OptimizerIndex = optimizer_index[additive_keep],
        Coordinate = as.character(adjacent$map$Coordinate[additive_keep]),
        Level = as.character(adjacent$map$Level[additive_keep]),
        Loading = direction[additive_keep],
        stringsAsFactors = FALSE
      )
    }
    expanded_direction <- rep(0, n_levels)
    expanded_direction[pair$PositiveIndex] <- 1
    expanded_direction[pair$NegativeIndex] <- -1
    expanded_keep <- which(expanded_direction != 0)
    loading_cursor <- loading_cursor + 1L
    loading_rows[[loading_cursor]] <- data.frame(
      PairId = pair$PairId,
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
        PairId = pair$PairId,
        CoordinateType = "optimizer_log_slope",
        OptimizerIndex = optimizer_log_slope_index[free_keep],
        Coordinate = paste0("log_slopes:", levels[free_keep]),
        Level = levels[free_keep],
        Loading = free_direction[free_keep],
        stringsAsFactors = FALSE
      )
    }
  }
  for (certificate_row in rate_certified_rows) {
    rate <- rate_certificates[certificate_row, , drop = FALSE]
    candidate <- rate_candidates[[certificate_row]]
    direction <- rate_directions[[certificate_row]]
    candidate_rows <- list()
    candidate_cursor <- 0L
    additive_keep <- which(abs(direction) > certificate_tolerance)
    if (length(additive_keep) > 0L) {
      candidate_cursor <- candidate_cursor + 1L
      candidate_rows[[candidate_cursor]] <- data.frame(
        PairId = rate$RateId,
        CoordinateType = "optimizer_additive",
        OptimizerIndex = optimizer_index[additive_keep],
        Coordinate = as.character(adjacent$map$Coordinate[additive_keep]),
        Level = as.character(adjacent$map$Level[additive_keep]),
        Loading = direction[additive_keep],
        stringsAsFactors = FALSE
      )
    }
    expanded_keep <- which(candidate$rates != 0)
    candidate_cursor <- candidate_cursor + 1L
    candidate_rows[[candidate_cursor]] <- data.frame(
      PairId = rate$RateId,
      CoordinateType = "expanded_log_slope",
      OptimizerIndex = rep(NA_integer_, length(expanded_keep)),
      Coordinate = paste0("expanded_log_slope:", levels[expanded_keep]),
      Level = levels[expanded_keep],
      Loading = candidate$rates[expanded_keep],
      stringsAsFactors = FALSE
    )
    free_direction <- candidate$rates[seq_len(n_levels - 1L)]
    free_keep <- which(free_direction != 0)
    if (length(free_keep) > 0L) {
      candidate_cursor <- candidate_cursor + 1L
      candidate_rows[[candidate_cursor]] <- data.frame(
        PairId = rate$RateId,
        CoordinateType = "optimizer_log_slope",
        OptimizerIndex = optimizer_log_slope_index[free_keep],
        Coordinate = paste0("log_slopes:", levels[free_keep]),
        Level = levels[free_keep],
        Loading = free_direction[free_keep],
        stringsAsFactors = FALSE
      )
    }
    candidate_loading <- do.call(rbind, candidate_rows)
    rate_loading_cursor <- rate_loading_cursor + 1L
    rate_loading_rows[[rate_loading_cursor]] <- candidate_loading
    loading_cursor <- loading_cursor + 1L
    loading_rows[[loading_cursor]] <- candidate_loading
  }
  direction_loadings <- if (length(loading_rows) > 0L) {
    do.call(rbind, loading_rows)
  } else {
    empty_loadings
  }
  rate_direction_loadings <- if (length(rate_loading_rows) > 0L) {
    do.call(rbind, rate_loading_rows)
  } else {
    empty_loadings
  }

  positive_count <- length(certified_rows) + length(rate_certified_rows)
  overall_complete <- isTRUE(rate_scope_complete)
  state <- if (positive_count > 0L) {
    "certified_competitive_joint_boundary_path"
  } else if (!overall_complete) {
    if (identical(rate_scope_state, "not_evaluated_size_limit")) {
      "not_evaluated_rate_size_limit"
    } else {
      "not_evaluated_rate_solver"
    }
  } else {
    "none_certified"
  }

  finish(
    state = state,
    complete = overall_complete,
    scope_complete = overall_complete,
    detail = if (positive_count > 0L) {
      paste0(
        "Certified ", length(certified_rows), " of ", nrow(certificates),
        " ordered slope pairs and ", length(rate_certified_rows), " of ",
        nrow(rate_certificates), " additional canonical constant-rate",
        " paths as competitive joint nonlinear boundary candidates."
      )
    } else if (!overall_complete) {
      paste0(
        "The ordered slope-pair family completed without a certificate, but",
        " the additional canonical constant-rate family was not completed",
        " (state: ", rate_scope_state, ")."
      )
    } else {
      paste0(
        "No competitive path was certified among ", nrow(certificates),
        " ordered joint additive/log-slope pairs or ",
        nrow(rate_certificates),
        " additional canonical constant-rate paths."
      )
    },
    target_status = targets,
    certificates = certificates,
    direction_loadings = direction_loadings,
    rate_certificates = rate_certificates,
    rate_direction_loadings = rate_direction_loadings,
    rate_scope_state = rate_scope_state,
    rate_scope_complete = rate_scope_complete,
    positive_certificate_present = positive_count > 0L,
    dimensions = dimensions,
    retained_log_likelihood = retained_log_likelihood,
    optimizer_log_likelihood = optimizer_log_likelihood,
    likelihood_difference = likelihood_difference
  )
}

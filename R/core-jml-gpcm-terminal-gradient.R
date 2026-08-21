# Fixed-objective JML GPCM terminal-gradient stability audit
# ==============================================================================
#
# A small gradient at a finite optimizer iterate is only a retained-point
# numerical property. In a separated or receding GPCM likelihood, finite
# iterates can have small gradients while the supremum remains on a boundary.
# This audit therefore binds the gradient to the exact unpenalized JML
# objective, verifies selected coordinates by central differences, reports
# optimizer blocks separately, and gives positive boundary certificates
# precedence over finite-point stationarity language.

mfrmr_jml_gpcm_terminal_gradient_contract_version <- function() {
  "mfrmr-jml-gpcm-fixed-objective-terminal-gradient-0.2.3-v1"
}

mfrmr_jml_gpcm_terminal_gradient_empty_block_table <- function() {
  data.frame(
    Block = character(0),
    FreeCoordinates = integer(0),
    GradientSupNorm = numeric(0),
    GradientRMS = numeric(0),
    WithinImplementationTolerance = logical(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_terminal_gradient_empty_probe_table <- function() {
  data.frame(
    OptimizerIndex = integer(0),
    Block = character(0),
    ParameterValue = numeric(0),
    Step = numeric(0),
    AnalyticGradient = numeric(0),
    NumericGradient = numeric(0),
    AbsoluteDifference = numeric(0),
    ScaledDifference = numeric(0),
    Evaluated = logical(0),
    stringsAsFactors = FALSE
  )
}

mfrmr_jml_gpcm_terminal_gradient_block_labels <- function(sizes) {
  unlist(lapply(names(sizes), function(block) {
    rep(as.character(block), as.integer(sizes[[block]] %||% 0L))
  }), use.names = FALSE)
}

mfrmr_jml_gpcm_terminal_gradient_probe_indices <- function(gradient,
                                                             sizes,
                                                             max_probes) {
  gradient <- as.numeric(gradient)
  max_probes <- as.integer(max_probes)
  if (length(gradient) == 0L || max_probes < 1L) return(integer(0))
  slices <- build_param_slices(sizes)
  selected <- integer(0)

  slope_indices <- as.integer(slices$log_slopes %||% integer(0))
  if (length(slope_indices) > 0L) {
    selected <- c(selected, slope_indices[seq_len(min(
      length(slope_indices), max_probes
    ))])
  }
  for (block in names(slices)) {
    indices <- as.integer(slices[[block]] %||% integer(0))
    if (length(indices) == 0L || length(selected) >= max_probes) next
    block_gradient <- abs(gradient[indices])
    if (all(!is.finite(block_gradient))) next
    selected <- c(selected, indices[which.max(block_gradient)])
    selected <- unique(selected)
  }
  if (length(selected) < max_probes && any(is.finite(gradient))) {
    order_all <- order(abs(gradient), decreasing = TRUE, na.last = NA)
    selected <- unique(c(selected, order_all))
  }
  as.integer(selected[seq_len(min(length(selected), max_probes))])
}

mfrmr_jml_gpcm_terminal_gradient_block_summary <- function(gradient,
                                                             sizes,
                                                             tolerance) {
  gradient <- as.numeric(gradient)
  labels <- mfrmr_jml_gpcm_terminal_gradient_block_labels(sizes)
  if (length(gradient) == 0L || length(labels) != length(gradient)) {
    return(mfrmr_jml_gpcm_terminal_gradient_empty_block_table())
  }
  rows <- lapply(unique(labels), function(block) {
    values <- gradient[labels == block]
    finite <- values[is.finite(values)]
    sup_norm <- if (length(finite) > 0L) max(abs(finite)) else NA_real_
    rms <- if (length(finite) > 0L) sqrt(mean(finite^2)) else NA_real_
    data.frame(
      Block = block,
      FreeCoordinates = as.integer(length(values)),
      GradientSupNorm = sup_norm,
      GradientRMS = rms,
      WithinImplementationTolerance = is.finite(sup_norm) &&
        is.finite(tolerance) && sup_norm <= tolerance,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

audit_mfrm_jml_gpcm_terminal_gradient <- function(
    idx,
    config,
    sizes,
    par,
    opt,
    boundary_classification =
      config$boundary_audit$gpcm_fixed_objective_classification %||% list(),
    max_free_dimension = 250L,
    max_numeric_probes = 8L,
    numeric_relative_step = 1e-5,
    derivative_tolerance = 1e-5,
    objective_tolerance = 1e-8,
    diagnostic_tolerance = 1e-10) {
  method <- as.character(config$method %||% "unknown")[1]
  model <- as.character(config$model %||% "unknown")[1]
  free_dimension <- sum(vapply(sizes, as.integer, integer(1)))
  empty_blocks <- mfrmr_jml_gpcm_terminal_gradient_empty_block_table()
  empty_probes <- mfrmr_jml_gpcm_terminal_gradient_empty_probe_table()

  finish <- function(state, evaluated, detail,
                     reason_codes,
                     retained_objective = NA_real_,
                     optimizer_objective = NA_real_,
                     objective_difference = NA_real_,
                     gradient_sup_norm = NA_real_,
                     gradient_rms = NA_real_,
                     stored_gradient_sup_norm = NA_real_,
                     stored_gradient_rms = NA_real_,
                     implementation_tolerance = NA_real_,
                     block_summary = empty_blocks,
                     numeric_probes = empty_probes,
                     objective_reconstruction_agrees = FALSE,
                     stored_gradient_agrees = FALSE,
                     numeric_gradient_agrees = FALSE,
                     tolerance_met = FALSE,
                     boundary_precedence = FALSE,
                     retained_point_stationarity_supported = FALSE,
                     optimizer_code_zero = FALSE,
                     polish_history_state = "not_evaluated") {
    list(
      contract_version =
        mfrmr_jml_gpcm_terminal_gradient_contract_version(),
      method = method,
      model = model,
      estimator_identity = "unpenalized_fixed_effects_jml_no_finite_box",
      objective_identity = "identified_conditional_joint_log_likelihood",
      coordinate_system = "identified_optimizer_free_coordinates",
      state = state,
      evaluated = isTRUE(evaluated),
      free_dimension = as.integer(free_dimension),
      retained_objective = as.numeric(retained_objective),
      optimizer_objective = as.numeric(optimizer_objective),
      objective_difference = as.numeric(objective_difference),
      gradient_sup_norm = as.numeric(gradient_sup_norm),
      gradient_rms = as.numeric(gradient_rms),
      stored_gradient_sup_norm = as.numeric(stored_gradient_sup_norm),
      stored_gradient_rms = as.numeric(stored_gradient_rms),
      implementation_tolerance = as.numeric(implementation_tolerance),
      implementation_tolerance_basis =
        "existing_common_terminal_gradient_gate_max_1e-4_or_10x_reltol",
      criterion_frozen_for_scientific_claims = FALSE,
      objective_reconstruction_agrees =
        isTRUE(objective_reconstruction_agrees),
      stored_gradient_agrees = isTRUE(stored_gradient_agrees),
      numeric_gradient_agrees = isTRUE(numeric_gradient_agrees),
      implementation_tolerance_met = isTRUE(tolerance_met),
      boundary_precedence = isTRUE(boundary_precedence),
      optimizer_code_zero = isTRUE(optimizer_code_zero),
      retained_point_stationarity_supported =
        isTRUE(retained_point_stationarity_supported),
      finite_interior_stationarity_certified = FALSE,
      global_finite_maximum_certified = FALSE,
      global_boundary_absence_certified = FALSE,
      standard_error_eligible = FALSE,
      confidence_interval_eligible = FALSE,
      external_comparison_eligible = FALSE,
      block_summary = block_summary,
      numeric_probes = numeric_probes,
      numeric_probe_relative_step = as.numeric(numeric_relative_step),
      numeric_probe_tolerance = as.numeric(derivative_tolerance),
      objective_tolerance = as.numeric(objective_tolerance),
      diagnostic_tolerance = as.numeric(diagnostic_tolerance),
      polish_history_state = polish_history_state,
      readiness_effect = "none_diagnostic_only",
      reason_codes = reason_codes,
      detail = detail,
      limitations = paste(
        "This audit classifies the analytic gradient only at the retained",
        "finite optimizer point under the identified unpenalized JML",
        "objective. Its implementation tolerance is not a frozen scientific",
        "criterion. A small gradient does not establish an interior finite",
        "maximum, global boundary absence, uncertainty, or external",
        "comparison eligibility, especially when a recession or competitive",
        "joint boundary path is certified."
      )
    )
  }

  controls <- c(
    max_free_dimension, max_numeric_probes, numeric_relative_step,
    derivative_tolerance, objective_tolerance, diagnostic_tolerance
  )
  valid_controls <- all(is.finite(controls)) && all(controls >= 0) &&
    max_free_dimension >= 1L && max_numeric_probes >= 1L &&
    numeric_relative_step > 0 && derivative_tolerance > 0 &&
    objective_tolerance > 0 && diagnostic_tolerance > 0
  if (!valid_controls) {
    return(finish(
      "not_evaluated_control", FALSE,
      "Terminal-gradient audit controls must be finite positive scalars.",
      "terminal_gradient_control_invalid"
    ))
  }
  if (!identical(model, "GPCM")) {
    out <- finish(
      "not_applicable_model", FALSE,
      "The fixed-objective terminal-gradient audit applies only to GPCM.",
      "terminal_gradient_not_applicable_model"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable"
    return(out)
  }
  if (!identical(method, "JML")) {
    out <- finish(
      "not_applicable_estimator", FALSE,
      "The conditional fixed-effects JML gradient audit is not reused for MML.",
      "terminal_gradient_not_applicable_estimator"
    )
    out$estimator_identity <- "not_applicable"
    out$objective_identity <- "not_applicable_marginal_estimator"
    return(out)
  }
  if (as.integer(config$gpcm_spec$n_params %||% 0L) == 0L) {
    return(finish(
      "not_required_unit_slope", FALSE,
      "The exact unit-slope reduction has no free GPCM log-slope coordinate.",
      "no_free_log_slope_coordinate"
    ))
  }
  if (free_dimension > as.integer(max_free_dimension)) {
    return(finish(
      "not_evaluated_size_limit", FALSE,
      paste0(
        "The fixed-objective numeric gradient oracle is limited to ",
        as.integer(max_free_dimension), " free coordinates; this fit has ",
        free_dimension, "."
      ),
      "terminal_gradient_free_dimension_limit"
    ))
  }
  if (is.null(par) || length(par) != free_dimension || any(!is.finite(par))) {
    return(finish(
      "not_evaluated_parameter_vector", FALSE,
      "The retained optimizer vector was unavailable, non-finite, or malformed.",
      "terminal_gradient_parameter_vector_invalid"
    ))
  }

  diagnostics <- opt$optimizer_diagnostics %||% list()
  optimizer_objective <- suppressWarnings(as.numeric(opt$value %||% NA_real_)[1])
  stored_sup <- suppressWarnings(as.numeric(
    diagnostics$TerminalGradientSupNorm %||% NA_real_
  )[1])
  stored_rms <- suppressWarnings(as.numeric(
    diagnostics$TerminalGradientRMS %||% NA_real_
  )[1])
  implementation_tolerance <- suppressWarnings(as.numeric(
    diagnostics$GradientReviewTolerance %||% NA_real_
  )[1])
  convergence_code <- suppressWarnings(as.integer(
    diagnostics$ConvergenceCode %||% opt$convergence %||% NA_integer_
  )[1])
  optimizer_code_zero <- is.finite(convergence_code) && convergence_code == 0L

  reconstructed <- tryCatch(
    list(
      objective = mfrm_loglik_jml(par, idx, config, sizes),
      gradient = mfrm_grad_jml(par, idx, config, sizes)
    ),
    error = function(e) e
  )
  if (inherits(reconstructed, "error")) {
    return(finish(
      "indeterminate_numerical", FALSE,
      paste0(
        "The fixed JML objective or analytic gradient could not be",
        " reconstructed: ", conditionMessage(reconstructed)
      ),
      "terminal_gradient_reconstruction_failed",
      optimizer_objective = optimizer_objective,
      stored_gradient_sup_norm = stored_sup,
      stored_gradient_rms = stored_rms,
      implementation_tolerance = implementation_tolerance,
      optimizer_code_zero = optimizer_code_zero
    ))
  }
  objective <- as.numeric(reconstructed$objective)[1]
  gradient <- as.numeric(reconstructed$gradient)
  metrics <- compute_gradient_metrics(gradient)
  objective_difference <- objective - optimizer_objective
  objective_scale <- max(1, abs(objective), abs(optimizer_objective))
  objective_agrees <- is.finite(objective) && is.finite(optimizer_objective) &&
    is.finite(objective_difference) &&
    abs(objective_difference) <= objective_tolerance * objective_scale
  gradient_valid <- length(gradient) == free_dimension &&
    all(is.finite(gradient)) &&
    is.finite(metrics$TerminalGradientSupNorm) &&
    is.finite(metrics$TerminalGradientRMS)
  diagnostic_scale <- max(
    1, abs(metrics$TerminalGradientSupNorm),
    abs(metrics$TerminalGradientRMS), abs(stored_sup), abs(stored_rms),
    na.rm = TRUE
  )
  stored_agrees <- gradient_valid && is.finite(stored_sup) &&
    is.finite(stored_rms) &&
    abs(metrics$TerminalGradientSupNorm - stored_sup) <=
      diagnostic_tolerance * diagnostic_scale &&
    abs(metrics$TerminalGradientRMS - stored_rms) <=
      diagnostic_tolerance * diagnostic_scale

  block_summary <- mfrmr_jml_gpcm_terminal_gradient_block_summary(
    gradient, sizes, implementation_tolerance
  )
  block_labels <- mfrmr_jml_gpcm_terminal_gradient_block_labels(sizes)
  probe_indices <- mfrmr_jml_gpcm_terminal_gradient_probe_indices(
    gradient, sizes, max_numeric_probes
  )
  probe_rows <- lapply(probe_indices, function(index) {
    step <- numeric_relative_step * max(1, abs(par[index]))
    plus <- par
    minus <- par
    plus[index] <- plus[index] + step
    minus[index] <- minus[index] - step
    values <- tryCatch(
      c(
        mfrm_loglik_jml(plus, idx, config, sizes),
        mfrm_loglik_jml(minus, idx, config, sizes)
      ),
      error = function(e) c(NA_real_, NA_real_)
    )
    numeric_gradient <- if (all(is.finite(values)) && is.finite(step) &&
                            step > 0) {
      (values[1] - values[2]) / (2 * step)
    } else {
      NA_real_
    }
    difference <- abs(numeric_gradient - gradient[index])
    scaled <- difference / max(
      1, abs(numeric_gradient), abs(gradient[index]), na.rm = TRUE
    )
    data.frame(
      OptimizerIndex = as.integer(index),
      Block = as.character(block_labels[index]),
      ParameterValue = as.numeric(par[index]),
      Step = as.numeric(step),
      AnalyticGradient = as.numeric(gradient[index]),
      NumericGradient = as.numeric(numeric_gradient),
      AbsoluteDifference = as.numeric(difference),
      ScaledDifference = as.numeric(scaled),
      Evaluated = is.finite(numeric_gradient) && is.finite(difference) &&
        is.finite(scaled),
      stringsAsFactors = FALSE
    )
  })
  numeric_probes <- if (length(probe_rows) > 0L) {
    do.call(rbind, probe_rows)
  } else {
    empty_probes
  }
  numeric_agrees <- nrow(numeric_probes) > 0L &&
    all(numeric_probes$Evaluated) &&
    all(numeric_probes$ScaledDifference <= derivative_tolerance)

  polish <- opt$optimizer_polish %||% list()
  stages <- as.data.frame(polish$Stages %||% data.frame(),
                          stringsAsFactors = FALSE)
  polish_history_state <- if (nrow(stages) == 0L) {
    "not_evaluated"
  } else if (!"Selected" %in% names(stages) ||
             sum(as.logical(stages$Selected), na.rm = TRUE) != 1L) {
    "inconsistent"
  } else {
    selected <- stages[which(as.logical(stages$Selected))[1], , drop = FALSE]
    selected_objective <- as.numeric(selected$Objective %||% NA_real_)[1]
    selected_gradient <- as.numeric(
      selected$TerminalGradientSupNorm %||% NA_real_
    )[1]
    if (is.finite(selected_objective) && is.finite(selected_gradient) &&
        abs(selected_objective - optimizer_objective) <=
          objective_tolerance * objective_scale &&
        abs(selected_gradient - stored_sup) <=
          diagnostic_tolerance * diagnostic_scale) {
      "consistent"
    } else {
      "inconsistent"
    }
  }

  coherent <- objective_agrees && stored_agrees && numeric_agrees &&
    gradient_valid && is.finite(implementation_tolerance) &&
    implementation_tolerance > 0 &&
    !identical(polish_history_state, "inconsistent")
  tolerance_met <- coherent &&
    metrics$TerminalGradientSupNorm <= implementation_tolerance
  boundary_contract_matches <- identical(
    as.character(boundary_classification$contract_version %||% "")[1],
    mfrmr_jml_gpcm_boundary_classification_contract_version()
  ) && identical(
    as.character(boundary_classification$objective_identity %||% "")[1],
    "identified_conditional_joint_log_likelihood"
  )
  boundary_positive <- boundary_contract_matches &&
    isTRUE(boundary_classification$positive_boundary_certificate_present)
  bounded_negative_complete <- boundary_contract_matches &&
    isTRUE(boundary_classification$no_path_certified_in_audited_families)

  common <- list(
    retained_objective = objective,
    optimizer_objective = optimizer_objective,
    objective_difference = objective_difference,
    gradient_sup_norm = metrics$TerminalGradientSupNorm,
    gradient_rms = metrics$TerminalGradientRMS,
    stored_gradient_sup_norm = stored_sup,
    stored_gradient_rms = stored_rms,
    implementation_tolerance = implementation_tolerance,
    block_summary = block_summary,
    numeric_probes = numeric_probes,
    objective_reconstruction_agrees = objective_agrees,
    stored_gradient_agrees = stored_agrees,
    numeric_gradient_agrees = numeric_agrees,
    tolerance_met = tolerance_met,
    boundary_precedence = boundary_positive,
    optimizer_code_zero = optimizer_code_zero,
    polish_history_state = polish_history_state
  )
  finalize <- function(state, evaluated, detail, reason_codes,
                       stationarity = FALSE) {
    do.call(finish, c(
      list(
        state = state,
        evaluated = evaluated,
        detail = detail,
        reason_codes = reason_codes,
        retained_point_stationarity_supported = stationarity
      ),
      common
    ))
  }

  if (!coherent) {
    if (boundary_positive) {
      return(finalize(
        "boundary_path_gradient_indeterminate", FALSE,
        paste(
          "A positive GPCM boundary-path certificate remains primary, but",
          "the retained-point objective/gradient record is numerically",
          "incoherent and cannot be interpreted."
        ),
        paste(
          "positive_boundary_certificate_primary",
          "terminal_gradient_fixed_objective_record_incoherent",
          sep = ";"
        )
      ))
    }
    return(finalize(
      "indeterminate_numerical", FALSE,
      paste(
        "The objective, stored analytic-gradient summary, selected-stage",
        "history, or fixed central-difference probes did not form a coherent",
        "fixed-objective terminal-gradient record."
      ),
      "terminal_gradient_fixed_objective_record_incoherent"
    ))
  }
  if (boundary_positive) {
    return(finalize(
      "boundary_path_gradient_nondecisive", TRUE,
      paste(
        "A positive GPCM boundary-path certificate takes precedence over",
        "the finite retained-point gradient. The gradient remains recorded",
        "but cannot establish a finite interior solution."
      ),
      "positive_boundary_certificate_overrides_finite_gradient"
    ))
  }
  if (!tolerance_met) {
    return(finalize(
      "terminal_gradient_exceeds_implementation_tolerance", TRUE,
      paste(
        "The reconstructed analytic terminal gradient exceeds the existing",
        "implementation review tolerance at the retained optimizer point."
      ),
      "terminal_gradient_review"
    ))
  }
  if (!optimizer_code_zero) {
    return(finalize(
      "small_gradient_optimizer_warning", TRUE,
      paste(
        "The retained-point gradient is within the implementation tolerance,",
        "but the optimizer did not return code zero."
      ),
      "small_terminal_gradient_with_nonzero_optimizer_code"
    ))
  }
  if (!bounded_negative_complete) {
    return(finalize(
      "retained_point_gradient_small_boundary_open", TRUE,
      paste(
        "The retained-point gradient is within the implementation tolerance,",
        "but the bounded GPCM path-family classification is incomplete or",
        "has a mismatched identity."
      ),
      "retained_point_gradient_small_boundary_classification_open"
    ))
  }
  finalize(
    "retained_point_gradient_within_implementation_tolerance", TRUE,
    paste(
      "The code-zero retained optimizer point has a coherent analytic",
      "gradient within the existing implementation tolerance, and both",
      "bounded GPCM path families completed without a positive certificate.",
      "This is retained-point first-order evidence only."
    ),
    paste(
      "retained_point_first_order_gradient_within_implementation_tolerance",
      "finite_global_maximum_not_certified",
      sep = ";"
    ),
    stationarity = TRUE
  )
}

# Accuracy-controlled ordinary direct MML. Only resolved fixed/adaptive modes
# reach model evaluation, information, scoring or saved calibration consumers.

mfrmr_mml_accuracy_orders <- function(max_points, initial = 31L, offset = 1L) {
  minimum <- if (offset == 1L) 7L else 3L
  if (!is.numeric(max_points) || is.complex(max_points) || length(max_points) != 1L ||
      !is.finite(max_points) || max_points != floor(max_points) ||
      max_points < minimum || max_points > .Machine$integer.max - 1) {
    stop("`max_points` must be one integer between ", minimum, " and ",
      .Machine$integer.max - 1, ".", call. = FALSE)
  }
  if (!is.numeric(initial) || is.complex(initial) || length(initial) != 1L ||
      !is.finite(initial) || initial < 2 || initial != floor(initial))
    stop("Automatic refinement requires integer `quad_points` >= 2.", call. = FALSE)
  limit <- floor((max_points - offset) / 2)
  order <- min(initial, limit)
  orders <- as.integer(order)
  while (order < limit) {
    order <- min(2 * order - 1, limit)
    orders <- c(orders, as.integer(order))
  }
  orders
}

mfrmr_mml_accuracy_control <- function(control = NULL) {
  defaults <- list(adaptive_points = c(31L, 61L, 121L),
    loglik_tol = 1e-7, posterior_tol = 1e-6, gradient_tol = 1e-4)
  if (is.null(control)) return(defaults)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) ||
       any(!names(control) %in% c(names(defaults), "max_points")) || anyDuplicated(names(control))))) {
    stop("`mml_accuracy_control` must be a named list of max_points or adaptive_points, loglik_tol, posterior_tol and/or gradient_tol.",
      call. = FALSE)
  }
  if (all(c("max_points", "adaptive_points") %in% names(control)))
    stop("Specify `max_points` or `adaptive_points`, not both.", call. = FALSE)
  capped <- "max_points" %in% names(control)
  if (capped) {
    points <- mfrmr_mml_accuracy_orders(control$max_points)
    defaults$adaptive_points <- NULL
  }
  defaults[names(control)] <- control
  if (!capped) points <- defaults$adaptive_points
  if (!is.numeric(points) || is.complex(points) || !length(points) || any(!is.finite(points)) ||
      any(points < 3 | points != floor(points) | points > (.Machine$integer.max - 1) / 2) ||
      is.unsorted(points, strictly = TRUE)) {
    stop("`adaptive_points` must contain strictly increasing integer orders >= 3.", call. = FALSE)
  }
  if (!capped) defaults$adaptive_points <- as.integer(points)
  for (name in c("loglik_tol", "posterior_tol", "gradient_tol")) {
    value <- defaults[[name]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop("`", name, "` must be a finite positive scalar.", call. = FALSE)
    }
  }
  # Preserve explicit-grid validation; a ceiling must not allocate unused rules.
  if (!capped) invisible(lapply(unique(c(points, 2L * points + 1L)), gauss_hermite_normal))
  defaults
}

mfrmr_check_mml_accuracy <- function(opt, idx, config, sizes, labels,
                                    adaptive_points, control) {
  nodes <- config$estimation_control$quad_points
  adaptive <- mfrmr_adaptive_integration(config)
  reference_nodes <- 2L * adaptive_points + 1L
  supnorm <- function(x) {
    if (any(!is.finite(x))) stop("Non-finite integration check.", call. = FALSE)
    if (length(x)) max(abs(x)) else 0
  }
  reference_config <- config
  reference_config$estimation_control$mml_integration <- "adaptive"
  # Each rule supplies moments and its moving-node gradient from one evaluation.
  evaluate <- function(order) mfrmr_make_adaptive_mml_evaluator(
    idx, reference_config, sizes, order, include_person_summaries = TRUE)(opt$par)
  low_value <- evaluate(adaptive_points)
  high_value <- evaluate(reference_nodes)
  low <- low_value$per_person
  high <- high_value$per_person
  candidate <- low
  retained_gradient <- low_value$gradient
  if (!adaptive) {
    params <- expand_params(opt$par, sizes, config)
    base <- compute_base_eta(idx, params, config)
    quad <- gauss_hermite_normal(nodes)
    bundle <- mfrm_mml_logprob_bundle(idx, config, quad, params, base,
      include_probs = TRUE, include_linear_part = identical(config$model, "GPCM"))
    posterior_bundle <- mfrm_mml_posterior_bundle(bundle)
    person <- posterior_bundle$person_bundle
    probability <- person$posterior
    center <- as.vector(probability %*% quad$nodes)
    basis <- bundle$quad_basis
    mu <- if (isTRUE(basis$transformed)) basis$mu[person$person_ids] else 0
    sigma <- if (isTRUE(basis$transformed)) basis$sigma else 1
    candidate <- data.frame(PersonIndex = person$person_ids,
      LogMarginal = person$log_marginal, EAP = mu + sigma * center,
      PosteriorSD = sigma * sqrt(rowSums(probability *
        (rep(quad$nodes, each = nrow(probability)) - center)^2)))
    candidate <- candidate[order(candidate$PersonIndex), , drop = FALSE]
    retained_gradient <- mfrm_grad_mml_core(params, base, idx, config, sizes, quad,
      logprob_bundle = bundle, posterior_bundle = posterior_bundle)
  }
  if (!identical(low$PersonIndex, seq_along(labels)) ||
      !identical(high$PersonIndex, low$PersonIndex) ||
      !identical(candidate$PersonIndex, low$PersonIndex)) {
    stop("Person integration check unavailable: incomplete Person alignment.", call. = FALSE)
  }
  # Include the retained rule and both adaptive rules; do not average away
  # inaccurate Persons or multiply their moment errors by pattern counts.
  difference <- function(field) pmax(abs(candidate[[field]] - high[[field]]),
    abs(low[[field]] - high[[field]]))
  per_person <- data.frame(Person = as.character(labels),
    LogMarginalDifference = difference("LogMarginal"),
    EAPDifference = difference("EAP"),
    PosteriorSDDifference = difference("PosteriorSD"), stringsAsFactors = FALSE)
  loglik_difference <- supnorm(per_person$LogMarginalDifference)
  eap_difference <- supnorm(per_person$EAPDifference)
  sd_difference <- supnorm(per_person$PosteriorSDDifference)
  gradient_difference <- max(supnorm(retained_gradient - high_value$gradient),
    supnorm(low_value$gradient - high_value$gradient))
  gradient_norm <- max(supnorm(retained_gradient), supnorm(high_value$gradient))
  passed <- isTRUE(opt$convergence == 0L) &&
    loglik_difference <= control$loglik_tol &&
    max(eap_difference, sd_difference) <= control$posterior_tol &&
    max(gradient_difference, gradient_norm) <= control$gradient_tol
  list(passed = passed, per_person = per_person,
    metrics = c(MaxLogMarginalDifference = loglik_difference,
      MaxEAPDifference = eap_difference, MaxPosteriorSDDifference = sd_difference,
      MaxGradientDifference = gradient_difference, MaxGradient = gradient_norm))
}

mfrmr_optimize_accuracy_mml <- function(start, idx, config, sizes, labels,
                                      maxit, reltol, optimizer, control) {
  history <- list()
  last_par <- start
  last_candidate <- NULL
  orders <- if (is.null(control$max_points)) control$adaptive_points else
    mfrmr_mml_accuracy_orders(control$max_points)
  initial_nodes <- config$estimation_control$quad_points
  fail <- function(detail) {
    attempted <- history[[length(history)]]
    origin <- if (is.null(last_candidate)) {
      "No optimizer candidate was returned; last_parameters contains initial values."
    } else {
      retained <- last_candidate$diagnostics
      paste0("Retained candidate: stage ", retained$Stage, " (", retained$Integration,
        ", ", retained$QuadPoints, " points); accuracy was not established.")
    }
    stop(structure(list(
      message = paste0("MML numerical accuracy was not established: ", detail,
        " ", origin, " Inspect last_candidate, last_attempt, history and control; no fit was returned."),
      call = NULL, history = do.call(rbind, history), control = control,
      provenance_version = 2L, initial_parameters = start,
      last_candidate = last_candidate, last_attempt = attempted,
      parameter_source = if (is.null(last_candidate)) "initial_values" else "optimizer_candidate",
      last_parameters = last_par,
      last_integration = if (is.null(last_candidate)) NA_character_ else last_candidate$diagnostics$Integration,
      last_quad_points = if (is.null(last_candidate)) NA_integer_ else last_candidate$diagnostics$QuadPoints),
      class = c("mfrmr_mml_accuracy_error", "error", "condition")))
  }
  for (stage in seq_len(length(orders) + 1L)) {
    adaptive <- stage > 1L
    nodes <- if (adaptive) orders[stage - 1L] else initial_nodes
    check_order <- if (adaptive) nodes else orders[1L]
    config$estimation_control$mml_integration <- if (adaptive) "adaptive" else "fixed"
    config$estimation_control$quad_points <- nodes
    row <- data.frame(Stage = stage, Integration = config$estimation_control$mml_integration,
      QuadPoints = nodes, ReferencePoints = 2L * check_order + 1L,
      OptimizerConvergence = NA_integer_, LogLik = NA_real_, ParameterChange = NA_real_,
      MaxLogMarginalDifference = NA_real_, MaxEAPDifference = NA_real_,
      MaxPosteriorSDDifference = NA_real_, MaxGradientDifference = NA_real_,
      MaxGradient = NA_real_, Status = "unavailable", FailurePhase = "", Detail = "", stringsAsFactors = FALSE)
    opt <- tryCatch(run_mfrm_optimization(last_par, "MML", idx, config, sizes,
      nodes, maxit, reltol, optimizer, suppress_convergence_warning = TRUE), error = identity)
    if (inherits(opt, "error")) {
      row$FailurePhase <- "optimization"
      row$Detail <- conditionMessage(opt)
      history[[stage]] <- row
      fail(row$Detail)
    }
    row$OptimizerConvergence <- opt$convergence
    row$LogLik <- -opt$value
    row$ParameterChange <- if (stage == 1L) NA_real_ else max(abs(opt$par - last_par), 0)
    last_par <- opt$par
    checked <- tryCatch(mfrmr_check_mml_accuracy(opt, idx, config, sizes, labels,
      check_order, control), error = identity)
    if (inherits(checked, "error")) {
      row$FailurePhase <- "accuracy_check"
      row$Detail <- conditionMessage(checked)
    } else {
      row[names(checked$metrics)] <- as.list(checked$metrics)
      row$Status <- if (checked$passed) "criteria_met" else "refinement_required"
    }
    history[[stage]] <- row
    last_candidate <- list(parameters = opt$par,
      parameter_names = mfrm_checkpoint_parameter_names(sizes),
      parameter_sizes = sizes, config = config, diagnostics = row)
    if (!inherits(checked, "error") && checked$passed) {
      config$mml_accuracy <- list(status = "criteria_met", control = control,
        history = do.call(rbind, history), per_person = checked$per_person,
        scope = "integration_and_stationarity_at_retained_calibration")
      return(list(opt = opt, config = config))
    }
  }
  fail(if (is.null(control$max_points))
    "the specified adaptive orders were exhausted before all criteria were met." else
    "max_points was exhausted before all numerical criteria were met.")
}

# Only the latest automatic fit owns reusable qualification in this R session.
# The full fit and numerical settings must agree; nothing is serialized in the fit.
mfrmr_gmfrm_inference_reuse <- local({
  source <- value <- settings <- NULL
  function(fit, qualified = NULL) {
    budget <- getOption("mfrmr.max_information_bytes", 256 * 1024^2)
    current <- list(budget = budget,
      cpp11 = getOption("mfrmr.use_cpp11_backend", TRUE))
    if (is.null(qualified) && !is.null(value) &&
        identical(fit, source) && identical(current, settings)) return(value)
    source <<- value <<- settings <<- NULL
    if (!is.null(qualified) && isTRUE(qualified$check$eligible) &&
        is.numeric(budget) && !is.complex(budget) && length(budget) == 1L &&
        is.finite(budget) && budget > 0 &&
        as.double(utils::object.size(fit)) + as.double(utils::object.size(qualified)) <= budget) {
      source <<- fit
      value <<- qualified
      settings <<- current
    }
    NULL
  }
})

# Preserve each two-family engine and start policy while refining its rule.
mfrmr_gmfrm_accuracy_control <- function(control = NULL) {
  defaults <- list(max_points = 241L, loglik_tol = 1e-7, posterior_tol = 1e-6)
  if (!is.null(control)) {
    if (!is.list(control) || (length(control) && (is.null(names(control)) ||
        anyNA(names(control)) || anyDuplicated(names(control)) ||
        any(!names(control) %in% names(defaults)))))
      stop("Two-family `mml_accuracy_control` accepts max_points, loglik_tol and posterior_tol; its stationarity and information checks are unchanged.", call. = FALSE)
    defaults[names(control)] <- control
  }
  mfrmr_mml_accuracy_orders(defaults$max_points, offset = -1L)
  for (name in c("loglik_tol", "posterior_tol")) {
    value <- defaults[[name]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0)
      stop("`", name, "` must be a finite positive scalar.", call. = FALSE)
  }
  defaults
}

mfrmr_gmfrm_rule_moments <- function(fit, order, setup) {
  if (mfrmr_adaptive_integration(fit$config)) {
    return(mfrmr_make_adaptive_mml_evaluator(setup$idx, fit$config, setup$sizes,
      order, include_person_summaries = TRUE)(fit$opt$par)$per_person)
  }
  params <- expand_params(fit$opt$par, setup$sizes, fit$config)
  quad <- gauss_hermite_normal(order)
  bundle <- mfrm_mml_logprob_bundle(setup$idx, fit$config, quad, params,
    compute_base_eta(setup$idx, params, fit$config), include_probs = TRUE,
    include_linear_part = TRUE)
  posterior <- mfrm_mml_posterior_bundle(bundle)$person_bundle
  center <- as.vector(posterior$posterior %*% quad$nodes)
  data.frame(PersonIndex = posterior$person_ids, LogMarginal = posterior$log_marginal,
    EAP = center, PosteriorSD = sqrt(rowSums(posterior$posterior *
      (rep(quad$nodes, each = nrow(posterior$posterior)) - center)^2)))
}

mfrmr_gmfrm_accuracy_check <- function(fit, control) {
  inference <- mfrm_gpcm_product_inference(fit)
  failed <- inference$checks$Check[!inference$checks$Passed]
  result <- list(passed = FALSE, retry = FALSE, inference = inference,
    per_person = NULL, metrics = c(MaxLogMarginalDifference = NA_real_,
      MaxEAPDifference = NA_real_, MaxPosteriorSDDifference = NA_real_),
    detail = inference$check$review)
  # Retry only the sensitivity or moment comparisons. An earlier failed
  # check does not identify its cause or establish a boundary.
  if (length(failed) && !identical(failed, "Quadrature sensitivity")) return(result)
  setup <- do.call(mfrm_gmfrm_common_setup,
    fit$gmfrm$specification[setdiff(names(fit$gmfrm$specification), "quadrature")])
  q <- fit$config$estimation_control$quad_points
  low <- mfrmr_gmfrm_rule_moments(fit, q, setup)
  high <- mfrmr_gmfrm_rule_moments(fit, 2L*q-1L, setup)
  low <- low[order(low$PersonIndex), ]; high <- high[order(high$PersonIndex), ]
  if (!identical(low$PersonIndex, high$PersonIndex) ||
      !identical(low$PersonIndex, seq_along(fit$prep$levels$Person)))
    stop("Incomplete Person alignment in two-family integration checks.", call. = FALSE)
  difference <- abs(low[-1L] - high[-1L])
  if (any(!is.finite(as.matrix(difference))))
    stop("Non-finite two-family integration checks.", call. = FALSE)
  names(difference) <- c("LogMarginalDifference", "EAPDifference", "PosteriorSDDifference")
  result$per_person <- data.frame(Person = fit$prep$levels$Person, difference)
  result$metrics[] <- vapply(difference, max, 0)
  moments_pass <- result$metrics[1] <= control$loglik_tol &&
    max(result$metrics[-1]) <= control$posterior_tol
  result$passed <- !length(failed) && moments_pass
  result$retry <- !result$passed
  result$detail <- if (result$passed) "Integration and local information checks passed; sampling coverage remains unqualified." else
    paste(if (length(failed)) inference$check$review else "",
      if (!moments_pass) "Person likelihood or posterior moments require a finer rule." else "")
  result
}

mfrm_fit_product_slopes_auto <- function(args, supplied) {
  mfrmr_gmfrm_inference_reuse(NULL)
  control <- args$mml_accuracy_control
  orders <- mfrmr_mml_accuracy_orders(control$max_points, args$quad_points, -1L)
  resolved <- if (args$mml_engine == "em") "fixed" else "adaptive"
  history <- list(); last_candidate <- NULL
  fail <- function(detail) {
    stop(structure(list(message = paste0("MML numerical accuracy was not established: ",
      detail, " Inspect history, control and last_candidate; no accuracy-qualified fit was returned."),
      call = NULL, history = do.call(rbind, history), control = control,
      provenance_version = 2L, last_candidate = last_candidate,
      last_attempt = history[[length(history)]],
      parameter_source = if (is.null(last_candidate)) "unavailable" else "optimizer_candidate",
      last_parameters = if (is.null(last_candidate)) NULL else last_candidate$parameters,
      last_integration = if (is.null(last_candidate)) NA_character_ else resolved,
      last_quad_points = if (is.null(last_candidate)) NA_integer_ else last_candidate$diagnostics$QuadPoints),
      class = c("mfrmr_mml_accuracy_error", "error", "condition")))
  }
  stage_args <- args
  stage_args$mml_integration <- resolved
  stage_args$mml_accuracy_control <- NULL
  stage_supplied <- setdiff(supplied, "mml_accuracy_control")
  for (stage in seq_along(orders)) {
    q <- orders[stage]; stage_args$quad_points <- q
    mfrmr_fit_progress(sprintf("Automatic stage %d/%d: fitting q=%d, reference=%d",
      stage, length(orders), q, 2L*q-1L))
    row <- data.frame(Stage = stage, Integration = resolved, QuadPoints = q,
      ReferencePoints = 2L*q-1L, LogLik = NA_real_,
      MaxLogMarginalDifference = NA_real_, MaxEAPDifference = NA_real_,
      MaxPosteriorSDDifference = NA_real_, Status = "unavailable",
      FailurePhase = "", Detail = "", stringsAsFactors = FALSE)
    fit <- tryCatch(mfrm_fit_product_slopes(stage_args, stage_supplied), error = identity)
    if (inherits(fit, "error")) {
      row$FailurePhase <- "optimization"; row$Detail <- conditionMessage(fit)
      history[[stage]] <- row; fail(row$Detail)
    }
    row$LogLik <- -fit$opt$value
    mfrmr_fit_progress(sprintf("Automatic stage %d/%d: information and precision checks, q=%d/%d",
      stage, length(orders), q, 2L*q-1L))
    check_warned <- FALSE
    checked <- tryCatch(withCallingHandlers(mfrmr_gmfrm_accuracy_check(fit, control),
      warning = function(w) check_warned <<- TRUE), error = identity)
    if (inherits(checked, "error")) {
      row$FailurePhase <- "accuracy_check"; row$Detail <- conditionMessage(checked)
    } else {
      row[names(checked$metrics)] <- as.list(checked$metrics)
      row$Status <- if (checked$passed) "criteria_met" else
        if (checked$retry) "refinement_required" else "check_refused"
      row$Detail <- checked$detail
    }
    history[[stage]] <- row
    sizes <- build_param_sizes(fit$config)
    last_candidate <- list(parameters = fit$opt$par,
      parameter_names = mfrm_checkpoint_parameter_names(sizes), parameter_sizes = sizes,
      config = fit$config, diagnostics = row,
      accuracy = if (inherits(checked, "error")) NULL else list(
        checks = checked$inference$checks, numerical_checks = checked$inference$numerical_checks,
        per_person = checked$per_person))
    if (inherits(checked, "error")) fail(row$Detail)
    if (checked$passed) {
      fit$config$mml_accuracy <- list(status = "criteria_met", control = control,
        history = do.call(rbind, history), per_person = checked$per_person,
        checks = checked$inference$checks, numerical_checks = checked$inference$numerical_checks,
        scope = "two_family_integration_and_local_information; coverage_unqualified")
      fit$config$estimation_control$mml_integration_requested <- "auto"
      fit$config$estimation_control$quad_points_requested <- as.integer(args$quad_points)
      fit$config$replay_inputs$mml_integration <- "auto"
      fit$config$replay_inputs$mml_accuracy_control <- control
      fit$config$replay_inputs$quad_points <- as.integer(args$quad_points)
      fit$summary$MMLAccuracyStatus <- "criteria_met"
      if (!check_warned) mfrmr_gmfrm_inference_reuse(fit, checked$inference)
      return(fit)
    }
    if (!checked$retry) fail(row$Detail)
  }
  fail("max_points was exhausted before all numerical criteria were met.")
}

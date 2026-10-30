# Same-data posterior replicate probabilities under the fitted two-family model.
# This checks the point calibration; it neither computes nor requires covariance.
mfrm_gmfrm_response_input <- function(fit, check_solution = TRUE,
                                      operation = "response diagnostics") {
  cfg <- fit$config
  adaptive <- mfrmr_adaptive_integration(cfg)
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !mfrm_has_product_slopes(fit) || !identical(cfg$model, "GPCM") ||
      !identical(cfg$method, "MML") || isTRUE(fit$population$active) ||
      isTRUE(cfg$population_spec$active) ||
      !identical(cfg$estimation_control$mml_engine_used, if (adaptive) "direct" else "em") ||
      !identical(cfg$posterior_basis, "fixed_standard_normal")) stop(
    paste0("Two-family ", operation, " require the original fixed-N(0,1) fit with fixed-grid EM or adaptive direct MML."), call. = FALSE)
  if (!isTRUE(fit$summary$Converged) || !identical(fit$opt$convergence, 0L) ||
      (adaptive && (!identical(fit$opt$optimizer_diagnostics$ConvergenceBasis, "optimizer_gradient") ||
        !identical(fit$opt$optimizer_diagnostics$ConvergenceSeverity, "pass")))) stop(
    paste0("Resolve the two-family convergence result before ", operation, "."), call. = FALSE)
  problem <- tryCatch(do.call(mfrm_gmfrm_problem, fit$gmfrm$specification), error = function(e) NULL)
  if (is.null(problem)) stop("Recover the original two-family fitting specification.", call. = FALSE)
  ref <- problem$common
  keys <- c("facet_names", "facet_levels", "facet_specs", "facet_signs", "interaction_specs",
    "slope_facet", "step_facet", "noncenter_facet", "gpcm_spec", "n_cat")
  q <- cfg$estimation_control$quad_points
  if (!identical(lapply(fit$prep$data, identity), lapply(ref$prep$data, identity)) ||
      !identical(fit$prep$levels, ref$prep$levels) ||
      !identical(fit$prep$source_columns, ref$prep$source_columns) ||
      !identical(fit$prep$score_map, ref$prep$score_map) ||
      !identical(cfg[keys], ref$config[keys]) ||
      !is.null(cfg$weight_col) || length(cfg$dummy_facets) ||
      !identical(cfg$facet_shrinkage %||% "none", "none") ||
      !isTRUE(all.equal(problem$specification$quadrature, gauss_hermite_normal(q), tolerance = 0)) ||
      !identical(cfg$estimation_control$quadrature, problem$specification$quadrature)) stop(
    "The saved two-family data, owner roles, population and integration settings must match the fitted specification.", call. = FALSE)
  if (check_solution) {
    if (adaptive) {
      state <- mfrmr_make_adaptive_mml_evaluator(ref$idx, cfg, ref$sizes, q)(fit$opt$par)
      reltol <- fit$opt$optimizer_diagnostics$EffectiveReltol
      if (length(reltol) != 1L || !is.finite(reltol) || reltol <= 0 ||
          !identical(fit$gmfrm$controls$reltol, cfg$estimation_control$reltol)) stop(
        "The saved adaptive optimizer controls must match the fitted specification.", call. = FALSE)
      fresh <- build_optimizer_diagnostics(fit$opt, gradient = state$gradient, reltol = reltol,
        gradient_tolerance = mfrm_optimizer_gradient_tolerance(cfg, reltol))
      if (!is.finite(state$value) || any(!is.finite(state$gradient)) ||
          !identical(fresh$ConvergenceSeverity, "pass") ||
          !identical(fresh$GradientReviewTolerance, fit$opt$optimizer_diagnostics$GradientReviewTolerance) ||
          !isTRUE(all.equal(state$value, fit$opt$value, tolerance = 1e-10))) stop(
        "The retained parameters must reproduce the adaptive likelihood and optimizer gradient check.", call. = FALSE)
    } else {
      state <- problem$marginal(fit$opt$par)
      tolerance <- fit$gmfrm$controls$score_tol
      if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
          !identical(tolerance, cfg$estimation_control$em_score_tolerance) ||
          !is.finite(state$logLik) || any(!is.finite(state$gradient)) ||
          max(abs(state$gradient)) > tolerance ||
          !isTRUE(all.equal(-state$logLik, fit$opt$value, tolerance = 1e-10))) stop(
        "The retained parameters must reproduce the EM likelihood and per-Person score stopping rule.", call. = FALSE)
    }
  }
  params <- expand_params(fit$opt$par, ref$sizes, cfg)
  equal_points <- function(actual, expected, keys) {
    fields <- c(keys, "Estimate")
    is.data.frame(actual) && all(fields %in% names(actual)) &&
      isTRUE(all.equal(actual[fields], expected[fields], tolerance = 1e-10, check.attributes = FALSE))
  }
  if (!equal_points(fit$facets$others, build_other_facet_table(cfg, ref$prep, params), c("Facet", "Level")) ||
      !equal_points(fit$steps, build_step_table(cfg, ref$prep, params), c("StepFacet", "Step")) ||
      !equal_points(fit$slopes, build_slope_table(cfg, ref$prep, params), c("SlopeOwner", "SlopeFacet"))) stop(
    "Finite saved facet, step and slope estimates must agree with the fitted two-family parameters.", call. = FALSE)
  data <- problem$specification$data
  columns <- ref$prep$source_columns
  data[c(columns$person, columns$facets)] <- lapply(data[c(columns$person, columns$facets)], as.character)
  rownames(data) <- NULL
  list(data = data, columns = columns, score_levels = 0:problem$specification$max_score,
    y = data[[columns$score]], person = data[[columns$person]],
    prep = ref$prep, config = cfg, parameters = fit$opt$par, indices = ref$idx,
    expanded_parameters = params, base_eta = compute_base_eta(ref$idx, params, cfg))
}

mfrm_gmfrm_response_probabilities <- function(input, rows, order) {
  rule <- gauss_hermite_normal(order)
  indices <- which(input$person == input$person[rows[1L]])
  n <- length(indices)
  nodes <- rule$nodes; log_weights <- log(rule$weights)
  if (mfrmr_adaptive_integration(input$config)) {
    # Adapt to the complete Person record, even when only some rows are returned.
    idx <- list(person = rep(1L, n), score_k = input$y[indices], weight = rep(1, n),
      step_idx = input$indices$step_idx[indices], slope_idx = input$indices$slope_idx[indices])
    basis <- mfrmr_adaptive_quadrature_basis(idx, input$config,
      input$expanded_parameters, rule, input$base_eta[indices], person_count = 1L)
    nodes <- basis$nodes[1L, ]; log_weights <- basis$log_weights[1L, ]
  }
  expanded <- rep(indices, length(rule$nodes))
  query <- input$data[expanded, input$columns$facets, drop = FALSE]
  query$Theta <- rep(nodes, each = n)
  bundle <- mfrm_gpcm_response_evaluator(input$prep, input$config, query)$evaluate(input$parameters)
  observed <- matrix(bundle$log_probabilities[cbind(seq_along(expanded), input$y[expanded] + 1L)], n)
  logw <- log_weights + colSums(observed)
  posterior <- exp(logw - max(logw)); posterior <- posterior / sum(posterior)
  result <- t(vapply(rows, function(i) {
    at <- match(i, indices) + n * (seq_along(rule$nodes) - 1L)
    colSums(bundle$probabilities[at, , drop = FALSE] * posterior)
  }, numeric(length(input$score_levels))))
  list(probabilities = result, normalization_error = rowSums(result) - 1)
}

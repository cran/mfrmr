# A representative of each observed-law equivalence class, never cell effects.
mfrmr_comparison_quotient <- function(adjacent, sizes, rank) {
  A <- as.matrix(adjacent$design)
  scales <- sqrt(colSums(A^2)); scales[scales == 0] <- 1
  decomposition <- svd(sweep(A, 2L, scales, "/"), nu = 0L, nv = ncol(A))
  V <- decomposition$v[, seq_len(rank), drop = FALSE]
  B <- V / scales
  indices <- adjacent$map$OptimizerIndex
  variance <- build_param_slices(sizes)[["log_sigma2"]] %||% integer()
  if (!setequal(c(indices, variance), seq_len(sum(unlist(sizes)))))
    stop("Comparison quotient does not cover the native parameter coordinates.", call. = FALSE)
  transform <- matrix(0, sum(unlist(sizes)), rank + length(variance))
  transform[indices, seq_len(rank)] <- B
  if (length(variance)) transform[variance, rank + 1L] <- 1
  list(transform = transform, linear_basis = B, scaled_basis = V, scales = scales,
    linear_indices = indices, variance_index = variance, rank = rank,
    nullity = ncol(A) - rank)
}

mfrmr_comparison_boundary <- function(adjacent, idx, cfg) {
  if (cfg$method == "MML") return(list(state = "finite_mle_not_certified"))
  mfrmr_additive_response_boundary(adjacent, idx, cfg)
}

# For fixed-slope additive logits, a recession direction favors the observed
# category over every alternative in every row. Under fixed-population MML
# this also increases the integrand at every ability value. Callers establish
# the model/scale scope; this is not a free-slope or estimated-population check.
mfrmr_additive_response_boundary <- function(adjacent, idx, cfg) {
  if (!requireNamespace("lpSolve", quietly = TRUE))
    return(list(state = "not_evaluated", reason = "lpSolve is needed to check additive response recession directions."))
  contrast <- mfrmr_jml_observed_contrast_design(adjacent$design, idx$score_k,
    length(idx$score_k), cfg$n_cat - 1L)
  base <- mfrmr_jml_recession_lp_base(contrast, representation = "sparse_triplet")
  check <- mfrmr_jml_recession_target_lp(base, Matrix::colSums(contrast),
    objective_tolerance = 1e-10, certificate_tolerance = 1e-7, timeout = 0L)
  finite <- isTRUE(check$evaluated) && is.finite(check$target_capacity) &&
    identical(check$reason, "no_target_recession_direction") && check$target_capacity <= 1e-9
  check$state <- if (isTRUE(check$certified)) "recession_detected" else if (finite)
    "no_recession_detected" else "not_evaluated"
  check
}

mfrmr_fit_comparisons <- function(a, supplied) {
  allowed <- c("data", "person", "facets", "score", "rating_min", "rating_max", "weight",
    "keep_original", "missing_codes", "model", "method", "step_facet", "facet_interactions",
    "min_obs_per_interaction", "interaction_policy", "noncenter_facet", "dummy_facets",
    "positive_facets", "quad_points", "maxit", "reltol", "optimizer", "population_formula",
    "person_data", "person_id", "category_policy", "repeated_ratings", "centering_scope",
    "fit_scope", "comparison_review")
  if (a$method == "CML") allowed <- setdiff(allowed,
    c("quad_points", "population_formula", "person_data", "person_id"))
  unsupported <- setdiff(supplied, allowed)
  if (length(unsupported)) stop("Comparison fitting does not accept: ",
    paste(unsupported, collapse = ", "), ". See fit_mfrm().", call. = FALSE)
  if (!a$model %in% c("RSM", "PCM") || !a$method %in% c("JML", "MML", "CML"))
    stop("Comparison fitting supports unadjusted RSM/PCM with JML, direct fixed-grid MML or CML.", call. = FALSE)
  prep <- prepare_mfrm_data(a$data, a$person, a$facets, a$score,
    a$rating_min, a$rating_max, a$weight, a$keep_original, a$missing_codes, a$repeated_ratings)
  population <- prepare_mfrm_population_scaffold(a$data, a$person,
    population_formula = a$population_formula, person_data = a$person_data, person_id = a$person_id)
  interactions <- build_facet_interaction_specs(prep, a$facet_interactions,
    min_obs_per_interaction = a$min_obs_per_interaction, interaction_policy = a$interaction_policy)
  built <- build_estimation_config(prep, a$model, a$method, a$step_facet, NULL,
    interaction_specs = interactions, weight_col = a$weight,
    facet_signs = build_facet_signs(prep$facet_names, a$positive_facets)$signs,
    positive_facets = a$positive_facets, noncenter_facet = a$noncenter_facet,
    dummy_facets = a$dummy_facets, anchor_df = NULL, group_anchor_df = NULL,
    population = population, centering_scope = a$centering_scope)
  cfg <- built$config; sizes <- built$sizes
  review <- mfrmr_match_comparison_model(prep, cfg, a$comparison_review)
  idx <- build_indices(prep, a$step_facet, interaction_specs = cfg$interaction_specs)
  if (a$method == "CML") return(mfrmr_fit_cml_comparisons(prep, idx, cfg, sizes, review, a))
  adjacent <- mfrmr_estimability_adjacent_design(prep, idx, cfg, sizes)
  table <- review$table
  table$EstimateStatus <- "unavailable"
  table$EstimateReason <- "The declared target is unidentified or its map was not evaluated."
  table$SE <- table$Lower <- table$Upper <- NA_real_
  table$NumericalState <- "not_evaluated"
  table$QuadratureState <- if (a$method == "MML") "not_evaluated" else "not_applicable"
  table$GradientMaxAbs <- table$FinerGradientMaxAbs <- table$QuadratureChangePerPerson <- NA_real_
  category <- audit_mfrm_category_support(prep, cfg, sizes)
  table$CategoryState <- category$readiness$CategoryState
  table$BoundaryState <- "not_evaluated"
  eligible <- table$ConditionalState %in% c("target_identified_only", "all_parameters_identified")
  checks <- list(category = category, boundary = NULL, starts = NULL, numerical = NULL,
    quadrature = if (a$method == "MML") c(a$quad_points, 2L * a$quad_points + 1L) else NULL)
  quotient <- opt <- NULL
  notes <- c("Only the declared interaction contrasts are reported; individual parameter representatives are not identified cell effects.",
    "The observed likelihood is unchanged: optimization removes conditional-map null directions and preserves every invariant target.",
    "Estimates are descriptive. SE and confidence limits are unavailable; FormalInferenceEligible is false.",
    "MML checks are local numerical evidence under the declared population model, not global marginal identification or finite-MLE certificates.")
  finish <- function() {
    out <- tibble::as_tibble(table)
    attr(out, "comparison_review") <- review
    attr(out, "fit_checks") <- checks
    attr(out, "notes") <- notes
    structure(list(comparisons = out, comparison_review = review, checks = checks,
      representation = quotient, optimization = opt, prep = prep, config = cfg,
      controls = a[c("quad_points", "maxit", "reltol", "optimizer")], notes = notes),
      class = "mfrm_comparison_fit")
  }
  if (!any(eligible) || any(table$ConditionalState == "not_evaluated")) return(finish())
  if (!table$CategoryState[1L] %in% c("adequate", "weak_information")) {
    table$EstimateReason[eligible] <- "Unsupported category/step information; no comparison estimate delivered."
    return(finish())
  }
  checks$boundary <- mfrmr_comparison_boundary(adjacent, idx, cfg)
  table$BoundaryState <- checks$boundary$state
  if (a$method == "JML" && table$BoundaryState[1L] != "no_recession_detected") {
    table$EstimateReason[eligible] <- "JML recession was detected or could not be checked; no finite comparison estimate delivered."
    return(finish())
  }
  quotient <- mfrmr_comparison_quotient(adjacent, sizes, table$Rank[1L])
  T <- quotient$transform
  mml <- a$method == "MML"
  evaluator <- function(q) make_mfrm_direct_evaluator(a$method,
    make_param_cache(sizes, cfg, idx, is_mml = mml), idx, cfg, sizes,
    if (mml) gauss_hermite_normal(q) else NULL)
  raw <- evaluator(a$quad_points)
  safe <- make_mfrm_boundary_safe_objective(raw)
  value <- function(z) safe$value(drop(T %*% z))
  gradient <- function(z) drop(crossprod(T, safe$gradient(drop(T %*% z))))
  start <- build_initial_param_vector(cfg, sizes)
  z <- drop(crossprod(quotient$scaled_basis,
    quotient$scales * start[quotient$linear_indices]))
  variance <- quotient$variance_index
  if (length(variance)) z <- c(z, start[variance])
  starts <- list(z)
  if (mml && length(variance)) for (delta in c(-log(2), log(2))) {
    other <- z; other[length(other)] <- other[length(other)] + delta
    starts[[length(starts) + 1L]] <- other
  }
  method <- resolve_mfrm_optimizer(a$optimizer, ncol(T))$Used
  tolerance <- max(1e-4, 10 * a$reltol)
  run <- function(z) {
    result <- stats::optim(z, value, gradient, method = method,
      control = build_mfrm_optim_control(method, a$maxit, a$reltol))
    # A function-tolerance stop alone does not establish a stationary point.
    for (tol in mfrm_optimizer_polish_tolerances(a$reltol)) {
      if (result$convergence != 0L ||
          max(abs(raw$gradient(drop(T %*% result$par)))) <= tolerance) break
      result <- stats::optim(result$par, value, gradient, method = method,
        control = build_mfrm_optim_control(method, a$maxit, tol))
    }
    par <- drop(T %*% result$par)
    list(reduced = result, representative = par, value = raw$value(par),
      gradient = max(abs(raw$gradient(par))))
  }
  attempts <- lapply(starts, function(z) tryCatch(run(z), error = identity))
  checks$starts <- do.call(rbind, lapply(attempts, function(x) {
    if (inherits(x, "error")) return(data.frame(Objective = NA_real_, GradientMaxAbs = NA_real_,
      Convergence = NA_integer_, Error = conditionMessage(x)))
    data.frame(Objective = x$value, GradientMaxAbs = x$gradient,
      Convergence = x$reduced$convergence, Error = "")
  }))
  ok <- with(checks$starts, is.finite(Objective) & is.finite(GradientMaxAbs) &
    GradientMaxAbs <= tolerance & !is.na(Convergence) & Convergence == 0L)
  if (!all(ok)) {
    table$NumericalState <- "review"
    table$EstimateReason[eligible] <- "One or more numerical starts did not yield a qualified stationary point."
    return(finish())
  }
  opt <- attempts[[which.min(checks$starts$Objective)]]
  estimates <- vapply(attempts, function(x) drop(review$target_map %*%
    x$representative[quotient$linear_indices]), numeric(nrow(table)))
  # Only invariant targets can be compared across likelihood representatives.
  estimates <- matrix(estimates, nrow = nrow(table))
  agreement <- max(apply(estimates[eligible, , drop = FALSE], 1L, function(x) diff(range(x))))
  objective_agreement <- diff(range(checks$starts$Objective)) / cfg$n_person
  numerical <- tryCatch({
    hessian <- stats::optimHess(opt$reduced$par, value, gradient)
    eigenvalues <- eigen((hessian + t(hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values
    fine <- if (mml) evaluator(2L * a$quad_points + 1L) else raw
    c(GradientMaxAbs = opt$gradient,
      FinerGradientMaxAbs = max(abs(fine$gradient(opt$representative))),
      QuadratureChangePerPerson = abs(fine$value(opt$representative) - opt$value) / cfg$n_person,
      CurvatureRatio = min(eigenvalues) / max(1, max(eigenvalues)),
      StartTargetDifference = agreement, StartObjectiveDifferencePerPerson = objective_agreement)
  }, error = identity)
  checks$numerical <- if (inherits(numerical, "error")) conditionMessage(numerical) else numerical
  table$NumericalState <- "review"
  table$EstimateReason[eligible] <- "Curvature, numerical-start agreement or integration checks did not pass."
  if (inherits(numerical, "error") || any(!is.finite(numerical))) return(finish())
  table$GradientMaxAbs <- numerical["GradientMaxAbs"]
  if (mml) {
    table$FinerGradientMaxAbs <- numerical["FinerGradientMaxAbs"]
    table$QuadratureChangePerPerson <- numerical["QuadratureChangePerPerson"]
    table$QuadratureState <- if (numerical["FinerGradientMaxAbs"] <= tolerance &&
      numerical["QuadratureChangePerPerson"] <= 1e-6) "stable_at_retained_point" else "review"
  }
  passed <- numerical["CurvatureRatio"] > 1e-8 && agreement <= 1e-4 &&
    objective_agreement <= 1e-6 && (!mml || table$QuadratureState[1L] == "stable_at_retained_point")
  if (passed) {
    table$NumericalState <- "ready"
    table$Estimate[eligible] <- drop(review$target_map %*% opt$representative[quotient$linear_indices])[eligible]
    table$EstimateStatus[eligible] <- "descriptive_only"
    table$EstimateReason[eligible] <- if (mml)
      "Invariant conditional-model target at a qualified local marginal fit; global marginal identification and formal uncertainty remain uncertified." else
      "Invariant conditional-model target; no JML recession detected and the reduced likelihood has positive local curvature. Formal uncertainty is not supplied."
  }
  finish()
}

#' @rdname fit_mfrm
#' @param ... Unused.
#' @export
print.mfrm_comparison_fit <- function(x, ...) {
  cat("Requested joint-fit interaction comparisons\n")
  columns <- c("Comparison", "ConditionalState", "Estimate", "EstimateStatus", "BoundaryState")
  if ("IdentificationBasis" %in% names(x$comparisons)) columns <- c(columns, "IdentificationBasis")
  if (x$config$method == "CML") columns <- c(columns, "SE", "Lower", "Upper")
  print(x$comparisons[columns])
  cat(if (x$config$method == "CML") "Identified contrasts: conditional-information SE and pointwise 95% Wald limits where available. Inspect EstimateReason for withheld uncertainty.\n" else
    "Descriptive estimates only; no standard errors or confidence limits.\n")
  invisible(x)
}

#' @rdname fit_mfrm
#' @export
summary.mfrm_comparison_fit <- function(object, ...) object$comparisons

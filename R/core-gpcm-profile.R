# Reuse the two-family EM cell aggregation for repeated marginal evaluations.
# This is the observed marginal score, not a frozen-posterior M-step objective.
# The EM kernel averages by Person; profiling requires the total likelihood.
mfrm_gmfrm_profile_evaluator <- function(problem, fallback) {
  cached_par <- NULL; state <- NULL
  evaluate <- function(par) {
    if (!identical(par, cached_par)) {
      cached_par <<- NULL
      candidate <- problem$marginal(par)
      state <<- list(value = candidate$value * problem$n_person,
        gradient = candidate$gradient * problem$n_person)
      # Preserve the original evaluator's typed boundary errors and rejection
      # behavior at invalid trials. Never replace a raw error by a finite value.
      if (any(!is.finite(c(state$value, state$gradient)))) state <<-
        list(value = fallback$value(par), gradient = fallback$gradient(par))
      # optim() may reuse its callback buffer; own the cache key, and publish
      # it only after both evaluations succeed.
      cached_par <<- par + 0
    }
    state
  }
  list(value = function(par) evaluate(par)$value,
    gradient = function(par) evaluate(par)$gradient,
    optimization_scale = problem$n_person)
}

# One linear log-slope constraint; all remaining calibration and population
# coordinates are nuisance parameters. No change to the likelihood or bounds.
mfrm_gpcm_profile_constraint <- function(start, contrast, value) {
  pivot <- which.max(abs(contrast))
  free <- setdiff(seq_along(start), pivot)
  basis <- diag(length(start))[, free, drop = FALSE]
  basis[pivot, ] <- -contrast[free] / contrast[pivot]
  offset <- rep(0, length(start)); offset[pivot] <- value / contrast[pivot]
  list(free = free, basis = basis,
    embed = function(z) offset + drop(basis %*% z))
}

mfrm_gpcm_profile_nuisance <- function(evaluator, reference, start, contrast,
                                       value, maxit) {
  map <- mfrm_gpcm_profile_constraint(start, contrast, value)
  safe <- make_mfrm_boundary_safe_objective(evaluator)
  fn <- function(z) safe$value(map$embed(z))
  gr <- function(z) drop(crossprod(map$basis, safe$gradient(map$embed(z))))
  raw_fn <- function(z) evaluator$value(map$embed(z))
  raw_gr <- function(z) drop(crossprod(map$basis, evaluator$gradient(map$embed(z))))
  stages <- list(); selected <- NULL; initial <- start[map$free]
  for (tol in c(1e-10, 1e-13)) {
    result <- tryCatch({
      # An invalid initial value is not a finite optimizer penalty.
      raw_fn(initial)
      opt <- optim(initial, fn, gr, method = "BFGS",
        control = c(build_mfrm_optim_control("BFGS", maxit, tol),
          # Match the EM search scale without rescaling LR values or checks.
          list(fnscale = evaluator$optimization_scale %||% 1)))
      list(par = opt$par, value = raw_fn(opt$par), gradient = raw_gr(opt$par),
        code = opt$convergence, error = "")
    }, error = function(e) list(error = conditionMessage(e)))
    stages[[length(stages) + 1L]] <- result
    if (!is.null(result$value) && all(is.finite(c(result$value, result$gradient))) &&
        (is.null(selected) || result$value < selected$value ||
          (result$value <= selected$value + 1e-9 &&
            max(abs(result$gradient)) < max(abs(selected$gradient))))) selected <- result
    if (!is.null(selected)) initial <- selected$par
    if (!is.null(selected) && selected$code == 0L && max(abs(selected$gradient)) <= 1e-4) break
  }
  if (!is.null(evaluator$optimization_scale) && !is.null(selected) && selected$code != 0L) {
    # A small nuisance slope can make location/step coordinates poorly scaled.
    # Retry the same constrained objective in fixed curvature coordinates;
    # the constraint, likelihood and final native-gradient checks are unchanged.
    scaled <- tryCatch({
      scale <- mfrm_optimizer_curvature_scale(selected$par, raw_fn, raw_gr)
      if (is.null(scale$transform)) stop(scale$error)
      center <- selected$par
      native <- function(z) center + drop(scale$transform %*% z)
      opt <- optim(rep(0, length(center)),
        function(z) fn(native(z)),
        function(z) drop(crossprod(scale$transform, gr(native(z)))), method = "BFGS",
        # The transform already accounts for the total-likelihood curvature.
        control = build_mfrm_optim_control("BFGS", maxit, 1e-13))
      par <- native(opt$par)
      list(par = par, value = raw_fn(par), gradient = raw_gr(par), code = opt$convergence,
        error = "", stage = "curvature_scaled", curvature = scale)
    }, error = function(e) list(error = conditionMessage(e), stage = "curvature_scaled"))
    stages[[length(stages) + 1L]] <- scaled
    if (!is.null(scaled$value) && all(is.finite(c(scaled$value, scaled$gradient))) &&
        (scaled$value < selected$value || (scaled$value <= selected$value + 1e-9 &&
          max(abs(scaled$gradient)) < max(abs(selected$gradient))))) selected <- scaled
  }
  reference_gradient <- if (!is.null(selected)) tryCatch(
    drop(crossprod(map$basis, reference$gradient(map$embed(selected$par)))),
    error = function(e) Inf) else Inf
  needs_polish <- !is.null(selected) &&
    (max(abs(selected$gradient)) > 1e-4 || any(!is.finite(reference_gradient)) ||
      max(abs(reference_gradient)) > 1e-4)
  if (needs_polish && selected$code == 0L) {
    # Reuse the fitter's one-step curvature polish after BFGS stops too early.
    # Check the higher-order score too: a passing low-order score can narrowly
    # miss the reference tolerance. Acceptance thresholds remain unchanged.
    # The helper accepts only a smaller gradient with no objective deterioration.
    proposal <- mfrm_optimizer_curvature_proposal(selected$par, raw_fn, raw_gr)
    polish <- tryCatch({
      if (is.null(proposal$par)) stop(proposal$error)
      list(par = proposal$par, value = raw_fn(proposal$par), gradient = raw_gr(proposal$par),
        code = selected$code, error = "", stage = "curvature_polish")
    }, error = function(e) list(error = conditionMessage(e), stage = "curvature_polish"))
    stages[[length(stages) + 1L]] <- polish
    if (!is.null(polish$value)) selected <- polish
  }
  row <- data.frame(NLL = NA_real_, MaxGradient = NA_real_, OptimizerCode = NA_integer_,
    NLLChange = NA_real_, GradientChange = NA_real_, ReferenceGradient = NA_real_,
    NuisanceMinEigenvalue = NA_real_, ReferenceMinEigenvalue = NA_real_,
    Passed = FALSE, Detail = paste(vapply(stages, `[[`, "", "error"), collapse = " "))
  parameters <- NULL
  if (!is.null(selected)) {
    parameters <- map$embed(selected$par)
    row$NLL <- selected$value; row$MaxGradient <- max(abs(selected$gradient))
    row$OptimizerCode <- selected$code
    checked <- tryCatch({
      g <- drop(crossprod(map$basis, reference$gradient(parameters)))
      c(abs(reference$value(parameters) - selected$value),
        max(abs(g - selected$gradient)), max(abs(g)))
    }, error = function(e) { row$Detail <<- conditionMessage(e); rep(NA_real_, 3L) })
    row[c("NLLChange", "GradientChange", "ReferenceGradient")] <- as.list(checked)
    row$Passed <- isTRUE(row$OptimizerCode == 0 && row$MaxGradient <= 1e-4 &&
      all(is.finite(checked)) && checked[1] <= 1e-5 && all(checked[-1] <= 1e-4))
    if (row$Passed) {
      # Two starts can coincide at a stationary nuisance maximum or saddle.
      # Review the constrained NLL Hessian at both integration orders; do not
      # floor eigenvalues or turn an unresolved solution into a profile point.
      curvature <- tryCatch({
        minimum <- function(objective) {
          h <- stats::optimHess(selected$par,
            function(z) objective$value(map$embed(z)),
            function(z) drop(crossprod(map$basis, objective$gradient(map$embed(z)))))
          min(eigen(h, symmetric = TRUE, only.values = TRUE)$values)
        }
        low <- minimum(evaluator)
        c(low, if (identical(evaluator, reference)) low else minimum(reference))
      }, error = function(e) {
        row$Detail <<- paste("Constrained nuisance curvature could not be verified:", conditionMessage(e))
        rep(NA_real_, 2L)
      })
      row[c("NuisanceMinEigenvalue", "ReferenceMinEigenvalue")] <- as.list(curvature)
      row$Passed <- all(is.finite(curvature)) && all(curvature > 0)
      if (!row$Passed && !nzchar(trimws(row$Detail))) row$Detail <-
        "Constrained nuisance curvature is not positive definite at both integration orders."
    }
    if (!row$Passed && !nzchar(trimws(row$Detail))) row$Detail <-
      "Constrained optimization or higher-order integration checks failed."
  }
  list(check = row, parameters = parameters, stages = stages)
}

mfrm_gpcm_profile_search <- function(evaluator, reference, start, contrast,
                                     level, maxit = 400L, max_steps = 8L,
                                     initial_step = .5) {
  center <- sum(contrast * start)
  baseline <- evaluator$value(start)
  cutoff <- stats::qchisq(level, 1)
  records <- list(); attempts <- list(); invalid_baseline <- FALSE
  evaluate <- function(value) {
    key <- sprintf("%.17g", value)
    if (!is.null(records[[key]])) return(records[[key]])
    results <- lapply(list(retained = start, neutral = rep(0, length(start))), function(p)
      mfrm_gpcm_profile_nuisance(evaluator, reference, p, contrast, value, maxit))
    checks <- do.call(rbind, lapply(results, `[[`, "check"))
    finite <- which(is.finite(checks$NLL))
    best <- if (length(finite)) finite[which.min(checks$NLL[finite])] else NA_integer_
    nll <- if (is.na(best)) NA_real_ else checks$NLL[best]
    spread <- if (length(finite) == 2L) diff(range(checks$NLL)) else NA_real_
    passed <- all(checks$Passed) && is.finite(spread) && spread <= 1e-5
    detail <- paste(unique(trimws(checks$Detail[nzchar(trimws(checks$Detail))])), collapse = " ")
    if (all(checks$Passed) && !passed) detail <- "The two constrained starts disagree."
    if (is.finite(nll) && nll < baseline - 1e-5) {
      invalid_baseline <<- TRUE; passed <- FALSE
      detail <- "A constrained search improved on the source likelihood; refit the source."
    }
    row <- data.frame(LogSlope = value, Slope = exp(value), NLL = nll,
      LR = 2 * (nll - baseline), StartDifference = spread, Passed = passed, Detail = detail)
    records[[key]] <<- row; attempts[[key]] <<- results
    row
  }
  origin <- evaluate(center)
  endpoint <- function(direction) {
    if (!origin$Passed) return(list(value = NA_real_, status = "source_profile_failed"))
    inner <- center; previous <- origin$LR
    failed_edge <- NULL
    for (i in seq_len(max_steps)) {
      # A failed far trial can lie beyond an accurately computable endpoint.
      # Contract toward the passing inner point, retaining every failed trial.
      # Never jump across an unresolved point or change integration settings.
      outer <- if (is.null(failed_edge)) center + direction * initial_step * 1.5^(i - 1L) else
        (inner + failed_edge) / 2
      row <- evaluate(outer)
      if (!row$Passed) {
        if (invalid_baseline) return(list(value = NA_real_, status = "source_likelihood_improved"))
        failed_edge <- outer
        next
      }
      if (row$LR < previous - 1e-4) return(list(value = NA_real_, status = "nonmonotone_profile"))
      if (row$LR >= cutoff) {
        root <- tryCatch(stats::uniroot(function(v) {
          z <- evaluate(v)
          if (!z$Passed) stop("Profile point failed.")
          z$LR - cutoff
        }, sort(c(inner, outer)), tol = 1e-6, maxiter = 100L)$root,
          error = function(e) NA_real_)
        if (!is.finite(root)) return(list(value = NA_real_, status = "endpoint_search_failed"))
        z <- evaluate(root)
        if (!z$Passed || abs(z$LR - cutoff) > 1e-4 || !is.finite(exp(root)) || exp(root) <= 0)
          return(list(value = NA_real_, status = "endpoint_check_failed"))
        return(list(value = exp(root), status = "computed"))
      }
      inner <- outer; previous <- row$LR
    }
    list(value = NA_real_, status = if (is.null(failed_edge)) "crossing_not_found" else "numerical_or_start_failure")
  }
  limits <- lapply(c(-1, 1), endpoint)
  endpoints <- data.frame(Side = c("Lower", "Upper"),
    Bound = vapply(limits, `[[`, 0, "value"), Status = vapply(limits, `[[`, "", "status"))
  if (invalid_baseline) {
    endpoints$Bound <- NA_real_; endpoints$Status <- "source_likelihood_improved"
  }
  profile <- do.call(rbind, records); rownames(profile) <- NULL
  profile <- profile[order(profile$LogSlope), , drop = FALSE]
  for (i in seq_len(2L)) {
    side <- profile[profile$Passed & (profile$LogSlope - center) * c(-1, 1)[i] >= 0, ]
    side <- side[order(abs(side$LogSlope - center)), ]
    if (any(diff(side$LR) < -1e-4) && !invalid_baseline) {
      endpoints$Bound[i] <- NA_real_; endpoints$Status[i] <- "nonmonotone_profile"
    }
  }
  # Finite scans cannot exclude unsampled/disconnected likelihood regions.
  list(profile = profile, attempts = attempts, endpoints = endpoints,
    baseline = baseline, cutoff = cutoff)
}

mfrm_gpcm_profile_interval <- function(fit, slope, level, control) {
  product <- mfrm_has_product_slopes(fit)
  adaptive <- mfrmr_adaptive_integration(fit$config)
  labels <- as.character(fit$config$gpcm_spec$levels)
  if (product) {
    owner <- names(slope)
    if (!is.character(slope) || length(slope) != 1L || is.na(slope) ||
        length(owner) != 1L || is.na(owner) || !owner %in% fit$config$slope_facet ||
        !slope %in% fit$config$gpcm_spec$component_levels[[owner]]) stop(
      "For two-family method = 'profile', identify one owner and level, e.g. slope = c(Task = 't1'); use your fitted column and level names.", call. = FALSE)
  } else if (!is.character(slope) || length(slope) != 1L || is.na(slope) || !slope %in% labels)
    stop("For method = 'profile', supply one exact slope-level name in `slope`.", call. = FALSE)
  defaults <- list(maxit = 400L, max_steps = 8L, initial_step = .5)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyDuplicated(names(control)) || any(!names(control) %in% names(defaults)))))
    stop("`profile_control` accepts only maxit, max_steps and initial_step.", call. = FALSE)
  controls <- utils::modifyList(defaults, control, keep.null = TRUE)
  for (name in names(controls)) {
    v <- controls[[name]]
    if (!is.numeric(v) || is.complex(v) || length(v) != 1L || !is.finite(v) || v <= 0 ||
        (name != "initial_step" && (v != floor(v) || v > .Machine$integer.max)))
      stop("Profile controls must be positive finite numbers; maxit and max_steps must be integers.", call. = FALSE)
  }
  population <- fit$config$population_spec
  if (!product && (length(fit$config$interaction_specs) || nrow(mfrmr_calibration_extract_anchors(fit)) ||
      !isTRUE(population$active) || !identical(population$design_columns, "(Intercept)")))
    stop("Slope profiling currently requires an estimated intercept-only normal population and no anchors or interactions.", call. = FALSE)
  info <- if (product) mfrm_gpcm_product_inference(fit) else mfrm_gpcm_inference(fit)
  if (!isTRUE(info$check$eligible)) stop(info$check$review, call. = FALSE)
  if (!is.null(info$check$caution)) warning(info$check$caution, call. = FALSE)
  config <- fit$config; sizes <- build_param_sizes(config)
  idx <- build_indices(fit$prep, config$step_facet, config$slope_facet, config$interaction_specs,
    gpcm_spec = config$gpcm_spec)
  order <- config$estimation_control$quad_points
  higher <- max(order + 10L, 2L * order - 1L)
  if (inherits(tryCatch(gauss_hermite_normal(higher), error = identity), "error")) higher <- order + 10L
  context <- function(q) {
    cfg <- config; cfg$estimation_control$quad_points <- q
    direct <- make_mfrm_direct_evaluator("MML", make_param_cache(sizes, cfg, idx, is_mml = TRUE),
      idx, cfg, sizes, gauss_hermite_normal(q))
    if (!product) return(direct)
    if (adaptive) {
      # Keep moving nodes, scales and Jacobians in every constrained evaluation.
      # Only the optimizer is scaled per Person; LR values stay on the total scale.
      direct$optimization_scale <- length(fit$prep$levels$Person)
      return(direct)
    }
    specification <- fit$gmfrm$specification
    specification$quadrature <- gauss_hermite_normal(q)
    mfrm_gmfrm_profile_evaluator(do.call(mfrm_gmfrm_problem, specification), direct)
  }
  low <- context(order); high <- context(higher); p <- fit$opt$par
  # Two-family stationarity is already checked in SE units by product_inference.
  # Profiling additionally checks absolute likelihood/gradient integration errors;
  # every constrained solve retains the same stricter nuisance-score checks.
  if (abs(low$value(p) - high$value(p)) > 1e-5 ||
      max(abs(low$gradient(p) - high$gradient(p))) > 1e-4 ||
      (!product && max(abs(high$gradient(p))) > 1e-4))
    stop("Source integration checks failed; refit at a finer quadrature before profiling.", call. = FALSE)
  target <- if (product) mfrm_gpcm_product_slope_targets(p, config)$components else
    mfrm_gpcm_slope_target(fit, "relative", NULL, "ratio")
  index <- if (product) which(target$metadata$Owner == owner & target$metadata$Level == slope) else match(slope, labels)
  contrast <- target$jacobian[index, ]
  started <- proc.time()[["elapsed"]]
  result <- do.call(mfrm_gpcm_profile_search, c(list(evaluator = low, reference = high,
    start = p, contrast = contrast, level = level), controls))
  result$elapsed <- proc.time()[["elapsed"]] - started
  se <- sqrt(drop(contrast %*% info$covariance %*% contrast))
  estimate <- target$estimate[index]
  eligible <- all(result$endpoints$Status == "computed")
  descriptions <- c(computed = "computed", source_profile_failed = "profile checks at the fitted slope failed",
    numerical_or_start_failure = "numerical checks or agreement between starts failed",
    nonmonotone_profile = "the sampled likelihood was not monotone away from the estimate",
    endpoint_search_failed = "an endpoint search evaluation failed",
    endpoint_check_failed = "endpoint accuracy was insufficient",
    crossing_not_found = "no cutoff crossing in the searched range",
    source_likelihood_improved = "a better source likelihood was found")
  result$endpoints$Detail <- unname(descriptions[result$endpoints$Status])
  review <- if (eligible) paste("Connected local profile-LR interval; chi-square(1) approximation.",
    "Two-start and integration checks do not prove global optimality or finite-sample coverage.") else
    paste("Profile interval unresolved:", paste(paste(result$endpoints$Side, result$endpoints$Detail), collapse = "; "),
      "Search limits do not establish an unbounded interval.")
  label <- if (product) paste(owner, slope, sep = " = ") else unname(slope)
  tab <- data.frame(SlopeFacet = label, Estimate = estimate, SE = NA_real_, LogSE = NA_real_,
    CI_Lower = result$endpoints$Bound[1], CI_Upper = result$endpoints$Bound[2],
    CI_Level = level, CIEligible = eligible, SEEligible = FALSE,
    CIUse = if (eligible) "approximate_pointwise_profile" else "unavailable", InferenceReview = review)
  description <- paste("Relative GPCM slope:", fit$config$slope_facet, slope, "(geometric mean one)")
  if (product) {
    tab$CoverageQualified <- FALSE
    tab$SlopeOwner <- owner; tab$SlopeLevel <- unname(slope)
    tab$ScaleReference <- target$metadata$ScaleReference[index]
    description <- paste("Two-family component slope:", owner, slope,
      if (tab$ScaleReference == "geometric_mean_one") "(first-family geometric mean one)" else "(free second-family slope on fixed N(0,1))")
  }
  out <- mfrm_gpcm_interval_result(tab, level, "Profile likelihood; chi-square(1) approximation", description)
  attr(out, "profile") <- result
  attr(out, "wald") <- data.frame(SlopeFacet = label, Estimate = estimate,
    Lower = exp(log(estimate) - stats::qnorm((1 - level) / 2, lower.tail = FALSE) * se),
    Upper = exp(log(estimate) + stats::qnorm((1 - level) / 2, lower.tail = FALSE) * se))
  attr(out, "source") <- mfrm_gpcm_inference_source(fit)
  attr(out, "settings") <- c(list(scale = if (product) "standardized" else "relative", method = "profile", slope = slope,
    quad_points = order, comparison_quad_points = higher,
    integration = if (adaptive) "adaptive" else "fixed",
    engine = config$estimation_control$mml_engine_used), controls)
  if (product) {
    attr(out, "checks") <- info$checks
    attr(out, "numerical_checks") <- info$numerical_checks
    attr(out, "score_rank") <- info$score_rank
    attr(out, "scale_note") <- "The first slope family has geometric mean one; second-family slopes are free on fixed N(0,1). No automatic rater-quality interpretation."
  }
  attr(out, "cautions") <- unique(c(info$check$caution,
    "A local connected profile is not a global likelihood or coverage certificate.", if (!eligible) review))
  if (!eligible) warning(review, call. = FALSE)
  out
}

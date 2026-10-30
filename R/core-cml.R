# Generalized ESF coefficients and conditional sufficient-statistic moments.
# Log coefficients avoid overflow; the covariance recursion avoids subtracting
# two large raw second moments.
mfrmr_cml_esf <- function(features, par, information = FALSE) {
  p <- length(par)
  logz <- 0
  means <- matrix(0, 1L, p)
  covs <- if (information) array(0, c(p, p, 1L)) else NULL
  for (x in features) {
    K <- nrow(x) - 1L
    eta <- drop(x %*% par)
    if (any(!is.finite(eta))) stop("Nonfinite CML category predictors.", call. = FALSE)
    nz <- length(logz) + K
    next_logz <- numeric(nz)
    next_means <- matrix(0, nz, p)
    next_covs <- if (information) array(0, c(p, p, nz)) else NULL
    for (s in seq_len(nz) - 1L) {
      k <- seq.int(max(0L, s - length(logz) + 1L), min(K, s))
      from <- s - k + 1L
      logs <- logz[from] + eta[k + 1L]
      peak <- max(logs)
      weights <- exp(logs - peak)
      next_logz[s + 1L] <- peak + log(sum(weights))
      weights <- weights / sum(weights)
      candidates <- means[from, , drop = FALSE] + x[k + 1L, , drop = FALSE]
      mu <- drop(crossprod(weights, candidates))
      next_means[s + 1L, ] <- mu
      if (information) {
        centered <- sweep(candidates, 2L, mu, "-") * sqrt(weights)
        V <- crossprod(centered)
        for (h in seq_along(from)) V <- V + weights[h] * covs[, , from[h]]
        next_covs[, , s + 1L] <- V
      }
    }
    logz <- next_logz; means <- next_means; covs <- next_covs
  }
  list(logz = logz, means = means, covariance = covs)
}

mfrmr_cml_problem <- function(prep, idx, cfg, sizes) {
  adjacent <- mfrmr_estimability_adjacent_design(prep, idx, cfg, sizes,
    include_person = FALSE, include_population_beta = FALSE)
  p <- sum(unlist(sizes)); K <- cfg$n_cat - 1L; n <- nrow(prep$data)
  if (!p || !setequal(adjacent$map$OptimizerIndex, seq_len(p)))
    stop("CML requires free structural coordinates covered by the native model map.", call. = FALSE)
  order <- order(adjacent$map$OptimizerIndex)
  A <- adjacent$design[, order, drop = FALSE]
  keys <- do.call(paste, c(idx$facets, sep = ":"))
  cells <- unique(keys)
  features <- lapply(match(cells, keys), function(i) {
    increments <- as.matrix(A[i + (seq_len(K) - 1L) * n, , drop = FALSE])
    out <- matrix(0, K + 1L, p)
    for (k in seq_len(K)) out[k + 1L, ] <- out[k, ] + increments[k, ]
    out
  })
  cell <- match(keys, cells)
  persons <- split(seq_len(n), idx$person)
  groups <- list(); group_index <- new.env(hash = TRUE, parent = emptyenv())
  observed <- numeric(p)
  totals <- vapply(persons, function(at) sum(idx$score_k[at]), numeric(1))
  maximum <- lengths(persons) * K
  extreme <- totals == 0 | totals == maximum
  for (i in which(!extreme)) {
    at <- persons[[i]]
    roster <- sort(cell[at])
    key <- paste(roster, collapse = ":")
    g <- group_index[[key]]
    if (is.null(g)) {
      g <- length(groups) + 1L
      group_index[[key]] <- g
      groups[[g]] <- list(features = features[roster],
        counts = numeric(length(at) * K + 1L))
    }
    groups[[g]]$counts[totals[i] + 1L] <- groups[[g]]$counts[totals[i] + 1L] + 1L
    for (j in at) observed <- observed + features[[cell[j]]][idx$score_k[j] + 1L, ]
  }
  list(groups = groups, observed = observed, parameters = p,
    map = adjacent$map[order, , drop = FALSE],
    person_totals = data.frame(Person = prep$levels$Person[as.integer(names(persons))], Total = totals,
      Maximum = maximum, Extreme = extreme, row.names = NULL))
}

mfrmr_cml_evaluate <- function(problem, par, information = FALSE) {
  value <- -sum(problem$observed * par)
  gradient <- -problem$observed
  hessian <- if (information) matrix(0, length(par), length(par)) else NULL
  for (g in problem$groups) {
    esf <- mfrmr_cml_esf(g$features, par, information)
    at <- which(g$counts > 0)
    value <- value + sum(g$counts[at] * esf$logz[at])
    gradient <- gradient + drop(crossprod(g$counts[at], esf$means[at, , drop = FALSE]))
    if (information) for (s in at)
      hessian <- hessian + g$counts[s] * esf$covariance[, , s]
  }
  list(value = value, gradient = gradient, information = hessian)
}

# Max-sum recursion independently checks the support certificate from the LP.
mfrmr_cml_support <- function(features, direction) {
  values <- 0
  for (x in features) {
    eta <- drop(x %*% direction)
    next_values <- rep(-Inf, length(values) + nrow(x) - 1L)
    for (k in seq_len(nrow(x))) {
      at <- seq_along(values) + k - 1L
      next_values[at] <- pmax(next_values[at], values + eta[k])
    }
    values <- next_values
  }
  values
}

# A DAG represents every pattern at the observed totals without enumerating
# their exponential number. Potentials bound each path's linear statistic.
mfrmr_cml_boundary <- function(problem, score_zero) {
  if (!requireNamespace("lpSolve", quietly = TRUE))
    stop("Native CML needs the suggested package 'lpSolve' for its finite-MLE check.", call. = FALSE)
  p <- problem$parameters
  rows <- list(); endpoint <- numeric(0); node_count <- 0L
  add <- function(columns, values) {
    keep <- values != 0
    rows[[length(rows) + 1L]] <<- list(columns = columns[keep], values = values[keep])
  }
  for (g in problem$groups) {
    totals <- which(g$counts > 0) - 1L
    m <- length(g$features); K <- nrow(g$features[[1L]]) - 1L
    nodes <- matrix(NA_integer_, m + 1L, max(totals) + 1L)
    nodes[1L, 1L] <- 0L
    for (j in seq_len(m)) for (s in 0:min(j * K, max(totals))) {
      if (!any(totals >= s & totals <= s + (m - j) * K)) next
      node_count <- node_count + 1L
      nodes[j + 1L, s + 1L] <- node_count
      for (k in 0:min(K, s)) {
        previous <- nodes[j, s - k + 1L]
        if (is.na(previous)) next
        x <- g$features[[j]][k + 1L, ]
        # Positive indices are direction coordinates; negative ones are nodes.
        columns <- c(seq_len(p), -node_count, if (previous > 0) -previous)
        add(columns, c(-x, 1, if (previous > 0) -1))
      }
    }
    at <- nodes[m + 1L, totals + 1L]
    endpoint <- c(endpoint, stats::setNames(g$counts[totals + 1L], at))
  }
  # Sum of maxima cannot exceed the observed statistic in a recession direction.
  add(c(seq_len(p), -as.integer(names(endpoint))), c(problem$observed, -endpoint))
  q <- p + node_count
  triplets <- lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    columns <- ifelse(r$columns > 0L, r$columns, p - r$columns)
    cbind(i, c(columns, columns + q), c(r$values, -r$values))
  })
  triplets[[length(triplets) + 1L]] <- cbind(length(rows) + 1L,
    c(seq_len(p), seq_len(p) + q), 1)
  objective <- c(-score_zero, rep(0, node_count), score_zero, rep(0, node_count))
  lp <- lpSolve::lp(direction = "max", objective.in = objective,
    dense.const = do.call(rbind, triplets),
    const.dir = c(rep(">=", length(rows)), "<="),
    const.rhs = c(rep(0, length(rows)), 1), timeout = 0L)
  if (lp$status != 0L || !is.finite(lp$objval))
    return(list(state = "not_evaluated", solver_status = lp$status))
  direction <- lp$solution[seq_len(p)] - lp$solution[seq_len(p) + q]
  gap <- -sum(problem$observed * direction)
  for (g in problem$groups) {
    at <- which(g$counts > 0)
    gap <- gap + sum(g$counts[at] * mfrmr_cml_support(g$features, direction)[at])
  }
  gain <- -sum(score_zero * direction)
  state <- if (lp$objval <= 1e-9) "no_recession_detected" else if (
    gain > 1e-7 && abs(gap) <= 1e-7 && sum(abs(direction)) <= 1 + 1e-7)
      "recession_detected" else "not_evaluated"
  list(state = state, solver_status = lp$status, objective = lp$objval,
    support_gap = gap, gain = gain, direction = direction)
}

mfrmr_cml_optimize <- function(problem, initial, a) {
  p <- problem$parameters
  saved_par <- saved_eval <- NULL
  evaluate <- function(par) {
    if (!identical(par, saved_par)) {
      saved_eval <<- mfrmr_cml_evaluate(problem, par)
      saved_par <<- par + 0
    }
    saved_eval
  }
  method <- resolve_mfrm_optimizer(a$optimizer, p)$Used
  value <- function(par) evaluate(par)$value
  gradient <- function(par) evaluate(par)$gradient
  opt <- stats::optim(initial, value, gradient,
    method = method, control = build_mfrm_optim_control(method, a$maxit, a$reltol))
  tolerance <- max(1e-6, min(1e-4, 10 * a$reltol))
  evaluate_solution <- function(par) tryCatch(mfrmr_cml_evaluate(problem, par, TRUE),
    error = function(e) {
      # An information failure must not discard a separately evaluable score.
      out <- mfrmr_cml_evaluate(problem, par)
      out$information_error <- conditionMessage(e)
      out
    })
  solution <- evaluate_solution(opt$par)
  # Function-tolerance stops can leave a visible score. Analytic Newton
  # refinement uses the same conditional likelihood and full information.
  for (iteration in seq_len(3L)) {
    if (opt$convergence != 0L || max(abs(solution$gradient)) <= tolerance) break
    step <- tryCatch(solve(solution$information, solution$gradient), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) break
    fraction <- 1
    repeat {
      candidate <- opt$par - fraction * step
      next_value <- value(candidate)
      if (next_value <= solution$value + 1e-12 * max(1, abs(solution$value)) || fraction < 1e-6) break
      fraction <- fraction / 2
    }
    if (fraction < 1e-6) break
    opt$par <- candidate
    solution <- evaluate_solution(opt$par)
  }
  opt$value <- solution$value
  gradient_max <- max(abs(solution$gradient))
  if (opt$convergence != 0L || any(!is.finite(opt$par)) || !is.finite(opt$value) ||
      length(solution$gradient) != p || !is.finite(gradient_max) || gradient_max > tolerance)
    stop("CML optimization did not reach a qualified stationary point (maximum absolute score = ",
      signif(gradient_max, 4), "); no Wald intervals returned.", call. = FALSE)
  H <- solution$information
  covariance <- matrix(NA_real_, p, p)
  ev <- rep(NA_real_, p); inverse_error <- NA_real_
  information_reason <- tryCatch({
    if (!is.null(solution$information_error)) stop(solution$information_error, call. = FALSE)
    if (!is.matrix(H) || !all(dim(H) == p) ||
        any(!is.finite(H)) || any(diag(H) <= 0))
      stop("CML conditional information requires finite entries and positive diagonal values.", call. = FALSE)
    scaled <- H / outer(sqrt(diag(H)), sqrt(diag(H)))
    ev <- eigen(scaled, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(ev)) || min(ev) <= 1e-8 * max(ev))
      stop("CML conditional information is singular or too ill-conditioned for unregularized Wald intervals.", call. = FALSE)
    candidate <- chol2inv(chol(H))
    inverse_error <- max(abs(H %*% candidate - diag(p)))
    if (!is.finite(inverse_error) || inverse_error > 1e-7 || any(!is.finite(candidate)))
      stop("CML information inversion failed its residual check.", call. = FALSE)
    covariance <- candidate
    ""
  }, error = conditionMessage)
  if (nzchar(information_reason)) warning("CML point estimates retained; uncertainty unavailable: ",
    information_reason, call. = FALSE)
  if (!is.matrix(H) || !all(dim(H) == p)) H <- matrix(NA_real_, p, p)
  list(opt = opt, covariance = covariance, information = H,
    checks = list(gradient_max = gradient_max, information_eigenvalues = ev, inverse_error = inverse_error,
      information_status = if (nzchar(information_reason)) "unavailable" else "available",
      information_reason = information_reason))
}

mfrmr_fit_cml <- function(a, supplied) {
  allowed <- c("data", "person", "facets", "score", "rating_min", "rating_max",
    "weight", "keep_original", "missing_codes", "model", "method", "step_facet",
    "facet_interactions", "min_obs_per_interaction", "interaction_policy",
    "positive_facets", "maxit", "reltol", "optimizer", "category_policy",
    "repeated_ratings", "centering_scope", "fit_scope")
  unsupported <- setdiff(supplied, allowed)
  if (length(unsupported)) stop("Native CML does not accept: ",
    paste(unsupported, collapse = ", "), ". See fit_mfrm().", call. = FALSE)
  if (!a$model %in% c("RSM", "PCM") || a$fit_scope != "full")
    stop("Native CML currently requires a full-rank RSM/PCM model and fit_scope = 'full'.", call. = FALSE)
  prep <- prepare_mfrm_data(a$data, a$person, a$facets, a$score,
    a$rating_min, a$rating_max, a$weight, a$keep_original, a$missing_codes, a$repeated_ratings,
    allow_single_person = TRUE)
  if (any(prep$data$Weight != 1) || !isTRUE(prep$keep_original))
    stop("Native CML requires unit weights and category_policy = 'preserve'.", call. = FALSE)
  if ("duplicate_person_facet_cells" %in% prep$preparation_notes$Condition)
    stop("CML intervals require resolution of repeated rating cells; correct duplicates or explicitly declare repeated_ratings = 'independent' when justified.", call. = FALSE)
  if (a$model == "RSM" && !is.null(a$step_facet))
    stop("RSM CML uses shared steps; omit step_facet.", call. = FALSE)
  interactions <- build_facet_interaction_specs(prep, a$facet_interactions,
    min_obs_per_interaction = a$min_obs_per_interaction, interaction_policy = a$interaction_policy)
  built <- build_estimation_config(prep, a$model, "CML", a$step_facet, NULL,
    interaction_specs = interactions, weight_col = a$weight,
    facet_signs = build_facet_signs(prep$facet_names, a$positive_facets)$signs,
    positive_facets = a$positive_facets, noncenter_facet = "Person",
    dummy_facets = NULL, anchor_df = NULL, group_anchor_df = NULL,
    centering_scope = a$centering_scope)
  cfg <- built$config; sizes <- built$sizes
  idx <- build_indices(prep, a$step_facet, interaction_specs = cfg$interaction_specs)
  problem <- mfrmr_cml_problem(prep, idx, cfg, sizes)
  p <- problem$parameters
  zero <- mfrmr_cml_evaluate(problem, numeric(p), TRUE)
  sd <- sqrt(pmax(diag(zero$information), 0))
  # Constant conditional statistics can leave round-off variance; scaling
  # such a column to unit variance would manufacture identification.
  sd[sd <= 1e-6 * max(sd)] <- 1
  eigen_zero <- eigen(zero$information / outer(sd, sd), symmetric = TRUE, only.values = TRUE)$values
  ranks <- vapply(c(1e-12, 1e-10, 1e-8), function(tol)
    sum(eigen_zero > tol * max(1, max(eigen_zero))), integer(1))
  if (any(ranks < p)) stop("The conditioned CML model is unidentified or numerically rank-sensitive (rank ",
    min(ranks), " of ", p, "). Person-constant effects disappear on conditioning; revise the design/model.", call. = FALSE)
  boundary <- mfrmr_cml_boundary(problem, zero$gradient)
  if (boundary$state != "no_recession_detected") stop("CML finite-MLE check: ",
    boundary$state, ". No finite estimate or Wald interval is returned.", call. = FALSE)
  fitted <- mfrmr_cml_optimize(problem, build_initial_param_vector(cfg, sizes), a)
  opt <- fitted$opt; H <- fitted$information; covariance <- fitted$covariance
  names(opt$par) <- problem$map$Coordinate
  dimnames(covariance) <- dimnames(H) <- list(names(opt$par), names(opt$par))
  params <- expand_params(opt$par, sizes, cfg)
  structure(list(opt = opt, prep = prep, config = cfg,
    facets = list(others = build_other_facet_table(cfg, prep, params)),
    steps = build_step_table(cfg, prep, params),
    interactions = list(effects = build_interaction_effect_table(cfg, prep, params)),
    covariance = covariance, information = H, person_totals = problem$person_totals,
    checks = c(list(conditional_ranks = ranks, boundary = boundary), fitted$checks),
    controls = a[c("maxit", "reltol", "optimizer")],
    notes = c("CML conditions on each Person total and observed rating roster; it estimates no Person ability or population distribution.",
      mfrm_single_person_caution(prep),
      "Extreme-total Persons contribute zero conditional information and remain recorded in person_totals.",
      "Wald intervals are asymptotic under the correctly specified response model, independent Persons, fixed structural dimension and information growing in every identified direction; finite-sample coverage is not guaranteed.",
      "More Persons alone does not establish the information-growth condition. Response-dependent assignment, stopping or missingness requires separate justification that the conditional response law is preserved.",
      if (nzchar(fitted$checks$information_reason)) paste("Point estimates retained; uncertainty unavailable:",
        fitted$checks$information_reason))),
    class = "mfrm_cml_fit")
}

mfrmr_cml_wald <- function(estimate, jacobian, covariance, level, unavailable_reason = NULL) {
  if (any(!is.finite(estimate)) || any(!is.finite(jacobian)))
    stop("CML targets and derivatives must be finite; rescale the contrast coefficients.", call. = FALSE)
  unavailable <- length(unavailable_reason) == 1L && !is.na(unavailable_reason) && nzchar(unavailable_reason)
  fixed <- rowSums(abs(jacobian)) == 0
  V <- if (unavailable) matrix(NA_real_, length(estimate), length(estimate),
    dimnames = list(rownames(jacobian), rownames(jacobian))) else
    symmetrize_matrix(jacobian %*% covariance %*% t(jacobian))
  V[fixed, ] <- 0; V[, fixed] <- 0
  if (!unavailable && any(!is.finite(V)))
    stop("CML contrast covariance must be finite; rescale the contrast coefficients.", call. = FALSE)
  if (!unavailable && any(diag(V)[!fixed] <= 0))
    stop("Non-fixed CML targets require positive representable variance; rescale the contrast coefficients.", call. = FALSE)
  se <- sqrt(diag(V))
  margin <- stats::qnorm((1 - level) / 2, lower.tail = FALSE) * se
  margin[fixed] <- NA_real_
  if (!unavailable && (any(!is.finite((estimate - margin)[!fixed])) || any(!is.finite((estimate + margin)[!fixed]))))
    stop("CML interval bounds must be finite; rescale the contrast coefficients.", call. = FALSE)
  list(table = data.frame(Estimate = as.numeric(estimate), SE = se,
    Lower = estimate - margin, Upper = estimate + margin,
    Status = ifelse(fixed, "fixed", if (unavailable) "information_unavailable" else "available"),
    CIEligible = !fixed & !unavailable,
    InferenceReview = ifelse(fixed,
      "Constraint-fixed target; SE is zero by construction and no inferential interval is supplied.",
      if (unavailable) paste("Point estimate retained; conditional uncertainty is unavailable:", unavailable_reason) else
        "Pointwise conditional-information normal approximation; numerical availability does not guarantee finite-sample coverage."),
    row.names = NULL), covariance = V)
}

mfrmr_cml_facet_intervals <- function(fit, facet, contrasts, method, clusters, adjust, level, scale) {
  if (match.arg(scale, c("native", "standardized")) != "native")
    stop("CML does not estimate a population SD; standardized contrasts require an MML fit that estimates that SD.", call. = FALSE)
  if (match.arg(method, c("model", "sandwich")) != "model" || !is.null(clusters) || !identical(adjust, FALSE))
    stop("CML facet intervals currently require method = 'model', clusters = NULL and adjust = FALSE.", call. = FALSE)
  if (!is.numeric(level) || is.complex(level) || is.object(level) || !is.null(dim(level)) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("Supply 0 < level < 1.", call. = FALSE)
  level <- unname(level)
  target <- mfrm_facet_contrasts(fit, facet, contrasts)
  sizes <- build_param_sizes(fit$config)
  params <- expand_params(fit$opt$par, sizes, fit$config)
  J <- matrix(0, nrow(target$jacobian), length(fit$opt$par),
    dimnames = list(rownames(target$contrasts), names(fit$opt$par)))
  J[, build_param_slices(sizes)[[facet]]] <- target$jacobian
  wald <- mfrmr_cml_wald(drop(target$contrasts %*% params$facets[[facet]]), J, fit$covariance, level,
    fit$checks$information_reason)
  tab <- cbind(Target = rownames(target$contrasts), wald$table)
  tab$ModelSE <- tab$SE; tab$ModelLower <- tab$Lower; tab$ModelUpper <- tab$Upper
  if (fit$prep$n_person == 1L) {
    tab$CoverageQualified <- FALSE
    tab$InferenceCaution <- mfrm_single_person_caution(fit$prep)
  }
  structure(list(table = tab, covariance = wald$covariance, model_covariance = wald$covariance,
    parameter_covariance = fit$covariance, model_parameter_covariance = fit$covariance,
    contrasts = target$contrasts, target_jacobian = J,
    settings = list(facet = facet, method = "model", level = level, adjust = FALSE,
      adjustment_factor = 1, persons = fit$prep$n_person, clusters = NA_integer_,
      cluster_score_rank = NA_integer_, free_parameters = length(fit$opt$par),
      reference = "pointwise normal", quad_points = NA_integer_,
      procedure = "cml_conditional_information", scale = "native_ability"),
    fit = fit, cautions = fit$notes), class = "mfrm_facet_intervals")
}

#' @rdname fit_mfrm
#' @param object,x A fitted CML result for its summary, print, covariance and
#'   log-likelihood methods, or a comparison-only result for its summary and
#'   print methods.
#' @param ... Unused by the CML methods.
#' @export
summary.mfrm_cml_fit <- function(object, ...) {
  locations <- do.call(rbind, lapply(object$config$facet_names, function(facet) {
    tab <- mfrm_facet_intervals(object, facet)$table
    cbind(Facet = facet, tab)
  }))
  sizes <- build_param_sizes(object$config)
  slices <- build_param_slices(sizes)
  J <- matrix(0, nrow(object$steps), length(object$opt$par))
  J[, slices$steps] <- as.matrix(mfrmr_step_jacobian_sparse(object$config, sizes)$jacobian)
  steps <- cbind(object$steps[setdiff(names(object$steps), "Estimate")],
    mfrmr_cml_wald(object$steps$Estimate, J, object$covariance, .95, object$checks$information_reason)$table)
  interactions <- object$interactions$effects
  if (nrow(interactions)) {
    offset <- 0L; row <- 0L
    J <- matrix(0, nrow(interactions), length(object$opt$par))
    for (spec in object$config$interaction_specs) {
      local <- as.matrix(mfrmr_interaction_jacobian_sparse(spec)$jacobian)
      J[row + seq_len(nrow(local)), slices[["interactions"]][offset + seq_len(ncol(local))]] <- local
      row <- row + nrow(local); offset <- offset + ncol(local)
    }
    interactions <- cbind(interactions[setdiff(names(interactions), "Estimate")],
      mfrmr_cml_wald(interactions$Estimate, J, object$covariance, .95, object$checks$information_reason)$table)
  }
  list(method = "CML", model = object$config$model, logLik = -object$opt$value,
    locations = locations, steps = steps, interactions = interactions,
    person_totals = object$person_totals, checks = object$checks, notes = object$notes)
}

#' @rdname fit_mfrm
#' @export
print.mfrm_cml_fit <- function(x, ...) {
  cat(x$config$model, "conditional maximum likelihood\n")
  cat(nrow(x$person_totals), "Persons;", sum(x$person_totals$Extreme), "extreme totals;",
    length(x$opt$par), "structural parameters\n")
  print(summary(x)$locations, row.names = FALSE)
  invisible(x)
}

#' @rdname mfrm_fit_extractors
#' @export
vcov.mfrm_cml_fit <- function(object, ...) {
  rlang::check_dots_empty()
  value <- object$covariance
  if (identical(object$checks$information_status, "unavailable")) {
    attr(value, "unavailable_reason") <- object$checks$information_reason
    attr(value, "checks") <- object$checks
    attr(value, "cautions") <- object$notes
  }
  if (object$prep$n_person == 1L) attr(value, "cautions") <- object$notes
  value
}

#' @rdname mfrm_fit_extractors
#' @export
logLik.mfrm_cml_fit <- function(object, ...) {
  rlang::check_dots_empty()
  structure(-object$opt$value, df = length(object$opt$par),
    nobs = nobs(object), informative_persons = sum(!object$person_totals$Extreme),
    parameter_count = length(object$opt$par), likelihood_basis = "conditional",
    ICReason = "Conditional likelihood only; its criteria must not be compared with marginal likelihood criteria.",
    class = "logLik")
}

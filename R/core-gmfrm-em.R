# Internal numerical implementation of the Uto--Ueno (2020), equation 9,
# two-slope model. The public route excludes anchors, weights, estimated
# population and general inferential outputs; component-slope intervals have
# separate experimental checks. The numerical caller supplies a
# fixed quadrature rule for N(0, 1). Muraki (1992) motivates MML--EM; sirt's
# rm.facets is an algorithm reference, not the same response model.
# slope_facets is ordered: the first supplies the GM1 slope reference and
# centered locations; the second supplies free slopes/locations and steps.
# Data column order and domain-specific labels do not determine those roles.
mfrm_gmfrm_problem <- function(data, max_score, quadrature,
                               slope_facets = c("Task", "Rater"),
                               person = "Person", score = "Score") {
  if (!is.character(slope_facets) || length(slope_facets) != 2L ||
      !is.character(person) || length(person) != 1L ||
      !is.character(score) || length(score) != 1L) {
    stop("Supply two slope-facet column names and one person and score column name.")
  }
  slope_facets <- unname(slope_facets); person <- unname(person); score <- unname(score)
  required <- c(person, slope_facets, score)
  reserved <- c("Person", "Score", "Weight", "score_k", "Theta", "InputRow",
    "SlopeOwner1", "SlopeOwner2", "Slope1", "Slope2", "EffectiveSlope",
    "ObservedContext", "ExpectedScore", "ScoreVariance", "ResponseSensitivity", "Information")
  if (anyNA(required) || any(!nzchar(trimws(required))) || anyDuplicated(required) ||
      any(slope_facets %in% reserved)) {
    stop("Model columns must be distinct nonempty names; slope facets cannot use reserved data or response columns.")
  }
  first <- slope_facets[1L]; second <- slope_facets[2L]
  if (!is.data.frame(data) || anyDuplicated(names(data)) || !all(required %in% names(data)) ||
      !nrow(data) || anyNA(data[required])) {
    stop("Supply observed rows for the declared person, slope facets and score without missing values.")
  }
  if (!is.numeric(max_score) || length(max_score) != 1L || !is.finite(max_score) ||
      max_score < 1 || max_score > .Machine$integer.max || max_score != floor(max_score) ||
      !is.numeric(data[[score]]) || is.complex(data[[score]]) || any(!is.finite(data[[score]])) ||
      any(data[[score]] != as.integer(data[[score]])) ||
      any(data[[score]] < 0 | data[[score]] > max_score)) {
    stop("Scores must be integers between zero and max_score.")
  }
  if (anyDuplicated(data[c(person, slope_facets)])) {
    stop("Repeated person-facet combinations require a separate observation model.")
  }
  nodes <- quadrature$nodes
  weights <- quadrature$weights
  if (!is.numeric(nodes) || !is.numeric(weights) || length(nodes) < 2L ||
      length(nodes) != length(weights) || any(!is.finite(c(nodes, weights))) ||
      any(weights <= 0) || abs(sum(weights) - 1) > 1e-10) {
    stop("Supply finite quadrature nodes and positive normalized weights.")
  }
  levels <- lapply(data[c(person, slope_facets)], function(x) sort(unique(as.character(x))))
  ids <- Map(function(x, lev) match(as.character(x), lev), data[names(levels)], levels)
  ni <- length(levels[[first]])
  nr <- length(levels[[second]])
  np <- length(levels[[person]])
  if (ni < 2L || nr < 2L) stop("This two-slope implementation requires at least two levels per slope facet.")
  contrast <- function(n) rbind(diag(n - 1L), rep(-1, n - 1L))
  ci <- contrast(ni)
  ct <- contrast(max_score)
  cells <- unique(as.data.frame(unname(ids[slope_facets]), check.names = FALSE))
  names(cells) <- slope_facets
  cell <- match(paste(ids[[first]], ids[[second]]), paste(cells[[first]], cells[[second]]))
  common <- mfrm_gmfrm_common_setup(data, max_score, slope_facets, person, score)
  slope_cells <- common$config$gpcm_spec$cells
  common_cell <- match(paste(cells[[first]], cells[[second]]),
                       paste(slope_cells[[first]], slope_cells[[second]]))
  nc <- nrow(cells)
  nq <- length(nodes)
  # The first facet has centered locations and geometric-mean-one slopes.
  # The second has free locations/slopes and owns the centered step sets.
  Z <- cbind(ci[cells[[first]], , drop = FALSE], diag(nr)[cells[[second]], , drop = FALSE])
  if (qr(Z)$rank < ncol(Z)) stop("Facet crossings do not identify both facet blocks.")
  nb <- ncol(Z) + nr * (max_score - 1L)
  D <- lapply(0:max_score, function(k) {
    step <- matrix(0, nc, nr * (max_score - 1L))
    if (k > 0L && max_score > 1L) {
      cum <- colSums(ct[seq_len(k), , drop = FALSE])
      for (r in seq_len(nr)) {
        cols <- (r - 1L) * (max_score - 1L) + seq_len(max_score - 1L)
        step[cells[[second]] == r, cols] <- rep(cum, each = sum(cells[[second]] == r))
      }
    }
    cbind(k * Z, step)
  })
  parameter_names <- c(paste0("location:", first, ":", levels[[first]][-ni]),
    paste0("location:", second, ":", levels[[second]]),
    if (max_score > 1L) unlist(lapply(levels[[second]], function(r)
      paste0("step:", second, ":", r, ":", seq_len(max_score - 1L)))),
    paste0("log_slope:", first, ":", levels[[first]][-ni]),
    paste0("log_slope:", second, ":", levels[[second]]))
  start <- setNames(numeric(nb + ncol(Z)), parameter_names)
  expanded_cell <- rep(seq_len(nc), nq)
  expanded_node <- rep(nodes, each = nc)
  Zq <- Z[expanded_cell, , drop = FALSE]
  Dq <- lapply(D, function(d) d[expanded_cell, , drop = FALSE])
  score_index <- data[[score]] + 1L
  obs_node <- outer(cell, nc * (seq_len(nq) - 1L), "+")
  quad_basis <- list(log_weights = matrix(log(weights), np, nq, byrow = TRUE))

  cached_par <- NULL; cached_kernel <- NULL
  kernel <- function(par) {
    if (length(par) != length(start) || any(!is.finite(par))) stop("Invalid parameter vector.")
    if (identical(par, cached_par)) return(cached_kernel)
    cached_par <<- NULL; cached_kernel <<- NULL
    params <- tryCatch(expand_params(par, common$sizes, common$config),
      mfrmr_gpcm_slope_numeric_boundary_error = function(e) NULL)
    if (is.null(params)) return(NULL)
    cell_slope <- params$slopes[common_cell]
    location <- unname(params$facets[[first]][cells[[first]]] + params$facets[[second]][cells[[second]]])
    step_cum <- unname(t(apply(params$steps_mat, 1L, function(s) c(0, cumsum(s)))))
    # This existing conditional GPCM kernel is shared with JML, prediction,
    # and diagnostics. Using it here does not change the MML objective to JML.
    bundle <- mfrm_jml_probability_bundle(
      eta = expanded_node - location[expanded_cell], score_k = rep(0L, nc * nq),
      model = "GPCM", step_cum = step_cum, criterion_idx = cells[[second]][expanded_cell],
      slopes = cell_slope, slope_idx = expanded_cell, include_log_probs = TRUE)
    slope <- cell_slope[expanded_cell]
    logits <- bundle$linear_part * slope
    if (any(!is.finite(logits))) return(NULL)
    if (any(!is.finite(bundle$log_probs))) return(NULL)
    cached_kernel <<- list(log_prob = bundle$log_probs, prob = bundle$probs,
      logits = logits, slope = slope)
    # Own the key: optim may reuse its callback parameter buffer.
    cached_par <<- par + 0
    cached_kernel
  }
  expected_counts <- function(posterior) {
    counts <- array(0, c(nc, nq, max_score + 1L))
    for (k in seq_len(max_score + 1L)) {
      selected <- which(score_index == k)
      if (length(selected)) {
        total <- rowsum(posterior[ids[[person]][selected], , drop = FALSE], cell[selected])
        counts[as.integer(rownames(total)), , k] <- total
      }
    }
    matrix(counts, nc * nq, max_score + 1L)
  }
  complete <- function(par, counts, derivatives = TRUE) {
    ker <- kernel(par)
    if (is.null(ker)) return(list(value = Inf, gradient = rep(NA_real_, length(par))))
    value <- -sum(counts * ker$log_prob) / np
    if (!derivatives) return(list(value = value))
    residual <- rowSums(counts) * ker$prob - counts
    gb <- Reduce(`+`, lapply(seq_along(Dq), function(k)
      -drop(crossprod(Dq[[k]], ker$slope * residual[, k]))))
    # A category-common offset has zero derivative contribution because the
    # residuals sum to zero. Use normalized log probabilities so rounding in
    # that sum is not multiplied by a large shared logit at a near tie.
    gg <- drop(crossprod(Zq, rowSums(residual * ker$log_prob)))
    list(value = value, gradient = c(gb, gg) / np)
  }
  marginal <- function(par, derivatives = TRUE) {
    ker <- kernel(par)
    if (is.null(ker)) return(list(value = Inf, gradient = rep(NA_real_, length(par))))
    log_observed <- matrix(ker$log_prob[cbind(as.vector(obs_node), rep(score_index, nq))], nrow(data), nq)
    # Center the Person likelihood before applying the quadrature weights, and
    # normalize in those shifted coordinates. Subtracting a restored marginal
    # log likelihood can otherwise erase the prior or inflate E-step counts.
    person_bundle <- mfrm_mml_person_bundle(log_observed, ids[[person]],
      quad_basis, include_posterior = TRUE)
    log_marginal <- person_bundle$log_marginal
    # The common helper keeps first-observed order; counts use sorted level IDs.
    posterior <- person_bundle$posterior[order(person_bundle$person_ids), , drop = FALSE]
    dimnames(posterior) <- list(levels[[person]], paste0("node", seq_len(nq)))
    counts <- expected_counts(posterior)
    list(value = -sum(log_marginal) / np, logLik = sum(log_marginal),
      posterior = posterior, counts = counts,
      gradient = if (derivatives) complete(par, counts)$gradient else NULL)
  }
  unpack <- function(par) {
    if (length(par) != length(start) || any(!is.finite(par))) stop("Invalid parameter vector.")
    params <- expand_params(par, common$sizes, common$config)
    dimnames(params$steps_mat) <- list(levels[[second]], seq_len(max_score))
    list(locations = setNames(lapply(slope_facets, function(owner)
        setNames(params$facets[[owner]], levels[[owner]])), slope_facets),
      slopes = setNames(lapply(slope_facets, function(owner)
        setNames(params$slope_components[[owner]]$Slope, levels[[owner]])), slope_facets),
      steps = params$steps_mat, step_facet = second)
  }
  response <- function(par, newdata) {
    evaluator <- mfrm_gpcm_response_evaluator(common$prep, common$config, newdata)
    bundle <- evaluator$evaluate(par)
    slope <- bundle$effective_slope
    variance <- bundle$variance
    summary <- data.frame(InputRow = seq_len(nrow(newdata)), Theta = newdata$Theta,
      SlopeOwner1 = first, SlopeOwner2 = second,
      Slope1 = bundle$slope_components[[first]], Slope2 = bundle$slope_components[[second]],
      EffectiveSlope = slope, ObservedContext = evaluator$contexts$ObservedContext,
      ExpectedScore = drop(bundle$probabilities %*% (0:max_score)), ScoreVariance = variance,
      ResponseSensitivity = slope * variance, Information = slope^2 * variance)
    for (owner in slope_facets) summary[[owner]] <- as.character(newdata[[owner]])
    list(summary = summary, probabilities = bundle$probabilities, log_probabilities = bundle$log_probabilities)
  }
  list(start = start, kernel = kernel, complete = complete, marginal = marginal,
    unpack = unpack, response = response, n_person = np, cells = cells, levels = levels,
    nodes = nodes, weights = weights, common = common,
    # Plain inputs reconstruct the same problem after saveRDS(), without
    # retaining closures/caches or silently choosing a new integration rule.
    specification = list(data = as.data.frame(data[required]), max_score = max_score,
      quadrature = list(nodes = nodes, weights = weights),
      slope_facets = slope_facets, person = person, score = score))
}

# Numerical M steps make this a generalized EM: every accepted step must
# increase Q and the same fixed-quadrature marginal likelihood. Convergence
# requires a small marginal score, not just a small change in likelihood.
mfrm_gmfrm_em <- function(problem, start = problem$start, maxit = 500L,
                          mstep_maxit = 100L, score_tol = 1e-6) {
  if (length(maxit) != 1L || !is.finite(maxit) || maxit < 1L || maxit != as.integer(maxit) ||
      length(mstep_maxit) != 1L || !is.finite(mstep_maxit) || mstep_maxit < 1L ||
      mstep_maxit != as.integer(mstep_maxit) || length(score_tol) != 1L ||
      !is.finite(score_tol) || score_tol <= 0) stop("Invalid EM controls.")
  par <- start
  state <- problem$marginal(par)
  if (!is.finite(state$value)) stop("Nonfinite initial likelihood.")
  trace <- list(data.frame(iteration = 0L, logLik = state$logLik,
    max_score = max(abs(state$gradient)), q_gain = NA_real_, mstep_code = NA_integer_, halvings = 0L,
    mstep_score = NA_real_, mstep_polished = FALSE, mstep_polish_error = ""))
  reason <- "iteration_limit"
  for (iteration in seq_len(maxit)) {
    if (max(abs(state$gradient)) <= score_tol) { reason <- "score_tolerance"; break }
    mfrmr_fit_progress(sprintf("Two-family EM: q=%d, iteration %d/%d",
      length(problem$nodes), iteration, maxit))
    counts <- state$counts
    old_q <- problem$complete(par, counts, FALSE)$value
    mstep <- stats::optim(par,
      function(p) problem$complete(p, counts, FALSE)$value,
      function(p) problem$complete(p, counts)$gradient, method = "BFGS",
      control = list(maxit = mstep_maxit, reltol = 1e-12))
    mstep_score <- max(abs(problem$complete(mstep$par, counts)$gradient))
    polished <- FALSE; polish_error <- ""
    if (mstep$convergence == 0L && is.finite(mstep_score) && mstep_score > score_tol) {
      # A small change in Q is not the requested score accuracy. Refine Q at
      # the SAME E-step counts, so this remains a generalized EM M step.
      # The shared helper requires positive usable curvature and an improved
      # gradient without objective deterioration. The ascent checks below
      # still apply to both Q and the observed marginal likelihood.
      proposal <- mfrm_optimizer_curvature_proposal(mstep$par,
        function(p) problem$complete(p, counts, FALSE)$value,
        function(p) problem$complete(p, counts)$gradient)
      polish_error <- proposal$error
      if (!is.null(proposal$par)) {
        mstep$par <- proposal$par
        mstep_score <- max(abs(problem$complete(mstep$par, counts)$gradient))
        polished <- TRUE
      }
    }
    accepted <- FALSE
    for (halvings in 0:25) {
      candidate <- par + (mstep$par - par) * 2^(-halvings)
      q <- problem$complete(candidate, counts, FALSE)$value
      proposed <- problem$marginal(candidate)
      # Tolerance permits floating-point roundoff only, per person.
      if (is.finite(q) && is.finite(proposed$value) && q <= old_q + 1e-12 &&
          proposed$value <= state$value + 1e-12) { accepted <- TRUE; break }
    }
    if (!accepted) { reason <- "no_ascent_step"; break }
    par <- candidate
    state <- proposed
    trace[[length(trace) + 1L]] <- data.frame(iteration = iteration, logLik = state$logLik,
      max_score = max(abs(state$gradient)), q_gain = old_q - q,
      mstep_code = mstep$convergence, halvings = halvings,
      mstep_score = mstep_score, mstep_polished = polished, mstep_polish_error = polish_error)
  }
  converged <- max(abs(state$gradient)) <= score_tol
  if (converged) reason <- "score_tolerance"
  common <- problem$common
  slopes <- build_slope_table(common$config, common$prep,
    expand_params(par, common$sizes, common$config))
  slope_readiness <- mfrmr_readiness_gpcm_slope_parameters(common$config, slopes)
  slopes <- apply_mfrm_slope_readiness(slopes, list(parameters = slope_readiness))
  list(slopes = slopes, slope_readiness = slope_readiness, par = par, logLik = state$logLik, posterior = state$posterior,
    slope_overview = mfrm_fit_slope_overview(common$config, slopes),
    scale_contract = mfrm_fit_scale_contract(list(config = common$config)),
    specification = problem$specification,
    controls = list(start = start, maxit = maxit, mstep_maxit = mstep_maxit,
      score_tol = score_tol),
    converged = converged, reason = reason, trace = do.call(rbind, trace),
    parameters = problem$unpack(par), max_score = max(abs(state$gradient)),
    parameter_map = problem$common$parameter_map,
    nonlinear_transformation = mfrmr_nonlinear_transformation_audit(
      par, problem$common$sizes, problem$common$config, "log_slopes"))
}

# Compare adaptive solutions, never fixed-grid and adaptive objective values.
# EM or supplied native vectors are starts only; inference keeps its own gates.
mfrm_gmfrm_adaptive_fit <- function(problem, config, quad_points, maxit, reltol,
                                    optimizer, policy) {
  capture <- function(call) {
    value <- NULL; error <- ""; warnings <- character()
    elapsed <- system.time(withCallingHandlers(tryCatch(value <- call(),
      error = function(e) error <<- conditionMessage(e)), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }))[["elapsed"]]
    list(value = value, error = error, warnings = unique(warnings), seconds = elapsed)
  }
  supplied <- is.list(policy)
  starts <- if (supplied) policy else list(neutral = problem$start)
  if (supplied) {
    if (is.object(starts) || length(starts) < 2L || is.null(names(starts)) ||
        anyNA(names(starts)) || any(!nzchar(trimws(names(starts)))) || anyDuplicated(names(starts)))
      stop("Supply a plain named list of at least two native starting vectors with distinct nonempty names.", call. = FALSE)
    for (start in starts) if (!is.numeric(start) || is.complex(start) ||
        !is.null(dim(start)) || length(start) != length(problem$start) || any(!is.finite(start)) ||
        (!is.null(names(start)) && !identical(names(start), names(problem$start))))
      stop("Each native start must be finite and match the model's free-parameter length and order; supplied parameter names must match exactly.", call. = FALSE)
    policy <- "explicit"
  }
  em <- NULL
  if (identical(policy, "neutral_em")) {
    mfrmr_fit_progress(sprintf("Adaptive start preparation: EM, q=%d", quad_points))
    em <- capture(function() mfrm_gmfrm_em(problem, maxit = maxit, score_tol = 1e-6))
    starts["em"] <- list(em$value$par)
    if (!is.null(em$value)) em$value <- em$value[c("par", "trace", "controls",
      "converged", "reason", "max_score")]
  }
  common <- problem$common
  evaluate <- mfrmr_make_adaptive_mml_evaluator(common$idx, config, common$sizes, quad_points)
  attempts <- lapply(names(starts), function(name) {
    mfrmr_fit_progress(sprintf("Adaptive start %s: q=%d", name, quad_points))
    start <- starts[[name]]; initial <- NA_real_
    trial <- capture(function() {
      if (is.null(start)) stop(paste("EM initialization unavailable:", em$error), call. = FALSE)
      initial <<- evaluate(start)$value
      if (!is.finite(initial)) stop("Nonfinite adaptive starting objective.", call. = FALSE)
      run_mfrm_optimization(start, "MML", common$idx, config, common$sizes,
        quad_points, maxit, reltol, optimizer)
    })
    c(list(start = start, initial_nll = initial), trial)
  })
  names(attempts) <- names(starts)
  table <- do.call(rbind, lapply(names(attempts), function(name) {
    x <- attempts[[name]]; opt <- x$value
    data.frame(Start = name, InitialNLL = x$initial_nll, FinalNLL = opt$value %||% NA_real_,
      ConvergenceCode = opt$convergence %||% NA_integer_,
      GradientCheck = opt$optimizer_diagnostics$ConvergenceSeverity %||% "unavailable",
      OptimizerStages = nrow(opt$optimizer_polish$Stages) %||% 0L,
      ElapsedSeconds = x$seconds, Error = x$error,
      Warnings = paste(x$warnings, collapse = " | "), Selected = FALSE)
  }))
  finite <- which(is.finite(table$FinalNLL))
  record <- list(policy = policy, selected = NA_character_, table = table,
    controls = list(quad_points = quad_points, maxit = maxit, reltol = reltol, optimizer = optimizer),
    attempts = attempts, em = em, total_seconds = sum(table$ElapsedSeconds) + (em$seconds %||% 0))
  if (!length(finite)) stop(structure(list(
    message = "No adaptive starting value produced a finite solution; inspect the initialization record.",
    call = NULL, initialization = record), class = c("mfrmr_gmfrm_start_error", "error", "condition")))
  # Prefer the lowest finite objective, even if it is not converged. Do not hide
  # a better unfinished solution by returning a worse one with a success flag.
  selected <- finite[which.min(table$FinalNLL[finite])]
  record$table$Selected[selected] <- TRUE
  record$selected <- table$Start[selected]
  opt <- attempts[[selected]]$value
  initial <- table$InitialNLL[is.finite(table$InitialNLL)]
  if (length(initial) && opt$value > min(initial) +
      64 * .Machine$double.eps * max(1, abs(opt$value), abs(initial))) {
    opt$convergence <- 2L
    opt$optimizer_diagnostics$ConvergenceCode <- 2L
    opt$optimizer_diagnostics$ConvergenceStatus <- "starting_value_review"
    opt$optimizer_diagnostics$ConvergenceSeverity <- "review"
    opt$optimizer_diagnostics$ConvergenceReason <- "better_starting_objective"
    opt$optimizer_diagnostics$ConvergenceDetail <-
      "A retained starting value has a lower adaptive objective than the selected terminal point."
    warning(opt$optimizer_diagnostics$ConvergenceDetail, call. = FALSE)
  }
  opt$mml_initialization <- record
  for (message in attempts[[selected]]$warnings) warning(message, call. = FALSE)
  if (any(nzchar(table$Error))) warning(
    "Some adaptive starting values failed; inspect fit$opt$mml_initialization.", call. = FALSE)
  opt
}

# Assemble an actual EM result with the shared fit tables/readiness contract.
# Public dispatch validates its input scope before using this adapter.
mfrm_gmfrm_fit_result <- function(problem, result) {
  if (!identical(problem$specification, result$specification)) {
    stop("The EM result and fitting specification must match.", call. = FALSE)
  }
  quad <- problem$specification$quadrature
  expected_quad <- gauss_hermite_normal(length(quad$nodes))
  # Eigensolver roundoff differs across platforms. Compare nodes at their
  # own scale and weights relatively, including the tiny tail weights. Keep
  # the supplied rule for every likelihood, score and retained result below.
  if (length(quad$weights) != length(expected_quad$weights) ||
      is.complex(quad$nodes) || is.complex(quad$weights) ||
      !isTRUE(all(abs(quad$nodes - expected_quad$nodes) <=
                  1e-12 * pmax(1, abs(expected_quad$nodes)))) ||
      !isTRUE(all(abs(quad$weights / expected_quad$weights - 1) <= 1e-10))) {
    stop("The shared fit result currently requires standard-normal Gauss-Hermite quadrature.", call. = FALSE)
  }
  common <- problem$common
  prep <- common$prep; config <- common$config; sizes <- common$sizes
  state <- problem$marginal(result$par)
  if (!is.finite(state$logLik) || !isTRUE(all.equal(state$logLik, result$logLik, tolerance = 1e-12))) {
    stop("The retained parameters do not reproduce the saved EM likelihood.", call. = FALSE)
  }
  score_tol <- result$controls$score_tol
  score_met <- max(abs(state$gradient)) <= score_tol
  if (!identical(score_met, result$converged)) {
    stop("The saved EM convergence flag does not match its marginal-score tolerance.", call. = FALSE)
  }
  config$estimation_control <- list(maxit = result$controls$maxit,
    quad_points = length(quad$nodes), quadrature = quad,
    mml_integration = "fixed", mml_engine_requested = "em", mml_engine_used = "em",
    optimizer_requested = "BFGS", optimizer_used = "BFGS",
    em_score_tolerance = score_tol, em_mstep_maxit = result$controls$mstep_maxit)
  opt <- list(par = result$par, value = -state$logLik,
    convergence = if (score_met) 0L else if (identical(result$reason, "iteration_limit")) 1L else 2L,
    counts = c("function" = NA_integer_, "gradient" = NA_integer_), message = result$reason,
    gradient = state$gradient * problem$n_person,
    mml_engine = list(Requested = "em", Used = "em", Fallback = FALSE,
      Detail = "Generalized EM with ascent-checked BFGS M steps and a marginal-score stopping rule per Person.",
      EMIterations = as.integer(tail(result$trace$iteration, 1)), EMConverged = score_met,
      EMRelativeChange = NA_real_), em_trace = result$trace)
  opt$optimizer_diagnostics <- build_optimizer_diagnostics(opt,
    gradient = state$gradient, reltol = NA_real_, optimizer_method = "BFGS",
    convergence_basis = "marginal_score_per_person", gradient_tolerance = score_tol)

  mfrm_gmfrm_assemble_fit(problem, config, opt, result$controls,
    result$nonlinear_transformation)
}

# The two numerical engines share model identity and output restrictions;
# the optimizer and integration metadata always describe the actual objective.
mfrm_gmfrm_assemble_fit <- function(problem, config, opt, controls, transformation) {
  common <- problem$common
  prep <- common$prep; sizes <- common$sizes
  quad <- problem$specification$quadrature
  config$method_input <- "MML"
  config$posterior_basis <- "fixed_standard_normal"
  config$step_facet_source <- "explicit"
  config$step_facet_note <- "The second slope facet owns the category steps."
  config$gpcm_mml_identification_requested <- "fixed_standard_normal"
  config$gpcm_estimator_family <- "marginal_maximum_likelihood"
  config$gpcm_slope_action <- "complete_adjacent_predictor"
  config$gpcm_latent_dimension_count <- 1L
  config$gpcm_statistical_penalty <- "none"
  config$gpcm_finite_parameter_box <- FALSE
  config$gpcm_extreme_person_policy <- "posterior_eap_under_population_model"
  config$estimability_audit <- audit_mfrm_estimability(prep, common$idx, config, sizes)
  config$estimability_audit$nonlinear_transformation <- transformation
  config$category_support_audit <- audit_mfrm_category_support(prep, config, sizes)
  data_review <- build_mfrm_data_review(prep,
    estimability_audit = config$estimability_audit,
    category_support_audit = config$category_support_audit)
  params <- expand_params(opt$par, sizes, config)
  mfrmr_fit_progress("Assembling two-family scores and diagnostics")
  person <- build_person_table("MML", common$idx, config, params, prep,
    quad_points = length(quad$nodes), quad = quad)
  config$boundary_audit <- audit_mfrm_person_boundary(person, config, prep)
  data_review$boundary <- config$boundary_audit
  person <- apply_mfrm_person_boundary(person, config$boundary_audit)
  slopes <- build_slope_table(config, prep, params)
  readiness <- build_mfrm_readiness_record(prep, data_review, config, opt, slopes)
  person <- apply_mfrm_person_source_readiness(person, readiness)
  slopes <- apply_mfrm_slope_readiness(slopes, readiness)
  structure(list(
    summary = build_estimation_summary("GPCM", "MML", prep, config, sizes, opt, readiness),
    facets = list(person = person, others = build_other_facet_table(config, prep, params)),
    steps = build_step_table(config, prep, params), slopes = slopes,
    interactions = list(effects = build_interaction_effect_table(config, prep, params),
      specs = config$interaction_specs),
    readiness = readiness, config = config, prep = prep, data_review = data_review, opt = opt,
    population = list(active = FALSE, posterior_basis = "fixed_standard_normal"),
    gmfrm = list(specification = problem$specification, controls = controls)),
    class = c("mfrm_fit", "list"))
}

# Scoped dispatch from fit_mfrm(); never silently drop a requested model option.
mfrm_fit_product_slopes <- function(args, supplied) {
  owners <- unname(args$slope_facet)
  if (!identical(args$model, "GPCM") || !identical(args$method, "MML") ||
      length(owners) != 2L || anyDuplicated(owners) ||
      length(args$facets) != 2L || !setequal(owners, args$facets)) {
    stop("Two slope families require model = 'GPCM', method = 'MML', and exactly two distinct slope_facet names matching facets.", call. = FALSE)
  }
  em <- identical(args$mml_engine, "em") && identical(args$mml_integration, "fixed")
  adaptive <- identical(args$mml_engine, "direct") && identical(args$mml_integration, "adaptive")
  if ((!em && !adaptive) ||
      !identical(args$gpcm_mml_identification, "fixed_standard_normal") ||
      !identical(args$step_facet, owners[2]) ||
      !identical(args$noncenter_facet, owners[2])) {
    stop(paste0("Two slope families require fixed integration with mml_engine = 'em', or adaptive integration with mml_engine = 'direct'; ",
      "gpcm_mml_identification = 'fixed_standard_normal', and both step_facet and noncenter_facet equal to the second slope facet ('",
      owners[2], "'). The first family has geometric mean one; the second is free on the fixed ability scale."), call. = FALSE)
  }
  allowed <- c("data", "person", "facets", "score", "rating_min", "rating_max",
    "keep_original", "category_policy", "model", "method", "step_facet", "slope_facet",
    "noncenter_facet", "quad_points", "maxit", "optimizer", "mml_engine",
    "gpcm_mml_identification", "mml_integration", if (em) "em_score_tol" else c("reltol", "gpcm_mml_start"))
  unsupported <- setdiff(supplied, allowed)
  if (length(unsupported)) {
    stop("These arguments are not supported by the two-family route: ",
      paste(unsupported, collapse = ", "),
      ". It uses observed unweighted ratings, no anchors or population covariates, and no automatic diagnostics. Use em_score_tol for fixed-grid EM or reltol for adaptive direct MML.", call. = FALSE)
  }
  if ((em && !args$optimizer %in% c("auto", "BFGS")) || args$quad_points < 2L) {
    stop("Two-family EM requires optimizer = 'auto' or 'BFGS'; both engines require at least two quadrature points.", call. = FALSE)
  }
  tol <- if (em) args$em_score_tol %||% 1e-6 else args$reltol
  if (!is.numeric(tol) || is.complex(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop(if (em) "`em_score_tol` must be a finite positive number." else
      "`reltol` must be a finite positive number.", call. = FALSE)
  }
  if (!is.null(args$rating_min) && (!is.numeric(args$rating_min) || is.complex(args$rating_min) ||
      length(args$rating_min) != 1L || !is.finite(args$rating_min) || args$rating_min != 0)) {
    stop("The two-family route currently requires rating_min = 0 and integer scores on that scale.", call. = FALSE)
  }
  score <- args$data[[args$score]]
  if (!is.numeric(score) || is.complex(score) || !length(score) || any(!is.finite(score))) {
    stop("The two-family route requires observed finite numeric scores; missing ratings must not be passed as scored observations.", call. = FALSE)
  }
  max_score <- args$rating_max %||% if (is.numeric(score) && length(score) && !anyNA(score)) max(score) else NA_real_
  if (!is.numeric(max_score) || is.complex(max_score) || length(max_score) != 1L ||
      !is.finite(max_score) || max_score < 1 || max_score > .Machine$integer.max || max_score != floor(max_score)) {
    stop("`rating_max` must declare a positive integer upper score for the two-family route.", call. = FALSE)
  }
  if (!isTRUE(args$keep_original) &&
      (!is.numeric(score) || anyNA(score) || !is.finite(max_score) ||
       !identical(sort(unique(as.numeric(score))), as.numeric(seq.int(0, max_score))))) {
    stop("The two-family route does not collapse score categories. Supply the complete observed 0:rating_max scale or category_policy = 'preserve' with the declared rating_max.", call. = FALSE)
  }
  for (name in c(args$person, owners)) {
    values <- as.character(args$data[[name]])
    if (anyNA(values) || any(values != trimws(values))) {
      stop("Two-family model IDs must be observed and have no surrounding whitespace: ", name, ".", call. = FALSE)
    }
  }
  problem <- mfrm_gmfrm_problem(args$data, max_score,
    gauss_hermite_normal(args$quad_points), owners, args$person, args$score)
  common <- problem$common
  category <- audit_mfrm_category_support(common$prep, common$config, common$sizes)
  if (identical(category$readiness$CategoryState[1], "unsupported_coordinate")) {
    mfrmr_stop_unsupported_category(category)
  }
  if (em) {
    result <- mfrm_gmfrm_em(problem, maxit = args$maxit, score_tol = tol)
    fit <- mfrm_gmfrm_fit_result(problem, result)
  } else {
    config <- common$config
    config$estimation_control <- list(maxit = args$maxit, reltol = tol,
      quad_points = args$quad_points, quadrature = problem$specification$quadrature,
      mml_integration = "adaptive", mml_engine_requested = "direct", mml_engine_used = "direct",
      gradient_review_per_person = 1e-6,
      optimizer_requested = args$optimizer)
    policy <- if (is.list(args$gpcm_mml_start)) args$gpcm_mml_start else
      match.arg(args$gpcm_mml_start %||% "neutral_em", c("neutral_em", "neutral"))
    config$estimation_control$gpcm_mml_start <- if (is.list(policy)) "explicit" else policy
    opt <- mfrm_gmfrm_adaptive_fit(problem, config, args$quad_points, args$maxit,
      tol, args$optimizer, policy)
    config$estimation_control$optimizer_used <- opt$optimizer_plan$Used
    fit <- mfrm_gmfrm_assemble_fit(problem, config, opt,
      list(maxit = args$maxit, reltol = tol),
      mfrmr_nonlinear_transformation_audit(opt$par, common$sizes, config, "log_slopes"))
  }
  # Retain a replayable public call; defaults in the single-family route do not
  # describe this likelihood and must not be substituted on reconstruction.
  replay <- args[allowed]
  replay$data <- NULL
  replay$rating_min <- 0
  replay$rating_max <- max_score
  if (em) replay$em_score_tol <- tol else replay$reltol <- tol
  if (adaptive) replay$gpcm_mml_start <- policy
  replay$package_version <- as.character(utils::packageVersion("mfrmr"))
  fit$config$replay_inputs <- replay
  fit$config$source_columns <- list(person = args$person, facets = args$facets, score = args$score, weight = NULL)
  fit$config$attached_diagnostics <- FALSE
  fit$config$estimation_control$optimizer_requested <- args$optimizer
  fit$config$public_product_slopes <- TRUE
  fit
}

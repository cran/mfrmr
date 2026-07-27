# ==============================================================================
# Optimizer dispatch and MML / EM scaffolding
# ==============================================================================
#
# Internal helpers for running the underlying optim() / nlminb() loops
# and the MML-EM hybrid scaffolding. Split out of `mfrm_core.R` so the
# engine-dispatch layer lives in a single file. The
# functions are internal (no @export); they are called from
# `mfrm_estimate()` once estimation configuration / parameter cache /
# initial values have been built upstream.

build_mfrm_mml_em_state <- function(par, idx, config, sizes, quad) {
  params <- expand_params(par, sizes, config)
  base_eta <- compute_base_eta(idx, params, config)
  logprob_bundle <- mfrm_mml_logprob_bundle(
    idx = idx,
    config = config,
    quad = quad,
    params = params,
    base_eta = base_eta
  )
  posterior_bundle <- mfrm_mml_posterior_bundle(logprob_bundle)

  list(
    params = params,
    base_eta = base_eta,
    logprob_bundle = logprob_bundle,
    posterior_bundle = posterior_bundle,
    marginal_loglik = sum(posterior_bundle$person_bundle$log_marginal)
  )
}

mfrm_grad_mml_complete_data_core <- function(params,
                                             base_eta,
                                             idx,
                                             config,
                                             sizes,
                                             quad,
                                             obs_posterior,
                                             step_cum = NULL,
                                             logprob_bundle = NULL) {
  if (identical(config$model, "GPCM")) {
    stop("Complete-data EM updates are currently implemented only for RSM/PCM.",
         call. = FALSE)
  }
  if (isTRUE(config$population_spec$active)) {
    stop("Complete-data EM updates are currently implemented only when `population = NULL`.",
         call. = FALSE)
  }
  if (facet_interactions_active(config)) {
    stop("Complete-data EM updates are currently implemented only for additive RSM/PCM fits without model-estimated facet interactions.",
         call. = FALSE)
  }

  n <- length(idx$score_k)
  if (n == 0L) return(rep(0, sum(unlist(sizes))))

  score_k <- idx$score_k
  weight <- idx$weight
  if (is.null(logprob_bundle)) {
    logprob_bundle <- mfrm_mml_logprob_bundle(
      idx = idx,
      config = config,
      quad = quad,
      params = params,
      base_eta = base_eta,
      step_cum = step_cum,
      include_probs = TRUE
    )
  }
  n_nodes <- ncol(obs_posterior)
  grad_facets_exp <- lapply(config$facet_names, function(f) numeric(length(params$facets[[f]])))
  names(grad_facets_exp) <- config$facet_names

  if (identical(config$model, "RSM")) {
    k_cat <- ncol(logprob_bundle$prob_list[[1]])
    n_steps <- k_cat - 1L
    k_vals <- 0:(k_cat - 1L)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    grad_step_centered <- numeric(n_steps)

    for (q in seq_len(n_nodes)) {
      probs_q <- logprob_bundle$prob_list[[q]]
      expected_q <- as.vector(probs_q %*% k_vals)
      residual_q <- score_k - expected_q
      if (!is.null(weight)) residual_q <- residual_q * weight

      obs_post_q <- obs_posterior[, q]
      w_residual <- residual_q * obs_post_q

      for (facet in config$facet_names) {
        sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
        rs <- rowsum(matrix(sign_f * w_residual, ncol = 1), idx$facets[[facet]], reorder = FALSE)
        f_ids <- as.integer(rownames(rs))
        grad_facets_exp[[facet]][f_ids] <- grad_facets_exp[[facet]][f_ids] + as.vector(rs)
      }

      step_resid <- (compute_P_geq(probs_q) - I_geq) * obs_post_q
      if (!is.null(weight)) step_resid <- step_resid * weight
      grad_step_centered <- grad_step_centered + colSums(step_resid)
    }

    grad_step_free <- project_sum_zero_gradient(grad_step_centered)
  } else {
    k_cat <- ncol(logprob_bundle$prob_list[[1]])
    n_steps <- k_cat - 1L
    n_criteria <- nrow(params$steps_mat)
    k_vals <- 0:(k_cat - 1L)
    I_geq <- outer(score_k, seq_len(n_steps), ">=") * 1.0
    grad_step_mat <- matrix(0, n_criteria, n_steps)

    for (q in seq_len(n_nodes)) {
      probs_q <- logprob_bundle$prob_list[[q]]
      expected_q <- as.vector(probs_q %*% k_vals)
      residual_q <- score_k - expected_q
      if (!is.null(weight)) residual_q <- residual_q * weight

      obs_post_q <- obs_posterior[, q]
      w_residual <- residual_q * obs_post_q

      for (facet in config$facet_names) {
        sign_f <- if (!is.null(config$facet_signs[[facet]])) config$facet_signs[[facet]] else -1
        rs <- rowsum(matrix(sign_f * w_residual, ncol = 1), idx$facets[[facet]], reorder = FALSE)
        f_ids <- as.integer(rownames(rs))
        grad_facets_exp[[facet]][f_ids] <- grad_facets_exp[[facet]][f_ids] + as.vector(rs)
      }

      step_resid <- (compute_P_geq(probs_q) - I_geq) * obs_post_q
      if (!is.null(weight)) step_resid <- step_resid * weight
      rs_step <- rowsum(step_resid, idx$step_idx, reorder = FALSE)
      rs_ids <- as.integer(rownames(rs_step))
      grad_step_mat[rs_ids, ] <- grad_step_mat[rs_ids, ] + rs_step
    }

    grad_step_free <- project_step_matrix_gradient(grad_step_mat)
  }

  grad_facet_free <- unlist(lapply(config$facet_names, function(f) {
    constraint_grad_project(grad_facets_exp[[f]], config$facet_specs[[f]])
  }))

  -c(grad_facet_free, grad_step_free)
}

normalize_mfrm_optimizer <- function(optimizer = "auto") {
  value <- as.character(optimizer %||% "auto")[1]
  aliases <- c(
    auto = "auto",
    bfgs = "BFGS",
    `l-bfgs-b` = "L-BFGS-B",
    lbfgsb = "L-BFGS-B",
    `l_bfgs_b` = "L-BFGS-B"
  )
  key <- tolower(value)
  resolved <- unname(aliases[key])
  if (length(resolved) == 0L || is.na(resolved)) {
    stop("`optimizer` must be one of 'auto', 'BFGS', or 'L-BFGS-B'.",
         call. = FALSE)
  }
  resolved
}

resolve_mfrm_optimizer <- function(optimizer, n_parameters,
                                   prefer_limited_memory = FALSE) {
  requested <- normalize_mfrm_optimizer(optimizer)
  # Fixed rather than option-driven so `optimizer = "auto"` is reproducible
  # across sessions and replayed analysis bundles.
  threshold <- 200L
  used <- if (identical(requested, "auto")) {
    if (isTRUE(prefer_limited_memory) ||
        as.integer(n_parameters) >= threshold) "L-BFGS-B" else "BFGS"
  } else {
    requested
  }
  list(
    Requested = requested,
    Used = used,
    ParameterCount = as.integer(n_parameters),
    AutoThreshold = threshold
  )
}

build_mfrm_optim_control <- function(method, maxit, reltol) {
  if (identical(method, "L-BFGS-B")) {
    return(list(
      maxit = as.integer(maxit),
      lmm = 20L,
      factr = max(1, as.numeric(reltol) / .Machine$double.eps),
      pgtol = max(as.numeric(reltol), sqrt(.Machine$double.eps))
    ))
  }
  list(maxit = as.integer(maxit), reltol = as.numeric(reltol))
}

mfrm_optimizer_polish_tolerances <- function(reltol) {
  reltol <- as.numeric(reltol)
  candidates <- c(1e-11, 1e-13, 3e-14, 1e-14)
  unique(candidates[is.finite(candidates) & candidates < reltol])
}

mfrm_optimizer_stage_row <- function(stage, selected = FALSE) {
  data.frame(
    Stage = as.integer(stage$index),
    StageLabel = as.character(stage$label),
    Method = as.character(stage$method),
    Reltol = as.numeric(stage$reltol),
    Factr = as.numeric(stage$control$factr %||% NA_real_),
    Pgtol = as.numeric(stage$control$pgtol %||% NA_real_),
    ConvergenceCode = as.integer(stage$diagnostics$ConvergenceCode %||% NA_integer_),
    ConvergenceStatus = as.character(stage$diagnostics$ConvergenceStatus %||% "unknown"),
    StopReason = as.character(stage$diagnostics$ConvergenceReason %||% "unknown"),
    ConvergenceSeverity = as.character(stage$diagnostics$ConvergenceSeverity %||% "review"),
    InferenceReady = identical(stage$diagnostics$ConvergenceSeverity, "pass"),
    TerminalGradientSupNorm = as.numeric(
      stage$diagnostics$TerminalGradientSupNorm %||% NA_real_
    ),
    TerminalGradientRMS = as.numeric(
      stage$diagnostics$TerminalGradientRMS %||% NA_real_
    ),
    FunctionEvaluations = as.integer(
      stage$diagnostics$FunctionEvaluations %||% NA_integer_
    ),
    GradientEvaluations = as.integer(
      stage$diagnostics$GradientEvaluations %||% NA_integer_
    ),
    Objective = as.numeric(stage$opt$value %||% NA_real_),
    MaxParameterChange = as.numeric(stage$max_parameter_change %||% NA_real_),
    ElapsedSeconds = as.numeric(stage$elapsed %||% NA_real_),
    Selected = isTRUE(selected),
    Error = as.character(stage$error %||% ""),
    stringsAsFactors = FALSE
  )
}

mfrm_optimizer_stage_is_better <- function(candidate, current) {
  candidate_value <- as.numeric(candidate$opt$value %||% Inf)
  current_value <- as.numeric(current$opt$value %||% Inf)
  if (!is.finite(candidate_value)) return(FALSE)
  if (is.finite(current_value)) {
    objective_slack <- 1e-8 * max(1, abs(current_value))
    if (candidate_value > current_value + objective_slack) return(FALSE)
  }

  candidate_ready <- identical(
    candidate$diagnostics$ConvergenceSeverity %||% "review", "pass"
  )
  current_ready <- identical(
    current$diagnostics$ConvergenceSeverity %||% "review", "pass"
  )
  if (candidate_ready != current_ready) return(candidate_ready)

  candidate_gradient <- as.numeric(
    candidate$diagnostics$TerminalGradientSupNorm %||% Inf
  )
  current_gradient <- as.numeric(
    current$diagnostics$TerminalGradientSupNorm %||% Inf
  )
  if (!is.finite(candidate_gradient)) candidate_gradient <- Inf
  if (!is.finite(current_gradient)) current_gradient <- Inf
  if (!isTRUE(all.equal(candidate_gradient, current_gradient, tolerance = 0))) {
    return(candidate_gradient < current_gradient)
  }

  is.finite(candidate_value) && candidate_value < current_value
}

# Lazily fused objective / gradient evaluator. optim() may request several
# objective values during a line search without requesting their gradients, so
# the reusable probability workspace is built for paired objective / gradient
# methods while the gradient itself is constructed only when gr() asks for it.
# At an identical parameter vector, marginal-likelihood and posterior work is
# reused; full probability matrices are retained only when the optimizer's
# call pattern makes that reuse worthwhile.
make_mfrm_direct_evaluator <- function(method, cache, idx, config, sizes, quad,
                                       reuse_probability_workspace = TRUE) {
  cached_par <- NULL
  cached_value <- NULL
  cached_gradient <- NULL
  probability_bundle <- NULL
  logprob_bundle <- NULL
  posterior_bundle <- NULL
  shared_evaluations <- 0L
  gradient_builds <- 0L
  value_hits <- 0L
  gradient_hits <- 0L

  ensure_shared <- function(par) {
    if (identical(par, cached_par)) return(invisible(NULL))
    cache$ensure(par)
    cached_par <<- par
    cached_gradient <<- NULL
    probability_bundle <<- NULL
    logprob_bundle <<- NULL
    posterior_bundle <<- NULL
    shared_evaluations <<- shared_evaluations + 1L

    if (identical(method, "JML")) {
      if (isTRUE(reuse_probability_workspace)) {
        probability_bundle <<- mfrm_jml_probability_bundle(
          eta = cache$eta(),
          score_k = idx$score_k,
          model = config$model,
          step_cum = cache$step_cum(),
          criterion_idx = idx$step_idx,
          slopes = cache$params()$slopes,
          slope_idx = idx$slope_idx
        )
        log_prob_obs <- probability_bundle$log_prob_obs
        if (!is.null(idx$weight)) log_prob_obs <- log_prob_obs * idx$weight
        cached_value <<- -sum(log_prob_obs)
      } else {
        cached_value <<- mfrm_loglik_jml_cached(cache, idx, config)
      }
    } else if (length(idx$score_k) == 0L) {
      cached_value <<- 0
    } else {
      logprob_bundle <<- mfrm_mml_logprob_bundle(
        idx = idx,
        config = config,
        quad = quad,
        params = cache$params(),
        base_eta = cache$base_eta(),
        step_cum = cache$step_cum(),
        include_probs = isTRUE(reuse_probability_workspace),
        include_linear_part = isTRUE(reuse_probability_workspace) &&
          identical(config$model, "GPCM")
      )
      posterior_bundle <<- mfrm_mml_posterior_bundle(logprob_bundle)
      cached_value <<- -sum(
        posterior_bundle$person_bundle$log_marginal
      )
    }
    invisible(NULL)
  }

  list(
    value = function(par) {
      same <- identical(par, cached_par)
      ensure_shared(par)
      if (same) value_hits <<- value_hits + 1L
      cached_value
    },
    gradient = function(par) {
      ensure_shared(par)
      if (!is.null(cached_gradient)) {
        gradient_hits <<- gradient_hits + 1L
        return(cached_gradient)
      }
      cached_gradient <<- if (identical(method, "JML")) {
        mfrm_grad_jml_cached(
          cache, idx, config, sizes,
          probability_bundle = probability_bundle
        )
      } else {
        mfrm_grad_mml_cached(
          cache, idx, config, sizes, quad,
          logprob_bundle = if (isTRUE(reuse_probability_workspace)) {
            logprob_bundle
          } else {
            NULL
          },
          posterior_bundle = posterior_bundle
        )
      }
      gradient_builds <<- gradient_builds + 1L
      cached_gradient
    },
    diagnostics = function() list(
      SharedEvaluations = shared_evaluations,
      GradientBuilds = gradient_builds,
      ValueCacheHits = value_hits,
      GradientCacheHits = gradient_hits
    )
  )
}

make_mfrm_em_mstep_evaluator <- function(cache, idx, config, sizes, quad,
                                         obs_posterior,
                                         reuse_probability_workspace = TRUE) {
  cached_par <- NULL
  cached_value <- NULL
  cached_gradient <- NULL
  logprob_bundle <- NULL

  ensure_shared <- function(par) {
    if (identical(par, cached_par)) return(invisible(NULL))
    cache$ensure(par)
    cached_par <<- par
    cached_gradient <<- NULL
    logprob_bundle <<- mfrm_mml_logprob_bundle(
      idx = idx,
      config = config,
      quad = quad,
      params = cache$params(),
      base_eta = cache$base_eta(),
      step_cum = cache$step_cum(),
      include_probs = isTRUE(reuse_probability_workspace)
    )
    cached_value <<- -sum(logprob_bundle$log_prob_mat * obs_posterior)
    invisible(NULL)
  }

  list(
    value = function(par) {
      ensure_shared(par)
      cached_value
    },
    gradient = function(par) {
      ensure_shared(par)
      if (is.null(cached_gradient)) {
        cached_gradient <<- mfrm_grad_mml_complete_data_core(
          params = cache$params(),
          base_eta = cache$base_eta(),
          idx = idx,
          config = config,
          sizes = sizes,
          quad = quad,
          obs_posterior = obs_posterior,
          step_cum = cache$step_cum(),
          logprob_bundle = if (isTRUE(reuse_probability_workspace)) {
            logprob_bundle
          } else {
            NULL
          }
        )
      }
      cached_gradient
    }
  )
}

run_mfrm_direct_optimization <- function(start,
                                         method,
                                         idx,
                                         config,
                                         sizes,
                                         quad_points,
                                         maxit,
                                         reltol,
                                         quad = NULL,
                                         optimizer = "auto",
                                         suppress_convergence_warning = FALSE) {
  optimizer_plan <- resolve_mfrm_optimizer(
    optimizer,
    length(start),
    prefer_limited_memory = identical(method, "MML")
  )
  if (method == "JML") {
    cache <- make_param_cache(sizes, config, idx, is_mml = FALSE)
  } else {
    quad <- quad %||% gauss_hermite_normal(quad_points)
    cache <- make_param_cache(sizes, config, idx, is_mml = TRUE)
  }

  evaluator <- make_mfrm_direct_evaluator(
    method = method,
    cache = cache,
    idx = idx,
    config = config,
    sizes = sizes,
    quad = quad,
    reuse_probability_workspace = identical(optimizer_plan$Used, "L-BFGS-B")
  )
  fn <- function(par, ...) evaluator$value(par)
  gr <- function(par, ...) evaluator$gradient(par)

  run_stage <- function(par, stage_method, stage_reltol, index, label,
                        fail_hard = FALSE) {
    started <- proc.time()[["elapsed"]]
    stage_control <- build_mfrm_optim_control(
      stage_method,
      maxit = maxit,
      reltol = stage_reltol
    )
    stage_opt <- tryCatch(
      optim(
        par = par,
        fn = fn,
        gr = gr,
        method = stage_method,
        control = stage_control
      ),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - started
    if (inherits(stage_opt, "error")) {
      if (isTRUE(fail_hard)) {
        stop("Model optimization failed: ", conditionMessage(stage_opt), ". ",
             "Possible causes: (1) insufficient data for the number of parameters, ",
             "(2) extreme score distributions, (3) near-constant responses. ",
             "Try reducing facets, increasing maxit, or checking data quality.",
             call. = FALSE)
      }
      placeholder <- list(
        par = par,
        value = NA_real_,
        convergence = NA_integer_,
        counts = c("function" = NA_integer_, "gradient" = NA_integer_),
        message = conditionMessage(stage_opt)
      )
      diagnostics <- build_optimizer_diagnostics(
        opt = placeholder,
        gradient = rep(NA_real_, length(par)),
        reltol = stage_reltol,
        maxit = maxit,
        optimizer_method = stage_method,
        convergence_basis = "optimizer_gradient"
      )
      return(list(
        index = index,
        label = label,
        method = stage_method,
        reltol = stage_reltol,
        control = stage_control,
        opt = placeholder,
        diagnostics = diagnostics,
        elapsed = elapsed,
        max_parameter_change = 0,
        error = conditionMessage(stage_opt)
      ))
    }

    final_gradient <- tryCatch(
      gr(stage_opt$par),
      error = function(e) rep(NA_real_, length(stage_opt$par))
    )
    diagnostics <- build_optimizer_diagnostics(
      opt = stage_opt,
      gradient = final_gradient,
      reltol = stage_reltol,
      maxit = maxit,
      optimizer_method = stage_method,
      convergence_basis = "optimizer_gradient"
    )
    list(
      index = index,
      label = label,
      method = stage_method,
      reltol = stage_reltol,
      control = stage_control,
      opt = stage_opt,
      diagnostics = diagnostics,
      elapsed = elapsed,
      max_parameter_change = if (length(stage_opt$par) > 0L) {
        max(abs(as.numeric(stage_opt$par) - as.numeric(par)), na.rm = TRUE)
      } else {
        0
      },
      error = ""
    )
  }

  initial <- run_stage(
    par = start,
    stage_method = optimizer_plan$Used,
    stage_reltol = reltol,
    index = 1L,
    label = "initial",
    fail_hard = TRUE
  )
  stages <- list(initial)
  selected <- initial
  polish_triggered <- identical(initial$opt$convergence, 0L) &&
    identical(initial$diagnostics$ConvergenceReason, "code_zero_large_gradient") &&
    is.finite(reltol) && reltol <= 1e-9

  if (isTRUE(polish_triggered)) {
    stage_index <- 1L
    for (stage_reltol in mfrm_optimizer_polish_tolerances(reltol)) {
      stage_index <- stage_index + 1L
      candidate <- run_stage(
        par = selected$opt$par,
        stage_method = optimizer_plan$Used,
        stage_reltol = stage_reltol,
        index = stage_index,
        label = "gradient_polish",
        fail_hard = FALSE
      )
      stages[[length(stages) + 1L]] <- candidate
      if (!nzchar(candidate$error) &&
          mfrm_optimizer_stage_is_better(candidate, selected)) {
        selected <- candidate
      }
      if (identical(selected$diagnostics$ConvergenceSeverity, "pass")) break
    }

    # A limited BFGS fallback can recover a small or moderate parameter vector
    # when L-BFGS-B's relative-objective stop remains ahead of the common
    # gradient gate. Avoid the quadratic-memory method for large JML vectors.
    fallback_limit <- 600L
    if (!identical(selected$diagnostics$ConvergenceSeverity, "pass") &&
        identical(optimizer_plan$Used, "L-BFGS-B") &&
        length(start) <= fallback_limit) {
      fallback_tolerances <- c(1e-13, 1e-14)
      for (stage_reltol in fallback_tolerances) {
        stage_index <- stage_index + 1L
        candidate <- run_stage(
          par = selected$opt$par,
          stage_method = "BFGS",
          stage_reltol = stage_reltol,
          index = stage_index,
          label = "bounded_bfgs_fallback",
          fail_hard = FALSE
        )
        stages[[length(stages) + 1L]] <- candidate
        if (!nzchar(candidate$error) &&
            mfrm_optimizer_stage_is_better(candidate, selected)) {
          selected <- candidate
        }
        if (identical(selected$diagnostics$ConvergenceSeverity, "pass")) break
      }
    }
  }

  opt <- selected$opt
  total_functions <- sum(vapply(stages, function(stage) {
    as.numeric(stage$diagnostics$FunctionEvaluations %||% 0)
  }, numeric(1)), na.rm = TRUE)
  total_gradients <- sum(vapply(stages, function(stage) {
    as.numeric(stage$diagnostics$GradientEvaluations %||% 0)
  }, numeric(1)), na.rm = TRUE)
  opt$optimizer_diagnostics <- selected$diagnostics
  opt$optimizer_diagnostics$RequestedReltol <- as.numeric(reltol)
  opt$optimizer_diagnostics$EffectiveReltol <- as.numeric(selected$reltol)
  opt$optimizer_diagnostics$OptimizerFactr <- as.numeric(
    selected$control$factr %||% NA_real_
  )
  opt$optimizer_diagnostics$OptimizerPgtol <- as.numeric(
    selected$control$pgtol %||% NA_real_
  )
  opt$optimizer_diagnostics$FunctionEvaluations <- as.integer(total_functions)
  opt$optimizer_diagnostics$GradientEvaluations <- as.integer(total_gradients)
  optimizer_plan$InitialUsed <- optimizer_plan$Used
  optimizer_plan$Used <- selected$method
  optimizer_plan$PolishTriggered <- isTRUE(polish_triggered)
  optimizer_plan$PolishSucceeded <- isTRUE(polish_triggered) &&
    identical(selected$diagnostics$ConvergenceSeverity, "pass")
  optimizer_plan$PolishStages <- length(stages) - 1L
  opt$optimizer_plan <- optimizer_plan
  stage_table <- do.call(rbind, lapply(stages, function(stage) {
    mfrm_optimizer_stage_row(
      stage,
      selected = identical(stage$index, selected$index)
    )
  }))
  opt$optimizer_polish <- list(
    Triggered = isTRUE(polish_triggered),
    Succeeded = isTRUE(optimizer_plan$PolishSucceeded),
    InitialMethod = as.character(initial$method),
    FinalMethod = as.character(selected$method),
    RequestedReltol = as.numeric(reltol),
    EffectiveReltol = as.numeric(selected$reltol),
    EffectiveFactr = as.numeric(selected$control$factr %||% NA_real_),
    EffectivePgtol = as.numeric(selected$control$pgtol %||% NA_real_),
    InitialTerminalGradientSupNorm = as.numeric(
      initial$diagnostics$TerminalGradientSupNorm
    ),
    FinalTerminalGradientSupNorm = as.numeric(
      selected$diagnostics$TerminalGradientSupNorm
    ),
    SelectedStage = as.integer(selected$index),
    Stages = stage_table
  )
  opt$evaluation_cache <- evaluator$diagnostics()

  if (!identical(opt$optimizer_diagnostics$ConvergenceSeverity, "pass") &&
      !isTRUE(suppress_convergence_warning)) {
    warning(
      "Optimization convergence review did not produce an inference-ready numerical solution ",
      "(code = ", opt$convergence,
      ", status = ", opt$optimizer_diagnostics$ConvergenceStatus, "). ",
      opt$optimizer_diagnostics$ConvergenceDetail, " ",
      if (isTRUE(polish_triggered)) {
        paste0("Bounded gradient polishing attempted ", length(stages) - 1L,
               " additional stage(s). ")
      } else {
        ""
      },
      "Inspect the model specification, data support, and starting values. ",
      "Do not interpret estimates until the review is resolved.",
      call. = FALSE
    )
  }

  opt
}

run_mfrm_mml_em_optimization <- function(start,
                                         idx,
                                         config,
                                         sizes,
                                         quad_points,
                                         maxit,
                                         reltol,
                                         m_step_maxit = NULL,
                                         m_step_reltol = NULL,
                                         optimizer = "auto",
                                         suppress_convergence_warning = FALSE,
                                         checkpoint = NULL) {
  quad <- gauss_hermite_normal(quad_points)
  par <- as.numeric(start)
  prev_loglik <- -Inf
  converged <- FALSE
  ll_trace <- numeric(0)
  rel_change <- NA_real_
  total_fn <- 0L
  total_gr <- 0L
  optimizer_plan <- resolve_mfrm_optimizer(
    optimizer,
    length(start),
    prefer_limited_memory = TRUE
  )

  # Resumable-fit checkpoint scaffolding. `checkpoint` is a list with
  # `file` (path) and `every_iter` (integer >= 1). When set, the
  # current EM state is `saveRDS()`-ed every `every_iter` outer EM
  # iterations so a long fit can resume after a crash via
  # `resume_mfrm_fit()`. We also try to load an existing checkpoint
  # before the first iteration; when found, `par` and `it` are
  # reseeded. The checkpoint format is intentionally tied to this
  # function's local state -- it should not be hand-edited.
  ckpt_file <- NULL
  ckpt_every <- 1L
  start_it <- 1L
  if (!is.null(checkpoint)) {
    if (!is.list(checkpoint) || is.null(checkpoint$file) ||
        !nzchar(as.character(checkpoint$file))) {
      stop("`checkpoint` must be a list with a non-empty `file` path.",
           call. = FALSE)
    }
    ckpt_file <- as.character(checkpoint$file)
    ckpt_every <- max(1L, as.integer(checkpoint$every_iter %||% 1L))
    if (file.exists(ckpt_file)) {
      saved <- tryCatch(readRDS(ckpt_file), error = function(e) NULL)
      if (is.list(saved) && identical(saved$.mfrm_checkpoint_kind %||% "",
                                       "mml_em")) {
        par <- as.numeric(saved$par)
        prev_loglik <- as.numeric(saved$prev_loglik %||% -Inf)
        ll_trace <- as.numeric(saved$ll_trace %||% numeric(0))
        total_fn <- as.integer(saved$total_fn %||% 0L)
        total_gr <- as.integer(saved$total_gr %||% 0L)
        start_it <- as.integer(saved$next_iter %||% 1L)
        message("Resumed MML EM from checkpoint at iteration ",
                start_it, " (", ckpt_file, ").")
      } else {
        warning("Existing checkpoint file '", ckpt_file, "' did not look ",
                "like an mfrmr MML EM checkpoint; starting from scratch.",
                call. = FALSE)
      }
    }
  }
  m_step_maxit <- if (is.null(m_step_maxit)) {
    max(5L, min(50L, as.integer(maxit)))
  } else {
    as.integer(m_step_maxit)
  }
  m_step_reltol <- if (is.null(m_step_reltol)) {
    max(as.numeric(reltol), 1e-5)
  } else {
    as.numeric(m_step_reltol)
  }

  for (it in seq.int(start_it, maxit)) {
    state <- build_mfrm_mml_em_state(par, idx, config, sizes, quad)
    ll_trace <- c(ll_trace, state$marginal_loglik)

    if (it > 1L) {
      rel_change <- abs(state$marginal_loglik - prev_loglik) / (abs(prev_loglik) + 1e-10)
      if (is.finite(rel_change) && rel_change < reltol) {
        converged <- TRUE
        break
      }
    }
    prev_loglik <- state$marginal_loglik

    cache <- make_param_cache(sizes, config, idx, is_mml = TRUE)
    obs_posterior_fixed <- state$posterior_bundle$obs_posterior

    evaluator <- make_mfrm_em_mstep_evaluator(
      cache = cache,
      idx = idx,
      config = config,
      sizes = sizes,
      quad = quad,
      obs_posterior = obs_posterior_fixed,
      reuse_probability_workspace = identical(optimizer_plan$Used, "L-BFGS-B")
    )
    fn <- function(par, ...) evaluator$value(par)
    gr <- function(par, ...) evaluator$gradient(par)

    m_opt <- tryCatch(
      optim(
        par = par,
        fn = fn,
        gr = gr,
        method = optimizer_plan$Used,
        control = build_mfrm_optim_control(
          optimizer_plan$Used,
          maxit = m_step_maxit,
          reltol = m_step_reltol
        )
      ),
      error = function(e) {
        stop("EM M-step optimization failed: ", conditionMessage(e), ". ",
             "Try increasing `maxit`, reducing model complexity, or using `mml_engine = 'direct'`.",
             call. = FALSE)
      }
    )

    par <- m_opt$par
    total_fn <- total_fn + as.integer(unname(m_opt$counts[["function"]] %||% 0L))
    total_gr <- total_gr + as.integer(unname(m_opt$counts[["gradient"]] %||% 0L))

    # Periodic checkpoint write. We snapshot the post-M-step state so
    # `resume_mfrm_fit()` continues at the next iteration with the
    # same `par` and accumulated trace. tryCatch shields the fit from
    # transient I/O failures (full disk, permission flap).
    if (!is.null(ckpt_file) && (it %% ckpt_every == 0L)) {
      tryCatch(
        saveRDS(
          list(
            .mfrm_checkpoint_kind = "mml_em",
            par = par,
            prev_loglik = prev_loglik,
            ll_trace = ll_trace,
            total_fn = total_fn,
            total_gr = total_gr,
            next_iter = it + 1L,
            quad_points = quad_points,
            maxit = maxit,
            reltol = reltol,
            timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
          ),
          file = ckpt_file
        ),
        error = function(e) {
          warning("MML EM checkpoint write to '", ckpt_file,
                  "' failed: ", conditionMessage(e),
                  ". Continuing without checkpoint.",
                  call. = FALSE)
        }
      )
    }
  }

  final_state <- build_mfrm_mml_em_state(par, idx, config, sizes, quad)
  final_loglik <- final_state$marginal_loglik
  if (length(ll_trace) == 0L ||
      !isTRUE(isTRUE(all.equal(tail(ll_trace, 1L), final_loglik, tolerance = 1e-12)))) {
    ll_trace <- c(ll_trace, final_loglik)
  }

  final_gradient <- tryCatch(
    mfrm_grad_mml_core(
      params = final_state$params,
      base_eta = final_state$base_eta,
      idx = idx,
      config = config,
      sizes = sizes,
      quad = quad
    ),
    error = function(e) rep(NA_real_, length(par))
  )

  opt <- list(
    par = par,
    value = -final_loglik,
    counts = stats::setNames(c(total_fn, total_gr), c("function", "gradient")),
    convergence = if (isTRUE(converged)) 0L else 1L,
    message = if (isTRUE(converged)) {
      "EM converged by relative log-likelihood change."
    } else {
      "EM reached max iterations before the relative log-likelihood change met the tolerance."
    },
    ll_trace = ll_trace,
    em_relative_change = rel_change,
    em_iterations = as.integer(length(ll_trace) - 1L)
  )
  opt$optimizer_plan <- optimizer_plan

  # `InferenceReady` has one meaning across direct, hybrid, and EM engines:
  # the engine-specific stopping rule must be satisfied and the common
  # terminal-gradient review must pass. Relative log-likelihood change alone
  # can be small while the score vector is still materially non-zero.
  opt$optimizer_diagnostics <- build_optimizer_diagnostics(
    opt = opt,
    gradient = final_gradient,
    reltol = reltol,
    maxit = maxit,
    optimizer_method = "EM",
    convergence_basis = "optimizer_gradient"
  )
  opt$optimizer_diagnostics$ConvergenceBasis <- "relative_loglik_and_gradient"
  opt$optimizer_diagnostics$RequestedReltol <- as.numeric(reltol)
  opt$optimizer_diagnostics$EffectiveReltol <- as.numeric(reltol)
  opt$optimizer_diagnostics$OptimizerFactr <- NA_real_
  opt$optimizer_diagnostics$OptimizerPgtol <- NA_real_
  if (identical(opt$optimizer_diagnostics$ConvergenceSeverity, "pass")) {
    opt$optimizer_diagnostics$ConvergenceReason <-
      "relative_loglik_and_gradient_tolerance_met"
    opt$optimizer_diagnostics$ConvergenceDetail <- paste(
      "EM relative log-likelihood change met the stopping tolerance and the",
      "terminal gradient met the common numerical review tolerance."
    )
  } else if (identical(opt$convergence, 0L) &&
             identical(opt$optimizer_diagnostics$ConvergenceReason,
                       "code_zero_large_gradient")) {
    opt$optimizer_diagnostics$ConvergenceReason <-
      "relative_loglik_met_large_gradient"
    opt$optimizer_diagnostics$ConvergenceDetail <- paste(
      "EM relative log-likelihood change met its stopping tolerance, but the",
      "terminal gradient exceeded the common numerical review tolerance."
    )
  }
  em_stage <- list(
    index = 1L,
    label = "em",
    method = "EM",
    reltol = reltol,
    control = list(reltol = reltol),
    opt = opt,
    diagnostics = opt$optimizer_diagnostics,
    elapsed = NA_real_,
    max_parameter_change = NA_real_,
    error = ""
  )
  opt$optimizer_polish <- list(
    Triggered = FALSE,
    Succeeded = FALSE,
    InitialMethod = "EM",
    FinalMethod = "EM",
    RequestedReltol = as.numeric(reltol),
    EffectiveReltol = as.numeric(reltol),
    EffectiveFactr = NA_real_,
    EffectivePgtol = NA_real_,
    InitialTerminalGradientSupNorm = as.numeric(
      opt$optimizer_diagnostics$TerminalGradientSupNorm
    ),
    FinalTerminalGradientSupNorm = as.numeric(
      opt$optimizer_diagnostics$TerminalGradientSupNorm
    ),
    SelectedStage = 1L,
    Stages = mfrm_optimizer_stage_row(em_stage, selected = TRUE)
  )
  opt$em_diagnostics <- list(
    EMIterations = as.integer(length(ll_trace) - 1L),
    EMConverged = isTRUE(converged),
    EMRelativeChange = rel_change,
    MStepMaxit = m_step_maxit,
    MStepReltol = m_step_reltol,
    MStepOptimizer = optimizer_plan$Used
  )

  if (!identical(opt$optimizer_diagnostics$ConvergenceSeverity, "pass") &&
      !isTRUE(suppress_convergence_warning)) {
    warning("EM did not produce an inference-ready numerical solution (status = ",
            opt$optimizer_diagnostics$ConvergenceStatus, "). ",
            opt$optimizer_diagnostics$ConvergenceDetail, " ",
            "Consider a tighter `reltol`, increasing maxit (current: ", maxit,
            "), or using `mml_engine = 'direct'` for bounded gradient polishing.",
            call. = FALSE)
  }

  opt
}

run_mfrm_optimization <- function(start,
                                  method,
                                  idx,
                                  config,
                                  sizes,
                                  quad_points,
                                  maxit,
                                  reltol,
                                  optimizer = "auto",
                                  suppress_convergence_warning = FALSE,
                                  checkpoint = NULL) {
  requested_engine <- normalize_mml_engine(config$estimation_control$mml_engine_requested %||% "direct")
  engine_plan <- resolve_mml_engine_plan(
    method = method,
    model = config$model,
    requested = requested_engine,
    population_active = isTRUE(config$population_spec$active),
    interaction_active = facet_interactions_active(config)
  )

  if (isTRUE(engine_plan$Fallback) &&
      identical(method, "MML") &&
      !isTRUE(suppress_convergence_warning)) {
    warning(engine_plan$Detail, call. = FALSE)
  }

  if (!identical(method, "MML") || identical(engine_plan$Used, "direct")) {
    opt <- run_mfrm_direct_optimization(
      start = start,
      method = method,
      idx = idx,
      config = config,
      sizes = sizes,
      quad_points = quad_points,
      maxit = maxit,
      reltol = reltol,
      optimizer = optimizer,
      suppress_convergence_warning = suppress_convergence_warning
    )
  } else if (identical(engine_plan$Used, "em")) {
    opt <- run_mfrm_mml_em_optimization(
      start = start,
      idx = idx,
      config = config,
      sizes = sizes,
      quad_points = quad_points,
      maxit = maxit,
      reltol = reltol,
      optimizer = optimizer,
      suppress_convergence_warning = suppress_convergence_warning,
      checkpoint = checkpoint
    )
  } else {
    em_maxit <- compute_hybrid_em_maxit(maxit)
    em_reltol <- compute_hybrid_em_reltol(reltol)
    em_opt <- run_mfrm_mml_em_optimization(
      start = start,
      idx = idx,
      config = config,
      sizes = sizes,
      quad_points = quad_points,
      maxit = em_maxit,
      reltol = em_reltol,
      m_step_maxit = compute_hybrid_em_mstep_maxit(em_maxit),
      m_step_reltol = em_reltol,
      optimizer = optimizer,
      suppress_convergence_warning = TRUE
    )
    opt <- run_mfrm_direct_optimization(
      start = em_opt$par,
      method = method,
      idx = idx,
      config = config,
      sizes = sizes,
      quad_points = quad_points,
      maxit = maxit,
      reltol = reltol,
      optimizer = optimizer,
      suppress_convergence_warning = suppress_convergence_warning
    )
    opt$em_diagnostics <- em_opt$em_diagnostics
    opt$em_warm_start_trace <- em_opt$ll_trace
  }

  if (identical(method, "MML")) {
    em_diag <- opt$em_diagnostics %||% list()
    opt$mml_engine <- list(
      Requested = engine_plan$Requested,
      Used = engine_plan$Used,
      Detail = engine_plan$Detail,
      Fallback = isTRUE(engine_plan$Fallback),
      EMIterations = as.integer(em_diag$EMIterations %||% NA_integer_),
      EMConverged = as.logical(em_diag$EMConverged %||% NA),
      EMRelativeChange = as.numeric(em_diag$EMRelativeChange %||% NA_real_)
    )
  }

  opt
}

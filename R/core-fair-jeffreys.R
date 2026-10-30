# Canonical category logits are offset_k + X_k beta. For centered rows Z_ek,
# I = sum p_ek Z_ek Z_ek' and d log|I| / d beta = sum p_ek q_ek Z_ek,
# where q_ek = Z_ek' I^{-1} Z_ek. Used by the explicit generating-model
# option and by the separate opt-in penalized estimator; it supplies no SE.
mfrm_fair_jeffreys_objective <- function(category, offset, y) {
  n <- nrow(offset); K <- ncol(offset); p <- ncol(category[[1]])
  cache <- new.env(parent = emptyenv())
  cache$par <- NULL; cache$value <- NULL
  evaluate <- function(beta, gradient = FALSE) {
    same <- identical(beta, cache$par)
    if (!same) {
      logits <- offset
      for (k in seq_len(K)) logits[, k] <- logits[, k] + as.numeric(category[[k]] %*% beta)
      logits <- logits - apply(logits, 1L, max)
      logp <- logits - log(rowSums(exp(logits)))
      probs <- exp(logp)
      nll <- -sum(logp[cbind(seq_len(n), y + 1L)])
      Z <- mfrm_fair_delta_center(category, probs)
      I <- mfrm_fair_delta_information(Z, probs)
      diagonal <- as.numeric(Matrix::diag(I))
      if (any(!is.finite(c(beta, nll, diagonal))) || any(diagonal <= 0))
        stop("Jeffreys information has a nonpositive/nonfinite diagonal.", call. = FALSE)
      D <- Matrix::Diagonal(x = 1 / sqrt(diagonal))
      # Keep the CHM-factor contract even when a coordinate transformation or
      # independent reference supplies dense category-design matrices.
      H <- Matrix::forceSymmetric(methods::as(D %*% I %*% D, "CsparseMatrix"))
      factor <- if (p) tryCatch(Matrix::Cholesky(H, LDL = FALSE, perm = TRUE),
        warning = function(w) NULL, error = function(e) NULL) else NULL
      if (p && is.null(factor))
        stop("Jeffreys information is not numerically positive definite; no regularization was applied.", call. = FALSE)
      # Older Matrix releases always return log|L| for a CHM factor. Request
      # that convention explicitly, then double it for log|H| on all versions.
      logdet <- if (p) 2 * as.numeric(Matrix::determinant(factor, logarithm = TRUE, sqrt = TRUE)$modulus) +
        sum(log(diagonal)) else 0
      value <- nll - logdet / 2
      if (!is.finite(value)) stop("Nonfinite Jeffreys penalized objective.", call. = FALSE)
      cache$value <- list(value = value, nll = nll, logdet = logdet, probabilities = probs,
        centered = Z, D = D, H = H, factor = factor, gradient = NULL)
      cache$par <- beta
    }
    z <- cache$value
    if (gradient && is.null(z$gradient)) {
      g <- adjustment <- numeric(p); residual <- 0
      for (k in seq_len(K)) {
        g <- g + as.numeric(Matrix::crossprod(category[[k]], z$probabilities[, k] - (y == k - 1L)))
        q <- numeric(n)
        if (p) for (first in seq.int(1L, n, by = 256L)) {
          rows <- seq.int(first, min(n, first + 255L))
          rhs <- z$D %*% Matrix::t(z$centered[[k]][rows, , drop = FALSE])
          solved <- Matrix::solve(z$factor, rhs, system = "A")
          error <- max(abs(z$H %*% solved - rhs)) / max(1, max(abs(rhs)))
          if (!is.finite(error) || error > 1e-8 || any(!is.finite(solved)))
            stop("Jeffreys information solve failed its residual check.", call. = FALSE)
          residual <- max(residual, error)
          q[rows] <- as.numeric(Matrix::colSums(rhs * solved))
        }
        adjustment <- adjustment + as.numeric(Matrix::crossprod(z$centered[[k]],
          z$probabilities[, k] * q)) / 2
      }
      z$ordinary_gradient <- g; z$adjustment <- adjustment
      z$gradient <- g - adjustment; z$solve_residual <- residual
      if (any(!is.finite(z$gradient))) stop("Nonfinite Jeffreys gradient.", call. = FALSE)
      cache$value <- z
    }
    z
  }
  list(evaluate = evaluate,
    fn = function(beta) tryCatch(evaluate(beta)$value, error = function(e) Inf),
    gr = function(beta) evaluate(beta, TRUE)$gradient)
}

mfrm_fair_jeffreys_design <- function(fit) {
  mfrm_fair_reference_fit(fit)
  cfg <- fit$config; prep <- fit$prep
  if (any(prep$data$Weight != 1) ||
      !identical(prep$rating_min_source, "declared") ||
      !identical(prep$rating_max_source, "declared") ||
      !identical(as.numeric(prep$score_map$OriginalScore), as.numeric(seq(prep$rating_min, prep$rating_max))) ||
      !identical(as.numeric(prep$score_map$InternalScore), as.numeric(seq(prep$rating_min, prep$rating_max))))
    stop("The Jeffreys generator requires unit weights and explicitly declared rubric limits with unchanged full score coding.", call. = FALSE)
  design <- mfrm_fair_delta_design(fit)
  idx <- design$idx; sizes <- build_param_sizes(cfg)
  zero <- numeric(length(fit$opt$par))
  pa <- expand_params(zero, sizes, cfg)
  eta <- pa$theta[idx$person] + compute_base_eta(idx, pa, cfg)
  steps <- if (cfg$model == "RSM") matrix(pa$steps, length(eta), cfg$n_cat - 1L, byrow = TRUE) else
    pa$steps_mat[idx$step_idx, , drop = FALSE]
  cumulative <- t(vapply(seq_len(nrow(steps)), function(j) c(0, cumsum(steps[j, ])), numeric(cfg$n_cat)))
  offset <- outer(eta, seq_len(cfg$n_cat) - 1L) - cumulative
  list(category = design$category, offset = offset, y = idx$score_k,
    parameter_map = design$parameter_map)
}

# Objective-change stopping may precede score-equation stationarity. Use only
# bounded local steps of the SAME penalized objective, with its exact Jacobian.
# This is not an alternative penalty, a relaxed gradient gate or a rank repair.
mfrm_fair_jeffreys_polish <- function(fit, par, objective, gradient_tol) {
  current <- objective$evaluate(par, TRUE)
  max_gradient <- function(x) if (length(x$gradient)) max(abs(x$gradient)) else 0
  diagnostics <- list(attempted = max_gradient(current) > gradient_tol, accepted_steps = 0L,
    initial_gradient = max_gradient(current), final_gradient = max_gradient(current),
    initial_objective = current$value, final_objective = current$value,
    max_step = 0, solve_residual = 0, reason = "not_needed")
  if (!diagnostics$attempted) return(list(par = par, diagnostics = diagnostics))
  design <- mfrm_fair_delta_design(fit)
  anchors <- data.frame(Facet = character(), Level = character(), Anchor = numeric())
  diagnostics$reason <- "step_limit"
  for (iteration in seq_len(3L)) {
    fit$opt$par <- par
    step <- tryCatch({
      J <- mfrm_fair_jeffreys_jacobian(fit, design, current$probabilities,
        anchors, mfrmr_empty_sparse_matrix(0L, 0L))
      mfrm_fair_delta_solve(J$information, matrix(current$gradient, ncol = 1L), "Jeffreys finishing")
    }, error = function(e) list(reason = conditionMessage(e)))
    if (!is.null(step$reason)) {diagnostics$reason <- step$reason; break}
    delta <- drop(step$solution)
    if (any(!is.finite(delta)) || max(abs(delta)) > 1) {
      diagnostics$reason <- "step_exceeds_local_bound"; break
    }
    accepted <- FALSE
    for (half in 0:8) {
      candidate_par <- par - delta * 2^(-half)
      candidate <- tryCatch(objective$evaluate(candidate_par, TRUE), error = function(e) NULL)
      if (!is.null(candidate) && is.finite(candidate$value) &&
          candidate$value <= current$value && max_gradient(candidate) < max_gradient(current) &&
          all(candidate$probabilities > 0 & candidate$probabilities < 1)) {
        accepted <- TRUE; break
      }
    }
    if (!accepted) {diagnostics$reason <- "no_decreasing_local_step"; break}
    diagnostics$accepted_steps <- diagnostics$accepted_steps + 1L
    diagnostics$max_step <- max(diagnostics$max_step, max(abs(candidate_par - par)))
    diagnostics$solve_residual <- max(diagnostics$solve_residual, step$residual)
    par <- candidate_par; current <- candidate
    if (max_gradient(current) <= gradient_tol) {diagnostics$reason <- "stationary"; break}
  }
  diagnostics$final_gradient <- max_gradient(current)
  diagnostics$final_objective <- current$value
  list(par = par, diagnostics = diagnostics)
}

mfrm_fair_jeffreys_generator <- function(fit, gradient_tol) {
  d <- mfrm_fair_jeffreys_design(fit)
  objective <- mfrm_fair_jeffreys_objective(d$category, d$offset, d$y)
  # Start at the canonical origin, not at an arbitrary extreme-score proxy.
  start <- numeric(ncol(d$category[[1]]))
  objective$evaluate(start, TRUE)
  started <- proc.time()[[3L]]
  reltol <- min(fit$config$replay_inputs$reltol, 1e-14)
  opt <- if (length(start)) stats::optim(start, objective$fn, objective$gr,
    method = "BFGS", control = list(maxit = fit$config$replay_inputs$maxit,
      reltol = reltol)) else
    list(par = start, convergence = 0L, counts = c("function" = 1L, gradient = 1L))
  final <- objective$evaluate(opt$par, TRUE)
  gradient <- if (length(final$gradient)) max(abs(final$gradient)) else 0
  polish <- NULL
  if (identical(as.integer(opt$convergence), 0L) && gradient > gradient_tol) {
    finished <- mfrm_fair_jeffreys_polish(fit, opt$par, objective, gradient_tol)
    opt$par <- finished$par; polish <- finished$diagnostics
    final <- objective$evaluate(opt$par, TRUE)
    gradient <- if (length(final$gradient)) max(abs(final$gradient)) else 0
  }
  if (!identical(as.integer(opt$convergence), 0L) || gradient > gradient_tol)
    stop("Jeffreys generating-model stationarity failed (maximum penalized gradient = ",
      signif(gradient, 5), ").", call. = FALSE)
  probs <- final$probabilities
  if (any(!is.finite(probs)) || any(probs <= 0 | probs >= 1) ||
      any(abs(rowSums(probs) - 1) > 1e-10))
    stop("Jeffreys generating probabilities are not numerically interior; no clipping was applied.", call. = FALSE)
  list(method = "jeffreys_canonical_joint", par = opt$par, probabilities = probs,
    parameter_map = d$parameter_map, objective = final$value, ordinary_nll = final$nll,
    information_logdet = final$logdet, max_gradient = gradient,
    ordinary_max_gradient = if (length(final$ordinary_gradient)) max(abs(final$ordinary_gradient)) else 0,
    solve_residual = final$solve_residual, convergence = opt$convergence, counts = opt$counts,
    polish = polish,
    seconds = proc.time()[[3L]] - started,
    control = list(optimizer = "BFGS", maxit = fit$config$replay_inputs$maxit,
      reltol = reltol, gradient_tol = gradient_tol),
    inference = "Canonical penalized point; no score-bias or coverage qualification.")
}

mfrm_fair_jeffreys_target <- function(fit, par, plan, parm) {
  cfg <- fit$config
  kernel <- mfrm_fair_reference_kernel(fit, plan, par)
  theta <- expand_params(par, build_param_sizes(cfg), cfg)$theta[match(parm, cfg$person_levels)]
  vapply(theta, function(t) {
    pr <- category_prob_pcm(t + kernel$offsets, kernel$cumulative, seq_len(nrow(plan)))
    sum(plan$NormalizedWeight * as.numeric(pr %*% (seq_len(cfg$n_cat) - 1L))) + fit$prep$rating_min
  }, 0)
}

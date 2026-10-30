# Exact ordinary additive JML curvature in the actual affine free coordinates.
# Used only to finish a stalled numerical solution, never to regularize a model.
mfrm_jml_polish_information <- function(par, idx, config, sizes) {
  if (!identical(config$method, "JML") ||
      !config$model %in% c("RSM", "PCM") || length(config$interaction_specs) ||
      !identical(config$facet_shrinkage %||% "none", "none"))
    stop("Information polishing requires ordinary additive RSM/PCM JML.")
  design <- mfrmr_estimability_adjacent_design(NULL, idx, config, sizes,
    include_person = TRUE, include_population_beta = FALSE)
  if (!identical(sort(design$map$OptimizerIndex), seq_along(par)))
    stop("JML information does not cover the optimizer coordinates.")
  order <- order(design$map$OptimizerIndex)
  n <- length(idx$score_k); p <- length(par)
  cache <- make_param_cache(sizes, config, idx, is_mml = FALSE)
  cache$ensure(par)
  probs <- mfrm_jml_probability_bundle(cache$eta(), idx$score_k,
    config$model, cache$step_cum(), criterion_idx = idx$step_idx)$probs
  if (any(!is.finite(probs)) || any(probs < 0) ||
      any(abs(rowSums(probs) - 1) > 1e-10)) stop("Invalid JML probabilities.")
  categories <- vector("list", ncol(probs))
  categories[[1]] <- mfrmr_empty_sparse_matrix(n, p)
  mean <- categories[[1]]
  for (k in seq_len(ncol(probs) - 1L)) {
    categories[[k + 1L]] <- categories[[k]] +
      design$design[seq_len(n) + (k - 1L) * n, order, drop = FALSE]
    mean <- mean + Matrix::Diagonal(x = probs[, k + 1L]) %*% categories[[k + 1L]]
  }
  weight <- idx$weight %||% rep(1, n)
  H <- mfrmr_empty_sparse_matrix(p, p)
  for (k in seq_along(categories)) {
    centered <- categories[[k]] - mean
    H <- H + Matrix::crossprod(centered,
      Matrix::Diagonal(x = weight * probs[, k]) %*% centered)
  }
  H
}

# At most three undamped Newton proposals, each checked on the original
# objective and gradient. A failed solve/proposal retains the preceding point.
# No finite replacement for an extreme Person is introduced here.
mfrm_jml_information_polish <- function(par, fn, gr, information,
                                        gradient_tol = 1e-4, active = seq_along(par)) {
  if (!is.numeric(active) || anyNA(active) || any(!is.finite(active)) ||
      any(active != floor(active)) || any(active < 1 | active > length(par)) || anyDuplicated(active))
    stop("Active JML coordinates must be distinct valid indices.")
  counts <- c("function" = 0L, "gradient" = 0L)
  value <- function(p) { counts[1L] <<- counts[1L] + 1L; fn(p) }
  gradient <- function(p) { counts[2L] <<- counts[2L] + 1L; gr(p) }
  steps <- 0L; residual <- NA_real_; error <- ""
  f <- value(par); g <- gradient(par)
  repeat {
    if (!length(par) || !all(is.finite(c(par, f, g))) || length(g) != length(par)) {
      error <- "JML polishing requires a finite objective and complete gradient."; break
    }
    if (max(abs(g)) <= gradient_tol || steps == 3L) break
    proposal <- tryCatch({
      if (!length(active)) stop("No active JML coordinates remain for finishing.")
      # A block update of the original objective: retain every observation and
      # test the COMPLETE gradient, including coordinates held at their input.
      H <- information(par)[active, active, drop = FALSE]
      diagonal <- as.numeric(Matrix::diag(H))
      if (any(!is.finite(diagonal)) || any(diagonal <= 0))
        stop("JML information has a nonpositive/nonfinite diagonal.")
      D <- Matrix::Diagonal(x = 1 / sqrt(diagonal))
      scaled <- Matrix::forceSymmetric(D %*% H %*% D)
      factor <- Matrix::Cholesky(scaled, LDL = FALSE, perm = TRUE)
      rhs <- D %*% g[active]
      z <- Matrix::solve(factor, rhs, system = "A")
      residual <- max(abs(scaled %*% z - rhs)) / max(1, max(abs(rhs)))
      delta <- as.numeric(D %*% z)
      if (!is.finite(residual) || residual > 1e-8 || any(!is.finite(delta)))
        stop("JML information solve failed its residual check.")
      # This is a local finishing step, not an unbounded search down a recession.
      if (max(abs(delta)) > 1) stop("JML information proposal is not local (step exceeds one logit).")
      candidate <- par
      candidate[active] <- par[active] - delta
      next_f <- value(candidate); next_g <- gradient(candidate)
      slack <- 64 * .Machine$double.eps * max(1, abs(f))
      if (!all(is.finite(c(candidate, next_f, next_g))) ||
          length(next_g) != length(g) || next_f > f + slack ||
          max(abs(next_g)) >= max(abs(g)))
        stop("JML information proposal did not improve stationarity without likelihood loss.")
      list(par = candidate, value = next_f, gradient = next_g)
    }, error = identity, warning = function(w) w)
    if (inherits(proposal, "condition")) { error <- conditionMessage(proposal); break }
    par <- proposal$par; f <- proposal$value; g <- proposal$gradient
    steps <- steps + 1L
  }
  list(par = par, value = f, gradient = g, counts = counts, steps = steps,
    residual = residual, error = error, fixed_coordinates = setdiff(seq_along(par), active),
    convergence = if (length(par) && length(g) == length(par) &&
      all(is.finite(c(par, f, g))) && max(abs(g)) <= gradient_tol) 0L else 1L)
}

# ==============================================================================
# Polytomous Rasch likelihoods and response-probability bundles
# ==============================================================================
#
# Pure-math helpers for the Rasch-family polytomous likelihoods used by
# the JML and MML estimators in `mfrm_core.R`. Split out of `mfrm_core.R`
# for 0.1.6 so the per-observation likelihood kernels and the cumulative
# response-probability helpers live in a single file. The functions are
# internal (no @export); they are called by the JML/MML likelihood
# builders and by the information / fit-statistic helpers.

# ---- likelihoods ----
# RSM log-likelihood: sum_i w_i log P(X_i = k_i | eta_i).
# Under the Rating Scale Model (Andrich, 1978):
#   P(X = k | eta) = exp(k*eta - tau_k) / sum_j exp(j*eta - tau_j)
# where tau_k = cumulative step parameters and eta = theta - sum(facets).
# Computation uses log-domain subtraction with logsumexp for stability.
loglik_rsm <- function(eta, score_k, step_cum, weight = NULL) {
  n <- length(eta)
  if (n == 0) return(0)
  k_cat <- length(step_cum)
  # Defensive guard: `m[cbind(i, score_k + 1)]` silently drops rows with
  # `score_k < 0` or `score_k >= k_cat`. prepare_mfrm_data() should catch
  # this earlier, but fail loudly here in case callers build score_k manually.
  if (any(!is.finite(score_k) | score_k < 0L | score_k >= k_cat)) {
    stop("`score_k` has values outside [0, ", k_cat - 1L,
         "]. This usually indicates `Score` fell outside `rating_min:rating_max`.",
         call. = FALSE)
  }
  eta_mat <- outer(eta, 0:(k_cat - 1))
  log_num <- eta_mat - matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
  row_max <- row_max_fast(log_num)
  log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
  log_num_obs <- log_num[cbind(seq_len(n), score_k + 1)]
  diff <- log_num_obs - log_denom
  if (is.null(weight)) {
    sum(diff)
  } else {
    sum(diff * weight)
  }
}

# PCM log-likelihood: same structure as RSM but with criterion-specific steps.
# Under the Partial Credit Model (Masters, 1982):
#   P(X = k | eta, criterion c) = exp(k*eta - tau_{c,k}) / sum_j exp(j*eta - tau_{c,j})
# step_cum_mat has one row per criterion level, columns = cumulative thresholds.
loglik_pcm <- function(eta, score_k, step_cum_mat, criterion_idx, weight = NULL,
                       criterion_splits = NULL) {
  n <- length(eta)
  if (n == 0) return(0)
  k_cat <- ncol(step_cum_mat)
  if (any(!is.finite(score_k) | score_k < 0L | score_k >= k_cat)) {
    stop("`score_k` has values outside [0, ", k_cat - 1L,
         "]. This usually indicates `Score` fell outside `rating_min:rating_max`.",
         call. = FALSE)
  }
  total <- 0
  splits <- criterion_splits %||% split(seq_len(n), criterion_idx)
  for (ci in seq_along(splits)) {
    rows <- splits[[ci]]
    if (length(rows) == 0) next
    c_idx <- as.integer(names(splits)[ci])
    eta_c <- eta[rows]
    step_cum <- step_cum_mat[c_idx, ]
    nr <- length(rows)
    eta_mat <- outer(eta_c, 0:(k_cat - 1))
    log_num <- eta_mat - matrix(step_cum, nrow = nr, ncol = k_cat, byrow = TRUE)
    row_max <- log_num[cbind(seq_len(nr), max.col(log_num))]
    log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
    log_num_obs <- log_num[cbind(seq_len(nr), score_k[rows] + 1)]
    diff <- log_num_obs - log_denom
    if (is.null(weight)) {
      total <- total + sum(diff)
    } else {
      total <- total + sum(diff * weight[rows])
    }
  }
  total
}

# GPCM log-likelihood: same adjacent-category structure as PCM but with a
# positive discrimination attached to each designated slope-facet level.
# Under the first-release target:
#   log(P_k / P_{k-1}) = a_c * (eta - tau_{c,k})
# so category k has kernel exp(a_c * (k * eta - tau_{c,k}^{cum})).
loglik_gpcm <- function(eta, score_k, step_cum_mat, criterion_idx, slopes,
                        slope_idx = criterion_idx, weight = NULL) {
  n <- length(eta)
  if (n == 0) return(0)
  k_cat <- ncol(step_cum_mat)
  if (any(!is.finite(score_k) | score_k < 0L | score_k >= k_cat)) {
    stop("`score_k` has values outside [0, ", k_cat - 1L,
         "]. This usually indicates `Score` fell outside `rating_min:rating_max`.",
         call. = FALSE)
  }
  if (length(criterion_idx) != n || length(slope_idx) != n) {
    stop("`criterion_idx` and `slope_idx` must have one entry per observation.",
         call. = FALSE)
  }

  criterion_idx <- as.integer(criterion_idx)
  slope_idx <- as.integer(slope_idx)
  if (any(!is.finite(criterion_idx)) || any(criterion_idx < 1L) ||
      any(criterion_idx > nrow(step_cum_mat))) {
    stop("`criterion_idx` must index valid rows of `step_cum_mat`.", call. = FALSE)
  }
  if (any(!is.finite(slope_idx)) || any(slope_idx < 1L) ||
      any(slope_idx > length(slopes))) {
    stop("`slope_idx` must index valid `slopes` entries.", call. = FALSE)
  }

  slope_obs <- as.numeric(slopes[slope_idx])
  if (any(!is.finite(slope_obs)) || any(slope_obs <= 0)) {
    stop("Observed GPCM slopes must be finite and strictly positive.", call. = FALSE)
  }

  step_cum_obs <- step_cum_mat[criterion_idx, , drop = FALSE]
  k_cat <- ncol(step_cum_obs)
  linear_part <- outer(eta, 0:(k_cat - 1)) - step_cum_obs
  log_num <- linear_part * matrix(slope_obs, nrow = n, ncol = k_cat)
  row_max <- row_max_fast(log_num)
  log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
  log_num_obs <- log_num[cbind(seq_len(n), score_k + 1)]
  diff <- log_num_obs - log_denom
  if (is.null(weight)) {
    sum(diff)
  } else {
    sum(diff * weight)
  }
}

# Category response probabilities under RSM.
# Returns an n x K matrix where K = number of categories.
# Each row sums to 1; probabilities are computed in log-domain for stability.
# category_prob_rsm / category_prob_pcm / category_prob_gpcm are
# defined in R/core-category-probabilities.R (split out of this file
# in 0.1.6 for clarity).

# Compute P(X >= s) matrix for s = 1,...,K-1 from category probabilities.
# Input: probs (n x K matrix of category probabilities, columns for k=0,...,K-1)
# Output: P_geq (n x (K-1) matrix), P_geq[i,s] = P(X_i >= s)
compute_P_geq_r <- function(probs) {
  k_cat <- ncol(probs)
  n_steps <- k_cat - 1
  if (n_steps == 0) return(matrix(0, nrow(probs), 0))
  n <- nrow(probs)
  P_geq <- matrix(0, n, n_steps)
  P_geq[, n_steps] <- probs[, k_cat]
  if (n_steps >= 2) {
    for (s in (n_steps - 1):1) {
      P_geq[, s] <- P_geq[, s + 1] + probs[, s + 1]
    }
  }
  P_geq
}

compute_P_geq <- function(probs) {
  if (isTRUE(getOption("mfrmr.use_cpp11_backend", FALSE)) &&
      mfrm_cpp11_backend_available() &&
      is.matrix(probs) &&
      is.numeric(probs)) {
    return(mfrm_cpp_compute_p_geq(probs))
  }
  compute_P_geq_r(probs)
}

compute_response_probability_bundle <- function(config, idx, params, eta) {
  n_obs <- length(eta)
  if (n_obs == 0L) {
    return(list(
      probs = matrix(0, nrow = 0L, ncol = max(config$n_cat %||% 0L, 0L)),
      expected_k = numeric(0),
      var_k = numeric(0),
      fourth_central_moment = numeric(0),
      score_information = numeric(0),
      slope_obs = numeric(0)
    ))
  }

  if (identical(config$model, "RSM")) {
    step_cum <- c(0, cumsum(params$steps))
    probs <- category_prob_rsm(eta, step_cum)
    slope_obs <- rep(1, n_obs)
  } else if (identical(config$model, "GPCM")) {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    slope_idx <- idx$slope_idx %||% idx$step_idx
    if (is.null(slope_idx)) {
      stop("GPCM response probabilities require a valid slope index.", call. = FALSE)
    }
    probs <- category_prob_gpcm(
      eta = eta,
      step_cum_mat = step_cum_mat,
      criterion_idx = idx$step_idx,
      slopes = params$slopes,
      slope_idx = slope_idx
    )
    slope_obs <- as.numeric(params$slopes[slope_idx])
  } else {
    step_cum_mat <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    probs <- category_prob_pcm(
      eta = eta,
      step_cum_mat = step_cum_mat,
      criterion_idx = idx$step_idx,
      criterion_splits = idx$criterion_splits
    )
    slope_obs <- rep(1, n_obs)
  }

  k_vals <- 0:(ncol(probs) - 1L)
  expected_k <- as.vector(probs %*% k_vals)
  var_k <- as.vector(probs %*% (k_vals^2)) - expected_k^2
  diff_k <- sweep(
    matrix(k_vals, nrow = nrow(probs), ncol = length(k_vals), byrow = TRUE),
    1L,
    expected_k,
    FUN = "-"
  )
  fourth_central_moment <- as.vector(rowSums(probs * diff_k^4))
  var_k <- ifelse(var_k <= 1e-10, NA_real_, var_k)
  # For bounded GPCM, the score information with respect to eta is
  # a^2 Var(X | eta); PCM/RSM are the a = 1 special case.
  score_information <- ifelse(
    is.finite(var_k) & is.finite(slope_obs),
    (slope_obs^2) * var_k,
    NA_real_
  )

  list(
    probs = probs,
    expected_k = expected_k,
    var_k = var_k,
    fourth_central_moment = fourth_central_moment,
    score_information = score_information,
    slope_obs = slope_obs
  )
}

# Convert mean-square fit statistic to a standardized z-score (ZSTD).
# Default uses the Wilson-Hilferty (1931) cube-root approximation:
#   ZSTD = (MnSq^(1/3) - (1 - 2/(9*df))) / sqrt(2/(9*df))
# When whexact = TRUE, uses the simpler linear approximation:
#   ZSTD = (MnSq - 1) * sqrt(df / 2)
# Values near 0 indicate expected fit; |ZSTD| > 2 flags potential misfit.
zstd_from_mnsq <- function(mnsq, df, whexact = FALSE) {
  mnsq <- as.numeric(mnsq)
  df <- as.numeric(df)

  if (length(df) == 1L && length(mnsq) > 1L) {
    df <- rep(df, length(mnsq))
  } else if (length(mnsq) == 1L && length(df) > 1L) {
    mnsq <- rep(mnsq, length(df))
  }

  n <- min(length(mnsq), length(df))
  if (n == 0L) return(numeric(0))

  out <- rep(NA_real_, n)
  m <- mnsq[seq_len(n)]
  d <- df[seq_len(n)]
  # Wilson-Hilferty approximation becomes numerically unstable and
  # sign-inverted for very small df: when d < 1 the term 2/(9*d)
  # dominates and MnSq ~ 0 can return large positive ZSTD.
  # Guard with d >= 1 and fall back to NA for ill-conditioned cells.
  ok <- is.finite(m) & is.finite(d) & (d >= 1)

  if (isTRUE(whexact)) {
    out[ok] <- (m[ok] - 1) * sqrt(d[ok] / 2)
  } else {
    out[ok] <- (m[ok]^(1 / 3) - (1 - 2 / (9 * d[ok]))) / sqrt(2 / (9 * d[ok]))
  }
  out
}

# FACETS/Winsteps-style standardized fit values use the same mean-square
# transform but a separate Wright-Masters df convention and a reported cap.
zstd_from_mnsq_facets <- function(mnsq, df, whexact = FALSE, cap = 9) {
  mnsq <- as.numeric(mnsq)
  df <- as.numeric(df)

  if (length(df) == 1L && length(mnsq) > 1L) {
    df <- rep(df, length(mnsq))
  } else if (length(mnsq) == 1L && length(df) > 1L) {
    mnsq <- rep(mnsq, length(df))
  }

  n <- min(length(mnsq), length(df))
  if (n == 0L) return(numeric(0))

  out <- rep(NA_real_, n)
  m <- mnsq[seq_len(n)]
  d <- df[seq_len(n)]
  ok <- is.finite(m) & is.finite(d) & (m > 0) & (d > 0)

  if (isTRUE(whexact)) {
    out[ok] <- (m[ok] - 1) * sqrt(d[ok] / 2)
  } else {
    out[ok] <- (m[ok]^(1 / 3) - (1 - 2 / (9 * d[ok]))) / sqrt(2 / (9 * d[ok]))
  }

  cap <- suppressWarnings(as.numeric(cap[1] %||% NA_real_))
  if (is.finite(cap) && cap > 0) {
    out <- pmin(pmax(out, -cap), cap)
  }
  out
}

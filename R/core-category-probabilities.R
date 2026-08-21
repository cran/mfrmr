# ==============================================================================
# Category response probabilities (RSM / PCM / GPCM)
# ==============================================================================
#
# Pure-math helpers for the Rasch-family polytomous response models.
# Split out of `mfrm_core.R` so the category probability
# kernels live in a single, browseable file. The functions are
# internal (no @export) and are called directly by the MML / JML
# likelihood builders, the information helpers, and the pathway /
# CCC plot helpers.
#
# All three helpers return an `n x K` matrix of P(X = k - 1 | ...)
# where K is the number of ordered categories. Each row sums to 1.
# The LogSumExp form is used so inputs on an extreme logit scale do
# not overflow exp().
# ==============================================================================


# Row-wise maximum with tie-break at the first column. Kept in this
# file because `category_prob_*` are the only callers that need it
# outside `mfrm_core.R`. It remains in the mfrmr namespace and is
# reachable from other files without modification.
row_max_fast <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) return(numeric(0))
  mat[cbind(seq_len(nrow(mat)), max.col(mat, ties.method = "first"))]
}

# Category response probabilities under the Rating Scale Model (RSM;
# Andrich, 1978). `eta` is the per-observation linear predictor minus
# the Rasch-Andrich cumulative thresholds (`step_cum`, of length
# `K = n_cat`). Shared thresholds across all items.
category_prob_rsm <- function(eta, step_cum) {
  n <- length(eta)
  if (n == 0) return(matrix(0, nrow = 0, ncol = length(step_cum)))
  k_cat <- length(step_cum)
  eta_mat <- outer(eta, 0:(k_cat - 1))
  log_num <- eta_mat - matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
  row_max <- row_max_fast(log_num)
  log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
  exp(log_num - matrix(log_denom, nrow = n, ncol = k_cat))
}

# Category response probabilities under the Partial Credit Model
# (PCM; Masters, 1982). `step_cum_mat` is a `K x n_cat` matrix of
# step-facet-level-specific cumulative thresholds; `criterion_idx`
# gives the row index for each observation. `criterion_splits` can be
# supplied by the caller to avoid re-splitting on every call.
category_prob_pcm <- function(eta, step_cum_mat, criterion_idx,
                              criterion_splits = NULL) {
  n <- length(eta)
  if (n == 0) return(matrix(0, nrow = 0, ncol = ncol(step_cum_mat)))
  k_cat <- ncol(step_cum_mat)
  probs <- matrix(0, nrow = n, ncol = k_cat)
  splits <- criterion_splits %||% split(seq_len(n), criterion_idx)
  for (ci in seq_along(splits)) {
    rows <- splits[[ci]]
    if (length(rows) == 0) next
    c_idx <- as.integer(names(splits)[ci])
    step_cum <- step_cum_mat[c_idx, ]
    eta_c <- eta[rows]
    eta_mat <- outer(eta_c, 0:(k_cat - 1))
    log_num <- eta_mat - matrix(step_cum, nrow = length(rows), ncol = k_cat, byrow = TRUE)
    row_max <- row_max_fast(log_num)
    log_denom <- row_max + log(rowSums(exp(log_num - row_max)))
    probs[rows, ] <- exp(log_num - matrix(log_denom, nrow = length(rows), ncol = k_cat))
  }
  probs
}

# Category response probabilities under the bounded Generalized
# Partial Credit Model (GPCM; Muraki, 1992). Discriminations
# (`slopes`) must be strictly positive; the internal convention is
# sum-to-zero on the log-slope scale with geometric-mean-one
# identification. `eta` already contains the person coordinate and every
# additive facet/interaction contribution. The slope multiplies that complete
# predictor and the owned cumulative step; it is not a trait-loading-only
# parameterization with unscaled facet intercepts.
category_prob_gpcm <- function(eta, step_cum_mat, criterion_idx, slopes,
                               slope_idx = criterion_idx) {
  n <- length(eta)
  if (n == 0) return(matrix(0, nrow = 0, ncol = ncol(step_cum_mat)))
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
  exp(log_num - matrix(log_denom, nrow = n, ncol = k_cat))
}

# Joint probability / observed-log-probability kernel used by the optimizer.
# Keeping both quantities in one bundle lets the objective and analytical
# gradient share the same normalization work at an identical parameter vector.
# The observed log probabilities are retained directly rather than recovered
# with log(probs), which preserves log-domain stability for extreme logits.
mfrm_jml_probability_bundle <- function(eta,
                                          score_k,
                                          model,
                                          step_cum,
                                          criterion_idx = NULL,
                                          slopes = NULL,
                                          slope_idx = criterion_idx) {
  n <- length(eta)
  k_cat <- if (identical(model, "RSM")) length(step_cum) else ncol(step_cum)
  if (n == 0L) {
    return(list(
      probs = matrix(0, nrow = 0L, ncol = k_cat),
      log_prob_obs = numeric(0),
      linear_part = NULL
    ))
  }
  if (any(!is.finite(score_k) | score_k < 0L | score_k >= k_cat)) {
    stop("`score_k` has values outside [0, ", k_cat - 1L,
         "]. This usually indicates `Score` fell outside `rating_min:rating_max`.",
         call. = FALSE)
  }

  k_vals <- 0:(k_cat - 1L)
  linear_part <- if (identical(model, "RSM")) {
    outer(eta, k_vals) -
      matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
  } else {
    if (length(criterion_idx) != n ||
        any(!is.finite(criterion_idx)) || any(criterion_idx < 1L) ||
        any(criterion_idx > nrow(step_cum))) {
      stop("`criterion_idx` must index a valid step row for every observation.",
           call. = FALSE)
    }
    outer(eta, k_vals) - step_cum[as.integer(criterion_idx), , drop = FALSE]
  }

  log_num <- linear_part
  if (identical(model, "GPCM")) {
    if (length(slope_idx) != n ||
        any(!is.finite(slope_idx)) || any(slope_idx < 1L) ||
        any(slope_idx > length(slopes))) {
      stop("`slope_idx` must index a valid slope for every observation.",
           call. = FALSE)
    }
    slope_obs <- as.numeric(slopes[as.integer(slope_idx)])
    if (any(!is.finite(slope_obs)) || any(slope_obs <= 0)) {
      stop("Observed GPCM slopes must be finite and strictly positive.",
           call. = FALSE)
    }
    log_num <- linear_part * slope_obs
  }

  # Reuse the normalized log-kernel object rather than retaining separate
  # shifted, denominator-expanded, and log-probability matrices.
  log_prob <- log_num - row_max_fast(log_num)
  log_prob <- log_prob - log(rowSums(exp(log_prob)))

  list(
    probs = exp(log_prob),
    log_prob_obs = log_prob[cbind(seq_len(n), score_k + 1L)],
    linear_part = if (identical(model, "GPCM")) linear_part else NULL
  )
}

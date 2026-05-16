# ==============================================================================
# Category response probabilities (RSM / PCM / GPCM)
# ==============================================================================
#
# Pure-math helpers for the Rasch-family polytomous response models.
# Split out of `mfrm_core.R` for 0.1.6 so the category probability
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
# identification.
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

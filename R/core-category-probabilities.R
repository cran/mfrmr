row_max_fast <- function(mat) {
    if (nrow(mat) == 0 || ncol(mat) == 0)
        return(numeric(0))
    mat[cbind(seq_len(nrow(mat)), max.col(mat, ties.method = "first"))]
}

mfrm_jml_probability_bundle <- function(eta, score_k, model, step_cum, criterion_idx = NULL, slopes = NULL,
    slope_idx = criterion_idx, include_log_probs = FALSE) {
    n <- length(eta)
    k_cat <- if (identical(model, "RSM"))
        length(step_cum)
    else ncol(step_cum)
    if (n == 0L) {
        return(list(probs = matrix(0, nrow = 0L, ncol = k_cat), log_probs = if (isTRUE(include_log_probs)) matrix(0,
            nrow = 0L, ncol = k_cat) else NULL, log_prob_obs = numeric(0), linear_part = NULL))
    }
    if (any(!is.finite(score_k) | score_k < 0L | score_k >= k_cat)) {
        stop("`score_k` has values outside [0, ", k_cat - 1L, "]. This usually indicates `Score` fell outside `rating_min:rating_max`.",
            call. = FALSE)
    }
    k_vals <- 0:(k_cat - 1L)
    linear_part <- if (identical(model, "RSM")) {
        outer(eta, k_vals) - matrix(step_cum, nrow = n, ncol = k_cat, byrow = TRUE)
    }
    else {
        if (length(criterion_idx) != n || any(!is.finite(criterion_idx)) || any(criterion_idx < 1L) ||
            any(criterion_idx > nrow(step_cum))) {
            stop("`criterion_idx` must index a valid step row for every observation.", call. = FALSE)
        }
        outer(eta, k_vals) - step_cum[as.integer(criterion_idx), , drop = FALSE]
    }
    log_num <- linear_part
    if (identical(model, "GPCM")) {
        if (length(slope_idx) != n || any(!is.finite(slope_idx)) || any(slope_idx < 1L) || any(slope_idx >
            length(slopes))) {
            stop("`slope_idx` must index a valid slope for every observation.", call. = FALSE)
        }
        slope_obs <- as.numeric(slopes[as.integer(slope_idx)])
        if (any(!is.finite(slope_obs)) || any(slope_obs <= 0)) {
            stop("Observed GPCM slopes must be finite and strictly positive.", call. = FALSE)
        }
        log_num <- linear_part * slope_obs
    }
    log_prob <- log_num - row_max_fast(log_num)
    log_prob <- log_prob - log(rowSums(exp(log_prob)))
    list(probs = exp(log_prob), log_probs = if (isTRUE(include_log_probs)) log_prob else NULL, log_prob_obs = log_prob[cbind(seq_len(n),
        score_k + 1L)], linear_part = if (identical(model, "GPCM")) linear_part else NULL)
}

# Conditional total-score extrapolation. Original parameters and observations
# are never modified; this is not an MLE for an observed extreme total.
.mfrm_unanchored_extremes <- function(res, obs) {
  out <- get_extreme_levels(obs[get_weights(obs) > 0, , drop = FALSE],
    c("Person", res$config$facet_names), res$prep$rating_min, res$prep$rating_max)
  for (f in names(out)) out[[f]] <- out[[f]][
    facet_anchor_status(f, out[[f]], res$config) != "A"]
  out
}

.mfrm_extreme_report <- function(res, obs, facet, levels, anchors, fraction, params) {
  n <- length(levels)
  out <- list(measure = rep(NA_real_, n), se = rep(NA_real_, n),
    adjustment = numeric(n), count = numeric(n))
  if (fraction == 0) return(out)
  config <- res$config; prep <- res$prep
  facets <- c("Person", config$facet_names)
  # A fixed finite anchor remains informative even if its observed scores
  # happen to be all minimum/maximum; do not remove it as an infinite MLE.
  extreme <- .mfrm_unanchored_extremes(res, obs)
  keep <- rep(TRUE, nrow(obs))
  for (f in setdiff(facets, facet))
    keep <- keep & !as.character(obs[[f]]) %in% extreme[[f]]
  offset <- numeric(nrow(obs))
  if (facet != "Person") offset <- res$facets$person$Estimate[
    match(as.character(obs$Person), as.character(res$facets$person$Person))]
  for (f in setdiff(config$facet_names, facet)) {
    sign <- config$facet_signs[[f]] %||% -1
    offset <- offset + sign * params$facets[[f]][match(as.character(obs[[f]]), prep$levels[[f]])]
  }
  direction <- if (facet == "Person") 1 else config$facet_signs[[facet]] %||% -1
  weight <- get_weights(obs)
  if (config$model == "RSM") steps <- c(0, cumsum(params$steps)) else {
    steps <- t(apply(params$steps_mat, 1, function(x) c(0, cumsum(x))))
    owner <- match(as.character(obs[[config$step_facet]]), prep$levels[[config$step_facet]])
  }
  for (j in seq_len(n)) {
    if (anchors[j] == "A" || !levels[j] %in% extreme[[facet]]) next
    rows <- which(as.character(obs[[facet]]) == levels[j] & keep & weight > 0)
    if (!length(rows) || any(!is.finite(offset[rows]))) next
    w <- weight[rows]; correction <- fraction * min(w)
    span <- prep$rating_max - prep$rating_min
    target <- if (all(obs$Observed[rows] == prep$rating_min)) correction else sum(w) * span - correction
    if (!is.finite(target) || target <= 0 || target >= span * sum(w)) next
    moments <- function(eta) {
      p <- if (config$model == "RSM") category_prob_rsm(offset[rows] + eta, steps) else
        category_prob_pcm(offset[rows] + eta, steps, owner[rows])
      k <- seq_len(ncol(p)) - 1L
      mean <- drop(p %*% k)
      variance <- rowSums(p * (matrix(k, nrow(p), length(k), byrow = TRUE) - mean)^2)
      c(score = sum(w * mean), information = sum(w * variance))
    }
    f <- function(eta) moments(eta)[1L] - target
    bound <- max(16, max(abs(offset[rows])) + max(abs(steps)) + 4)
    while (bound < 1e6 && (f(-bound) >= 0 || f(bound) <= 0)) bound <- bound * 2
    if (!is.finite(bound) || f(-bound) >= 0 || f(bound) <= 0) next
    eta <- stats::uniroot(f, c(-bound, bound), tol = 1e-11)$root
    info <- moments(eta)[2L]
    if (!is.finite(info) || info <= 0) next
    out$measure[j] <- eta / direction
    out$se[j] <- 1 / sqrt(info)
    out$adjustment[j] <- correction
    out$count[j] <- sum(w)
  }
  out
}

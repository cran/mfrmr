# Fixed assignment: remove facet incidence, then diagonalize Person incidence.
.mfrm_mvgt_spectral_ratio <- function(g, weights, level) {
  block <- .mfrm_mvgt_residual_block(g)
  ans <- list(log_bounds = rep(NA_real_, 2L), status = block$status,
    residual_projection = block, valid = FALSE)
  if (!block$valid) return(ans)
  n <- nrow(g$data)
  indicator <- function(id) {
    index <- match(id, unique(id))
    z <- matrix(0, n, max(index)); z[cbind(seq_len(n), index)] <- 1
    z
  }
  person <- indicator(g$data[[g$design$person]])
  facet <- indicator(g$data[[unname(g$design$facets)]])
  projected <- qr.resid(qr(facet), person)
  decomposition <- svd(projected, nv = 0L)
  rank <- sum(decomposition$d > max(decomposition$d)*1e-10)
  ans$person_df <- rank; ans$residual_df <- block$df
  ans$facet_rank <- ncol(facet)
  if (rank == 0L || rank != n-ncol(facet)-block$df) {
    ans$status <- "No resolved Person contrasts after removing the facet"; return(ans)
  }
  u <- decomposition$u[,seq_len(rank),drop = FALSE]
  lambda <- decomposition$d[seq_len(rank)]^2
  ans$eigenvalues <- lambda
  ans$orthogonality_error <- max(abs(crossprod(facet, u)))
  if (ans$orthogonality_error > 1e-10) {
    ans$status <- "Person projection failed numerical checks"; return(ans)
  }
  a <- weights*block$score_scale; scale <- max(abs(a))
  if (!is.finite(scale) || scale == 0)
    stop("Composite units exceed numeric range; rescale scores or weights.", call. = FALSE)
  a <- a/scale
  y <- sweep(as.matrix(g$data[g$design$scores]), 2L, block$score_scale, "/")
  y <- drop(sweep(y, 2L, colMeans(y), "-") %*% a)
  z <- drop(crossprod(u, y))
  residual_ms <- drop(crossprod(a, block$mean_product %*% a))
  energy <- sum(abs(a)*sqrt(block$score_energy))^2
  ans$person_squares <- z*z; ans$residual_ms <- residual_ms
  ans$normalized_score_scale <- scale
  if (residual_ms*block$df <= block$relative_resolution*energy) {
    ans$status <- "Residual information at numerical resolution"; return(ans)
  }
  log_terms <- 2*log(abs(z))
  log_sum <- function(x) {
    m <- max(x)
    if (m == -Inf) -Inf else m+log(sum(exp(x-m)))
  }
  log_f <- function(log_ratio) {
    v <- log(lambda)+log_ratio
    softplus <- pmax(v, 0)+log1p(exp(-abs(v)))
    log_sum(log_terms-softplus)-log(rank)-log(residual_ms)
  }
  log_f0 <- log_f(-Inf)
  tail <- (1-level)/2
  critical <- log(c(stats::qf(tail, rank, block$df, lower.tail = FALSE),
    stats::qf(tail, rank, block$df)))
  ans$log_f_zero <- log_f0; ans$log_critical <- critical
  ans$status <- "Numerical F quantile failure"
  if (any(!is.finite(critical))) return(ans)
  tolerance <- min(1e-8, 1e-5*(critical[1L]-critical[2L]))
  resolution <- 64*.Machine$double.eps*max(1, abs(critical))
  ans$endpoint_tolerance <- tolerance; ans$quantile_resolution <- resolution
  if (tolerance < resolution) {
    ans$status <- "F quantile interval below numeric resolution"; return(ans)
  }
  if (log_f0 < critical[2L]) {
    ans$status <- "Empty parameter-space interval"; return(ans)
  }
  # F(0)/(1+lambda_max*t) <= F(t) <= F(0)/(1+lambda_min*t).
  # These analytic root brackets avoid an arbitrary upper search limit.
  endpoint <- function(target) {
    if (log_f0 <= target) return(-Inf)
    v <- log_f0-target
    bracket <- v+log(-expm1(-v))-log(range(lambda)[2:1])
    if (diff(bracket) < 1e-12) return(mean(bracket))
    stats::uniroot(function(x) log_f(x)-target, bracket,
      tol = min(1e-11, tolerance), check.conv = TRUE)$root
  }
  limits <- tryCatch(vapply(critical, endpoint, numeric(1)), error = function(e) NULL)
  ans$status <- "Spectral F endpoint calculation failed"
  if (is.null(limits)) return(ans)
  finite <- is.finite(limits)
  error <- if (any(finite)) max(abs(vapply(limits[finite], log_f, numeric(1))-critical[finite])) else 0
  ans$endpoint_error <- error
  if (anyNA(limits) || limits[1L] > limits[2L] || error > tolerance) return(ans)
  if (is.finite(limits[2L]) && limits[1L] == limits[2L]) {
    ans$status <- "Variance-ratio interval below numeric resolution"; return(ans)
  }
  ans$log_bounds <- limits; ans$valid <- TRUE
  ans$status <- if (limits[2L] == -Inf) "Upper endpoint at zero" else
    if (limits[1L] == -Inf) "Includes zero" else "Available"
  ans
}

.mfrm_mvgt_spectral_table <- function(d, tab, weights, level, reference) {
  g <- d$gstudy
  pivot <- .mfrm_mvgt_spectral_ratio(g, weights, level)
  counts <- d$design_grid[[1L]]
  out <- tab[c("Scenario", names(d$design_grid), "Kind", "Score")]
  out$Metric <- "G"; out$Estimate <- tab$G; out$PointStatus <- tab$GStatus
  out$Lower <- stats::plogis(log(counts)+pivot$log_bounds[1L])
  out$Upper <- stats::plogis(log(counts)+pivot$log_bounds[2L])
  out$LogVarianceRatioLower <- pivot$log_bounds[1L]
  out$LogVarianceRatioUpper <- pivot$log_bounds[2L]
  out$Level <- level; out$Method <- "Exact spectral F under normal effects"
  out$Status <- pivot$status
  out$Validation <- if (pivot$valid) "Exact normal-model coverage for positive residual variance" else "No interval supplied"
  if (pivot$valid && is.null(reference) && pivot$log_bounds[2L] > pivot$log_bounds[1L]) {
    unresolved <- out$Lower == out$Upper
    out$Lower[unresolved] <- out$Upper[unresolved] <- NA_real_
    out$Status[unresolved] <- "G interval below numeric resolution"
    out$Validation[unresolved] <- "No interval supplied"
  }
  if (!is.null(reference)) {
    bounds <- .mfrm_mvgt_count_difference_range(pivot$log_bounds, counts, counts[reference], pivot$status)
    out$Reference <- as.integer(reference); out$ReferenceCount <- counts[reference]
    out$ReferenceValue <- tab$G[reference]; out$Value <- tab$G
    out$Metric <- "G difference"; out$Estimate <- out$Value-out$ReferenceValue
    out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
    out$Lower <- bounds$lower; out$Upper <- bounds$upper; out$Status <- bounds$status
    out$Method <- "Spectral F confidence-set image"
    if (pivot$valid) out$Validation <- "At least nominal normal-model coverage; set image may be conservative"
    attr(out, "reference") <- list(scenario = as.integer(reference), count = counts[reference],
      value = tab$G[reference], point_status = tab$GStatus[reference])
    out <- out[-reference, , drop = FALSE]
  }
  out$Validation[!is.finite(out$Lower) | !is.finite(out$Upper)] <- "No interval supplied"
  rownames(out) <- NULL
  attr(out, "weights") <- weights; attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- g$design[c("facets", "counts", "method", "balanced", "complete", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "source_data") <- g$data
  attr(out, "pivot") <- pivot
  attr(out, "notes") <- c(
    "One-facet G intervals condition on the observed assignment and remove the entire facet incidence. Independent normal Person effects, independent homoskedastic normal residuals and positive residual variance are required. Missingness or assignment must preserve these assumptions. Future plans remain complete balanced ratings.",
    "The spectral F pivot is a specialization of the variance-ratio family in E, Hannig and Iyer (2008), equations 3-7. It separates all positive Person eigenvalues from pure residual directions. It is not their fiducial method or the optimal Burch-Iyer partition. No likelihood fitting or simulation is used.",
    "Removing facet effects leaves the population G target unchanged but discards directions that may carry information. Exact normal-model coverage does not imply a short, most efficient or distribution-robust interval. Balanced one-facet limits agree with the usual exact F method.",
    "Estimate and PointStatus retain the original ANOVA/MINQUE output, including unavailable or negative-component points. No alternative point estimate or sampling SE is supplied. LogVarianceRatioLower/Upper concern log(Person variance / residual variance), with -Inf meaning zero; they are not mean-square ratios or G endpoints.",
    "Every plan for this prespecified score maps one nonnegative variance-ratio confidence set. The same event protects all returned G plans or differences jointly; no guarantee extends to other scores, metrics or adaptively chosen weights. An empty set, absent residual contrasts or unresolved numerical information remains unavailable, not a zero-width interval.",
    "Unresolved F quantiles, ratio bounds or mapped interval widths remain unavailable with an explicit Status. The pivot retains endpoint_tolerance and quantile_resolution; these numerical checks do not test normality or missingness assumptions. Identical-count differences remain exactly zero only when the base set is available.",
    "Residual combines the Person-by-facet interaction and observation error. This method does not provide Phi, SEM, nested or multifacet intervals. Save the full object, including pivot diagnostics, for source-matched summaries and plots without recalculation.",
    if (!is.null(reference)) "Differences are comparison minus reference and map the entire ratio set, including its interior extremum. They do not subtract marginal endpoints. Positive differences favor relative reliability, not necessarily cost-effectiveness.")
  p <- vapply(g$components, function(m) drop(crossprod(weights, m %*% weights)), numeric(1))
  .mfrm_mvgt_interval_result(out, d, list(estimate = p), "spectral_f")
}

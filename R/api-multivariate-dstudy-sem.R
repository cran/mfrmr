# One common residual SD determines every one-facet relative SEM and difference.
.mfrm_mvgt_residual_sem <- function(d, tab, weights, level, reference) {
  g <- d$gstudy
  components <- .mfrm_mvgt_residual_intervals(g, NULL, NULL, weights, level)
  base <- components[components$Source == "Residual", , drop = FALSE]
  counts <- d$design_grid[[1L]]
  multiplier <- 1/sqrt(counts)
  if (!is.null(reference)) {
    a <- counts[reference]
    # Rationalize the subtraction to retain differences between large counts.
    multiplier <- ((a-counts)/(sqrt(a)+sqrt(counts)))/sqrt(a)/sqrt(counts)
  }
  center <- if (is.finite(base$IntervalCenter) && base$IntervalCenter >= 0)
    sqrt(base$IntervalCenter) else NA_real_
  lower <- sqrt(base$Lower); upper <- sqrt(base$Upper)
  out <- tab[c("Scenario", names(d$design_grid), "Kind", "Score")]
  out$Metric <- "RelativeSEM"; out$Estimate <- tab$RelativeSEM
  out$PointStatus <- tab$RelativeSEMStatus
  out$IntervalCenter <- multiplier*center
  out$Lower <- pmin(multiplier*lower, multiplier*upper)
  out$Upper <- pmax(multiplier*lower, multiplier*upper)
  out$Level <- level; out$Method <- "Residual-contrast chi-square SEM interval"
  out$Status <- base$Status
  out$Validation <- base$Validation
  if (!is.null(reference)) {
    out$Reference <- as.integer(reference); out$ReferenceCount <- a
    out$ReferenceValue <- tab$RelativeSEM[reference]; out$Value <- tab$RelativeSEM
    out$Metric <- "RelativeSEM difference"
    out$Estimate <- out$Value-out$ReferenceValue
    out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
    if (base$Available) out$Status <- ifelse(counts == a, "Identical counts", "Exact normal-model interval")
    out$Method <- "Residual-contrast SEM difference set image"
    attr(out, "reference") <- list(scenario = as.integer(reference), count = a,
      value = tab$RelativeSEM[reference], point_status = tab$RelativeSEMStatus[reference])
    out <- out[-reference, , drop = FALSE]; rownames(out) <- NULL
  }
  if (base$Available) {
    collapsed <- base$Lower < base$Upper & (lower == upper |
      (out$Lower == out$Upper & multiplier[out$Scenario] != 0))
    out$Lower[collapsed] <- out$Upper[collapsed] <- NA_real_
    out$Status[collapsed] <- "SEM interval below numeric resolution"
  }
  attr(out, "weights") <- weights; attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- attr(base, "source_design")
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "source_data") <- g$data
  attr(out, "base_interval") <- base
  attr(out, "residual_projection") <- attr(base, "residual_projection")
  keep <- if (is.null(reference)) seq_along(counts) else seq_along(counts)[-reference]
  attr(out, "diagnostics") <- setNames(lapply(keep, function(i)
    list(count = counts[i], multiplier = multiplier[i], residual_df = attr(base, "residual_projection")$df)),
    as.character(out$Scenario))
  attr(out, "notes") <- c(
    "One-facet relative SEM equals the residual SD divided by sqrt(future count). The source may be incomplete; future plans remain complete balanced assignments. Independent homoskedastic normal residuals and omissions preserving this model are required.",
    "All plan intervals transform one residual chi-square event for this prespecified score/composite. The same event protects all returned plans or differences jointly under the model, but gives no joint guarantee across scores, metrics or adaptively selected weights. More future ratings add no source information.",
    "Estimate and PointStatus retain the original ANOVA/MINQUE SEM or difference, including unavailable points. IntervalCenter is the separate projected estimate. Both and their limits are in score units; an original point can lie outside the interval. No sampling SE of the original SEM estimator is supplied.",
    "base_interval retains the residual variance interval and original/projected variance estimates in squared score units. It is not a SEM interval. residual_projection retains degrees of freedom and numerical checks. No likelihood fitting or pivotal simulation is used.",
    "No residual contrasts, numerically unresolved residual variation or a failed base interval leaves bounds unavailable. Identical-count differences are [0,0] only when the base interval is available; failed base calculations are not converted into successful zero-width intervals.",
    "This method supplies one-facet RelativeSEM only, not AbsoluteSEM, G/Phi or multifacet SEM. Save the full result with saveRDS(); summaries and curves can reuse its source-matched intervals without recalculation.",
    if (!is.null(reference)) "Differences are comparison minus reference. Multiply the common residual-SD interval by 1/sqrt(comparison count)-1/sqrt(reference count), reversing endpoints for a negative multiplier. Negative SEM differences favor the comparison plan. Marginal interval endpoints are not subtracted.")
  .mfrm_mvgt_interval_result(out, d, list(estimate = components$Estimate), "residual_contrast")
}

# Error-variance intervals use the independent ANOVA mean squares directly.
# Taking their square roots preserves the variance interval's coverage event.
.mfrm_mvgt_sem_table <- function(d, tab, p, metric, level, reference) {
  g <- d$gstudy
  absolute <- identical(metric, "AbsoluteSEM")
  variance <- if (absolute) "AbsoluteErrorVariance" else "RelativeErrorVariance"
  rows <- diagnostics <- vector("list", nrow(tab))
  for (i in seq_len(nrow(tab))) {
    counts <- as.numeric(d$design_grid[i, ])
    persons <- g$design$counts[1L]
    if (!is.null(g$design$nesting)) {
      # Component order follows parent then child, independently of the
      # user's facet-column order; child counts mean children per parent.
      child <- counts[match(g$design$count_columns[[names(g$design$nesting)]], names(d$design_grid))]
      parent <- counts[match(g$design$count_columns[[unname(g$design$nesting)]], names(d$design_grid))]
      w <- c(0, 0, 0, 1/parent, 1/(child*parent))
      if (absolute) w <- w + c(0, 1/parent, 1/(child*parent), 0, 0)
      observed_child <- g$design$child_counts[1L]
      main <- 1/(observed_child*parent)
      remainder <- (observed_child-child)/(observed_child*child*parent)
      a <- if (absolute) c(0, main/persons, remainder/persons,
        (persons-1)*main/persons, (persons-1)*remainder/persons) else
          c(0, 0, 0, main, remainder)
    } else if (length(counts) == 1L) {
      w <- c(0, as.integer(absolute), 1)/counts[1L]
      a <- if (absolute) c(0, 1/persons, (persons-1)/persons)/counts[1L] else w
    } else {
      w <- c(0, 0, 0, 1/counts[1L], 1/counts[2L], 0, 1/prod(counts))
      if (absolute) w <- w + c(0, 1/counts[1L], 1/counts[2L], 0, 0, 1/prod(counts), 0)
      observed <- as.numeric(g$design$counts[-1L])
      first <- 1/(observed[2L]*counts[1L]); second <- 1/(observed[1L]*counts[2L])
      remainder <- (prod(observed)-observed[1L]*counts[2L]-observed[2L]*counts[1L]) /
        (prod(observed)*prod(counts))
      a <- if (absolute) c(0, first/persons, second/persons,
        (persons-1)*first/persons, (persons-1)*second/persons,
        remainder/persons, (persons-1)*remainder/persons) else
          c(0, 0, 0, first, second, 0, remainder)
    }
    # These algebraic forms of w' A use integer count differences before
    # division: a rounding residual must not select a signed MLS formula.
    ans <- .mfrm_mvgt_linear_variance_interval(p$ms, p$design$df, a, level)
    row <- tab[i, c("Scenario", names(d$design_grid), "Kind", "Score"), drop = FALSE]
    row$Metric <- metric; row$Estimate <- tab[[metric]][i]
    row$PointStatus <- tab[[paste0(metric, "Status")]][i]
    row$ErrorVariance <- tab[[variance]][i]
    row$ErrorVarianceSE <- ans$se * p$restore
    row$RawVarianceLower <- ans$raw[1L] * p$restore
    row$RawVarianceUpper <- ans$raw[2L] * p$restore
    # Restore the SD scale after the square root to avoid squaring large
    # scores when the SEM endpoints themselves remain representable.
    row$Lower <- sqrt(ans$bounds[1L]) * sqrt(p$restore)
    row$Upper <- sqrt(ans$bounds[2L]) * sqrt(p$restore)
    row$Status <- ans$status
    if (ans$status %in% c("Available", "Includes zero")) {
      limits <- c(row$Lower, row$Upper)
      if (any(!is.finite(limits)) || any(limits == 0 & ans$bounds > 0)) {
        row$Status <- "Unrepresentable SEM endpoints"
      } else if (limits[1L] == limits[2L] && ans$bounds[1L] < ans$bounds[2L]) {
        row$Status <- "SEM interval below numeric resolution"
      }
      if (!row$Status %in% c("Available", "Includes zero")) row$Lower <- row$Upper <- NA_real_
    }
    row$Level <- level
    row$Method <- paste(ans$method, "SEM interval")
    row$Validation <- if (ans$method == "Chi-square")
      "Exact normal-model interval for positive mean-square expectation" else
      "MLS approximation; pointwise coverage is not guaranteed"
    rows[[i]] <- row
    diagnostics[[i]] <- list(component_weights = stats::setNames(w, names(g$components)),
      mean_square_coefficients = stats::setNames(a, names(g$components)),
      normalized_error_variance = ans$estimate, normalized_raw_bounds = ans$raw,
      method = ans$method, status = ans$status)
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "weights") <- p$weights; attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- g$design[c("facets", "counts", "nesting", "child_counts",
    "method", "balanced", "complete", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  names(diagnostics) <- as.character(tab$Scenario)
  attr(out, "diagnostics") <- diagnostics
  attr(out, "mean_squares") <- list(normalized = p$ms, scale = p$restore,
    df = p$design$df, component_map = p$design$A)
  attr(out, "notes") <- c(
    "Intervals quantify sampling uncertainty in population measurement SEM, not intervals for Person scores. Source variance components are combined before forming an interval; component endpoints are not added.",
    "One contributing mean square uses a chi-square interval; positive sums use Graybill-Wang MLS and signed sums use Ting MLS. Square-root transformation preserves the underlying variance-interval coverage; it does not make MLS exact or distribution robust.",
    "Intervals are pointwise under independent normal random effects and fixed score weights. They are not generally simultaneous across plans or scores. Missing/empty intervals and negative component estimates remain explicit.",
    "ErrorVariance, ErrorVarianceSE and RawVarianceLower/Upper have squared score units; Lower/Upper and Estimate have score units. A negative error-variance estimate can leave a point unavailable while an interval exists.",
    "The nonnegative variance-space intersection is transformed before drawing the SEM interval. If all contributing mean squares are zero, the interval is unavailable; a zero SE is not proof of exact population knowledge. Positive intervals whose width or endpoints cannot be represented remain unavailable with an explicit Status.",
    "In a limited normal-model screen, nominal 95% coverage for positive error variances ranged from 93.45% to 99.00% (2,000 replications per source condition); approximate coverage is not assured. More future ratings do not increase source-study information.")
  if (!is.null(reference)) {
    # All one-facet plans share one single-rating error variance. Map that
    # one interval; marginal endpoint subtraction would discard dependence.
    counts <- if (is.null(g$design$nesting)) d$design_grid[[1L]] else
      d$design_grid[[g$design$count_columns[[unname(g$design$nesting)]]]]
    factor <- (counts[reference]-counts) /
      (sqrt(counts) * (sqrt(counts[reference])+sqrt(counts)))
    base <- c(out$Lower[reference], out$Upper[reference])
    bounds <- t(vapply(factor, function(k) if (all(is.finite(base)))
      sort(k * base) else c(NA_real_, NA_real_), numeric(2)))
    out$Reference <- as.integer(reference); out$ReferenceCount <- counts[reference]
    out$ReferenceValue <- tab[[metric]][reference]; out$Value <- tab[[metric]]
    out$Metric <- paste(metric, "difference"); out$Estimate <- out$Value - out$ReferenceValue
    out$PointStatus <- ifelse(is.finite(out$Estimate), "Available", "Point difference unavailable")
    out$Lower <- bounds[,1L]; out$Upper <- bounds[,2L]
    out$Status <- if (all(is.finite(base))) ifelse(counts == counts[reference],
      "Identical counts", ifelse(out$Lower <= 0 & out$Upper >= 0, "Includes zero", "Available")) else
        rep(out$Status[reference], nrow(out))
    # Retain the common source interval instead of presenting individual-plan
    # variance diagnostics as though they were variances of differences.
    attr(out, "reference") <- list(scenario = as.integer(reference), count = counts[reference],
      value = tab[[metric]][reference], point_status = tab[[paste0(metric,"Status")]][reference],
      sem_interval = base, variance_estimate = out$ErrorVariance[reference],
      variance_se = out$ErrorVarianceSE[reference],
      variance_raw_bounds = c(out$RawVarianceLower[reference], out$RawVarianceUpper[reference]))
    attr(out, "notes")[4L] <- "Reference SEM and difference intervals have score units. Reference error-variance estimates, SEs and raw bounds have squared score units and are retained in the reference attribute."
    out[c("ErrorVariance", "ErrorVarianceSE", "RawVarianceLower", "RawVarianceUpper")] <- NULL
    out$Method <- sub("SEM interval$", "SEM difference set image", out$Method)
    out <- out[-reference, , drop = FALSE]; rownames(out) <- NULL
    attr(out, "notes") <- c(attr(out, "notes"),
      "Differences are comparison minus reference. Each compared SEM equals a common error SD divided by sqrt(future count), using the parent count for nested plans with a common child count. The whole common interval is multiplied by the count difference factor. Negative SEM differences favor the comparison plan.",
      "The same base coverage event protects all returned differences for this score jointly. Exact chi-square or approximate MLS coverage is retained according to Method and Validation; no extra independent-plan assumption is used.",
      "For differences, reference retains the base SEM interval and error-variance diagnostics; no variance of the SEM difference or delta-method SE is implied.")
  }
  .mfrm_mvgt_interval_result(out, d, p, "mls")
}

.mfrm_mvgt_linear_variance_interval <- function(ms, df, a, level) {
  z <- abs(a) * ms; scale <- max(z)
  estimate <- sum(a * ms)
  se <- if (scale == 0) 0 else sqrt(sum((z/scale)^2 * 2/(df+2))) * scale
  method <- if (sum(a != 0) == 1L && all(a >= 0)) "Chi-square" else
    if (all(a >= 0)) "Graybill-Wang MLS" else "Ting MLS"
  raw <- bounds <- c(NA_real_, NA_real_)
  status <- "Degenerate mean squares"
  if (scale > 0) {
    tail <- (1-level)/2
    raw <- if (method == "Chi-square") {
      j <- which(a != 0)
      estimate * (df[j]/c(stats::qchisq(tail, df[j], lower.tail = FALSE), stats::qchisq(tail, df[j])))
    } else estimate + c(-.mfrm_mvgt_mls_width(ms, df, a, tail),
      .mfrm_mvgt_mls_width(ms, df, -a, tail))
    status <- if (any(!is.finite(raw))) "Numerical interval failure" else
      if (method == "Chi-square" && any(raw <= 0)) "Unrepresentable chi-square endpoints" else
      if (raw[1L] == raw[2L]) "Variance interval below numeric resolution" else
      if (raw[2L] < 0) "Empty nonnegative interval" else
        if (raw[1L] <= 0) "Includes zero" else "Available"
    if (status %in% c("Available", "Includes zero")) bounds <- c(max(0, raw[1L]), raw[2L])
  }
  list(estimate = estimate, se = se, raw = raw, bounds = bounds, status = status, method = method)
}

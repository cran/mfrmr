# First-order coefficient uncertainty from the original component estimator.
.mfrm_mvgt_delta_table <- function(d, tab, weights, metric, level) {
  g <- d$gstudy
  sampling <- vcov(g, weights = weights/max(abs(weights)), assumption = "normal")
  info <- attr(sampling, "diagnostics")
  theta <- attr(sampling, "summary")$Estimate/info$component_scale
  C <- info$normalized_covariance
  critical <- stats::qnorm((1-level)/2, lower.tail = FALSE)
  out <- tab[c("Scenario", names(d$design_grid), "Kind", "Score")]
  out$Metric <- metric; out$Estimate <- tab[[metric]]
  out$PointStatus <- tab[[paste0(metric, "Status")]]
  out[c("SE", "SamplingVariance", "Lower", "Upper")] <- NA_real_
  rownames(out) <- NULL
  out$SEStatus <- out$Status <- "Point projection unavailable"
  out$AtBoundary <- FALSE
  diagnostics <- vector("list", nrow(tab))
  for (i in seq_len(nrow(tab))) {
    error <- .mfrm_mvgt_profile_error(d, i, metric)
    u <- theta[1L]; e <- sum(error*theta)
    gradient <- rep(NA_real_, length(theta))
    if (is.finite(out$Estimate[i]) && is.finite(u+e) && u+e > 0) {
      gradient <- -out$Estimate[i]*error/(u+e)
      gradient[1L] <- (1-out$Estimate[i])/(u+e)
      out$AtBoundary[i] <- u == 0 || e == 0
      # Unidentified covariance entries with zero weights do not contribute.
      active <- which(gradient != 0)
      variance <- if (all(is.finite(gradient)))
        drop(crossprod(gradient[active], C[active, active, drop=FALSE] %*% gradient[active])) else NA_real_
      out$SamplingVariance[i] <- variance
      out$SEStatus[i] <- if (!is.finite(variance)) "Sampling variance unavailable" else
        if (variance < 0) "Negative sampling-variance estimate" else
          if (variance == 0) "Zero first-order variance" else "Available"
      out$Status[i] <- out$SEStatus[i]
      if (is.finite(variance) && variance >= 0) out$SE[i] <- sqrt(variance)
      if (is.finite(variance) && variance > 0) {
        bounds <- out$Estimate[i] + c(-1, 1)*critical*out$SE[i]
        if (all(is.finite(bounds)) && bounds[1L] == bounds[2L]) {
          out$Status[i] <- "Delta interval below numeric resolution"
        } else if (all(is.finite(bounds))) {
          out$Lower[i] <- bounds[1L]; out$Upper[i] <- bounds[2L]
          out$Status[i] <- "Computed delta approximation"
        } else out$Status[i] <- "Interval exceeds numeric range"
      }
    }
    diagnostics[[i]] <- list(error_coefficients=error, normalized_gradient=gradient,
      component_scale=info$component_scale)
  }
  out$Level <- level; out$Method <- "Delta normal approximation"
  out$SamplingCovariancePSD <- info$positive_semidefinite
  out$OutsideParameterSpace <- out$Lower < 0 | out$Upper > 1
  out$Validation <- "First-order normal approximation; small-sample and boundary coverage not guaranteed"
  attr(out, "sampling_covariance") <- sampling
  attr(out, "diagnostics") <- diagnostics
  attr(out, "weights") <- weights; attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- g$design[c("facets", "counts", "method", "balanced", "complete", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  if (identical(g$design$method, "Multivariate MINQUE(0)")) attr(out, "source_data") <- g$data
  attr(out, "notes") <- c(
    "Explicit delta-method sampling SEs and normal-reference intervals for the original ANOVA/MINQUE G or Phi estimator, under independent normal random effects. SE is uncertainty in the coefficient, not measurement SEM or the SE of a REML estimate.",
    "The full component sampling covariance comes from vcov(g), with its Gaussian moment correction. Balanced ANOVA uses 2*MS^2/(df+2), as in the random crossed derivation of Almehrizi (2021, equations 21-25). Nested and incomplete sources use the same gradient with their own component covariance; the paper does not validate those extensions.",
    "The covariance correction does not make a nonlinear coefficient SE unbiased or the normal-reference interval exact. Few facet levels, weak source information, component boundaries or nonnormal effects can impair coverage. Assignment and omissions must preserve the stated model. Future plans remain complete balanced ratings.",
    "More source Persons alone need not justify a normal reference: with a fixed random-rater count and positive rater variance, the rater-population variance entering Phi is not consistently estimated as Persons increase. Future plan counts also add no source information.",
    "Lower and Upper are Estimate +/- the normal critical value times SE, without clipping to [0,1]. OutsideParameterSpace and AtBoundary flag limits of the approximation. A wide interval remains visible; no method is chosen because its interval is narrower.",
    "Negative component and sampling-variance estimates are retained. An indefinite component sampling covariance is flagged, but a finite positive projected variance still supplies an explicit approximation. Undefined points or nonfinite/negative projected variances cannot supply a Wald interval. A zero first-order variance retains SE=0 but withholds a degenerate normal-reference interval. A positive SE is retained when the requested width rounds away, with unavailable limits and a numeric-resolution status.",
    "Intervals are pointwise for a prespecified score/composite and plan. Comparing methods on one dataset is sensitivity analysis, not a coverage test. This method supplies neither SEM nor reference-plan differences; mfrm_multivariate_d_compare() uses its separately documented plug-in covariance for paired differences.",
    "No refitting or sampling is performed. The retained sampling_covariance uses weights divided by their largest absolute value; the gradient uses the matching normalized component scale, while weights retains the original score definition. Save the complete object for summaries and plots without recalculation.")
  .mfrm_mvgt_interval_result(out, d, list(estimate=theta), "delta")
}

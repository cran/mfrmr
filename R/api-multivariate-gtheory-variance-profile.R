# Absolute targets require the residual scale that coefficient profiles remove.
.mfrm_mvgt_variance_at <- function(source, target, weights) {
  k <- length(weights); df <- length(source$y)-1L
  active <- which(weights > 0)
  scalar <- length(active) == 1L && active < k
  if (target == 0 && weights[k] > 0)
    return(list(target = target, weights = weights, valid = FALSE, lr = NA_real_,
      check = list(kkt = NA_real_, error = NA_real_)))
  if (target == 0) {
    free <- setdiff(seq_len(k-1L), active)
    decode <- function(par) { v <- numeric(k-1L); v[free] <- expm1(par); v }
    objective <- function(par) source$objective(decode(par))
    starts <- lapply(list(source$fit$par[free], rep(1, length(free)),
      rep(.1, length(free))), log1p)
    lower <- rep(0, length(free))
  } else {
    index <- if (scalar) active else integer()
    decode <- function(par) {
      v <- expm1(par)
      if (scalar) v[index] <- exp(par[index])
      v
    }
    objective <- function(par) {
      v <- decode(par)
      if (any(!is.finite(v)) || (scalar && v[index] == 0)) return(Inf)
      at <- source$objective(v, details = TRUE)
      logscale <- log(target)-(if (scalar) log(weights[index])+par[index] else
        log(sum(weights*c(v, 1))))
      logratio <- logscale-log(at$sigma2)
      at$objective+df*(logratio+expm1(-logratio))
    }
    starts <- lapply(list(source$fit$par, rep(1, k-1L), rep(.1, k-1L)), log1p)
    lower <- rep(0, k-1L)
    if (scalar) {
      lower[index] <- -Inf
      for (j in seq_along(starts))
        starts[[j]][index] <- log(target)-log(weights[index])-log(source$sigma2)
    }
  }
  fitted <- .mfrm_mvgt_profile_minimize(starts, objective, lower)
  fit <- fitted$best; v <- decode(fit$par)
  at <- source$objective(v, details = TRUE)
  sigma2 <- if (target == 0) at$sigma2 else target/sum(weights*c(v, 1))
  check <- if (is.finite(fit$objective))
    .mfrm_mvgt_profile_gradient(fit$par, objective, lower, refinements = 4L) else
      list(kkt = NA_real_, error = NA_real_, valid = FALSE)
  lr <- fit$objective-source$fit$objective
  list(target = target, weights = weights, relative_variances = c(v, 1),
    sigma2 = sigma2, components = sigma2*c(v, 1), lr = lr, fit = fit,
    attempts = fitted$attempts, check = check,
    coordinates = "log1p variance ratios; log ratio for a fixed positive nonresidual component",
    valid = source$valid && fit$convergence == 0 && check$valid &&
      is.finite(lr) && lr >= -1e-6)
}

.mfrm_mvgt_variance_set <- function(source, weights, level) {
  ans <- list(bounds = c(NA_real_, NA_real_), candidate_bounds = c(NA_real_, NA_real_),
    center = NA_real_, status = source$status, error_coefficients = weights)
  if (is.null(source$components)) return(ans)
  center <- sum(weights*source$components)
  ans$center <- center
  anchor <- if (center > 0) center else sum(weights)*source$sigma2
  if (!source$valid || !is.finite(anchor) || anchor <= 0) return(ans)
  cutoff <- .mfrm_mvgt_profile_cutoff(level, source$fit$objective)
  ans$cutoff <- cutoff; ans$critical <- critical <- cutoff$critical
  if (!cutoff$valid) {
    ans$status <- "Likelihood-ratio cutoff below numeric resolution"
    return(ans)
  }
  zero_allowed <- tail(weights, 1) == 0
  cache <- new.env(parent = emptyenv()); failures <- character()
  transform <- function(x) anchor*expm1(x)
  xcenter <- log1p(center/anchor)
  evaluate <- function(x) {
    key <- sprintf("%.17g", x)
    if (!exists(key, cache, inherits = FALSE)) {
      p <- .mfrm_mvgt_variance_at(source, transform(x), weights)
      p$parameter <- x; assign(key, p, cache)
    }
    get(key, cache, inherits = FALSE)
  }
  statistic <- function(x) { p <- evaluate(x); if (p$valid) p$lr-critical else NA_real_ }
  atcenter <- statistic(xcenter)
  endpoint <- function(path, side) {
    previous <- xcenter
    if (!is.finite(atcenter) || atcenter > 1e-6) return(NA_real_)
    if (!length(path)) return(0)
    for (x in path) {
      value <- statistic(x)
      if (!is.finite(value)) {
        failures <<- c(failures, paste(side, "profile evaluation unresolved"))
        return(NA_real_)
      }
      if (value >= 0) {
        root <- tryCatch(stats::uniroot(statistic, sort(c(previous, x)),
          tol = min(1e-10, abs(x-previous)*1e-7), check.conv = TRUE)$root,
          error = function(e) { failures <<- c(failures, conditionMessage(e)); NA_real_ })
        if (!is.finite(root)) return(NA_real_)
        residual <- statistic(root)
        if (!is.finite(residual) || abs(residual) > cutoff$tolerance) {
          failures <<- c(failures, paste(side, "endpoint does not recover cutoff"))
          return(NA_real_)
        }
        return(transform(root))
      }
      previous <- x
    }
    if (side == "Lower" && zero_allowed) return(0)
    failures <<- c(failures, paste(side, "tail remains unresolved"))
    NA_real_
  }
  lower_path <- if (zero_allowed) {
    if (xcenter == 0) numeric() else
      unique(c(pmax(0, xcenter-.25*seq_len(ceiling(xcenter/.25))), 0))
  } else log1p(center/anchor*exp(-.5*seq_len(72L)))
  lower <- endpoint(lower_path, "Lower")
  upper <- endpoint(xcenter+.25*seq_len(144L), "Upper")
  profiles <- as.list(cache)
  scan <- do.call(rbind, lapply(profiles, function(p) data.frame(Target = p$target,
    Parameter = p$parameter, LR = p$lr, KKT = p$check$kkt,
    GradientError = p$check$error, Valid = p$valid)))
  scan <- scan[order(scan$Target), ]; rownames(scan) <- NULL
  ans$anchor <- anchor; ans$critical <- critical
  ans$profiles <- profiles; ans$scan <- scan; ans$failures <- failures
  ans$candidate_bounds <- c(lower, upper)
  valid <- all(scan$Valid) && !length(failures) &&
    all(is.finite(ans$candidate_bounds)) && lower <= center && upper >= center && lower < upper
  ans$status <- "Unresolved numerical profile; inspect candidate bounds and diagnostics"
  if (valid) {
    ans$bounds <- ans$candidate_bounds
    ans$status <- "Computed profile approximation around REML center"
  }
  ans
}

.mfrm_mvgt_component_profile <- function(g, parm, score, weights, level) {
  scores <- g$design$scores
  facets <- names(g$design$facets)
  if (!inherits(g, "mfrm_multivariate_gstudy") ||
      !.mfrm_mvgt_valid_facets(facets) || !is.character(scores) || !length(scores) ||
      anyNA(scores) || anyDuplicated(scores) ||
      !identical(names(g$components), .mfrm_mvgt_sources(c("Person", facets), g$design$nesting)) ||
      !all(vapply(g$components, function(m) is.matrix(m) && is.numeric(m) &&
        !is.complex(m) && identical(dimnames(m), list(scores, scores)) &&
        all(is.finite(m)) && identical(m, t(m)), logical(1))))
    stop("The G-study components or design identities are incomplete or altered; recreate the result.", call. = FALSE)
  selection <- .mfrm_mvgt_interval_score(scores, score, weights)
  weights <- selection$weights
  sources <- names(g$components)
  selected <- .mfrm_mvgt_component_selection(parm, sources)
  source <- .mfrm_mvgt_profile_source(g, weights)
  estimate <- source$raw_estimate
  if (any(!is.finite(estimate)))
    stop("Component units exceed numeric range; rescale the scores or weights.", call. = FALSE)
  if (source$valid && (!is.finite(source$score_scale) || source$score_scale <= 0))
    stop("Component units exceed numeric range; rescale the scores or weights.", call. = FALSE)
  restore <- function(x) (x*source$score_scale)*source$score_scale
  rows <- diagnostics <- vector("list", length(selected))
  for (j in seq_along(selected)) {
    i <- selected[j]
    a <- setNames(as.numeric(seq_along(sources) == i), sources)
    ans <- .mfrm_mvgt_variance_set(source, a, level)
    reported <- .mfrm_mvgt_profile_restore(ans, restore)
    values <- reported$values
    rows[[j]] <- data.frame(Source = sources[i], Score = selection$score,
      Estimate = unname(estimate[i]), PointStatus = if (estimate[i] < 0) "Negative estimate" else "Nonnegative estimate",
      IntervalCenter = values[1L], SE = NA_real_, SEStatus = "Not supplied by profile inversion",
      Lower = values[2L], Upper = values[3L], RawLower = NA_real_, RawUpper = NA_real_,
      CandidateLower = values[4L], CandidateUpper = values[5L], Level = level,
      Method = "REML component profile (approximate)", MethodCode = "reml_profile",
      Status = reported$status, Available = reported$available,
      Validation = "Experimental chi-square calibration; boundary and finite-sample coverage not guaranteed")
    diagnostics[[j]] <- ans
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  names(diagnostics) <- sources[selected]
  source$objective <- NULL
  attr(out, "reml_source") <- source
  attr(out, "source_data") <- g$data
  attr(out, "source_summary") <- g[c("design", "components", "mean_products", "degrees_of_freedom")]
  attr(out, "diagnostics") <- diagnostics
  attr(out, "weights") <- weights; attr(out, "kind") <- selection$kind
  attr(out, "assumption") <- "normal"
  attr(out, "source_design") <- g$design[c("facets", "counts", "nesting", "child_counts",
    "balanced", "complete", "method", "score_convention")]
  attr(out, "source_rows") <- g$data_usage$counts
  attr(out, "notes") <- c(
    "Explicit experimental component intervals under independent normal random effects, a common mean and positive residual variance. Source omissions must be noninformative conditional on the declared fixed assignment. Raw used ratings are needed for fitting.",
    "Estimate retains the original ANOVA/MINQUE component, including negative estimates. IntervalCenter is the separate REML variance estimate; both target the same population component. The original point may lie outside the interval. Component matrices are not repaired.",
    "Estimate, IntervalCenter and interval endpoints have squared score units. SE is unavailable because this profile inversion does not calculate a sampling SE; it is not zero. RawLower/RawUpper are inapplicable MLS columns and remain NA, not unconstrained REML limits.",
    "Fix each component in the Gaussian restricted likelihood, optimize the remaining variances and invert a chi-square(1) cutoff. Nonresidual components may be exactly zero; residual variance is strictly positive. Component intervals cannot be added to obtain error-variance or coefficient intervals.",
    "Search outward from the REML center to the first checked cutoff crossing in each direction. The returned numerical interval contains that center but does not enumerate remote disconnected accepted regions. Multiple starts, a second optimizer and boundary-gradient checks do not certify global optima or rule out undetected crossings.",
    "Available denotes numerical completion only. An unresolved source, required evaluation, root or tail leaves accepted bounds unavailable; resolved CandidateLower/Upper remain exploratory. A search cap is not an infinite endpoint. The chi-square calibration is approximate, not an exact boundary mixture or finite-sample guarantee.",
    "Limited normal-model component checks with 40 samples per condition gave true-value inclusion of 92.5-100% in positive-component populations and 90-100% for checked zero components, counting numerical failures as misses. This coarse feasibility screen is not full interval-endpoint coverage or a generally recommended method.",
    "Intervals are pointwise for prespecified scores/weights and components. They are not simultaneous or distribution robust; scalar projection does not provide off-diagonal covariance intervals. The method is not selected automatically.",
    "Save the complete result with saveRDS() to retain used ratings, weights, diagnostics and assumptions. Diagnostics use standardized-score variance units; multiply by reml_source$score_scale twice to restore squared score units. Plot original estimates and preserve the interval notes when sharing figures.",
    "A likelihood-ratio cutoff below numerical resolution, or a nonzero interval lost when restoring score units, leaves bounds unavailable with a specific Status. Estimates and representable REML centers are retained; diagnostics$cutoff records the criterion tolerance and resolution.")
  out
}

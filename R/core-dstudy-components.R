# D-study conventions never overwrite the original G-study component estimates.
.mfrm_mvds_components <- function(g, policy) {
  if (identical(policy, "raw")) return(g$components)
  if (length(g$design$scores) != 1L ||
      !identical(g$design$method, "Balanced multivariate ANOVA") ||
      !isTRUE(g$design$balanced) || !isTRUE(g$design$complete))
    stop("GENOVA component treatments require a complete balanced univariate ANOVA G-study.", call. = FALSE)
  p <- .mfrm_mvgt_anova_projection(g)
  if (identical(policy, "genova_algorithm")) {
    effective <- pmax(p$estimate, 0)
  } else {
    # EMS = L * components, using SOURCE counts, before future-plan scaling.
    # Back substitution uses already zero-replaced higher-order components.
    L <- solve(p$design$A)
    effective <- numeric(length(p$ms))
    for (i in rev(seq_along(effective))) {
      higher <- which(seq_along(effective) > i)
      effective[i] <- max(0, (p$ms[i] - sum(L[i, higher] * effective[higher])) / L[i, i])
    }
  }
  effective <- effective * p$restore
  if (any(!is.finite(effective)))
    stop("Adjusted components exceed numerical range; rescale the scores.", call. = FALSE)
  out <- g$components
  for (i in seq_along(out)) out[[i]][] <- effective[i]
  out
}

.mfrm_mvds_policy <- function(x) {
  policy <- x$negative_components
  if (is.null(policy)) {
    recorded <- x$coefficients$ComponentTreatment
    if (!is.null(x$effective_components) ||
        (!is.null(recorded) && (anyNA(recorded) || any(recorded != "raw"))))
      stop("The saved component treatment is missing; recreate the D-study.", call. = FALSE)
    return("raw") # Before this option, all saved D-studies retained raw components.
  }
  if (!is.character(policy) || length(policy) != 1L || is.na(policy) ||
      !policy %in% c("raw", "genova_algorithm", "genova_ems"))
    stop("The saved component treatment is invalid; recreate the D-study.", call. = FALSE)
  recorded <- x$coefficients$ComponentTreatment
  if (is.null(recorded) || anyNA(recorded) || any(recorded != policy))
    stop("D-study results do not match the stored component treatment; recreate the result.", call. = FALSE)
  policy
}

.mfrm_mvds_component_note <- function(x) {
  switch(.mfrm_mvds_policy(x),
    raw = "Component treatment: raw estimates retained.",
    genova_algorithm = "Component treatment: GENOVA ALGORITHM; negative variance components replaced by zero for D-study projections. Original G-study estimates are retained.",
    genova_ems = "Component treatment: GENOVA EMS; components recalculated with zero replacements using source-design mean-square equations. Original G-study estimates are retained.")
}

.mfrm_mvds_require_raw <- function(x) {
  if (!identical(.mfrm_mvds_policy(x), "raw"))
    stop("Adjusted-estimator delta comparisons are not available for GENOVA component treatments. Use confint(d, reference = ...) for supported population plan-difference intervals, or recreate the raw D-study for original-estimator delta comparisons.", call. = FALSE)
  invisible(NULL)
}

# Point reporting does not change the population target or its source study.
# Keep the original interval procedure intact, including failures and any
# separately fitted center. Never attach its SE to the adjusted estimator.
.mfrm_mvds_adjusted_intervals <- function(ci, d) {
  metric <- sub(" difference$", "", unique(ci$Metric))
  tab <- d$coefficients
  tab <- tab[tab$Kind == ci$Kind[1L] & tab$Score == ci$Score[1L], , drop = FALSE]
  idx <- match(ci$Scenario, tab$Scenario)
  ci$IntervalEstimate <- ci$Estimate
  ci$IntervalPointStatus <- ci$PointStatus
  ci$ComponentTreatment <- .mfrm_mvds_policy(d)
  ci$IntervalBasis <- "Original source study; population parameter"
  ref <- attr(ci, "reference")
  if (is.null(ref)) {
    ci$Estimate <- tab[[metric]][idx]
    ci$PointStatus <- tab[[paste0(metric, "Status")]][idx]
  } else {
    ci$IntervalValue <- ci$Value
    ci$IntervalReferenceValue <- ci$ReferenceValue
    ci$Value <- tab[[metric]][idx]
    ci$ReferenceValue <- tab[[metric]][match(ref$scenario, tab$Scenario)]
    ci$Estimate <- ci$Value - ci$ReferenceValue
    ci$PointStatus <- ifelse(is.finite(ci$Estimate), "Available", "Point difference unavailable")
    attr(ci, "interval_reference") <- ref
    ref$value <- ci$ReferenceValue[1L]
    ref$point_status <- tab[[paste0(metric, "Status")]][match(ref$scenario, tab$Scenario)]
    attr(ci, "reference") <- ref
  }
  # Rename SE and SEStatus together: retaining only SEStatus would make R's
  # partial `$SE` matching return a character status in place of a numeric SE.
  se_columns <- grep("^SE($|[A-Z])", names(ci), value = TRUE)
  names(ci)[match(se_columns, names(ci))] <- paste0("Interval", se_columns)
  attr(ci, "negative_components") <- .mfrm_mvds_policy(d)
  attr(ci, "effective_components") <- d$effective_components
  attr(ci, "notes") <- c(.mfrm_mvds_component_note(d),
    "Estimate uses the selected component treatment. IntervalEstimate and IntervalPointStatus retain the original source-study projection. Bounds and their availability are unchanged population-parameter intervals; they are not recentered on the adjusted point and need not contain it. Existing IntervalCenter values retain their method-specific meaning.",
    "IntervalSE, SamplingVariance, error-variance diagnostics and other sampling calculations, when supplied, retain their original interval-method basis. They are not sampling uncertainty of the GENOVA-adjusted estimator. Zero replacement is nonsmooth at a boundary; unchanged GENOVA standard errors do not qualify adjusted-estimator inference. In the method notes below, references to the original Estimate or its SE mean IntervalEstimate or IntervalSE for this adjusted-point result.",
    attr(ci, "notes"))
  ci
}

.mfrm_mvds_interval_policy <- function(ci, object) {
  policy <- .mfrm_mvds_policy(object)
  recorded <- attr(ci, "negative_components")
  if (is.null(recorded)) recorded <- "raw" # Legacy original-source intervals.
  if (!identical(recorded, policy))
    stop("Interval component treatment does not match this D-study; recreate intervals with confint(d).", call. = FALSE)
  if (identical(policy, "raw")) {
    if (!is.null(ci$ComponentTreatment) && anyNA(match(ci$ComponentTreatment, "raw")))
      stop("Invalid interval component treatment.", call. = FALSE)
    return(invisible(NULL))
  }
  required <- c("ComponentTreatment", "IntervalBasis", "IntervalEstimate", "IntervalPointStatus")
  if (!all(required %in% names(ci)) ||
      anyNA(ci[c("ComponentTreatment", "IntervalBasis")]) ||
      any(ci$ComponentTreatment != policy) ||
      any(ci$IntervalBasis != "Original source study; population parameter") ||
      !isTRUE(all.equal(attr(ci, "effective_components"), object$effective_components, tolerance = 0)))
    stop("Incomplete or inconsistent adjusted-point interval basis; recreate intervals with confint(d).", call. = FALSE)
  # Check the unadjusted reference without recomputing an interval or refitting.
  raw <- mfrm_d_study(object$gstudy, object$design_grid, object$weights)
  original <- ci
  original$Estimate <- ci$IntervalEstimate
  original$PointStatus <- ci$IntervalPointStatus
  if (!is.null(attr(ci, "reference"))) {
    if (!all(c("IntervalValue", "IntervalReferenceValue") %in% names(ci)) ||
        is.null(attr(ci, "interval_reference")))
      stop("Missing original interval reference values.", call. = FALSE)
    original$Value <- ci$IntervalValue
    original$ReferenceValue <- ci$IntervalReferenceValue
    attr(original, "reference") <- attr(ci, "interval_reference")
  }
  original$ComponentTreatment <- NULL
  attr(original, "negative_components") <- NULL
  .mfrm_mvds_summary_check(original, raw)
  invisible(NULL)
}

# Reconstruct from the saved policy and source; refuse silently stale consumers.
.mfrm_mvds_replay <- function(x) {
  if (!inherits(x, "mfrm_multivariate_d_study"))
    stop("Expected a multivariate D-study result.", call. = FALSE)
  d <- mfrm_d_study(x$gstudy, x$design_grid, x$weights,
    negative_components = .mfrm_mvds_policy(x))
  expected <- d$coefficients
  if (is.null(x$negative_components) && !"ComponentTreatment" %in% names(x$coefficients))
    expected$ComponentTreatment <- NULL
  same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 128 * .Machine$double.eps))
  if (!same(expected, x$coefficients) || !same(d$covariances, x$covariances) ||
      !same(d$component_diagnostics, x$component_diagnostics) ||
      (!is.null(x$negative_components) &&
       (!same(d$effective_components, x$effective_components) ||
        !same(d$effective_component_diagnostics, x$effective_component_diagnostics))))
    stop("D-study results do not match the stored component treatment, components, weights and plans; recreate the result.", call. = FALSE)
  d
}

# Keep the reason with each saved row, including when normal intervals were
# not requested. A retained optimizer estimate is not certified by this table.
mfrm_extended_normal_intervals <- function(object, estimate, se, level,
                                           requested = TRUE, supported = rep(TRUE, length(estimate))) {
  n <- length(estimate)
  if (is.null(se)) se <- rep(NA_real_, n)
  status <- rep(if (requested) "available" else "not_requested", n)
  review <- rep(if (requested)
    "Explicit pointwise normal approximation; numerical availability does not establish nominal coverage." else
    "Normal intervals were not requested. Use the documented confint() target to request this approximation.", n)
  reasons <- c(
    if (!isTRUE(object$checks$NumericalReady)) "Numerical readiness is unresolved; inspect the saved optimization and integration checks.",
    if (!isTRUE(object$checks$InformationPositive)) "Positive information has not been established; inspect the saved information checks.",
    if (isTRUE(object$checks$EstimatedVarianceBoundary)) "The estimated random-effect variance is at zero; regular normal intervals are withheld. More iterations alone do not remove this boundary restriction.",
    if (isTRUE(object$checks$EstimatedPersonVarianceBoundary)) "The estimated ability variance is at zero; regular normal intervals are withheld. More iterations alone do not remove this boundary restriction.")
  if (length(reasons)) {
    status[] <- if (!isTRUE(object$checks$NumericalReady) || !isTRUE(object$checks$InformationPositive))
      "numerical_review_failed" else "estimated_variance_boundary"
    review[] <- paste(c(reasons, "Saved point estimates are retained for review; inspect fit$checks before interpretation."), collapse = " ")
  } else {
    invalid <- !is.finite(estimate) | !is.finite(se) | se < 0
    status[invalid] <- "invalid_estimate_or_se"
    review[invalid] <- "A finite point estimate and finite nonnegative SE are required. Inspect the saved estimate and covariance; missing bounds are not zero uncertainty."
  }
  status[!supported] <- "unsupported_target"
  review[!supported] <- "This normal-interval route does not support population SDs or variances. A retained point or SE does not authorize a regular variance interval."
  lower <- upper <- rep(NA_real_, n)
  eligible <- status == "available"
  margin <- stats::qnorm((1 - level) / 2, lower.tail = FALSE) * se[eligible]
  lower[eligible] <- estimate[eligible] - margin
  upper[eligible] <- estimate[eligible] + margin
  invalid <- eligible & (!is.finite(lower) | !is.finite(upper))
  status[invalid] <- "nonfinite_bounds"
  review[invalid] <- "The normal bounds are outside finite numerical range; the saved point remains, but this interval cannot be reported. Review the estimate and SE scale."
  lower[invalid] <- upper[invalid] <- NA_real_
  data.frame(Lower = lower, Upper = upper, IntervalStatus = status,
    CIEligible = status == "available", InferenceReview = review, row.names = NULL)
}

# One presentation rule for current fits, older fits and saved-result collectors.
mfrm_extended_calibration_table <- function(object, intervals = c("none", "normal"), level = .95) {
  intervals <- match.arg(intervals)
  mfrm_random_rater_interval_level(level)
  tab <- object$calibration_table
  intervals <- mfrm_extended_normal_intervals(object, tab$Estimate, tab[["SE"]], level,
    requested = intervals == "normal", supported = tab$Parameter %in% c("Fixed facet", "Step"))
  tab[names(intervals)] <- intervals
  tab
}

mfrm_extended_calibration_confint <- function(object, level) {
  tab <- mfrm_extended_calibration_table(object, "normal", level)
  tab <- tab[tab$Parameter %in% c("Fixed facet", "Step"), , drop = FALSE]
  bounds <- as.matrix(tab[c("Lower", "Upper")])
  rownames(bounds) <- paste(tab$Parameter, tab$Facet, tab$Level, sep = ": ")
  structure(bounds, level = level, method = "Observed-information normal approximation",
    target = "Fixed-facet and step calibration parameters",
    diagnostics = tab, checks = object$checks,
    note = "Pointwise approximation only; nominal coverage is not established. Missing bounds retain numerical, information, boundary or SE restrictions. No variance, Person-difference or rater-quality inference.")
}

#' Explicit normal-approximation intervals for testlet calibration
#'
#' @param object A result from [fit_mfrm_testlet()].
#' @param parm `"calibration"` returns all fixed-facet and step intervals.
#'   Variance components and Person abilities are not included.
#' @param level Nominal confidence level between zero and one; default 0.95.
#' @param ... Unused.
#' @return A matrix with `Lower` and `Upper`, one row per fixed-facet level or
#'   step, and attributes `level`, `method`, `target`, `note`, `checks` and
#'   `diagnostics`. `attr(ci, "diagnostics")` retains each estimate, SE,
#'   `IntervalStatus`, `CIEligible` and plain-language `InferenceReview`; `checks`
#'   retains the source checks without recomputing them. Save the whole object
#'   with [saveRDS()]; for CSV, export the diagnostics table to retain reasons.
#' @section Reading missing bounds:
#'   `estimated_variance_boundary` means an estimated variance is at zero;
#'   regular normal inference is withheld even when a point is retained.
#'   `numerical_review_failed` means numerical readiness or positive information
#'   has not been established; `InferenceReview` retains all detected reasons,
#'   including any simultaneous boundary. Increasing iterations alone does not
#'   repair a variance-boundary restriction. Inspect `attr(ci, "checks")` first.
#'   `invalid_estimate_or_se` and `nonfinite_bounds` identify row-specific
#'   calculation failures. Missing bounds are neither zero uncertainty nor
#'   infinite endpoints. Retaining a point does not certify its reliability.
#'   Default calibration tables also distinguish `not_requested` from
#'   `unsupported_target` (population SDs/variances). `CIEligible` is true only
#'   for a returned normal interval; it is not a coverage guarantee.
#' @details The explicit request computes estimate plus/minus a normal quantile
#'   times the saved observed-information SE. It needs no fitting, scoring or
#'   live optimizer and does not change the source fit. Nominal finite-sample
#'   coverage is not established. Estimated variance boundaries, unresolved
#'   numerical/information checks and invalid SEs retain missing bounds.
#'   No regular variance interval, simultaneous comparison or automatic rater
#'   classification is supplied. Default fit tables, summaries and plots omit
#'   these bounds. Use `summary(fit, calibration_intervals = "normal")`,
#'   `plot(fit, intervals = "normal")`, or
#'   `mfrm_results(fit, calibration_intervals = "normal")` to select the same
#'   approximation in those outputs; their level arguments retain its nominal
#'   interpretation. [predict.mfrm_testlet()] supplies the separate conditional
#'   Person-scoring intervals.
#' @seealso [confint.mfrm_random_rater()], [mfrmr_interval_guide()]
#' @export
confint.mfrm_testlet <- function(object, parm = "calibration", level = .95, ...) {
  rlang::check_dots_empty()
  if (!identical(parm, "calibration")) stop("Use parm = 'calibration'; testlet-variance and ability-SD intervals are not supplied.", call. = FALSE)
  mfrm_extended_calibration_confint(object, level)
}

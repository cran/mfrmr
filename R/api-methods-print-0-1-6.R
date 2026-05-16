# Lightweight print delegates added in 0.1.6 to close the first-use surface.
#
# NEWS 0.1.5 committed to `Status`, `Key warnings`, and `Next actions`
# leading every first-use print, but 13 classes shipped only `summary()`
# methods. Without a matching `print()` S3, typing the object at the R
# console fell back to the default list printer and dumped hundreds of
# rows of raw tables.
#
# The following delegates route `print(x)` through `summary(x)`, so the
# curated first-use output defined in the summary methods becomes the
# console default. Callers that want the raw list can still use
# `unclass(x)` or direct subsetting.

#' @export
print.mfrm_apa_outputs <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_bias <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_bundle <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_design_evaluation <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_diagnostics <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_facet_dashboard <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_future_branch_active_branch <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_plausible_values <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_population_prediction <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_reporting_checklist <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_signal_detection <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_threshold_profiles <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @export
print.mfrm_unit_prediction <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

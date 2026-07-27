# Lightweight print delegates for concise first-use output.
#
# These classes provide curated `summary()` methods. Matching `print()`
# delegates prevent the default list printer from dumping large raw tables
# when a user types an object at the R console.
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
print.mfrm_structural_design_review <- function(x, ...) {
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

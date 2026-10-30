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

#' @rawNamespace S3method(print,mfrm_apa_outputs,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_bias,print.mfrm_bundle)
NULL

#' @export
print.mfrm_bundle <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @rawNamespace S3method(print,mfrm_design_evaluation,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_diagnostics,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_facet_dashboard,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_structural_design_review,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_plausible_values,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_population_prediction,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_reporting_checklist,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_signal_detection,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_threshold_profiles,print.mfrm_bundle)
NULL

#' @rawNamespace S3method(print,mfrm_unit_prediction,print.mfrm_bundle)
NULL

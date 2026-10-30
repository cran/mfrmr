# Saved Person tables are keyed by identifier, never by their display order.
mfrm_saved_person_estimates <- function(object) {
  person <- object$facets$person
  if (!is.data.frame(person) || !all(c("Person", "Estimate") %in% names(person)) ||
      anyNA(person$Person) || anyDuplicated(person$Person))
    stop("Saved primary Person estimates need distinct, nonmissing identifiers.", call. = FALSE)
  person$Estimate[match(object$prep$levels$Person, person$Person)]
}

# Legacy diagnostics deliberately retain optimizer-coordinate JML calculations
# and their variance screening threshold. Keep that basis visible in row exports.
mfrm_annotate_diagnostic_responses <- function(out, res, idx, bundle) {
  primary <- mfrm_saved_person_estimates(res)[idx$person]
  variance <- bundle$conditional_var_k
  available <- is.finite(out$StdResidual)
  valid <- is.finite(out$Expected) & is.finite(variance) & variance >= 0
  small <- valid & variance > 0 & variance <= 1e-10
  zero <- valid & variance == 0
  out$PrimaryPersonMeasure <- primary
  out$PersonMeasureBasis <- if (res$config$method == "JML") "JML optimizer" else "saved EAP"
  out$ConditionalVariance <- variance
  out$PearsonAvailable <- available
  out$ResponseStatus <- ifelse(available, "available_plugin", "unavailable")
  out$ResponseStatus[small] <- "small_variance"
  out$ResponseStatus[zero] <- "zero_variance"
  out$ResponseReason <- ifelse(available,
    "Conditional point calculation; calibration and Person uncertainty are excluded.",
    "Conditional probabilities or the standardized residual are unavailable at the diagnostic estimates.")
  out$ResponseReason[small] <- paste(
    "Positive conditional variance is at or below the legacy diagnostic cutoff (1e-10).",
    "Expected ratings and raw residuals remain available; standardized diagnostics are withheld.")
  out$ResponseReason[zero] <- paste(
    "Conditional variance is zero; expected ratings and raw residuals remain available,",
    "but the Pearson residual is undefined.")
  if (res$config$method == "JML") {
    proxy <- !is.finite(primary) | abs(primary - out$PersonMeasure) > 1e-10
    out$ResponseReason[proxy] <- paste(
      "JML diagnostic probabilities use optimizer coordinates, which differ from the primary Person estimate.",
      out$ResponseReason[proxy])
  }
  out
}

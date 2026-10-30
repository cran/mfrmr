#' Extract coefficients, covariance, likelihood and observation counts
#'
#' Standard R accessors for native many-facet fits. They retain the estimator's
#' coordinates and limitations; extracting a number does not qualify inference.
#'
#' @param object A native result from [fit_mfrm()]. Imported and extended-model
#'   results use their documented summaries and model-specific methods.
#' @param type For `coef()`, `"free"` (default) returns the independent fitted
#'   coordinates; `"locations"`, `"steps"`, `"slopes"` or `"persons"` returns
#'   the corresponding saved primary estimates. For `nobs()`, `"likelihood"`
#'   (default) counts Persons for MML/CML and retained response rows for JML;
#'   `"persons"` and `"ratings"` request those counts explicitly.
#' @param k Penalty per free parameter for single-fit `AIC()`; default 2.
#' @param ... Unused for coefficient, covariance, likelihood and count
#'   extraction; misspelled arguments are rejected. Multiple-model calls to
#'   `AIC()`/`BIC()` are refused: use [compare_mfrm()] for supported comparisons
#'   with matched observations, estimation and integration settings.
#'
#' @details
#' **Coordinates.** The default `coef(fit)` is a named numeric vector aligned
#' with both dimensions of `vcov(fit)`. These are free coordinates, not every
#' displayed facet level: anchors and dependent centered coordinates are omitted.
#' MML/JML names such as `Rater[1]` identify positions within parameter blocks,
#' not literal level names. Steps are free step offsets, slopes are log slopes,
#' and estimated population variance is on the log-variance scale. Population
#' regression coefficients remain on their fitted design scale. CML retains its
#' existing conditional-coordinate names. Use `type = "locations"` for all
#' saved facet locations with their level labels, and `summary(fit)` for tables,
#' readiness, constraints and interpretation. Location signs are the fitted
#' model's signs; this accessor does not reorient them as severity.
#' Default [confint()][confint.mfrm_fit] intervals use these same coordinates;
#' GPCM positive-slope intervals require explicit `parm = "slopes"`.
#'
#' `type = "persons"` returns saved primary Person estimates: MML EAPs and
#' JML point estimates are different targets. Person posterior SDs are not
#' structural standard errors and do not occur in `vcov(fit)`. CML has no
#' Person estimates. Free-coordinate extraction refuses corrected or penalized
#' JML and any fit whose primary estimates are nonfinite: finite optimizer
#' proxies must not replace primary boundary values. Explicit component
#' extraction preserves the saved infinite or missing primary estimates.
#'
#' **Covariance.** MML `vcov()` reuses the same qualified full joint covariance
#' as [mfrm_facet_intervals()], including its model/design restrictions and
#' numerical checks. It may reevaluate information and integration, but does
#' not refit. An ineligible source remains an error. For an otherwise eligible
#' ordinary fixed-population RSM/PCM fit, information failure returns an `NA`
#' matrix with `unavailable_reason`, `checks` and `cautions`; default `confint()`
#' retains coefficient estimates and reasons in its `diagnostics` attribute
#' while bounds are missing. No regularized inverse is substituted. If even
#' one dense `NA` matrix exceeds `mfrmr.max_information_bytes`, `vcov()` errors
#' with that reason; use the smaller target-level [mfrm_facet_intervals()] result.
#' Other MML routes retain their separate qualification checks and error rules.
#' RSM/PCM singleton-category fits can pass the interval-specific information
#' review while retaining `InferenceReady = FALSE`; returned `cautions` and
#' checks explain the small-sample limitation. No count label alone certifies
#' or invalidates this covariance.
#' No JML structural-truth covariance is supplied. CML `vcov()` reuses saved
#' conditional information. If CML's finite numerical point passed its
#' identification, finite-MLE and stationarity checks but information failed,
#' `vcov()` returns an `NA` matrix with `unavailable_reason`, `checks` and
#' `cautions` attributes. `coef()` and `logLik()` still expose the retained
#' values; default `confint()` supplies missing bounds with diagnostic reasons.
#' The warning also limits confidence in point accuracy along weak directions;
#' no regularized covariance or formal inference is substituted.
#' Free-coordinate covariance is not the covariance
#' of a selected component vector or of all displayed levels; use
#' [mfrm_facet_intervals()] or [confint.mfrm_fit()] for supported transformations.
#' Comparison-only `coef()` returns the declared contrast estimates, retaining
#' unavailable rows as `NA`; its `vcov()` is available only for native CML
#' conditional comparisons, with the same row/column order.
#'
#' **Likelihood and counts.** `logLik()` returns the retained MML marginal or
#' ordinary JML joint log likelihood, with `likelihood_basis` and `ICReason`
#' attributes. Nonunit weights give a weighted objective, not an ordinary
#' sampling likelihood. The `df` attribute is the free parameter count only
#' when the current information-criterion contract and solution checks permit
#' its use; otherwise it is `NA`, so `AIC()` and `BIC()` are also `NA`.
#' `parameter_count` retains the actual dimension separately. JML and two-family
#' GPCM do not gain information-criterion eligibility through these accessors.
#' Corrected JML, penalized objectives and comparison-only fits do not supply
#' a generic `logLik()` here. For estimated-population/one-family GPCM MML,
#' likelihood qualification may reevaluate solution information.
#'
#' CML retains its conditional log likelihood and free dimension. Its
#' `nobs` attribute counts all retained Persons, including extreme totals
#' contributing zero conditional information; `informative_persons` records
#' nonextreme totals separately. This preserves the existing CML method and
#' is not a new validation of conditional information criteria. Conditional
#' and marginal likelihoods must not be compared. For MML, BIC uses Persons,
#' not rating rows or the sum of rating weights. Counts describe the retained
#' analysis, not omitted input records or an effective sample size under
#' dependence. Use `nobs(fit, type = "ratings")` to count retained ratings.
#'
#' @return `coef()` returns a named numeric vector; `vcov()` a matrix in the
#'   matching default coordinates. MML covariance retains `checks`,
#'   `numerical_checks`, `information_review` and `cautions` attributes when
#'   supplied by the interval route. `logLik()` returns a `logLik` object; `nobs()` an
#'   integer count; single-fit `AIC()`/`BIC()` a numeric value or `NA` when
#'   unsupported. Existing CML covariance coordinates and numeric values are
#'   unchanged. No accessor silently computes new Person scores.
#' @seealso [fit_mfrm()], [summary.mfrm_fit()], [mfrm_fit_responses], [mfrm_facet_intervals()],
#'   [compare_mfrm()], [mfrmr_compatibility_layer]
#' @examples
#' ratings <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", quad_points = 31, maxit = 500)
#' coef(fit, type = "locations") # Labelled facet locations, including constraints
#' coef(fit)                     # Independent coordinates, including steps
#' nobs(fit)                     # Persons contributing to marginal likelihood
#' nobs(fit, type = "ratings")    # Retained response rows
#' ll <- logLik(fit)
#' attr(ll, "ICReason")          # Read before interpreting AIC/BIC
#' AIC(fit)
#' # When covariance is qualified, its names match the free coefficients:
#' # V <- vcov(fit); identical(rownames(V), names(coef(fit)))
#' @name mfrm_fit_extractors
#' @importFrom stats coef vcov logLik nobs AIC BIC
NULL

mfrm_extractor_native <- function(object) {
  if (inherits(object, "mfrm_imported_fit") ||
      !inherits(object, c("mfrm_fit", "mfrm_cml_fit", "mfrm_comparison_fit")))
    stop("These extractors require a native fit_mfrm() result; use the source model's summary and methods.", call. = FALSE)
}

#' @rdname mfrm_fit_extractors
#' @export
coef.mfrm_fit <- function(object, type = c("free", "locations", "steps", "slopes", "persons"), ...) {
  rlang::check_dots_empty()
  mfrm_extractor_native(object)
  type <- match.arg(type)
  if (type == "free") {
    stop_if_jml_adjustment(object, "coef(type = 'free')")
    if (!identical(object$config$facet_shrinkage %||% "none", "none"))
      stop("Penalized coefficients require their model-specific summary; use explicit primary component extraction.", call. = FALSE)
    primary <- c(object$facets$others[["Estimate"]], object$facets$person[["Estimate"]],
      object$steps[["Estimate"]], object$slopes[["Estimate"]])
    if (any(!is.finite(primary)))
      stop("Nonfinite primary estimates cannot be replaced with optimizer proxies. Use coef(type = 'locations'), 'persons', 'steps' or 'slopes', and inspect summary().", call. = FALSE)
    values <- object$opt$par
    labels <- if (inherits(object, "mfrm_cml_fit")) names(values) else
      mfrm_checkpoint_parameter_names(build_param_sizes(object$config))
  } else {
    tab <- switch(type, locations = object$facets$others, steps = object$steps,
      slopes = object$slopes, persons = object$facets$person)
    if (!is.data.frame(tab) || !"Estimate" %in% names(tab))
      stop("This fit has no saved ", type, " estimates; inspect its model-specific summary.", call. = FALSE)
    keys <- switch(type, locations = c("Facet", "Level"), steps = c("StepFacet", "Step"),
      slopes = intersect(c("SlopeOwner", "SlopeFacet"), names(tab)), persons = "Person")
    if (!length(keys) || !all(keys %in% names(tab))) stop("Saved coefficient labels are incomplete.", call. = FALSE)
    labels <- do.call(paste, c(lapply(tab[keys], function(x) encodeString(as.character(x), quote = '"')), sep = ":"))
    values <- tab$Estimate
  }
  if (!is.numeric(values) || length(labels) != length(values) || anyNA(labels) ||
      anyDuplicated(labels) || any(!nzchar(labels)))
    stop("Saved coefficient coordinates are incomplete or ambiguous.", call. = FALSE)
  stats::setNames(as.numeric(values), labels)
}

#' @rdname mfrm_fit_extractors
#' @name coef.mfrm_cml_fit
#' @aliases coef.mfrm_cml_fit
#' @usage \method{coef}{mfrm_cml_fit}(object, type = c("free", "locations", "steps", "slopes", "persons"), ...)
#' @rawNamespace S3method(coef,mfrm_cml_fit,coef.mfrm_fit)
NULL

#' @rdname mfrm_fit_extractors
#' @export
coef.mfrm_comparison_fit <- function(object, ...) {
  rlang::check_dots_empty()
  stats::setNames(as.numeric(object$comparisons$Estimate), as.character(object$comparisons$Comparison))
}

#' @rdname mfrm_fit_extractors
#' @export
vcov.mfrm_fit <- function(object, ...) {
  rlang::check_dots_empty()
  mfrm_extractor_native(object)
  if (!identical(object$config$method, "MML"))
    stop("Structural covariance requires MML or native CML; JML point estimates do not supply structural-truth Wald inference.", call. = FALSE)
  coefficients <- coef(object)
  facet <- object$config$facet_names[1L]
  ci <- mfrm_facet_intervals(object, facet)
  value <- ci$parameter_covariance
  if (is.null(value) && is.character(ci$unavailable_reason) &&
      length(ci$unavailable_reason) == 1L && !is.na(ci$unavailable_reason) &&
      nzchar(ci$unavailable_reason)) {
    # Even an NA matrix needs quadratic storage; respect the information budget.
    if (8 * as.double(length(coefficients))^2 > getOption("mfrmr.max_information_bytes", 256 * 1024^2))
      stop(paste("The requested NA covariance matrix exceeds option mfrmr.max_information_bytes.",
        ci$unavailable_reason, "Use mfrm_facet_intervals() for target-level results."), call. = FALSE)
    value <- matrix(NA_real_, length(coefficients), length(coefficients))
    attr(value, "unavailable_reason") <- ci$unavailable_reason
  } else if (!is.matrix(value) || !identical(dim(value), rep(length(coefficients), 2L)) ||
      any(!is.finite(value)) || (!is.null(ci$checks) && !all(ci$checks$Passed %in% TRUE))) {
    reasons <- ci$checks$Detail[!ci$checks$Passed %in% TRUE]
    stop(paste("Joint covariance is unavailable.", paste(unique(reasons), collapse = " "),
      "Inspect mfrm_facet_intervals() checks and inference reasons."), call. = FALSE)
  }
  dimnames(value) <- list(names(coefficients), names(coefficients))
  attr(value, "checks") <- ci$checks
  attr(value, "numerical_checks") <- ci$numerical_checks
  attr(value, "information_review") <- ci$information_review
  attr(value, "cautions") <- ci$cautions
  value
}

#' @rdname mfrm_fit_extractors
#' @export
vcov.mfrm_comparison_fit <- function(object, ...) {
  rlang::check_dots_empty()
  if (!identical(object$config$method, "CML") || !is.matrix(object$covariance))
    stop("Comparison-only covariance is available for native CML targets; other comparison-only fits provide descriptive estimates.", call. = FALSE)
  value <- object$covariance
  if (identical(object$optimization$checks$information_status, "unavailable")) {
    attr(value, "unavailable_reason") <- object$optimization$checks$information_reason
    attr(value, "checks") <- object$checks
    attr(value, "cautions") <- object$notes
  }
  value
}

mfrm_coefficient_intervals <- function(object, parm, level) {
  mfrm_extractor_native(object)
  if (mfrm_has_product_slopes(object))
    stop("Free-coordinate intervals are not supplied for two-family GPCM; use parm = 'slopes' for its separately checked experimental component intervals.", call. = FALSE)
  b <- coef(object)
  labels <- names(b)
  if (!length(b) || length(labels) != length(b) || anyNA(labels) ||
      anyDuplicated(labels) || any(!nzchar(labels)))
    stop("Saved coefficient coordinates are incomplete or ambiguous.", call. = FALSE)
  selected <- seq_along(b)
  if (!is.null(parm)) {
    if (is.object(parm) || !is.null(dim(parm)))
      stop("`parm` must select coefficient names or integer positions.", call. = FALSE)
    if (is.character(parm)) selected <- match(parm, labels) else {
      if (!is.numeric(parm) || is.complex(parm) || any(!is.finite(parm)) || any(parm != floor(parm)))
        stop("`parm` must select coefficient names or integer positions.", call. = FALSE)
      selected <- parm
    }
    if (!length(selected) || anyNA(selected) || anyDuplicated(selected) ||
        any(selected < 1 | selected > length(b)))
      stop("`parm` must select distinct existing coefficients; use names(coef(object)). For GPCM slope intervals, use parm = 'slopes'.", call. = FALSE)
  }
  V <- vcov(object)
  if (!is.matrix(V) || !is.numeric(V) || !identical(dim(V), rep(length(b), 2L)) ||
      !setequal(rownames(V), labels) || !setequal(colnames(V), labels) ||
      anyDuplicated(rownames(V)) || anyDuplicated(colnames(V)))
    stop("Covariance dimensions and names must match the saved coefficient coordinates.", call. = FALSE)
  covariance_notes <- attributes(V)
  covariance_notes <- covariance_notes[intersect(c("checks", "numerical_checks", "information_review", "cautions", "unavailable_reason"), names(covariance_notes))]
  V <- V[labels, labels, drop = FALSE]
  variance <- diag(V)[selected]
  b <- b[selected]
  available <- is.finite(b) & is.finite(variance) & variance > 0
  reason <- rep("Pointwise normal approximation in the coef() coordinates; numerical availability does not guarantee finite-sample coverage.", length(b))
  reason[!available] <- "A finite coefficient and positive finite sampling variance are required; no interval is supplied."
  if (length(covariance_notes$unavailable_reason) == 1L &&
      !is.na(covariance_notes$unavailable_reason) && nzchar(covariance_notes$unavailable_reason))
    reason[!available] <- paste("Point estimate retained; uncertainty is unavailable:",
      covariance_notes$unavailable_reason)
  comparison <- inherits(object, "mfrm_comparison_fit")
  if (comparison) {
    flags <- object$comparisons[["FormalInferenceEligible"]]
    allowed <- if (is.logical(flags) && length(flags) == length(labels))
      flags[selected] %in% TRUE else rep(FALSE, length(selected))
    review <- rep("Saved conditional-target eligibility is unavailable; recreate the comparison fit before interval reporting.", length(selected))
    saved_review <- object$comparisons[["EstimateReason"]]
    if (is.character(saved_review) && length(saved_review) == length(labels)) {
      stated <- !is.na(saved_review[selected]) & nzchar(saved_review[selected])
      if (is.logical(flags) && length(flags) == length(labels)) review[stated] <- saved_review[selected][stated]
    }
    reason[!allowed] <- review[!allowed]
    available <- available & allowed
  }
  se <- sqrt(ifelse(is.finite(variance) & variance >= 0, variance, NA_real_))
  margin <- stats::qnorm((1 - level) / 2, lower.tail = FALSE) * se
  out <- cbind(b - margin, b + margin)
  representable <- rowSums(!is.finite(out)) == 0
  reason[available & !representable] <- "Normal interval bounds are not finite on these coefficient coordinates."
  available <- available & representable
  out[!available, ] <- NA_real_
  dimnames(out) <- list(names(b), paste0(format(100 * c((1 - level) / 2, (1 + level) / 2), trim = TRUE), "%"))
  attr(out, "level") <- level
  attr(out, "target") <- if (comparison) "Saved conditional comparisons" else "Free model coordinates, matching coef()"
  attr(out, "method") <- if (object$config$method == "CML") "Conditional-information pointwise Wald" else "Joint MML information pointwise Wald"
  attr(out, "diagnostics") <- data.frame(Parameter = names(b), Estimate = unname(b), SE = unname(se),
    CIEligible = available, Status = ifelse(available, "available", "unavailable"),
    InferenceReview = reason, row.names = NULL)
  for (key in names(covariance_notes)) if (!is.null(covariance_notes[[key]]))
    attr(out, key) <- covariance_notes[[key]]
  if (object$config$method == "CML") {
    attr(out, "checks") <- object$checks
    attr(out, "cautions") <- object$notes
  }
  out
}

#' @rdname mfrm_fit_extractors
#' @export
nobs.mfrm_fit <- function(object, type = c("likelihood", "persons", "ratings"), ...) {
  rlang::check_dots_empty()
  mfrm_extractor_native(object)
  type <- match.arg(type)
  if (type == "likelihood") type <- if (object$config$method %in% c("MML", "CML")) "persons" else "ratings"
  data <- object$prep$data
  if (!is.data.frame(data) || !"Person" %in% names(data) || anyNA(data$Person))
    stop("The retained response roster is unavailable or incomplete.", call. = FALSE)
  if (type == "persons") length(unique(as.character(data$Person))) else nrow(data)
}

#' @rdname mfrm_fit_extractors
#' @name nobs.mfrm_cml_fit
#' @aliases nobs.mfrm_cml_fit
#' @usage \method{nobs}{mfrm_cml_fit}(object, type = c("likelihood", "persons", "ratings"), ...)
#' @rawNamespace S3method(nobs,mfrm_cml_fit,nobs.mfrm_fit)
NULL

#' @rdname mfrm_fit_extractors
#' @name nobs.mfrm_comparison_fit
#' @aliases nobs.mfrm_comparison_fit
#' @usage \method{nobs}{mfrm_comparison_fit}(object, type = c("likelihood", "persons", "ratings"), ...)
#' @rawNamespace S3method(nobs,mfrm_comparison_fit,nobs.mfrm_fit)
NULL

#' @rdname mfrm_fit_extractors
#' @export
logLik.mfrm_fit <- function(object, ...) {
  rlang::check_dots_empty()
  mfrm_extractor_native(object)
  stop_if_jml_adjustment(object, "logLik()")
  if (!identical(object$config$facet_shrinkage %||% "none", "none"))
    stop("A penalized objective is not exposed as an ordinary log likelihood; use the model-specific summary.", call. = FALSE)
  contract <- mfrm_extract_fit_ic_contract(object)
  if (!isTRUE(contract$StoredICConsistent) || !is.finite(contract$LogLik))
    stop("The retained likelihood metadata are missing or inconsistent; inspect the source fit rather than reconstructing a likelihood from legacy summaries.", call. = FALSE)
  eligible <- isTRUE(contract$ICSelectable)
  reason <- switch(as.character(contract$ICStatus),
    descriptive_jml = "JML retains a joint objective for descriptive use; structural model selection by AIC/BIC is unsupported.",
    unsupported_product_slopes = "Information criteria are not qualified for two-family GPCM.",
    suppressed_nonunit_weight = "Nonunit rating weights define a weighted objective; ordinary likelihood information criteria are unavailable.",
    "The retained likelihood contract does not permit information criteria.")
  if (eligible) {
    check <- mfrm_ic_fit_check(object, contract)
    eligible <- isTRUE(check$eligible)
    reason <- if (eligible) "Current MML likelihood and solution checks permit information criteria; model comparisons require matched sources." else check$review
  } else if (isTRUE(contract$ICEligible)) reason <- "The retained quadrature grid is not qualified for information-criterion selection."
  basis <- if (identical(object$config$method, "MML")) "marginal" else "joint"
  if (!contract$WeightPolicy %in% c("unweighted", "explicit_unit")) basis <- paste0("weighted_", basis, "_objective")
  structure(as.numeric(contract$LogLik), df = if (eligible) as.integer(contract$Npar) else NA_integer_,
    nobs = nobs(object), parameter_count = as.integer(contract$Npar), likelihood_basis = basis,
    ICEligible = eligible, ICReason = reason, class = "logLik")
}

#' @rdname mfrm_fit_extractors
#' @export
logLik.mfrm_comparison_fit <- function(object, ...) {
  rlang::check_dots_empty()
  stop("Comparison-only fits do not supply a generic full-model log likelihood or information criteria; use summary().", call. = FALSE)
}

#' @rdname mfrm_fit_extractors
#' @export
AIC.mfrm_fit <- function(object, ..., k = 2) {
  if (length(list(...))) stop("Use compare_mfrm() for supported model comparisons; conditional and marginal likelihoods must not be mixed.", call. = FALSE)
  if (!is.numeric(k) || is.complex(k) || length(k) != 1L || !is.finite(k))
    stop("`k` must be one finite numeric penalty.", call. = FALSE)
  stats::AIC(logLik(object), k = k)
}

#' @rdname mfrm_fit_extractors
#' @export
BIC.mfrm_fit <- function(object, ...) {
  if (length(list(...))) stop("Use compare_mfrm() for supported model comparisons; conditional and marginal likelihoods must not be mixed.", call. = FALSE)
  stats::BIC(logLik(object))
}

#' @rdname mfrm_fit_extractors
#' @name AIC.mfrm_cml_fit
#' @aliases AIC.mfrm_cml_fit
#' @usage \method{AIC}{mfrm_cml_fit}(object, ..., k = 2)
#' @rawNamespace S3method(AIC,mfrm_cml_fit,AIC.mfrm_fit)
NULL

#' @rdname mfrm_fit_extractors
#' @name BIC.mfrm_cml_fit
#' @aliases BIC.mfrm_cml_fit
#' @usage \method{BIC}{mfrm_cml_fit}(object, ...)
#' @rawNamespace S3method(BIC,mfrm_cml_fit,BIC.mfrm_fit)
NULL

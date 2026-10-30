#' Extract fitted rating values and residuals
#'
#' Describe the retained ratings at the fitted calibration and saved Person
#' estimates, without fitting, scoring Persons again or computing covariance.
#'
#' @param object A native ordinary RSM/PCM or one-family GPCM fit from
#'   [fit_mfrm()], estimated by MML or uncorrected JML.
#' @param type For `fitted()` and `predict()`, `"response"` returns expected
#'   ratings and `"prob"` category probabilities. For `residuals()`,
#'   `"response"` returns observed minus expected ratings and `"pearson"`
#'   divides that difference by the conditional rating SD.
#' @param newdata Must be `NULL`. These methods describe retained source
#'   responses only. For ability scoring from supplied responses use
#'   [predict_mfrm_units()] or [score_mfrm_calibration()]; those return Person
#'   abilities, not predicted rating values. [score_mfrm_persons()] supplies
#'   conditional source-Person scoring for its documented model scope.
#' @param ... Unused; additional arguments are rejected. Prediction intervals,
#'   standard errors, new rating designs and missing-score imputation are not
#'   supplied by these methods.
#'
#' @details
#' **Target.** Probabilities are evaluated at saved primary Person estimates:
#' MML uses the saved EAP, while JML uses its primary point estimate. Thus
#' `fitted(fit)` is a conditional plug-in expectation. Evaluating probabilities
#' at an EAP is different from averaging probabilities over a posterior.
#' For the latter target, use [mfrm_response_diagnostics()] within its stated
#' model/design scope. These methods do not change [diagnose_mfrm()] or its
#' legacy diagnostic conventions. Legacy diagnostic helpers use finite JML
#' optimizer proxies and omit variances at or below `1e-10`; the primary
#' boundary limits and small positive variances here can therefore give
#' different residuals. Their observation tables retain `PersonMeasureBasis`,
#' `PrimaryPersonMeasure`, `ConditionalVariance` and `ResponseReason` to explain
#' these differences. Do not assume the routes agree at these limits.
#'
#' For model-score categories \eqn{k}, let \eqn{p_{ik}} be the conditional
#' probability. The mean is \eqn{\mu_i=\sum_k k p_{ik}}, the variance is
#' \eqn{V_i=\sum_k (k-\mu_i)^2p_{ik}}, and the Pearson residual is
#' \eqn{(y_i-\mu_i)/\sqrt{V_i}}. Rating weights affect estimation but are
#' not multiplied into these per-rating residuals or probabilities. The
#' variance has no artificial positive floor. A zero variance leaves the
#' mean/raw residual available but the Pearson residual missing, with a reason.
#'
#' A positive/negative infinite JML Person estimate with finite structural
#' calibration has limiting probability one at the maximum/minimum category.
#' These limits use the primary estimate, not the finite optimizer proxy.
#' Other unavailable Person estimates leave the corresponding rows missing.
#' Nonfinite structural estimates, corrected or penalized fits, two-family
#' GPCM, imported fits, CML and comparison-only fits are explicitly refused.
#' CML eliminates Person abilities and does not supply this plug-in target.
#' The supported posterior or corrected-JML routes remain available through
#' [mfrm_response_diagnostics()].
#'
#' **Rows and scale.** Output order is `fit$prep$data` order, including repeated
#' rating events, and excludes records omitted before fitting. Vector names
#' and probability-matrix row names are retained-row positions. The `rows`
#' attribute records these positions and original input-row
#' positions where retained metadata permit recovery; unknown positions are
#' `NA`. No response is silently added or dropped during extraction.
#' The separate `identifiers` attribute retains Person/facet columns in the
#' same order, including user-supplied facet names that coincide with `Row`.
#'
#' Values use the modeled score scale (`InternalScore` in `score_map`), as in
#' ordinary plug-in diagnostics. With collapsed categories this differs from
#' the original score labels. The saved `score_map` attribute retains that
#' mapping; probability columns name the modeled categories. Inspect the map
#' before comparing fitted values with an original, unrecoded data column.
#'
#' **Interpretation.** These same-data summaries hold calibration and Person
#' estimates fixed. They omit their estimation uncertainty and are not held-out
#' forecasts, calibrated goodness-of-fit tests or evidence of rater quality.
#' Point calculation does not establish convergence or inference readiness.
#' Read `source_readiness` and row-specific `diagnostics` before interpreting
#' values, and keep these attributes when saving the result.
#'
#' @return A named numeric vector, or for `type = "prob"` a numeric matrix
#'   with one column per modeled category. Attributes include `rows`, `identifiers`,
#'   `score_map`, `diagnostics` (availability, conditional variance and reason),
#'   `settings` (quantity, target, Person basis and limitations) and `source_readiness`.
#'   `predict(fit)` and `fitted(fit)` return the same source-rating target.
#' @seealso [mfrm_fit_extractors], [diagnose_mfrm()], [mfrm_response_diagnostics()],
#'   [predict_mfrm_units()], [mfrmr_compatibility_layer]
#' @examples
#' ratings <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'   method = "MML", quad_points = 31, maxit = 500)
#' expected <- fitted(fit) # Expected ratings at saved EAP abilities
#' head(data.frame(Observed = fit$prep$data$Score, Expected = as.numeric(expected),
#'   Residual = as.numeric(residuals(fit))))
#' head(attr(expected, "diagnostics")) # Check unavailable and zero-variance rows
#' attr(expected, "settings")         # Read the target and interpretation limits
#' head(fitted(fit, type = "prob"))    # One column per modeled score category
#' @name mfrm_fit_responses
#' @importFrom stats fitted residuals predict
NULL

mfrm_fit_response_values <- function(object, type, residual = FALSE) {
  if (inherits(object, c("mfrm_cml_fit", "mfrm_comparison_fit")))
    stop("CML/comparison-only fits do not supply primary Person estimates for plug-in rating predictions; use their documented conditional summaries.", call. = FALSE)
  mfrm_extractor_native(object)
  if (mfrm_has_jml_adjustment(object) || mfrm_has_product_slopes(object))
    stop("Use mfrm_response_diagnostics() for supported corrected-JML or two-family response calculations; their probability definitions differ from ordinary plug-in extraction.", call. = FALSE)
  cfg <- object$config; prep <- object$prep
  if (!cfg$method %in% c("MML", "JML") || !cfg$model %in% c("RSM", "PCM", "GPCM") ||
      !identical(cfg$facet_shrinkage %||% "none", "none"))
    stop("Plug-in rating extraction requires an ordinary unpenalized MML or JML fit.", call. = FALSE)
  structural <- c(object$facets$others[["Estimate"]], object$steps[["Estimate"]],
    object$slopes[["Estimate"]])
  if (any(!is.finite(structural)) || any(!is.finite(object$opt$par)))
    stop("Finite structural calibration is required; nonfinite primary estimates cannot be replaced by optimizer proxies.", call. = FALSE)
  if (cfg$model == "GPCM" && (!length(object$slopes[["Estimate"]]) ||
      any(object$slopes[["Estimate"]] <= 0)))
    stop("GPCM plug-in responses require strictly positive primary slopes; joint slope/ability boundary limits are not supplied.", call. = FALSE)
  theta <- mfrm_saved_person_estimates(object)
  idx <- build_indices(prep, step_facet = cfg$step_facet, slope_facet = cfg$slope_facet,
    interaction_specs = cfg$interaction_specs)
  params <- expand_params(object$opt$par, build_param_sizes(cfg), cfg)
  # Only finite stand-ins reach the kernel; their rows are replaced below.
  finite_theta <- theta; finite_theta[!is.finite(finite_theta)] <- 0
  eta <- compute_eta(idx, params, cfg, theta_override = finite_theta)
  p <- compute_response_probability_bundle(cfg, idx, params, eta)$probs
  row_theta <- theta[idx$person]
  boundary <- cfg$method == "JML" & is.infinite(row_theta)
  unavailable <- !is.finite(row_theta) & !boundary
  p[unavailable, ] <- NA_real_
  if (any(boundary)) {
    p[boundary, ] <- 0
    p[cbind(which(boundary), ifelse(row_theta[boundary] > 0, ncol(p), 1L))] <- 1
  }
  valid <- rowSums(!is.finite(p) | p < 0 | p > 1) == 0 & abs(rowSums(p) - 1) <= 1e-10
  valid[is.na(valid)] <- FALSE
  p[!valid, ] <- NA_real_
  categories <- prep$rating_min + seq_len(ncol(p)) - 1L
  mu <- drop(p %*% categories)
  variance <- mfrm_category_variance(p, categories, mu)
  raw <- prep$data$Score - mu
  pearson <- rep(NA_real_, length(mu))
  positive <- valid & is.finite(variance) & variance > 0
  pearson[positive] <- raw[positive] / sqrt(variance[positive])
  pearson[!is.finite(pearson)] <- NA_real_
  diagnostics <- data.frame(Row = seq_along(mu), ProbabilityAvailable = valid,
    PearsonAvailable = is.finite(pearson), Variance = variance,
    Status = ifelse(valid, "available_plugin", "unavailable"),
    Reason = ifelse(valid, "Conditional point calculation; calibration and Person uncertainty are excluded.",
      "Conditional probabilities are unavailable at the saved primary estimates."))
  diagnostics$Reason[unavailable] <- "A finite primary Person estimate is unavailable; no optimizer or display proxy was substituted."
  zero <- valid & is.finite(variance) & variance == 0
  diagnostics$Status[zero] <- "zero_variance"
  diagnostics$Reason[zero] <- "Conditional variance is zero; the mean and raw residual remain available, but the Pearson residual is undefined."
  diagnostics$Reason[boundary] <- paste("Limiting probabilities at an infinite JML Person estimate.", diagnostics$Reason[boundary])
  ids <- as.character(seq_along(mu))
  dimnames(p) <- list(ids, as.character(categories))
  out <- if (residual) if (type == "pearson") pearson else raw else if (type == "prob") p else mu
  if (!is.matrix(out)) names(out) <- ids
  rows <- data.frame(Row = seq_along(mu), InputRow = NA_integer_)
  n_input <- prep$row_retention$Rows[match("input_selected_columns", prep$row_retention$Stage)]
  omitted <- prep$omitted_input_rows %||% integer()
  if (length(n_input) == 1L && is.finite(n_input) && n_input >= nrow(rows) &&
      is.numeric(omitted) && !is.complex(omitted) && !anyNA(omitted) && !anyDuplicated(omitted) &&
      all(omitted >= 1 & omitted <= n_input & omitted == floor(omitted)) &&
      n_input - length(omitted) == nrow(rows))
    rows$InputRow <- setdiff(seq_len(n_input), omitted)
  attr(out, "rows") <- rows
  attr(out, "identifiers") <- as.data.frame(prep$data[c("Person", cfg$facet_names)])
  attr(out, "score_map") <- prep$score_map
  attr(out, "diagnostics") <- diagnostics
  attr(out, "settings") <- list(target = "Conditional plug-in source ratings",
    quantity = if (residual) if (type == "pearson") "Pearson residual" else "raw residual" else
      if (type == "prob") "category probabilities" else "expected rating",
    person_basis = if (cfg$method == "MML") "saved EAP" else "primary JML estimate, including boundary limits",
    score_scale = "modeled categories", calibration_uncertainty = FALSE,
    limitation = "Same-data descriptive values; no convergence, goodness-of-fit, forecast or coverage guarantee.")
  attr(out, "source_readiness") <- mfrmr_get_readiness_record(object)$fit
  out
}

#' @rdname mfrm_fit_responses
#' @export
fitted.mfrm_fit <- function(object, type = c("response", "prob"), ...) {
  rlang::check_dots_empty()
  mfrm_fit_response_values(object, match.arg(type))
}

#' @rdname mfrm_fit_responses
#' @export
residuals.mfrm_fit <- function(object, type = c("response", "pearson"), ...) {
  rlang::check_dots_empty()
  mfrm_fit_response_values(object, match.arg(type), residual = TRUE)
}

#' @rdname mfrm_fit_responses
#' @export
predict.mfrm_fit <- function(object, newdata = NULL, type = c("response", "prob"), ...) {
  rlang::check_dots_empty()
  if (!is.null(newdata)) stop("predict() currently describes retained source ratings only. For supported native MML/JML fits, use predict_mfrm_units() or score_mfrm_calibration() to estimate abilities from supplied responses; those have a different target.", call. = FALSE)
  fitted(object, type = match.arg(type))
}

#' @rdname mfrm_fit_responses
#' @name fitted.mfrm_cml_fit
#' @aliases fitted.mfrm_cml_fit
#' @usage \method{fitted}{mfrm_cml_fit}(object, type = c("response", "prob"), ...)
#' @rawNamespace S3method(fitted,mfrm_cml_fit,fitted.mfrm_fit)
NULL
#' @rdname mfrm_fit_responses
#' @name fitted.mfrm_comparison_fit
#' @aliases fitted.mfrm_comparison_fit
#' @usage \method{fitted}{mfrm_comparison_fit}(object, type = c("response", "prob"), ...)
#' @rawNamespace S3method(fitted,mfrm_comparison_fit,fitted.mfrm_fit)
NULL
#' @rdname mfrm_fit_responses
#' @name residuals.mfrm_cml_fit
#' @aliases residuals.mfrm_cml_fit
#' @usage \method{residuals}{mfrm_cml_fit}(object, type = c("response", "pearson"), ...)
#' @rawNamespace S3method(residuals,mfrm_cml_fit,residuals.mfrm_fit)
NULL
#' @rdname mfrm_fit_responses
#' @name residuals.mfrm_comparison_fit
#' @aliases residuals.mfrm_comparison_fit
#' @usage \method{residuals}{mfrm_comparison_fit}(object, type = c("response", "pearson"), ...)
#' @rawNamespace S3method(residuals,mfrm_comparison_fit,residuals.mfrm_fit)
NULL
#' @rdname mfrm_fit_responses
#' @name predict.mfrm_cml_fit
#' @aliases predict.mfrm_cml_fit
#' @usage \method{predict}{mfrm_cml_fit}(object, newdata = NULL, type = c("response", "prob"), ...)
#' @rawNamespace S3method(predict,mfrm_cml_fit,predict.mfrm_fit)
NULL
#' @rdname mfrm_fit_responses
#' @name predict.mfrm_comparison_fit
#' @aliases predict.mfrm_comparison_fit
#' @usage \method{predict}{mfrm_comparison_fit}(object, newdata = NULL, type = c("response", "prob"), ...)
#' @rawNamespace S3method(predict,mfrm_comparison_fit,predict.mfrm_fit)
NULL

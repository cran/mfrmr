#' Pointwise intervals for a conditional-normal population model
#'
#' Report coefficients, prespecified linear contrasts and residual variance from
#' an estimated RSM/PCM MML population model. Coefficient intervals use a normal
#' approximation; residual
#' variance intervals exponentiate a normal interval for log residual variance.
#' These experimental intervals have unqualified repeated-sample coverage.
#' Method qualification is incomplete even when numerical checks pass.
#' Saved cautions are also retained as an `interpretation_notes` table in
#' [mfrm_results()], reports and table exports.
#'
#' @param object A native RSM/PCM MML fit with an active population model, or
#'   saved native-location [mfrm_facet_intervals()] output from that fit. The
#'   latter reuses its numerically checked full joint covariance without differentiation.
#'   For `summary()` and `vcov()`, supply the population-interval result.
#' @param level Pointwise confidence level, strictly between zero and one.
#' @param contrasts Optional numeric matrix of prespecified linear combinations
#'   of population coefficients. Give each row a distinct nonempty target name
#'   and name every column with a fitted coefficient name, including zero-weight
#'   columns, as in `names(fit$population$coefficients)`. Columns are matched by
#'   name. Each row must be finite and nonzero; dependent rows are allowed for
#'   pointwise intervals. `NULL` retains the coefficient/variance-only output.
#' @inheritParams mfrm_facet_intervals
#' @details
#' For the bounded fixed-SD anchor-origin route, `anchor_covariance` with
#' `independent_anchors = TRUE` propagates an independently estimated anchor
#' covariance into the population coefficients and contrasts. Alternatively,
#' pass a saved facet-interval result already containing source propagation and
#' omit both arguments; its sensitivity and conditional information are reused
#' once. See [mfrm_facet_intervals()] for the formula, required independence,
#' supported scope and complete bundled-data example. Point estimates are
#' unchanged; conditional SE/bounds remain in separate `Conditional*` columns.
#' `vcov()` then returns source-propagated coefficient covariance. The fixed
#' population variance still has no SE or interval. Drift, source-estimate bias,
#' dependent samples and uncertainty in the supplied SD are not corrected.
#'
#' The fitted model is \eqn{\theta_i\mid X_i\sim N(X_i\beta,\sigma_e^2)}.
#' Residual variance is conditional on the fitted predictors; it is not the
#' marginal variance of ability. Predictor transformations and coefficient
#' coding are those retained by the fit. The covariance is taken from the
#' inverse of the full information, including all structural nuisance terms.
#' A contrast matrix \eqn{L} gives estimates \eqn{L\widehat\beta} and covariance
#' \eqn{L V_\beta L^{\mathsf T}}, using the coefficient block of that inverse.
#' Contrast rows are added to the interval table. These describe conditional
#' population means or their linear combinations, not differences between two
#' individual abilities or causal effects. For example, in `~ X`, weights
#' `(0, 2)` compare the means at `X = 1` and `X = -1`. Interactions, factor
#' reference levels and transformations must use the fitted coding. The helper
#' accepts weights in coefficient coordinates; it does not infer raw covariate
#' profiles or rebuild a transformation from a new sample.
#'
#' Use direct MML with fixed or adaptive integration, unit observation weights,
#' full-rank population coding, centered additive facets and centered steps,
#' without constraint groups or
#' fitted interactions. The fixed-SD anchor-origin extension described below
#' also admits supplied anchors. Source consistency, category
#' support, local solution, full joint information and finer-integration checks
#' must pass as in [mfrm_facet_intervals()]. The procedure can return eligible
#' intervals while the source fit retains `InferenceReady = FALSE`; inspect
#' the interval's `CIEligible`, `Status`, checks and cautions. This does not
#' change overall fit readiness or qualify other inference procedures.
#' A recorded identification or boundary failure prevents intervals. Failed numerical checks retain point
#' estimates and reasons with unavailable uncertainty. This is an interior,
#' model-based approximation: it does not handle a zero-variance boundary,
#' misspecified populations, dependent Persons or informative assignment.
#' Passing numerical checks does not prove global optimality or coverage.
#' When residual variance is estimated, a zero-residual-variance candidate at the retained nuisance estimates is also
#' compared with the fitted likelihood at both integration orders. A higher
#' boundary likelihood beyond the numerical allowance prevents intervals.
#' This inexpensive comparison does not optimize the boundary or rule out other
#' better boundary candidates. Requests using saved facet inference also perform
#' this likelihood comparison, without refitting or differentiating again.
#' A fit with explicit `population_sd` conditions on that positive fixed SD;
#' zero variance is not a feasible candidate. Its residual-variance row retains
#' the supplied value with `Status = "fixed"`, missing SE/bounds and
#' `CIEligible = FALSE`. Coefficients and their contrasts still have model-based
#' intervals when numerical checks pass. With `population_formula = ~1`, direct
#' anchors on the selected `noncenter_facet` can identify the population mean
#' and free facet locations. By default, these intervals condition on the
#' supplied anchors and SD. The optional `anchor_covariance` extension adds
#' independent source-anchor uncertainty while still conditioning on the SD.
#'
#' \code{SE} is on the reported estimate scale (delta-method SE for variance).
#' \code{WorkingEstimate} and \code{WorkingSE} use coefficient/log-variance
#' coordinates. \code{WorkingLower} and \code{WorkingUpper} retain the interval
#' on that scale. Exponentiation can give displayed variance bounds of zero or
#' infinity outside floating-point range; \code{NaturalScaleStatus} identifies
#' that representation limit while the finite log-scale interval is retained.
#' \code{vcov()} returns their covariance, with explicit coordinate
#' names; its final coordinate is log residual variance, not natural variance,
#' when variance is estimated. A fixed variance has no free covariance coordinate.
#' It remains the coefficient/log-variance covariance when contrasts are supplied;
#' `contrast_covariance` separately contains the covariance of the requested rows.
#' No multiplicity adjustment or automatic significance decision is made.
#'
#' Save the result directly, or attach it through
#' \code{mfrm_results(fit, intervals = list(population = ci), compute = "never")}.
#' Tables, reports and exports reuse saved intervals and numerical checks.
#' Supplied contrast weights and their covariance are also exported. These
#' tables prefix coordinate columns with `beta:`/`contrast:`
#' to distinguish them from the `Target` identifier column.
#' The saved result retains the unchanged fit, including its population formula, trained
#' terms, factor levels and coding; save the raw reference profiles separately
#' if you used them to construct the weights.
#' The ordinary population summary remains a point-estimate table.
#'
#' @section Reading the result:
#' Start with `summary(ci)`, then inspect `ci$checks`, `ci$numerical_checks`
#' and `ci$cautions`. Keep these with the interval table when saving or sharing.
#' [plot(ci)][plot.mfrm_population_intervals] displays saved coefficients;
#' `plot(ci, type = "contrasts")` displays the supplied linear contrasts.
#' Unavailable intervals remain visible, and [as_ggplot()] reuses the same data.
#' - `Status = "available"` and `CIEligible = TRUE` mean that the numerical
#'   checks permit the experimental interval. They do not certify nominal
#'   coverage. `InferenceReview` records this distinction in each row.
#' - `Status = "fixed"` retains a supplied residual variance without an
#'   inferential interval; missing uncertainty is not an estimated zero SE.
#' - `Status = "numerical_review_failed"` retains the point estimate and
#'   failed-check reason, with unavailable SEs and bounds. `nonpositive_variance`
#'   or `nonfinite_bounds` instead identifies an unusable target-specific
#'   uncertainty calculation. Missing bounds are not infinite endpoints.
#' - `NaturalScaleStatus = "outside_floating_point_range"` means the finite
#'   log-variance interval cannot be fully represented after exponentiation.
#'   Inspect `WorkingLower` and `WorkingUpper`; a displayed zero is not an
#'   interval qualified for a zero-variance boundary.
#'
#' Coverage evidence is specific to the target, design and complete estimation
#' procedure. In the retained RSM checks, mean-at-\eqn{X=1} and log-variance
#' intervals at 20 Persons remain inconclusive under the prespecified practical
#' criteria. The independent 240-Person slope/mean-difference confirmation
#' supports those targets in its tested design; it does not establish a minimum
#' sample size or qualify arbitrary regressions, PCM designs or all contrasts.
#' All outputs from this helper therefore retain their experimental status.
#'
#' These intervals do not add calibration uncertainty to individual scoring.
#' [predict_mfrm_units()] has separate source and scoring checks; its conditional
#' posterior intervals hold the fitted calibration and population parameters
#' fixed. Availability here does not authorize an otherwise refused scoring
#' result. See [mfrmr_interval_guide()] for the different uncertainty targets.
#' @return An \code{mfrm_population_intervals} object with the interval table,
#'   working-coordinate covariance, full parameter covariance, numerical checks,
#'   source identity and unchanged source fit. With `contrasts`, also contains
#'   the column-aligned `contrasts` matrix and `contrast_covariance`; targets in
#'   the table have type `coefficient_contrast` and prefix `contrast:`.
#' @seealso [mfrm_facet_intervals()], [predict_mfrm_units()], [mfrm_results()],
#'   [plot.mfrm_population_intervals()], [apa_table()], [export_mfrm_results()]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("mfrm_population_intervals", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_core")
#' persons <- data.frame(Person = unique(ratings$Person))
#' # Fictional predictor for this example; use your observed covariates.
#' persons$X <- seq(-1, 1, length.out = nrow(persons))
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'   model = "RSM", method = "MML", population_formula = ~ X,
#'   person_data = persons, quad_points = 61)
#' L <- matrix(c(0, 2), nrow = 1,
#'   dimnames = list("X = 1 minus X = -1", c("(Intercept)", "X")))
#' ci <- mfrm_population_intervals(fit, contrasts = L)
#' subset(summary(ci), Type == "coefficient_contrast")
#' plot(ci)
#' p <- plot(ci, type = "contrasts")
#' plot_data(p, "contrast_weights")
#' plot_data(p, "notes")
#' results <- mfrm_results(fit, intervals = list(population = ci),
#'   include = "fit", compute = "never")
#' apa_table(results$population_intervals$population)
#' }
#' @export
#' @order 1
mfrm_population_intervals <- function(object, level = .95, contrasts = NULL,
                                      anchor_covariance = NULL, independent_anchors = FALSE) {
  if (inherits(object, "mfrm_facet_intervals") && !is.null(object$anchor_uncertainty)) {
    if (!is.null(anchor_covariance) || !identical(independent_anchors, FALSE))
      stop("Saved anchor-source intervals already contain the covariance and independence declaration; omit both arguments when reusing them.", call. = FALSE)
    mfrmr_anchor_uncertainty_validate_saved(object)
    request <- object$anchor_uncertainty
    object$parameter_covariance <- request$conditional_parameter_covariance
  } else {
    fit <- if (inherits(object, "mfrm_facet_intervals")) object$fit else object
    request <- mfrmr_anchor_uncertainty_request(fit, anchor_covariance, independent_anchors)
  }
  mfrm_population_intervals_impl(object, level, contrasts, anchor_request = request)
}

mfrm_population_intervals_impl <- function(object, level, contrasts, covariance = NULL, anchor_request = NULL) {
  if (!is.numeric(level) || is.complex(level) || is.object(level) || !is.null(dim(level)) ||
      length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("Supply 0 < level < 1.", call. = FALSE)
  level <- unname(level)
  reuse <- inherits(object, "mfrm_facet_intervals")
  fit <- if (reuse) object$fit else object
  stop_if_jml_adjustment(fit, "mfrm_population_intervals()")
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !fit$config$model %in% c("RSM", "PCM") || !isTRUE(fit$config$population_spec$active))
    stop("Supply a native RSM/PCM fit with an estimated normal population.", call. = FALSE)
  params <- mfrm_native_location_source(fit)
  identity <- mfrm_native_location_identity(fit)
  columns <- fit$config$population_spec$design_columns
  if (!is.null(contrasts)) {
    if (!is.matrix(contrasts) || !is.numeric(contrasts) || is.complex(contrasts) ||
        !nrow(contrasts) || ncol(contrasts) != length(columns) || any(!is.finite(contrasts)) ||
        is.null(colnames(contrasts)) || anyDuplicated(colnames(contrasts)) ||
        !setequal(colnames(contrasts), columns) || is.null(rownames(contrasts)) ||
        anyNA(rownames(contrasts)) || any(!nzchar(trimws(rownames(contrasts)))) ||
        anyDuplicated(rownames(contrasts)) || any(rowSums(abs(contrasts)) == 0))
      stop("`contrasts` must have distinct target row names, all fitted coefficient column names, and finite nonzero rows.", call. = FALSE)
    contrasts <- contrasts[, columns, drop = FALSE]
  }
  if (reuse) {
    if (!identical(object$settings$procedure, "mml_native_location_model_v1") ||
        !identical(identity, object$source) || !is.data.frame(object$checks) || !nrow(object$checks))
      stop("Saved native-location inference must match its source fit and qualification metadata.", call. = FALSE)
    eligible <- all(object$checks$Passed %in% TRUE) && !is.null(object$parameter_covariance)
    required <- c("Identification and boundary", "Category support", "Joint information", "Local solution",
      "Full gradient", "Stationarity and inverse", "Comparison information", "Quadrature comparison")
    if (eligible && (!identical(object$checks$Check, required) ||
        !is.data.frame(object$numerical_checks) || nrow(object$numerical_checks) != 1L ||
        !all(c("CurvatureScaledGradient", "InverseResidual", "ComparisonInverseResidual",
          "QuadratureScoreShift", "QuadratureCovarianceChange") %in% names(object$numerical_checks))))
      stop("Saved native-location inference needs its complete numerical qualification record.", call. = FALSE)
    inference <- list(covariance = if (eligible) object$parameter_covariance else NULL,
      checks = object$checks, numerical_checks = object$numerical_checks,
      comparison_information_review = object$comparison_information_review,
      check = list(eligible = eligible,
        review = paste(object$checks$Detail[!object$checks$Passed], collapse = " "),
        caution = object$cautions[!startsWith(object$cautions, "Experimental native location intervals:") &
          !startsWith(object$cautions, "Native location intervals:") &
          !startsWith(object$cautions, "Population-SD-standardized contrasts use ")]))
    inference <- mfrm_native_location_boundary_candidate(fit, inference)
    information_review <- object$information_review
  } else {
    inference <- if (is.null(covariance)) mfrm_native_location_inference(fit) else
      mfrm_native_location_inference(fit, covariance = covariance)
    information_review <- inference$information$solution_information$inverse_review
  }
  sizes <- build_param_sizes(fit$config); slices <- build_param_slices(sizes)
  fixed_variance <- !is.null(fit$config$population_spec$fixed_sd)
  at <- c(slices[["beta"]], slices[["log_sigma2"]])
  working <- c(unname(params$population$coefficients), log(params$population$sigma2))
  targets <- c(paste0("beta:", columns), "log_residual_variance")
  V <- inference$covariance
  if (!is.null(V) && (!is.matrix(V) || !identical(dim(V), rep(length(fit$opt$par), 2L)) ||
      any(!is.finite(V)) || !isTRUE(all.equal(V, t(V))) ||
      is.null(tryCatch(chol(V), error = function(e) NULL))))
    stop("Saved joint covariance must be finite, positive definite and match all free coordinates.", call. = FALSE)
  covariance <- if (is.null(V)) matrix(NA_real_, length(at), length(at)) else V[at, at, drop = FALSE]
  covariance_targets <- if (fixed_variance) head(targets, -1L) else targets
  dimnames(covariance) <- list(covariance_targets, covariance_targets)
  se <- covariance_diag_se(covariance)
  if (fixed_variance) se <- c(se, NA_real_)
  terms <- columns; types <- rep("coefficient", length(columns))
  contrast_covariance <- NULL
  if (!is.null(contrasts)) {
    beta_at <- seq_along(columns)
    contrast_covariance <- contrasts %*% covariance[beta_at, beta_at, drop = FALSE] %*% t(contrasts)
    working <- c(working[beta_at], drop(contrasts %*% working[beta_at]), tail(working, 1L))
    se <- c(se[beta_at], covariance_diag_se(contrast_covariance), tail(se, 1L))
    targets <- c(targets[beta_at], paste0("contrast:", rownames(contrasts)), tail(targets, 1L))
    terms <- c(terms, rownames(contrasts))
    types <- c(types, rep("coefficient_contrast", nrow(contrasts)))
  }
  status <- rep(if (isTRUE(inference$check$eligible)) "available" else "numerical_review_failed", length(working))
  status[status == "available" & (!is.finite(se) | se <= 0)] <- "nonpositive_variance"
  if (fixed_variance) status[length(status)] <- "fixed"
  se[status != "available"] <- NA_real_
  margin <- stats::qnorm((1-level)/2, lower.tail=FALSE)*se
  working_lower <- working-margin; working_upper <- working+margin
  lower <- working_lower; upper <- working_upper
  estimate <- working; natural_se <- se; last <- length(working)
  estimate[last] <- params$population$sigma2
  lower[last] <- exp(lower[last]); upper[last] <- exp(upper[last])
  natural_se[last] <- estimate[last]*se[last]
  valid <- is.finite(working_lower) & is.finite(working_upper) & is.finite(se)
  natural_status <- rep("finite", length(working))
  if (valid[last] && (lower[last] == 0 || !is.finite(upper[last]) || !is.finite(natural_se[last])))
    natural_status[last] <- "outside_floating_point_range"
  status[status == "available" & !valid] <- "nonfinite_bounds"
  unavailable <- status != "available"
  se[unavailable] <- natural_se[unavailable] <- lower[unavailable] <- upper[unavailable] <- NA_real_
  working_lower[unavailable] <- working_upper[unavailable] <- NA_real_
  natural_status[unavailable] <- "unavailable"
  if (fixed_variance) natural_status[last] <- "finite"
  review <- if (isTRUE(inference$check$eligible))
    "Experimental population normal/log-variance intervals; numerical checks passed, sampling coverage unqualified." else
    inference$check$review
  table <- data.frame(Target = targets, Type = c(types, "residual_variance"),
    Term = c(terms, "Residual variance"), Estimate = estimate, SE = natural_se,
    Lower = lower, Upper = upper, WorkingEstimate = working, WorkingSE = se,
    WorkingLower = working_lower, WorkingUpper = working_upper,
    WorkingScale = c(rep("coefficient", length(terms)), "log_residual_variance"),
    NaturalScaleStatus = natural_status,
    Status = status, CIEligible = !unavailable, InferenceReview = review, row.names = NULL)
  table$InferenceReview[status %in% c("nonpositive_variance", "nonfinite_bounds")] <-
    "Finite positive uncertainty and finite bounds on the coefficient/log-variance scale are required."
  table$InferenceReview[status == "fixed"] <-
    "Residual variance is fixed by population_sd; no estimated variance uncertainty or inferential interval is supplied."
  caution <- "Experimental pointwise population intervals: sampling coverage is unqualified and method qualification is incomplete. Residual variance is conditional on the fitted predictors; these intervals do not include population misspecification or dependence."
  out <- list(table = table, covariance = covariance, parameter_covariance = V,
    checks = inference$checks, numerical_checks = inference$numerical_checks,
    information_review = information_review,
    comparison_information_review = inference$comparison_information_review,
    settings = list(procedure = "mml_normal_population_model_v1", method = "model", level = level,
      model = fit$config$model, population = "conditional_normal", quad_points = fit$config$estimation_control$quad_points,
      integration = fit$config$estimation_control$mml_integration %||% "fixed",
      inference_reused = reuse, adjustment = "pointwise"),
    source = identity, fit = fit, cautions = unique(c(caution, inference$check$caution)))
  if (fixed_variance) {
    out$settings$population_sd <- fit$config$population_spec$fixed_sd
    out$settings$population_sd_estimated <- FALSE
    out$cautions <- c(out$cautions,
      "Population coefficient intervals condition on the supplied population SD; its uncertainty is not propagated.")
  }
  anchor_note <- mfrm_gpcm_anchor_note(fit)
  if (!is.null(anchor_note)) out$cautions <- c(out$cautions, anchor_note)
  if (!is.null(contrasts)) {
    out$contrasts <- contrasts
    out$contrast_covariance <- contrast_covariance
    out$cautions <- c(out$cautions, "Contrasts use the fitted population coding and describe conditional population means, not individual ability differences or causal effects.")
  }
  class(out) <- "mfrm_population_intervals"
  out <- mfrmr_anchor_uncertainty_intervals(out, anchor_request)
  warning(paste(out$cautions, collapse = " "), call. = FALSE)
  out
}

#' @rdname mfrm_population_intervals
#' @param x For `print.mfrm_population_intervals()`: An interval result.
#'
#' For `plot()`: A saved result from [mfrm_population_intervals()].
#' @param ... For `print.mfrm_population_intervals()`: Unused.
#'
#' For `plot()`: Unused; unsupported arguments are rejected.
#' @export
print.mfrm_population_intervals <- function(x, ...) {
  cat("Conditional-normal population intervals\n")
  print(x$table, row.names=FALSE)
  for (caution in x$cautions) print_wrapped_line(paste0("Caution: ", caution))
  invisible(x)
}

#' @rdname mfrm_population_intervals
#' @name summary.mfrm_population_intervals
#' @aliases summary.mfrm_population_intervals
#' @usage \method{summary}{mfrm_population_intervals}(object, ...)
#' @rawNamespace S3method(summary,mfrm_population_intervals,summary.mfrm_pooled)
NULL

#' @rdname mfrm_population_intervals
#' @importFrom stats vcov
#' @name vcov.mfrm_population_intervals
#' @aliases vcov.mfrm_population_intervals
#' @usage \method{vcov}{mfrm_population_intervals}(object, ...)
#' @rawNamespace S3method(vcov,mfrm_population_intervals,vcov.mfrm_facet_intervals)
NULL

mfrm_population_interval_tables <- function(x) {
  tab <- x$table; tab$ConfidenceLevel <- x$settings$level
  tables <- list(intervals=tab, settings=as.data.frame(x$settings),
    interpretation_notes=data.frame(Note=as.character(x$cautions)),
    working_covariance=data.frame(Target=rownames(x$covariance), x$covariance,
      check.names=FALSE, row.names=NULL))
  if (!is.null(x$contrasts)) {
    tables$contrast_weights <- data.frame(Target = paste0("contrast:", rownames(x$contrasts)),
      x$contrasts, check.names = FALSE, row.names = NULL)
    tables$contrast_covariance <- data.frame(Target = paste0("contrast:", rownames(x$contrast_covariance)),
      x$contrast_covariance, check.names = FALSE, row.names = NULL)
    names(tables$contrast_weights)[-1L] <- paste0("beta:", colnames(x$contrasts))
    names(tables$contrast_covariance)[-1L] <- paste0("contrast:", colnames(x$contrast_covariance))
    coding <- population_coding_summary_table(x$fit$population)
    if (nrow(coding)) tables$population_coding <- coding
  }
  tables <- c(tables, mfrmr_anchor_uncertainty_tables(x))
  for (name in c("checks", "numerical_checks", "information_review", "comparison_information_review"))
    if (!is.null(x[[name]])) tables[[name]] <- x[[name]]
  tables
}

mfrm_population_results_inputs <- function(fit, intervals) {
  if (inherits(intervals, "mfrm_population_intervals")) intervals <- list(inference=intervals)
  selected <- if (is.list(intervals)) vapply(intervals, inherits, logical(1), what="mfrm_population_intervals") else logical()
  if (!any(selected)) return(list(other=intervals, population=NULL))
  if (is.null(names(intervals)) || anyNA(names(intervals)) || anyDuplicated(tolower(names(intervals))) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", names(intervals))))
    stop("Name saved intervals with distinct letter/number/underscore names.", call. = FALSE)
  mfrm_native_location_source(fit)
  identity <- mfrm_native_location_identity(fit)
  population <- intervals[selected]; names(population) <- tolower(names(population))
  for (x in population) {
    if (!identical(x$settings$procedure, "mml_normal_population_model_v1") || !identical(x$source, identity))
      stop("Saved population intervals must match the source fit and qualification metadata.", call. = FALSE)
  }
  list(other=if (all(selected)) NULL else intervals[!selected], population=population)
}

mfrm_population_results_attach <- function(out, inputs) {
  if (is.null(inputs)) return(out)
  out$population_intervals <- inputs
  for (name in names(inputs)) {
    key <- paste0("population_", name)
    tables <- mfrm_population_interval_tables(inputs[[name]])
    names(tables) <- paste(key, names(tables), sep="_")
    out$tables <- c(out$tables, tables); out$components[[key]] <- inputs[[name]]
    out$status <- rbind(out$status, mfrm_results_status_row(key, "review",
      "Saved experimental population intervals; conditional residual variance, numerical checks and unavailable targets retained."))
  }
  out$table_index <- mfrm_results_table_index(out$tables)
  out$notes <- unique(c(out$notes, unlist(lapply(inputs, `[[`, "cautions"), use.names=FALSE)))
  out
}

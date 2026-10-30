#' Read an exported Bayesian measurement summary
#'
#' Bring a saved Stan/brms summary back into mfrmr for figures and reports,
#' without loading either backend or running an estimator.
#' @param path An export directory containing `posterior-summary.rds` and
#'   its matching `model-data.rds`, or the path to that summary file.
#' @return An `mfrm_posterior_summary` containing unrounded parameter and
#'   diagnostic tables, source settings, priors, mappings and constraints.
#' @details Run the exported bundle's `Rscript run.R --summarize` first.
#' The reader checks the recorded model-data checksum, parameter mappings,
#' fixed-value semantics and consistency of the saved sampling checks. Older
#' summaries without a checksum must be regenerated with the current bundled
#' `run.R` and `summarize.R`; no new sampling is needed. Keep the original
#' `model-data.rds` bytes with the fit when moving an export directory.
#'
#' These checks detect mismatched files and inconsistent summaries; they do
#' not independently verify draws, convergence, prior adequacy or model fit.
#' The object contains marginal summaries, not joint posterior draws or a
#' native fitted model. It supports [mfrm_report()], [apa_table()],
#' [build_summary_table_bundle()] and [export_summary_appendix()].
#' Use the backend's saved fit for new credible levels, contrasts, predictions,
#' trace plots or posterior predictive checks. Native scoring and MML/JML
#' diagnostics are not available from this object.
#' Estimated-population one-family Stan summaries retain the population mean
#' and SD as reported parameters. Their schema-2 files require this newer reader;
#' existing fixed-population schema-1 summaries remain supported. Use
#' `plot(x, quantity = "population_sd")` to inspect between-Person variation
#' on the native metric. This SD is not an individual Person standard error.
#'
#' `confint()` retrieves the stored equal-tailed 95% *posterior credible*
#' intervals, retaining fixed/unavailable rows and review flags. It does not
#' calculate frequentist confidence intervals or change the saved level.
#' Fixed anchors and constraint-determined values have no estimated SD or
#' interval. External calibration uncertainty has not been propagated.
#' `summary()` returns the complete validated summary; printing is compact.
#' The bundled synthetic example uses short chains and requires sampling
#' review. Its numbers are workflow illustrations, not assessment findings.
#' @seealso [export_mfrm_model()], [plot.mfrm_posterior_summary()]
#' @examples
#' directory <- system.file("extdata", "posterior-example", package = "mfrmr")
#' posterior <- read_mfrm_posterior(directory)
#' summary(posterior) # Read the decision and next action before the estimates.
#' intervals <- confint(posterior)
#' head(intervals[, c("facet", "level", "value_role", "lower", "upper", "interval_status")])
#' report <- mfrm_report(posterior, style = "apa", study_context = c(
#'   score_direction = "Higher coded scores indicate more favorable performance"))
#' summary(report, view = "reader")
#' @export
#' @order 1
read_mfrm_posterior <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("`path` must identify one saved export summary or directory.", call. = FALSE)
  file <- if (dir.exists(path)) file.path(path, "posterior-summary.rds") else path
  model_file <- file.path(dirname(file), "model-data.rds")
  if (!file.exists(file) || !file.exists(model_file))
    stop("Keep posterior-summary.rds and its matching model-data.rds together.", call. = FALSE)
  x <- readRDS(file)
  if (!is.list(x) || length(x$schema_version) != 1L || !x$schema_version %in% c(1L, 2L) ||
      !identical(x$model_data_md5, unname(tools::md5sum(model_file))))
    stop("The summary has no matching model-data checksum. Regenerate it with the current run.R and summarize.R using --summarize; no new sampling is needed.", call. = FALSE)
  bundle <- readRDS(model_file)
  if (!is.list(bundle) || !identical(bundle$source_method, "MML") ||
      !is.character(bundle$model) || length(bundle$model) != 1L ||
      !bundle$model %in% c("RSM", "PCM", "GPCM"))
    stop("Expected a model-data bundle from export_mfrm_model().", call. = FALSE)
  helpers <- .mfrm_posterior_helpers()
  expected_version <- if (is.null(bundle$population)) 1L else 2L
  if (!identical(x$schema_version, expected_version) ||
      (!is.null(bundle$schema_version) && !identical(bundle$schema_version, expected_version)) ||
      (expected_version == 2L && (!identical(bundle$schema_version, 2L) ||
        !identical(x$backend, "Stan") || !identical(bundle$model, "GPCM"))))
    stop("Incompatible population bundle and summary versions; regenerate the matching summary.", call. = FALSE)
  x$parameter_map <- helpers$mfrm_parameter_map(bundle)
  x$population <- bundle$population
  d <- bundle$stan_data
  x$source <- data.frame(Model = bundle$model, Backend = x$backend,
    Estimation = "External Bayesian sampling", ExportSourceMethod = bundle$source_method,
    PackageVersion = bundle$package_version,
    Ratings = d$N, Persons = d$J, Categories = d$C,
    Population = helpers$mfrm_population_label(x$population), ResponseMetric = "Adjacent-category logit (D = 1)",
    Interval = x$interval, Level = x$level, Chains = x$chains, ModelDataMD5 = x$model_data_md5)
  x$priors <- bundle$priors
  x$constraints <- bundle$constraints
  x$categories <- bundle$categories
  class(x) <- "mfrm_posterior_summary"
  .mfrm_validate_posterior(x)
  x
}

.mfrm_posterior_helpers <- function() {
  env <- new.env(parent = baseenv())
  sys.source(system.file("model-export", "summarize.R", package = "mfrmr"), envir = env)
  env
}

.mfrm_validate_posterior <- function(x) {
  fail <- function() stop("Inconsistent exported posterior summary; regenerate it from the matching saved fit and bundle.", call. = FALSE)
  if (!inherits(x, "mfrm_posterior_summary") ||
      !identical(x$schema_version, if (is.null(x$population)) 1L else 2L) ||
      !is.character(x$backend) || length(x$backend) != 1L || !x$backend %in% c("Stan", "brms") ||
      !identical(x$level, .95) || !identical(x$interval, "equal-tailed posterior credible interval for varying parameters") ||
      length(x$chains) != 1L || !is.numeric(x$chains) || !is.finite(x$chains) || x$chains < 1 || x$chains != as.integer(x$chains)) fail()
  p <- x$parameters; map <- x$parameter_map
  numeric <- c("mean", "median", "sd", "lower", "upper", "rhat", "ess_bulk", "ess_tail", "mcse_lower", "mcse_upper", "fixed_value")
  required <- c("variable", "quantity", "facet", "level", "transition", "value_role", "role", "reference", "direction",
    "sampling_review", "interval_status", "interval_computed", numeric)
  if (!is.data.frame(p) || !nrow(p) || !all(required %in% names(p)) ||
      !is.character(p$variable) || anyNA(p$variable) || any(!nzchar(p$variable)) || anyDuplicated(p$variable) ||
      !all(vapply(p[numeric], function(v) is.numeric(v) && all(is.na(v) | is.finite(v)), logical(1))) ||
      !is.data.frame(map) || !all(c("variable", "value_role", "fixed_value") %in% names(map)) ||
      anyNA(map$variable) || anyDuplicated(map$variable) || !all(map$variable %in% p$variable)) fail()
  if (any(as.matrix(p[c("sd", "ess_bulk", "ess_tail", "mcse_lower", "mcse_upper")]) < 0, na.rm = TRUE) ||
      any(p$rhat <= 0, na.rm = TRUE) ||
      any(as.matrix(p[p$quantity %in% c("slope", "population_sd"), c("mean", "median", "lower", "upper"), drop = FALSE]) <= 0, na.rm = TRUE)) fail()
  at <- match(map$variable, p$variable)
  for (name in names(map)) if (!identical(p[[name]][at], map[[name]])) fail()
  internal <- !p$variable %in% map$variable
  if (anyNA(p$value_role) || !all(p$value_role %in% c("estimated", "fixed_anchor", "fixed_by_constraint", "sampling_coordinate")) ||
      !identical(p$role, ifelse(p$variable %in% map$variable, "reported_parameter", "free_coordinate")) ||
      any(p$value_role[internal] != "sampling_coordinate") || anyNA(p$quantity) ||
      any(p$quantity[internal] != "internal_coordinate") ||
      any(!is.na(as.matrix(p[internal, c("facet", "level", "transition", "direction", "fixed_value"), drop = FALSE])))) fail()
  fixed <- p$value_role %in% c("fixed_anchor", "fixed_by_constraint")
  if (any(!is.finite(p$mean[fixed])) || any(abs(p$mean[fixed] - p$fixed_value[fixed]) > 1e-7 * pmax(1, abs(p$fixed_value[fixed]))) ||
      any(!is.na(as.matrix(p[fixed, setdiff(numeric, c("mean", "median", "fixed_value")), drop = FALSE]))) ||
      any(p$lower > p$upper, na.rm = TRUE) ||
      !identical(p$interval_computed, is.finite(p$lower) & is.finite(p$upper) & !fixed)) fail()
  if (!is.data.frame(x$sampler) || !is.data.frame(x$checks) || !is.data.frame(x$decision)) fail()
  expected <- .mfrm_posterior_helpers()$mfrm_sampling_review(p, x$sampler, x$chains)
  fields <- c("Check", "Value", "Criterion", "Status")
  if (!all(fields %in% names(x$checks)) ||
      !isTRUE(all.equal(x$checks[fields], expected[fields], check.attributes = FALSE, tolerance = 0))) fail()
  status <- if (any(expected$Status != "no_flag")) "sampling_review_required" else "no_sampling_flags_detected"
  if (nrow(x$decision) != 1L || !identical(x$decision$Status, status) ||
      !all(c("Interpretation", "Why", "NextAction") %in% names(x$decision)) ||
      anyNA(p$sampling_review) || any(p$sampling_review != status)) fail()
  interval_status <- ifelse(fixed, "fixed_value_no_interval", ifelse(!p$interval_computed, "unavailable",
    if (status == "sampling_review_required") "computed_sampling_review_required" else "computed_model_review_required"))
  if (!identical(p$interval_status, interval_status) || !is.data.frame(x$source) || nrow(x$source) != 1L ||
      !identical(x$source$ModelDataMD5, x$model_data_md5) || !identical(x$source$Backend, x$backend) ||
      !identical(x$source$Level, x$level) || !identical(x$source$Population,
        .mfrm_posterior_helpers()$mfrm_population_label(x$population))) fail()
  invisible(x)
}

#' @rdname read_mfrm_posterior
#' @param object A saved `mfrm_posterior_summary`.
#' @param ... Unused.
#' @export
print.mfrm_posterior_summary <- function(x, ...) {
  rlang::check_dots_empty()
  .mfrm_validate_posterior(x)
  cat("External Bayesian measurement summary\n")
  cat(x$source$Model, "|", x$backend, "|", x$source$Population, "\n")
  cat(x$decision$Interpretation, "\nWhy: ", x$decision$Why,
    "\nNext: inspect $checks, then model fit, prior sensitivity and target-specific MCSE.\n", sep = "")
  cat(sum(x$parameters$role == "reported_parameter"), "mapped parameters; fixed values have no estimated interval.\n")
  main <- x$parameters[x$parameters$role == "reported_parameter" & x$parameters$quantity != "log_slope", ]
  tab <- utils::head(main[!duplicated(main[c("quantity", "value_role")]), ], 8L)
  label <- paste(tab$facet, tab$level, sep = ": ")
  steps <- !is.na(tab$transition)
  label[steps] <- paste(label[steps], tab$transition[steps], sep = " | ")
  display <- data.frame(Parameter = label, Role = tab$value_role,
    Mean = tab$mean, Lower = tab$lower, Upper = tab$upper,
    Interval = ifelse(!is.na(tab$fixed_value), "Fixed value", ifelse(!tab$interval_computed, "Unavailable",
      ifelse(tab$sampling_review == "sampling_review_required", "Sampling review", "Model review remains"))))
  print(round_numeric_df(display, digits = 3L), row.names = FALSE)
  cat("Use $parameters, confint(x), plot(x), or mfrm_report(x). Sampling checks do not establish model validity.\n")
  invisible(x)
}

#' @rdname read_mfrm_posterior
#' @export
summary.mfrm_posterior_summary <- function(object, ...) {
  rlang::check_dots_empty()
  .mfrm_validate_posterior(object)
  object
}

#' @rdname read_mfrm_posterior
#' @param parm Optional saved parameter names (for example `"locations[1]"`)
#'   or row indices in `$parameters`. `NULL` returns mapped model parameters.
#' @param level The stored level, currently `0.95`. Other levels require
#'   posterior draws and are rejected rather than approximated from an SD.
#' @export
confint.mfrm_posterior_summary <- function(object, parm = NULL, level = .95, ...) {
  rlang::check_dots_empty()
  .mfrm_validate_posterior(object)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || abs(level - object$level) > 1e-12)
    stop("Only the stored 95% posterior credible interval is available. Use the backend's draws for another level.", call. = FALSE)
  p <- object$parameters
  if (is.null(parm)) at <- which(p$role == "reported_parameter") else if (is.character(parm)) {
    at <- match(parm, p$variable)
  } else if (is.numeric(parm) && all(is.finite(parm)) && all(parm == as.integer(parm))) {
    at <- as.integer(parm)
  } else stop("`parm` must select saved parameter names or row indices.", call. = FALSE)
  if (!length(at) || anyNA(at) || any(at < 1L | at > nrow(p)) || anyDuplicated(at))
    stop("`parm` must select distinct existing parameters.", call. = FALSE)
  out <- p[at, , drop = FALSE]
  out$interval_type <- object$interval
  out$credible_level <- object$level
  attr(out, "notes") <- object$notes
  out
}

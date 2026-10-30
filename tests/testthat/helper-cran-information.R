# Only unchanged test fixtures reuse a freshly computed covariance in CRAN runs.
# Production code and NOT_CRAN=true runs always perform the ordinary calculation.
local_cran_fixture_information <- function(fit, key, .env = parent.frame()) {
  if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) return(invisible(NULL))
  original <- get("compute_mml_parameter_covariance", asNamespace("mfrmr"))
  original_review <- get("prediction_local_calibration_review", asNamespace("mfrmr"))
  own_information <- NULL
  source_fit <- fit
  cache_key <- paste0("cran-information:", key)
  testthat::local_mocked_bindings(
    compute_mml_parameter_covariance = function(res, retained_matrices = 0L,
                                                solution_tolerance = NULL) {
      # Options include the workspace budget and backend selection. Changed
      # data, parameters, options or evaluator bindings must be checked afresh.
      relevant_options <- options()
      relevant_options <- relevant_options[grepl("^mfrmr[.]", names(relevant_options)) |
        names(relevant_options) %in% c("warn", "digits", "OutDec", "scipen")]
      settings <- list(options = relevant_options, retained = retained_matrices,
        solution_tolerance = solution_tolerance,
        original = original,
        evaluator = mfrmr:::make_mfrm_direct_evaluator,
        adaptive = mfrmr:::mfrmr_make_adaptive_mml_evaluator)
      previous <- .mfrmr_test_cache[[cache_key]]
      same <- identical(res, source_fit)
      if (same && !is.null(previous) && identical(res, previous$fit) &&
          identical(settings, previous$settings)) return(previous$value)
      warned <- FALSE
      value <- withCallingHandlers(original(res, retained_matrices, solution_tolerance = solution_tolerance),
        warning = function(w) warned <<- TRUE)
      evaluated <- identical(value$status, "ok") ||
        identical(value$solution_information$status, "evaluated_diagnostic_only")
      if (same && !warned && evaluated)
        .mfrmr_test_cache[[cache_key]] <- list(fit = res, settings = settings, value = value)
      value
    },
    prediction_local_calibration_review = function(fit, covariance = NULL) {
      # The downstream scoring/export tests share an unchanged calibration.
      # Still evaluate every new fit, caller-supplied covariance and modified
      # review/evaluator binding. Failed or warned reviews never enter the cache.
      relevant_options <- options()
      relevant_options <- relevant_options[grepl("^mfrmr[.]", names(relevant_options)) |
        names(relevant_options) %in% c("warn", "digits", "OutDec", "scipen")]
      information <- mfrmr:::compute_mml_parameter_covariance
      if (identical(information, own_information)) information <- original
      settings <- list(options = relevant_options, information = information,
        original = original_review,
        check = mfrmr:::mfrm_ic_fit_check,
        contract = mfrmr:::mfrm_extract_fit_ic_contract,
        native = mfrmr:::mfrm_native_location_source,
        indices = mfrmr:::build_indices, sizes = mfrmr:::build_param_sizes,
        parameters = mfrmr:::make_param_cache, rule = mfrmr:::gauss_hermite_normal,
        evaluator = mfrmr:::make_mfrm_direct_evaluator,
        adaptive = mfrmr:::mfrmr_make_adaptive_mml_evaluator)
      review_key <- paste0("cran-local-review:", key)
      previous <- .mfrmr_test_cache[[review_key]]
      same <- is.null(covariance) && identical(fit, source_fit)
      if (same && !is.null(previous) && identical(fit, previous$fit) &&
          identical(settings, previous$settings)) return(previous$value)
      warned <- FALSE
      value <- withCallingHandlers(original_review(fit, covariance),
        warning = function(w) warned <<- TRUE)
      if (same && !warned && isTRUE(value$eligible))
        .mfrmr_test_cache[[review_key]] <- list(fit = fit, settings = settings, value = value)
      value
    }, .package = "mfrmr", .env = .env)
  own_information <- get("compute_mml_parameter_covariance", asNamespace("mfrmr"))
  invisible(NULL)
}

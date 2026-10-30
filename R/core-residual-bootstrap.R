# Same-design fitted-model reference for exploratory residual PCA.
validate_residual_bootstrap_source <- function(fit) {
  if (!inherits(fit, "mfrm_fit") ||
      inherits(fit, c("mfrm_testlet", "mfrm_random_rater"))) {
    stop("The model bootstrap requires the original fit_mfrm() fit, not detached diagnostics or an extended-model fit.", call. = FALSE)
  }
  config <- fit$config
  replay <- config$replay_inputs
  if (!identical(config$method, "MML") || !config$model %in% c("RSM", "PCM")) {
    stop("The residual model bootstrap currently supports RSM/PCM fitted by MML.", call. = FALSE)
  }
  if (is.null(replay) || is.null(config$estimation_control$quad_points)) {
    stop("The source fit lacks the settings needed to reproduce estimation. Refit it first.", call. = FALSE)
  }
  if (isTRUE(config$population_spec$active) ||
      !identical(config$estimation_control$mml_integration, "fixed") ||
      !identical(config$facet_shrinkage, "none") || isTRUE(replay$shrink_person) ||
      length(config$interaction_specs) > 0L ||
      NROW(replay$anchors) > 0L || NROW(replay$group_anchors) > 0L ||
      length(replay$dummy_facets) > 0L || length(replay$positive_facets) > 0L ||
      any(!is.finite(fit$prep$data$Weight) | fit$prep$data$Weight != 1)) {
    stop(paste("The residual model bootstrap requires a fixed standard-normal population,",
      "fixed quadrature, additive facets, unit weights, and no anchors, shrinkage,",
      "dummy facets or positive facets."), call. = FALSE)
  }
  if (!identical(residual_bootstrap_numerical_state(fit), "ready")) {
    stop("Resolve the source fit's numerical convergence before generating a model reference.", call. = FALSE)
  }
  invisible(TRUE)
}

residual_bootstrap_numerical_state <- function(fit) {
  if (any(!is.finite(fit$opt$par)) || !is.finite(fit$opt$value)) return("failed")
  as.character(mfrmr_readiness_numerical_component(fit$opt)$State)
}

residual_model_bootstrap <- function(fit, overall, by_facet, mode, facets,
                                     max_factors, reps, quantile, seed) {
  bundles <- c(if (!is.null(overall)) list(overall), unname(by_facet))
  scope <- c(if (!is.null(overall)) "overall", rep("facet", length(by_facet)))
  facet <- c(if (!is.null(overall)) NA_character_, names(by_facet))
  observed <- lapply(bundles, extract_pca_eigenvalues)
  draws <- lapply(observed, function(x) matrix(NA_real_, reps, length(x)))
  seeds <- with_preserved_rng_seed(seed, sample.int(.Machine$integer.max, reps))
  trials <- vector("list", reps * length(bundles))
  for (i in seq_len(reps)) {
    warnings <- character()
    numerical <- NA_character_
    readiness <- NULL
    refit_ok <- FALSE
    result <- tryCatch(withCallingHandlers(with_preserved_rng_seed(seeds[i], {
      data <- mfrm_gpcm_bootstrap_generate(fit, seeds[i])
      refit <- mfrm_gpcm_bootstrap_refit(fit, data)
      readiness <- mfrmr_get_readiness_record(refit)$fit
      numerical <- residual_bootstrap_numerical_state(refit)
      if (!identical(numerical, "ready")) {
        stop("The refit did not pass the numerical convergence checks.", call. = FALSE)
      }
      refit_ok <- TRUE
      analyze_residual_pca(refit, mode = mode, facets = facets,
        pca_max_factors = max_factors, parallel = FALSE)
    }), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    simulated <- if (!inherits(result, "error")) {
      c(if (!is.null(result$overall)) list(result$overall), unname(result$by_facet))
    } else NULL
    for (j in seq_along(bundles)) {
      error <- if (inherits(result, "error")) conditionMessage(result) else ""
      original <- bundles[[j]]
      available <- FALSE
      if (!nzchar(error)) {
        candidate <- simulated[[j]]
        eigenvalues <- extract_pca_eigenvalues(candidate)
        error <- original$error %||% candidate$error %||% ""
        if (!nzchar(error) &&
            (!identical(dimnames(original$residual_matrix), dimnames(candidate$residual_matrix)) ||
             !identical(is.na(original$residual_matrix), is.na(candidate$residual_matrix)) ||
             length(eigenvalues) != length(observed[[j]]) || !length(eigenvalues) ||
             any(!is.finite(eigenvalues)))) {
          error <- "The refit did not reproduce the residual columns, observation pattern or complete eigenvalues."
        }
        available <- !nzchar(error)
        if (available) draws[[j]][i, ] <- eigenvalues
      }
      trials[[(i - 1L) * length(bundles) + j]] <- data.frame(
        Replicate = i, Seed = seeds[i], Scope = scope[j], Facet = facet[j],
        RefitAvailable = refit_ok, PCAAvailable = available,
        NumericalState = numerical, Error = error,
        CategoryState = readiness$CategoryState %||% NA_character_,
        BoundaryState = readiness$BoundaryState %||% NA_character_,
        ReviewReasons = readiness$ReasonCodes %||% "",
        Warning = paste(unique(c(warnings,
          if (!is.null(simulated)) simulated[[j]]$warning)), collapse = "; "),
        stringsAsFactors = FALSE)
    }
  }
  trials <- dplyr::bind_rows(trials)
  for (j in seq_along(bundles)) {
    successful <- sum(trials$PCAAvailable[seq.int(j, nrow(trials), length(bundles))])
    bundles[[j]]$parallel <- list(
      table = if (successful == reps) residual_parallel_table(observed[[j]],
        draws[[j]], reps, quantile, "model_bootstrap") else data.frame(),
      successful_reps = successful,
      error = if (nzchar(bundles[[j]]$error %||% "")) {
        paste0("The observed residual PCA is unavailable: ", bundles[[j]]$error)
      } else if (successful < reps) paste0(successful, " of ", reps,
        " model-bootstrap replicates produced usable refits and residual correlations. ",
        "The comparison is unavailable because retaining only successful replicates would change the reference distribution.") else NULL,
      warning = NULL, draws = draws[[j]],
      settings = list(method = "model_bootstrap", reps = reps,
        quantile = quantile, seed = seed %||% NA_integer_))
  }
  list(overall = if (!is.null(overall)) bundles[[1L]],
    by_facet = stats::setNames(bundles[which(scope == "facet")], names(by_facet)),
    trials = trials,
    settings = data.frame(Model = fit$config$model, Method = "MML",
      Population = "Fixed standard normal", Integration = "Fixed quadrature",
      QuadraturePoints = fit$config$estimation_control$quad_points,
      EngineRequested = fit$config$estimation_control$mml_engine_requested,
      ObservationPattern = "Retained response rows; assignment held fixed",
      PersonGeneration = "One independent normal draw per Person, shared across their rows",
      ParameterGeneration = "Fitted facet and step estimates held fixed",
      Refitting = "Original identification and estimation settings",
      ReferenceUse = "Exploratory; no calibrated dimensionality decision",
      stringsAsFactors = FALSE))
}

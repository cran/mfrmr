# The ordinary reference deliberately shares the predictive target, not the
# plug-in diagnostics returned by diagnose_mfrm(). No fitting takes place.
mfrm_ordinary_response_input <- function(fit) {
  cfg <- fit$config; prep <- fit$prep
  if (!inherits(fit, "mfrm_fit") || inherits(fit, "mfrm_imported_fit") ||
      !identical(cfg$model, "RSM") || !identical(cfg$method, "MML")) stop(
    "Ordinary posterior predictive residuals require a native RSM MML fit.", call. = FALSE)
  if (!identical(as.character(fit$summary$NumericalState), "ready") ||
      !isTRUE(mfrm_convergence_state(fit)$code_converged) ||
      identical(fit$readiness$fit$FitReadiness, "blocked")) stop(
    "Resolve the ordinary fit's numerical and blocked source checks before response diagnostics.", call. = FALSE)
  columns <- prep$source_columns %||% cfg$source_columns
  facets <- columns$facets; anchors <- extract_anchor_tables(cfg)
  if (is.null(columns$person) || is.null(columns$score) ||
      length(cfg$facet_interactions) || nrow(as.data.frame(anchors$anchors)) ||
      nrow(as.data.frame(anchors$groups)) || length(cfg$dummy_facets) ||
      !identical(cfg$noncenter_facet, "Person") ||
      !identical(cfg$facet_shrinkage %||% "none", "none") ||
      length(cfg$facet_signs[facets]) != length(facets) ||
      anyNA(cfg$facet_signs[facets]) || any(cfg$facet_signs[facets] != -1) ||
      !is.data.frame(prep$data) || !all(is.finite(prep$data$Weight) & prep$data$Weight == 1)) stop(
    "Use unit weights, additive unanchored severity facets, Person as the noncentered facet, and no shrinkage or interactions.", call. = FALSE)
  pop <- fit$population %||% cfg$population_spec
  estimated <- isTRUE(pop$active)
  if (estimated && (!identical(pop$design_columns, "(Intercept)") ||
      !is.matrix(cfg$population_spec$design_matrix) || ncol(cfg$population_spec$design_matrix) != 1L ||
      anyNA(cfg$population_spec$design_matrix) || any(cfg$population_spec$design_matrix != 1))) stop(
    "Use an intercept-only normal population or the known N(0,1) default; latent regressions are not supported by this response-diagnostic route.", call. = FALSE)
  if (isTRUE((pop$response_rows_omitted %||% 0) > 0)) stop("Population-data omissions are not supported by this response-diagnostic route.", call. = FALSE)
  mu <- if (estimated) as.numeric(pop$coefficients) else 0
  variance <- if (estimated) pop$sigma2 else 1
  if (length(mu) != 1L || !is.finite(mu) || length(variance) != 1L ||
      !is.finite(variance) || variance < 0) stop("The fitted normal population is unavailable.", call. = FALSE)
  map <- prep$score_map
  if (!is.data.frame(map) || !all(c("OriginalScore", "InternalScore") %in% names(map))) stop("The saved fit lacks its category mapping.", call. = FALSE)
  map <- map[order(map$InternalScore), , drop = FALSE]; categories <- as.numeric(map$OriginalScore)
  if (length(categories) < 2L || any(!is.finite(categories)) || any(categories != floor(categories)) ||
      any(diff(categories) != 1) || any(diff(map$InternalScore) != 1) ||
      anyNA(prep$data$score_k) || !all(prep$data$score_k == match(prep$data$Score, map$InternalScore) - 1L)) stop(
    "Retain consecutive integer source categories and their adjacent-category mapping.", call. = FALSE)
  steps <- fit$steps$Estimate
  offset <- rep(0, nrow(prep$data))
  for (f in facets) {
    tab <- fit$facets$others[fit$facets$others$Facet == f, , drop = FALSE]
    if (anyDuplicated(tab$Level)) stop("Fixed-facet levels are duplicated.", call. = FALSE)
    at <- match(as.character(prep$data[[f]]), as.character(tab$Level))
    if ("ParameterStatus" %in% names(tab) && anyNA(match(tab$ParameterStatus[at], c("estimable", "fixed", "weak_information")))) stop("A fixed facet has unavailable parameter status.", call. = FALSE)
    offset <- offset + tab$Estimate[at]
  }
  if (length(steps) != length(categories) - 1L || any(!is.finite(steps)) || any(!is.finite(offset))) stop("Finite fixed-facet and step estimates are required.", call. = FALSE)
  data <- as.data.frame(prep$data[c("Person", facets, "Score")])
  names(data) <- c(columns$person, facets, columns$score)
  data[[columns$score]] <- categories[match(prep$data$Score, map$InternalScore)]
  data[setdiff(names(data), columns$score)] <- lapply(data[setdiff(names(data), columns$score)], as.character)
  rownames(data) <- NULL
  n <- prep$row_retention$Rows[match("input_selected_columns", prep$row_retention$Stage)]
  if (length(n) != 1L || !is.finite(n) || n < nrow(data)) stop("Original input-row accounting is unavailable.", call. = FALSE)
  omitted <- prep$omitted_data; positions <- prep$omitted_input_rows
  source <- data[rep(NA_integer_, n), , drop = FALSE]; rownames(source) <- NULL
  if (n > nrow(data)) {
    if (!is.data.frame(omitted) || nrow(omitted) != n - nrow(data) ||
        length(positions) != nrow(omitted) || anyNA(positions) || anyDuplicated(positions) ||
        any(positions != floor(positions) | positions < 1 | positions > n)) stop(
      "The older fit lacks omitted input-row identities; refit from the complete assigned-score roster.", call. = FALSE)
    if (anyNA(omitted[c("Person", facets)]) || any(!is.na(omitted$Score)) ||
        anyNA(omitted$Weight) || any(omitted$Weight != 1)) stop(
      "Only missing assigned scores are supported; omitted identifiers or excluded weights need a different data review.", call. = FALSE)
    missing_data <- omitted[c("Person", facets, "Score")]; names(missing_data) <- names(data)
    missing_data[setdiff(names(data), columns$score)] <- lapply(missing_data[setdiff(names(data), columns$score)], as.character)
    source[positions, ] <- missing_data
  } else positions <- integer()
  source[setdiff(seq_len(n), positions), ] <- data
  list(data = data, assigned_data = source, input_rows = n, omitted_rows = positions,
    columns = columns, score_levels = categories, y = prep$data$score_k,
    person = as.character(data[[columns$person]]), offset = offset, steps = steps,
    mean = mu, sd = sqrt(variance))
}

mfrm_ordinary_response_probabilities <- function(input, rows, order) {
  rule <- gauss_hermite_normal(order); theta <- input$mean + input$sd * rule$nodes
  indices <- which(input$person == input$person[rows[1L]])
  logp <- lapply(indices, function(i) {
    logw <- outer(theta - input$offset[i], 0:length(input$steps)) -
      matrix(c(0, cumsum(input$steps)), length(theta), length(input$steps) + 1L, byrow = TRUE)
    logw <- logw - row_max_fast(logw)
    logw - mfrm_testlet_logsum_rows(logw)
  })
  loglik <- Reduce(`+`, lapply(seq_along(indices), function(j) logp[[j]][, input$y[indices[j]] + 1L]))
  logw <- (loglik - max(loglik)) + log(rule$weights)
  posterior <- exp(logw - max(logw)); posterior <- posterior / sum(posterior)
  result <- t(vapply(rows, function(i) colSums(exp(logp[[match(i, indices)]]) * posterior), numeric(length(input$score_levels))))
  list(probabilities = result, normalization_error = rowSums(result) - 1)
}

mfrm_response_source <- function(fit) {
  if (mfrm_has_jml_adjustment(fit)) return(list(class=class(fit),config=fit$config,
    facets=fit$facets,steps=fit$steps,slopes=fit$slopes,estimator=fit$jml_adjustment$estimator,
    point_available=fit$jml_adjustment$point$available))
  if (mfrm_extended_fit(fit)) return(mfrm_extended_prediction_source(fit))
  if (mfrm_has_product_slopes(fit)) return(list(class = class(fit), config = fit$config,
    parameters = fit$opt$par, objective = fit$opt$value, convergence = fit$opt$convergence,
    facets = fit$facets$others, steps = fit$steps, slopes = fit$slopes,
    population = fit$population, specification = fit$gmfrm, summary = fit$summary))
  list(class = class(fit), config = fit$config[setdiff(names(fit$config), "attached_diagnostics")],
    facets = fit$facets$others, steps = fit$steps, population = fit$population,
    readiness = fit$readiness, summary = fit$summary, convergence = fit$opt$convergence)
}

mfrm_validate_response_diagnostics <- function(fit, diagnostics) {
  input <- if (mfrm_has_jml_adjustment(fit)) mfrm_jml_response_input(fit) else
    if (mfrm_has_product_slopes(fit)) mfrm_gmfrm_response_input(fit, check_solution = FALSE) else
    if (mfrm_extended_fit(fit)) fit$input else mfrm_ordinary_response_input(fit)
  if (!inherits(diagnostics, "mfrm_response_diagnostics") ||
      !identical(diagnostics$source, mfrm_response_source(fit)) ||
      (mfrm_has_jml_adjustment(fit) && !identical(diagnostics$settings$probability_method,"corrected_jml_plugin")) ||
      !identical(diagnostics$source_data, input$assigned_data %||% input$data) ||
      !identical(diagnostics$source_observed, input$data)) stop(
    "Supply mfrm_response_diagnostics() output with matching calibration and the exact source roster.", call. = FALSE)
  invisible(TRUE)
}

mfrm_response_diagnostic_overview <- function(diagnostics) {
  rows <- diagnostics$rows
  available <- sum(rows$Status == "available_conditional")
  missing <- sum(rows$Status == "missing_score")
  zero <- sum(rows$Status == "zero_variance")
  unresolved <- nrow(rows) - available - missing - zero
  omitted <- nrow(diagnostics$source_data) - nrow(rows)
  status <- if (available + zero == 0L) "not_available" else
    if (unresolved + missing + zero + omitted > 0L) "review" else "available"
  detail <- paste0(available, " of ", nrow(rows),
    " selected rows have standardized residuals; ", unresolved, " unresolved, ",
    missing, " missing scores, ", zero, " zero-variance rows. ", omitted,
    " source rows were not selected. Same-data description with calibration fixed; no calibrated fit cutoffs.")
  data.frame(SourceRows = nrow(diagnostics$source_data), Selected = nrow(rows),
    Available = available, Unresolved = unresolved, Missing = missing,
    ZeroVariance = zero, NotIncluded = omitted, Status = status, Detail = detail)
}

# Older results retain the diagnostic rows but lack their availability overview.
# Reconstruct display metadata on a local copy; never refit or reintegrate.
mfrm_results_restore_response_overview <- function(x) {
  if (!is.null(x$tables$response_overview)) return(x)
  diagnostics <- x$response_diagnostics
  if (is.null(diagnostics) && inherits(x$diagnostics, "mfrm_response_diagnostics"))
    diagnostics <- x$diagnostics
  if (is.null(diagnostics)) return(x)
  overview <- mfrm_response_diagnostic_overview(diagnostics)
  x$tables$response_overview <- overview
  row <- mfrm_results_status_row("response_diagnostics", overview$Status, overview$Detail)
  at <- match("response_diagnostics", x$status$Section)
  if (is.na(at)) x$status <- rbind(x$status, row) else x$status[at, names(row)] <- row
  if (!is.null(x$tables$section_status)) x$tables$section_status <- x$status
  x$table_index <- mfrm_results_table_index(x$tables)
  x
}

mfrm_response_diagnostic_tables <- function(diagnostics) {
  p <- diagnostics$probabilities
  settings <- diagnostics$settings
  if (identical(settings$probability_method,"corrected_jml_plugin"))
    settings$probability_method <- "Conditional fitted probabilities at corrected calibration and Person profiles"
  list(response_overview = mfrm_response_diagnostic_overview(diagnostics),
    response_residuals = diagnostics$rows,
    response_measures = diagnostics$measures,
    response_probabilities = data.frame(
      InputRow = rep(diagnostics$rows$InputRow, ncol(p)),
      Score = rep(as.numeric(colnames(p)), each = nrow(p)), Probability = as.vector(p)),
    response_diagnostic_settings = mfrm_extended_settings_table(settings))
}

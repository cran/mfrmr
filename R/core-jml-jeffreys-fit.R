# A separate point-estimator contract. Never inherit mfrm_fit: its SEs,
# ordinary-MLE boundary replacements and FACETS diagnostics do not apply.
mfrm_is_jeffreys_fit <- function(fit) inherits(fit, "mfrm_jeffreys_fit")

mfrm_fit_jeffreys <- function(a, supplied) {
  allowed <- c("data", "person", "facets", "score", "rating_min", "rating_max",
    "weight", "keep_original", "category_policy", "missing_codes", "model", "method",
    "step_facet", "anchors", "noncenter_facet", "dummy_facets", "positive_facets",
    "centering_scope", "maxit", "reltol", "optimizer", "repeated_ratings",
    "fit_scope", "attach_diagnostics", "jml_penalty")
  unsupported <- setdiff(supplied, allowed)
  if (length(unsupported)) stop("Jeffreys JML does not accept: ",
    paste(unsupported, collapse = ", "), ". See fit_mfrm().", call. = FALSE)
  if (a$method != "JML" || !a$model %in% c("RSM", "PCM") || a$fit_scope != "full" ||
      isTRUE(a$attach_diagnostics) || !a$optimizer %in% c("auto", "BFGS"))
    stop("Jeffreys estimation requires additive RSM/PCM JML, fit_scope = 'full', optimizer = 'auto' or 'BFGS', and attach_diagnostics = FALSE.", call. = FALSE)
  if (is.null(a$rating_min) || is.null(a$rating_max) || !isTRUE(a$keep_original))
    stop("Jeffreys estimation requires explicitly declared rubric limits and category_policy = 'preserve'.", call. = FALSE)
  prep <- prepare_mfrm_data(a$data, a$person, a$facets, a$score,
    a$rating_min, a$rating_max, a$weight, a$keep_original, a$missing_codes, a$repeated_ratings,
    allow_single_person = TRUE, allow_constant_scores = TRUE)
  if (any(prep$data$Weight != 1))
    stop("Jeffreys estimation currently requires unit weights.", call. = FALSE)
  if ("duplicate_person_facet_cells" %in% prep$preparation_notes$Condition)
    stop("Resolve repeated rating cells or explicitly declare repeated_ratings = 'independent' when justified.", call. = FALSE)
  families <- resolve_step_and_slope_facets(a$model, a$step_facet, NULL, prep$facet_names, "JML")
  noncenter <- sanitize_noncenter_facet(a$noncenter_facet, prep$facet_names)
  dummy <- sanitize_dummy_facets(a$dummy_facets, prep$facet_names)
  signs <- build_facet_signs(prep$facet_names, a$positive_facets)
  review <- audit_anchor_tables(prep, a$anchors, noncenter_facet = noncenter,
    dummy_facets = dummy, centering_scope = a$centering_scope)
  invalid <- c("anchor_schema_mismatch", "unknown_anchor_facets", "unknown_anchor_levels",
    "invalid_anchor_values", "duplicate_anchors")
  if (any(vapply(review$issues[invalid], nrow, integer(1)) > 0L) ||
      (!is.null(a$anchors) && nrow(normalize_anchor_df(a$anchors)) != nrow(a$anchors)))
    stop("Jeffreys estimation requires unique, valid, finite direct anchors; no input anchor may be discarded.", call. = FALSE)
  anchors <- review$anchors
  built <- build_estimation_config(prep, a$model, "JML", families$step_facet, NULL,
    weight_col = a$weight, facet_signs = signs$signs, positive_facets = signs$positive_facets,
    noncenter_facet = noncenter, dummy_facets = dummy, anchor_df = anchors,
    group_anchor_df = NULL, centering_scope = a$centering_scope)
  cfg <- built$config
  cfg$jml_penalty <- "jeffreys"
  cfg$estimator <- "jeffreys_canonical_joint"
  cfg$facet_shrinkage <- "none"
  cfg$step_facet_source <- families$step_facet_source
  idx <- build_indices(prep, cfg$step_facet)
  audit <- audit_mfrm_estimability(prep, idx, cfg, built$sizes)
  if (!isTRUE(audit$complete) || length(audit$zero_coordinates) || nrow(audit$null_blocks) ||
      !identical(audit$design$state, "identified") || isTRUE(audit$design$tolerance_sensitive))
    stop("Jeffreys estimation requires an identified, numerically rank-stable additive design; a penalty cannot repair unidentified coordinates.", call. = FALSE)
  cfg$estimability_audit <- audit
  # This audit concerns ordinary likelihood recession, not nonidentification
  # of the penalized estimate. Retain evidence without invoking ordinary guards.
  cfg$category_support_audit <- audit_mfrm_category_support(prep, cfg, built$sizes)
  replay <- a[intersect(allowed, names(a))]
  replay[c("data", "category_policy")] <- NULL
  replay$anchors <- anchors
  replay$step_facet <- cfg$step_facet
  replay$package_version <- as.character(utils::packageVersion("mfrmr"))
  cfg$replay_inputs <- replay
  fit <- structure(list(prep = prep, config = cfg,
    opt = list(par = numeric(sum(unlist(built$sizes))))), class = "mfrm_jeffreys_fit")
  solution <- mfrm_fair_jeffreys_generator(fit, gradient_tol = 1e-4)
  fit$opt <- list(par = solution$par, value = solution$objective,
    convergence = solution$convergence, counts = solution$counts)
  fit$checks <- solution[setdiff(names(solution), c("par", "probabilities", "parameter_map"))]
  fit$parameter_map <- solution$parameter_map
  fit$probabilities <- solution$probabilities
  pa <- expand_params(solution$par, built$sizes, cfg)
  fit$facets <- list(person = data.frame(Person = cfg$person_levels,
    Estimate = unname(pa$theta), PrimaryEstimate = unname(pa$theta)),
    others = build_other_facet_table(cfg, prep, pa))
  fit$steps <- build_step_table(cfg, prep, pa)
  fit$notes <- c(
    "Jeffreys-penalized canonical joint estimate; this is not ordinary FACETS JML or the corrected-GPCM estimator.",
    "Local penalized stationarity and positive information do not establish a global optimum, unbiased expected scores or interval coverage.",
    "No parameter SE, Wald interval, likelihood-ratio test or ordinary FACETS fit diagnostic is supplied.",
    "Use fair_average_table() with an explicit reference_plan for common-plan scores and experimental refit or analytic delta intervals.")
  fit
}

#' Inspect a Jeffreys-penalized fit
#'
#' These methods inspect the separate result returned by
#' [fit_mfrm()] with `method = "JML", jml_penalty = "jeffreys"`.
#' The result does not inherit ordinary JML uncertainty or FACETS diagnostics.
#'
#' @param object,x A `mfrm_jeffreys_fit` result.
#' @param ... Unused; additional arguments are rejected.
#' @return `summary()` returns Person, facet and step point-estimate tables,
#'   checks, category-support evidence and notes. `coef()` returns the three
#'   point-estimate tables. `fitted()` returns category probabilities in prepared
#'   response-row order, with columns from `rating_min` to `rating_max`.
#'   `print()` displays the estimator and numerical status.
#'   `vcov()`, `confint()` and `logLik()` on this fit stop explicitly: parameter
#'   uncertainty and likelihood-based model comparison are not implemented.
#'   To estimate common-plan expected-score intervals, use
#'   [fair_average_table()] followed by [confint.mfrm_fair_average()].
#' @details See [fit_mfrm()] for the penalty, assumptions and restrictions.
#'   `checks$ordinary_nll` is the ordinary negative log likelihood evaluated
#'   at the penalized point; `checks$objective` also includes the information
#'   penalty. Neither is an ordinary maximized likelihood for AIC or testing.
#'   `category_support` describes observed category usage, including missing
#'   categories; its ordinary-MLE warnings are not penalty-rank diagnostics.
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("summary.mfrm_jeffreys_fit", run.dontrun = TRUE)
#' \dontrun{
#' ratings <- load_mfrmr_data("example_core")
#' # Each row is one Person's rating by a Rater on a Criterion.
#' fit <- fit_mfrm(ratings, "Person", c("Rater", "Criterion"), "Score",
#'   method = "JML", jml_penalty = "jeffreys", rating_min = 1, rating_max = 4,
#'   category_policy = "preserve")
#' summary(fit)$persons # Finite penalized Person points, without SEs
#' summary(fit)$checks$max_gradient # Penalized score-equation diagnostic
#' plan <- unique(ratings[c("Rater", "Criterion")])
#' scores <- fair_average_table(fit, reference_plan = plan)
#' scores$reference_scores # Expected averages under the same plan for all Persons
#' intervals <- confint(scores, parm = c("P001", "P002"), method = "delta")
#' intervals[, c("Target", "Estimate", "SE", "Available", "CIEligible")]
#' # Numerical availability is not a qualification of small-sample coverage.
#' }
#' @export
summary.mfrm_jeffreys_fit <- function(object, ...) {
  rlang::check_dots_empty()
  list(estimator = object$config$estimator, model = object$config$model,
    persons = object$facets$person, facets = object$facets$others, steps = object$steps,
    checks = object$checks, category_support = object$config$category_support_audit,
    notes = object$notes)
}

#' @rdname summary.mfrm_jeffreys_fit
#' @export
print.mfrm_jeffreys_fit <- function(x, ...) {
  rlang::check_dots_empty()
  cat(x$config$model, "Jeffreys-penalized joint estimation (experimental)\n")
  cat(x$config$n_person, "Persons; maximum penalized gradient:", signif(x$checks$max_gradient, 4), "\n")
  cat("Parameter uncertainty and score-bias/coverage qualification are unavailable.\n")
  invisible(x)
}

#' @rdname summary.mfrm_jeffreys_fit
#' @export
coef.mfrm_jeffreys_fit <- function(object, ...) {
  rlang::check_dots_empty()
  list(persons = object$facets$person, facets = object$facets$others, steps = object$steps)
}

#' @rdname summary.mfrm_jeffreys_fit
#' @export
fitted.mfrm_jeffreys_fit <- function(object, ...) {
  rlang::check_dots_empty()
  p <- object$probabilities
  colnames(p) <- seq(object$prep$rating_min, object$prep$rating_max)
  p
}

#' @rdname summary.mfrm_jeffreys_fit
#' @export
vcov.mfrm_jeffreys_fit <- function(object, ...) {
  stop("Parameter covariance for the Jeffreys estimator is unavailable; ordinary or penalized inverse curvature is not supplied as sampling covariance.", call. = FALSE)
}

#' @rdname summary.mfrm_jeffreys_fit
#' @export
confint.mfrm_jeffreys_fit <- function(object, ...) {
  stop("Parameter intervals for the Jeffreys estimator are unavailable. For common-plan scores use confint(fair_average_table(object, reference_plan = ...), ...).", call. = FALSE)
}

#' @rdname summary.mfrm_jeffreys_fit
#' @export
logLik.mfrm_jeffreys_fit <- function(object, ...) {
  stop("The penalized objective is not a maximized ordinary likelihood for AIC or likelihood-ratio inference; inspect summary(object)$checks instead.", call. = FALSE)
}

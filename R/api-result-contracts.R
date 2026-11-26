# The public class denotes the statistical result, independently of its implementation.
.mfrm_validate_fit <- function(x) {
  if (inherits(x, "mfrm_drift")) stop("Ordinary time-invariant diagnostics do not apply to dynamic `mfrm_drift` fits. Use summary(), confint(), and plot() on the drift result.", call. = FALSE)
  if (!inherits(x, "mfrm_fit") || !identical(x$schema_version, 1L) ||
      !is.data.frame(x$person) || !is.data.frame(x$facets) || !is.list(x$encoding)) {
    stop("This operation requires a current fit_mfrm() calibration. Refit older saved objects from their original data and model specification; changing their class does not convert their contents.", call. = FALSE)
  }
  invisible(x)
}

.mfrm_validate_gstudy <- function(x) {
  if (!inherits(x, "mfrm_gtheory") || !identical(x$schema_version, 1L) ||
      !x$method_code %in% c("anova", "minque0", "urgenova") ||
      !is.list(x$components) || !is.list(x$design) || !is.data.frame(x$variance_components))
    stop("Recreate this G-Study with mfrm_g_study(); its saved result schema is unavailable or incompatible.", call. = FALSE)
  invisible(x)
}

.mfrm_new_gstudy <- function(x, method, data) {
  x$schema_version <- 1L
  x$method_code <- method
  class(x) <- "mfrm_gtheory"
  if (method == "urgenova") {
    x$components <- x$vcov_components
    if (is.null(x$components)) x$components <- stats::setNames(lapply(x$variance_components$RawEstimate,
      function(v) matrix(v, 1L, 1L, dimnames = list(x$scores, x$scores))), x$variance_components$Component)
    x$design <- list(person = x$person, facets = stats::setNames(x$facets, x$facets), scores = x$scores,
      counts = c(Person = x$n_person, stats::setNames(x$facet_levels, x$facets)),
      nesting = x$nesting, method = x$method, rows = x$n_obs)
    x$data <- data
  } else {
    x$person <- x$design$person
    x$facets <- unname(x$design$facets)
    x$scores <- x$score <- x$design$scores
    x$weights <- stats::setNames(rep(1/length(x$scores), length(x$scores)), x$scores)
    x$nesting <- x$design$nesting
    x$n_obs <- x$design$rows
    x$n_person <- unname(x$design$counts[1L])
    x$facet_levels <- unname(x$design$counts[-1L])
    x$method <- x$design$method
    x$truncate_negative <- FALSE
    x$vcov_components <- x$components
    estimate <- vapply(x$components, function(v) as.numeric(crossprod(x$weights, v %*% x$weights)), 0)
    x$variance_components <- data.frame(Component = names(estimate), RawEstimate = unname(estimate),
      Estimate = unname(estimate), SE = NA_real_, df = NA_real_,
      PercentVariance = if (sum(estimate) > 0) 100*estimate/sum(estimate) else NA_real_,
      ContainsPerson = grepl("Person|Residual", names(estimate)), row.names = NULL)
    x$covariance_basis <- "Score- or composite-specific Gaussian moment covariance; request vcov() explicitly."
    x$d_study <- data.frame()
    x$univariate <- NULL
    x$g_study_summary <- .mfrm_gstudy_metrics(NA_real_, NA_real_, NA_real_)
    if (isTRUE(x$design$balanced) && isTRUE(x$design$complete)) {
      ds <- .mfrm_d_study_moments(x, weights = x$weights)
      x$d_study <- ds$coefficients
      x$univariate <- ds$coefficients[ds$coefficients$Kind == "Score", , drop = FALSE]
      cc <- ds$coefficients[ds$coefficients$Kind == "Composite", , drop = FALSE]
      x$g_study_summary <- .mfrm_gstudy_metrics(cc$UniverseVariance[1L], cc$RelativeErrorVariance[1L], cc$AbsoluteErrorVariance[1L])
    }
    x$elapsed_sec <- NA_real_
  }
  # Optional method-specific numerical ingredients have fixed, explicit slots.
  optional <- c("mean_products", "degrees_of_freedom", "estimation", "data_usage", "component_diagnostics",
    "score_scale", "vc_vcov", "obs_n_per_facet", "non_person_subsets")
  for (nm in setdiff(optional, names(x))) x[nm] <- list(NULL)
  .mfrm_validate_gstudy(x)
  x
}

vcov.mfrm_gtheory <- function(object, ...) {
  .mfrm_validate_gstudy(object)
  if (object$method_code == "urgenova") .mfrm_gstudy_vcov_henderson(object, ...)
  else .mfrm_gstudy_vcov_moments(object, ...)
}
confint.mfrm_gtheory <- function(object, parm = NULL, level = 0.95, ...) {
  .mfrm_validate_gstudy(object)
  if (object$method_code == "urgenova") .mfrm_gstudy_confint_henderson(object, parm, level, ...)
  else .mfrm_gstudy_confint_moments(object, parm, level, ...)
}
mfrm_d_study.mfrm_gtheory <- function(x, design_grid = NULL, ...) {
  .mfrm_validate_gstudy(x)
  if (x$method_code == "urgenova") .mfrm_d_study_henderson(x, design_grid, ...)
  else .mfrm_d_study_moments(x, design_grid, ...)
}
print.mfrm_gtheory <- function(x, ...) {
  .mfrm_validate_gstudy(x)
  if (x$method_code == "urgenova") .mfrm_gstudy_print_henderson(x, ...)
  else .mfrm_gstudy_print_moments(x, ...)
}
summary.mfrm_gtheory <- function(object, ...) {
  .mfrm_validate_gstudy(object)
  if (object$method_code == "urgenova") object <- .mfrm_refresh_gstudy(object)
  list(variance_components = object$variance_components, components = object$components,
    component_diagnostics = object$component_diagnostics, univariate = object$univariate,
    g_study_summary = object$g_study_summary, d_study = object$d_study, method = object$method)
}
plot.mfrm_gtheory <- function(x, type = NULL, palette = c("okabe_ito", "viridis", "monochrome"), base_size = 11, ...) {
 .mfrm_validate_gstudy(x)
 palette <- match.arg(palette)
 if (is.null(type)) type <- if(x$method_code == "urgenova") "dstudy" else "variance"
 type <- match.arg(type, c("dstudy", "variance"))
 if (x$method_code == "urgenova") return(.mfrm_gstudy_plot_henderson(x, type=type, palette=palette, base_size=base_size, ...))
 if (type == "variance") return(.mfrm_gstudy_plot_moments(x, palette=palette, base_size=base_size, ...))
 plot(mfrm_d_study(x, ...), palette=palette, base_size=base_size)
}

# A local quadrature sensitivity check, not a bound on integration error or coverage.
.mfrm_refine_quadrature <- function(fit) {
  q <- fit$control$quad_points
  if (!is.finite(q) || q < 3L || q > (.Machine$integer.max - 1)/2)
    stop("Automatic quadrature requires a valid initial order of at least three.", call. = FALSE)
  control <- fit$control
  control$quad_points <- as.integer(2*q + 1)
  control$warm_init <- as.double(fit$xi)
  t0 <- proc.time()[[3L]]
  raw <- mfrm_cpp_fast_fit_spec(fit$encoding, control)
  fine <- .mfrm_format_result(fit$encoding, raw, fit$method, fit$score_method,
    control$xtreme, control, fit$elapsed_sec + proc.time()[[3L]] - t0)
  metrics <- c(parameter = max(abs(fine$xi - fit$xi)/pmax(1, abs(fine$xi))),
    expected_score = max(abs(fitted(fine) - fitted(fit))),
    loglik_per_rating = abs(fine$loglik - fit$loglik)/fit$n_obs)
  tolerance <- c(parameter = 1e-4, expected_score = 1e-4, loglik_per_rating = 1e-6)
  stable <- isTRUE(fit$converged) && isTRUE(fine$converged) &&
    all(is.finite(metrics)) && all(metrics <= tolerance)
  fine$quadrature_check <- list(orders = c(q, control$quad_points),
    selected_order = control$quad_points, differences = metrics, tolerance = tolerance,
    stable = stable, status = if (stable) "stable_at_checked_orders" else "quadrature_review_required")
  fine$readiness <- .mfrm_readiness(fine)
  fine$summary$InferenceReady <- fine$readiness$ready
  fine$summary$InferenceStatus <- fine$readiness$status
  if (!stable) warning("MML estimates or predictions changed across quadrature orders; increase quad_points and review quadrature_check before inference.", call. = FALSE)
  fine
}

.mfrm_threads <- function(threads) {
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
      threads < 1 || threads > .Machine$integer.max || threads != floor(threads))
    stop("threads must be a positive integer.", call. = FALSE)
  limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  as.integer(if (nzchar(limit) && !limit %in% c("false", "no", "0")) min(threads, 2L) else threads)
}

.mfrm_describe_numeric <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  center <- if (n) mean(x) else NA_real_
  spread <- stats::sd(x)
  z <- if (is.finite(spread) && spread > 0) (x - center) / spread else NA_real_
  data.frame(vars = 1L, n = n, mean = center, sd = spread,
    median = if (n) stats::median(x) else NA_real_, min = if (n) min(x) else NA_real_,
    max = if (n) max(x) else NA_real_, range = if (n) diff(range(x)) else NA_real_,
    skew = mean(z^3), kurtosis = mean(z^4) - 3, se = spread / sqrt(n), row.names = "X1")
}

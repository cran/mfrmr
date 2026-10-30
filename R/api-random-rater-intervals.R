#' Bootstrap prediction intervals for observed random raters
#'
#' Regenerate persons, shared raters and scores on the observed assignment,
#' refit the same RSM, and retain prediction errors for every planned replicate.
#' @param object A numerically ready [fit_mfrm_random_rater()] result with
#'   positive information, positive population SD and finite positive rater
#'   prediction SEs. Estimated zero-variance sources are not supported.
#' @param nsim Number of planned bootstrap datasets; default 499, minimum 2.
#'   Each requires a full refit. Small values are useful for checking a workflow,
#'   not for stable tail quantiles.
#' @param seed Required nonnegative integer simulation seed. The caller's random
#'   number state is restored. Per-replicate seeds and RNG kind are retained.
#' @param level For `mfrm_random_rater_intervals()`: Pointwise nominal coverage; default 0.95.
#'
#' For `plot()`: Nominal pointwise coverage; defaults to the saved level.
#' @param parm Rater IDs for `confint()`; `NULL` returns all fitted raters.
#' @param method For `mfrm_random_rater_intervals()`: `"studentized"` (default) uses generated-minus-estimated rater
#'   effects divided by each refit's calibration-adjusted `PredictionSE`.
#'   `"error"` uses unscaled prediction errors as a comparison. Changing method
#'   or level through `confint()` reuses saved draws without new fitting.
#'
#' For `plot()`: `"studentized"` or `"error"`; see [mfrm_random_rater_intervals()].
#' @param x For `mfrm_random_rater_intervals()`: A bootstrap result.
#'
#' For `plot()`: A result from [mfrm_random_rater_intervals()].
#' @param ... Unused.
#' @details Every generated dataset draws one ability per person and one
#'   severity per rater, shared across their observed rows. Generating fixed
#'   facets, steps and both population SDs are the source estimates. Each
#'   Person ability is drawn from the fitted normal population. The refit
#'   re-estimates calibration and each population SD unless originally fixed.
#'   Earlier saved fits retain their known N(0,1) ability population. Per-trial
#'   records include the refitted ability SD and its zero-boundary status.
#'   The target is a realized effect relative to the population mean, not a
#'   centered effect, new-rater score or fixed rater coefficient. Coverage is
#'   marginal over new persons, rater effects and scores at this assignment;
#'   it is not conditional coverage for each fixed true severity.
#'
#'   Studentized endpoints are the source estimate plus its `PredictionSE`
#'   times the empirical tail quantiles (type 1) of bootstrap errors divided
#'   by refitted `PredictionSE`. Unscaled endpoints add error quantiles directly.
#'   All planned replicates remain in the denominator. Failed fits and
#'   unavailable studentizers are unresolved, including studentizers at
#'   estimated variance boundaries. Lower-tail calculations place unresolved
#'   roots at minus infinity; upper-tail calculations place them at plus
#'   infinity. The resulting limits enclose empirical bootstrap limits for any
#'   completion of those roots. They may be unbounded; failures are never
#'   silently removed or replaced. A numerically ready boundary refit still
#'   supplies an unscaled prediction error.
#'
#'   This is a model-based bootstrap candidate, not a general finite-sample
#'   coverage guarantee. Linear mixed-model bootstrap theory motivates the
#'   construction but does not establish its accuracy for this crossed ordinal
#'   RSM, few raters, variance boundaries or Laplace approximation. Basic error
#'   intervals and studentized intervals need separate studies of their coverage.
#'   The fit's normal-population and assignment assumptions remain essential.
#'   Omitted scores stay omitted; this conditions on analyzed rows and does not
#'   simulate a missingness mechanism or impute assigned scores. Rater contrasts,
#'   familywise intervals and simultaneous rater classification are not provided.
#' @section Review unresolved refits:
#' Inspect `$trials` before interpreting the intervals. `FitReady` combines
#' numerical and information checks; it does not certify interval coverage.
#' New results also retain `OptimizerCode`, `NumericalReady`,
#' `InformationPositive`, `PersonQuadratureStable`, `QuadraturePoints`, `CheckPoints`,
#' `LogLikDifference`, `GradientDifference`, `EstimatedVarianceBoundary` and
#' `PersonVarianceUpperBoundary` from each refit. Together with the recorded
#' gradient and Person-variance boundary, these distinguish unresolved numerical
#' calculations from variance boundaries and missing regular studentizers.
#' A failed refit has missing additional check values and its recorded error; missing
#' checks are not passes. Older saved trial tables lack these additional fields
#' and cannot acquire them without rerunning the corresponding refits.
#' Failed checks must not be removed, selectively retried until successful or
#' relabeled as adequate interval coverage. Increasing Person quadrature can
#' address integration precision but does not qualify the rater Laplace
#' approximation, the interval method or the population assumptions.
#'
#' Increasing `nsim` alone does not resolve unavailable prediction errors.
#' For example, with 499 planned refits at 95%, 13 unresolved studentized
#' errors for one rater make both of that rater's limits infinite under the
#' type-1 completion rule. Twelve unresolved errors leave finite empirical
#' limits when all other errors and the source estimate/SE are finite. This
#' describes the calculation, not a threshold establishing accurate coverage.
#' If the unresolved fraction remains above the nominal tail probability,
#' a larger run still has unbounded limits. Inspect causes in `$trials` before
#' committing to more refits; do not discard unresolved draws.
#' @return An `mfrm_random_rater_intervals` object containing the source rater
#'   table, `intervals`, aligned error/studentized matrices, generated effects,
#'   refitted estimates/SEs, per-trial checks, warnings/errors, seeds, model
#'   settings and analysis data. `confint()` returns a two-column matrix with
#'   availability, method and Monte Carlo resolution attributes. No native
#'   pointers are stored; save with `saveRDS()`.
#' @references Chatterjee, S., Lahiri, P. and Li, H. (2008). Parametric bootstrap
#'   approximation to the distribution of EBLUP and related prediction intervals
#'   in linear mixed models. *Annals of Statistics*, 36, 1221--1245.
#'   \doi{10.1214/07-AOS512}.
#' @seealso [fit_mfrm_random_rater()], [plot.mfrm_random_rater_intervals()]
#' @examples example <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr"))
#' fit <- example$random_rater$fit
#' intervals <- example$random_rater$intervals
#' # To regenerate (requires RTMB >= 2.0 and repeated model fitting):
#' # intervals <- mfrm_random_rater_intervals(fit, nsim = 19, seed = 923701)
#' # These 19 saved trials illustrate mechanics, not accurate 2.5% tails.
#' # Every trial is retained; unresolved trials can give unbounded intervals.
#' summary(intervals)
#' confint(intervals, method = "error")
#' # Changing the level reuses the same draws, with no further fitting:
#' confint(intervals, level = .90)
#' @export
#' @order 1
mfrm_random_rater_intervals <- function(object, nsim = 499L, seed, level = .95) {
  if (!inherits(object, "mfrm_random_rater") || !isTRUE(object$checks$NumericalReady) ||
      !isTRUE(object$checks$InformationPositive)) stop("Supply a shared-rater fit with resolved numerical and information checks.", call. = FALSE)
  if (object$calibration$rater_sd <= 0 || any(!is.finite(object$raters$PredictionSE)) ||
      any(object$raters$PredictionSE <= 0)) stop("Bootstrap intervals require positive source rater SD and prediction SEs; estimated boundaries do not supply regular studentizers. Rater-SD profiling is a separate target and is unavailable at an estimated Person-variance boundary.", call. = FALSE)
  if (!is.numeric(nsim) || is.complex(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim != floor(nsim) || nsim < 2 || nsim > 100000) stop("`nsim` must be an integer from 2 to 100000.", call. = FALSE)
  if (missing(seed) || !is.numeric(seed) || is.complex(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) stop("Supply a nonnegative integer `seed`.", call. = FALSE)
  mfrm_random_rater_interval_level(level)
  with_preserved_rng_seed(seed, {
    if (!requireNamespace("RTMB", quietly = TRUE) || utils::packageVersion("RTMB") < "2.0") {
      stop("Install optional RTMB version 2.0 or later for bootstrap refitting.", call. = FALSE)
    }
    seeds <- sample.int(.Machine$integer.max, nsim)
    draws <- lapply(seeds, function(s) mfrm_random_rater_bootstrap_one(object, s))
    get_matrix <- function(name) {
      z <- do.call(rbind, lapply(draws, `[[`, name))
      dimnames(z) <- list(as.character(seq_len(nsim)), object$raters$Rater); z
    }
    out <- list(source = object, error = get_matrix("error"), studentized = get_matrix("studentized"),
      truth = get_matrix("truth"), estimates = get_matrix("estimate"), prediction_se = get_matrix("se"),
      trials = do.call(rbind, lapply(draws, `[[`, "trial")),
      settings = list(nsim = nsim, seed = seed, replicate_seeds = seeds, rng_kind = RNGkind(),
        person_sd = object$calibration$person_sd %||% 1,
        fixed_person_sd = mfrm_random_rater_fixed_person_sd(object),
        level = level, target = "Realized rater effects relative to the population mean",
        sampling = "New persons, shared raters and responses on analyzed rows", quantile_type = 1L))
    out$trials$Replicate <- seq_len(nsim)
    class(out) <- "mfrm_random_rater_intervals"
    out$intervals <- confint(out, level = level)
    out
  })
}

mfrm_random_rater_interval_level <- function(level) {
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) stop("Supply 0 < level < 1.", call. = FALSE)
}

mfrm_random_rater_generate <- function(object, seed) {
  with_preserved_rng_seed(seed, {
    input <- object$input; n <- length(input$y)
    theta <- stats::rnorm(max(input$person), sd = object$calibration$person_sd %||% 1)
    severity <- stats::rnorm(max(input$rater), sd = object$calibration$rater_sd)
    eta <- theta[input$person] - severity[input$rater]
    if (length(object$calibration$beta)) eta <- eta - as.vector(input$X %*% object$calibration$beta)
    steps <- object$calibration$steps
    logw <- outer(eta, 0:length(steps)) - matrix(c(0, cumsum(steps)), n, length(steps) + 1L, byrow = TRUE)
    w <- exp(logw - apply(logw, 1L, max)); w <- w / rowSums(w)
    u <- stats::runif(n); cumulative <- w[, 1L]; category <- integer(n)
    for (k in seq_along(steps)) {
      category <- category + as.integer(u > cumulative)
      cumulative <- cumulative + w[, k + 1L]
    }
    data <- input$data; data[[input$columns$score]] <- input$score_levels[category + 1L]
    list(data = data, truth = setNames(severity, input$levels[[input$columns$rater]]), theta = theta)
  })
}

mfrm_random_rater_bootstrap_one <- function(object, seed) {
  generated <- mfrm_random_rater_generate(object, seed)
  warnings <- character(); error_message <- ""; columns <- object$input$columns
  fit <- tryCatch(withCallingHandlers(fit_mfrm_random_rater(generated$data,
    columns$person, columns$rater, columns$score, columns$facets, object$input$score_levels,
    rater_sd = object$settings$fixed_rater_sd, quad_points = object$settings$quad_points,
    maxit = object$settings$maxit, person_sd = mfrm_random_rater_fixed_person_sd(object),
    person_variance_max = object$settings$person_variance_max %||% 16), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) { error_message <<- conditionMessage(e); NULL })
  labels <- object$raters$Rater
  estimate <- se <- error <- studentized <- rep(NA_real_, length(labels))
  ready <- !is.null(fit) && isTRUE(fit$checks$NumericalReady) && isTRUE(fit$checks$InformationPositive)
  boundary <- !is.null(fit) && (isTRUE(fit$checks$EstimatedVarianceBoundary) ||
    isTRUE(fit$checks$EstimatedPersonVarianceBoundary))
  if (ready) {
    at <- match(labels, fit$raters$Rater)
    if (anyNA(at)) stop("A bootstrap refit changed the rater identities.", call. = FALSE)
    estimate <- fit$raters$Estimate[at]; se <- fit$raters$PredictionSE[at]
    error <- unname(generated$truth[labels]) - estimate
    usable <- is.finite(se) & se > 0
    studentized[usable] <- error[usable] / se[usable]
  }
  check <- function(name) if (is.null(fit)) NA else fit$checks[[name]] %||% NA
  list(error = error, studentized = studentized, truth = unname(generated$truth[labels]),
    estimate = estimate, se = se, trial = data.frame(Seed = seed, FitReady = ready,
      EstimatedBoundary = boundary, RaterSD = if (!is.null(fit)) fit$calibration$rater_sd else NA_real_,
      PersonSD = if (!is.null(fit)) fit$calibration$person_sd %||% 1 else NA_real_,
      EstimatedPersonVarianceBoundary = !is.null(fit) && isTRUE(fit$checks$EstimatedPersonVarianceBoundary),
      MaxGradient = if (!is.null(fit)) fit$checks$MaxGradient else NA_real_,
      OptimizerCode = check("OptimizerCode"), NumericalReady = check("NumericalReady"),
      InformationPositive = check("InformationPositive"),
      PersonQuadratureStable = check("PersonQuadratureStable"),
      QuadraturePoints = check("QuadraturePoints"), CheckPoints = check("CheckPoints"),
      LogLikDifference = check("LogLikDifference"), GradientDifference = check("GradientDifference"),
      EstimatedVarianceBoundary = check("EstimatedVarianceBoundary"),
      PersonVarianceUpperBoundary = check("PersonVarianceUpperBoundary"),
      Error = error_message, Warnings = paste(unique(warnings), collapse = " | ")))
}

#' @rdname mfrm_random_rater_intervals
#' @export
confint.mfrm_random_rater_intervals <- function(object, parm = NULL, level = .95,
                                                method = c("studentized", "error"), ...) {
  rlang::check_dots_empty(); mfrm_random_rater_interval_level(level); method <- match.arg(method)
  labels <- object$source$raters$Rater
  if (is.null(parm)) parm <- labels
  if (!is.character(parm) || !length(parm) || anyNA(parm) || anyDuplicated(parm) ||
      !all(parm %in% labels)) stop("`parm` must contain distinct fitted rater IDs.", call. = FALSE)
  at <- match(parm, labels); roots <- object[[method]][, at, drop = FALSE]
  source_se <- if (method == "studentized") object$source$raters$PredictionSE[at] else rep(1, length(at))
  lower <- upper <- roots; missing <- !is.finite(roots)
  lower[missing] <- -Inf; upper[missing] <- Inf
  tails <- (1 - level) / 2
  bounds <- cbind(Lower = apply(lower, 2L, stats::quantile, probs = tails, type = 1, names = FALSE),
    Upper = apply(upper, 2L, stats::quantile, probs = 1 - tails, type = 1, names = FALSE))
  bounds <- object$source$raters$Estimate[at] + source_se * bounds
  rownames(bounds) <- parm
  attr(bounds, "availability") <- data.frame(Rater = parm, Planned = nrow(roots),
    Known = colSums(!missing), Unresolved = colSums(missing), Finite = apply(is.finite(bounds), 1L, all))
  attr(bounds, "method") <- method; attr(bounds, "level") <- level
  attr(bounds, "expected_tail_draws") <- nrow(roots) * tails
  attr(bounds, "note") <- "Pointwise model-based bootstrap; unresolved roots widen limits, possibly to infinity. Coverage and Monte Carlo tail accuracy require separate assessment."
  bounds
}

#' @rdname mfrm_random_rater_intervals
#' @export
summary.mfrm_random_rater_intervals <- function(object, ...) {
  list(intervals = object$intervals, availability = attr(object$intervals, "availability"),
    trials = c(Planned = nrow(object$trials), FitReady = sum(object$trials$FitReady),
      EstimatedBoundary = sum(object$trials$EstimatedBoundary)), settings = object$settings)
}

#' @rdname mfrm_random_rater_intervals
#' @export
print.mfrm_random_rater_intervals <- function(x, ...) {
  cat("Pointwise bootstrap prediction intervals for observed random raters\n")
  print(x$intervals[, , drop = FALSE])
  cat("Planned refits:", nrow(x$trials), " Ready:", sum(x$trials$FitReady),
    " Variance boundaries:", sum(x$trials$EstimatedBoundary), "\n")
  cat("Expected draws per nominal tail:", attr(x$intervals, "expected_tail_draws"), "\n")
  cat("Model-based candidate; finite-sample coverage is not guaranteed.\n")
  invisible(x)
}

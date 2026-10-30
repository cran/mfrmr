#' Predict category probabilities for observed or replacement raters
#'
#' Evaluate a shared-rater RSM at explicitly supplied ability values.
#'
#' @param object A result from [fit_mfrm_random_rater()].
#' @param newdata Data frame with the fitted rater column and every fixed-facet
#'   column. Observed-rater labels must be known; replacement-rater labels must
#'   be new. Fixed-facet levels must always be known.
#' @param ability A finite numeric value for all rows, one value per row, or
#'   the name of a numeric column of `newdata`. These are specified abilities
#'   on the fitted logit scale (mean ability zero, Rasch slope one), not estimated
#'   person scores. A value of one is one logit, not one population SD.
#' @param rater `"observed"` uses the approximate conditional distribution of
#'   an observed rater given the calibration responses. `"new"` integrates a
#'   replacement rater from the estimated population distribution. The target
#'   is explicit; an unfamiliar label is never silently treated as observed.
#' @param quad_points Normal quadrature points for the prediction integral,
#'   from 7 to 121; default 41. Compare orders when needed.
#' @param ... Unused.
#'
#' @return A list with `probabilities` (rows by categories), `expected_scores`,
#'   `newdata`, `ability`, `settings` and matching `source` metadata for
#'   [mfrm_results()]. No model is refitted. Save this list
#'   together with the fitted model and prediction inputs. The settings record
#'   the effect-sharing unit and whether integration uses the observed rater's
#'   conditional distribution or the replacement-rater population. These are
#'   preserved in the prediction settings table of [mfrm_results()].
#' @details Predictions condition on estimated calibration and supplied ability.
#'   They do not propagate calibration-estimation uncertainty or uncertainty
#'   about ability. Observed-rater integration uses the conditional Laplace
#'   mean/mode and covariance, not the calibration-adjusted prediction SE.
#'   New-rater integration uses the population SD: substituting severity zero
#'   generally produces different probabilities.
#'
#'   Rows are marginal probabilities for individual future ratings, not a joint
#'   distribution or an interval for their average. Ratings with the same new
#'   rater ID share one random severity in the model; multiplying these row
#'   probabilities would discard that dependence. This function does not score
#'   latent abilities from responses or impute missing assigned scores.
#'   Use [score_mfrm_random_rater()] for conditional Person scoring from a
#'   complete joint response roster.
#'   Supplying a score column in `newdata` does not update rater distributions
#'   or estimate abilities. Use the scoring function and its complete-roster
#'   rules when the input is a new batch of observed responses.
#' @export
predict.mfrm_random_rater <- function(object, newdata, ability,
                                      rater = c("observed", "new"), quad_points = 41L, ...) {
  rlang::check_dots_empty()
  rater <- match.arg(rater)
  if (!isTRUE(object$checks$NumericalReady) || !isTRUE(object$checks$InformationPositive)) {
    stop("Resolve the fit's numerical and information checks before prediction.", call. = FALSE)
  }
  if (!is.numeric(quad_points) || is.complex(quad_points) || length(quad_points) != 1L ||
      !is.finite(quad_points) || quad_points != floor(quad_points) || quad_points < 7 || quad_points > 121) {
    stop("Use integer quad_points from 7 to 121.", call. = FALSE)
  }
  columns <- object$input$columns
  required <- c(columns$rater, columns$facets)
  if (!is.data.frame(newdata) || !nrow(newdata) || anyDuplicated(names(newdata)) || !all(required %in% names(newdata))) {
    stop("`newdata` must contain the fitted rater and fixed-facet columns.", call. = FALSE)
  }
  n <- nrow(newdata)
  if (is.character(ability) && length(ability) == 1L && !is.na(ability) && ability %in% names(newdata)) ability <- newdata[[ability]]
  if (!is.numeric(ability) || is.complex(ability) || !is.null(dim(ability)) ||
      !length(ability) %in% c(1L, n) || any(!is.finite(ability))) stop("Supply finite numeric abilities, one per row or one common value.", call. = FALSE)
  ability <- rep(ability, length.out = n)
  ids <- lapply(newdata[required], function(x) {
    if (!is.atomic(x) || !is.null(dim(x)) || anyNA(x) || any(!nzchar(trimws(as.character(x))))) {
      stop("Prediction identifiers must be complete and nonempty.", call. = FALSE)
    }
    as.character(x)
  })
  known <- object$input$levels[[columns$rater]]
  observed <- match(ids[[columns$rater]], known)
  if (rater == "observed" && anyNA(observed)) stop("Observed-rater prediction requires known rater IDs.", call. = FALSE)
  if (rater == "new" && any(!is.na(observed))) stop("Replacement-rater prediction requires new rater IDs.", call. = FALSE)
  location <- rep(0, n); at <- 0L
  for (facet in columns$facets) {
    basis <- object$input$basis[[facet]]
    index <- match(ids[[facet]], rownames(basis))
    if (anyNA(index)) stop("Prediction contains an unknown fixed-facet level: ", facet, call. = FALSE)
    beta <- object$calibration$beta[at + seq_len(ncol(basis))]; at <- at + ncol(basis)
    location <- location + as.vector(basis[index, , drop = FALSE] %*% beta)
  }
  mean <- if (rater == "observed") object$rater_mode[observed] else rep(0, n)
  sd <- if (rater == "observed") sqrt(diag(object$conditional_rater_covariance))[observed] else rep(object$calibration$rater_sd, n)
  rule <- gauss_hermite_normal(quad_points)
  steps <- object$calibration$steps; categories <- object$input$score_levels
  probabilities <- matrix(0, n, length(categories), dimnames = list(NULL, as.character(categories)))
  for (q in seq_along(rule$nodes)) {
    eta <- ability - location - mean - sd * rule$nodes[q]
    logw <- outer(eta, 0:length(steps)) - matrix(c(0, cumsum(steps)), n, length(categories), byrow = TRUE)
    w <- exp(logw - apply(logw, 1L, max))
    probabilities <- probabilities + rule$weights[q] * w / rowSums(w)
  }
  list(probabilities = probabilities, source = mfrm_extended_prediction_source(object),
    expected_scores = data.frame(Row = seq_len(n), Rater = ids[[columns$rater]],
      Ability = ability, ExpectedScore = as.vector(probabilities %*% categories)),
    newdata = newdata, ability = ability, settings = list(rater = rater,
      effect_sharing = "One severity effect per rater ID, shared across Persons; row probabilities are not independent joint predictions",
      rater_integration = if (rater == "observed")
        "Conditional Laplace distribution given calibration responses; conditional rater covariance" else
        "Normal replacement-rater population; fitted rater SD",
      calibration_uncertainty = "Not included; fitted calibration held fixed",
      ability_uncertainty = "Not included; specified ability values",
      row_target = "Marginal probabilities, not a joint rating distribution", quad_points = quad_points))
}

#' @rdname fit_mfrm_random_rater
#' @details `plot()` method: Plot observed-rater severity under a shared-rater RSM
#' @inheritParams plot.mfrm_testlet_scores
#' @seealso [plot.mfrm_testlet_scores()] for view selection and accessible
#'   output. The distribution view summarizes the observed conditional modes,
#'   not the assumed rater population.
#' @param x A result from [fit_mfrm_random_rater()].
#'
#' Also for `summary()`: A fitted `mfrm_random_rater` object.
#' @param draw Logical; `FALSE` returns plot data without opening a device.
#' @param intervals `"none"` (default) displays point estimates without
#'   automatic prediction intervals. `"normal"` explicitly requests the
#'   first-order approximation from `confint(x, parm = "raters")`; nominal
#'   coverage is not established. Required for the interval-width precision
#'   view and sorting by interval width (`sort = "uncertainty"`).
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object with the rater table and settings.
#'   Use [plot_data()] for custom graphics.
#' @details For `plot()`: Points are conditional modes of severity relative to the rater
#'   population mean zero; positive values mean stricter ratings. No individual
#'   rater intervals are supplied automatically, including for earlier saved
#'   fits. Explicit normal whiskers are approximate 95% prediction intervals,
#'   not qualified rater classifications. At an estimated zero variance or with
#'   unresolved numerical checks, even requested whiskers are unavailable.
#' @export
plot.mfrm_random_rater <- function(x, draw = TRUE,
    style = c("interval", "precision", "distribution"),
    sort = c("input", "estimate", "uncertainty"), palette = c("accessible", "mono"),
    title = NULL, caption = NULL, show_title = FALSE, show_notes = FALSE,
    show_labels = TRUE, reference = 0, text_scale = 1, point_size = 2.5,
    intervals = c("none", "normal"), ...) {
  rlang::check_dots_empty()
  intervals <- match.arg(intervals); style <- match.arg(style); sort <- match.arg(sort)
  if ((style == "precision" || sort == "uncertainty") && intervals == "none") stop("Interval-width display or sorting requires an explicit intervals = 'normal' approximation, or plot a saved bootstrap interval result.", call. = FALSE)
  note <- if (isTRUE(x$checks$EstimatedVarianceBoundary) ||
      isTRUE(x$checks$EstimatedPersonVarianceBoundary)) "Estimated variance is zero; regular intervals unavailable" else
    if (!x$checks$NumericalReady || !x$checks$InformationPositive) "Numerical review required; intervals unavailable" else
      "Explicit first-order 95% normal approximation; nominal coverage is not established"
  if (intervals == "none") note <- if (!isTRUE(x$checks$NumericalReady) || !isTRUE(x$checks$InformationPositive))
    "Point estimates only; numerical review required" else
      "Point estimates only; individual-rater interval coverage is not established"
  settings <- x$settings; settings$rater_intervals <- intervals
  mfrm_extended_estimate_plot(mfrm_random_rater_table(x, intervals), x$raters$Rater, "random_rater_severity",
    "Shared-rater severity", "Rater severity (logits)", note, draw, settings,
    style = style, sort = sort, palette = palette, title = title, caption = caption,
    show_title = show_title, show_notes = show_notes, show_labels = show_labels,
    reference = reference, text_scale = text_scale, point_size = point_size, extra = list(checks = x$checks))
}

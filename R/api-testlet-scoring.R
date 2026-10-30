#' Score Persons using a fixed testlet calibration
#'
#' Recompute Person scores from all supplied rows, integrating local effects
#' within each Person/testlet block and conditioning on fitted calibration.
#' @param object A ready [fit_mfrm_testlet()] result.
#' @param newdata Long-format ratings with the fitted column names. `NULL`
#'   reuses the source assignment, including missing-score rows. New Person and
#'   testlet labels are allowed; every fixed-facet level must be known.
#'   If the testlet column is also a fixed facet, its labels must therefore
#'   be known. For example, with `testlet = "Rater"` and `"Rater"` in `facets`,
#'   this function cannot score a new rater whose fixed severity was not fitted.
#' @param level Conditional equal-tail posterior interval probability; default
#'   0.95. These are not calibration-adjusted frequentist confidence intervals.
#' @param quad_points Local-effect quadrature order, 7 to 121. Default uses the
#'   fitting order; each Person is also checked at `2 * quad_points + 1`.
#' @param missing `"fail"` or `"omit"`. Default reuses the fit's missing-score
#'   policy for source replay and requires explicit omission for new data.
#' @param ... Unused.
#' @details All supplied rows for a Person are scored jointly. This does not
#'   append to cached responses or condition on a stored local-effect mode.
#'   The saved settings retain this integration, effect-sharing and roster
#'   interpretation. [mfrm_results()] retains the complete `scoring_data` as
#'   its `scoring_roster` table, including explicitly omitted scores, alongside
#'   observed block counts. Column names and row order are preserved.
#'   To add ratings for an existing Person, supply that Person's complete set
#'   of ratings once, with memberships that correctly identify effects shared
#'   within that set. Scores from separate calls are not automatically linked
#'   within a Person or testlet. No fixed effects are re-estimated.
#'   Unequal block sizes do not imply equal block weights. The influence of an
#'   additional response depends on its block and the fitted model; integrating
#'   local effects is not a general correction for assignment or content bias.
#'   See `vignette("mfrmr-testlet-applications")` for a complete-roster
#'   comparison of one-point changes in five- versus two-criterion tasks.
#'
#'   EAP, posterior SD and intervals use continuous ability integration. Interval
#'   endpoints invert the continuous CDF rather than a finite quadrature-grid
#'   CDF. Local-effect integration is compared at two orders; errors or unresolved
#'   differences retain a row with unavailable scores and an explanation.
#'   An explicitly omitted, entirely missing Person returns the fitted or
#'   specified mean-zero normal ability population
#'   with status `"prior_only"`, not evidence of measured average ability.
#'
#'   Calibration is held fixed, including ability and testlet variances. If the
#'   fitted ability variance is zero, all requested Persons retain unavailable
#'   rows with a reason; no zero-width ability interval is reported. At an estimated
#'   zero variance this is scoring under that fitted submodel, not evidence that
#'   dependence is absent. Prior/posterior terminology refers to latent Person
#'   scoring; calibration remains frequentist MML. A 480-dataset same-source
#'   simulation estimating both variances, replayed after correcting numerical
#'   start selection, found Person-interval coverage of 94.4% with 120 Persons
#'   and 91.7% with 24 Persons under balanced positive local dependence;
#'   coverage was 93.3% in one
#'   sparse, unequal-block condition with 120 Persons. The original numerical
#'   failures and their repair remain separately documented. Accounting for
#'   dependence improved coverage relative to ordinary RSM, but did not consistently improve EAP
#'   mean squared error. See `vignette("mfrmr-testlets")` for the comparison,
#'   Monte Carlo uncertainty and numerical-selection correction.
#'   These intervals neither correct calibration-estimation uncertainty nor
#'   guarantee 95% coverage for each fixed ability. Person contrasts and
#'   simultaneous decisions are not provided by this route.
#' @return An `mfrm_testlet_scores` object with one table row per requested
#'   Person, omitted-row accounting, block sizes, settings, calibration and
#'   prediction-source metadata for [mfrm_results()].
#'   `summary()` reports all statuses; `plot()` returns [plot_data()] when
#'   `draw = FALSE`. Saved scoring results need no live optimizer.
#'   Keep `table$Status` and `table$Reason` with exported scores: `prior_only`
#'   describes the assumed population without observed responses for that Person;
#'   `unavailable` withholds the estimate, SD and bounds with a recorded cause.
#' @param x For `predict.mfrm_testlet()`: A scoring result.
#'
#' For `plot()`: A result from [predict.mfrm_testlet()].
#' @param persons Distinct Person IDs to return, or `NULL` for all. This
#'   selects outputs while retaining the complete supplied scoring roster.
#' @seealso [fit_mfrm_testlet()], [plot.mfrm_testlet_scores()], [mfrm_results()]
#' @examples # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("predict.mfrm_testlet", run.dontrun = TRUE)
#' \dontrun{
#' example <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr"))
#' fit <- example$testlet$fit
#' scores <- example$testlet$scores
#' # To recompute these conditional scores:
#' # ratings <- load_mfrmr_data("example_core")
#' # scores <- predict(fit, persons = as.character(unique(ratings$Person)[1:4]))
#' # Selecting outputs retains the complete source roster for model maps.
#' scores$table
#' plot(scores)
#' }
#' @export
#' @order 1
predict.mfrm_testlet <- function(object, newdata = NULL, level = .95,
    quad_points = object$settings$quad_points,
    missing = if (is.null(newdata)) object$settings$missing else "fail", persons = NULL, ...) {
  rlang::check_dots_empty()
  if (!isTRUE(object$checks$NumericalReady) || !isTRUE(object$checks$InformationPositive)) {
    stop("Resolve the testlet fit's numerical and information checks before scoring.", call. = FALSE)
  }
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) stop("Supply 0 < level < 1.", call. = FALSE)
  if (!is.numeric(quad_points) || is.complex(quad_points) || length(quad_points) != 1L ||
      !is.finite(quad_points) || quad_points != floor(quad_points) || quad_points < 7 || quad_points > 121) stop("Use integer quad_points from 7 to 121.", call. = FALSE)
  missing <- match.arg(missing, c("fail", "omit")); columns <- object$input$columns
  data <- if (is.null(newdata)) object$input$assigned_data %||% object$input$data else newdata
  if (is.null(newdata) && is.null(object$input$assigned_data) && length(object$input$omitted_rows)) stop(
    "The saved fit lacks its missing-score roster; supply the complete roster as newdata.", call. = FALSE)
  input <- mfrm_testlet_data(data, columns$person, columns$score, columns$testlet,
    columns$facets, object$input$score_levels, missing, reference = object$input)
  person_variance <- mfrm_testlet_person_variance(input, object$parameters)
  persons <- persons %||% input$persons
  if (!is.character(persons) || !length(persons) || anyNA(persons) || anyDuplicated(persons) || !all(persons %in% input$persons)) stop("Choose distinct Person IDs in the scoring roster.", call. = FALSE)
  person_sd <- sqrt(person_variance)
  rules <- lapply(c(quad_points, 2 * quad_points + 1), gauss_hermite_normal)
  table <- lapply(match(persons, input$persons), function(p) {
    observed <- sum(lengths(input$groups[[p]])); reason <- ""; difference <- NA_real_
    if (person_variance == 0) {
      value <- c(Estimate = NA_real_, ConditionalSD = NA_real_, Lower = NA_real_, Upper = NA_real_)
      status <- "unavailable"
      reason <- "The fitted ability variance is zero; individual ability scores are withheld."
    } else if (observed == 0) {
      alpha <- (1 - level) / 2
      value <- c(Estimate = 0, ConditionalSD = person_sd, Lower = person_sd * stats::qnorm(alpha),
        Upper = person_sd * stats::qnorm(alpha, lower.tail = FALSE))
      status <- "prior_only"; difference <- 0
      reason <- "No observed responses for this Person; the estimate, SD and interval describe only the assumed population, not measured individual ability."
    } else {
      value <- tryCatch(withCallingHandlers({
        low <- mfrm_testlet_person_interval(input, object$parameters, p, rules[[1]], level)
        high <- mfrm_testlet_person_interval(input, object$parameters, p, rules[[2]], level)
        delta <- abs(high - low); difference <- max(delta)
        if (any(!is.finite(delta)) || any(delta > c(1e-7, 1e-7, 1e-5, 1e-5, 1e-6))) {
          stop("Local-effect quadrature is unresolved; increase quad_points and review the scoring pattern.")
        }
        high[1:4]
      }, warning = function(w) stop(conditionMessage(w), call. = FALSE)),
      error = function(e) { reason <<- conditionMessage(e); NULL })
      status <- if (is.null(value)) "unavailable" else "available_conditional"
      if (is.null(value)) value <- setNames(rep(NA_real_, 4), c("Estimate", "ConditionalSD", "Lower", "Upper"))
    }
    data.frame(Person = input$persons[p], Observed = observed, Testlets = length(input$groups[[p]]),
      as.list(value), Status = status, IntegrationDifference = difference, Reason = reason, row.names = NULL)
  })
  structure(list(table = do.call(rbind, table), blocks = input$blocks, scoring_data = input$assigned_data,
    data_usage = c(Input = input$input_rows, Observed = length(input$y), Omitted = length(input$omitted_rows)),
    omitted_rows = input$omitted_rows, calibration = object$calibration,
    source = mfrm_extended_prediction_source(object),
    settings = list(level = level, quad_points = quad_points, check_points = 2 * quad_points + 1,
      missing = missing, calibration_uncertainty = FALSE,
      effect_sharing = "One local effect per Person/testlet pair; repeated labels across Persons are different effects",
      local_effect_integration = "Integrated within each Person/block, then over ability; no stored local-effect mode is reused",
      roster = "Supplied table replaces source; output selection retains the full scoring roster",
      estimated_variance_boundary = object$checks$EstimatedVarianceBoundary,
      estimated_person_variance_boundary = isTRUE(object$checks$EstimatedPersonVarianceBoundary),
      person_variance = person_variance,
      target = "Person ability conditional on calibration and the fitted or specified normal population")), class = "mfrm_testlet_scores")
}

mfrm_testlet_person_interval <- function(input, par, person, rule, level) {
  person_sd <- sqrt(mfrm_testlet_person_variance(input, par))
  if (person_sd == 0) stop("Individual scoring requires a positive ability variance.")
  # Integrate in standard-normal coordinates for both narrow and broad priors.
  # Moments and CDF inversion revisit identical points; reuse only exact values
  # within this Person/calibration/rule, without interpolating or rounding nodes.
  cache <- new.env(parent = emptyenv())
  loglikelihood <- function(z) {
    keys <- sprintf("%a", z)
    fresh <- !vapply(keys, exists, logical(1), envir = cache, inherits = FALSE) & !duplicated(keys)
    if (any(fresh)) {
      values <- as.vector(mfrm_testlet_kernel(input, par, person_sd * z[fresh], rule,
        persons = person)$loglik)
      new_keys <- keys[fresh]
      for (i in seq_along(values)) assign(new_keys[i], values[i], envir = cache)
    }
    vapply(keys, get, numeric(1), envir = cache, inherits = FALSE, USE.NAMES = FALSE)
  }
  mfrm_person_posterior_interval(loglikelihood, person_sd, level)
}

mfrm_person_posterior_interval <- function(loglikelihood, person_sd, level, cut_points = NULL) {
  # Center the likelihood before adding the standard-normal prior, including
  # during mode finding; a common likelihood constant must not quantize it.
  offset <- loglikelihood(0)
  if (length(offset) != 1L || !is.finite(offset)) stop("Conditional log likelihood is nonfinite.")
  logdensity <- function(z) (loglikelihood(z) - offset) + stats::dnorm(z, log = TRUE)
  mode <- stats::optimize(function(x) -logdensity(x), c(-20, 20), tol = 1e-8)$minimum
  if (abs(mode) > 19.9) stop("The conditional posterior mode is outside the scoring search range.")
  peak_likelihood <- loglikelihood(mode)
  peak_prior <- stats::dnorm(mode, log = TRUE)
  density <- function(theta) exp((loglikelihood(theta) - peak_likelihood) +
    (stats::dnorm(theta, log = TRUE) - peak_prior))
  integral <- function(f, lower, upper) {
    z <- stats::integrate(f, lower, upper, rel.tol = 1e-9, abs.tol = 1e-12, subdivisions = 500L)
    if (!is.finite(z$value) || z$message != "OK") stop("Continuous posterior integration failed.")
    z$value
  }
  full <- function(f) integral(f, -Inf, mode) + integral(f, mode, Inf)
  left_mass <- integral(density, -Inf, mode)
  right_mass <- integral(density, mode, Inf)
  mass <- left_mass + right_mass
  if (mass <= 0) stop("Continuous posterior mass is unavailable.")
  mean_shift <- full(function(x) (x - mode) * density(x)) / mass
  variance <- full(function(x) (x - mode)^2 * density(x)) / mass - mean_shift^2
  if (!is.finite(variance) || variance <= 0) stop("Continuous posterior variance is unavailable.")
  alpha <- (1 - level) / 2
  # Anchor ordinary CDF evaluations at the mode, reusing the expensive tails.
  # Very small tail requests retain direct integration to avoid cancellation.
  lower_tail <- function(x) {
    if (alpha < 1e-5) return(integral(density, -Inf, x) / mass)
    if (x < mode) (left_mass - integral(density, x, mode)) / mass else
      (left_mass + integral(density, mode, x)) / mass
  }
  upper_tail <- function(x) {
    if (alpha < 1e-5) return(integral(density, x, Inf) / mass)
    if (x > mode) (right_mass - integral(density, mode, x)) / mass else
      (right_mass + integral(density, x, mode)) / mass
  }
  span <- max(8, 10 * sqrt(variance))
  for (i in 1:6) {
    bounds <- mode + c(-span, span)
    if (lower_tail(bounds[1]) < alpha && upper_tail(bounds[2]) < alpha) break
    span <- 2 * span
  }
  lower <- stats::uniroot(function(x) lower_tail(x) - alpha, bounds, tol = 1e-8, check.conv = TRUE)$root
  upper <- stats::uniroot(function(x) upper_tail(x) - alpha, bounds, tol = 1e-8, check.conv = TRUE)$root
  if (max(abs(c(lower_tail(lower), upper_tail(upper)) - alpha)) > min(1e-7, alpha * 1e-4)) stop("Continuous posterior tail accuracy is unresolved.")
  answer <- c(Estimate = person_sd * (mode + mean_shift), ConditionalSD = person_sd * sqrt(variance),
    Lower = person_sd * lower, Upper = person_sd * upper,
    LogNormalizer = peak_likelihood + (log(mass) + peak_prior))
  if (!is.null(cut_points)) answer <- c(answer,
    mfrm_cut_probabilities(function(a,b) integral(density,a,b), cut_points, mass, mode=mode))
  answer
}

#' @rdname predict.mfrm_testlet
#' @export
summary.mfrm_testlet_scores <- function(object, ...) {
  list(scores = object$table, status_counts = table(factor(object$table$Status,
    levels = c("available_conditional", "prior_only", "unavailable"))),
    data_usage = object$data_usage, settings = object$settings)
}

#' @rdname predict.mfrm_testlet
#' @export
print.mfrm_testlet_scores <- function(x, ...) {
  cat("Conditional Person scores from a fixed testlet calibration\n")
  print(x$table, row.names = FALSE)
  cat("Intervals exclude calibration-estimation uncertainty.\n")
  invisible(x)
}

#' @rdname fit_mfrm_testlet
#' @details `plot()` method: Plot testlet-model fixed-facet estimates
#' @inheritParams plot.mfrm_testlet_scores
#' @param x An object from [fit_mfrm_testlet()].
#'
#' Also for `fit_mfrm_testlet()`: An `mfrm_testlet` result.
#' @param facet A fixed facet to display; defaults to the first fixed facet.
#' @param draw `FALSE` returns plot data without opening a device.
#' @param intervals `"none"` (default) shows points only; `"normal"`
#'   explicitly requests observed-information normal-approximation bounds.
#'   Interval-width views or sorting require this explicit selection.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object; use [plot_data()].
#' @details For `plot()`: Requested whiskers are observed-information normal approximations
#'   for calibration parameters; finite-sample coverage is not established. Missing
#'   intervals remain missing; points do not classify rater quality.
#'   [plot_data()] and [plot_data_components()] provide the same table without
#'   drawing. [as_ggplot()] preserves the interval meaning and missing whiskers.
#' @seealso [mfrmr_workflow_methods], [mfrmr_output_guide()],
#'   [mfrmr_interval_guide()], [plot.mfrm_testlet_scores()]
#' @export
plot.mfrm_testlet <- function(x, facet = x$input$columns$facets[1], draw = TRUE,
    style = c("interval", "precision", "distribution"),
    sort = c("input", "estimate", "uncertainty"), palette = c("accessible", "mono"),
    title = NULL, caption = NULL, show_title = FALSE, show_notes = FALSE,
    show_labels = TRUE, reference = 0, text_scale = 1, point_size = 2.5,
    intervals = c("none", "normal"), level = .95, ...) {
  rlang::check_dots_empty()
  intervals <- match.arg(intervals); style <- match.arg(style); sort <- match.arg(sort)
  if ((style == "precision" || sort == "uncertainty") && intervals == "none") stop("Interval-width display or sorting requires explicit intervals = 'normal'.", call. = FALSE)
  if (!is.character(facet) || length(facet) != 1L || is.na(facet) || !facet %in% x$input$columns$facets) stop("Choose one fitted fixed facet.", call. = FALSE)
  tab <- mfrm_extended_calibration_table(x, intervals, level)
  tab <- tab[tab$Parameter == "Fixed facet" & tab$Facet == facet, , drop = FALSE]
  note <- if (intervals == "none") "Point estimates only; calibration bounds not requested" else
    paste0("Explicit ", format(100 * level, trim = TRUE, digits = 15), "% normal approximation; nominal coverage not established | missing whiskers: unavailable")
  mfrm_extended_estimate_plot(tab, tab$Level, "testlet_calibration", paste("Fixed", facet, "effects"),
    "Facet severity (logits)", note, draw,
    list(facet = facet, checks = x$checks, calibration_intervals = intervals, level = level), style = style, sort = sort, palette = palette, title = title, caption = caption,
    show_title = show_title, show_notes = show_notes, show_labels = show_labels,
    reference = reference, text_scale = text_scale, point_size = point_size)
}

#' @rdname predict.mfrm_testlet
#' @details `plot()` method: Plot conditional Person scores from a testlet model
#' @param style `"interval"` keeps every labeled row; `"precision"` plots
#'   estimates against interval width (upper minus lower bound); `"distribution"`
#'   shows the empirical cumulative distribution of finite point estimates,
#'   excluding prior-only and unavailable scores. No density is reconstructed.
#' @param sort Row order for interval plots: input order, increasing estimate,
#'   or increasing interval width. Missing values come last; ties retain input
#'   order. Sorting is descriptive, not a test of differences.
#' @param palette `"accessible"` uses blue and dark orange with filled/open
#'   symbols; `"mono"` uses dark grey with the same symbol distinctions.
#' @param title,caption `NULL` uses the default text; a string replaces it;
#'   `""` removes it. Interpretation notes remain in the saved plot data.
#' @param show_title,show_notes Show the title and caption, respectively.
#'   `FALSE` hides text without deleting its metadata.
#' @param show_labels Show IDs on interval and precision plots. Distribution
#'   plots summarize estimates and do not label individual IDs. For crowded
#'   precision plots use `FALSE` and consult the retained table.
#' @param reference Vertical reference line in logits; `NULL` omits it.
#'   This is an orientation aid, not a quality threshold.
#' @param text_scale Positive multiplier for text size.
#' @param point_size Positive point size in ggplot millimetres; base graphics
#'   use the corresponding relative size (2.5 is the default).
#' @section Choosing a view:
#' For `plot()`:
#'
#' Interval plots answer where estimates lie and how uncertain they are.
#'   Precision plots help identify estimates with wider intervals; width is
#'   neither a fit statistic nor a reliability coefficient. Finite intervals
#'   are required and excluded rows remain in `display_data` with reasons.
#'   Distribution plots summarize the selected point estimates, which are
#'   affected by shrinkage and the selected roster. They do not estimate the
#'   latent population distribution, show posterior densities, or imply
#'   independent observations.
#' @section Accessible and reusable output:
#' For `plot()`:
#'
#' Colours are supplemented by point shapes; open circles mean prior only.
#'   Use `palette = "mono"` for monochrome reproduction. [plot_data()] retains
#'   the complete source `table`, `display_data` with inclusion reasons,
#'   `alt_text`, interpretation `notes`, and display settings. Supply the text
#'   alternative and table alongside exported figures; image files alone do
#'   not automatically expose that information to screen readers.
#'   [as_ggplot()] preserves the view and display controls, supplies
#'   `labs(alt = ...)`, and allows further `labs()`, `theme()` and scale edits.
#'   Rendering does not fit, score or resample. Use the same controls through
#'   `plot(results, type = "scores", style = "precision")`.
#' @param draw `FALSE` returns plot data without opening a device.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object retaining all requested Persons,
#'   their statuses, interval endpoints and settings; use [plot_data()].
#' @details For `plot()`: In the default interval view, empty rows remain labeled. Open circles mark prior-only results;
#'   filled points are response-based conditional scores. The intervals exclude
#'   calibration-estimation uncertainty and are not tests of Person differences.
#'   [as_ggplot()] preserves all labeled rows, prior-only symbols and the
#'   conditional-interval note. [plot_data()] retains statuses and reasons for
#'   custom graphics or table export; conversion does not rerun scoring.
#' @seealso [plot.mfrm_testlet()], [plot_data_components()], [mfrmr_interval_guide()]
#' @export
plot.mfrm_testlet_scores <- function(x, draw = TRUE,
    style = c("interval", "precision", "distribution"),
    sort = c("input", "estimate", "uncertainty"), palette = c("accessible", "mono"),
    title = NULL, caption = NULL, show_title = FALSE, show_notes = FALSE,
    show_labels = TRUE, reference = 0, text_scale = 1, point_size = 2.5, ...) {
  rlang::check_dots_empty()
  mfrm_extended_estimate_plot(x$table, x$table$Person, "testlet_scores", "Conditional Person scores",
    "Person ability (logits)", sprintf("%.0f%% conditional intervals | calibration fixed | open: prior only", 100 * x$settings$level),
    draw, x$settings, style = style, sort = sort, palette = palette, title = title, caption = caption,
    show_title = show_title, show_notes = show_notes, show_labels = show_labels,
    reference = reference, text_scale = text_scale, point_size = point_size)
}

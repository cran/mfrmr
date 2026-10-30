#' Score Persons while integrating shared random raters
#'
#' Compute conditional Person scores using the complete joint scoring roster,
#' with calibration held fixed and shared raters integrated jointly.
#' @param object A numerically ready [fit_mfrm_random_rater()] result.
#' @param newdata Complete long-format scoring roster, with the fitted columns.
#'   `NULL` uses the source roster. A supplied table replaces the source roster;
#'   it is not appended to it. New Person and rater IDs are allowed, but fixed
#'   facet levels must be known. Include all responses intended to inform the
#'   shared raters, including responses of Persons not selected for display.
#' @param persons Character vector of distinct Person IDs to return; `NULL`
#'   returns all Persons in the scoring roster. This selects outputs, not data.
#' @param level Conditional equal-tail interval probability; default 0.95.
#' @param quad_points Other-Person quadrature order, 7 to 241; default uses the
#'   fit's order. Scores are checked at `2 * quad_points + 1` points.
#' @param missing `"fail"` or `"omit"`. Source replay uses the fit's policy;
#'   omission in a supplied roster must be requested explicitly.
#' @details For each requested Person, the other Persons' abilities are
#'   integrated by quadrature. At each value of the requested ability, the
#'   full shared-rater vector is integrated by a conditional Laplace
#'   approximation. Multiplying this likelihood by the normal ability density
#'   and normalizing gives an approximate continuous marginal posterior.
#'   EAP, posterior SD and equal-tail endpoints are calculated from that
#'   density. Neither independent rater marginals nor rater modes substituted
#'   as known values define this calculation.
#'
#'   Writing \eqn{L_j(u)} for Person j's likelihood integrated over that
#'   Person's normal ability population, the target for Person p is
#'   \deqn{q_p(t) \propto \phi_p(t) \int \phi_r(u)
#'     P(y_p\mid t,u) \prod_{j\ne p} L_j(u)\,du.}
#'   Calibration arguments are suppressed in this expression. The rater
#'   integral is approximated separately at each t, and the resulting density
#'   is normalized over t; a rater mode is not held fixed across t.
#'
#'   `InterpolationDifference` checks refinement of the continuous density,
#'   and `IntegrationDifference` compares other-Person quadrature orders.
#'   Small values do not bound the separate rater Laplace approximation error
#'   or establish interval coverage. Endpoints invert a continuous CDF, not the steps of
#'   a quadrature-grid CDF. Computation grows with roster size and requested
#'   Persons. Start with a few IDs when checking a large roster; `persons = NULL`
#'   computes every Person's posterior and can be slow. Selecting outputs keeps
#'   all scoring responses. No fixed effects or population variances are re-estimated.
#'
#'   With `newdata`, reusing a source rater ID does not import that rater's
#'   saved conditional distribution. Shared-rater distributions are calculated
#'   from the supplied responses and fitted population parameters. To retain
#'   earlier responses, include them in the complete scoring table, with each
#'   intended rating event appearing once. By contrast,
#'   [predict.mfrm_random_rater()] with `rater = "observed"` uses the rater's
#'   conditional distribution from the calibration responses at abilities you
#'   supply; it does not estimate abilities from new responses.
#'
#'   Intervals condition on all supplied observed responses, fitted calibration
#'   and the fitted or specified normal ability population. They exclude
#'   calibration/population estimation uncertainty, are not fixed-ability
#'   frequentist confidence intervals and do not test Person differences.
#'   They do not qualify few-rater uncertainty or the normality and assignment
#'   assumptions. With an estimated zero rater variance they describe that
#'   fitted submodel, not certainty that rater differences are absent.
#'
#'   An entirely missing Person explicitly retained with `missing = "omit"`
#'   returns the normal prior with status `"prior_only"`, not a measured average
#'   ability. Zero fitted ability variance gives `"unavailable"` scores rather
#'   than zero-width intervals. Errors or unresolved numerical checks preserve
#'   the requested Person's row and reason. No missing scores are imputed.
#'   Older saved fits keep their original known N(0,1) ability population.
#'   If an older fit omitted scores without saving the original roster, supply
#'   that roster explicitly to retain its unobserved Persons and row accounting.
#' @section Numerical evidence:
#' An independent joint-posterior comparison at fixed calibration covered
#' three-category RSMs with six or 24 raters and rotating or weakly linked
#' assignments. Eight scoring rosters contained 240 Persons and 1,440 responses;
#' four reduced rosters contained 48 Persons and 96 responses. The 48 selected
#' EAPs and posterior SDs met a 0.05-logit numerical tolerance, including
#' Monte Carlo uncertainty. Conditional-CDF calculations using the saved joint
#' rater draws supported the stated 0.10-logit tolerance for all 96 interval
#' endpoints, allowing for Monte Carlo and numerical error.
#'
#' These are bounded numerical comparisons, not universal accuracy bounds or
#' performance cutoffs. They do not establish repeated-sampling 95 percent
#' coverage when calibration is estimated, qualify the calibration likelihood,
#' or test Person differences. See `vignette("mfrmr-random-raters")` for the
#' reference-precision limitations and interpretation.
#' @return An `mfrm_random_rater_scores` object with `table`, complete scoring
#'   roster, data accounting, settings and source metadata for [mfrm_results()].
#'   Save it with `saveRDS()`. Its summary and plots need no live optimizer.
#'   Keep `table$Status` and `table$Reason` with exported scores: `prior_only`
#'   describes the assumed population without observed responses for that Person;
#'   `unavailable` withholds the estimate, SD and bounds with a recorded cause.
#' @param x For `score_mfrm_random_rater()`: A scoring result.
#'
#' For `plot()`: A result from [score_mfrm_random_rater()].
#' @param ... Unused.
#' @seealso [predict.mfrm_random_rater()], [plot.mfrm_random_rater_scores()]
#' @examples example <- readRDS(system.file("examples", "extended-models.rds", package = "mfrmr"))
#' fit <- example$random_rater$fit
#' scores <- example$random_rater$scores
#' # To recompute (requires RTMB >= 2.0 and can take several minutes):
#' # scores <- score_mfrm_random_rater(fit, persons = c("P001", "P002"))
#' # The complete source roster still informs both selected Persons.
#' scores$table
#' plot(scores)
#' mfrm_results(fit, scores = scores, compute = "never")
#' @export
#' @order 1
score_mfrm_random_rater <- function(object, newdata = NULL, persons = NULL,
    level = .95, quad_points = object$settings$quad_points,
    missing = if (is.null(newdata)) object$settings$missing else "fail") {
  if (!inherits(object, "mfrm_random_rater") || !isTRUE(object$checks$NumericalReady) ||
      !isTRUE(object$checks$InformationPositive)) stop("Resolve the fit's numerical and information checks before scoring.", call. = FALSE)
  mfrm_random_rater_interval_level(level)
  if (!is.numeric(quad_points) || is.complex(quad_points) || length(quad_points) != 1L ||
      !is.finite(quad_points) || quad_points != floor(quad_points) || quad_points < 7 || quad_points > 241) {
    stop("Use integer quad_points from 7 to 241.", call. = FALSE)
  }
  missing <- match.arg(missing, c("fail", "omit")); columns <- object$input$columns
  if (is.null(newdata) && is.null(object$input$assigned_data) && length(object$input$omitted_rows)) {
    stop("This saved fit lacks its original missing-score roster; supply the complete scoring roster as newdata.", call. = FALSE)
  }
  data <- if (is.null(newdata)) object$input$assigned_data %||% object$input$data else newdata
  input <- mfrm_random_rater_data(data, columns$person, columns$rater, columns$facets,
    columns$score, object$input$score_levels, missing, reference = object$input)
  ids <- input$levels[[columns$person]]
  if (is.null(persons)) persons <- ids
  if (!is.character(persons) || !length(persons) || anyNA(persons) || anyDuplicated(persons) ||
      !all(persons %in% ids)) stop("`persons` must name distinct Person IDs in the complete scoring roster.", call. = FALSE)
  ps <- object$calibration$person_sd %||% 1
  result <- lapply(persons, function(id) {
    person <- match(id, ids); observed <- sum(input$person == person)
    reason <- ""; difference <- interpolation <- NA_real_
    value <- setNames(rep(NA_real_, 4), c("Estimate", "ConditionalSD", "Lower", "Upper"))
    if (ps == 0) {
      status <- "unavailable"; reason <- "The fitted ability variance is zero; individual ability scores are withheld."
    } else if (observed == 0) {
      alpha <- (1 - level) / 2
      value <- c(Estimate = 0, ConditionalSD = ps, Lower = ps * stats::qnorm(alpha),
        Upper = ps * stats::qnorm(alpha, lower.tail = FALSE))
      status <- "prior_only"; difference <- interpolation <- 0
      reason <- "No observed responses for this Person; the estimate, SD and interval describe only the assumed population, not measured individual ability."
    } else {
      value <- tryCatch(withCallingHandlers({
        low <- mfrm_random_rater_person_interval(input, object$calibration, person, quad_points, level)
        high <- mfrm_random_rater_person_interval(input, object$calibration, person, 2 * quad_points + 1, level)
        delta <- abs(high$value - low$value); difference <- max(delta)
        interpolation <- max(low$difference, high$difference)
        if (any(!is.finite(delta)) || any(delta > c(1e-7, 1e-7, 1e-5, 1e-5, 1e-6))) {
          stop("Other-Person quadrature is unresolved; increase quad_points and review the scoring roster.")
        }
        high$value[1:4]
      }, warning = function(w) stop(conditionMessage(w), call. = FALSE)),
      error = function(e) { reason <<- conditionMessage(e); NULL })
      status <- if (is.null(value)) "unavailable" else "available_conditional"
      if (is.null(value)) value <- setNames(rep(NA_real_, 4), c("Estimate", "ConditionalSD", "Lower", "Upper"))
    }
    data.frame(Person = id, Observed = observed,
      Raters = length(unique(input$rater[input$person == person])), as.list(value), Status = status,
      IntegrationDifference = difference, InterpolationDifference = interpolation, Reason = reason, row.names = NULL)
  })
  structure(list(table = do.call(rbind, result), scoring_data = input$assigned_data,
    data_usage = c(Input = input$input_rows, Observed = length(input$y), Omitted = length(input$omitted_rows),
      RosterPersons = length(ids), RequestedPersons = length(persons)),
    omitted_rows = input$omitted_rows, calibration = object$calibration,
    source = mfrm_extended_prediction_source(object),
    settings = list(level = level, quad_points = quad_points, check_points = 2 * quad_points + 1,
      missing = missing, calibration_uncertainty = FALSE,
      person_variance = ps^2, rater_variance = object$calibration$rater_sd^2,
      estimated_rater_variance_boundary = isTRUE(object$checks$EstimatedVarianceBoundary),
      estimated_person_variance_boundary = isTRUE(object$checks$EstimatedPersonVarianceBoundary),
      target = "Person marginal posterior conditional on calibration and the complete scoring roster",
      rater_integration = "Joint conditional Laplace; approximation accuracy is not certified by numerical checks",
      roster = "Supplied table replaces source; output selection retains all observed responses")),
    class = "mfrm_random_rater_scores")
}

mfrm_random_rater_person_interval <- function(input, calibration, person, quad_points, level) {
  ps <- calibration$person_sd %||% 1
  obj <- mfrm_random_rater_objective(input, quad_points, fixed_sd = calibration$rater_sd,
    fixed_person_sd = ps, focal_person = person)
  par <- c(calibration$beta, calibration$steps)
  cache <- new.env(parent = emptyenv())
  loglikelihood <- function(z) vapply(z, function(x) {
    key <- sprintf("%a", x)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      value <- -as.numeric(obj$fn(c(par, ps * x)))
      if (!is.finite(value)) stop("Conditional shared-rater integration is nonfinite.")
      assign(key, value, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }, numeric(1))
  offset <- loglikelihood(0)
  logdensity <- function(z) (loglikelihood(z) - offset) + stats::dnorm(z, log = TRUE)
  mode <- stats::optimize(function(z) -logdensity(z), c(-20, 20), tol = 1e-8)$minimum
  if (abs(mode) > 19.9) stop("The conditional posterior mode is outside the scoring search range.")
  h <- .01
  curvature <- -(logdensity(mode + h) - 2 * logdensity(mode) + logdensity(mode - h)) / h^2
  span <- if (is.finite(curvature) && curvature > 0) min(12, max(1, 8 / sqrt(curvature))) else 12
  previous <- NULL
  # Reuse nested density knots; true likelihood evaluations are retained in
  # the tails. Refinement checks the continuous integrals, not grid quantiles.
  for (n in c(65L, 129L, 257L, 513L, 1025L)) {
    knots <- seq(mode - span, mode + span, length.out = n)
    spline <- stats::splinefun(knots, loglikelihood(knots) - offset, method = "natural")
    interpolated <- function(z) {
      inside <- z >= knots[1] & z <= tail(knots, 1)
      value <- numeric(length(z)); value[inside] <- spline(z[inside])
      value[!inside] <- loglikelihood(z[!inside]) - offset; value
    }
    value <- mfrm_person_posterior_interval(interpolated, ps, level)
    if (!is.null(previous)) {
      delta <- abs(value - previous)
      if (all(is.finite(delta)) && all(delta < c(1e-8, 1e-8, 1e-6, 1e-6, 1e-7))) {
        value["LogNormalizer"] <- value["LogNormalizer"] + offset
        return(list(value = value, difference = max(delta), density_points = n))
      }
    }
    previous <- value
  }
  stop("Continuous posterior density interpolation is unresolved; scores are unavailable.")
}

#' @rdname score_mfrm_random_rater
#' @name summary.mfrm_random_rater_scores
#' @aliases summary.mfrm_random_rater_scores
#' @usage \method{summary}{mfrm_random_rater_scores}(object, ...)
#' @rawNamespace S3method(summary,mfrm_random_rater_scores,summary.mfrm_testlet_scores)
NULL

#' @rdname score_mfrm_random_rater
#' @export
print.mfrm_random_rater_scores <- function(x, ...) {
  cat("Conditional Person scores with jointly integrated shared raters\n")
  print(x$table, row.names = FALSE)
  cat("Intervals exclude calibration-estimation uncertainty and use a rater Laplace approximation.\n")
  invisible(x)
}

#' @rdname score_mfrm_random_rater
#' @details `plot()` method: Plot Person scores from a shared-rater model
#' @inheritParams plot.mfrm_testlet_scores
#' @seealso [plot.mfrm_testlet_scores()] for view selection, accessible output
#'   and the meaning of the empirical distribution.
#' @param draw `FALSE` returns plot data without opening a device.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object; use [plot_data()].
#' @details For `plot()`: In the default interval view, all requested Persons remain labeled. Open circles indicate
#'   prior-only results; unavailable scores have no point or interval.
#'   Intervals are conditional on calibration and use joint rater Laplace
#'   integration. They are not tests of Person differences. [as_ggplot()]
#'   preserves these meanings without rerunning scoring.
#' @export
plot.mfrm_random_rater_scores <- function(x, draw = TRUE,
    style = c("interval", "precision", "distribution"),
    sort = c("input", "estimate", "uncertainty"), palette = c("accessible", "mono"),
    title = NULL, caption = NULL, show_title = FALSE, show_notes = FALSE,
    show_labels = TRUE, reference = 0, text_scale = 1, point_size = 2.5, ...) {
  rlang::check_dots_empty()
  mfrm_extended_estimate_plot(x$table, x$table$Person, "random_rater_scores", "Conditional Person scores",
    "Person ability (logits)", sprintf("%.0f%% conditional intervals | calibration fixed | rater Laplace | open: prior only", 100 * x$settings$level),
    draw, x$settings, style = style, sort = sort, palette = palette, title = title, caption = caption,
    show_title = show_title, show_notes = show_notes, show_labels = show_labels,
    reference = reference, text_scale = text_scale, point_size = point_size)
}

#' Evaluate screening outcomes against a planned simulation roster
#'
#' In a simulation where the true problem status is known, summarize how often
#' a warning rule detects affected raters or flags unaffected raters. Supply the
#' planned trials and observed flags; this function does not generate or fit
#' simulated ratings. Unknown truth in real ratings cannot supply these rates.
#' Unavailable screens remain unavailable, and correlated raters are not counted
#' as independent simulation trials.
#'
#' @param roster Data frame with `Condition`, `Replicate`, `Target` identifiers
#'   and logical `Affected`: `TRUE` denotes the prespecified departure that the
#'   screen is intended to detect. Include every planned target and replication,
#'   including failed runs. Within a condition, use the same targets and truth
#'   labels in every replication. Put different designs or methods in separate
#'   conditions. Identifiers are matched exactly as character labels.
#' @param results Data frame with the same three identifier columns and logical
#'   `Flag` (`TRUE`, `FALSE` or `NA`). Omitted result rows remain unavailable.
#'   Extra or duplicate keys are refused. Additional columns are retained in
#'   the source results, for example numerical status and failure messages.
#' @param rule Nonempty description of the prespecified screen, including its
#'   thresholds, comparison family and any selection/refitting procedure.
#' @param level Confidence level for exact binomial Monte Carlo intervals;
#'   default 0.95. These describe simulation uncertainty, not an interval for
#'   rater severity or a guarantee of screening accuracy.
#'
#' @return An `mfrm_screening_performance` object with `by_target`, `by_family`,
#'   aligned `outcomes`, replication-level `family_outcomes`, the original
#'   `roster` and `results`, and `settings`. Tables retain full precision.
#'   `Planned`, `Available`, `Unavailable` and `Positive` give denominators and
#'   counts. `Rate` and `MCSE` condition on available outcomes; `MCLower` and
#'   `MCUpper` are exact binomial bounds. `AllTrialsLower` and `AllTrialsUpper`
#'   bound the realized proportion across all planned trials by assigning every
#'   unknown outcome negative or positive. These latter bounds are not confidence
#'   intervals. A family with no relevant targets has zero planned trials and
#'   undefined rates. `CompleteScreens` counts fully observed families.
#'
#' @details A target's rate is sensitivity when `Affected = TRUE` and a false
#'   flag rate when `Affected = FALSE`. The family summaries count whether any
#'   unaffected target was flagged, and whether any affected target was flagged,
#'   once per replication. They do not pool all raters into a binomial sample.
#'   An observed positive makes an any-flag event known even if another target
#'   is unavailable. Without a positive, any unavailable member leaves the event
#'   unknown; only a complete negative family establishes no flag.
#'
#'   Replications must be independent and identically generated within each
#'   condition for the binomial uncertainty calculation. Pairing conditions does
#'   not invalidate their separate intervals, but a comparison between paired
#'   conditions needs its own paired analysis. Zero observed false flags does
#'   not imply a zero population probability: the exact upper bound stays positive.
#'   Conditional rates may be biased for the intended all-trial rate if failures
#'   depend on outcomes. Inspect availability and all-trial bounds together.
#'
#'   This helper neither creates flags nor establishes their validity. A severity
#'   difference alone is not rater misfit. Known simulation truth is required;
#'   labeling unknown real-world rater quality does not create a validation study.
#'   Define computable screening output separately from formal inference eligibility.
#'   The recorded rule cannot prove that it was selected before seeing outcomes.
#'
#' @references Morris, T. P., White, I. R. and Crowther, M. J. (2019).
#'   Using simulation studies to evaluate statistical methods.
#'   *Statistics in Medicine*, 38, 2074--2102. \doi{10.1002/sim.8086}.
#' @seealso [mfrm_screening_sensitivity()] for several threshold bands,
#'   [evaluate_mfrm_signal_detection()], [facet_quality_dashboard()]
#' @examples
#' roster <- expand.grid(Condition = "Null", Replicate = 1:4,
#'                       Target = c("R1", "R2"), stringsAsFactors = FALSE)
#' roster$Affected <- FALSE
#' results <- roster[c("Condition", "Replicate", "Target")]
#' results$Flag <- c(FALSE, TRUE, FALSE, NA, FALSE, NA, FALSE, NA)
#' performance <- mfrm_screening_performance(roster, results,
#'   rule = "Infit or Outfit outside [0.5, 1.5]; both raters; no exclusion/refit")
#' summary(performance)
#' plot(performance)
#' @export
#' @order 1
mfrm_screening_performance <- function(roster, results, rule, level = .95) {
  ids <- c("Condition", "Replicate", "Target")
  validate <- function(data, outcome, name) {
    if (!is.data.frame(data) || anyDuplicated(names(data)) ||
        !all(c(ids, outcome) %in% names(data)) ||
        !is.logical(data[[outcome]]) || !is.null(dim(data[[outcome]])) ||
        (outcome == "Affected" && anyNA(data[[outcome]]))) {
      stop("`", name, "` needs Condition, Replicate, Target and logical ", outcome,
        if (outcome == "Affected") " without missing truth labels." else ".", call. = FALSE)
    }
    for (id in ids) {
      value <- data[[id]]
      if (!is.atomic(value) || !is.null(dim(value)) || anyNA(value) ||
          any(!nzchar(trimws(as.character(value))))) {
        stop("Use complete nonempty identifiers in `", name, "`.", call. = FALSE)
      }
      data[[id]] <- as.character(value)
    }
    if (anyDuplicated(data[ids])) stop("Duplicate condition/replication/target in `", name, "`.", call. = FALSE)
    data[c(ids, outcome)]
  }
  planned <- validate(roster, "Affected", "roster")
  observed <- validate(results, "Flag", "results")
  if (!nrow(planned)) stop("`roster` must contain the planned trials.", call. = FALSE)
  if (!is.character(rule) || length(rule) != 1L || is.na(rule) || !nzchar(trimws(rule))) {
    stop("Describe the prespecified screening rule in `rule`.", call. = FALSE)
  }
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) stop("Supply 0 < level < 1.", call. = FALSE)
  for (condition in unique(planned$Condition)) {
    d <- planned[planned$Condition == condition, ]
    if (nrow(d) != length(unique(d$Replicate)) * length(unique(d$Target)) ||
        any(vapply(split(d$Affected, d$Target), function(a) length(unique(a)) != 1L, logical(1)))) {
      stop("Each condition needs the same planned targets and truth labels in every replication.", call. = FALSE)
    }
  }
  levels <- lapply(ids, function(id) unique(c(planned[[id]], observed[[id]])))
  key <- function(data) do.call(paste, Map(match, data[ids], levels))
  matched <- match(key(observed), key(planned))
  if (anyNA(matched)) stop("Results contain targets or trials outside the planned roster.", call. = FALSE)
  aligned <- planned
  aligned$Flag <- rep(NA, nrow(planned)); aligned$Reported <- FALSE
  aligned$Flag[matched] <- observed$Flag
  aligned$Reported[matched] <- TRUE
  rownames(aligned) <- NULL
  target_rows <- family_rows <- event_rows <- list()
  for (condition in unique(aligned$Condition)) {
    d <- aligned[aligned$Condition == condition, ]
    for (target in unique(d$Target)) {
      a <- d[d$Target == target, ]
      target_rows[[length(target_rows) + 1L]] <- cbind(data.frame(
        Condition = condition, Target = target, Affected = a$Affected[1],
        Metric = if (a$Affected[1]) "Sensitivity" else "False flag rate"),
        mfrm_screening_rate(a$Flag, level))
    }
    for (affected in c(FALSE, TRUE)) {
      a <- d[d$Affected == affected, ]
      metric <- if (affected) "Any affected target flagged" else "Any unaffected target flagged"
      events <- lapply(unique(a$Replicate), function(rep) {
        flags <- a$Flag[a$Replicate == rep]
        data.frame(Condition = condition, Replicate = rep, Metric = metric,
          Flag = if (any(flags %in% TRUE)) TRUE else if (anyNA(flags)) NA else FALSE,
          Complete = !anyNA(flags))
      })
      events <- if (length(events)) do.call(rbind, events) else data.frame(
        Condition = character(), Replicate = character(), Metric = character(),
        Flag = logical(), Complete = logical())
      event_rows[[length(event_rows) + 1L]] <- events
      family_rows[[length(family_rows) + 1L]] <- cbind(data.frame(
        Condition = condition, Metric = metric, Targets = length(unique(a$Target)),
        CompleteScreens = sum(events$Complete)), mfrm_screening_rate(events$Flag, level))
    }
  }
  out <- list(by_target = do.call(rbind, target_rows), by_family = do.call(rbind, family_rows),
    outcomes = aligned, family_outcomes = do.call(rbind, event_rows), roster = roster,
    results = results, settings = list(rule = rule, level = level,
      uncertainty_unit = "Independent simulation replication within condition"))
  class(out) <- "mfrm_screening_performance"
  out
}

mfrm_screening_rate <- function(flag, level = .95) {
  planned <- length(flag)
  available <- sum(!is.na(flag))
  positive <- sum(flag %in% TRUE)
  rate <- if (available) positive / available else NA_real_
  bounds <- if (available) stats::binom.test(positive, available, conf.level = level)$conf.int else c(NA_real_, NA_real_)
  data.frame(Planned = planned, Available = available, Unavailable = planned - available,
    Positive = positive, Rate = rate, MCSE = if (available) sqrt(rate * (1 - rate) / available) else NA_real_,
    MCLower = bounds[1], MCUpper = bounds[2],
    AllTrialsLower = if (planned) positive / planned else NA_real_,
    AllTrialsUpper = if (planned) (positive + planned - available) / planned else NA_real_)
}

#' @rdname mfrm_screening_performance
#' @param object An object returned by `mfrm_screening_performance()`.
#' @param ... For `summary.mfrm_screening_performance()`: Unused by print and summary methods.
#'
#' For `plot()`: Unused.
#' @export
summary.mfrm_screening_performance <- function(object, ...) {
  list(by_target = object$by_target, by_family = object$by_family, settings = object$settings)
}

#' @rdname mfrm_screening_performance
#' @export
print.mfrm_screening_performance <- function(x, ...) {
  cat("Screening performance against planned simulation targets\n")
  cat("Rule:", x$settings$rule, "\n")
  print(x$by_family, row.names = FALSE)
  cat("Rates condition on known outcomes; all-trial bounds retain unresolved outcomes.\n")
  invisible(x)
}

#' @rdname mfrm_screening_performance
#' @details `plot()` method: Plot screening performance and unresolved outcomes
#' @param x Output from [mfrm_screening_performance()].
#'
#' Also for `summary()`: An object returned by `mfrm_screening_performance()`.
#' @param scope `"family"` (default) shows any-flag events once per replication;
#'   `"target"` shows each individual target separately.
#' @param metric `"false_flags"` (default) or `"detection"`.
#' @param condition Optional condition labels to include, matched exactly.
#' @param draw Logical; `FALSE` returns exact plot data without opening a device.
#' @return For `plot()`: Invisibly, an `mfrm_plot_data` object. Use [plot_data()] for custom
#'   graphics; automatic [as_ggplot()] conversion is unavailable.
#' @details For `plot()`: Blue points and whiskers show conditional rates and their exact
#'   binomial Monte Carlo intervals. Gray segments bound the realized all-trial
#'   proportion by assigning unresolved outcomes negative or positive; these
#'   are not confidence intervals. A row with no known outcomes has no point
#'   estimate. Right-hand counts show available/planned events. No model is
#'   fitted and no threshold is selected by this plot.
#' @export
plot.mfrm_screening_performance <- function(x, scope = c("family", "target"),
                                            metric = c("false_flags", "detection"),
                                            condition = NULL, draw = TRUE, ...) {
  rlang::check_dots_empty()
  scope <- match.arg(scope); metric <- match.arg(metric)
  if (!is.logical(draw) || length(draw) != 1L || is.na(draw)) stop("`draw` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(condition) && (!is.character(condition) || anyNA(condition) ||
      !length(condition) || any(!condition %in% x$by_family$Condition))) {
    stop("Use recorded condition labels in `condition`.", call. = FALSE)
  }
  tab <- if (scope == "family") x$by_family else x$by_target
  label <- if (metric == "false_flags") "Any unaffected target flagged" else "Any affected target flagged"
  keep <- if (scope == "family") tab$Metric == label else tab$Affected == (metric == "detection")
  if (!is.null(condition)) keep <- keep & tab$Condition %in% condition
  tab <- tab[keep & tab$Planned > 0, , drop = FALSE]
  if (!nrow(tab)) stop("No planned targets for the selected metric and condition.", call. = FALSE)
  labels <- if (scope == "family") tab$Condition else paste(tab$Condition, tab$Target, sep = " / ")
  title <- if (scope == "family") label else if (metric == "false_flags") "False flags by target" else "Detection by target"
  out <- new_mfrm_plot_data("screening_performance", list(table = tab, labels = labels,
    title = title, scope = scope, metric = metric, settings = x$settings,
    subtitle = sprintf("Conditional rates and %.0f%% Monte Carlo intervals", 100 * x$settings$level),
    caption = "Gray segments: bounds from unresolved outcomes (not confidence intervals)"))
  if (!draw) return(invisible(out))
  old <- graphics::par(mar = c(7, 12, 5, 5))
  on.exit(graphics::par(old), add = TRUE)
  y <- rev(seq_len(nrow(tab)))
  graphics::plot(NA_real_, NA_real_, xlim = c(0, 1), ylim = c(.5, nrow(tab) + .5),
    xlab = "Proportion", ylab = "", yaxt = "n", main = "")
  graphics::segments(tab$AllTrialsLower, y, tab$AllTrialsUpper, y, col = "grey75", lwd = 6)
  graphics::segments(tab$MCLower, y, tab$MCUpper, y, col = "#2274A5", lwd = 2)
  graphics::points(tab$Rate, y, pch = 16, col = "#2274A5")
  axis_labels <- vapply(labels, function(label) paste(strwrap(label, width = 32), collapse = "\n"), character(1))
  graphics::axis(2, at = y, labels = axis_labels, las = 1, cex.axis = .8)
  graphics::axis(4, at = y, labels = paste0(tab$Available, "/", tab$Planned), las = 1, tick = FALSE, cex.axis = .8)
  graphics::text(1.08, nrow(tab) + .8, "Known /\nplanned", xpd = NA, cex = .7)
  invisible(out)
}

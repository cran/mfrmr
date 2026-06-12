# Response-time diagnostics are intentionally descriptive. They do not alter
# the ordered-response likelihood fitted by fit_mfrm().

.response_time_check_column <- function(data, column, role) {
  if (!is.character(column) || length(column) != 1L || !nzchar(column)) {
    stop("`", role, "` must be a single column name.", call. = FALSE)
  }
  if (!column %in% names(data)) {
    stop("`", role, "` column '", column, "' was not found in `data`.",
         call. = FALSE)
  }
  column
}

.response_time_check_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop("`", name, "` must be a single number in (0, 1).",
         call. = FALSE)
  }
  as.numeric(x)
}

.response_time_empty_summary <- function(group_cols = character()) {
  out <- as.data.frame(
    setNames(replicate(length(group_cols), character(0), simplify = FALSE),
             group_cols),
    stringsAsFactors = FALSE
  )
  out$N <- integer(0)
  out$MeanTime <- numeric(0)
  out$MedianTime <- numeric(0)
  out$SDTime <- numeric(0)
  out$MeanLogTime <- numeric(0)
  out$MedianLogTime <- numeric(0)
  out$RapidResponses <- integer(0)
  out$RapidRate <- numeric(0)
  out$SlowResponses <- integer(0)
  out$SlowRate <- numeric(0)
  out
}

.response_time_group_summary <- function(obs, group_cols) {
  if (!is.data.frame(obs) || nrow(obs) == 0L) {
    return(.response_time_empty_summary(group_cols))
  }
  keep <- stats::complete.cases(obs[, group_cols, drop = FALSE])
  obs <- obs[keep, , drop = FALSE]
  if (nrow(obs) == 0L) {
    return(.response_time_empty_summary(group_cols))
  }
  key <- interaction(obs[, group_cols, drop = FALSE], drop = TRUE, sep = "\r")
  idx <- split(seq_len(nrow(obs)), key)
  rows <- lapply(idx, function(i) {
    grp <- obs[i[1L], group_cols, drop = FALSE]
    t <- obs$Time[i]
    lt <- obs$LogTime[i]
    data.frame(
      grp,
      N = length(i),
      MeanTime = mean(t, na.rm = TRUE),
      MedianTime = stats::median(t, na.rm = TRUE),
      SDTime = if (length(i) > 1L) stats::sd(t, na.rm = TRUE) else NA_real_,
      MeanLogTime = mean(lt, na.rm = TRUE),
      MedianLogTime = stats::median(lt, na.rm = TRUE),
      RapidResponses = sum(obs$RapidFlag[i], na.rm = TRUE),
      RapidRate = mean(obs$RapidFlag[i], na.rm = TRUE),
      SlowResponses = sum(obs$SlowFlag[i], na.rm = TRUE),
      SlowRate = mean(obs$SlowFlag[i], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.response_time_flag_rows <- function(summary_tbl, source, label_cols,
                                     min_n_flag, rapid_rate_warn,
                                     slow_rate_warn) {
  empty <- data.frame(
    Source = character(0), Group = character(0), Flag = character(0),
    Rate = numeric(0), N = integer(0), ThresholdRate = numeric(0),
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(summary_tbl) || nrow(summary_tbl) == 0L) return(empty)
  label <- apply(summary_tbl[, label_cols, drop = FALSE], 1L, paste,
                 collapse = " / ")
  rapid_idx <- which(summary_tbl$N >= min_n_flag &
                       is.finite(summary_tbl$RapidRate) &
                       summary_tbl$RapidRate >= rapid_rate_warn)
  slow_idx <- which(summary_tbl$N >= min_n_flag &
                     is.finite(summary_tbl$SlowRate) &
                     summary_tbl$SlowRate >= slow_rate_warn)
  out <- list()
  if (length(rapid_idx) > 0L) {
    out[[length(out) + 1L]] <- data.frame(
      Source = source,
      Group = label[rapid_idx],
      Flag = "high_rapid_response_rate",
      Rate = summary_tbl$RapidRate[rapid_idx],
      N = summary_tbl$N[rapid_idx],
      ThresholdRate = rapid_rate_warn,
      stringsAsFactors = FALSE
    )
  }
  if (length(slow_idx) > 0L) {
    out[[length(out) + 1L]] <- data.frame(
      Source = source,
      Group = label[slow_idx],
      Flag = "high_slow_response_rate",
      Rate = summary_tbl$SlowRate[slow_idx],
      N = summary_tbl$N[slow_idx],
      ThresholdRate = slow_rate_warn,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0L) empty else do.call(rbind, out)
}

#' Review response-time patterns outside the MFRM likelihood
#'
#' @description
#' Build a descriptive response-time review table from the same long-format
#' rating-event data used by [fit_mfrm()]. This helper does not fit a joint
#' response-time model and does not change MFRM estimates. It summarizes
#' response-time distributions, distributional rapid/slow flags, and person /
#' facet / score-level response-time patterns for screening and reporting
#' context.
#'
#' @param data A data.frame in long format with one row per observed rating
#'   event.
#' @param person Column name for the person identifier.
#' @param facets Optional character vector of facet columns to summarize.
#' @param time Column name containing positive response times.
#' @param score Optional ordered-score column. When supplied, score-level
#'   response-time summaries are returned.
#' @param time_unit Label for the response-time unit, such as `"seconds"`.
#' @param min_time Minimum valid response time. Values must be strictly greater
#'   than this threshold; default 0.
#' @param rapid_threshold Optional numeric response-time cutoff for rapid
#'   responses. When `NULL`, it is estimated from `rapid_quantile`.
#' @param slow_threshold Optional numeric response-time cutoff for slow
#'   responses. When `NULL`, it is estimated from `slow_quantile`.
#' @param rapid_quantile Quantile used when `rapid_threshold = NULL`.
#' @param slow_quantile Quantile used when `slow_threshold = NULL`.
#' @param rapid_rate_warn Group-level rapid-response rate that creates a
#'   descriptive flag; default 0.25.
#' @param slow_rate_warn Group-level slow-response rate that creates a
#'   descriptive flag; default 0.25.
#' @param min_n_flag Minimum group size before rapid/slow rates are flagged;
#'   default 3.
#'
#' @return An object of class `mfrm_response_time_review`, a list with
#'   `overview`, `thresholds`, `observations`, `person_summary`,
#'   `facet_summary`, `score_summary`, `flags`, and `notes`.
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' toy$ResponseTime <- 12 + as.numeric(factor(toy$Person)) * 0.4 +
#'   as.numeric(toy$Score)
#' rt <- response_time_review(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   time = "ResponseTime"
#' )
#' summary(rt)
#' plot_response_time_review(rt, draw = FALSE)
#' @seealso [plot_response_time_review()], [fit_mfrm()],
#'   [mfrmr_visual_diagnostics], [mfrmr_output_guide()]
#' @concept response time
#' @concept visual diagnostics
#' @concept quality control
#' @export
response_time_review <- function(data,
                                 person,
                                 facets = NULL,
                                 time,
                                 score = NULL,
                                 time_unit = "seconds",
                                 min_time = 0,
                                 rapid_threshold = NULL,
                                 slow_threshold = NULL,
                                 rapid_quantile = 0.05,
                                 slow_quantile = 0.95,
                                 rapid_rate_warn = 0.25,
                                 slow_rate_warn = 0.25,
                                 min_n_flag = 3L) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("`data` must be a non-empty data.frame.", call. = FALSE)
  }
  person <- .response_time_check_column(data, person, "person")
  time <- .response_time_check_column(data, time, "time")
  if (!is.null(score)) {
    score <- .response_time_check_column(data, score, "score")
  }
  if (!is.null(facets)) {
    if (!is.character(facets) || length(facets) == 0L ||
        any(!nzchar(facets))) {
      stop("`facets` must be a character vector of column names.",
           call. = FALSE)
    }
    missing_facets <- setdiff(facets, names(data))
    if (length(missing_facets) > 0L) {
      stop("`facets` column(s) not found in `data`: ",
           paste(missing_facets, collapse = ", "), ".", call. = FALSE)
    }
  } else {
    facets <- character(0)
  }
  if (!is.numeric(min_time) || length(min_time) != 1L ||
      !is.finite(min_time) || min_time < 0) {
    stop("`min_time` must be a single non-negative number.", call. = FALSE)
  }
  rapid_quantile <- .response_time_check_probability(
    rapid_quantile, "rapid_quantile"
  )
  slow_quantile <- .response_time_check_probability(
    slow_quantile, "slow_quantile"
  )
  if (rapid_quantile >= slow_quantile) {
    stop("`rapid_quantile` must be smaller than `slow_quantile`.",
         call. = FALSE)
  }
  rapid_rate_warn <- .response_time_check_probability(
    rapid_rate_warn, "rapid_rate_warn"
  )
  slow_rate_warn <- .response_time_check_probability(
    slow_rate_warn, "slow_rate_warn"
  )
  min_n_flag <- max(1L, as.integer(min_n_flag[1L]))

  time_values <- suppressWarnings(as.numeric(data[[time]]))
  person_values <- as.character(data[[person]])
  valid <- is.finite(time_values) & time_values > min_time &
    !is.na(person_values) & nzchar(person_values)
  if (!any(valid)) {
    stop("No rows contain finite response times greater than `min_time` ",
         "and non-missing person IDs.", call. = FALSE)
  }
  valid_times <- time_values[valid]

  if (is.null(rapid_threshold)) {
    rapid_threshold <- as.numeric(stats::quantile(
      valid_times, probs = rapid_quantile, na.rm = TRUE, names = FALSE,
      type = 7
    ))
  } else {
    rapid_threshold <- as.numeric(rapid_threshold[1L])
  }
  if (is.null(slow_threshold)) {
    slow_threshold <- as.numeric(stats::quantile(
      valid_times, probs = slow_quantile, na.rm = TRUE, names = FALSE,
      type = 7
    ))
  } else {
    slow_threshold <- as.numeric(slow_threshold[1L])
  }
  if (!is.finite(rapid_threshold) || !is.finite(slow_threshold) ||
      rapid_threshold <= min_time || slow_threshold <= min_time) {
    stop("Rapid and slow response-time thresholds must be finite and ",
         "greater than `min_time`.", call. = FALSE)
  }

  obs <- data.frame(
    Row = which(valid),
    Person = person_values[valid],
    Time = valid_times,
    LogTime = log(valid_times),
    RapidFlag = valid_times <= rapid_threshold,
    SlowFlag = valid_times >= slow_threshold,
    stringsAsFactors = FALSE
  )
  if (!is.null(score)) {
    obs$Score <- as.character(data[[score]][valid])
  }

  person_summary <- .response_time_group_summary(obs, "Person")

  facet_summary <- lapply(facets, function(f) {
    facet_obs <- obs
    facet_obs$Facet <- f
    facet_obs$Level <- as.character(data[[f]][valid])
    facet_obs <- facet_obs[!is.na(facet_obs$Level) & nzchar(facet_obs$Level),
                           , drop = FALSE]
    .response_time_group_summary(facet_obs, c("Facet", "Level"))
  })
  facet_summary <- if (length(facet_summary) > 0L) {
    do.call(rbind, facet_summary)
  } else {
    .response_time_empty_summary(c("Facet", "Level"))
  }
  rownames(facet_summary) <- NULL

  score_summary <- if (!is.null(score)) {
    score_obs <- obs[!is.na(obs$Score) & nzchar(obs$Score), , drop = FALSE]
    .response_time_group_summary(score_obs, "Score")
  } else {
    .response_time_empty_summary("Score")
  }

  flags <- do.call(rbind, list(
    .response_time_flag_rows(person_summary, "person", "Person",
                             min_n_flag, rapid_rate_warn, slow_rate_warn),
    .response_time_flag_rows(facet_summary, "facet", c("Facet", "Level"),
                             min_n_flag, rapid_rate_warn, slow_rate_warn),
    .response_time_flag_rows(score_summary, "score", "Score",
                             min_n_flag, rapid_rate_warn, slow_rate_warn)
  ))
  rownames(flags) <- NULL

  notes <- c(
    "Response-time review is descriptive; it does not change fit_mfrm estimates.",
    if (is.null(score)) {
      "No score column was supplied, so score-level response-time summaries are empty."
    } else {
      "Score-level summaries are descriptive and should not be read as response-time model parameters."
    },
    if (nrow(data) - nrow(obs) > 0L) {
      sprintf("%d row(s) were excluded because response time or person ID was invalid.",
              nrow(data) - nrow(obs))
    },
    if (rapid_threshold >= slow_threshold) {
      "Rapid and slow thresholds overlap because the response-time distribution is narrow; review flags cautiously."
    }
  )
  notes <- notes[nzchar(notes)]

  overview <- data.frame(
    Rows = nrow(data),
    ValidRows = nrow(obs),
    DroppedRows = nrow(data) - nrow(obs),
    Persons = length(unique(obs$Person)),
    Facets = length(facets),
    TimeColumn = time,
    ScoreColumn = as.character(score %||% ""),
    TimeUnit = as.character(time_unit[1L] %||% "seconds"),
    MedianTime = stats::median(obs$Time, na.rm = TRUE),
    MeanLogTime = mean(obs$LogTime, na.rm = TRUE),
    RapidThreshold = rapid_threshold,
    SlowThreshold = slow_threshold,
    RapidRate = mean(obs$RapidFlag, na.rm = TRUE),
    SlowRate = mean(obs$SlowFlag, na.rm = TRUE),
    FlaggedGroups = nrow(flags),
    InterpretationBoundary = paste(
      "Descriptive response-time screening; not a joint speed-accuracy model",
      "and not a fit/pass-fail rule."
    ),
    stringsAsFactors = FALSE
  )
  thresholds <- data.frame(
    Threshold = c("rapid", "slow"),
    Value = c(rapid_threshold, slow_threshold),
    Basis = c(
      if (is.null(match.call()$rapid_threshold)) {
        paste0("quantile_", rapid_quantile)
      } else {
        "user_supplied"
      },
      if (is.null(match.call()$slow_threshold)) {
        paste0("quantile_", slow_quantile)
      } else {
        "user_supplied"
      }
    ),
    TimeUnit = overview$TimeUnit,
    stringsAsFactors = FALSE
  )

  out <- list(
    overview = overview,
    thresholds = thresholds,
    observations = obs,
    person_summary = person_summary,
    facet_summary = facet_summary,
    score_summary = score_summary,
    flags = flags,
    notes = notes,
    config = list(
      person = person,
      facets = facets,
      score = score,
      time = time,
      time_unit = overview$TimeUnit,
      min_time = min_time,
      rapid_quantile = rapid_quantile,
      slow_quantile = slow_quantile,
      rapid_rate_warn = rapid_rate_warn,
      slow_rate_warn = slow_rate_warn,
      min_n_flag = min_n_flag
    )
  )
  class(out) <- c("mfrm_response_time_review", class(out))
  out
}

#' @export
print.mfrm_response_time_review <- function(x, ...) {
  cat("<mfrm_response_time_review>\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    ov <- x$overview[1L, , drop = FALSE]
    cat("  rows       : ", ov$ValidRows, " valid / ", ov$Rows, " total\n",
        sep = "")
    cat("  persons    : ", ov$Persons, "\n", sep = "")
    cat("  median time: ", signif(ov$MedianTime, 4), " ", ov$TimeUnit, "\n",
        sep = "")
    cat("  thresholds : rapid <= ", signif(ov$RapidThreshold, 4),
        ", slow >= ", signif(ov$SlowThreshold, 4), "\n", sep = "")
    cat("  flags      : ", ov$FlaggedGroups, " group(s)\n", sep = "")
  }
  cat("Use `summary(x)` for top flagged groups and ",
      "`plot_response_time_review(x)` for draw-free plot data or graphics.\n",
      sep = "")
  invisible(x)
}

#' Summarize a response-time review
#'
#' @param object An object returned by [response_time_review()].
#' @param top_n Number of top groups to retain in summary previews.
#' @param ... Unused.
#' @return A list of class `summary.mfrm_response_time_review`.
#' @export
summary.mfrm_response_time_review <- function(object, top_n = 10L, ...) {
  if (!inherits(object, "mfrm_response_time_review")) {
    stop("`object` must be output from response_time_review().",
         call. = FALSE)
  }
  top_n <- max(1L, as.integer(top_n[1L]))
  order_top <- function(tbl, rate_col) {
    if (!is.data.frame(tbl) || nrow(tbl) == 0L) return(tbl)
    ord <- order(tbl[[rate_col]], tbl$N, decreasing = TRUE, na.last = TRUE)
    tbl[utils::head(ord, top_n), , drop = FALSE]
  }
  out <- list(
    overview = object$overview,
    thresholds = object$thresholds,
    top_rapid_persons = order_top(object$person_summary, "RapidRate"),
    top_slow_persons = order_top(object$person_summary, "SlowRate"),
    top_rapid_facets = order_top(object$facet_summary, "RapidRate"),
    top_slow_facets = order_top(object$facet_summary, "SlowRate"),
    flags = if (nrow(object$flags) > top_n) {
      object$flags[seq_len(top_n), , drop = FALSE]
    } else {
      object$flags
    },
    notes = object$notes
  )
  class(out) <- c("summary.mfrm_response_time_review", class(out))
  out
}

#' @export
print.summary.mfrm_response_time_review <- function(x, ...) {
  cat("mfrmr response-time review\n\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    print(x$overview, row.names = FALSE)
  }
  if (!is.null(x$thresholds) && nrow(x$thresholds) > 0L) {
    cat("\nThresholds:\n")
    print(x$thresholds, row.names = FALSE)
  }
  if (!is.null(x$flags) && nrow(x$flags) > 0L) {
    cat("\nFlagged groups:\n")
    print(x$flags, row.names = FALSE)
  }
  if (!is.null(x$notes) && length(x$notes) > 0L) {
    cat("\nNotes:\n")
    for (note in x$notes) cat("- ", note, "\n", sep = "")
  }
  invisible(x)
}

.response_time_plot_table <- function(x, type, facet = NULL, top_n = 25L) {
  top_n <- max(1L, as.integer(top_n[1L]))
  if (identical(type, "distribution")) {
    return(x$observations)
  }
  if (identical(type, "person")) {
    tbl <- x$person_summary
    if (nrow(tbl) == 0L) stop("No person response-time rows are available.",
                              call. = FALSE)
    tbl <- tbl[order(tbl$MedianTime, decreasing = TRUE), , drop = FALSE]
    return(utils::head(tbl, top_n))
  }
  if (identical(type, "facet")) {
    tbl <- x$facet_summary
    if (nrow(tbl) == 0L) stop("No facet response-time rows are available.",
                              call. = FALSE)
    if (is.null(facet)) facet <- as.character(tbl$Facet[1L])
    tbl <- tbl[as.character(tbl$Facet) == as.character(facet[1L]), ,
               drop = FALSE]
    if (nrow(tbl) == 0L) {
      stop("No response-time rows are available for facet '", facet, "'.",
           call. = FALSE)
    }
    tbl <- tbl[order(tbl$MedianTime, decreasing = TRUE), , drop = FALSE]
    return(utils::head(tbl, top_n))
  }
  if (identical(type, "score")) {
    tbl <- x$score_summary
    if (nrow(tbl) == 0L) {
      stop("No score-level response-time rows are available. ",
           "Supply `score` to response_time_review().", call. = FALSE)
    }
    score_num <- suppressWarnings(as.numeric(tbl$Score))
    tbl <- tbl[order(ifelse(is.finite(score_num), score_num, seq_len(nrow(tbl)))),
               , drop = FALSE]
    return(tbl)
  }
  stop("Unsupported response-time plot type.", call. = FALSE)
}

#' Plot response-time review summaries
#'
#' @description
#' Draw or return reusable plot data for a [response_time_review()] object.
#' Plot types are descriptive screening views and do not represent a joint
#' response-time model.
#'
#' @param x A `mfrm_response_time_review` object.
#' @param type Plot type: `"distribution"`, `"person"`, `"facet"`, or
#'   `"score"`.
#' @param facet Optional facet name when `type = "facet"`.
#' @param top_n Maximum number of person or facet rows to plot.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics. If `FALSE`, return only an
#'   `mfrm_plot_data` object.
#' @param ... Unused.
#'
#' @return Invisibly, an `mfrm_plot_data` object containing the plot table,
#'   thresholds, overview, and interpretation notes.
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' toy$ResponseTime <- 10 + seq_len(nrow(toy)) %% 6 + as.numeric(toy$Score)
#' rt <- response_time_review(
#'   toy, person = "Person", facets = c("Rater", "Criterion"),
#'   score = "Score", time = "ResponseTime"
#' )
#' plot_response_time_review(rt, type = "distribution", draw = FALSE)
#' plot_response_time_review(rt, type = "person", draw = FALSE)
#' @seealso [response_time_review()], [plot_data_components()],
#'   [mfrmr_output_guide()]
#' @concept response time
#' @concept visual diagnostics
#' @export
plot_response_time_review <- function(x,
                                      type = c("distribution", "person",
                                               "facet", "score"),
                                      facet = NULL,
                                      top_n = 25L,
                                      preset = c("standard", "publication",
                                                 "compact", "monochrome"),
                                      draw = TRUE,
                                      ...) {
  if (!inherits(x, "mfrm_response_time_review")) {
    stop("`x` must be output from response_time_review().", call. = FALSE)
  }
  type <- match.arg(type)
  style <- resolve_plot_preset(preset)
  tbl <- .response_time_plot_table(x, type = type, facet = facet,
                                   top_n = top_n)
  thresholds <- x$thresholds
  unit <- as.character(x$overview$TimeUnit[1L] %||% "seconds")
  title <- switch(
    type,
    distribution = "Response-time distribution",
    person = "Response-time review by person",
    facet = sprintf("Response-time review by %s",
                    as.character(tbl$Facet[1L] %||% "facet")),
    score = "Response-time review by score"
  )
  subtitle <- "Descriptive screening view; not a joint speed-accuracy model"

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    if (identical(type, "distribution")) {
      graphics::hist(
        tbl$LogTime,
        breaks = "FD",
        col = style$fill_soft,
        border = style$axis,
        main = title,
        xlab = paste0("log response time (", unit, ")")
      )
      graphics::abline(v = log(thresholds$Value), lty = c(2, 3),
                       col = c(style$warn, style$accent_secondary), lwd = 1.2)
      graphics::legend(
        "topright",
        legend = c("Rapid threshold", "Slow threshold"),
        lty = c(2, 3),
        col = c(style$warn, style$accent_secondary),
        bty = "n", cex = 0.85
      )
    } else if (identical(type, "score")) {
      y <- tbl$MedianTime
      x_pos <- seq_len(nrow(tbl))
      graphics::plot(
        x_pos, y, type = "b", pch = 19, col = style$accent_primary,
        xaxt = "n", xlab = "Score", ylab = paste0("Median time (", unit, ")"),
        main = title
      )
      graphics::axis(1, at = x_pos, labels = tbl$Score)
      graphics::abline(h = thresholds$Value, lty = c(2, 3),
                       col = c(style$warn, style$accent_secondary))
    } else {
      label_col <- if (identical(type, "person")) "Person" else "Level"
      labels <- as.character(tbl[[label_col]])
      y <- seq_len(nrow(tbl))
      graphics::par(mar = c(5, max(6, min(14, max(nchar(labels)) * 0.55)),
                            3, 2))
      graphics::plot(
        tbl$MedianTime, y,
        yaxt = "n", pch = 19, col = style$accent_primary,
        xlab = paste0("Median time (", unit, ")"), ylab = "",
        main = title
      )
      graphics::axis(2, at = y, labels = labels, las = 1, cex.axis = 0.75)
      graphics::abline(v = thresholds$Value, lty = c(2, 3),
                       col = c(style$warn, style$accent_secondary))
    }
  }

  ref_values <- if (identical(type, "distribution")) {
    log(thresholds$Value)
  } else {
    thresholds$Value
  }
  legend_primary <- if (identical(type, "distribution")) {
    "Response-time histogram"
  } else {
    "Median response time"
  }

  invisible(new_mfrm_plot_data(
    "response_time_review",
    list(
      table = tbl,
      thresholds = thresholds,
      overview = x$overview,
      notes = x$notes,
      type = type,
      facet = as.character(facet %||% ""),
      top_n = if (type %in% c("person", "facet")) {
        max(1L, as.integer(top_n[1L]))
      } else {
        NA_integer_
      },
      title = title,
      subtitle = subtitle,
      preset = style$name,
      legend = new_plot_legend(
        label = c(legend_primary, "Rapid threshold", "Slow threshold"),
        role = c("location", "threshold", "threshold"),
        aesthetic = c("point", "line", "line"),
        value = c(style$accent_primary, style$warn,
                  style$accent_secondary)
      ),
      reference_lines = new_reference_lines(
        axis = rep(if (identical(type, "score")) "h" else "v", 2L),
        value = ref_values,
        label = paste(thresholds$Threshold, "threshold"),
        linetype = c("dashed", "dotted"),
        role = c("rapid_response", "slow_response")
      )
    )
  ))
}

#' @export
plot.mfrm_response_time_review <- function(x,
                                           type = c("distribution", "person",
                                                    "facet", "score"),
                                           ...) {
  plot_response_time_review(x, type = type, ...)
}

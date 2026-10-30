# ==============================================================================
# Additional visualization helpers
# ==============================================================================
#
# Each function in this file follows the established mfrmr plot conventions:
#  * accepts an mfrm_fit (or related class) plus ergonomic options;
#  * resolves a `preset = c("standard", "publication", "compact", "monochrome")` style via
#    `resolve_plot_preset()` and `apply_plot_preset()`;
#  * returns an `mfrm_plot_data` object with a stable data contract so
#    downstream pipelines can re-render or export the underlying tables.
# ==============================================================================


#' Plot RSM/PCM threshold ladders with disorder highlighting
#'
#' Renders the Rasch-Andrich threshold structure as a vertical ladder per
#' step-facet level. Each tick is a `tau_k`; lines connecting adjacent
#' thresholds are coloured to make disordered crossings (`tau_{k+1} <
#' tau_k`) visually obvious. For RSM there is one ladder; for PCM (and
#' GPCM) there is one ladder per `step_facet` level.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param highlight_disorder Logical. When `TRUE` (default), draw
#'   disordered segments with the preset's `fail` colour. The number of
#'   affected groups is returned in the notes, outside the figure.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`,
#'   or `"monochrome"`).
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object with a `data` slot containing
#'   columns `Group`, `Step`, `Threshold`, `Disordered` for each ladder
#'   row.
#'
#' @section Interpreting output:
#' Each ladder shows centered step offsets for one step-owner level, including
#' when GPCM slopes belong to a different facet. Slopes do not duplicate the
#' threshold families. Disorder describes point estimates; inspect owner-level
#' counts, curves and the rubric before considering category changes.
#' Missing steps remain unavailable and are not bridged by connecting lines.
#' No intervals or full ability-axis crossings are calculated. Save the result
#' of `plot()` to retain its model context and interpretation notes.
#' Automatic ggplot conversion of this ladder is unavailable; use base graphics
#' or [plot_data()] for custom graphics. The saved [rating_scale_table()]
#' threshold view has its own ggplot conversion.
#'
#' @seealso [category_structure_report()], [category_curves_report()],
#'   [plot.mfrm_fit()] (`type = "ccc"`).
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_threshold_ladder", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_operational")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "MML", quad_points = 7, maxit = 30)
#' p <- plot(fit, type = "thresholds", draw = FALSE)
#' head(p$data$data)
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_threshold_ladder <- function(fit,
                                  highlight_disorder = TRUE,
                                  preset = c("standard", "publication", "compact", "monochrome"),
                                  draw = TRUE) {
  stop_if_product_slopes(fit, "plot(x, type = \"thresholds\")")
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  if (inherits(fit, "mfrm_imported_fit")) {
    stop("Imported thresholds do not supply the native centered-step contract.", call. = FALSE)
  }
  steps <- calc_step_order(fit$steps, expected_steps = fit$config$n_cat - 1L,
    expected_facets = fit$config$facet_levels[[fit$config$step_facet %||% ""]])
  if (!nrow(steps)) {
    stop("Step estimates are not available for this fit.", call. = FALSE)
  }
  context <- category_model_context(fit)
  rows <- data.frame(Group = steps$StepFacet, Step = steps$Step,
    StepOwner = context$StepOwner, StepIndex = steps$StepIndex,
    Threshold = steps$Estimate, Disordered = steps$Spacing < 0,
    stringsAsFactors = FALSE)
  n_disorder_groups <- length(unique(rows$Group[rows$Disordered %in% TRUE]))
  plot_title <- "Threshold ladder"
  plot_subtitle <- if (isTRUE(highlight_disorder) && n_disorder_groups > 0L) {
    sprintf("%d ladder(s) with at least one disordered step",
            n_disorder_groups)
  } else {
    "Centered step offsets (logit)"
  }

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_mar <- graphics::par("mar")
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    graphics::par(mar = c(max(old_mar[1], 5), old_mar[2], 1.5, old_mar[4]))
    groups <- unique(rows$Group)
    n_grp <- length(groups)
    x_pos <- seq_along(groups)
    y_lim <- if (any(is.finite(rows$Threshold))) range(rows$Threshold, finite = TRUE) else c(-1, 1)
    y_lim <- y_lim + c(-0.15, 0.15) * max(diff(y_lim), 1)
    plot(
      x = numeric(0), y = numeric(0),
      xlim = c(0.5, n_grp + 0.5), ylim = y_lim,
      xaxt = "n", xlab = context$StepOwner, ylab = "Step offset (logit)",
      main = ""
    )
    axis(1, at = x_pos, labels = groups, las = 2, cex.axis = 0.85)
    abline(h = 0, lty = 2, col = style$neutral)
    for (i in seq_along(groups)) {
      sub <- rows[rows$Group == groups[i], , drop = FALSE]
      sub <- sub[order(seq_len(nrow(sub))), , drop = FALSE]
      tau_vec <- sub$Threshold
      x <- x_pos[i]
      if (length(tau_vec) >= 2L) {
        for (k in seq_len(length(tau_vec) - 1L)) {
          if (!all(is.finite(tau_vec[c(k, k + 1L)]))) next
          if (!isTRUE(sub$StepIndex[k + 1L] - sub$StepIndex[k] == 1L)) next
          col_seg <- if (isTRUE(highlight_disorder) && isTRUE(sub$Disordered[k + 1L])) {
            style$fail
          } else {
            style$accent_primary
          }
          segments(x, tau_vec[k], x, tau_vec[k + 1L],
                             col = col_seg, lwd = 2)
        }
      }
      points(rep(x, length(tau_vec)), tau_vec,
                       pch = 19, col = style$accent_primary)
      text(x + 0.15, tau_vec, labels = sub$Step,
                     cex = 0.7, adj = 0)
      if (any(!is.finite(tau_vec))) text(x, mean(y_lim),
        paste(sub$Step[!is.finite(tau_vec)], "unavailable", collapse = "; "), cex = 0.65)
    }
  }

  out <- new_mfrm_plot_data(
    "threshold_ladder",
    list(
      data = rows,
      model_context = context,
      threshold_coverage = summarize_threshold_order(steps),
      notes = c(category_owner_notes(context),
        "Ladders show one set per step owner, independently of slope ownership. Only consecutive available steps are connected; disorder is not an automatic category-merging rule."),
      n_disorder_groups = n_disorder_groups,
      title = plot_title,
      subtitle = plot_subtitle,
      legend = new_plot_legend(
        label = c("Ordered step", "Disordered crossing"),
        role = c("location", "alert"),
        aesthetic = c("segment", "segment"),
        value = c(style$accent_primary, style$fail)
      ),
      reference_lines = new_reference_lines(
        "h", 0, "Centred logit reference", "dashed", "reference"
      ),
      preset = style$name
    )
  )
  invisible(out)
}


#' Plot per-person fit
#'
#' Per-person diagnostic bubble plot inspired by FACETS Table 6 / KIDMAP
#' summaries. Each bubble represents one person at the intersection of
#' Infit (x) and Outfit (y), sized by total observations and coloured by
#' the standard 0.5/1.5 fit envelope: green when both Infit and Outfit
#' fall in `[lower, upper]`, amber when one statistic is outside, red
#' when both are outside. Set `fit_index = "loglik"` for a ranked view of
#' the report-ready `lz_star` / `lz` index instead.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output. When omitted,
#'   `diagnose_mfrm(fit, residual_pca = "none")` is run internally.
#' @param lower Lower fit threshold (default `0.5`, Linacre 2002).
#' @param upper Upper fit threshold (default `1.5`).
#' @param top_n_label Maximum number of persons whose label is drawn.
#'   The default mean-square view uses largest `|Infit - 1| + |Outfit - 1|`;
#'   `fit_index = "loglik"` uses largest absolute report index. Default `12`.
#' @param fit_index Plot focus. `"meansquare"` keeps the Infit/Outfit bubble
#'   plot. `"loglik"` draws the report index selected by
#'   [compute_person_fit_indices()] (`lz_star` when available, otherwise
#'   `lz` with a caveat).
#' @param preset Visual preset, including `"monochrome"`.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object whose reusable plot data include
#'   `data` with one row per person, `plot_long` for custom R graphics,
#'   `person_fit_indices` from [compute_person_fit_indices()], and compact
#'   flag/status summaries.
#'
#' @section Interpreting output:
#' The default 0.5-1.5 envelope follows Linacre (2002) Rasch
#' Measurement Transactions. Persons in the green centre are
#' within the screening band; amber and red corners are candidates for
#' review (overfit / underfit) using
#' [unexpected_response_table()] for follow-up.
#'
#' @seealso [diagnose_mfrm()], [unexpected_response_table()],
#'   [build_misfit_casebook()].
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_person_fit", run.dontrun = TRUE)
#' \dontrun{
#' library(mfrmr)
#'
#' # Saved synthetic MML/RSM analysis of example_operational.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "operational-mml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#'
#' # Use the diagnostics saved with this fit
#' diagnostics <- example_fit$diagnostics
#'
#' # Compare each person's Infit and Outfit with the reference value of 1
#' plot(fit, type = "person_fit", diagnostics = diagnostics)
#' # Values above 1 indicate more response variation than expected, below 1 less
#' # The reference bands flag patterns for review, not automatic exclusion
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_person_fit <- function(fit,
                            diagnostics = NULL,
                            lower = 0.5,
                            upper = 1.5,
                            top_n_label = 12L,
                            preset = c("standard", "publication", "compact", "monochrome"),
                            draw = TRUE,
                            fit_index = c("meansquare", "loglik")) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  fit_index <- match.arg(fit_index)
  style <- resolve_plot_preset(preset)
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(suppressWarnings(
      diagnose_mfrm(fit, residual_pca = "none")
    ))
  }
  mfrm_results_validate_diagnostics_identity(fit, diagnostics, helper = "plot(x, type = \"person_fit\")")
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower) || lower < 0 ||
      !is.numeric(upper) || length(upper) != 1L || !is.finite(upper) || upper <= lower) {
    stop("Fit thresholds must be finite numbers with 0 <= lower < upper.", call. = FALSE)
  }
  m <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  if (!all(c("Facet", "Level", "Infit", "Outfit") %in% names(m))) {
    stop("Diagnostics measures are missing required columns.",
         call. = FALSE)
  }
  m <- m[as.character(m$Facet) == "Person", , drop = FALSE]
  if (nrow(m) == 0L) {
    stop("No Person rows are available in the diagnostics.",
         call. = FALSE)
  }
  n_col <- if ("N" %in% names(m)) "N" else if ("N.x" %in% names(m)) "N.x" else NA_character_
  m$N <- if (!is.na(n_col)) suppressWarnings(as.numeric(m[[n_col]])) else 1
  in_band <- m$Infit >= lower & m$Infit <= upper
  out_band <- m$Outfit >= lower & m$Outfit <= upper
  available_fit <- is.finite(m$Infit) & is.finite(m$Outfit)
  m$Status <- ifelse(!available_fit, "not_available", ifelse(in_band & out_band, "in_band",
                     ifelse(in_band | out_band, "one_outside", "both_outside")))
  status_color <- c(in_band = style$success,
                    one_outside = style$warn,
                    both_outside = style$fail, not_available = style$neutral)
  m$Color <- unname(status_color[m$Status])
  m$Score <- abs(m$Infit - 1) + abs(m$Outfit - 1)
  m <- m[order(-m$Score), , drop = FALSE]

  person_fit_indices <- tryCatch(
    compute_person_fit_indices(diagnostics, fit = fit),
    error = function(e) e
  )
  person_fit_status <- data.frame(
    Available = !inherits(person_fit_indices, "error"),
    Status = if (inherits(person_fit_indices, "error")) "error" else "available",
    Message = if (inherits(person_fit_indices, "error")) {
      conditionMessage(person_fit_indices)
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )
  if (inherits(person_fit_indices, "error")) {
    person_fit_indices <- data.frame(Person = character(), stringsAsFactors = FALSE)
  } else {
    person_fit_indices <- as.data.frame(person_fit_indices, stringsAsFactors = FALSE)
  }
  if (nrow(person_fit_indices) > 0L) {
    merge_tbl <- person_fit_indices
    names(merge_tbl)[names(merge_tbl) == "N"] <- "PersonFitN"
    m$.row_id <- seq_len(nrow(m))
    m <- merge(m, merge_tbl,
               by.x = "Level", by.y = "Person",
               all.x = TRUE, sort = FALSE)
    m <- m[order(m$.row_id), , drop = FALSE]
    m$.row_id <- NULL
  }
  person_cols <- c(
    "LogLik", "lz", "lz_star", "lz_star_status", "lz_star_c",
    "lz_star_variance", "lz_flag_5pct", "lz_flag_1pct",
    "lz_star_flag_5pct", "lz_star_flag_1pct", "ReportIndex",
    "ReportValue", "ReportFlagLevel", "ReportFlag", "ReviewStatus",
    "ReviewReason", "ReportCaveat"
  )
  for (nm in person_cols) {
    if (!nm %in% names(m)) {
      m[[nm]] <- if (nm %in% c("lz_flag_5pct", "lz_flag_1pct",
                               "lz_star_flag_5pct", "lz_star_flag_1pct",
                               "ReportFlag")) {
        NA
      } else if (nm %in% c("lz_star_status", "ReportIndex",
                           "ReportFlagLevel", "ReviewStatus",
                           "ReviewReason", "ReportCaveat")) {
        NA_character_
      } else {
        NA_real_
      }
    }
  }
  review_color <- c(
    not_flagged = style$success,
    review_5pct = style$warn,
    review_1pct = style$fail,
    not_available = style$neutral
  )
  m$ReviewColor <- unname(review_color[as.character(m$ReviewStatus)])
  m$ReviewColor[is.na(m$ReviewColor)] <- style$neutral

  plot_title <- "Person fit"
  plot_subtitle <- if (identical(fit_index, "meansquare")) {
    sprintf(
      "Infit and Outfit per person (heuristic review band [%g, %g], Linacre 2002)",
      lower, upper
    )
  } else {
    "Conditional lz* or uncorrected lz; normal cutoffs are uncalibrated screening references"
  }

  plot_available <- if (fit_index == "meansquare") is.finite(m$Infit) & is.finite(m$Outfit) else is.finite(m$ReportValue)
  plot_subtitle <- paste0(plot_subtitle, sprintf("; %d of %d persons available", sum(plot_available), nrow(m)))

  z_5pct <- stats::qnorm(0.975)
  z_1pct <- stats::qnorm(0.995)

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_mar <- graphics::par("mar")
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    graphics::par(mar = c(max(old_mar[1], 5), old_mar[2], 1.5, old_mar[4]))
    if (identical(fit_index, "meansquare") && !any(plot_available)) {
      plot.new()
      text(.5, .5, "No complete Infit/Outfit pairs are available")
    } else if (identical(fit_index, "meansquare")) {
      cex_size <- 0.6 + 1.6 * sqrt(m$N / max(m$N, na.rm = TRUE))
      x_rng <- range(c(m$Infit, lower, upper), finite = TRUE)
      y_rng <- range(c(m$Outfit, lower, upper), finite = TRUE)
      x_pad <- if (diff(x_rng) > 0) c(-0.05, 0.05) * diff(x_rng) else c(-0.1, 0.1)
      y_pad <- if (diff(y_rng) > 0) c(-0.05, 0.05) * diff(y_rng) else c(-0.1, 0.1)
      plot(
        x = m$Infit, y = m$Outfit,
        xlab = "Infit MnSq", ylab = "Outfit MnSq",
        main = "",
        pch = 21, bg = m$Color, col = "white",
        cex = cex_size,
        xlim = x_rng + x_pad, ylim = y_rng + y_pad
      )
      abline(h = 1, v = 1, lty = 3, col = style$neutral)
      abline(h = c(lower, upper), v = c(lower, upper),
                       lty = 2, col = style$grid)
      n_lbl <- min(top_n_label, nrow(m))
      if (n_lbl > 0L) {
        text(
          x = m$Infit[seq_len(n_lbl)],
          y = m$Outfit[seq_len(n_lbl)],
          labels = as.character(m$Level[seq_len(n_lbl)]),
          cex = 0.7, pos = 4, offset = 0.5
        )
      }
    } else {
      ranked <- m[is.finite(m$ReportValue), , drop = FALSE]
      ranked <- ranked[order(ranked$ReportValue), , drop = FALSE]
      if (nrow(ranked) == 0L) {
        plot.new()
        text(.5, .5, "No finite lz/lz* report index available")
      } else {
        y <- seq_len(nrow(ranked))
        x_rng <- range(c(ranked$ReportValue, -z_1pct, -z_5pct, z_5pct, z_1pct),
                       finite = TRUE)
        x_pad <- if (diff(x_rng) > 0) c(-0.08, 0.08) * diff(x_rng) else c(-0.5, 0.5)
        plot(
          x = ranked$ReportValue, y = y,
          xlab = "Person-fit report index",
          ylab = "Persons ordered by report index",
          yaxt = "n",
          main = "",
          pch = 21, bg = ranked$ReviewColor, col = "white",
          cex = 0.85 + 1.2 * sqrt(ranked$N / max(ranked$N, na.rm = TRUE)),
          xlim = x_rng + x_pad
        )
        abline(v = 0, lty = 3, col = style$neutral)
        abline(v = c(-z_5pct, z_5pct), lty = 2, col = style$grid)
        abline(v = c(-z_1pct, z_1pct), lty = 3, col = style$grid)
        label_order <- order(abs(ranked$ReportValue), decreasing = TRUE)
        n_lbl <- min(top_n_label, length(label_order))
        if (n_lbl > 0L) {
          idx <- label_order[seq_len(n_lbl)]
          text(
            x = ranked$ReportValue[idx],
            y = y[idx],
            labels = as.character(ranked$Level[idx]),
            cex = 0.7, pos = 4, offset = 0.5
          )
        }
      }
    }
  }

  payload <- data.frame(
    Person = as.character(m$Level),
    Infit = m$Infit,
    Outfit = m$Outfit,
    N = m$N,
    PersonFitN = m$PersonFitN %||% rep(NA_integer_, nrow(m)),
    Status = m$Status,
    LogLik = suppressWarnings(as.numeric(m$LogLik)),
    lz = suppressWarnings(as.numeric(m$lz)),
    lz_star = suppressWarnings(as.numeric(m$lz_star)),
    lz_star_status = as.character(m$lz_star_status),
    ReportIndex = as.character(m$ReportIndex),
    ReportValue = suppressWarnings(as.numeric(m$ReportValue)),
    ReportFlagLevel = as.character(m$ReportFlagLevel),
    ReportFlag = as.logical(m$ReportFlag),
    ReviewStatus = as.character(m$ReviewStatus),
    ReviewReason = as.character(m$ReviewReason),
    ReportCaveat = as.character(m$ReportCaveat),
    stringsAsFactors = FALSE
  )
  plot_long <- rbind(
    data.frame(Person = payload$Person, Metric = "Infit", Value = payload$Infit,
               MetricFamily = "mean_square", DisplayedByDefault = identical(fit_index, "meansquare"),
               Status = payload$Status, ReviewStatus = payload$ReviewStatus,
               stringsAsFactors = FALSE),
    data.frame(Person = payload$Person, Metric = "Outfit", Value = payload$Outfit,
               MetricFamily = "mean_square", DisplayedByDefault = identical(fit_index, "meansquare"),
               Status = payload$Status, ReviewStatus = payload$ReviewStatus,
               stringsAsFactors = FALSE),
    data.frame(Person = payload$Person, Metric = "lz", Value = payload$lz,
               MetricFamily = "log_likelihood", DisplayedByDefault = FALSE,
               Status = payload$Status, ReviewStatus = payload$ReviewStatus,
               stringsAsFactors = FALSE),
    data.frame(Person = payload$Person, Metric = "lz_star", Value = payload$lz_star,
               MetricFamily = "log_likelihood", DisplayedByDefault = FALSE,
               Status = payload$Status, ReviewStatus = payload$ReviewStatus,
               stringsAsFactors = FALSE),
    data.frame(Person = payload$Person, Metric = "ReportValue", Value = payload$ReportValue,
               MetricFamily = "log_likelihood", DisplayedByDefault = identical(fit_index, "loglik"),
               Status = payload$Status, ReviewStatus = payload$ReviewStatus,
               stringsAsFactors = FALSE)
  )
  flag_summary <- data.frame(
    Area = c("mean_square", "report_index"),
    Rows = c(nrow(payload), nrow(payload)),
    AvailableRows = c(sum(is.finite(payload$Infit) & is.finite(payload$Outfit)), sum(is.finite(payload$ReportValue))),
    UnavailableRows = c(sum(!is.finite(payload$Infit) | !is.finite(payload$Outfit)), sum(!is.finite(payload$ReportValue))),
    FlaggedRows = c(
      if (any(payload$Status != "not_available")) sum(payload$Status %in% c("one_outside", "both_outside")) else NA_integer_,
      if (any(is.finite(payload$ReportValue))) sum(payload$ReportFlag %in% TRUE) else NA_integer_
    ),
    Review1PctRows = c(
      NA_integer_,
      sum(payload$ReviewStatus == "review_1pct", na.rm = TRUE)
    ),
    Review5PctRows = c(
      NA_integer_,
      sum(payload$ReviewStatus == "review_5pct", na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  reference_lines <- if (identical(fit_index, "meansquare")) {
    new_reference_lines(
      axis = c("h", "v", "h", "v"),
      value = c(lower, lower, upper, upper),
      label = rep("Fit envelope", 4),
      linetype = rep("dashed", 4),
      role = rep("threshold", 4)
    )
  } else {
    new_reference_lines(
      axis = rep("v", 4),
      value = c(-z_1pct, -z_5pct, z_5pct, z_1pct),
      label = c("Normal reference -2.576", "Normal reference -1.960", "Normal reference 1.960", "Normal reference 2.576"),
      linetype = c("dotted", "dashed", "dashed", "dotted"),
      role = rep("person-fit z threshold", 4)
    )
  }
  out <- new_mfrm_plot_data(
    "person_fit",
    list(
      data = payload,
      plot_long = plot_long,
      person_fit_indices = person_fit_indices,
      person_fit_status = person_fit_status,
      flag_summary = flag_summary,
      plot_settings = data.frame(
        FitIndex = fit_index,
        Lower = lower,
        Upper = upper,
        Z5Pct = z_5pct,
        Z1Pct = z_1pct,
        Preset = style$name,
        stringsAsFactors = FALSE
      ),
      lower = lower,
      upper = upper,
      fit_index = fit_index,
      title = plot_title,
      subtitle = plot_subtitle,
      legend = new_plot_legend(
        label = c("Within band", "One outside", "Both outside", "Unavailable",
                  "No threshold exceeded", "Above 1.960", "Above 2.576"),
        role = rep("status", 7),
        aesthetic = rep("point", 7),
        value = c(unname(status_color),
                  unname(review_color[c("not_flagged", "review_5pct", "review_1pct")]))
      ),
      reference_lines = reference_lines,
      preset = style$name
    )
  )
  invisible(out)
}


#' Plot facet locations with qualified normal uncertainty bands
#'
#' Ranks the levels of a chosen facet by fitted location and
#' draws available normal uncertainty bounds around each point
#' estimate. Optional descriptive bands mark absolute distances of `0.5`
#' and `1.0` logit from zero; they are display aids, not calibration rules.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output. When omitted,
#'   `diagnose_mfrm(fit, residual_pca = "none")` is run internally.
#'   Supplied diagnostics must match the fitted analysis and current precision
#'   contract; they are reused without covariance recomputation.
#' @param facet Facet name to plot (default `"Rater"`). Any non-Person
#'   facet name is accepted.
#' @param ci_level Confidence level used for the whiskers (default
#'   `0.95`). Bounds use `+/- z * ModelSE` (or `SE` only when the
#'   `ModelSE` column is absent). A missing `ModelSE` is not replaced rowwise.
#' @param show_bands Logical, default `FALSE`. Explicit `TRUE` draws neutral
#'   gray `+/-0.5` and `+/-1.0` logit guides and describes them in the
#'   returned annotations. These distances are not equivalence bounds.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object. Its `data$data` table contains
#'   columns `Level`, `Estimate`, `SE`, `CI_Lower`, `CI_Upper`, and `Band`,
#'   where `Band` contains numerical absolute-location ranges, not judgments
#'   about severity or quality. The ranges remain in the table when hidden.
#'   Other columns include `Fixed`, `IntervalStatus`, `SEColumn`, `SE_Method`,
#'   `PrecisionTier`, `CIUse`, `CIEligible`, `SupportsFormalInference` and
#'   the normal interval method/level. Source precision labels are retained;
#'   missing precision metadata cannot qualify an interval for inference.
#'   The enclosing `data` list retains the facet, confidence level and plot
#'   annotations and diagnostics identity check (`source_status`); keep these
#'   with the table when reporting the results. Fit readiness and interpretation
#'   notes are returned. Titles and notes are not automatically drawn.
#'
#' @section Interpreting output:
#' Zero is the sum-to-zero reference for the default centered facet; describe
#' any different constraints or anchors used in the fit. With the default
#' negative facet orientation, higher estimates reduce expected coded scores.
#' This means stricter scoring when higher scores represent better performance;
#' check the rubric direction before using that substantive interpretation.
#' The axis names the selected facet's location rather than assuming that
#' lower coded scores are less favorable. Returned notes state its fitted sign.
#' The optional guides and numerical `Band` ranges describe distance from zero,
#' not the sign of severity, operational interchangeability, or a need for
#' training. No guide bands are drawn by default. Pairwise claims require the
#' uncertainty of the contrast; check the SE basis in the supplied diagnostics.
#' Two locations at -0.5 and +0.5 are one logit apart despite lying in the same
#' range. A changed origin can change these ranges without changing pairwise
#' differences. Use [analyze_facet_equivalence()] with a justified practical
#' bound for equivalence; a non-significant difference does not establish it.
#'
#' Here "profile" means an ordered display of facet locations, not a
#' likelihood-profile confidence interval. Ordinary JML bands remain
#' exploratory; MML fallback, regularized or review-only bands retain their
#' restrictions. Plotting does not validate normal-approximation coverage.
#' Fixed values use diamonds without sampling intervals. Crosses show finite
#' points with unavailable or nonpositive SEs; unavailable estimates retain
#' labeled rows. The horizontal range includes finite points even without SEs.
#' [as_ggplot()] preserves these states, saved endpoints, guide bands and notes
#' without fitting a model or recomputing an interval. Older saved plot payloads
#' can reproduce their stored geometry but cannot acquire missing source or
#' uncertainty qualifications through conversion.
#' Older saved plots retain their original colored guides and `Band` labels
#' (`gentle`, `moderate`, `strict`). Those labels never established equivalence.
#' Regenerate the plot from the saved fit and matching diagnostics to obtain
#' the neutral display; this does not require a new fit.
#'
#' @seealso [diagnose_mfrm()], [analyze_facet_equivalence()],
#'   [plot(x, type = "facet_equivalence")][plot.mfrm_fit].
#'
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept rater severity
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_rater_severity_profile", run.dontrun = TRUE)
#' \dontrun{
#' library(mfrmr)
#'
#' # Saved synthetic MML/RSM analysis of example_operational.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "operational-mml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#'
#' # Use the diagnostics saved with this fit
#' diagnostics <- example_fit$diagnostics
#'
#' # Compare signed severity estimates and their intervals
#' severity <- plot(fit, type = "rater_severity", diagnostics = diagnostics, show_bands = FALSE)
#' severity$data$data[, c("Level", "Estimate", "SE", "CI_Lower", "CI_Upper")]
#' # Higher estimates mean stricter ratings with this example's default orientation
#' # The optional magnitude bands are omitted from this first comparison
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_rater_severity_profile <- function(fit,
                                        diagnostics = NULL,
                                        facet = "Rater",
                                        ci_level = 0.95,
                                        show_bands = FALSE,
                                        preset = c("standard", "publication", "compact", "monochrome"),
                                        draw = TRUE) {
  stop_if_product_slopes(fit, "plot(x, type = \"rater_severity\")")
  if (inherits(fit, "mfrm_imported_fit")) {
    stop("Severity profiles require a native, supported mfrm_fit object.", call. = FALSE)
  }

  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.character(facet) || length(facet) != 1L || is.na(facet) ||
      !nzchar(facet) || facet == "Person") {
    stop("`facet` must name one non-Person facet.", call. = FALSE)
  }
  if (!is.logical(show_bands) || length(show_bands) != 1L || is.na(show_bands)) {
    stop("`show_bands` must be TRUE or FALSE.", call. = FALSE)
  }
  readiness <- .mfrm_fit_plot_readiness(fit)
  style <- resolve_plot_preset(preset)
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(suppressWarnings(
      diagnose_mfrm(fit, residual_pca = "none")
    ))
  }
  source_status <- mfrm_results_validate_diagnostics_identity(
    fit, diagnostics, helper = 'plot(x, type = "rater_severity")')
  m <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  m <- m[as.character(m$Facet) == as.character(facet), , drop = FALSE]
  if (nrow(m) == 0L) {
    stop(sprintf("No rows for facet '%s' in diagnostics$measures.", facet),
         call. = FALSE)
  }
  m$Estimate <- suppressWarnings(as.numeric(m$Estimate))
  se_column <- if ("ModelSE" %in% names(m)) "ModelSE" else if ("SE" %in% names(m)) "SE" else NA_character_
  m$SE <- if (is.na(se_column)) rep(NA_real_, nrow(m)) else suppressWarnings(as.numeric(m[[se_column]]))
  defaults <- list(SE_Method = "Not recorded", PrecisionTier = "unknown",
    SupportsFormalInference = FALSE, SEUse = "review_before_reporting",
    CIBasis = "Interval basis not recorded", CIUse = "review_before_reporting",
    CIEligible = FALSE, CILabel = "Interval basis not recorded; review source diagnostics")
  complete_basis <- all(names(defaults) %in% names(m))
  for (nm in names(defaults)) if (!nm %in% names(m)) m[[nm]] <- rep(defaults[[nm]], nrow(m))
  if (!complete_basis) for (nm in c("PrecisionTier", "SupportsFormalInference", "SEUse", "CIUse", "CIEligible", "CILabel")) {
    m[[nm]] <- rep(defaults[[nm]], nrow(m))
  }
  m$Fixed <- if ("Fixed" %in% names(m)) m$Fixed %in% TRUE else rep(FALSE, nrow(m))
  m$IntervalStatus <- ifelse(!is.finite(m$Estimate), "estimate_unavailable",
    ifelse(m$Fixed, "fixed", ifelse(!is.finite(m$SE), "se_unavailable",
      ifelse(m$SE <= 0, "se_nonpositive", "normal_approximation"))))
  available <- m$IntervalStatus == "normal_approximation"
  z_ci <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
  m$CI_Lower <- ifelse(available, m$Estimate - z_ci * m$SE, NA_real_)
  m$CI_Upper <- ifelse(available, m$Estimate + z_ci * m$SE, NA_real_)
  m$CIEligible <- available & m$CIEligible %in% TRUE
  m$SupportsFormalInference <- available & m$SupportsFormalInference %in% TRUE
  m$CIUse[!available] <- ifelse(m$Fixed[!available], "not_applicable", "unavailable")
  m$CI_Method <- ifelse(m$Fixed, "Not applicable to fixed values", "Normal approximation")
  m$CI_Level <- ci_level
  m$SEColumn <- se_column
  m$Band <- ifelse(!is.finite(m$Estimate), NA_character_,
    ifelse(abs(m$Estimate) <= 0.5, "|location| <= 0.5",
      ifelse(abs(m$Estimate) <= 1.0, "0.5 < |location| <= 1.0", "|location| > 1.0")))
  m <- m[order(ifelse(is.finite(m$Estimate), m$Estimate, NA_real_), na.last = TRUE), , drop = FALSE]
  payload <- as.data.frame(m[, c("Level", "Estimate", "SE", "CI_Lower", "CI_Upper", "Band",
    "Fixed", "IntervalStatus", "SEColumn", names(defaults), "CI_Level", "CI_Method"), drop = FALSE])
  payload$Level <- as.character(payload$Level)
  geometry <- .mfrm_severity_geometry(payload)
  band_colors <- rep(style$neutral, 3L)
  facet_sign <- unname(fit$config$facet_signs[facet])
  axis_label <- paste0(facet, " location (logit)")
  plot_title <- sprintf("%s severity profile", facet)
  if (!isTRUE(readiness$ready)) plot_title <- paste("REVIEW ONLY -", plot_title)
  plot_subtitle <- sprintf(
    "Zero reference; %g%% normal uncertainty; inspect CIUse and SE_Method%s",
    round(100 * ci_level),
    if (isTRUE(show_bands)) "; shaded guides at +/-0.5 and +/-1.0 logits" else ""
  )

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_mar <- graphics::par("mar")
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    label_width <- max(graphics::strwidth(as.character(m$Level), units = "inches", cex = 0.85))
    label_margin <- label_width / (graphics::par("csi") * graphics::par("mex")) + 1.2
    if (label_width + graphics::par("mai")[4] + 1.25 >= graphics::par("fin")[1]) {
      stop("Facet labels need more horizontal space; use a wider graphics device or inspect the table with `draw = FALSE`.", call. = FALSE)
    }
    graphics::par(mar = c(max(old_mar[1], 5),
                          max(old_mar[2], label_margin), old_mar[3:4]))
    y <- seq_len(nrow(m))
    xrange <- geometry$xrange
    plot(
      x = m$Estimate, y = y,
      type = "n",
      xlim = xrange,
      yaxt = "n",
      xlab = axis_label,
      ylab = "",
      main = ""
    )
    if (isTRUE(show_bands)) {
      usr <- graphics::par("usr")
      rect(-0.5, usr[3], 0.5, usr[4], border = NA,
                     col = grDevices::adjustcolor(band_colors[1L], alpha.f = 0.10))
      rect(-1.0, usr[3], -0.5, usr[4], border = NA,
                     col = grDevices::adjustcolor(band_colors[2L], alpha.f = 0.08))
      rect(0.5, usr[3], 1.0, usr[4], border = NA,
                     col = grDevices::adjustcolor(band_colors[3L], alpha.f = 0.08))
    }
    abline(v = 0, lty = 2, col = style$neutral)
    valid <- geometry$interval
    if (any(valid)) {
      segments(m$CI_Lower[valid], y[valid],
                         m$CI_Upper[valid], y[valid],
                         col = style$accent_primary, lwd = 2)
    }
    points(m$Estimate[geometry$shown], y[geometry$shown],
      pch = geometry$shape[geometry$shown], col = style$accent_primary)
    if (any(!geometry$shown)) text(mean(xrange), y[!geometry$shown],
      "Unavailable", col = style$neutral, cex = 0.85)
    axis(2, at = y, labels = as.character(m$Level), las = 1,
                   cex.axis = 0.85)
  }

  out <- new_mfrm_plot_data(
    "rater_severity_profile",
    list(
      data = payload,
      facet = facet,
      ci_level = ci_level,
      show_bands = show_bands,
      band_colors = band_colors,
      source_status = source_status,
      facet_sign = facet_sign,
      axis_label = axis_label,
      notes = data.frame(Type = "Interpretation", Text = c(
        readiness$detail,
        "Whiskers use the stored SE and a normal approximation, not likelihood-profile intervals. CIUse and SE_Method qualify their use; exploratory JML and fallback/regularized bands do not become formal intervals.",
        "Diamonds mark fixed values; crosses mark finite estimates without intervals. Unavailable estimates remain labeled rows. These states do not indicate zero uncertainty.",
        .mfrm_score_direction_notes(fit, facets = facet, include_person = FALSE)$Text,
        "Zero refers to the fitted scale constraints/anchors, not an external standard of correct scoring. Band ranges and optional gray guides describe distance from zero, not homogeneity, practical equivalence or rater quality. Changing the origin can change a Band without changing pairwise differences.",
        "Locations at -0.5 and +0.5 are one logit apart despite sharing a Band. Use analyze_facet_equivalence() with a justified bound and eligible contrast uncertainty; a non-significant difference does not establish equivalence.",
        "Individual intervals are not pairwise tests. Review common ratings and assignment before comparing raters."
      )),
      title = plot_title,
      subtitle = plot_subtitle,
      legend = new_plot_legend(
        label = c("Estimate", "Normal interval/band", "Fixed value", "No interval", if (isTRUE(show_bands)) c("+/-0.5 guide", "+/-1.0 guide")),
        role = c("location", "uncertainty", "fixed", "unavailable", if (isTRUE(show_bands)) c("band", "band")),
        aesthetic = c("point", "segment", "shape", "shape", if (isTRUE(show_bands)) c("fill", "fill")),
        value = c(style$accent_primary, style$accent_primary, "5", "4",
                  if (isTRUE(show_bands)) band_colors[1:2])
      ),
      reference_lines = new_reference_lines(
        "v", 0, "Zero reference", "dashed", "reference"
      ),
      preset = style$name
    )
  )
  out <- .mfrm_attach_plot_readiness(out, readiness)
  invisible(out)
}

.mfrm_severity_geometry <- function(table) {
  fixed <- if ("Fixed" %in% names(table)) table$Fixed %in% TRUE else rep(FALSE, nrow(table))
  shown <- is.finite(table$Estimate)
  interval <- shown & !fixed & is.finite(table$CI_Lower) & is.finite(table$CI_Upper) &
    table$CI_Lower <= table$CI_Upper
  list(shown = shown, interval = interval,
    shape = ifelse(fixed, 5, ifelse(interval, 19, 4)),
    xrange = range(c(table$Estimate[shown], table$CI_Lower[interval],
      table$CI_Upper[interval], -1.05, 1.05), finite = TRUE))
}


#' Summary plot of differential functioning effect sizes
#'
#' Compact effect-size summary for a [analyze_dff()] / [analyze_dif()]
#' result. Shows each contrast's signed effect size as a horizontal bar
#' with a vertical reference at zero and a neutral colour. Residual comparisons
#' do not provide differential-functioning tests or classifications.
#'
#' @param x Output from [analyze_dff()] or [analyze_dif()].
#' @param top_n Maximum rows shown (default `30`).
#' @param sort_by `"abs_effect"` (default), `"effect"`, or
#'   `"classification"`.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#' @param ci_level Optional confidence level for approximate normal
#'   intervals drawn from `Effect +/- z * SE` when finite standard errors are
#'   available for linked refits. Residual comparisons require `NULL` (default)
#'   because their interval uncertainty is not established.
#' @param effect_thresholds Optional numeric vector of absolute effect-size
#'   guide lines to draw at `+/- threshold`. These are display aids, not ETS
#'   classification boundaries.
#' @param effect_axis_label Optional x-axis label override. When `NULL`, the
#'   label is chosen from the DFF method.
#'
#' @return An `mfrm_plot_data` object whose `data` slot contains
#'   columns `Pair`, `Effect`, `SE`, `Classification`, `Color`.
#'
#' @section Interpreting output:
#' Bars are anchored at zero. Width corresponds to effect size on the
#' contrast's native scale. For `method = "residual"`, this is the
#' observed-minus-expected average difference between groups, in score units. For
#' `method = "refit"`, this is the subgroup parameter difference on the
#' fitted logit scale when linking support allows a comparable contrast.
#' Residual differences do not isolate differential functioning. Linked refit
#' intervals omit estimated-anchor uncertainty and cross-refit covariance.
#'
#' @seealso [analyze_dff()], [analyze_dif()], [plot(x, type = "heatmap")][plot.mfrm_dif].
#'
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept DFF DIF
#'
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_dif_summary <- function(x,
                             top_n = 30L,
                             sort_by = c("abs_effect", "effect", "classification"),
                             preset = c("standard", "publication", "compact", "monochrome"),
                             draw = TRUE,
                             ci_level = NULL,
                             effect_thresholds = NULL,
                             effect_axis_label = NULL) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(x, c("mfrm_dff", "mfrm_dif"))) {
    stop("`x` must be output from analyze_dff() or analyze_dif().",
         call. = FALSE)
  }
  validate_dff_residual_output(x)
  sort_by <- match.arg(sort_by)
  top_n <- .validate_dff_count_arg(top_n, "top_n")
  ci_level <- .validate_dff_probability(ci_level, "ci_level")
  effect_thresholds <- .validate_dff_threshold_vector(effect_thresholds,
                                                      "effect_thresholds")
  if (!is.null(effect_axis_label) &&
      (!is.character(effect_axis_label) || length(effect_axis_label) != 1L ||
       is.na(effect_axis_label) || !nzchar(effect_axis_label))) {
    stop("`effect_axis_label` must be a single non-empty character string.",
         call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  tbl <- as.data.frame(x$dif_table, stringsAsFactors = FALSE)
  if (nrow(tbl) == 0L) {
    stop("Empty differential-functioning table.", call. = FALSE)
  }
  effect_col <- if ("Contrast" %in% names(tbl)) "Contrast" else "Effect"
  if (!effect_col %in% names(tbl)) {
    stop("DFF/DIF table is missing a `Contrast` (or `Effect`) column.",
         call. = FALSE)
  }
  tbl$Effect <- suppressWarnings(as.numeric(tbl[[effect_col]]))
  tbl$SE <- if ("SE" %in% names(tbl)) suppressWarnings(as.numeric(tbl$SE)) else NA_real_
  pairs <- as.data.frame(x$dif_table, stringsAsFactors = FALSE)
  pairs$PlotStatus <- ifelse(is.finite(tbl$Effect), "Not selected by top_n", "Contrast unavailable")
  tbl$.PlotRow <- seq_len(nrow(tbl))
  n_fixed <- sum(tbl[["Classification"]] == "Fixed target; no comparison", na.rm = TRUE)
  tbl <- tbl[is.finite(tbl$Effect), , drop = FALSE]
  if (nrow(tbl) == 0L) {
    stop("No finite contrast values to plot.", call. = FALSE)
  }
  method <- x$config$method %||% unique(as.character(tbl$Method))[1] %||% NA_character_
  classification_system <- unique(as.character(tbl$ClassificationSystem %||% NA_character_))
  classification_system <- classification_system[!is.na(classification_system)]
  classification_system <- classification_system[1] %||% NA_character_
  axis_label <- effect_axis_label %||% .dff_effect_axis_label(method)
  pair_cols <- intersect(c("Level", "Group1", "Group2"), names(tbl))
  if (length(pair_cols) > 0L) {
    tbl$Pair <- do.call(paste, c(tbl[pair_cols], sep = " | "))
  } else {
    tbl$Pair <- as.character(seq_len(nrow(tbl)))
  }
  tbl$Classification <- if ("Classification" %in% names(tbl)) {
    as.character(tbl$Classification)
  } else if ("ETS" %in% names(tbl)) {
    as.character(tbl$ETS)
  } else {
    NA_character_
  }
  classification_color <- c(
    A = style$success, B = style$warn, C = style$fail,
    Negligible = style$success, Moderate = style$warn, Large = style$fail
  )
  tbl$Color <- ifelse(
    is.na(tbl$Classification),
    style$neutral,
    unname(classification_color[as.character(tbl$Classification)])
  )
  tbl$Color[is.na(tbl$Color)] <- style$neutral
  ord <- switch(
    sort_by,
    abs_effect = order(-abs(tbl$Effect), na.last = TRUE),
    effect = order(tbl$Effect, na.last = TRUE),
    classification = order(tbl$Classification, -abs(tbl$Effect), na.last = TRUE)
  )
  tbl <- tbl[ord, , drop = FALSE]
  tbl <- tbl[seq_len(min(nrow(tbl), top_n)), , drop = FALSE]
  pairs$PlotStatus[tbl$.PlotRow] <- "Displayed"
  if (!is.null(ci_level) && identical(method, "residual")) {
    stop("Confidence intervals are unavailable for residual differences. Omit `ci_level` to plot the differences.", call. = FALSE)
  }
  if (!is.null(ci_level)) {
    z <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
    finite_se <- is.finite(tbl$SE)
    tbl$CI_Lower <- ifelse(finite_se, tbl$Effect - z * tbl$SE, NA_real_)
    tbl$CI_Upper <- ifelse(finite_se, tbl$Effect + z * tbl$SE, NA_real_)
  } else {
    tbl$CI_Lower <- NA_real_
    tbl$CI_Upper <- NA_real_
  }
  plot_title <- if (identical(method, "residual")) "Group residual differences" else "Subgroup location differences"
  plot_subtitle <- sprintf(
    "%d comparisons shown",
    nrow(tbl)
  )
  if (identical(method, "residual")) {
    plot_subtitle <- paste(plot_subtitle, "Residual differences do not isolate differential functioning.", sep = ". ")
  }

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_mar <- graphics::par("mar")
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    label_width <- max(graphics::strwidth(tbl$Pair, units = "inches", cex = 0.7))
    label_margin <- label_width / (graphics::par("csi") * graphics::par("mex")) + 1.2
    if (label_width + graphics::par("mai")[4] + 1.25 >= graphics::par("fin")[1]) {
      stop("Comparison labels need more horizontal space; use a wider graphics device or inspect the table with `draw = FALSE`.", call. = FALSE)
    }
    graphics::par(mar = c(max(old_mar[1], 5),
                          max(old_mar[2], label_margin), 1.5, old_mar[4]))
    y <- rev(seq_len(nrow(tbl)))
    x_values <- c(tbl$Effect, tbl$CI_Lower, tbl$CI_Upper, 0,
                  effect_thresholds, -effect_thresholds)
    x_range <- range(x_values, finite = TRUE)
    if (!all(is.finite(x_range)) || diff(x_range) == 0) {
      x_range <- x_range[1] + c(-1, 1)
    }
    xlim <- x_range + c(-0.05, 0.05) * diff(x_range)
    plot(
      x = tbl$Effect, y = y, type = "n",
      xlim = xlim, yaxt = "n",
      xlab = axis_label, ylab = "",
      main = ""
    )
    abline(v = 0, lty = 2, col = style$neutral)
    if (length(effect_thresholds) > 0L) {
      for (thr in effect_thresholds) {
        abline(v = c(-thr, thr), lty = 3, col = style$warn)
      }
    }
    if (!is.null(ci_level) && any(is.finite(tbl$CI_Lower) & is.finite(tbl$CI_Upper))) {
      segments(tbl$CI_Lower, y, tbl$CI_Upper, y,
                         col = style$neutral, lwd = 1)
    }
    segments(0, y, tbl$Effect, y,
                       col = tbl$Color, lwd = 5)
    points(tbl$Effect, y, pch = 19, col = tbl$Color)
    axis(2, at = y, labels = tbl$Pair, las = 1, cex.axis = 0.7)
  }

  payload <- data.frame(
    Pair = as.character(tbl$Pair),
    Effect = tbl$Effect,
    SE = tbl$SE,
    CI_Lower = tbl$CI_Lower,
    CI_Upper = tbl$CI_Upper,
    Classification = tbl$Classification,
    ClassificationSystem = if ("ClassificationSystem" %in% names(tbl)) {
      as.character(tbl$ClassificationSystem)
    } else {
      NA_character_
    },
    Color = tbl$Color,
    stringsAsFactors = FALSE
  )
  threshold_lines <- if (length(effect_thresholds) > 0L) {
    labs <- names(effect_thresholds)
    labs[is.na(labs) | !nzchar(labs)] <- paste0("Effect threshold ",
                                                effect_thresholds[is.na(labs) | !nzchar(labs)])
    new_reference_lines(
      "v",
      rep(c(-1, 1), each = length(effect_thresholds)) *
        rep(effect_thresholds, times = 2),
      rep(labs, times = 2),
      "dotted",
      "threshold"
    )
  } else {
    new_reference_lines()
  }
  out <- new_mfrm_plot_data(
    "dif_summary",
    list(
      data = payload,
      pairs = pairs,
      cases = x$cases,
      analysis_config = x$config,
      comparison_review = x$comparison_review,
      retention_note = sprintf(
        "%d of %d source comparisons displayed; %d unavailable and %d not selected by top_n. See pairs$PlotStatus and the source comparison conditions; omission does not mean no difference.",
        nrow(tbl), nrow(pairs), sum(pairs$PlotStatus == "Contrast unavailable"),
        sum(pairs$PlotStatus == "Not selected by top_n")),
      notes = data.frame(Type = "Interpretation", Text = if (identical(method, "residual")) {
        "Residual differences do not isolate differential functioning; no test or confidence interval is provided."
      } else paste0("Linked differences condition on the estimated anchors.",
        if (n_fixed > 0L) paste0(" ", n_fixed,
          " fixed-target row(s) have no contrast and are omitted; equal anchors do not establish invariance.")
        else ""), stringsAsFactors = FALSE),
      sort_by = sort_by,
      title = plot_title,
      subtitle = plot_subtitle,
      legend = new_plot_legend(
        label = "Group difference",
        role = "estimate", aesthetic = "bar", value = style$neutral
      ),
      reference_lines = rbind(
        new_reference_lines("v", 0, "Zero contrast", "dashed", "reference"),
        threshold_lines
      ),
      interpretation_guide = .dff_interpretation_guide(
        metric = "summary",
        method = method,
        classification_system = classification_system,
        effect_thresholds = effect_thresholds
      ),
      gpcm_boundary = x$gpcm_boundary %||% data.frame(),
      settings = list(
        ci_level = ci_level,
        effect_thresholds = effect_thresholds,
        effect_axis_label = axis_label
      ),
      preset = style$name
    )
  )
  guide <- out$data$interpretation_guide
  out$data$notes <- rbind(.mfrm_plot_annotation_notes(out$data),
    data.frame(Type = guide$Item, Text = paste(guide$Meaning, guide$ReportingNote)))
  if (!is.null(x$comparison_review)) out$data$notes <- rbind(out$data$notes,
    data.frame(Type = "Comparison review", Text = "A separate prospective interaction review is retained in comparison_review; its states do not classify the displayed residual/refit effects."))
  invisible(out)
}


#' Manuscript-oriented four-panel draft (Wright + severity + threshold + summary)
#'
#' Builds a 2x2 draft composite for an `mfrm_fit`, suitable for reviewing a
#' possible "Figure 1" in the Rasch-family `RSM`/`PCM` manuscript route. Panels: (1)
#' single-panel FACETS-style Wright map (without CI whiskers), (2)
#' rater severity profile with CI whiskers, (3) threshold ladder, (4)
#' a one-line reliability / separation summary block. Each panel reuses
#' the standalone plot helper so the visual language is consistent
#' with the rest of the package.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output.
#' @param rater_facet Facet name to use as the "rater" axis (default
#'   `"Rater"`).
#' @param ci_level Confidence level for the rater severity panel.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw the composite immediately with
#'   `graphics::layout()`.
#'
#' @return Invisibly, an `mfrm_plot_data` object whose `data` slot
#'   bundles the four panel data objects under `wright`, `severity`,
#'   `threshold`, and `summary`. Fit readiness is retained in
#'   `data$fit_readiness`, `data$interpretation_status`, and
#'   `data$interpretation_note`.
#'
#' @section Interpreting output:
#' Designed as a single-figure Methods or Results draft. The
#' summary panel prints the model class, sample size, log-likelihood, the
#' canonical MML IC panel or an explicit ineligibility/legacy label, and the
#' largest non-Person facet's separation /
#' reliability if available. A fit that has not passed its numerical, data,
#' design, and stability checks produce one warning and a visible
#' `"REVIEW ONLY"` label. Resolve that review before treating the composite as
#' report-ready evidence.
#'
#' @seealso [plot.mfrm_fit()] (`type = "wright"`),
#'   [plot(x, type = "rater_severity")][plot.mfrm_fit], [plot(x, type = "thresholds")][plot.mfrm_fit],
#'   [build_apa_outputs()].
#'
#' @concept confidence intervals
#' @concept reporting workflow
#' @concept visual diagnostics
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_apa_figure_one", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' p <- plot(fit, type = "apa_summary", draw = FALSE)
#' names(p$data)
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_apa_figure_one <- function(fit,
                                diagnostics = NULL,
                                rater_facet = "Rater",
                                ci_level = 0.95,
                                preset = c("standard", "publication", "compact", "monochrome"),
                                draw = TRUE) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  readiness <- .mfrm_fit_plot_readiness(fit)
  muffle_nested_readiness <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        if (identical(conditionMessage(w), readiness$detail)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(suppressWarnings(
      diagnose_mfrm(fit, residual_pca = "none")
    ))
  }
  wright <- muffle_nested_readiness(
    plot(fit, type = "wright", renderer = "facets", show_ci = FALSE,
         preset = preset, draw = FALSE)
  )
  severity <- plot_rater_severity_profile(
    fit, diagnostics = diagnostics, facet = rater_facet,
    ci_level = ci_level, preset = preset, draw = FALSE
  )
  threshold <- plot_threshold_ladder(fit, preset = preset, draw = FALSE)
  summary_lines <- character(0)
  s <- fit$summary
  first_summary_value <- function(tbl, candidates, fallback = NA) {
    if (!is.data.frame(tbl) || nrow(tbl) < 1L) {
      return(fallback)
    }
    for (nm in candidates) {
      if (nm %in% names(tbl)) {
        val <- tbl[[nm]][1]
        if (length(val) == 1L && !is.null(val)) {
          return(val)
        }
      }
    }
    fallback
  }
  if (is.data.frame(s) && nrow(s) >= 1L) {
    n_obs <- first_summary_value(
      s,
      c("N_Obs", "N", "Observations", "TotalObservations"),
      fallback = fit$prep$n_obs %||% NA_integer_
    )
    n_person <- first_summary_value(
      s,
      c("N_Person", "Persons", "NPersons", "PersonN"),
      fallback = fit$prep$n_person %||%
        if (!is.null(fit$facets$person)) nrow(fit$facets$person) else NA_integer_
    )
    summary_lines <- c(
      sprintf("Model: %s | Method: %s",
              as.character(first_summary_value(s, "Model", NA_character_)),
              as.character(first_summary_value(s, "Method", NA_character_))),
      sprintf("N obs = %s | Persons = %s",
              format(n_obs, big.mark = ","),
              format(n_person, big.mark = ",")),
      mfrm_ic_console_lines(s, digits = 2L)
    )
  }
  rel <- diagnostics$reliability
  if (is.data.frame(rel) && nrow(rel) >= 1L) {
    rel_non_person <- rel[as.character(rel$Facet) != "Person", , drop = FALSE]
    if (nrow(rel_non_person) >= 1L) {
      pick <- which.max(suppressWarnings(as.numeric(rel_non_person$Separation)))
      if (length(pick) == 1L && is.finite(pick)) {
        summary_lines <- c(
          summary_lines,
          sprintf("%s separation = %.2f | reliability = %.2f",
                  as.character(rel_non_person$Facet[pick]),
                  suppressWarnings(as.numeric(rel_non_person$Separation[pick])),
                  suppressWarnings(as.numeric(rel_non_person$Reliability[pick])))
        )
      }
    }
  }
  if (length(summary_lines) == 0L) {
    summary_lines <- "Summary table unavailable."
  }

  if (isTRUE(draw)) {
    old_par <- graphics::par()[c("mfrow", "cex", "mex")]
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE))
    muffle_nested_readiness(
      plot(fit, type = "wright", renderer = "facets", show_ci = FALSE,
           preset = preset, draw = TRUE)
    )
    plot_rater_severity_profile(fit, diagnostics = diagnostics,
                                facet = rater_facet,
                                ci_level = ci_level,
                                preset = preset, draw = TRUE)
    plot_threshold_ladder(fit, preset = preset, draw = TRUE)
    apply_plot_preset(style)
    plot.new()
    display_lines <- unlist(lapply(summary_lines, strwrap,
      width = max(25L, floor(graphics::par("pin")[1] * 15))), use.names = FALSE)
    text(
      x = 0, y = seq(0.95, 0.05, length.out = length(display_lines)),
      labels = display_lines, adj = 0, cex = 0.9
    )
  }

  out <- new_mfrm_plot_data(
    "apa_figure_one",
    list(
      data = list(
        wright = wright,
        severity = severity,
        threshold = threshold,
        summary = summary_lines
      ),
      title = "APA Figure 1: Wright + severity + threshold + summary",
      subtitle = "2x2 manuscript-oriented draft composite",
      preset = style$name
    )
  )
  out <- .mfrm_attach_plot_readiness(out, readiness)
  if (!isTRUE(readiness$ready)) {
    warning(readiness$detail, call. = FALSE)
  }
  invisible(out)
}

# ==============================================================================
# Additional visualization helpers (added in 0.1.6)
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
#' bounded GPCM) there is one ladder per `step_facet` level.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param highlight_disorder Logical. When `TRUE` (default), draw
#'   disordered segments with the preset's `fail` colour and add a
#'   subtitle counting the disordered groups.
#' @param preset Visual preset (`"standard"`, `"publication"`, `"compact"`,
#'   or `"monochrome"`).
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object with a `data` slot containing
#'   columns `Group`, `Step`, `Threshold`, `Disordered` for each ladder
#'   row.
#'
#' @section Interpreting output:
#' Within each ladder, thresholds should ascend monotonically. A
#' disordered crossing (highlighted in the fail colour) suggests that
#' the corresponding category is rarely the most likely response over
#' any logit interval, and is a common trigger for category-collapsing
#' decisions.
#'
#' @seealso [category_structure_report()], [category_curves_report()],
#'   [plot.mfrm_fit()] (`type = "ccc"`).
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' p <- plot_threshold_ladder(fit, draw = FALSE)
#' head(p$data$data)
#' @export
plot_threshold_ladder <- function(fit,
                                  highlight_disorder = TRUE,
                                  preset = c("standard", "publication", "compact", "monochrome"),
                                  draw = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  spec <- tryCatch(build_step_curve_spec(fit), error = function(e) NULL)
  if (is.null(spec) || length(spec$groups) == 0L) {
    stop("Step estimates are not available for this fit.", call. = FALSE)
  }
  rows <- do.call(rbind, lapply(names(spec$groups), function(grp) {
    g <- spec$groups[[grp]]
    tau <- as.numeric(g$tau)
    if (length(tau) == 0L) return(NULL)
    disordered <- c(FALSE, diff(tau) < 0)
    data.frame(
      Group = grp,
      Step = paste0("tau_", seq_along(tau)),
      Threshold = tau,
      Disordered = disordered,
      stringsAsFactors = FALSE
    )
  }))
  rownames(rows) <- NULL
  if (is.null(rows) || nrow(rows) == 0L) {
    stop("No threshold rows to plot.", call. = FALSE)
  }
  n_disorder_groups <- length(unique(rows$Group[rows$Disordered]))
  plot_title <- "Threshold ladder"
  plot_subtitle <- if (isTRUE(highlight_disorder) && n_disorder_groups > 0L) {
    sprintf("%d ladder(s) with at least one disordered step",
            n_disorder_groups)
  } else {
    "Adjacent-category Rasch thresholds (logit)"
  }

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    groups <- unique(rows$Group)
    n_grp <- length(groups)
    x_pos <- seq_along(groups)
    y_lim <- range(rows$Threshold, finite = TRUE) +
      c(-0.15, 0.15) * diff(range(rows$Threshold, finite = TRUE))
    graphics::plot(
      x = numeric(0), y = numeric(0),
      xlim = c(0.5, n_grp + 0.5), ylim = y_lim,
      xaxt = "n", xlab = "", ylab = "Threshold (logit)",
      main = plot_title
    )
    graphics::title(sub = plot_subtitle, line = 2.2, cex.sub = 0.9)
    graphics::axis(1, at = x_pos, labels = groups, las = 2, cex.axis = 0.85)
    graphics::abline(h = 0, lty = 2, col = style$neutral)
    for (i in seq_along(groups)) {
      sub <- rows[rows$Group == groups[i], , drop = FALSE]
      sub <- sub[order(seq_len(nrow(sub))), , drop = FALSE]
      tau_vec <- sub$Threshold
      x <- x_pos[i]
      if (length(tau_vec) >= 2L) {
        for (k in seq_len(length(tau_vec) - 1L)) {
          col_seg <- if (isTRUE(highlight_disorder) && tau_vec[k + 1L] < tau_vec[k]) {
            style$fail
          } else {
            style$accent_primary
          }
          graphics::segments(x, tau_vec[k], x, tau_vec[k + 1L],
                             col = col_seg, lwd = 2)
        }
      }
      graphics::points(rep(x, length(tau_vec)), tau_vec,
                       pch = 19, col = style$accent_primary)
      graphics::text(x + 0.15, tau_vec, labels = sub$Step,
                     cex = 0.7, adj = 0)
    }
  }

  out <- new_mfrm_plot_data(
    "threshold_ladder",
    list(
      data = rows,
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
#' fit-acceptable; amber and red corners are candidates for misfit
#' review (overfit / underfit) using
#' [unexpected_response_table()] for follow-up.
#'
#' @seealso [diagnose_mfrm()], [unexpected_response_table()],
#'   [build_misfit_casebook()].
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' p <- plot_person_fit(fit, draw = FALSE)
#' head(p$data$data)
#' @export
plot_person_fit <- function(fit,
                            diagnostics = NULL,
                            lower = 0.5,
                            upper = 1.5,
                            top_n_label = 12L,
                            preset = c("standard", "publication", "compact", "monochrome"),
                            draw = TRUE,
                            fit_index = c("meansquare", "loglik")) {
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
  m <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  if (!all(c("Facet", "Level", "Infit", "Outfit") %in% names(m))) {
    stop("Diagnostics measures are missing required columns.",
         call. = FALSE)
  }
  m <- m[as.character(m$Facet) == "Person", , drop = FALSE]
  m <- m[is.finite(m$Infit) & is.finite(m$Outfit), , drop = FALSE]
  if (nrow(m) == 0L) {
    stop("No finite Infit/Outfit values for any Person row.",
         call. = FALSE)
  }
  n_col <- if ("N" %in% names(m)) "N" else if ("N.x" %in% names(m)) "N.x" else NA_character_
  m$N <- if (!is.na(n_col)) suppressWarnings(as.numeric(m[[n_col]])) else 1
  in_band <- m$Infit >= lower & m$Infit <= upper
  out_band <- m$Outfit >= lower & m$Outfit <= upper
  m$Status <- ifelse(in_band & out_band, "in_band",
                     ifelse(in_band | out_band, "one_outside", "both_outside"))
  status_color <- c(in_band = style$success,
                    one_outside = style$warn,
                    both_outside = style$fail)
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
    merge_tbl <- person_fit_indices[, setdiff(names(person_fit_indices), "N"), drop = FALSE]
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
        FALSE
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
      "Infit and Outfit per person (acceptance band [%g, %g], Linacre 2002)",
      lower, upper
    )
  } else {
    "Report index per person: lz* when available, otherwise lz with caveat"
  }

  z_5pct <- stats::qnorm(0.975)
  z_1pct <- stats::qnorm(0.995)

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    if (identical(fit_index, "meansquare")) {
      cex_size <- 0.6 + 1.6 * sqrt(m$N / max(m$N, na.rm = TRUE))
      x_rng <- range(c(m$Infit, lower, upper), finite = TRUE)
      y_rng <- range(c(m$Outfit, lower, upper), finite = TRUE)
      x_pad <- if (diff(x_rng) > 0) c(-0.05, 0.05) * diff(x_rng) else c(-0.1, 0.1)
      y_pad <- if (diff(y_rng) > 0) c(-0.05, 0.05) * diff(y_rng) else c(-0.1, 0.1)
      graphics::plot(
        x = m$Infit, y = m$Outfit,
        xlab = "Infit MnSq", ylab = "Outfit MnSq",
        main = plot_title,
        pch = 21, bg = m$Color, col = "white",
        cex = cex_size,
        xlim = x_rng + x_pad, ylim = y_rng + y_pad
      )
      graphics::title(sub = plot_subtitle, line = 2.2, cex.sub = 0.9)
      graphics::abline(h = 1, v = 1, lty = 3, col = style$neutral)
      graphics::abline(h = c(lower, upper), v = c(lower, upper),
                       lty = 2, col = style$grid)
      n_lbl <- min(top_n_label, nrow(m))
      if (n_lbl > 0L) {
        graphics::text(
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
        graphics::plot.new()
        graphics::title(main = plot_title, sub = "No finite lz/lz* report index available")
      } else {
        y <- seq_len(nrow(ranked))
        x_rng <- range(c(ranked$ReportValue, -z_1pct, -z_5pct, z_5pct, z_1pct),
                       finite = TRUE)
        x_pad <- if (diff(x_rng) > 0) c(-0.08, 0.08) * diff(x_rng) else c(-0.5, 0.5)
        graphics::plot(
          x = ranked$ReportValue, y = y,
          xlab = "Person-fit report index",
          ylab = "Persons ordered by report index",
          yaxt = "n",
          main = plot_title,
          pch = 21, bg = ranked$ReviewColor, col = "white",
          cex = 0.85 + 1.2 * sqrt(ranked$N / max(ranked$N, na.rm = TRUE)),
          xlim = x_rng + x_pad
        )
        graphics::title(sub = plot_subtitle, line = 2.2, cex.sub = 0.9)
        graphics::abline(v = 0, lty = 3, col = style$neutral)
        graphics::abline(v = c(-z_5pct, z_5pct), lty = 2, col = style$grid)
        graphics::abline(v = c(-z_1pct, z_1pct), lty = 3, col = style$grid)
        label_order <- order(abs(ranked$ReportValue), decreasing = TRUE)
        n_lbl <- min(top_n_label, length(label_order))
        if (n_lbl > 0L) {
          idx <- label_order[seq_len(n_lbl)]
          graphics::text(
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
    FlaggedRows = c(
      sum(payload$Status != "in_band", na.rm = TRUE),
      sum(payload$ReportFlag %in% TRUE, na.rm = TRUE)
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
      label = c("-1% threshold", "-5% threshold", "5% threshold", "1% threshold"),
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
        label = c("In band", "One outside", "Both outside",
                  "Not flagged", "Review 5%", "Review 1%"),
        role = rep("status", 6),
        aesthetic = rep("point", 6),
        value = c(unname(status_color),
                  unname(review_color[c("not_flagged", "review_5pct", "review_1pct")]))
      ),
      reference_lines = reference_lines,
      preset = style$name
    )
  )
  invisible(out)
}


#' Plot per-rater severity ranking with confidence interval whiskers
#'
#' Ranks the levels of a chosen rater facet by estimated severity and
#' draws each level as a horizontal CI whisker around the point
#' estimate. Optional gentle / strict guidance bands at `+/-0.5` and
#' `+/-1.0` logit relative to the centred mean make rater calibration
#' easy to read for training feedback.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output. When omitted,
#'   `diagnose_mfrm(fit, residual_pca = "none")` is run internally.
#' @param facet Facet name to plot (default `"Rater"`). Any non-Person
#'   facet name is accepted.
#' @param ci_level Confidence level used for the whiskers (default
#'   `0.95`). Bounds use `+/- z * ModelSE`.
#' @param show_bands Logical. When `TRUE` (default) draw shaded
#'   `+/-0.5` (gentle) and `+/-1.0` (strict) logit guidance bands.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object whose `data` slot contains
#'   columns `Level`, `Estimate`, `SE`, `CI_Lower`, `CI_Upper`,
#'   `Band`.
#'
#' @section Interpreting output:
#' The vertical reference line at zero is the sum-to-zero centring
#' point. Levels well within `+/- 0.5 logit` (gentle band) are
#' typically interchangeable in operational scoring; levels outside
#' `+/- 1.0 logit` (strict band) deserve targeted training or
#' anchoring.
#'
#' @seealso [diagnose_mfrm()], [analyze_facet_equivalence()],
#'   [plot_facet_equivalence()].
#'
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept rater severity
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' p <- plot_rater_severity_profile(fit, draw = FALSE)
#' head(p$data$data)
#' @export
plot_rater_severity_profile <- function(fit,
                                        diagnostics = NULL,
                                        facet = "Rater",
                                        ci_level = 0.95,
                                        show_bands = TRUE,
                                        preset = c("standard", "publication", "compact", "monochrome"),
                                        draw = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(suppressWarnings(
      diagnose_mfrm(fit, residual_pca = "none")
    ))
  }
  m <- as.data.frame(diagnostics$measures, stringsAsFactors = FALSE)
  m <- m[as.character(m$Facet) == as.character(facet), , drop = FALSE]
  if (nrow(m) == 0L) {
    stop(sprintf("No rows for facet '%s' in diagnostics$measures.", facet),
         call. = FALSE)
  }
  m$Estimate <- suppressWarnings(as.numeric(m$Estimate))
  m$SE <- suppressWarnings(as.numeric(m$ModelSE %||% m$SE))
  m <- m[is.finite(m$Estimate), , drop = FALSE]
  z_ci <- stats::qnorm(1 - (1 - ci_level) / 2)
  m$CI_Lower <- m$Estimate - z_ci * m$SE
  m$CI_Upper <- m$Estimate + z_ci * m$SE
  m$Band <- ifelse(abs(m$Estimate) <= 0.5, "gentle",
                   ifelse(abs(m$Estimate) <= 1.0, "moderate", "strict"))
  m <- m[order(m$Estimate), , drop = FALSE]
  plot_title <- sprintf("%s severity profile", facet)
  plot_subtitle <- sprintf(
    "Sum-to-zero centred; +/-0.5 gentle band, +/-1.0 strict band; %g%% CI from ModelSE",
    round(100 * ci_level)
  )

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    y <- seq_len(nrow(m))
    xrange <- range(c(m$CI_Lower, m$CI_Upper, -1.05, 1.05), finite = TRUE,
                    na.rm = TRUE)
    graphics::plot(
      x = m$Estimate, y = y,
      type = "n",
      xlim = xrange,
      yaxt = "n",
      xlab = "Severity (logit)",
      ylab = "",
      main = plot_title
    )
    graphics::title(sub = plot_subtitle, line = 2.2, cex.sub = 0.9)
    if (isTRUE(show_bands)) {
      usr <- graphics::par("usr")
      graphics::rect(-0.5, usr[3], 0.5, usr[4], border = NA,
                     col = grDevices::adjustcolor(style$success, alpha.f = 0.10))
      graphics::rect(-1.0, usr[3], -0.5, usr[4], border = NA,
                     col = grDevices::adjustcolor(style$warn, alpha.f = 0.08))
      graphics::rect(0.5, usr[3], 1.0, usr[4], border = NA,
                     col = grDevices::adjustcolor(style$warn, alpha.f = 0.08))
    }
    graphics::abline(v = 0, lty = 2, col = style$neutral)
    valid <- is.finite(m$CI_Lower) & is.finite(m$CI_Upper)
    if (any(valid)) {
      graphics::segments(m$CI_Lower[valid], y[valid],
                         m$CI_Upper[valid], y[valid],
                         col = style$accent_primary, lwd = 2)
    }
    graphics::points(m$Estimate, y, pch = 19, col = style$accent_primary)
    graphics::axis(2, at = y, labels = as.character(m$Level), las = 1,
                   cex.axis = 0.85)
  }

  payload <- data.frame(
    Level = as.character(m$Level),
    Estimate = m$Estimate,
    SE = m$SE,
    CI_Lower = m$CI_Lower,
    CI_Upper = m$CI_Upper,
    Band = m$Band,
    stringsAsFactors = FALSE
  )
  out <- new_mfrm_plot_data(
    "rater_severity_profile",
    list(
      data = payload,
      facet = facet,
      ci_level = ci_level,
      title = plot_title,
      subtitle = plot_subtitle,
      legend = new_plot_legend(
        label = c("Estimate", "CI whisker", "+/-0.5 gentle", "+/-1.0 strict"),
        role = c("location", "uncertainty", "band", "band"),
        aesthetic = c("point", "segment", "fill", "fill"),
        value = c(style$accent_primary, style$accent_primary,
                  style$success, style$warn)
      ),
      reference_lines = new_reference_lines(
        "v", 0, "Sum-to-zero centred", "dashed", "reference"
      ),
      preset = style$name
    )
  )
  invisible(out)
}


#' Summary plot of differential functioning effect sizes
#'
#' Compact effect-size summary for a [analyze_dff()] / [analyze_dif()]
#' result. Shows each contrast's signed effect size as a horizontal bar
#' with a vertical reference at zero, coloured by the method-appropriate
#' classification. ETS-style A / B / C colours are used only when they
#' are actually available; residual-method screening labels otherwise use
#' the neutral colour.
#'
#' @param x Output from [analyze_dff()] or [analyze_dif()].
#' @param top_n Maximum rows shown (default `30`).
#' @param sort_by `"abs_effect"` (default), `"effect"`, or
#'   `"classification"`.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#' @param ci_level Optional confidence level for approximate normal
#'   intervals drawn from `Effect +/- z * SE` when finite standard errors are
#'   available. Use `NULL` (default) to omit intervals.
#' @param effect_thresholds Optional numeric vector of absolute effect-size
#'   guide lines to draw at `+/- threshold`. These are display aids; only use
#'   ETS-like values when the source rows support ETS interpretation.
#' @param effect_axis_label Optional x-axis label override. When `NULL`, the
#'   label is chosen from the DFF method.
#'
#' @return An `mfrm_plot_data` object whose `data` slot contains
#'   columns `Pair`, `Effect`, `SE`, `Classification`, `Color`.
#'
#' @section Interpreting output:
#' Bars are anchored at zero. Width corresponds to effect size on the
#' contrast's native scale. For `method = "residual"`, this is the
#' observed-minus-expected average screening contrast between groups. For
#' `method = "refit"`, this is the subgroup parameter difference on the
#' fitted logit scale when linking support allows a comparable contrast.
#' The ETS classification (A negligible, B moderate, C large) drives bar
#' colour only when `ClassificationSystem == "ETS"`; otherwise the bar
#' uses the preset's neutral.
#'
#' @seealso [analyze_dff()], [analyze_dif()], [plot_dif_heatmap()].
#'
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept DFF DIF
#'
#' @examples
#' toy <- load_mfrmr_data("example_bias")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' dff <- analyze_dff(fit, diagnostics = diag,
#'                    facet = "Rater", group = "Group", data = toy)
#' unique(dff$dif_table$ClassificationSystem)
#' p <- plot_dif_summary(dff, draw = FALSE)
#' head(p$data$data)
#' @export
plot_dif_summary <- function(x,
                             top_n = 30L,
                             sort_by = c("abs_effect", "effect", "classification"),
                             preset = c("standard", "publication", "compact", "monochrome"),
                             draw = TRUE,
                             ci_level = NULL,
                             effect_thresholds = NULL,
                             effect_axis_label = NULL) {
  if (!inherits(x, c("mfrm_dff", "mfrm_dif"))) {
    stop("`x` must be output from analyze_dff() or analyze_dif().",
         call. = FALSE)
  }
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
  if (!is.null(ci_level)) {
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    finite_se <- is.finite(tbl$SE)
    tbl$CI_Lower <- ifelse(finite_se, tbl$Effect - z * tbl$SE, NA_real_)
    tbl$CI_Upper <- ifelse(finite_se, tbl$Effect + z * tbl$SE, NA_real_)
  } else {
    tbl$CI_Lower <- NA_real_
    tbl$CI_Upper <- NA_real_
  }
  plot_title <- "Differential functioning summary"
  plot_subtitle <- sprintf(
    "%d row(s) shown; sorted by %s",
    nrow(tbl), sort_by
  )

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    y <- rev(seq_len(nrow(tbl)))
    x_values <- c(tbl$Effect, tbl$CI_Lower, tbl$CI_Upper, 0,
                  effect_thresholds, -effect_thresholds)
    x_range <- range(x_values, finite = TRUE)
    if (!all(is.finite(x_range)) || diff(x_range) == 0) {
      x_range <- x_range[1] + c(-1, 1)
    }
    xlim <- x_range + c(-0.05, 0.05) * diff(x_range)
    graphics::plot(
      x = tbl$Effect, y = y, type = "n",
      xlim = xlim, yaxt = "n",
      xlab = axis_label, ylab = "",
      main = plot_title
    )
    graphics::title(sub = plot_subtitle, line = 2.2, cex.sub = 0.9)
    graphics::abline(v = 0, lty = 2, col = style$neutral)
    if (length(effect_thresholds) > 0L) {
      for (thr in effect_thresholds) {
        graphics::abline(v = c(-thr, thr), lty = 3, col = style$warn)
      }
    }
    if (!is.null(ci_level) && any(is.finite(tbl$CI_Lower) & is.finite(tbl$CI_Upper))) {
      graphics::segments(tbl$CI_Lower, y, tbl$CI_Upper, y,
                         col = style$neutral, lwd = 1)
    }
    graphics::segments(0, y, tbl$Effect, y,
                       col = tbl$Color, lwd = 5)
    graphics::points(tbl$Effect, y, pch = 19, col = tbl$Color)
    graphics::axis(2, at = y, labels = tbl$Pair, las = 1, cex.axis = 0.7)
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
      sort_by = sort_by,
      title = plot_title,
      subtitle = plot_subtitle,
      legend = new_plot_legend(
        label = c("A: negligible", "B: moderate", "C: large", "unclassified"),
        role = rep("classification", 4),
        aesthetic = rep("bar", 4),
        value = c(style$success, style$warn, style$fail, style$neutral)
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
  invisible(out)
}


#' Manuscript-ready four-panel composite (Wright + severity + threshold + summary)
#'
#' Builds a 2x2 publication composite for an `mfrm_fit`, suitable for a
#' "Figure 1" in the Rasch-family `RSM`/`PCM` manuscript route. Panels: (1)
#' Wright map, (2)
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
#'   `threshold`, and `summary`.
#'
#' @section Interpreting output:
#' Designed for a single-figure Methods or Results overview. The
#' summary panel prints the model class, sample size, log-likelihood,
#' AIC/BIC, and the largest non-Person facet's separation /
#' reliability if available.
#'
#' @seealso [plot.mfrm_fit()] (`type = "wright"`),
#'   [plot_rater_severity_profile()], [plot_threshold_ladder()],
#'   [build_apa_outputs()].
#'
#' @concept confidence intervals
#' @concept reporting workflow
#' @concept visual diagnostics
#'
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                 method = "JML", maxit = 30)
#' p <- plot_apa_figure_one(fit, draw = FALSE)
#' names(p$data)
#' @export
plot_apa_figure_one <- function(fit,
                                diagnostics = NULL,
                                rater_facet = "Rater",
                                ci_level = 0.95,
                                preset = c("standard", "publication", "compact", "monochrome"),
                                draw = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(suppressWarnings(
      diagnose_mfrm(fit, residual_pca = "none")
    ))
  }
  wright <- plot(fit, type = "wright", draw = FALSE)
  severity <- plot_rater_severity_profile(
    fit, diagnostics = diagnostics, facet = rater_facet,
    ci_level = ci_level, draw = FALSE
  )
  threshold <- plot_threshold_ladder(fit, draw = FALSE)
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
    loglik <- suppressWarnings(as.numeric(first_summary_value(s, "LogLik", NA_real_)))
    aic <- suppressWarnings(as.numeric(first_summary_value(s, "AIC", NA_real_)))
    bic <- suppressWarnings(as.numeric(first_summary_value(s, "BIC", NA_real_)))
    summary_lines <- c(
      sprintf("Model: %s | Method: %s",
              as.character(first_summary_value(s, "Model", NA_character_)),
              as.character(first_summary_value(s, "Method", NA_character_))),
      sprintf("N obs = %s | Persons = %s",
              format(n_obs, big.mark = ","),
              format(n_person, big.mark = ",")),
      sprintf("LogLik = %.2f | AIC = %.2f | BIC = %.2f",
              loglik,
              aic,
              bic)
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
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE))
    plot(fit, type = "wright", preset = preset, draw = TRUE)
    plot_rater_severity_profile(fit, diagnostics = diagnostics,
                                facet = rater_facet,
                                ci_level = ci_level,
                                preset = preset, draw = TRUE)
    plot_threshold_ladder(fit, preset = preset, draw = TRUE)
    apply_plot_preset(style)
    graphics::plot.new()
    graphics::title(main = "Fit summary", line = 1)
    graphics::text(
      x = 0, y = seq(1, 0, length.out = length(summary_lines) + 2L)[-c(1, length(summary_lines) + 2L)],
      labels = summary_lines, adj = 0, cex = 0.9
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
      subtitle = "2x2 publication composite",
      preset = style$name
    )
  )
  invisible(out)
}

# ==============================================================================
# Screening / case-level visualization helpers
# ==============================================================================
#
# Four helpers grouped here for case-level screening: a Guttman scalogram,
# a Q-Q plot of person residual aggregates, a per-rater severity
# trajectory across waves, and a pairwise rater-agreement heatmap. They
# all follow the mfrmr plot contract: accept an `mfrm_fit` (or a related
# class), resolve a preset via resolve_plot_preset() + apply_plot_preset(),
# and return an `mfrm_plot_data` object that downstream code can
# re-render or export.
# ==============================================================================


#' Guttman-style scalogram of person x item observed responses
#'
#' Draws a person x item (or person x facet-level) matrix coloured by
#' observed category, with rows ordered by person measure and columns
#' ordered by location measure. Unexpected responses (those that fall
#' far from the expected category at a given theta) are highlighted
#' with a heavy border so the visual reads as a Rasch-convention
#' Guttman scalogram.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output; used to pick
#'   up unexpected-response flags when available.
#' @param column_facet Facet name used for the columns. Default
#'   `"Criterion"` when the fit contains it, otherwise the last entry
#'   of `fit$config$facet_names`.
#' @param top_n_persons Maximum number of persons shown (default
#'   `40`). Persons closest to the median measure are retained when
#'   the population exceeds this cap.
#' @param highlight_unexpected Logical. When `TRUE` (default), draw a
#'   heavy border around cells flagged as unexpected by
#'   [unexpected_response_table()].
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object whose `data` slot bundles the
#'   scalogram matrix and the optional unexpected-response overlay.
#'
#' @seealso [unexpected_response_table()] for the case-level review of
#'   the cells flagged in the overlay;
#'   [plot(x, type = "agreement", view = "heatmap")][plot.mfrm_fit] for a complementary rater-pair
#'   view of the same residual structure;
#'   [diagnose_mfrm()] for the underlying diagnostics bundle.
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_guttman_scalogram", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' p <- plot(fit, type = "guttman", draw = FALSE)
#' dim(p$data$matrix)
#' # Look for: a clean monotone "staircase" of higher scores in the
#' #   upper-right triangle and lower scores in the lower-left, once
#' #   rows are sorted by person ability. Cells circled by the
#' #   unexpected-response overlay break the staircase and warrant
#' #   case-level review with `unexpected_response_table()`.
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_guttman_scalogram <- function(fit,
                                   diagnostics = NULL,
                                   column_facet = NULL,
                                   top_n_persons = 40L,
                                   highlight_unexpected = TRUE,
                                   preset = c("standard", "publication", "compact", "monochrome"),
                                   draw = TRUE) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  top_n_persons <- max(2L, as.integer(top_n_persons))

  obs <- as.data.frame(fit$prep$data %||% NULL, stringsAsFactors = FALSE)
  if (nrow(obs) == 0L) {
    stop("Prepared observation data are not available on this fit.",
         call. = FALSE)
  }
  if (is.null(column_facet)) {
    facet_names <- as.character(fit$config$facet_names %||% character())
    column_facet <- if ("Criterion" %in% facet_names) "Criterion" else
      utils::tail(facet_names, 1L)
  }
  if (!column_facet %in% names(obs)) {
    stop("Column facet '", column_facet,
         "' not found in prepared data.", call. = FALSE)
  }
  column_facet <- as.character(column_facet)

  # Collapse observations to (Person, column_facet) -> mean Score when
  # multiple rows exist for the same cell.
  cell <- stats::aggregate(
    obs$Score, by = list(Person = obs$Person,
                          ColVar = obs[[column_facet]]),
    FUN = function(v) round(mean(v, na.rm = TRUE))
  )
  names(cell)[names(cell) == "x"] <- "Score"
  persons <- as.character(unique(cell$Person))
  cols <- as.character(unique(cell$ColVar))

  # Order rows by person measure.
  person_tbl <- as.data.frame(fit$facets$person, stringsAsFactors = FALSE)
  theta_lookup <- stats::setNames(
    suppressWarnings(as.numeric(person_tbl$Estimate)),
    as.character(person_tbl$Person)
  )
  persons <- persons[order(theta_lookup[persons], na.last = TRUE)]

  # Order columns by facet-level measure.
  others <- as.data.frame(fit$facets$others, stringsAsFactors = FALSE)
  col_sub <- others[as.character(others$Facet) == column_facet, , drop = FALSE]
  col_lookup <- stats::setNames(
    suppressWarnings(as.numeric(col_sub$Estimate)),
    as.character(col_sub$Level)
  )
  cols <- cols[order(col_lookup[cols], na.last = TRUE)]

  # Cap persons to top_n_persons around the median.
  if (length(persons) > top_n_persons) {
    mid <- floor(length(persons) / 2L)
    keep <- seq(
      max(1L, mid - top_n_persons %/% 2L),
      min(length(persons), mid - top_n_persons %/% 2L + top_n_persons - 1L)
    )
    persons <- persons[keep]
  }
  cell <- cell[as.character(cell$Person) %in% persons &
                  as.character(cell$ColVar) %in% cols, , drop = FALSE]
  mat <- matrix(NA_integer_, nrow = length(persons), ncol = length(cols),
                dimnames = list(persons, cols))
  for (k in seq_len(nrow(cell))) {
    mat[as.character(cell$Person[k]),
        as.character(cell$ColVar[k])] <- as.integer(cell$Score[k])
  }

  # Optional unexpected overlay.
  unexpected_overlay <- tryCatch({
    if (isTRUE(highlight_unexpected)) {
      u_tbl <- unexpected_response_table(fit, diagnostics = diagnostics)
      u_tbl <- as.data.frame(u_tbl$table %||% u_tbl, stringsAsFactors = FALSE)
      if (nrow(u_tbl) > 0L && all(c("Person", column_facet) %in% names(u_tbl))) {
        data.frame(
          Person = as.character(u_tbl$Person),
          Column = as.character(u_tbl[[column_facet]]),
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(Person = character(0), Column = character(0),
                   stringsAsFactors = FALSE)
      }
    } else {
      data.frame(Person = character(0), Column = character(0),
                 stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    data.frame(Person = character(0), Column = character(0),
               stringsAsFactors = FALSE)
  })

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    n_cat <- length(unique(stats::na.omit(as.vector(mat))))
    cat_palette <- grDevices::hcl.colors(max(3L, n_cat), "YlGnBu", rev = TRUE)
    image(
      x = seq_len(ncol(mat)), y = seq_len(nrow(mat)),
      z = t(mat), col = cat_palette,
      xaxt = "n", yaxt = "n",
      xlab = column_facet, ylab = "Person (sorted by measure)",
      main = ""
    )
    axis(1, at = seq_len(ncol(mat)), labels = cols,
                   las = 2, cex.axis = 0.75)
    if (nrow(mat) <= 25L) {
      axis(2, at = seq_len(nrow(mat)), labels = rownames(mat),
                     las = 1, cex.axis = 0.7)
    }
    if (nrow(unexpected_overlay) > 0L) {
      for (k in seq_len(nrow(unexpected_overlay))) {
        r <- match(unexpected_overlay$Person[k], rownames(mat))
        c <- match(unexpected_overlay$Column[k], colnames(mat))
        if (is.finite(r) && is.finite(c)) {
          rect(c - 0.5, r - 0.5, c + 0.5, r + 0.5,
                         border = style$fail, lwd = 2)
        }
      }
    }
  }

  out <- new_mfrm_plot_data(
    "guttman_scalogram",
    list(
      matrix = mat,
      unexpected = unexpected_overlay,
      persons = persons,
      column_facet = column_facet,
      title = "Guttman scalogram",
      subtitle = sprintf(
        "%d person(s) x %d %s level(s); unexpected cells highlighted",
        length(persons), length(cols), column_facet
      ),
      preset = style$name
    )
  )
  invisible(out)
}


#' Normal quantile-quantile plot of person standardized residuals
#'
#' Produces a Q-Q plot of per-person standardized residuals. Under the
#' fitted Rasch-family model the residuals are approximately N(0, 1),
#' so deviations from the reference line diagnose distributional
#' misfit that mean-square summaries may miss.
#'
#' @param fit An `mfrm_fit`.
#' @param diagnostics Optional [diagnose_mfrm()] output; required
#'   entries are generated internally when absent.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object with a `data` slot containing
#'   `Person`, `Theoretical`, `Sample` columns.
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_residual_qq", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' p <- plot(fit, type = "residual_qq", draw = FALSE)
#' head(p$data$data)
#' # Look for: points hugging the y = x reference line. Heavy upper-
#' #   right tails indicate persons whose residual aggregates exceed
#' #   the standard normal expectation; pair with `plot()`
#' #   for case-level follow-up. This is an exploratory screen; do
#' #   not treat tail behaviour as a definitive normality test.
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_residual_qq <- function(fit,
                             diagnostics = NULL,
                             preset = c("standard", "publication", "compact", "monochrome"),
                             draw = TRUE) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  obs <- if (is.null(diagnostics)) {
    compute_obs_table(.mfrm_validate_diagnostic_fit(fit))
  } else diagnostics$obs
  obs <- as.data.frame(obs %||% NULL, stringsAsFactors = FALSE)
  if (nrow(obs) == 0L || !"StdResidual" %in% names(obs)) {
    stop("Standardized residuals are not available on this fit.",
         call. = FALSE)
  }
  obs <- obs[is.finite(obs$StdResidual), , drop = FALSE]
  person_sum <- stats::aggregate(obs$StdResidual,
                                  by = list(Person = obs$Person),
                                  FUN = function(v)
                                    sum(v, na.rm = TRUE) /
                                      sqrt(max(1L, length(v))))
  names(person_sum)[2] <- "StdResidSum"
  person_sum <- person_sum[is.finite(person_sum$StdResidSum), ,
                            drop = FALSE]
  q_sample <- sort(person_sum$StdResidSum)
  q_theory <- stats::qnorm(stats::ppoints(length(q_sample)))
  payload <- data.frame(
    Person = as.character(person_sum$Person[order(person_sum$StdResidSum)]),
    Theoretical = q_theory,
    Sample = q_sample,
    stringsAsFactors = FALSE
  )

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    plot(
      x = q_theory, y = q_sample,
      xlab = "Theoretical quantile (N(0, 1))",
      ylab = "Sample quantile (person standardized residual)",
      main = "",
      pch = 19, col = style$accent_primary
    )
    abline(a = 0, b = 1, lty = 2, col = style$neutral)
  }

  out <- new_mfrm_plot_data(
    "residual_qq",
    list(
      data = payload,
      title = "Normal Q-Q of person residuals",
      subtitle = sprintf("%d person(s) in the Q-Q plot", nrow(payload)),
      preset = style$name,
      reference_lines = new_reference_lines(
        "diag", 1, "Identity (N(0,1))", "dashed", "reference"
      )
    )
  )
  invisible(out)
}


#' Rater-severity trajectory across an ordered wave / occasion variable
#'
#' Plots each rater's severity estimate across a user-supplied
#' ordering variable (e.g. `Session`, `Wave`, `AdminDate`), producing
#' one line per rater. When the ordering column is time-like (numeric
#' or date), the x-axis is drawn on that scale; otherwise the values
#' are rendered as discrete ordered categories. Useful for rater
#' training / drift feedback loops.
#'
#' @section Anchor-linking caveat:
#' Each wave is fit independently under its own sum-to-zero
#' identification, so the per-wave severity logits live on separate
#' scales unless you actively link them. Before interpreting movement
#' across waves as rater drift, link the waves by either (i) holding
#' common anchors fixed across fits (see
#' [mfrmr_linking_and_dff] for the supported linking route), or
#' (ii) harmonizing the scale post-hoc with a Stocking-Lord type
#' transformation and reviewing the result via [plot(x)][plot.mfrm_anchor_drift].
#' The trajectory plot itself does not perform linking; it only
#' visualizes the supplied fits on their as-fit scales.
#'
#' @param fits A named list of `mfrm_fit` objects, one per wave. Names
#'   become the x-axis labels in their supplied order. Fits are
#'   assumed to have been placed on a common scale via anchor-linking
#'   or an equivalent post-hoc transformation (see the caveat above).
#' @param facet Facet whose levels are tracked (default `"Rater"`).
#' @param ci_level Confidence level for the per-wave CI ribbons
#'   drawn around each trajectory (default `0.95`).
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object whose `data` slot is a long
#'   data.frame with `Wave`, `Level`, `Estimate`, `SE`, `CI_Lower`,
#'   `CI_Upper` columns.
#'
#' @seealso [plot(x)][plot.mfrm_anchor_drift], [mfrmr_linking_and_dff]
#'
#' @concept confidence intervals
#' @concept visual diagnostics
#' @concept linking
#'
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("plot_rater_trajectory", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' fit_a <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                   method = "JML", maxit = 300)
#' fit_b <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                   method = "JML", maxit = 300)
#' p <- plot(fit_a, type = "trajectory", comparison = list(T2 = fit_b),
#'   labels = c("T1", "T2"), draw = FALSE)
#' head(p$data$data)
#' # Look for: stable trajectories (small wave-to-wave shifts within
#' #   each rater's CI ribbon) once the waves are anchor-linked. A
#' #   rater whose line drifts >0.5 logits across waves is the typical
#' #   "calibration drift" signal. Without anchor linking the per-wave
#' #   logits are on different scales and the picture cannot be read
#' #   as drift; see the Anchor-linking caveat in the docstring.
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_rater_trajectory <- function(fits,
                                  facet = "Rater",
                                  ci_level = 0.95,
                                  preset = c("standard", "publication", "compact", "monochrome"),
                                  draw = TRUE) {
  for (fit in fits) stop_if_product_slopes(fit, "plot((fits)[[1L]], type = \"trajectory\", comparison = (fits)[-1L], labels = names(fits))")

  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!is.list(fits) || length(fits) < 2L) {
    stop("`fits` must be a named list of at least two mfrm_fit objects.",
         call. = FALSE)
  }
  if (is.null(names(fits)) || any(!nzchar(names(fits)))) {
    stop("`fits` must have distinct non-empty names.", call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  z_ci <- stats::qnorm((1 - ci_level) / 2, lower.tail = FALSE)
  rows <- list()
  for (w in names(fits)) {
    fit <- fits[[w]]
    if (!inherits(fit, "mfrm_fit")) next
    others <- as.data.frame(fit$facets$others, stringsAsFactors = FALSE)
    sub <- others[as.character(others$Facet) == facet, , drop = FALSE]
    if (nrow(sub) == 0L) next
    # Per-level SE: prefer attached SE, else approximate via diagnose.
    se_vec <- if ("SE" %in% names(sub)) {
      suppressWarnings(as.numeric(sub$SE))
    } else {
      rep(NA_real_, nrow(sub))
    }
    if (all(!is.finite(se_vec))) {
      diag <- tryCatch(
        suppressMessages(suppressWarnings(
          diagnose_mfrm(fit, residual_pca = "none",
                        diagnostic_mode = "legacy")
        )),
        error = function(e) NULL
      )
      if (!is.null(diag) && !is.null(diag$measures)) {
        m <- as.data.frame(diag$measures, stringsAsFactors = FALSE)
        m <- m[as.character(m$Facet) == facet, c("Level", "ModelSE"),
               drop = FALSE]
        se_vec <- suppressWarnings(as.numeric(
          m$ModelSE[match(as.character(sub$Level), as.character(m$Level))]
        ))
      }
    }
    est <- suppressWarnings(as.numeric(sub$Estimate))
    rows[[w]] <- data.frame(
      Wave = w,
      Level = as.character(sub$Level),
      Estimate = est,
      SE = se_vec,
      CI_Lower = est - z_ci * se_vec,
      CI_Upper = est + z_ci * se_vec,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    stop("No fits contained the requested facet '", facet, "'.",
         call. = FALSE)
  }
  long <- do.call(rbind, rows)
  long$WaveOrder <- match(long$Wave, names(fits))
  long <- long[order(long$Level, long$WaveOrder), , drop = FALSE]

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    levels_kept <- unique(long$Level)
    pal <- grDevices::hcl.colors(length(levels_kept), "Dark 3")
    xlim <- c(0.8, length(fits) + 0.2)
    y_range <- range(c(long$Estimate, long$CI_Lower, long$CI_Upper),
                     finite = TRUE, na.rm = TRUE)
    plot(
      x = numeric(0), y = numeric(0),
      xlim = xlim, ylim = y_range,
      xaxt = "n", xlab = "Wave", ylab = "Severity (logit)",
      main = ""
    )
    axis(1, at = seq_along(fits), labels = names(fits),
                   las = 1)
    abline(h = 0, lty = 2, col = style$neutral)
    for (k in seq_along(levels_kept)) {
      lvl <- levels_kept[k]
      sub <- long[long$Level == lvl, , drop = FALSE]
      lines(sub$WaveOrder, sub$Estimate,
                      col = pal[k], lwd = 2)
      points(sub$WaveOrder, sub$Estimate,
                       pch = 19, col = pal[k])
      valid <- is.finite(sub$CI_Lower) & is.finite(sub$CI_Upper)
      if (any(valid)) {
        segments(sub$WaveOrder[valid], sub$CI_Lower[valid],
                           sub$WaveOrder[valid], sub$CI_Upper[valid],
                           col = pal[k], lwd = 1)
      }
    }
    legend("topright", legend = levels_kept, col = pal,
                     lwd = 2, bty = "n", cex = 0.75, inset = 0.02)
  }

  out <- new_mfrm_plot_data(
    "rater_trajectory",
    list(
      data = long[, setdiff(names(long), "WaveOrder"), drop = FALSE],
      facet = facet,
      ci_level = ci_level,
      title = sprintf("%s severity trajectory", facet),
      subtitle = sprintf("%d wave(s), %d level(s); %g%% CI whiskers",
                         length(fits), length(unique(long$Level)),
                         round(100 * ci_level)),
      preset = style$name,
      reference_lines = new_reference_lines(
        "h", 0, "Sum-to-zero reference", "dashed", "reference"
      )
    )
  )
  invisible(out)
}


#' Pairwise rater-agreement heatmap
#'
#' Summarizes inter-rater agreement as a symmetric rater x rater
#' heatmap. Cells are coloured by the chosen agreement metric: exact
#' agreement proportion by default, or the Pearson-style `Corr` column
#' from [interrater_agreement_table()] when `metric = "correlation"`.
#' The plot is a compact alternative to [plot(x, type = "agreement")][plot.mfrm_fit]'s
#' bar chart when the rater count exceeds ~6 pairs.
#'
#' @param fit An `mfrm_fit`.
#' @param diagnostics Optional [diagnose_mfrm()] output; piped through
#'   to [interrater_agreement_table()] when supplied.
#' @param rater_facet Name of the rater facet (default `"Rater"`).
#' @param metric Column to colour by: `"exact"` (default) or
#'   `"correlation"`. Quadratic-weighted kappa is not currently
#'   computed by [interrater_agreement_table()] and is therefore not
#'   offered here.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` object whose `data` slot bundles the
#'   rater x rater matrix, raw pairwise rows, matching settings and palette.
#'   [as_ggplot()] preserves the selected metric, facet labels and unavailable
#'   cells, with a fixed 0--1 exact-agreement or -1--1 correlation scale.
#'   The `"monochrome"` preset applies to both renderers. Self comparisons
#'   are unavailable, not estimated perfect agreement. Save the plot payload
#'   with [base::saveRDS()] to convert it without recomputing diagnostics;
#'   interpretation notes remain in `plot_data(p, "notes")`.
#'
#' @seealso [interrater_agreement_table()] for the underlying numeric
#'   table; [plot(x, type = "guttman")][plot.mfrm_fit] for a complementary
#'   person-by-element view of residual structure;
#'   [diagnose_mfrm()] for the diagnostics bundle the heatmap
#'   reads from.
#'
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("plot_rater_agreement_heatmap", run.dontrun = TRUE)
#' \dontrun{
#' # Saved synthetic JML/RSM analysis of example_core.
#' # For your own data, fit a model first; see ?fit_mfrm.
#' example_fit <- readRDS(system.file("examples", "core-jml.rds", package = "mfrmr"))
#' fit <- example_fit$fit
#' p <- plot(fit, type = "agreement", view = "heatmap", draw = FALSE)
#' dim(p$data$matrix)
#' # Look for (default `metric = "exact"`):
#' # - Off-diagonal cells close to the corresponding entry of
#' #   `summary(diag)$interrater$ExactAgreement` indicate consistent
#' #   pair behaviour; cells well below the average mark a pair
#' #   that disagrees more than the rest.
#' # - With `metric = "correlation"` the colour scale switches to
#' #   `[-1, 1]`; positive cells = pairs agree on relative ordering,
#' #   negative cells = pairs systematically rank persons in opposite
#' #   directions. Review context, exposure and sample size before deciding
#' #   whether any rater needs follow-up.
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_rater_agreement_heatmap <- function(fit,
                                         diagnostics = NULL,
                                         rater_facet = "Rater",
                                         metric = c("exact", "correlation"),
                                         preset = c("standard", "publication", "compact", "monochrome"),
                                         draw = TRUE) {
  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  metric <- match.arg(metric)
  style <- resolve_plot_preset(preset)
  agree <- interrater_agreement_table(fit, diagnostics = diagnostics,
                                       rater_facet = rater_facet)
  if (is.null(agree) || is.null(agree$pairs)) {
    stop("interrater_agreement_table() did not return pairwise rows.",
         call. = FALSE)
  }
  pairs <- as.data.frame(agree$pairs, stringsAsFactors = FALSE)
  if (!all(c("Rater1", "Rater2") %in% names(pairs))) {
    stop("Pairwise agreement table requires Rater1 and Rater2 columns.",
         call. = FALSE)
  }
  candidate_cols <- switch(
    metric,
    exact = c("Exact", "ExactAgreement", "Exact_Agreement"),
    correlation = c("Corr", "Correlation", "Pearson")
  )
  hit <- which(candidate_cols %in% names(pairs))
  value_col <- if (length(hit) > 0L) candidate_cols[hit[1]] else NA_character_
  if (is.na(value_col)) {
    stop("No column matching metric '", metric,
         "' was found in the pairwise table (available: ",
         paste(names(pairs), collapse = ", "), ").", call. = FALSE)
  }
  pairs$Value <- suppressWarnings(as.numeric(pairs[[value_col]]))

  raters <- sort(unique(c(pairs$Rater1, pairs$Rater2)))
  mat <- matrix(NA_real_, nrow = length(raters), ncol = length(raters),
                dimnames = list(raters, raters))
  # Self comparisons are not estimated agreement evidence.
  diag(mat) <- NA_real_
  for (k in seq_len(nrow(pairs))) {
    i <- match(pairs$Rater1[k], raters)
    j <- match(pairs$Rater2[k], raters)
    if (is.finite(i) && is.finite(j)) {
      mat[i, j] <- pairs$Value[k]
      mat[j, i] <- pairs$Value[k]
    }
  }

  zlim <- if (identical(metric, "correlation")) c(-1, 1) else c(0, 1)
  palette_name <- if (identical(metric, "correlation")) "RdBu" else "YlGnBu"
  cols <- if (identical(style$name, "monochrome")) {
    grDevices::gray.colors(20L, start = 0.95, end = 0.15)
  } else grDevices::hcl.colors(20L, palette_name, rev = TRUE)
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_par <- graphics::par()["mar"]
    on.exit(graphics::par(old_par), add = TRUE)
    label_lines <- max(graphics::strwidth(raters, units = "inches", cex = .8)) /
      (graphics::par("csi") * graphics::par("mex"))
    graphics::par(mar = c(rep(max(4, label_lines + 2.5), 2), 1, 1))
    image(
      x = seq_len(ncol(mat)), y = seq_len(nrow(mat)),
      z = mat, col = cols, zlim = zlim,
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "",
      main = ""
    )
    axis(1, at = seq_along(raters), labels = raters,
                   las = 2, cex.axis = 0.8)
    axis(2, at = seq_along(raters), labels = raters,
                   las = 1, cex.axis = 0.8)
    mtext(agree$settings$rater_facet, side = 1, line = label_lines + 1.4)
    mtext(agree$settings$rater_facet, side = 2, line = label_lines + 1.4)
    for (i in seq_along(raters)) for (j in seq_along(raters)) {
      v <- mat[i, j]
      if (is.finite(v)) {
        color_index <- findInterval(v, seq(zlim[1], zlim[2], length.out = length(cols) + 1L),
          all.inside = TRUE)
        luminance <- sum(grDevices::col2rgb(cols[color_index]) / 255 * c(.299, .587, .114))
        text(j, i, sprintf("%.2f", v), cex = 0.7,
                       col = if (luminance < .5) "white" else "black")
      }
    }
  }

  out <- new_mfrm_plot_data(
    "rater_agreement_heatmap",
    list(
      matrix = mat,
      pairs = pairs,
      metric = metric,
      settings = agree$settings,
      colors = cols,
      title = sprintf("Pairwise rater agreement (%s)", metric),
      subtitle = sprintf("%d raters; %d / %d pair comparisons available. Blank cells are unassessed; self-comparisons are omitted.",
                         length(raters), sum(is.finite(pairs$Value)), nrow(pairs)),
      notes = agree$notes,
      preset = style$name
    )
  )
  invisible(out)
}

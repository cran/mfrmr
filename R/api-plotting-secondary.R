# ==============================================================================
# Secondary visualization helpers
# ==============================================================================
#
# Four helpers grouped here for follow-up review after the primary
# Wright / pathway / CCC / fit / dashboard layer:
#
# - plot_local_dependence_heatmap: pairwise standardized residuals
#   between facet levels (Q3-style heatmap).
# - plot_reliability_snapshot: facet x separation / strata /
#   reliability bar overview.
# - plot_residual_matrix: person x facet-level standardized residuals
#   for follow-up after the Guttman scalogram.
# - plot_shrinkage_funnel: empirical-Bayes shrinkage caterpillar /
#   funnel display.
#
# Each helper follows the mfrmr plot contract: accept an `mfrm_fit`
# (or related class), resolve a preset via resolve_plot_preset() +
# apply_plot_preset(), and return an mfrm_plot_data object that
# downstream code can re-render or export.
# ==============================================================================


#' Pairwise standardized-residual heatmap for local-dependence review
#'
#' Builds an N x N heatmap of pairwise standardized residuals between
#' facet levels, computed from the diagnostics observation table.
#' Cells with large absolute values flag pairs of facet elements (e.g.
#' two raters, two items) whose residuals co-move more than the
#' main-effects MFRM expects, which is the standard Yen Q3-style
#' indicator of local response dependence.
#'
#' This helper complements [plot_marginal_pairwise()]: the marginal
#' version uses posterior-integrated agreement residuals on a
#' top-N pair list, while this view shows every pair on a shared color
#' scale so an analyst can scan for diagonal blocks or hotspots.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output. Computed
#'   on demand when omitted.
#' @param facet Facet whose levels are placed on both axes (default
#'   `"Rater"`).
#' @param min_pairs Minimum number of shared response opportunities
#'   required to retain a pair. Pairs below the threshold are shown
#'   as `NA`.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` whose `data` slot bundles the symmetric
#'   residual `matrix`, the long-form `pairs` table, and the threshold
#'   used.
#' @seealso [plot_marginal_pairwise()], [plot_qc_dashboard()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
#' p <- plot_local_dependence_heatmap(fit, draw = FALSE)
#' dim(p$data$matrix)
#' # Look for: |off-diagonal correlation| < 0.2 is the typical
#' #   acceptable regime; values >= 0.3 (Yen 1984 / Marais 2013
#' #   guideline) flag pairs that may share dependence beyond the
#' #   main-effects MFRM. Inspect those cells in `diag$obs`.
#' @export
plot_local_dependence_heatmap <- function(fit,
                                          diagnostics = NULL,
                                          facet = "Rater",
                                          min_pairs = 5L,
                                          preset = c("standard", "publication", "compact", "monochrome"),
                                          draw = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  facet <- as.character(facet[1])
  facet_names <- as.character(fit$config$facet_names %||% character(0))
  if (!facet %in% facet_names) {
    stop("`facet` must be one of: ", paste(facet_names, collapse = ", "), ".",
         call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(diagnose_mfrm(fit, residual_pca = "none",
                                                    diagnostic_mode = "legacy"))
  }
  obs <- as.data.frame(diagnostics$obs %||% data.frame(),
                        stringsAsFactors = FALSE)
  needed <- c(facet, "Person", "StdResidual")
  missing <- setdiff(needed, names(obs))
  if (length(missing) > 0L) {
    stop("`diagnostics$obs` is missing required columns: ",
         paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  obs[[facet]] <- as.character(obs[[facet]])
  obs$Person <- as.character(obs$Person)
  obs$StdResidual <- suppressWarnings(as.numeric(obs$StdResidual))
  obs <- obs[is.finite(obs$StdResidual), , drop = FALSE]

  wide <- tryCatch(
    tidyr::pivot_wider(
      obs[, c("Person", facet, "StdResidual"), drop = FALSE],
      id_cols = "Person",
      names_from = !!rlang::sym(facet),
      values_from = "StdResidual",
      values_fn = mean
    ),
    error = function(e) NULL
  )
  if (is.null(wide)) {
    stop("Could not pivot residuals to a person x level matrix.",
         call. = FALSE)
  }
  rater_cols <- setdiff(names(wide), "Person")
  if (length(rater_cols) < 2L) {
    stop("At least two distinct levels are required for the heatmap.",
         call. = FALSE)
  }

  rater_cols <- sort(rater_cols)
  mat <- matrix(NA_real_, nrow = length(rater_cols), ncol = length(rater_cols),
                 dimnames = list(rater_cols, rater_cols))
  pairs_rows <- list()
  for (i in seq_along(rater_cols)) {
    for (j in seq_along(rater_cols)) {
      if (i == j) {
        mat[i, j] <- 1
        next
      }
      a <- suppressWarnings(as.numeric(wide[[rater_cols[i]]]))
      b <- suppressWarnings(as.numeric(wide[[rater_cols[j]]]))
      ok <- is.finite(a) & is.finite(b)
      n_pair <- sum(ok)
      if (n_pair < as.integer(min_pairs)) {
        next
      }
      r <- suppressWarnings(stats::cor(a[ok], b[ok], use = "complete.obs"))
      if (!is.finite(r)) next
      mat[i, j] <- r
      pairs_rows[[length(pairs_rows) + 1L]] <- data.frame(
        Level1 = rater_cols[i],
        Level2 = rater_cols[j],
        ResidualCor = r,
        N = n_pair,
        stringsAsFactors = FALSE
      )
    }
  }
  pairs_df <- if (length(pairs_rows) > 0L) {
    do.call(rbind, pairs_rows)
  } else {
    data.frame(Level1 = character(0), Level2 = character(0),
               ResidualCor = numeric(0), N = integer(0),
               stringsAsFactors = FALSE)
  }

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    cols <- grDevices::hcl.colors(20L, "RdBu", rev = TRUE)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(5, 5, 3, 2))
    graphics::image(
      x = seq_len(ncol(mat)), y = seq_len(nrow(mat)),
      z = mat, col = cols, zlim = c(-1, 1),
      xaxt = "n", yaxt = "n", xlab = facet, ylab = facet,
      main = sprintf("Pairwise residual correlation (Q3 view): %s", facet)
    )
    graphics::axis(1, at = seq_along(rater_cols), labels = rater_cols,
                    las = 2, cex.axis = 0.8)
    graphics::axis(2, at = seq_along(rater_cols), labels = rater_cols,
                    las = 1, cex.axis = 0.8)
    for (i in seq_along(rater_cols)) for (j in seq_along(rater_cols)) {
      v <- mat[i, j]
      if (is.finite(v)) {
        graphics::text(j, i, sprintf("%.2f", v), cex = 0.7,
                        col = if (abs(v) > 0.6) "white" else "black")
      }
    }
  }

  invisible(new_mfrm_plot_data(
    "local_dependence_heatmap",
    list(
      matrix = mat,
      pairs = pairs_df,
      facet = facet,
      min_pairs = as.integer(min_pairs),
      title = sprintf("Pairwise residual correlation (Q3 view): %s", facet),
      subtitle = sprintf("%d level(s); pairs with N >= %d retained",
                          length(rater_cols), as.integer(min_pairs)),
      preset = style$name
    )
  ))
}


#' Facet reliability and separation snapshot bar plot
#'
#' Compact facet-level visual of the Wright & Masters (1982)
#' separation, strata, and reliability indices that
#' [diagnose_mfrm()] computes. Helpful as a single small figure for
#' "are persons / raters / criteria distinguishable?" review.
#' These are Rasch/FACETS-style separation indices on the fitted logit
#' scale, not ICCs; use [compute_facet_icc()] for the complementary
#' observed-score variance-share view.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output. Computed on
#'   demand when omitted.
#' @param metric `"reliability"` (default), `"separation"`, or
#'   `"strata"`.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` whose `data` slot bundles a tidy
#'   `Facet`, `Metric`, `Value` data frame.
#' @seealso [diagnose_mfrm()], [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
#' p <- plot_reliability_snapshot(fit, draw = FALSE)
#' p$data$table
#' # Look for (default `metric = "reliability"`):
#' # - >= 0.9 strong, 0.7-0.9 adequate, < 0.7 weak (Wright & Masters 1982).
#' # - The Person row is the operative reliability for ability scores.
#' # - Non-Person rows (Rater / Criterion) report the same index but
#' #   should be read as "are facet elements distinguishable?"; values
#' #   close to 1 mean facet means differ reliably from each other.
#' @export
plot_reliability_snapshot <- function(fit,
                                      diagnostics = NULL,
                                      metric = c("reliability", "separation", "strata"),
                                      preset = c("standard", "publication", "compact", "monochrome"),
                                      draw = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  metric <- match.arg(metric)
  style <- resolve_plot_preset(preset)
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(diagnose_mfrm(fit, residual_pca = "none",
                                                   diagnostic_mode = "legacy"))
  }
  rel <- as.data.frame(diagnostics$reliability %||% data.frame(),
                        stringsAsFactors = FALSE)
  if (nrow(rel) == 0L || !"Facet" %in% names(rel)) {
    stop("`diagnostics$reliability` is empty or missing the Facet column.",
         call. = FALSE)
  }
  metric_col <- switch(metric,
    reliability = if ("Reliability" %in% names(rel)) "Reliability" else NA_character_,
    separation = if ("Separation" %in% names(rel)) "Separation" else NA_character_,
    strata = if ("Strata" %in% names(rel)) "Strata" else NA_character_
  )
  if (is.na(metric_col)) {
    stop("Diagnostics$reliability does not include the requested metric column: ",
         metric, ".", call. = FALSE)
  }
  vals <- suppressWarnings(as.numeric(rel[[metric_col]]))
  ok <- is.finite(vals)
  rel <- rel[ok, , drop = FALSE]
  vals <- vals[ok]
  if (length(vals) == 0L) {
    stop("No finite values for the requested metric.", call. = FALSE)
  }
  ord <- order(vals, decreasing = TRUE)
  rel <- rel[ord, , drop = FALSE]
  vals <- vals[ord]
  payload <- data.frame(
    Facet = as.character(rel$Facet),
    Metric = metric,
    Value = vals,
    stringsAsFactors = FALSE
  )

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(5, 6, 3, 2))
    bar_max <- if (identical(metric, "reliability")) {
      max(1.0, max(payload$Value, na.rm = TRUE) * 1.05)
    } else {
      max(payload$Value, na.rm = TRUE) * 1.10
    }
    bp <- graphics::barplot(
      payload$Value, names.arg = payload$Facet,
      horiz = TRUE, las = 1,
      xlim = c(0, bar_max),
      col = grDevices::adjustcolor(style$accent_primary, alpha.f = 0.7),
      border = style$accent_primary,
      xlab = paste0(metric, " (", metric_col, ")"),
      main = sprintf("Facet %s snapshot", metric)
    )
    graphics::text(payload$Value, bp,
                    labels = sprintf("%.2f", payload$Value),
                    pos = 4, cex = 0.85, col = style$neutral)
    if (identical(metric, "reliability")) {
      graphics::abline(v = c(0.7, 0.9), lty = c(3, 2), col = style$grid)
      graphics::mtext("0.7 = adequate, 0.9 = strong (Wright & Masters, 1982)",
                       side = 1, line = 3.5, cex = 0.7)
    }
  }

  invisible(new_mfrm_plot_data(
    "reliability_snapshot",
    list(
      table = payload,
      metric = metric,
      title = sprintf("Facet %s snapshot", metric),
      subtitle = paste0("Sorted by ", metric_col, "; N facet rows = ",
                         nrow(payload)),
      preset = style$name,
      reference_lines = if (identical(metric, "reliability")) {
        new_reference_lines(
          axis = c("v", "v"),
          value = c(0.7, 0.9),
          label = c("Adequate (0.7)", "Strong (0.9)"),
          linetype = c("dotted", "dashed"),
          role = c("threshold", "threshold")
        )
      } else NULL
    )
  ))
}


#' Person x facet-level standardized-residual matrix
#'
#' Visualizes the person x element matrix of standardized residuals
#' from [diagnose_mfrm()] as a heatmap. Complements
#' [plot_guttman_scalogram()] (which shows raw responses) by exposing
#' the residual structure directly: large positive cells show
#' under-prediction, negative cells over-prediction.
#'
#' @param fit An `mfrm_fit` from [fit_mfrm()].
#' @param diagnostics Optional [diagnose_mfrm()] output. Computed on
#'   demand when omitted.
#' @param facet Facet whose levels become the column axis (default
#'   `"Rater"`).
#' @param top_n_persons Cap on the number of rows. Defaults to 40
#'   to keep the figure legible; persons are kept by largest absolute
#'   residual mean.
#' @param preset Visual preset.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` whose `data` slot bundles the residual
#'   `matrix` (rows = Person, columns = facet level) and the long-form
#'   `obs` table.
#' @seealso [plot_guttman_scalogram()], [plot_unexpected()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
#' p <- plot_residual_matrix(fit, top_n_persons = 12, draw = FALSE)
#' dim(p$data$matrix)
#' # Look for: cell values within ~|2| are routine; |residual| > 2 is
#' #   misfit at the 5% level and |residual| > 3 at the 1% level
#' #   (Wright & Linacre 1994). Persons with multiple high-magnitude
#' #   cells across the same facet level point at scoring drift.
#' @export
plot_residual_matrix <- function(fit,
                                 diagnostics = NULL,
                                 facet = "Rater",
                                 top_n_persons = 40L,
                                 preset = c("standard", "publication", "compact", "monochrome"),
                                 draw = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  facet <- as.character(facet[1])
  facet_names <- as.character(fit$config$facet_names %||% character(0))
  if (!facet %in% facet_names) {
    stop("`facet` must be one of: ", paste(facet_names, collapse = ", "), ".",
         call. = FALSE)
  }
  style <- resolve_plot_preset(preset)
  if (is.null(diagnostics)) {
    diagnostics <- suppressMessages(diagnose_mfrm(fit, residual_pca = "none",
                                                    diagnostic_mode = "legacy"))
  }
  obs <- as.data.frame(diagnostics$obs %||% data.frame(),
                        stringsAsFactors = FALSE)
  needed <- c("Person", facet, "StdResidual")
  missing <- setdiff(needed, names(obs))
  if (length(missing) > 0L) {
    stop("`diagnostics$obs` is missing required columns: ",
         paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  obs$Person <- as.character(obs$Person)
  obs[[facet]] <- as.character(obs[[facet]])
  obs$StdResidual <- suppressWarnings(as.numeric(obs$StdResidual))

  agg <- stats::aggregate(
    obs$StdResidual,
    by = list(Person = obs$Person, Level = obs[[facet]]),
    FUN = function(z) mean(z, na.rm = TRUE)
  )
  names(agg)[3] <- "StdResidual"

  person_priority <- stats::aggregate(
    abs(agg$StdResidual),
    by = list(Person = agg$Person),
    FUN = function(z) mean(z, na.rm = TRUE)
  )
  names(person_priority)[2] <- "AbsMean"
  person_priority <- person_priority[order(person_priority$AbsMean,
                                            decreasing = TRUE), , drop = FALSE]
  cap <- max(1L, as.integer(top_n_persons))
  if (nrow(person_priority) > cap) {
    person_priority <- person_priority[seq_len(cap), , drop = FALSE]
  }
  agg <- agg[agg$Person %in% person_priority$Person, , drop = FALSE]

  persons <- as.character(person_priority$Person)
  levels_x <- sort(unique(as.character(agg$Level)))
  mat <- matrix(NA_real_, nrow = length(persons), ncol = length(levels_x),
                 dimnames = list(persons, levels_x))
  ix <- match(agg$Person, persons)
  jx <- match(agg$Level, levels_x)
  ok <- is.finite(ix) & is.finite(jx)
  mat[cbind(ix[ok], jx[ok])] <- agg$StdResidual[ok]

  if (isTRUE(draw)) {
    apply_plot_preset(style)
    cols <- grDevices::hcl.colors(40L, "RdBu", rev = TRUE)
    zmax <- max(3, max(abs(mat), na.rm = TRUE))
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(5, 6, 3, 2))
    graphics::image(
      x = seq_along(levels_x), y = seq_along(persons),
      z = t(mat), col = cols, zlim = c(-zmax, zmax),
      xaxt = "n", yaxt = "n", xlab = facet, ylab = "Person",
      main = sprintf("Standardized residual matrix (Person x %s)", facet)
    )
    graphics::axis(1, at = seq_along(levels_x), labels = levels_x,
                    las = 2, cex.axis = 0.7)
    graphics::axis(2, at = seq_along(persons),
                    labels = persons, las = 1, cex.axis = 0.7)
  }

  invisible(new_mfrm_plot_data(
    "residual_matrix",
    list(
      matrix = mat,
      obs = agg,
      facet = facet,
      title = sprintf("Standardized residual matrix (Person x %s)", facet),
      subtitle = sprintf("Top %d person(s) by mean |StdResidual|",
                          length(persons)),
      preset = style$name,
      reference_lines = new_reference_lines(
        "h", 0, "Zero residual", "dashed", "reference"
      )
    )
  ))
}


#' Empirical-Bayes shrinkage funnel / caterpillar
#'
#' Visualizes empirical-Bayes shrinkage by drawing one row per facet
#' level with the raw (pre-shrinkage) and shrunken estimates plus the
#' shrinkage factor. Rows are ordered by absolute shrinkage so the
#' levels that move most under the prior appear at the top.
#'
#' Requires a fit produced via [apply_empirical_bayes_shrinkage()] or
#' a `fit_mfrm(..., facet_shrinkage = "empirical_bayes")` run, so that
#' `fit$facets$others` carries `Estimate`, `ShrunkEstimate`, and
#' `ShrinkageFactor` columns.
#'
#' @param fit An `mfrm_fit` augmented with empirical-Bayes shrinkage.
#' @param facet Facet to draw (default: first non-person facet with
#'   shrinkage columns present).
#' @param top_n Maximum number of rows to draw (default 30).
#' @param preset Visual preset.
#' @param show_ci Logical. When `TRUE`, draw approximate confidence-interval
#'   whiskers for raw and shrunken estimates when `SE` / `ShrunkSE` evidence is
#'   available.
#' @param ci_level Confidence level used when `show_ci = TRUE`; default 0.95.
#' @param draw If `TRUE`, draw with base graphics.
#'
#' @return An `mfrm_plot_data` whose `data` slot bundles the long
#'   `Level`, `RawEstimate`, `ShrunkEstimate`, `ShrinkageFactor` table.
#'   When `show_ci = TRUE`, the table also includes `RawCI_Lower`,
#'   `RawCI_Upper`, `ShrunkCI_Lower`, `ShrunkCI_Upper`, and `CI_Level`.
#' @seealso [apply_empirical_bayes_shrinkage()],
#'   [mfrmr_visual_diagnostics]
#' @examples
#' toy <- load_mfrmr_data("example_core")
#' fit <- fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
#'                  method = "JML", maxit = 30)
#' fit_eb <- apply_empirical_bayes_shrinkage(fit)
#' p <- plot_shrinkage_funnel(fit_eb, draw = FALSE)
#' head(p$data$table)
#' # Look for: short segments (Raw and Shrunken close together) =
#' #   little pooling. Long segments fanning toward the centre = the
#' #   prior pulled the estimate strongly; this is most pronounced for
#' #   small-N levels. ShrinkageFactor near 1 means most of the
#' #   movement was driven by the prior rather than the data.
#' @export
plot_shrinkage_funnel <- function(fit,
                                  facet = NULL,
                                  top_n = 30L,
                                  preset = c("standard", "publication", "compact", "monochrome"),
                                  show_ci = FALSE,
                                  ci_level = 0.95,
                                  draw = TRUE) {
  if (!inherits(fit, "mfrm_fit")) {
    stop("`fit` must be an mfrm_fit object from fit_mfrm().", call. = FALSE)
  }
  if (isTRUE(show_ci)) {
    ci_level <- .validate_shrinkage_ci_level(ci_level)
  }
  others <- as.data.frame(fit$facets$others %||% data.frame(),
                           stringsAsFactors = FALSE)
  needed <- c("Facet", "Level", "Estimate", "ShrunkEstimate", "ShrinkageFactor")
  missing <- setdiff(needed, names(others))
  if (length(missing) > 0L) {
    stop("This fit does not carry empirical-Bayes shrinkage columns. ",
         "Run apply_empirical_bayes_shrinkage(fit) first. Missing: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(facet)) {
    facet <- as.character(others$Facet[1])
  }
  facet <- as.character(facet[1])
  others <- others[as.character(others$Facet) == facet, , drop = FALSE]
  if (nrow(others) == 0L) {
    stop("No rows in fit$facets$others for facet '", facet, "'.",
         call. = FALSE)
  }
  others$Estimate <- suppressWarnings(as.numeric(others$Estimate))
  others$ShrunkEstimate <- suppressWarnings(as.numeric(others$ShrunkEstimate))
  others$ShrinkageFactor <- suppressWarnings(as.numeric(others$ShrinkageFactor))
  se_col <- if ("ModelSE" %in% names(others)) {
    "ModelSE"
  } else if ("SE" %in% names(others)) {
    "SE"
  } else {
    NA_character_
  }
  others$SE <- if (!is.na(se_col)) {
    suppressWarnings(as.numeric(others[[se_col]]))
  } else {
    rep(NA_real_, nrow(others))
  }
  others$ShrunkSE <- if ("ShrunkSE" %in% names(others)) {
    suppressWarnings(as.numeric(others$ShrunkSE))
  } else {
    rep(NA_real_, nrow(others))
  }
  others <- others[is.finite(others$Estimate) &
                     is.finite(others$ShrunkEstimate), , drop = FALSE]
  if (nrow(others) == 0L) {
    stop("No rows had finite Estimate / ShrunkEstimate for facet '", facet, "'.",
         call. = FALSE)
  }
  ord <- order(abs(others$Estimate - others$ShrunkEstimate),
                decreasing = TRUE)
  others <- others[ord, , drop = FALSE]
  cap <- max(1L, as.integer(top_n))
  if (nrow(others) > cap) {
    others <- others[seq_len(cap), , drop = FALSE]
  }
  payload <- data.frame(
    Facet = as.character(others$Facet),
    Level = as.character(others$Level),
    RawEstimate = others$Estimate,
    SE = others$SE,
    ShrunkEstimate = others$ShrunkEstimate,
    ShrunkSE = others$ShrunkSE,
    ShrinkageFactor = others$ShrinkageFactor,
    Movement = others$ShrunkEstimate - others$Estimate,
    stringsAsFactors = FALSE
  )
  payload <- payload[order(payload$RawEstimate), , drop = FALSE]
  payload$RowOrder <- seq_len(nrow(payload))
  if (isTRUE(show_ci)) {
    payload <- .add_shrinkage_ci_columns(
      payload,
      ci_level = ci_level,
      raw_estimate_col = "RawEstimate",
      raw_se_col = "SE",
      raw_prefix = "Raw",
      shrunk_estimate_col = "ShrunkEstimate",
      shrunk_se_col = "ShrunkSE",
      shrunk_prefix = "Shrunk"
    )
  }

  style <- resolve_plot_preset(preset)
  if (isTRUE(draw)) {
    apply_plot_preset(style)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(5, 6, 3, 2))
    x_values <- c(payload$RawEstimate, payload$ShrunkEstimate)
    if (isTRUE(show_ci)) {
      x_values <- c(x_values, payload$RawCI_Lower, payload$RawCI_Upper,
                    payload$ShrunkCI_Lower, payload$ShrunkCI_Upper)
    }
    xr <- range(x_values[is.finite(x_values)], na.rm = TRUE)
    if (!all(is.finite(xr))) xr <- c(-1, 1)
    if (diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
    xr <- xr + c(-0.06, 0.06) * diff(xr)
    graphics::plot(NA, xlim = xr, ylim = c(0.5, nrow(payload) + 0.5),
                    yaxt = "n",
                    xlab = "Estimate (logits)", ylab = facet,
                    main = sprintf("Empirical-Bayes shrinkage funnel: %s",
                                    facet))
    graphics::abline(v = 0, lty = 3, col = style$grid)
    graphics::axis(2, at = payload$RowOrder, labels = payload$Level,
                    las = 1, cex.axis = 0.7)
    graphics::segments(payload$RawEstimate, payload$RowOrder,
                        payload$ShrunkEstimate, payload$RowOrder,
                        col = style$neutral, lwd = 1.2)
    if (isTRUE(show_ci)) {
      graphics::segments(payload$RawCI_Lower, payload$RowOrder - 0.12,
                         payload$RawCI_Upper, payload$RowOrder - 0.12,
                         col = style$accent_primary, lwd = 1)
      graphics::segments(payload$ShrunkCI_Lower, payload$RowOrder + 0.12,
                         payload$ShrunkCI_Upper, payload$RowOrder + 0.12,
                         col = style$accent_tertiary, lwd = 1)
    }
    graphics::points(payload$RawEstimate, payload$RowOrder,
                      pch = 1, col = style$accent_primary, cex = 1.0)
    graphics::points(payload$ShrunkEstimate, payload$RowOrder,
                      pch = 19, col = style$accent_tertiary, cex = 1.0)
    legend_labels <- c("Raw", "Shrunken")
    legend_pch <- c(1, 19)
    legend_lty <- c(NA, NA)
    legend_col <- c(style$accent_primary, style$accent_tertiary)
    if (isTRUE(show_ci)) {
      legend_labels <- c(legend_labels,
                         sprintf("%g%% CI whisker", round(100 * ci_level)))
      legend_pch <- c(legend_pch, NA)
      legend_lty <- c(legend_lty, 1)
      legend_col <- c(legend_col, style$accent_primary)
    }
    graphics::legend("topright",
                      legend = legend_labels,
                      pch = legend_pch, lty = legend_lty, col = legend_col,
                      bty = "n", cex = 0.85, inset = 0.02)
  }

  legend_label <- c("Raw estimate", "Shrunken estimate",
                    "Shrinkage movement", "Sum-to-zero reference")
  legend_role <- c("location", "location", "movement", "reference")
  legend_aesthetic <- c("point", "point", "line", "line")
  legend_value <- c(style$accent_primary, style$accent_tertiary,
                    style$neutral, style$grid)
  if (isTRUE(show_ci)) {
    legend_label <- c(legend_label, sprintf("%g%% CI whisker",
                                            round(100 * ci_level)))
    legend_role <- c(legend_role, "interval")
    legend_aesthetic <- c(legend_aesthetic, "line")
    legend_value <- c(legend_value, style$accent_primary)
  }
  invisible(new_mfrm_plot_data(
    "shrinkage_funnel",
    list(
      table = payload,
      facet = facet,
      show_ci = isTRUE(show_ci),
      ci_level = if (isTRUE(show_ci)) ci_level else NA_real_,
      title = sprintf("Empirical-Bayes shrinkage funnel: %s", facet),
      subtitle = sprintf("%d level(s); rows ordered by raw estimate",
                          nrow(payload)),
      preset = style$name,
      legend = new_plot_legend(
        label = legend_label,
        role = legend_role,
        aesthetic = legend_aesthetic,
        value = legend_value
      ),
      reference_lines = new_reference_lines(
        "v", 0, "Sum-to-zero reference", "dashed", "reference"
      )
    )
  ))
}

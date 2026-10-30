#' Uncertainty in GPCM category probabilities and information curves
#'
#' Evaluate specified rating contexts at known ability values and propagate the
#' full calibration covariance to category probabilities or per-rating Fisher
#' information. This is uncertainty in a fitted curve, not a Person score interval.
#'
#' @param fit A native GPCM MML fit. Two-family fits supply provisional point
#'   curves only; intervals and covariance adjustments are unavailable.
#' @param newdata Data frame with each fitted non-Person facet and a numeric
#'   `Theta` column. Each row is one rating context at one known ability value.
#'   Use original facet labels. Unknown levels and missing inputs are refused.
#'   Previously unobserved combinations of known levels are allowed and labelled
#'   in the returned `contexts` table.
#' @param type `"probability"` (one row per category) or `"information"`
#'   (one row per supplied context/ability). Information is per rating, not a
#'   sum over all observed exposures as in [compute_information()].
#' @param method,clusters,adjust,level,simultaneous As in [confint.mfrm_fit()].
#' @param object,x A saved `mfrm_curve_intervals` result.
#' @param title,subtitle Optional returned text; NULL removes it. Figure text
#'   is not drawn. Use `plot_data(p, "title")` and `plot_data(p, "notes")`.
#' @param palette Optional vector of colors, one per category/series. Line
#'   types distinguish series in addition to color. The default is colorblind friendly.
#' @param draw Draw the ggplot immediately? Default TRUE; FALSE returns it only.
#' @param caption Optional returned caption. If omitted, unavailable intervals are
#'   counted and explained, along with unobserved combinations of fitted levels.
#'   Two-family point-only curves instead state that calibration intervals
#'   are unavailable for the entire plot, without marking every curve point.
#'   NULL removes the caption without removing markers
#'   or the reasons in the saved table.
#' @param ... Unused.
#'
#' @details The complete GPCM predictor, including fitted facet interactions,
#'   is evaluated from the free parameter vector. A central-difference Jacobian
#'   propagates the full joint covariance. Probability limits use a logit delta
#'   approximation; information limits use a log delta approximation for
#'   `slope^2 * Var(category | Theta, context)`. Bounds respect their support.
#'   Rounded probabilities zero/one and nonpositive information keep the point
#'   but have unavailable intervals. The supplied Theta values are fixed on
#'   the fitted native scale; their estimation uncertainty is not included.
#'   Estimated population parameters enter the joint covariance as nuisance
#'   parameters, not as a request to transform or integrate the Theta grid.
#'
#'   Bonferroni applies to the finite collection of all output rows in this
#'   call, not to the continuous curve between them. Ribbons connect grid-point
#'   intervals for display. Sandwich interpretation and independence assumptions
#'   are those of [confint.mfrm_fit()]. These are asymptotic approximations;
#'   numerical agreement does not establish sampling coverage.
#'   A reanalysis of saved correctly specified GPCM fits found substantial
#'   undercoverage for probability intervals in small incomplete designs,
#'   including Bonferroni-adjusted finite-grid families. Adjustment cannot
#'   repair inaccurate marginal approximations. Refitting every dataset in the
#'   affected unequal-rater-slope condition with the updated optimizer left
#'   its interval results unchanged. See the GPCM scope vignette
#'   for the evaluated grid, denominators and source limitations; an available
#'   interval is not a finite-sample coverage certification.
#'   Printing and the returned plot description identify the approximation.
#'   Custom plot text may omit that description; retain the method and its
#'   limitations in the figure legend or accompanying report.
#'   For one-family curves, crosses mark retained estimates whose intervals
#'   are unavailable. Two-family curves have no calibration intervals; the
#'   returned subtitle and caption state this without covering curves in crosses.
#'   A context with only one supplied ability value retains colored points.
#'   Panel labels put each facet on its own line. With many rating contexts,
#'   supply the comparisons of interest in `newdata` or enlarge the exported
#'   figure; plotting does not select or discard contexts automatically. Ribbons
#'   stop at unavailable grid points; a missing ribbon does not mean zero
#'   uncertainty. Consult the table's `InferenceReview` for each reason.
#'
#' @section Two slope families:
#' For a two-family fit, the two component slopes are multiplied in each
#' specified context. The shared evaluator returns provisional fitted values,
#' including after numerical nonconvergence; it does not certify their
#' reliability. `SE`, `Lower` and `Upper` remain missing and `CIEligible` is
#' false. Printing and default plot text identify this scope without claiming
#' nominal intervals. Sandwich and simultaneous adjustments are not available.
#' The existing plot, saved-data, report and export routes preserve these values
#' and missing intervals. A slope value of one is not a rater-quality threshold.
#'
#' @return A `mfrm_curve_intervals` list with `table`, `settings`, and the exact
#'   `newdata`, plus `contexts` identifying each `InputRow` as an observed or
#'   unobserved combination of fitted facet levels (`ObservedContext`). This
#'   records the retained rating design, not statistical identification or
#'   interval reliability. Its plot method returns a ggplot with the plotted data available
#'   in `plot$data`; descriptions and cautions are available through
#'   `plot_data(p, "notes")`, separately from the figure.
#'   `summary()` returns settings, the complete interval table, contexts,
#'   numerical cautions and the anchor note without recalculating information.
#' @seealso [confint.mfrm_fit()], [category_curves_report()], [compute_information()]
#' @details A numerically verified but ill-conditioned information inverse can
#'   supply intervals with a warning, as described in [confint.mfrm_fit()].
#'   The warning is saved in `cautions` and the table's `InferenceReview`, and
#'   appears in printing and the returned plot description. A custom subtitle,
#'   including NULL, changes the display without removing saved diagnostics.
#'   In a one-family fit with fixed facet anchors or group means, model and
#'   sandwich intervals condition on those supplied values. Their source
#'   calibration uncertainty is not propagated. Fixing a facet location does
#'   not fix its slope, thresholds or whole response curve. Availability still
#'   depends on the joint-information and numerical checks described above.
#' @examples
#' # Supplementary workflow: run explicitly to keep routine checks brief.
#' # example("mfrm_curve_intervals", run.dontrun = TRUE)
#' \dontrun{
#' # How do category probabilities vary with ability in one rating context?
#' ratings <- load_mfrmr_data("example_core")
#' # Reuse the synthetic GPCM fit; estimation uses 61 quadrature points.
#' # Full regeneration recipe: examples/measurement-fits.R in the package.
#' fit <- readRDS(system.file("examples", "measurement-fits.rds",
#'   package = "mfrmr", mustWork = TRUE))$gpcm
#' # Select actual labels from the example; replace them with your context.
#' grid <- expand.grid(Theta = seq(-2, 2, length.out = 21),
#'   Rater = unique(ratings$Rater)[1],
#'   Criterion = unique(ratings$Criterion)[1])
#' curves <- mfrm_curve_intervals(fit, grid)
#' # At ability zero, read the probabilities of the original score categories.
#' subset(curves$table, Theta == 0,
#'   select = c(Category, Estimate, Lower, Upper, CIEligible, InferenceReview))
#' plot(curves)
#' # Ribbons describe calibration uncertainty at fixed ability, not uncertainty
#' # about an individual's ability. They are pointwise, not simultaneous bands.
#' file <- tempfile(fileext = ".rds")
#' saveRDS(curves, file) # Keep the complete result, including context and notes
#' saved <- readRDS(file)
#' saved$settings
#' }
#' @export
mfrm_curve_intervals <- function(fit, newdata, type = c("probability", "information"),
    method = c("model", "sandwich"), clusters = NULL, adjust = FALSE,
    level = .95, simultaneous = c("none", "bonferroni")) {
  stop_if_jml_adjustment(fit, "mfrm_curve_intervals")
  type <- match.arg(type); method <- match.arg(method); simultaneous <- match.arg(simultaneous)
  mfrm_random_rater_interval_level(level)
  point_only <- mfrm_has_product_slopes(fit)
  if (point_only && simultaneous != "none") {
    stop("Two-family curve intervals are unavailable; simultaneous adjustment is not supported.", call. = FALSE)
  }
  inference <- mfrm_gpcm_inference(fit, method, clusters, adjust)
  if (!is.null(inference$check$caution)) warning(inference$check$caution, call. = FALSE)
  facets <- fit$config$facet_names
  response <- mfrm_gpcm_response_evaluator(fit$prep, fit$config, newdata)
  newdata <- response$newdata
  evaluate <- function(par) {
    value <- response$evaluate(par)
    if (type == "probability") return(as.vector(t(value$probabilities)))
    value$information
  }
  estimate <- evaluate(fit$opt$par)
  n <- length(estimate)
  se <- rep(NA_real_, n)
  reason <- rep(if (isTRUE(inference$check$eligible))
    paste("Approximate calibration uncertainty at fixed ability and rating context; finite-sample coverage is not guaranteed.",
      inference$check$caution %||% "")
    else inference$check$review, n)
  if (isTRUE(inference$check$eligible) && length(inference$check$anchor_note)) {
    reason <- paste(reason, inference$check$anchor_note)
  }
  if (isTRUE(inference$check$eligible)) {
    jac <- mfrmr_numeric_transformation_jacobian(evaluate, fit$opt$par, relative_step = 1e-5)
    if (jac$valid) {
      variance <- rowSums((jac$jacobian %*% inference$covariance) * jac$jacobian)
      se <- ifelse(variance >= 0, sqrt(pmax(0, variance)), NA_real_)
    } else reason[] <- "The curve derivative could not be evaluated."
  }
  transform <- if (type == "probability") stats::qlogis else log
  inverse <- if (type == "probability") stats::plogis else exp
  interior <- is.finite(estimate) & estimate > 0 & (type != "probability" | estimate < 1)
  link_se <- se / if (type == "probability") (estimate * (1-estimate)) else estimate
  critical <- stats::qnorm((1-level)/(2 * if (simultaneous == "bonferroni") n else 1), lower.tail = FALSE)
  low <- high <- rep(NA_real_, n)
  usable <- interior & is.finite(link_se) & link_se > 0 & isTRUE(inference$check$eligible)
  low[usable] <- inverse(transform(estimate[usable]) - critical*link_se[usable])
  high[usable] <- inverse(transform(estimate[usable]) + critical*link_se[usable])
  usable <- usable & is.finite(low) & is.finite(high)
  low[!usable] <- high[!usable] <- NA_real_
  reason[!interior] <- "The fitted curve is at a numerical support boundary; an interval is unavailable."
  reason[interior & !usable & isTRUE(inference$check$eligible)] <- "Nonpositive variance, unavailable derivative or nonrepresentable interval."
  k <- if (type == "probability") fit$config$n_cat else 1L
  tab <- newdata[rep(seq_len(nrow(newdata)), each = k), , drop = FALSE]
  tab$InputRow <- rep(seq_len(nrow(newdata)), each = k)
  tab$Category <- if (type == "probability") rep(fit$prep$score_map$OriginalScore[
    match(seq_len(k)-1L+fit$prep$rating_min, fit$prep$score_map$InternalScore)], nrow(newdata)) else NA_real_
  tab$Estimate <- estimate; tab$SE <- se; tab$Lower <- low; tab$Upper <- high
  tab$CIEligible <- usable; tab$InferenceReview <- reason
  rownames(tab) <- NULL
  structure(list(table = tab, newdata = newdata, contexts = response$contexts,
    cautions = inference$check$caution,
    anchor_note = inference$check$anchor_note,
    source = mfrm_gpcm_inference_source(fit),
    settings = list(type = type, method = method, clusters = inference$clusters,
      adjust = adjust, level = level, simultaneous = simultaneous, facets = facets,
      point_only = point_only,
      target = if (point_only) "Provisional fitted curves at fixed ability and rating context; intervals unavailable" else
        "Calibration uncertainty at fixed native-scale ability and rating context")),
    class = "mfrm_curve_intervals")
}

#' @rdname mfrm_curve_intervals
#' @export
summary.mfrm_curve_intervals <- function(object, ...) {
  rlang::check_dots_empty()
  object[c("settings", "table", "contexts", "cautions", "anchor_note")]
}

#' @rdname mfrm_curve_intervals
#' @export
print.mfrm_curve_intervals <- function(x, ...) {
  if (isTRUE(x$settings$point_only)) {
    cat("Two-family GPCM", x$settings$type, "curves: provisional fitted values; intervals unavailable\n")
  } else cat("GPCM", x$settings$type, "curve intervals:", x$settings$method, "covariance\n")
  print(x$table[setdiff(names(x$table), "InferenceReview")], row.names = FALSE)
  cat("Available intervals:", sum(x$table$CIEligible), "of", nrow(x$table),
    if (isTRUE(x$settings$point_only)) "; fitted values at fixed ability values only.\n" else
      "; fixed ability values, calibration uncertainty only.\n")
  if (!is.null(x$contexts) && any(!x$contexts$ObservedContext)) {
    cat("Unobserved combinations of fitted levels:", sum(!x$contexts$ObservedContext),
      "of", nrow(x$contexts), "input rows; see contexts.\n")
  }
  if (!isTRUE(x$settings$point_only)) cat("Approximate intervals; numerical availability does not establish nominal coverage.\n")
  if (length(x$anchor_note)) cat(x$anchor_note, "\n")
  for (reason in unique(x$table$InferenceReview[!x$table$CIEligible])) {
    cat("  ", reason, "\n", sep = "")
  }
  for (caution in x$cautions) cat("Caution:", caution, "\n")
  invisible(x)
}

#' @rdname mfrm_curve_intervals
#' @export
plot.mfrm_curve_intervals <- function(x, title = "GPCM curve uncertainty",
    subtitle = paste("Approximate calibration intervals; fixed ability", x$settings$method,
      paste0(100*x$settings$level, "%"),
      if (x$settings$simultaneous == "none") "pointwise" else "Bonferroni grid points", sep = " | "),
    palette = NULL, draw = TRUE, caption = NULL, ...) {
  rlang::check_dots_empty()
  point_only <- isTRUE(x$settings$point_only)
  if (point_only) {
    if (missing(title)) title <- "Two-family GPCM fitted curves"
    if (missing(subtitle)) subtitle <- "Provisional fitted values; calibration intervals are unavailable"
  }
  if (missing(subtitle) && length(x$cautions)) subtitle <- paste(subtitle,
    paste(x$cautions, collapse = " "), sep = "\n")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 to plot curve intervals.", call. = FALSE)
  tab <- x$table
  automatic_caption <- missing(caption)
  if (automatic_caption && point_only) {
    caption <- "Fitted curves only; calibration intervals are unavailable for two-family GPCM."
  } else if (automatic_caption && any(!tab$CIEligible)) caption <- paste(
      "Intervals unavailable at", sum(!tab$CIEligible), "of", nrow(tab),
      "points. Crosses mark retained estimates; see the result table for reasons.")
  if (automatic_caption && !is.null(x$contexts) && any(!x$contexts$ObservedContext)) {
    caption <- paste(c(caption, paste(sum(!x$contexts$ObservedContext),
      "input rows combine fitted levels not observed together; see contexts.")), collapse = " ")
  }
  if (!is.null(subtitle)) subtitle <- paste(strwrap(subtitle, width = 75), collapse = "\n")
  if (!is.null(caption)) caption <- paste(strwrap(caption, width = 85), collapse = "\n")
  tab$Context <- apply(tab[x$settings$facets], 1L, function(z)
    paste(paste(x$settings$facets, z, sep = " = "), collapse = ", "))
  contexts <- tab[!duplicated(tab$Context), , drop = FALSE]
  context_labels <- setNames(apply(contexts[x$settings$facets], 1L, function(z)
    paste(paste(x$settings$facets, z, sep = " = "), collapse = "\n")), contexts$Context)
  tab$Series <- if (x$settings$type == "probability") factor(tab$Category) else factor("Information")
  isolated <- if (point_only) stats::ave(tab$Theta, tab$Context,
    FUN = function(z) length(unique(z))) == 1L else rep(FALSE, nrow(tab))
  # Keep a ribbon from bridging over an unavailable point in a series.
  index <- order(tab$Context, tab$Series, tab$Theta, tab$InputRow)
  tab$IntervalGroup <- integer(nrow(tab))
  tab$IntervalGroup[index] <- cumsum(!tab$CIEligible[index])
  tab$IntervalGroup <- interaction(tab$Context, tab$Series, tab$IntervalGroup, drop = TRUE)
  if (is.null(palette)) palette <- rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#000000"),
    length.out = nlevels(tab$Series))
  if (!is.character(palette) || length(palette) < nlevels(tab$Series) || anyNA(palette)) stop("Supply one color per series.", call. = FALSE)
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$Theta, y = .data$Estimate,
    color = .data$Series, fill = .data$Series, linetype = .data$Series, group = .data$Series)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Lower, ymax = .data$Upper,
      group = .data$IntervalGroup), alpha = .15, color = NA, na.rm = TRUE,
      show.legend = any(tab$CIEligible)) +
    ggplot2::geom_line(data = tab[!isolated, , drop = FALSE], linewidth = .7) +
    ggplot2::facet_wrap(~Context, labeller = ggplot2::as_labeller(context_labels)) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_linetype_manual(values = rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), length.out = nlevels(tab$Series))) +
    ggplot2::labs(x = "Ability (native scale)",
      y = if (x$settings$type == "probability") "Category probability" else "Information per rating",
      color = if (x$settings$type == "probability") "Category" else "Series",
      fill = if (x$settings$type == "probability") "Category" else "Series",
      linetype = if (x$settings$type == "probability") "Category" else "Series",
      alt = if (point_only) paste("Two-family GPCM", x$settings$type,
        "provisional fitted curves. Calibration intervals are unavailable for the entire plot.",
        "Contexts with only one supplied ability value retain colored points.") else
        paste("GPCM", x$settings$type, "curves with", 100*x$settings$level,
        "percent", if (x$settings$simultaneous == "none") "pointwise" else "Bonferroni-adjusted",
        "intervals;", sum(!tab$CIEligible), "unavailable grid-point intervals.",
        "Crosses mark estimates without intervals; ribbons do not bridge unavailable points.")) +
    ggplot2::theme_minimal() + ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9))
  if (any(isolated)) p <- p +
    ggplot2::geom_point(data = tab[isolated, , drop = FALSE], size = 2.5)
  if (!point_only && any(!tab$CIEligible)) p <- p +
    ggplot2::geom_point(data = tab[!tab$CIEligible, , drop = FALSE],
      ggplot2::aes(shape = "Unavailable"), color = "#333333", size = 2.5,
      show.legend = c(shape = TRUE, color = FALSE, fill = FALSE, linetype = FALSE)) +
    ggplot2::scale_shape_manual(values = c(Unavailable = 4), name = "Interval")
  attr(p, "mfrmr_plot_data") <- new_mfrm_plot_data("gpcm_curve_intervals",
    list(table = tab, settings = x$settings, title = title, subtitle = subtitle, caption = caption,
      cautions = x$cautions, contexts = x$contexts, uncertainty_note = x$anchor_note))
  p <- .mfrmr_gg_labs(p, attr(p, "mfrmr_plot_data")$data,
    x = "Ability (native scale)",
    y = if (x$settings$type == "probability") "Category probability" else "Information per rating")
  if (isTRUE(draw)) print(p)
  invisible(p)
}

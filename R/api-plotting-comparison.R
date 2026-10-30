# Paired diagnostic displays reuse the native coordinates and curve engine.

#' Compare Wright maps or category curves from two fitted models
#'
#' @param reference,comparison Two `mfrm_fit` objects. All signed differences
#'   are comparison minus reference, on the fitted coordinate scale.
#' @param type `"wright"` for distributions/locations or `"ccc"` for category
#'   probabilities at the zero additive-facet reference profile.
#' @param view `"comparison"` for paired displays or `"difference"` for
#'   matched location differences versus means, or probability differences.
#' @param labels Two distinct, nonempty display labels, reference first.
#' @param curve_groups Optional step-facet levels to compare. By default all
#'   groups are required in both fits. An RSM common scale is explicitly paired
#'   with each selected PCM/GPCM group. Non-RSM fits must share a step owner.
#' @param panel CCC layout: `"auto"` uses category panels in monochrome or
#'   with more than five categories, otherwise `"group"`. `"category"` always
#'   separates categories; `"group"` overlays categories within each group.
#' @param theta_range Finite increasing length-two predictor range for CCCs.
#' @param theta_points Integer number of grid points, at least two.
#' @param preset Existing visual preset; monochrome uses shapes/line types
#'   as well as grey tones.
#' @param show_title,show_notes Display title and explanatory subtitle.
#'   Notes and source readiness remain in the returned object.
#' @param draw Draw with the optional ggplot2 renderer. `FALSE` returns data
#'   without requiring ggplot2. Use [as_ggplot()] to customize/export a view.
#'
#' @details
#' This is a descriptive paired-model display, not an information-criterion
#' comparison or an external-software importer. Recorded centering, anchors,
#' orientation, category coding, estimator method and coordinate basis must
#' agree. These checks do not establish scale equivalence: no origin or unit
#' transformation is estimated, and no observations or model parameters are
#' refitted. Population SDs, full scale contracts and compared settings are
#' returned for inspection. Different data/assignments can affect differences.
#'
#' Wright maps use all source levels, with persons displayed as violins only
#' when at least two distinct eligible estimates exist; otherwise points are
#' used. Facet and step points are offset horizontally only. Step panels pair
#' corresponding adjacent transitions; RSM steps are repeated as references
#' across selected PCM/GPCM groups, not estimated separately. Excluded,
#' non-finite, boundary-separated and unmatched levels remain in the tables.
#' Differences are unavailable for those rows, never set to zero. No SE or CI
#' for a difference is calculated. Person distributions are distributions of
#' fitted point estimates, not posterior or population density estimates.
#'
#' CCCs reuse [plot.mfrm_fit()] probabilities, retaining GPCM slopes. Additive
#' facet effects and interactions are fixed at zero in both fits. Equal
#' original category mappings and a common predictor grid are required.
#' The maximum absolute probability difference is a grid diagnostic, with
#' its signed value, category and grid location from the first maximizing row returned;
#' it is not a continuous-domain supremum or a significance test. Use explicit
#' group selection or a larger export when many panels are needed.
#' `Category` retains the native internal code; `OriginalCategory` records the
#' score label shown in the plot. `ExpectedScoreDifference` uses the native
#' internal score coding, which need not equal the original score increments.
#'
#' @return An `mfrm_plot_data` object. `data$locations` and `data$differences`
#'   describe Wright coordinates and matches; `data$probabilities`,
#'   `data$differences` and `data$summary` describe CCCs. Both types return
#'   `group_selection`, `category_labels`, `basis`, `scale_contracts`,
#'   `fit_readiness`, `notes` and display `settings`. `source_plots` preserves
#'   the native draw-free payloads and their exclusions/interpretation metadata.
#' @seealso [plot.mfrm_fit()], [as_ggplot()], [compare_mfrm()]
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("plot_compare_mfrm", run.dontrun = TRUE)
#' \dontrun{
#' # Load one dataset for both models
#' library(mfrmr)
#' toy <- load_mfrmr_data("example_operational")
#'
#' # RSM: shared category thresholds
#' fit_rsm <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   model = "RSM"
#' )
#'
#' # PCM: separate category thresholds for each criterion
#' fit_pcm <- fit_mfrm(
#'   data = toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "MML",
#'   model = "PCM",
#'   step_facet = "Criterion"
#' )
#'
#' # Compare fitted locations; drawing requires the optional ggplot2 package
#' # If needed, install it once with install.packages("ggplot2")
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(fit_rsm, comparison = fit_pcm, labels = c("RSM", "PCM"))
#'
#'   # Run separately: compare score-category probabilities for one criterion
#'   plot(fit_rsm, comparison = fit_pcm, type = "ccc",
#'     curve_groups = "Content", labels = c("RSM", "PCM"))
#' }
#'
#' # Plot data are available even without ggplot2
#' paired <- plot(fit_rsm, comparison = fit_pcm, draw = FALSE)
#' head(paired$data$differences)
#' # Differences are PCM minus RSM; no difference SE or significance test is computed
#' }
#' @inheritSection mfrmr_visual_diagnostics Session plot defaults
#' @keywords internal
#' @section Migration in 0.2.4:
#' This renderer is internal. Use the corresponding `plot()` route in
#' [mfrmr_plot_migration]; the former exported name is no longer available.
plot_compare_mfrm <- function(reference, comparison, type = c("wright", "ccc"),
                              view = c("comparison", "difference"),
                              labels = c("Reference", "Comparison"),
                              curve_groups = NULL, panel = c("auto", "category", "group"),
                              theta_range = c(-6, 6), theta_points = 241L,
                              preset = c("standard", "publication", "compact", "monochrome"),
                              show_title = FALSE, show_notes = FALSE, draw = TRUE) {
  stop_if_product_slopes(reference, "plot(reference, comparison = comparison)")
  stop_if_product_slopes(comparison, "plot(reference, comparison = comparison)")

  if (missing(preset)) preset <- .mfrm_default_plot_preset()
  type <- match.arg(type); view <- match.arg(view); panel <- match.arg(panel)
  preset <- match.arg(preset)
  fits <- list(reference, comparison)
  if (!all(vapply(fits, inherits, logical(1), "mfrm_fit"))) stop(
    "Both inputs must be mfrm_fit objects.", call. = FALSE)
  if (!is.character(labels) || length(labels) != 2L || anyNA(labels) ||
      any(!nzchar(trimws(labels))) || anyDuplicated(labels)) stop(
    "`labels` must contain two distinct nonempty strings.", call. = FALSE)
  for (flag in list(show_title, show_notes, draw)) if (!is.logical(flag) ||
      length(flag) != 1L || is.na(flag)) stop("Display flags must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(theta_range) || length(theta_range) != 2L ||
      any(!is.finite(theta_range)) || diff(theta_range) <= 0 ||
      !is.numeric(theta_points) || length(theta_points) != 1L ||
      !is.finite(theta_points) || theta_points < 2 || theta_points != floor(theta_points)) stop(
    "Use an increasing finite theta_range and integer theta_points >= 2.", call. = FALSE)

  signatures <- lapply(fits, normalize_compare_signature)
  fields <- c("method", "facets", "rating_min", "rating_max", "score_map",
    "noncenter_facet", "dummy_facets", "positive_facets", "anchors", "group_anchors")
  basis <- data.frame(Component = fields, Matches = vapply(fields, function(k)
    same_signature_component(signatures[[1]][[k]], signatures[[2]][[k]]), logical(1)))
  scales <- lapply(fits, mfrm_fit_scale_contract)
  basis <- rbind(basis, data.frame(Component = "coordinate_basis",
    Matches = identical(scales[[1]]$CoordinateBasis, scales[[2]]$CoordinateBasis)))
  basis <- rbind(basis, data.frame(Component = "prepared_rating_range",
    Matches = identical(c(reference$prep$rating_min, reference$prep$rating_max),
      c(comparison$prep$rating_min, comparison$prep$rating_max))))
  if (!all(basis$Matches)) stop("Comparison basis differs: ",
    paste(basis$Component[!basis$Matches], collapse = ", "),
    ". No automatic scale alignment is performed.", call. = FALSE)
  if (any(vapply(signatures, function(s) length(s$facets) == 0L ||
      anyNA(c(s$rating_min, s$rating_max, s$method)), logical(1)))) stop(
    "Recorded facet, score-range and estimator metadata are required.", call. = FALSE)
  specs <- lapply(fits, build_step_curve_spec)
  if (!all(vapply(specs, function(s) s$model %in% c("RSM", "PCM", "GPCM"), logical(1)))) stop(
    "Supported models are RSM, PCM and GPCM.", call. = FALSE)
  if (any(vapply(specs, function(s) any(vapply(s$groups, function(g)
      length(g$tau) != length(s$categories) - 1L || any(!is.finite(g$tau)), logical(1))), logical(1)))) stop(
    "Each selected model requires complete finite step estimates for its retained categories.", call. = FALSE)
  category_labels <- lapply(seq_along(fits), function(i)
    wright_score_labels(fits[[i]], specs[[i]]$categories))
  if (!identical(specs[[1]]$categories, specs[[2]]$categories) ||
      !identical(category_labels[[1]]$OriginalScore, category_labels[[2]]$OriginalScore)) stop(
    "Both fits must have identical retained categories and original score mappings.", call. = FALSE)
  shared <- vapply(specs, function(s) s$model == "RSM", logical(1))
  if (!any(shared) && !identical(specs[[1]]$step_facet, specs[[2]]$step_facet)) stop(
    "Non-RSM curves must have the same step-facet owner.", call. = FALSE)
  groups <- if (all(shared)) "Common" else unique(unlist(lapply(specs[!shared], function(s) names(s$groups))))
  if (is.null(curve_groups)) curve_groups <- groups
  if (!is.character(curve_groups) || !length(curve_groups) || anyNA(curve_groups) ||
      anyDuplicated(curve_groups) || any(!curve_groups %in% groups)) stop(
    "`curve_groups` must select known, unique step groups.", call. = FALSE)
  selection <- data.frame(CurveGroup = groups, Selected = groups %in% curve_groups,
    ReferenceGroup = if (shared[1]) rep("Common", length(groups)) else groups,
    ComparisonGroup = if (shared[2]) rep("Common", length(groups)) else groups)
  selection$ReferenceAvailable <- selection$ReferenceGroup %in% names(specs[[1]]$groups)
  selection$ComparisonAvailable <- selection$ComparisonGroup %in% names(specs[[2]]$groups)
  if (any(selection$Selected & (!selection$ReferenceAvailable | !selection$ComparisonAvailable))) stop(
    "A selected curve group is missing in one fit; explicitly select groups present in both fits.", call. = FALSE)
  selected <- selection[match(curve_groups, selection$CurveGroup), , drop = FALSE]
  source <- lapply(fits, function(f) plot(f, type = type, draw = FALSE,
    top_n = Inf, show_ci = FALSE, theta_range = theta_range, theta_points = theta_points,
    preset = preset))
  tagged <- function(tables) dplyr::bind_rows(lapply(seq_along(tables), function(i) {
    t <- tables[[i]]; t$Fit <- rep(labels[i], nrow(t)); t
  }))
  notes <- tagged(lapply(source, function(p) p$data$notes))
  notes <- dplyr::bind_rows(notes, data.frame(Type = c("basis", "uncertainty", "selection"), Text = c(
      "Descriptive comparison on fitted coordinates; no alignment or refitting. Matching recorded settings does not prove scale equivalence or identical observations. Differences are comparison minus reference.",
    "No difference SE, covariance, confidence interval or significance test is computed. Source readiness is retained; overlap does not establish equivalence.",
    sprintf("Selected %d of %d step groups. RSM shared steps/curves are reused for each selected non-RSM group, not independently estimated.", nrow(selected), nrow(selection)))))
  data <- list(title = if (type == "wright") "Paired Wright maps" else "Paired category curves",
    subtitle = if (type == "wright") "Fitted locations; no scale alignment or difference intervals" else
      "Reference profile: additive facet effects and interactions fixed at zero",
    display = list(show_title = show_title, show_notes = show_notes), preset = preset,
    settings = list(type = type, view = view, labels = labels, panel = panel,
      theta_range = theta_range, theta_points = theta_points, alignment = "none"),
    group_selection = selection, category_labels = category_labels[[1]], basis = basis,
    scale_contracts = tagged(scales), fit_readiness = tagged(lapply(source, function(p) p$data$fit_readiness)),
    notes = notes, source_plots = stats::setNames(source, labels))
  if (any(vapply(source, function(p) identical(p$data$interpretation_status, "review_only"), logical(1)))) {
    data$interpretation_status <- "review_only"
    data$subtitle <- paste("REVIEW ONLY -", data$subtitle)
  }
  if (type == "ccc") {
    probabilities <- lapply(seq_along(source), function(i) {
      p <- source[[i]]$data$probabilities
      dplyr::bind_rows(lapply(seq_len(nrow(selected)), function(j) {
        t <- p[p$CurveGroup == selected[[if (i == 1L) "ReferenceGroup" else "ComparisonGroup"]][j], ]
        t$SourceCurveGroup <- t$CurveGroup; t$CurveGroup <- selected$CurveGroup[j]; t
      }))
    })
    keys <- c("CurveGroup", "Theta", "Category")
    a <- probabilities[[1]][, c(keys, "Probability", "ExpectedScore")]
    b <- probabilities[[2]][, c(keys, "Probability", "ExpectedScore")]
    differences <- dplyr::full_join(a, b, by = keys, suffix = c("_Reference", "_Comparison"))
    differences$Difference <- differences$Probability_Comparison - differences$Probability_Reference
    differences$ExpectedScoreDifference <- differences$ExpectedScore_Comparison - differences$ExpectedScore_Reference
    if (any(!is.finite(differences$Difference))) stop("CCC comparison contains unavailable probability coordinates.", call. = FALSE)
    data$probabilities <- tagged(probabilities)
    data$differences <- differences
    data$summary <- dplyr::bind_rows(lapply(curve_groups, function(g) {
      t <- differences[differences$CurveGroup == g, ]
      i <- which.max(abs(t$Difference))
      data.frame(CurveGroup = g, MaxAbsProbabilityDifference = abs(t$Difference[i]),
        Difference = t$Difference[i], Theta = t$Theta[i], Category = t$Category[i])
    }))
    for (nm in c("probabilities", "differences", "summary")) data[[nm]]$OriginalCategory <-
      data$category_labels$OriginalScore[match(as.character(data[[nm]]$Category), data$category_labels$InternalScore)]
    data$notes <- dplyr::bind_rows(data$notes, data.frame(Type = "score_coding", Text =
      "Category is the internal code; OriginalCategory is the displayed score label. ExpectedScoreDifference uses native internal score coding, not necessarily original score increments."))
  } else {
    data$locations <- tagged(lapply(seq_along(fits), function(i)
      .comparison_wright_locations(fits[[i]], source[[i]]$data, specs[[i]], selected,
        if (i == 1L) "ReferenceGroup" else "ComparisonGroup", category_labels[[i]])))
    keys <- c("Kind", "Facet", "Level")
    t <- data$locations
    keep <- c(keys, "Estimate", "SourceEstimate", "Status")
    differences <- dplyr::full_join(t[t$Fit == labels[1], keep], t[t$Fit == labels[2], keep],
      by = keys, suffix = c("_Reference", "_Comparison"))
    for (nm in c("Status_Reference", "Status_Comparison")) differences[[nm]][is.na(differences[[nm]])] <- "unmatched"
    differences$Difference <- differences$Estimate_Comparison - differences$Estimate_Reference
    differences$Mean <- (differences$Estimate_Comparison + differences$Estimate_Reference) / 2
    data$differences <- differences
    data$notes <- dplyr::bind_rows(data$notes, data.frame(Type = "availability", Text = sprintf(
      "%d of %d matched/unmatched location rows have finite differences. Boundary, non-finite and source-excluded values remain unavailable. Violins summarize eligible person point estimates only; individual identities remain in the tables.",
      sum(is.finite(differences$Difference)), nrow(differences))))
  }
  out <- new_mfrm_plot_data("paired_model_comparison", data)
  if (draw) print(as_ggplot(out))
  invisible(out)
}

.comparison_wright_locations <- function(fit, source, spec, selection, group_column, categories) {
  persons <- as.data.frame(fit$facets$person)
  facets <- as.data.frame(fit$facets$others)
  rows <- dplyr::bind_rows(
    data.frame(Kind = "Person", Facet = "Person", Level = as.character(persons$Person),
      SourceEstimate = persons$Estimate),
    data.frame(Kind = "Facet", Facet = as.character(facets$Facet), Level = as.character(facets$Level),
      SourceEstimate = facets$Estimate))
  allowed <- dplyr::bind_rows(
    data.frame(Kind = "Person", Facet = "Person", Level = as.character(source$person$Person)),
    data.frame(Kind = "Facet", Facet = source$locations$Group[source$locations$PlotType == "Facet level"],
      Level = source$locations$Label[source$locations$PlotType == "Facet level"]))
  allowed$SourceIncluded <- TRUE
  rows <- dplyr::left_join(rows, allowed, by = c("Kind", "Facet", "Level"))
  rows$Status <- ifelse(is.finite(rows$SourceEstimate),
    ifelse(is.na(rows$SourceIncluded), "source_plot_excluded", "available"), "nonfinite")
  boundary <- source$locations[source$locations$PlotType == "Facet level" & source$locations$BoundarySeparated,
    c("Group", "Label"), drop = FALSE]
  names(boundary) <- c("Facet", "Level"); boundary$Kind <- rep("Facet", nrow(boundary))
  boundary$Boundary <- rep(TRUE, nrow(boundary))
  rows <- dplyr::left_join(rows, boundary, by = c("Kind", "Facet", "Level"))
  rows$Status[!is.na(rows$Boundary)] <- "boundary_separated"
  if ("ParameterStatus" %in% names(persons)) {
    excluded_persons <- persons$Person[persons$ParameterStatus %in% c("unbounded_low", "unbounded_high")]
    rows$Status[rows$Kind == "Person" & rows$Level %in% excluded_persons] <- "boundary_separated"
  }
  rows <- rows[, c("Kind", "Facet", "Level", "SourceEstimate", "Status")]
  steps <- dplyr::bind_rows(lapply(seq_len(nrow(selection)), function(j) {
    s <- spec$step_points[spec$step_points$CurveGroup == selection[[group_column]][j], ]
    data.frame(Kind = "Step", Facet = selection$CurveGroup[j],
      Level = paste(categories$OriginalScore[s$StepIndex], categories$OriginalScore[s$StepIndex + 1L], sep = " -> "),
      SourceEstimate = s$Threshold, Status = ifelse(is.finite(s$Threshold), "available", "nonfinite"))
  }))
  rows <- dplyr::bind_rows(rows, steps)
  if (anyNA(rows[, c("Kind", "Facet", "Level")]) || anyDuplicated(rows[, c("Kind", "Facet", "Level")])) stop(
    "Wright comparison requires unique non-missing facet/level identities within each kind.", call. = FALSE)
  rows$Estimate <- ifelse(rows$Status == "available", rows$SourceEstimate, NA_real_)
  rows
}

.mfrmr_gg_comparison <- function(payload) {
  settings <- payload$settings
  labels <- settings$labels
  colors <- .plot_series_colors(labels, payload$preset)
  difference <- settings$view == "difference"
  if (settings$type == "ccc") {
    t <- if (difference) payload$differences else payload$probabilities
    categories <- as.character(payload$category_labels$InternalScore)
    t$Category <- factor(as.character(t$Category), levels = categories,
      labels = payload$category_labels$OriginalScore)
    t$CurveGroup <- factor(t$CurveGroup, levels = payload$group_selection$CurveGroup[payload$group_selection$Selected])
    panel <- settings$panel
    if (panel == "auto") panel <- if (payload$preset == "monochrome" || length(categories) > 5L) "category" else "group"
    if (panel == "group" && payload$preset == "monochrome" && length(categories) > 1L) warning(
      "Overlaid categories may be indistinguishable in monochrome; use panel = 'category'.", call. = FALSE)
    if (difference) {
      t$Value <- t$Difference
      p <- ggplot2::ggplot(t, ggplot2::aes(.data$Theta, .data$Value, colour = .data$Category)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.35) +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::coord_cartesian(ylim = c(-1, 1) * max(0.001, abs(t$Difference)))
    } else {
      t$Fit <- factor(t$Fit, levels = labels)
      p <- ggplot2::ggplot(t, ggplot2::aes(.data$Theta, .data$Probability,
        colour = .data$Category, linetype = .data$Fit)) + ggplot2::geom_line(linewidth = 0.75) +
        ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "dashed"), labels)) +
        ggplot2::coord_cartesian(ylim = c(0, 1))
    }
    p <- p + ggplot2::scale_colour_manual(values = .plot_series_colors(levels(t$Category), payload$preset))
    if (panel == "category") {
      panel_count <- nlevels(t$CurveGroup) * nlevels(t$Category)
      p <- p + ggplot2::facet_wrap(~ CurveGroup + Category,
        labeller = ggplot2::labeller(CurveGroup = ggplot2::label_wrap_gen(24),
          Category = function(x) paste("Category", x)), ncol = min(3L, ceiling(sqrt(panel_count)))) +
        ggplot2::guides(colour = "none")
    } else p <- p + ggplot2::facet_wrap(~ CurveGroup, labeller = ggplot2::label_wrap_gen(24))
    panel_count <- nlevels(t$CurveGroup) * if (panel == "category") nlevels(t$Category) else 1L
    xlab <- "Relative response propensity (logits)"
    ylab <- if (difference) paste("Probability difference:", labels[2], "-", labels[1]) else "Category probability"
  } else {
    t <- if (difference) payload$differences else payload$locations
    panels <- unique(t[, c("Kind", "Facet"), drop = FALSE])
    panels$Panel <- seq_len(nrow(panels))
    t <- dplyr::left_join(t, panels, by = c("Kind", "Facet"))
    panel_labels <- ifelse(panels$Kind == "Step", paste0("Steps: ", panels$Facet), panels$Facet)
    panel_labels <- make.unique(panel_labels)
    t$Panel <- factor(t$Panel, levels = panels$Panel, labels = panel_labels)
    if (difference) {
      t <- t[is.finite(t$Difference), ]
      p <- ggplot2::ggplot(t, ggplot2::aes(.data$Mean, .data$Difference)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
        ggplot2::geom_point(shape = 21, fill = colors[2], size = 1.8, alpha = 0.8) +
        ggplot2::facet_wrap(~ Panel, scales = "free_x", drop = FALSE, labeller = ggplot2::label_wrap_gen(24))
      empty <- panels$Panel[!panels$Panel %in% as.integer(t$Panel)]
      if (length(empty)) p <- p + ggplot2::geom_label(data = data.frame(
        Panel = factor(empty, levels = panels$Panel, labels = panel_labels)),
        ggplot2::aes(x = 0, y = 0, label = "No finite matched values"),
        inherit.aes = FALSE, size = 3)
      xlab <- "Mean of fitted locations (logits)"
      ylab <- paste("Location difference (logits):", labels[2], "-", labels[1])
    } else {
      t <- t[is.finite(t$Estimate), ]
      t$Fit <- factor(t$Fit, levels = labels)
      # ponytail: deterministic horizontal spreading; use custom ggplot labels for dense individual audits.
      t <- t[order(t$Panel, t$Fit, t$Level), ]
      t$X <- as.numeric(t$Fit) + stats::ave(seq_len(nrow(t)), interaction(t$Panel, t$Fit, drop = TRUE),
        FUN = function(z) if (length(z) == 1L) 0 else seq(-0.08, 0.08, length.out = length(z)))
      persons <- t[t$Kind == "Person", ]
      if (nrow(persons)) persons <- persons[stats::ave(persons$Estimate, persons$Fit,
        FUN = function(z) length(unique(z))) >= 2L, ]
      p <- ggplot2::ggplot(t, ggplot2::aes(.data$X, .data$Estimate, colour = .data$Fit, shape = .data$Fit)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35)
      if (nrow(persons)) p <- p + ggplot2::geom_violin(data = persons,
        ggplot2::aes(x = as.numeric(.data$Fit), group = .data$Fit, fill = .data$Fit),
        width = 0.65, alpha = 0.2, colour = "grey45", linewidth = 0.4, trim = TRUE)
      steps <- t[t$Kind == "Step", ]
      if (nrow(steps)) {
        centers <- stats::aggregate(Estimate ~ Panel + Level, steps, mean)
        yr <- range(t$Estimate)
        if (diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
        centers$LabelY <- stats::ave(centers$Estimate, centers$Panel, FUN = function(z)
          .spread_wright_label_positions(z, yr[1], yr[2]))
        steps <- dplyr::left_join(steps, centers[, c("Panel", "Level", "LabelY")], by = c("Panel", "Level"))
        p <- p + ggplot2::geom_segment(data = steps,
          ggplot2::aes(xend = 1.5, yend = .data$LabelY), colour = "grey70", linewidth = 0.3) +
          ggplot2::geom_label(data = centers, ggplot2::aes(x = 1.5, y = .data$LabelY, label = .data$Level),
            inherit.aes = FALSE, size = 2.5)
      }
      p <- p + ggplot2::geom_point(size = 1.5, alpha = 0.7) +
        ggplot2::scale_x_continuous(breaks = 1:2,
          labels = vapply(labels, function(s) paste(strwrap(s, width = 14), collapse = "\n"), character(1)),
          limits = c(0.5, 2.5)) +
        ggplot2::scale_colour_manual(values = colors) + ggplot2::scale_fill_manual(values = colors) +
        ggplot2::scale_shape_manual(values = stats::setNames(c(1, 17), labels)) +
        ggplot2::facet_wrap(~ Panel, nrow = 1, drop = FALSE, labeller = ggplot2::label_wrap_gen(24)) +
        ggplot2::guides(colour = "none", fill = "none", shape = "none")
      xlab <- NULL; ylab <- "Fitted location (logits)"
    }
    panel_count <- nrow(panels)
  }
  if (panel_count > 8L) warning("Many comparison panels; select curve_groups or enlarge the export. All selected panels are retained.", call. = FALSE)
  p <- p + .mfrmr_gg_theme() + ggplot2::theme(legend.position = "bottom",
    strip.text = ggplot2::element_text(size = 9))
  .mfrmr_gg_labs(p, payload, x = xlab, y = paste(strwrap(ylab, width = 48L), collapse = "\n"))
}

# Optional ggplot2 rendering for draw-free mfrmr plot payloads.

.require_mfrmr_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "`as_ggplot()` requires the optional `ggplot2` package. ",
      "Install ggplot2 or use the base-R renderer.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.mfrmr_gg_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
}

.mfrmr_gg_labs <- function(p, payload, x = NULL, y = NULL, fallback = NULL) {
  title <- as.character(payload$title %||% fallback %||% "")
  subtitle <- as.character(payload$subtitle %||% "")
  title <- title[!is.na(title) & nzchar(title)][1] %||% NULL
  subtitle <- subtitle[!is.na(subtitle) & nzchar(subtitle)][1] %||% NULL
  p + ggplot2::labs(title = title, subtitle = subtitle, x = x, y = y)
}

.mfrmr_gg_add_references <- function(p, reference_lines) {
  refs <- normalize_reference_lines(reference_lines)
  if (nrow(refs) == 0L) return(p)
  for (i in seq_len(nrow(refs))) {
    value <- suppressWarnings(as.numeric(refs$value[i]))
    if (!is.finite(value)) next
    lty <- as.character(refs$linetype[i] %||% "dashed")
    if (!nzchar(lty)) lty <- "dashed"
    axis <- tolower(as.character(refs$axis[i]))
    if (axis %in% c("v", "x", "vertical")) {
      p <- p + ggplot2::geom_vline(
        xintercept = value, linetype = lty, colour = "grey45", linewidth = 0.35
      )
    } else if (axis %in% c("h", "y", "horizontal")) {
      p <- p + ggplot2::geom_hline(
        yintercept = value, linetype = lty, colour = "grey45", linewidth = 0.35
      )
    }
  }
  p
}

.mfrmr_gg_component <- function(payload, component = NULL) {
  if (!is.null(component)) {
    if (length(component) != 1L || is.na(component) || !nzchar(component)) {
      stop("`component` must be a single non-empty name.", call. = FALSE)
    }
    if (!component %in% names(payload)) {
      stop(
        "`component` must be one of: ", paste(names(payload), collapse = ", "),
        call. = FALSE
      )
    }
    return(list(name = component, value = payload[[component]]))
  }
  priority <- c(
    "plot_long", "pathway_long", "information_long", "probabilities",
    "expected", "data", "locations", "table", "matrix"
  )
  for (nm in c(priority, setdiff(names(payload), priority))) {
    value <- payload[[nm]]
    if ((is.data.frame(value) && nrow(value) > 0L) ||
        (is.matrix(value) && length(value) > 0L)) {
      return(list(name = nm, value = value))
    }
  }
  stop(
    "No tabular component is available for automatic ggplot conversion. ",
    "Use `plot_data_components()` to inspect the payload.",
    call. = FALSE
  )
}

.mfrmr_gg_matrix <- function(mat, payload, flag_matrix = NULL,
                             show_values = FALSE, value_digits = 2L,
                             flag_color = "black") {
  long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(long) <- c("Row", "Column", "Value")
  long$Flag <- FALSE
  if (is.matrix(flag_matrix) && identical(dim(flag_matrix), dim(mat))) {
    long$Flag <- as.logical(as.data.frame(as.table(flag_matrix))$Freq)
    long$Flag[is.na(long$Flag)] <- FALSE
  }
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$Column, y = .data$Row, fill = .data$Value)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "firebrick", na.value = "grey90"
    ) +
    ggplot2::coord_equal() +
    .mfrmr_gg_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major = ggplot2::element_blank()
    )
  if (isTRUE(show_values)) {
    long$Label <- ifelse(
      is.finite(long$Value),
      formatC(long$Value, format = "f", digits = as.integer(value_digits)),
      ""
    )
    p <- p + ggplot2::geom_text(
      data = long, ggplot2::aes(label = .data$Label), size = 3
    )
  }
  flagged <- long[long$Flag, , drop = FALSE]
  if (nrow(flagged) > 0L) {
    p <- p + ggplot2::geom_tile(
      data = flagged, fill = NA, colour = flag_color, linewidth = 0.8
    )
  }
  .mfrmr_gg_labs(
    p, payload, x = "Group or contrast", y = "Facet level",
    fallback = "DIF/DFF heatmap"
  )
}

.mfrmr_gg_effect_summary <- function(df, payload) {
  if (!all(c("Pair", "Effect") %in% names(df))) {
    stop("Effect-summary conversion requires `Pair` and `Effect` columns.", call. = FALSE)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$Effect <- suppressWarnings(as.numeric(df$Effect))
  df <- df[is.finite(df$Effect), , drop = FALSE]
  df$PairFactor <- factor(as.character(df$Pair), levels = rev(unique(as.character(df$Pair))))
  if (!"Color" %in% names(df)) df$Color <- "#4D4D4D"
  p <- ggplot2::ggplot(df, ggplot2::aes(y = .data$PairFactor)) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0, xend = .data$Effect, yend = .data$PairFactor,
        colour = .data$Color
      ),
      linewidth = 1.1, lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$Effect, colour = .data$Color), size = 2.2
    ) +
    ggplot2::scale_colour_identity() +
    .mfrmr_gg_theme()
  if (all(c("CI_Lower", "CI_Upper") %in% names(df))) {
    ci <- df[is.finite(df$CI_Lower) & is.finite(df$CI_Upper), , drop = FALSE]
    if (nrow(ci) > 0L) {
      p <- p + ggplot2::geom_segment(
        data = ci,
        ggplot2::aes(
          x = .data$CI_Lower, xend = .data$CI_Upper,
          y = .data$PairFactor, yend = .data$PairFactor
        ),
        inherit.aes = FALSE, colour = "grey35", linewidth = 0.45
      )
    }
  }
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  .mfrmr_gg_labs(
    p, payload,
    x = payload$settings$effect_axis_label %||% "Effect", y = NULL,
    fallback = "Differential functioning summary"
  )
}

.mfrmr_gg_wright <- function(payload) {
  loc <- as.data.frame(payload$locations %||% data.frame(), stringsAsFactors = FALSE)
  required <- c("PlotType", "Group", "Label", "Estimate", "X")
  if (!all(required %in% names(loc))) {
    stop("Wright-map conversion requires locations with PlotType, Group, Label, Estimate, and X.", call. = FALSE)
  }
  loc$Estimate <- suppressWarnings(as.numeric(loc$Estimate))
  loc$X <- suppressWarnings(as.numeric(loc$X))
  loc <- loc[is.finite(loc$Estimate) & is.finite(loc$X), , drop = FALSE]
  groups <- as.character(payload$group_levels %||% unique(loc$Group))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  hist_obj <- payload$person_hist
  hist_df <- data.frame()
  if (is.list(hist_obj) && length(hist_obj$counts) > 0L &&
      length(hist_obj$breaks) == length(hist_obj$counts) + 1L) {
    counts <- suppressWarnings(as.numeric(hist_obj$counts))
    denom <- max(counts, na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- 1
    hist_df <- data.frame(
      xmin = -0.85 * counts / denom,
      xmax = 0,
      ymin = utils::head(hist_obj$breaks, -1L),
      ymax = utils::tail(hist_obj$breaks, -1L)
    )
  }
  p <- ggplot2::ggplot()
  y_limits <- suppressWarnings(as.numeric(
    payload$y_range %||% range(loc$Estimate, finite = TRUE)
  ))
  if (length(y_limits) != 2L || !all(is.finite(y_limits)) || diff(y_limits) <= 0) {
    y_limits <- range(loc$Estimate, finite = TRUE) + c(-0.5, 0.5)
  }
  odd_groups <- seq_along(groups)[seq_along(groups) %% 2L == 1L]
  if (length(odd_groups) > 0L) {
    shade <- data.frame(
      xmin = odd_groups - 0.5,
      xmax = odd_groups + 0.5,
      ymin = y_limits[1],
      ymax = y_limits[2]
    )
    p <- p + ggplot2::geom_rect(
      data = shade,
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax
      ),
      inherit.aes = FALSE, fill = "#EEF2F7", colour = NA, alpha = 0.55
    )
  }
  if (nrow(hist_df) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = hist_df,
      ggplot2::aes(
        xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax
      ),
      inherit.aes = FALSE, fill = "#94A3B8", colour = "white", alpha = 0.7
    )
  }
  density <- as.data.frame(payload$group %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(density) > 0L && all(c("Group", "Theta", "Density") %in% names(density))) {
    density$Theta <- suppressWarnings(as.numeric(density$Theta))
    density$Density <- suppressWarnings(as.numeric(density$Density))
    density_max <- max(density$Density, na.rm = TRUE)
    if (is.finite(density_max) && density_max > 0) {
      density$DensityX <- -0.85 * density$Density / density_max
      p <- p + ggplot2::geom_path(
        data = density,
        ggplot2::aes(
          x = .data$DensityX, y = .data$Theta,
          colour = .data$Group, group = .data$Group
        ),
        inherit.aes = FALSE, linewidth = 0.7
      )
    }
  }
  summary_df <- as.data.frame(payload$group_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(summary_df) > 0L &&
      all(c("XBase", "Min", "Max", "Q1", "Q3", "Median", "PlotType") %in% names(summary_df))) {
    p <- p +
      ggplot2::geom_segment(
        data = summary_df,
        ggplot2::aes(
          x = .data$XBase, xend = .data$XBase,
          y = .data$Min, yend = .data$Max, colour = .data$PlotType
        ),
        inherit.aes = FALSE, linewidth = 0.55, alpha = 0.5
      ) +
      ggplot2::geom_segment(
        data = summary_df,
        ggplot2::aes(
          x = .data$XBase, xend = .data$XBase,
          y = .data$Q1, yend = .data$Q3, colour = .data$PlotType
        ),
        inherit.aes = FALSE, linewidth = 2.2, alpha = 0.8
      ) +
      ggplot2::geom_segment(
        data = summary_df,
        ggplot2::aes(
          x = .data$XBase - 0.12, xend = .data$XBase + 0.12,
          y = .data$Median, yend = .data$Median
        ),
        inherit.aes = FALSE, colour = "grey20", linewidth = 0.45
      )
  }
  person_stats <- as.data.frame(payload$person_stats %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(person_stats) > 0L && all(c("Mean", "Median") %in% names(person_stats))) {
    stat_lines <- data.frame(
      Statistic = c("Person mean", "Person median"),
      Value = c(person_stats$Mean[1], person_stats$Median[1])
    )
    stat_lines <- stat_lines[is.finite(stat_lines$Value), , drop = FALSE]
    if (nrow(stat_lines) > 0L) {
      p <- p + ggplot2::geom_segment(
        data = stat_lines,
        ggplot2::aes(
          x = -0.85, xend = 0, y = .data$Value, yend = .data$Value,
          linetype = .data$Statistic
        ),
        inherit.aes = FALSE, colour = "grey25", linewidth = 0.4
      )
    }
  }
  if (all(c("CI_Lower", "CI_Upper") %in% names(loc))) {
    ci <- loc[is.finite(loc$CI_Lower) & is.finite(loc$CI_Upper), , drop = FALSE]
    if (nrow(ci) > 0L) {
      p <- p + ggplot2::geom_segment(
        data = ci,
        ggplot2::aes(
          x = .data$X, xend = .data$X,
          y = .data$CI_Lower, yend = .data$CI_Upper,
          colour = .data$PlotType
        ),
        inherit.aes = FALSE, linewidth = 0.4, alpha = 0.65
      )
    }
  }
  p <- p +
    ggplot2::geom_point(
      data = loc,
      ggplot2::aes(
        x = .data$X, y = .data$Estimate,
        colour = .data$PlotType, shape = .data$PlotType
      ),
      size = 2.1, alpha = 0.9
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::scale_shape_manual(values = c(`Facet level` = 16, `Step threshold` = 17)) +
    ggplot2::scale_x_continuous(
      breaks = c(-0.425, seq_along(groups)),
      labels = c("Persons", truncate_axis_label(groups, width = 16L)),
      limits = c(-0.95, max(1, length(groups)) + 0.75)
    ) +
    .mfrmr_gg_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    )
  labels <- as.data.frame(payload$label_points %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(labels) > 0L && all(c("X", "Estimate", "Label") %in% names(labels))) {
    labels$LabelShort <- truncate_axis_label(labels$Label, width = 14L)
    p <- p + ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = .data$X, y = .data$Estimate, label = .data$LabelShort),
      inherit.aes = FALSE, hjust = -0.1, size = 2.5, check_overlap = TRUE
    )
  }
  .mfrmr_gg_labs(p, payload, x = NULL, y = "Logit scale", fallback = "Wright map")
}

.mfrmr_gg_expected_pathway <- function(payload) {
  expected <- as.data.frame(payload$expected %||% data.frame(), stringsAsFactors = FALSE)
  if (!all(c("Theta", "ExpectedScore", "CurveGroup") %in% names(expected))) {
    stop("Expected-score pathway conversion requires Theta, ExpectedScore, and CurveGroup.", call. = FALSE)
  }
  expected$Theta <- suppressWarnings(as.numeric(expected$Theta))
  expected$ExpectedScore <- suppressWarnings(as.numeric(expected$ExpectedScore))
  expected <- expected[is.finite(expected$Theta) & is.finite(expected$ExpectedScore), , drop = FALSE]
  groups <- unique(as.character(expected$CurveGroup))
  y_rng <- range(expected$ExpectedScore, finite = TRUE)
  dominance <- as.data.frame(payload$dominance_regions %||% data.frame(), stringsAsFactors = FALSE)
  use_strips <- nrow(dominance) > 0L && length(groups) <= 8L && all(
    c("CurveGroup", "Category", "ThetaStart", "ThetaEnd", "ThetaMid") %in% names(dominance)
  )
  p <- ggplot2::ggplot()
  y_limits <- y_rng
  if (use_strips) {
    band_h <- max(diff(y_rng) * 0.035, 0.08)
    band_gap <- band_h * 0.35
    index <- stats::setNames(seq_along(groups), groups)
    dominance$Band <- unname(index[as.character(dominance$CurveGroup)])
    dominance$ymax <- y_rng[1] - band_h * 0.45 - (dominance$Band - 1L) * (band_h + band_gap)
    dominance$ymin <- dominance$ymax - band_h
    dominance$ymid <- (dominance$ymin + dominance$ymax) / 2
    y_limits[1] <- min(dominance$ymin, na.rm = TRUE) - band_h
    p <- p +
      ggplot2::geom_rect(
        data = dominance,
        ggplot2::aes(
          xmin = .data$ThetaStart, xmax = .data$ThetaEnd,
          ymin = .data$ymin, ymax = .data$ymax, fill = .data$CurveGroup
        ),
        inherit.aes = FALSE, alpha = 0.16, colour = "white", linewidth = 0.2
      ) +
      ggplot2::geom_text(
        data = dominance,
        ggplot2::aes(
          x = .data$ThetaMid, y = .data$ymid,
          label = .data$Category, colour = .data$CurveGroup
        ),
        inherit.aes = FALSE, size = 2.3, show.legend = FALSE,
        check_overlap = TRUE
      )
  }
  steps <- as.data.frame(payload$steps %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(steps) > 0L && all(c("Threshold", "PathY", "ThresholdLabel", "CurveGroup") %in% names(steps))) {
    steps$LabelShort <- truncate_axis_label(steps$ThresholdLabel, width = 6L)
    p <- p +
      ggplot2::geom_vline(
        data = steps,
        ggplot2::aes(xintercept = .data$Threshold, colour = .data$CurveGroup),
        inherit.aes = FALSE, linetype = "dotted", alpha = 0.25,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = steps,
        ggplot2::aes(x = .data$Threshold, y = .data$PathY, colour = .data$CurveGroup),
        inherit.aes = FALSE, shape = 18, size = 2, show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = steps,
        ggplot2::aes(
          x = .data$Threshold, y = .data$PathY,
          label = .data$LabelShort, colour = .data$CurveGroup
        ),
        inherit.aes = FALSE, vjust = -0.7, size = 2.3,
        show.legend = FALSE, check_overlap = TRUE
      )
  }
  p <- p +
    ggplot2::geom_line(
      data = expected,
      ggplot2::aes(
        x = .data$Theta, y = .data$ExpectedScore,
        colour = .data$CurveGroup, group = .data$CurveGroup
      ),
      linewidth = 0.85
    ) +
    ggplot2::coord_cartesian(ylim = y_limits, clip = "off") +
    .mfrmr_gg_theme() +
    ggplot2::theme(legend.position = "bottom")
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  .mfrmr_gg_labs(
    p, payload, x = "Theta / logit", y = "Expected score",
    fallback = "Expected-score pathway"
  )
}

.mfrmr_gg_fit_pathway <- function(payload) {
  tbl <- as.data.frame(payload$table %||% data.frame(), stringsAsFactors = FALSE)
  required <- c("Facet", "Level", "FitValue", "Measure", "ElementType")
  if (!all(required %in% names(tbl))) {
    stop("Fit-pathway conversion requires Facet, Level, FitValue, Measure, and ElementType.", call. = FALSE)
  }
  tbl <- tbl[is.finite(tbl$FitValue) & is.finite(tbl$Measure), , drop = FALSE]
  if (!"LabelText" %in% names(tbl)) tbl$LabelText <- ""
  p <- ggplot2::ggplot()
  if (identical(payload$fit_scale, "mnsq")) {
    p <- p + ggplot2::annotate(
      "rect", xmin = payload$fit_range[1], xmax = payload$fit_range[2],
      ymin = -Inf, ymax = Inf, fill = "#E8F0FE", alpha = 0.45
    )
  } else {
    p <- p + ggplot2::annotate(
      "rect", xmin = -payload$zstd_cut, xmax = payload$zstd_cut,
      ymin = -Inf, ymax = Inf, fill = "#E8F0FE", alpha = 0.45
    )
  }
  if (isTRUE(payload$show_ci) && all(c("CI_Lower", "CI_Upper") %in% names(tbl))) {
    ci <- tbl[is.finite(tbl$CI_Lower) & is.finite(tbl$CI_Upper), , drop = FALSE]
    if (nrow(ci) > 0L) {
      p <- p + ggplot2::geom_errorbar(
        data = ci,
        ggplot2::aes(
          x = .data$FitValue, ymin = .data$CI_Lower, ymax = .data$CI_Upper,
          colour = .data$Facet
        ),
        inherit.aes = FALSE, width = 0.02, alpha = 0.6
      )
    }
  }
  p <- p +
    ggplot2::geom_point(
      data = tbl,
      ggplot2::aes(
        x = .data$FitValue, y = .data$Measure,
        colour = .data$Facet, shape = .data$ElementType
      ),
      size = 2.2, alpha = 0.85
    ) +
    ggplot2::scale_shape_manual(values = c(`Facet level` = 16, Person = 15)) +
    .mfrmr_gg_theme()
  labels <- tbl[nzchar(as.character(tbl$LabelText)), , drop = FALSE]
  if (nrow(labels) > 0L) {
    p <- p + ggplot2::geom_text(
      data = labels,
      ggplot2::aes(
        x = .data$FitValue, y = .data$Measure,
        label = .data$LabelText, colour = .data$Facet
      ),
      inherit.aes = FALSE, hjust = -0.1, size = 2.5,
      show.legend = FALSE, check_overlap = TRUE
    )
  }
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  if (identical(payload$panel, "facet") && "Panel" %in% names(tbl)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula("~ Panel"), scales = "free_y")
  }
  x_label <- if (identical(payload$fit_scale, "mnsq")) {
    paste0(payload$fit_stat, " MnSq")
  } else {
    paste0(payload$fit_stat, " ZSTD (", payload$zstd_method, " df)")
  }
  .mfrmr_gg_labs(p, payload, x = x_label, y = "Measure (logits)", fallback = "Fit pathway")
}

.mfrmr_gg_bubble <- function(payload) {
  tbl <- as.data.frame(payload$table %||% data.frame(), stringsAsFactors = FALSE)
  view <- as.character(payload$view %||% "measure")
  stat <- as.character(payload$fit_stat %||% "Infit")
  if (identical(view, "infit_outfit")) {
    tbl$.x <- tbl$Infit
    tbl$.y <- tbl$Outfit
    xlab <- "Infit MnSq"
    ylab <- "Outfit MnSq"
  } else {
    tbl$.x <- tbl$Estimate
    tbl$.y <- tbl[[stat]]
    xlab <- "Measure (logits)"
    ylab <- paste0(stat, " MnSq")
  }
  tbl <- tbl[is.finite(tbl$.x) & is.finite(tbl$.y), , drop = FALSE]
  p <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = .data$.x, y = .data$.y, colour = .data$Facet)
  ) + ggplot2::geom_point(size = 2.6, alpha = 0.7) + .mfrmr_gg_theme()
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  .mfrmr_gg_labs(p, payload, x = xlab, y = ylab, fallback = "Bubble chart")
}

.mfrmr_gg_ccc <- function(payload,
                          slope_aes = c("none", "linewidth", "alpha", "colour"),
                          facet_by = c("auto", "none", "curve_group"),
                          show_overlay = TRUE) {
  slope_aes <- match.arg(slope_aes)
  facet_by <- match.arg(facet_by)
  prob <- as.data.frame(payload$probabilities %||% data.frame(), stringsAsFactors = FALSE)
  required <- c("Theta", "Probability", "Category", "CurveGroup")
  if (!all(required %in% names(prob))) {
    stop("CCC conversion requires Theta, Probability, Category, and CurveGroup.", call. = FALSE)
  }
  if (!"Slope" %in% names(prob)) prob$Slope <- 1
  prob$Slope <- suppressWarnings(as.numeric(prob$Slope))
  prob$Trace <- paste(prob$CurveGroup, prob$Category, sep = " | ")
  mapping <- switch(
    slope_aes,
    linewidth = ggplot2::aes(
      x = .data$Theta, y = .data$Probability, colour = .data$Category,
      group = .data$Trace, linewidth = .data$Slope
    ),
    alpha = ggplot2::aes(
      x = .data$Theta, y = .data$Probability, colour = .data$Category,
      group = .data$Trace, alpha = .data$Slope
    ),
    colour = ggplot2::aes(
      x = .data$Theta, y = .data$Probability, colour = .data$Slope,
      linetype = .data$Category, group = .data$Trace
    ),
    none = ggplot2::aes(
      x = .data$Theta, y = .data$Probability, colour = .data$Category,
      group = .data$Trace
    )
  )
  p <- ggplot2::ggplot(prob, mapping)
  p <- if (identical(slope_aes, "linewidth")) {
    p + ggplot2::geom_line()
  } else {
    p + ggplot2::geom_line(linewidth = 0.85)
  }
  p <- p + ggplot2::coord_cartesian(ylim = c(0, 1)) +
    .mfrmr_gg_theme() +
    ggplot2::theme(legend.position = "bottom")
  group_n <- length(unique(as.character(prob$CurveGroup)))
  if (identical(facet_by, "auto")) facet_by <- if (group_n > 1L) "curve_group" else "none"
  if (identical(facet_by, "curve_group") && group_n > 1L) {
    p <- p + ggplot2::facet_wrap(stats::as.formula("~ CurveGroup"))
  }
  overlay <- as.data.frame(payload$overlay %||% data.frame(), stringsAsFactors = FALSE)
  if (isTRUE(show_overlay) && nrow(overlay) > 0L &&
      all(c("Theta", "Proportion", "Category") %in% names(overlay))) {
    p <- p + ggplot2::geom_point(
      data = overlay,
      ggplot2::aes(
        x = .data$Theta, y = .data$Proportion, fill = .data$Category
      ),
      inherit.aes = FALSE, shape = 21, alpha = 0.75
    )
  }
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  .mfrmr_gg_labs(
    p, payload, x = "Theta / logit", y = "Category probability",
    fallback = "Category characteristic curves"
  )
}

.mfrmr_gg_tabular <- function(component, payload) {
  value <- component$value
  if (is.matrix(value)) return(.mfrmr_gg_matrix(value, payload))
  df <- as.data.frame(value, stringsAsFactors = FALSE)
  if (all(c("Pair", "Effect") %in% names(df))) {
    return(.mfrmr_gg_effect_summary(df, payload))
  }
  if (all(c("Theta", "Value") %in% names(df))) {
    df$Series <- if ("Series" %in% names(df)) as.character(df$Series) else "Series"
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$Theta, y = .data$Value,
        colour = .data$Series, group = .data$Series
      )
    ) + ggplot2::geom_line() + .mfrmr_gg_theme()
    return(.mfrmr_gg_labs(p, payload, x = "Theta / logit", y = "Value"))
  }
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  label_cols <- names(df)[vapply(df, function(z) is.character(z) || is.factor(z), logical(1))]
  if (length(numeric_cols) > 0L && length(label_cols) > 0L) {
    dat <- data.frame(
      Value = suppressWarnings(as.numeric(df[[numeric_cols[1]]])),
      Label = as.character(df[[label_cols[1]]])
    )
    dat <- dat[is.finite(dat$Value) & nzchar(dat$Label), , drop = FALSE]
    dat$Label <- factor(dat$Label, levels = rev(unique(dat$Label)))
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Value, y = .data$Label)) +
      ggplot2::geom_col(fill = "#2B6CB0") + .mfrmr_gg_theme()
    return(.mfrmr_gg_labs(p, payload, x = numeric_cols[1], y = NULL))
  }
  if (length(numeric_cols) >= 2L) {
    dat <- data.frame(X = df[[numeric_cols[1]]], Y = df[[numeric_cols[2]]])
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$X, y = .data$Y)) +
      ggplot2::geom_point() + .mfrmr_gg_theme()
    return(.mfrmr_gg_labs(p, payload, x = numeric_cols[1], y = numeric_cols[2]))
  }
  stop("The selected component is not suitable for automatic ggplot conversion.", call. = FALSE)
}

.mfrmr_gg_simulation_scan <- function(payload) {
  .require_mfrmr_ggplot2()
  df <- as.data.frame(payload$data %||% data.frame(), stringsAsFactors = FALSE)
  x_var <- as.character(payload$x_var %||% "")
  if (!nzchar(x_var) || !x_var %in% names(df) || !"y" %in% names(df)) {
    stop("Simulation plot conversion requires `data`, `x_var`, and a `y` column.", call. = FALSE)
  }
  df$.mfrmr_x <- suppressWarnings(as.numeric(df[[x_var]]))
  df$.mfrmr_y <- suppressWarnings(as.numeric(df$y))
  df$.mfrmr_group <- if ("group" %in% names(df)) as.character(df$group) else "All designs"
  df <- df[is.finite(df$.mfrmr_x) & is.finite(df$.mfrmr_y), , drop = FALSE]
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$.mfrmr_x, y = .data$.mfrmr_y,
      colour = .data$.mfrmr_group, group = .data$.mfrmr_group
    )
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2) +
    .mfrmr_gg_theme() +
    ggplot2::theme(legend.position = "bottom")
  payload$title <- payload$title %||% if (identical(payload$plot, "signal_detection")) {
    "Signal-detection simulation"
  } else {
    "Design simulation"
  }
  payload$subtitle <- payload$subtitle %||% payload$interpretation_note %||% NULL
  .mfrmr_gg_labs(
    p,
    payload,
    x = payload$x_label %||% x_var,
    y = payload$display_metric %||% payload$metric_col %||% payload$metric %||% "Value"
  )
}

#' Convert draw-free mfrmr plot data to ggplot2
#'
#' `as_ggplot()` is an optional renderer for an `mfrm_plot_data` object or an
#' object whose `plot()` method supports `draw = FALSE`. Base graphics remain
#' the default; the returned `ggplot` can be restyled or composed downstream.
#'
#' Dedicated conversions are provided for Wright maps, theta-to-expected-score
#' pathways, fit-statistic-to-measure pathways, category characteristic curves,
#' bubble charts, and DIF/DFF summaries and heatmaps. Other draw-free payloads
#' use a conservative tabular fallback; inspect [plot_data_components()] when
#' automatic inference is not appropriate.
#'
#' @param x An `mfrm_plot_data` object, or an mfrmr object with a draw-free
#'   plot method.
#' @param type Optional plot type passed to `plot()` for a non-plot-data input.
#' @param component Optional tabular payload component to convert.
#' @param ... Arguments passed to the draw-free plot method. CCC conversion
#'   additionally accepts `slope_aes`, `facet_by`, and `show_overlay`.
#'
#' @return A `ggplot2` plot object.
#' @examples
#' \donttest{
#' fit <- fit_mfrm(load_mfrmr_data("example_core"), "Person",
#'                 c("Rater", "Criterion"), "Score", maxit = 30)
#' as_ggplot(fit, type = "wright")
#' as_ggplot(fit, type = "fit_pathway", include_person = TRUE)
#' as_ggplot(plot(fit, type = "ccc", draw = FALSE))
#' }
#' @export
as_ggplot <- function(x, type = NULL, component = NULL, ...) {
  UseMethod("as_ggplot")
}

#' @rdname as_ggplot
#' @export
as_ggplot.default <- function(x, type = NULL, component = NULL, ...) {
  dots <- list(...)
  conversion_names <- c("slope_aes", "facet_by", "show_overlay")
  dot_names <- names(dots) %||% rep("", length(dots))
  conversion <- nzchar(dot_names) & dot_names %in% conversion_names
  plot_obj <- do.call(
    as_mfrm_plot_data_object,
    c(list(x = x, type = type), dots[!conversion])
  )
  do.call(
    as_ggplot,
    c(list(x = plot_obj, component = component), dots[conversion])
  )
}

#' @rdname as_ggplot
#' @export
as_ggplot.mfrm_design_evaluation <- function(x, type = NULL, component = NULL, ...) {
  payload <- plot(x, draw = FALSE, ...)
  .mfrmr_gg_simulation_scan(payload)
}

#' @rdname as_ggplot
#' @export
as_ggplot.mfrm_design_evaluation_plot_data <- function(x, type = NULL,
                                                       component = NULL, ...) {
  payload <- as.list(x)
  if (!is.null(component)) {
    .require_mfrmr_ggplot2()
    return(.mfrmr_gg_tabular(.mfrmr_gg_component(payload, component), payload))
  }
  .mfrmr_gg_simulation_scan(payload)
}

#' @rdname as_ggplot
#' @export
as_ggplot.mfrm_signal_detection <- function(x, type = NULL, component = NULL, ...) {
  payload <- plot(x, draw = FALSE, ...)
  .mfrmr_gg_simulation_scan(payload)
}

#' @rdname as_ggplot
#' @export
as_ggplot.mfrm_signal_detection_plot_data <- function(x, type = NULL,
                                                      component = NULL, ...) {
  payload <- as.list(x)
  if (!is.null(component)) {
    .require_mfrmr_ggplot2()
    return(.mfrmr_gg_tabular(.mfrmr_gg_component(payload, component), payload))
  }
  .mfrmr_gg_simulation_scan(payload)
}

#' @rdname as_ggplot
#' @export
as_ggplot.mfrm_plot_data <- function(x, type = NULL, component = NULL, ...) {
  .require_mfrmr_ggplot2()
  payload <- x$data %||% list()
  dots <- list(...)
  if (is.null(component) && identical(x$name, "wright_map")) {
    return(.mfrmr_gg_wright(payload))
  }
  if (is.null(component) && identical(x$name, "pathway_map")) {
    return(.mfrmr_gg_expected_pathway(payload))
  }
  if (is.null(component) && identical(x$name, "fit_pathway")) {
    return(.mfrmr_gg_fit_pathway(payload))
  }
  if (is.null(component) && identical(x$name, "bubble")) {
    return(.mfrmr_gg_bubble(payload))
  }
  if (is.null(component) && x$name %in% c(
    "category_characteristic_curves", "category_characteristic_curves_overlay"
  )) {
    slope_aes <- match.arg(
      as.character(dots$slope_aes %||% "none"),
      c("none", "linewidth", "alpha", "colour", "color")
    )
    if (identical(slope_aes, "color")) slope_aes <- "colour"
    facet_by <- match.arg(
      as.character(dots$facet_by %||% "auto"),
      c("auto", "none", "curve_group")
    )
    return(.mfrmr_gg_ccc(
      payload,
      slope_aes = slope_aes,
      facet_by = facet_by,
      show_overlay = isTRUE(dots$show_overlay %||% TRUE)
    ))
  }
  if (is.null(component) && identical(x$name, "dif_summary") && is.data.frame(payload$data)) {
    return(.mfrmr_gg_effect_summary(payload$data, payload))
  }
  if (is.null(component) && identical(x$name, "dif_heatmap") && is.matrix(payload$matrix)) {
    settings <- payload$settings %||% list()
    return(.mfrmr_gg_matrix(
      payload$matrix,
      payload,
      flag_matrix = payload$flag_matrix,
      show_values = isTRUE(settings$show_values),
      value_digits = settings$value_digits %||% 2L,
      flag_color = settings$flag_color %||% "black"
    ))
  }
  .mfrmr_gg_tabular(.mfrmr_gg_component(payload, component), payload)
}

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
  labels <- lapply(list(title = payload$title %||% fallback,
                       subtitle = payload$subtitle, caption = payload$caption), function(value) {
    value <- as.character(value)
    value <- value[!is.na(value) & nzchar(value)]
    if (!length(value)) return(NULL)
    # ponytail: wrap at 72 columns; use labs() for unusually narrow export sizes.
    paste(unlist(lapply(strsplit(value[1L], "\n", fixed = TRUE)[[1L]],
                        strwrap, width = 72L)), collapse = "\n")
  })
  display <- payload[["display"]]
  if (!is.list(display) || is.data.frame(display)) display <- list()
  if (!isTRUE(display$show_title)) labels$title <- NULL
  if (!isTRUE(display$show_notes)) {
    labels$subtitle <- NULL
    labels$caption <- NULL
  }
  p <- p + ggplot2::labs(
    title = labels$title, subtitle = labels$subtitle, caption = labels$caption, x = x, y = y
  )
  attr(p, "mfrmr_notes") <- .mfrm_plot_annotation_notes(payload)
  attr(p, "mfrmr_plot_data") <- new_mfrm_plot_data(payload$plot_name %||% "plot", payload)
  p
}

.mfrmr_gg_fit_measure_bands <- function(x) {
  payload <- x$data
  df <- payload$table
  style <- resolve_plot_preset(payload$preset %||% "standard")
  palette <- c(underfit = style$warn, overfit = style$accent_secondary,
    mixed = style$fail, within_band = style$success, not_available = style$neutral)
  fixed <- df[["Fixed"]] %in% TRUE
  if (length(fixed) != nrow(df)) fixed <- rep(FALSE, nrow(df))
  shown <- is.finite(df$Measure)
  interval <- shown & !fixed & is.finite(df$CI_Lower) & is.finite(df$CI_Upper)
  df$Row <- seq_len(nrow(df))
  df$Shape <- ifelse(fixed, 5, ifelse(interval, 21, 4))
  df$Colour <- unname(palette[match(df$FitStatus, names(palette))])
  df$Colour[is.na(df$Colour)] <- style$accent_primary
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Measure, y = .data$Row)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = style$neutral) +
    ggplot2::geom_segment(data = df[interval, , drop = FALSE],
      ggplot2::aes(x = .data$CI_Lower, xend = .data$CI_Upper, yend = .data$Row),
      colour = style$neutral, linewidth = .6) +
    ggplot2::geom_point(data = df[shown, , drop = FALSE],
      ggplot2::aes(shape = .data$Shape, colour = .data$Colour, fill = .data$Colour), size = 2) +
    ggplot2::scale_shape_identity() + ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_reverse(breaks = df$Row[shown],
      labels = paste(df$Facet[shown], df$Level[shown], sep = ": ")) +
    .mfrmr_gg_theme()
  p <- .mfrmr_gg_labs(p, payload, x = "Measure (logits)", y = NULL)
  attr(p, "mfrmr_plot_data") <- x
  p
}

.mfrmr_gg_severity <- function(x) {
  payload <- x$data
  df <- payload$data
  geometry <- .mfrm_severity_geometry(df)
  style <- resolve_plot_preset(payload$preset %||% "standard")
  df$Row <- seq_len(nrow(df))
  df$Shape <- geometry$shape
  show_bands <- payload$show_bands %||% any(payload$legend$role == "band")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Estimate, y = .data$Row))
  if (isTRUE(show_bands)) {
    # Older saved plots keep their original guide colors and data labels.
    band_colors <- payload$band_colors %||% c(style$success, style$warn, style$warn)
    p <- p + ggplot2::annotate("rect", xmin = c(-.5, -1, .5), xmax = c(.5, -.5, 1),
      ymin = -Inf, ymax = Inf, fill = band_colors,
      alpha = c(.10, .08, .08))
  }
  p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = style$neutral) +
    ggplot2::geom_segment(data = df[geometry$interval, , drop = FALSE],
      ggplot2::aes(x = .data$CI_Lower, xend = .data$CI_Upper, yend = .data$Row),
      colour = style$accent_primary, linewidth = .6) +
    ggplot2::geom_point(data = df[geometry$shown, , drop = FALSE],
      ggplot2::aes(shape = .data$Shape), colour = style$accent_primary, size = 2) +
    ggplot2::geom_text(data = df[!geometry$shown, , drop = FALSE],
      ggplot2::aes(x = mean(geometry$xrange)), label = "Unavailable",
      colour = style$neutral, size = 3) +
    ggplot2::scale_shape_identity() +
    ggplot2::scale_y_continuous(breaks = df$Row, labels = df$Level) +
    ggplot2::scale_x_continuous(limits = geometry$xrange) + .mfrmr_gg_theme()
  p <- .mfrmr_gg_labs(p, payload, x = payload$axis_label %||% "Severity (logit)", y = NULL)
  attr(p, "mfrmr_plot_data") <- x
  p
}

.mfrmr_gg_shrinkage <- function(x) {
  payload <- x$data
  validate_shrinkage_output(list(shrinkage_report = payload$shrinkage_report))
  df <- payload$data
  style <- resolve_plot_preset(payload$preset %||% "standard")
  df <- df[order(df$Facet, -abs(df$ShrinkageFactor)), , drop = FALSE]
  df$Row <- seq_len(nrow(df))
  if (nrow(df) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No shrinkage applied.") +
      ggplot2::theme_void()
  } else {
    points <- rbind(
      data.frame(Row = df$Row, Value = df$Estimate, Role = "Original"),
      data.frame(Row = df$Row, Value = df$ShrunkEstimate, Role = "Shrunk"))
    arrows <- df[is.finite(df$Estimate) & is.finite(df$ShrunkEstimate), , drop = FALSE]
    p <- ggplot2::ggplot() +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = style$neutral) +
      ggplot2::geom_segment(data = arrows,
        ggplot2::aes(x = .data$Estimate, xend = .data$ShrunkEstimate,
          y = .data$Row, yend = .data$Row), colour = style$neutral,
        arrow = grid::arrow(length = grid::unit(.06, "inches")))
    if (isTRUE(payload$show_ci)) {
      # Render the saved endpoints; conversion does not recompute uncertainty.
      bands <- rbind(
        data.frame(Row = df$Row - .12, Lower = df$CI_Lower,
          Upper = df$CI_Upper, Role = "Original"),
        data.frame(Row = df$Row + .12, Lower = df$ShrunkCI_Lower,
          Upper = df$ShrunkCI_Upper, Role = "Shrunk"))
      bands <- bands[is.finite(bands$Lower) & is.finite(bands$Upper), , drop = FALSE]
      p <- p + ggplot2::geom_segment(data = bands,
        ggplot2::aes(x = .data$Lower, xend = .data$Upper, y = .data$Row,
          yend = .data$Row, colour = .data$Role), show.legend = FALSE)
    }
    p <- p +
      ggplot2::geom_point(data = points[is.finite(points$Value), , drop = FALSE],
        ggplot2::aes(x = .data$Value, y = .data$Row,
          colour = .data$Role, shape = .data$Role, fill = .data$Role), size = 2.5) +
      ggplot2::scale_colour_manual(values = c(Original = style$accent_primary,
        Shrunk = style$accent_tertiary), name = NULL) +
      ggplot2::scale_shape_manual(values = c(Original = 16, Shrunk = 21), name = NULL) +
      ggplot2::scale_fill_manual(values = c(Original = style$accent_primary,
        Shrunk = "white"), name = NULL) +
      ggplot2::scale_y_reverse(breaks = df$Row,
        labels = paste(df$Facet, df$Level, sep = " / ")) +
      .mfrmr_gg_theme()
  }
  p <- .mfrmr_gg_labs(p, payload, x = "logit", y = NULL)
  attr(p, "mfrmr_plot_data") <- x
  p
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
  residual <- identical(payload$plot_name, "local_dependence_heatmap")
  halo <- identical(payload$plot_name, "halo_network")
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
      low = if (halo) payload$palette[["non_halo"]] %||% "steelblue" else "steelblue",
      mid = "white", high = if (halo) payload$palette[["halo"]] %||% "firebrick" else "firebrick",
      na.value = "grey90", limits = if (residual || halo) c(-1, 1) else NULL,
      name = if (halo) "Retained correlation" else if (residual) "Residual correlation" else ggplot2::waiver()
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
    p, payload, x = if (halo) "Rater x criterion node" else if (residual) payload$facet else "Group or contrast",
    y = if (halo) "Rater x criterion node" else if (residual) payload$facet else "Facet level",
    fallback = "DIF/DFF heatmap"
  )
}

.mfrmr_gg_interrater <- function(x) {
  payload <- x$data
  df <- payload$table
  view <- payload$plot
  if (length(view) != 1L || is.na(view) || !view %in% c("exact", "corr", "difference") || !is.data.frame(df) ||
      !nrow(df) || !all(c("Rater1", "Rater2", "Exact", "Corr", "MeanDiff", "MAD") %in% names(df))) {
    stop("Recreate the agreement plot from its saved result to retain the selected view.", call. = FALSE)
  }
  style <- resolve_plot_preset(payload$preset %||% "standard")
  pal <- payload$palette %||% c(ok = style$accent_primary, flag = style$warn,
    expected = style$accent_secondary)
  df$Colour <- if ("Flag" %in% names(df)) {
    ifelse(is.na(df$Flag), style$neutral, ifelse(df$Flag, pal["flag"], pal["ok"]))
  } else pal["ok"]
  if (view == "difference") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$MeanDiff, y = .data$MAD)) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$Colour), na.rm = TRUE) +
      ggplot2::scale_colour_identity()
    if (!any(is.finite(df$MeanDiff) & is.finite(df$MAD)))
      p <- p + ggplot2::annotate("text", x = 0, y = 0,
        label = "No matched-score differences are available.")
    xlab <- "Mean score difference (Rater1 - Rater2)"
    ylab <- "Mean absolute difference"
  } else {
    df$Index <- seq_len(nrow(df))
    df$Value <- if (view == "exact") df$Exact else df$Corr
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Index, y = .data$Value)) +
      ggplot2::geom_col(ggplot2::aes(fill = .data$Colour), na.rm = TRUE) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_x_continuous(breaks = df$Index,
        labels = truncate_axis_label(paste0(df$Rater1, " | ", df$Rater2), width = 28L),
        limits = c(0.5, nrow(df) + 0.5)) +
      ggplot2::scale_y_continuous(limits = if (view == "exact") c(0, 1) else c(-1, 1))
    if (view == "exact" && "ExpectedExact" %in% names(df) && any(is.finite(df$ExpectedExact))) {
      if (sum(is.finite(df$ExpectedExact)) > 1L) p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$ExpectedExact,
        colour = "Model expected", group = 1), na.rm = TRUE)
      p <- p + ggplot2::geom_point(ggplot2::aes(y = .data$ExpectedExact,
          colour = "Model expected"), shape = 21, fill = "white", na.rm = TRUE) +
        ggplot2::scale_colour_manual(values = c("Model expected" = unname(pal["expected"])), name = NULL)
    }
    xlab <- "Rater pair"
    ylab <- if (view == "exact") "Exact agreement" else "Correlation"
  }
  p <- p + .mfrmr_gg_theme() + ggplot2::theme(legend.position = "bottom")
  if (view != "difference") p <- p + ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = payload$label_angle %||% 45, hjust = 1))
  p <- .mfrmr_gg_labs(.mfrmr_gg_add_references(p, payload$reference_lines), payload,
    x = xlab, y = ylab)
  attr(p, "mfrmr_plot_data") <- x
  p
}

.mfrmr_gg_agreement_heatmap <- function(x) {
  payload <- x$data
  mat <- payload$matrix
  correlation <- identical(payload$metric, "correlation")
  if (!is.matrix(mat) || !length(mat) || length(payload$metric) != 1L ||
      is.na(payload$metric) || !payload$metric %in% c("exact", "correlation"))
    stop("Recreate the agreement heatmap to retain its matrix and metric.", call. = FALSE)
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("Row", "Column", "Value")
  df$Label <- ifelse(is.finite(df$Value), sprintf("%.2f", df$Value), "")
  colors <- payload$colors %||% grDevices::hcl.colors(20L,
    if (correlation) "RdBu" else "YlGnBu", rev = TRUE)
  limits <- if (correlation) c(-1, 1) else c(0, 1)
  # Use palette luminance, including monochrome, for readable cell text.
  index <- round((df$Value - limits[1]) / diff(limits) * (length(colors) - 1)) + 1L
  rgb <- grDevices::col2rgb(colors[pmax(1L, pmin(length(colors), index))]) / 255
  df$TextColour <- ifelse(colSums(rgb * c(.299, .587, .114)) < .5, "white", "black")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Column, y = .data$Row, fill = .data$Value)) +
    ggplot2::geom_tile(colour = "white", linewidth = .25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Label, colour = .data$TextColour), size = 3) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradientn(colours = colors, limits = limits,
      na.value = "grey90", name = if (correlation) "Correlation" else "Exact agreement") +
    ggplot2::coord_equal() + .mfrmr_gg_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major = ggplot2::element_blank())
  facet <- payload$settings$rater_facet %||% "Rater"
  p <- .mfrmr_gg_labs(p, payload, x = facet, y = facet)
  attr(p, "mfrmr_plot_data") <- x
  p
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
  location_colors <- resolve_palette(payload[["palette"]],
    .plot_series_colors(c("facet_level", "step_threshold"), payload$preset %||% "standard"))
  names(location_colors) <- c("Facet level", "Step threshold")
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
          colour = .data$Group, linetype = .data$Group, group = .data$Group
        ),
        inherit.aes = FALSE, linewidth = 0.7
      )
    }
  }
  summary_df <- as.data.frame(payload$group_summary %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(summary_df) > 0L &&
      all(c("XBase", "Min", "Max", "Q1", "Q3", "Median", "PlotType") %in% names(summary_df))) {
    summary_df <- summary_df[
      as.character(summary_df$PlotType) == "Facet level",
      ,
      drop = FALSE
    ]
  }
  if (nrow(summary_df) > 0L) {
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
  step_df <- loc[as.character(loc$PlotType) == "Step threshold", , drop = FALSE]
  if (nrow(step_df) > 1L && all(c("Group", "XBase") %in% names(step_df))) {
    step_ladder <- do.call(
      rbind,
      lapply(split(step_df, as.character(step_df$Group)), function(tbl) {
        if (nrow(tbl) < 2L) return(NULL)
        data.frame(
          XBase = tbl$XBase[1L],
          Min = min(tbl$Estimate, na.rm = TRUE),
          Max = max(tbl$Estimate, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )
    if (!is.null(step_ladder) && nrow(step_ladder) > 0L) {
      p <- p + ggplot2::geom_segment(
        data = step_ladder,
        ggplot2::aes(
          x = .data$XBase, xend = .data$XBase,
          y = .data$Min, yend = .data$Max
        ),
        inherit.aes = FALSE, colour = location_colors[["Step threshold"]],
        linewidth = 0.55, alpha = 0.5
      )
    }
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
    labels$LabelX <- suppressWarnings(as.numeric(labels$LabelX %||% labels$X))
    labels$LabelY <- suppressWarnings(as.numeric(labels$LabelY %||% labels$Estimate))
    labels$PointY <- suppressWarnings(as.numeric(
      labels$DisplayEstimate %||% labels$Estimate
    ))
    labels$LabelHjust <- suppressWarnings(as.numeric(
      labels$LabelHjust %||% rep(0, nrow(labels))
    ))
    labels$LabelShort <- truncate_axis_label(
      labels$LabelText %||% labels$Label,
      width = 22L
    )
    p <- p +
      ggplot2::geom_segment(
        data = labels,
        ggplot2::aes(
          x = .data$X, y = .data$PointY,
          xend = .data$LabelX, yend = .data$LabelY
        ),
        inherit.aes = FALSE, colour = "grey45",
        linewidth = 0.25, alpha = 0.65
      ) +
      ggplot2::geom_text(
        data = labels,
        ggplot2::aes(
          x = .data$LabelX, y = .data$LabelY,
          label = .data$LabelShort, hjust = .data$LabelHjust
        ),
        inherit.aes = FALSE, size = 2.5, check_overlap = FALSE
      )
  }
  density_groups <- unique(as.character(density$Group))
  p <- p + ggplot2::scale_colour_manual(values = c(location_colors,
      .plot_series_colors(density_groups, payload$preset %||% "standard"))) +
    ggplot2::scale_linetype_manual(values = c(.plot_series_linetypes(density_groups),
      `Person mean` = "dotted", `Person median` = "dashed"))
  .mfrmr_gg_labs(p, payload, x = NULL, y = payload$axis_label %||% "Logit scale", fallback = "Wright map")
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
  colors <- resolve_palette(payload[["palette"]],
    .plot_series_colors(groups, payload$preset %||% "standard"))
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
          label = .data$Category
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
          label = .data$LabelShort
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
        colour = .data$CurveGroup, linetype = .data$CurveGroup, group = .data$CurveGroup
      ),
      linewidth = 0.85
    ) +
    ggplot2::scale_colour_manual(values = colors, limits = groups) +
    ggplot2::scale_fill_manual(values = colors, limits = groups) +
    ggplot2::scale_linetype_manual(values = .plot_series_linetypes(groups), limits = groups) +
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
    ggplot2::scale_colour_manual(values = resolve_palette(payload[["palette"]],
      .plot_series_colors(unique(as.character(tbl$Facet)), payload$preset %||% "standard"))) +
    .mfrmr_gg_theme()
  labels <- tbl[nzchar(as.character(tbl$LabelText)), , drop = FALSE]
  if (nrow(labels) > 0L) {
    p <- p + ggplot2::geom_text(
      data = labels,
      ggplot2::aes(
        x = .data$FitValue, y = .data$Measure,
        label = .data$LabelText
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
  radius <- as.numeric(payload$radius)
  if (length(radius) != nrow(tbl) || any(!is.finite(radius) | radius < 0)) {
    stop("Saved bubble radii must match the table. Recreate the payload with plot(..., type = \"bubble\", draw = FALSE).", call. = FALSE)
  }
  # ggplot point sizes use physical units, while base circles use plot units.
  # Preserve the saved radius ratios; a native maximum of .15 maps to 6 mm.
  tbl$.diameter <- 40 * radius
  tbl <- tbl[is.finite(tbl$.x) & is.finite(tbl$.y), , drop = FALSE]
  facets <- unique(as.character(tbl$Facet))
  legend <- payload$legend
  if (!is.data.frame(legend) || !all(c("role", "label", "value") %in% names(legend))) {
    stop("Saved bubble facet colours are unavailable. Recreate the payload with plot(..., type = \"bubble\", draw = FALSE).", call. = FALSE)
  }
  legend <- legend[legend$role == "facet", , drop = FALSE]
  at <- match(facets, legend$label)
  if (anyNA(at)) stop("Saved bubble facet colours do not cover the table.", call. = FALSE)
  colours <- stats::setNames(as.character(legend$value[at]), facets)
  tbl$Facet <- factor(as.character(tbl$Facet), levels = facets)
  p <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = .data$.x, y = .data$.y, colour = .data$Facet)
  ) + ggplot2::geom_point(ggplot2::aes(size = .data$.diameter), alpha = 0.7) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_colour_manual(values = colours, breaks = facets) + .mfrmr_gg_theme()
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
  categories <- unique(as.character(prob$Category))
  prob$Category <- factor(as.character(prob$Category), levels = categories)
  colors <- resolve_palette(payload[["palette"]],
    .plot_series_colors(categories, payload$preset %||% "standard"))
  if (!"Slope" %in% names(prob)) prob$Slope <- 1
  prob$Slope <- suppressWarnings(as.numeric(prob$Slope))
  prob$Trace <- paste(prob$CurveGroup, prob$Category, sep = " | ")
  mapping <- switch(
    slope_aes,
    linewidth = ggplot2::aes(
      x = .data$Theta, y = .data$Probability, colour = .data$Category,
      group = .data$Trace, linetype = .data$Category, linewidth = .data$Slope
    ),
    alpha = ggplot2::aes(
      x = .data$Theta, y = .data$Probability, colour = .data$Category,
      group = .data$Trace, linetype = .data$Category, alpha = .data$Slope
    ),
    colour = ggplot2::aes(
      x = .data$Theta, y = .data$Probability, colour = .data$Slope,
      linetype = .data$Category, group = .data$Trace
    ),
    none = ggplot2::aes(
      x = .data$Theta, y = .data$Probability, colour = .data$Category,
      group = .data$Trace, linetype = .data$Category
    )
  )
  p <- ggplot2::ggplot(prob, mapping)
  p <- if (identical(slope_aes, "linewidth")) {
    p +
      ggplot2::geom_line() +
      ggplot2::scale_linewidth_continuous(range = c(0.45, 1.35))
  } else {
    p + ggplot2::geom_line(linewidth = 0.85)
  }
  p <- p + ggplot2::scale_linetype_manual(
    values = .plot_series_linetypes(categories), limits = categories) +
    ggplot2::scale_fill_manual(values = colors, limits = categories)
  p <- p + if (identical(slope_aes, "colour")) {
    if (identical(payload$preset, "monochrome")) {
      ggplot2::scale_colour_gradient(low = "gray15", high = "gray55")
    } else ggplot2::scale_colour_viridis_c(end = 0.8)
  } else ggplot2::scale_colour_manual(values = colors, limits = categories)
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
    overlay$Category <- factor(as.character(overlay$Category), levels = categories)
    p <- p + ggplot2::geom_point(
      data = overlay,
      ggplot2::aes(
          x = .data$Theta, y = .data$Proportion, fill = .data$Category,
          shape = .data$Category
      ),
      inherit.aes = FALSE, alpha = 0.75
    ) + ggplot2::scale_shape_manual(values = stats::setNames(
      rep(c(21, 22, 24, 23, 25), length.out = length(categories)), categories))
  }
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  .mfrmr_gg_labs(
    p, payload, x = "Theta / logit", y = "Category probability",
    fallback = "Category characteristic curves"
  )
}

.mfrmr_gg_calibration_score <- function(payload, type) {
  df <- as.data.frame(payload$data %||% data.frame(), stringsAsFactors = FALSE)
  required <- switch(
    type,
    interval = c("Person", "Estimate", "Lower", "Upper", "ReviewFlag"),
    precision = c("Person", "Observations", "SD", "ReviewFlag"),
    edge_mass = c(
      "Person", "Estimate", "QuadratureEdgeMass", "ReviewFlag"
    )
  )
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop(
      "Portable-score conversion is missing field(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  df$ReviewStatus <- factor(
    ifelse(df$ReviewFlag, "Review", "Scored"),
    levels = c("Scored", "Review")
  )
  encoding <- as.data.frame(
    payload$display_encoding %||% data.frame(), stringsAsFactors = FALSE
  )
  required_encoding <- c("Status", "Colour", "Shape")
  if (nrow(encoding) >= 2L &&
      all(required_encoding %in% names(encoding)) &&
      all(c("Scored", "Review") %in% encoding$Status)) {
    encoding <- encoding[match(c("Scored", "Review"), encoding$Status), ]
    colors <- stats::setNames(as.character(encoding$Colour), encoding$Status)
    shapes <- stats::setNames(as.numeric(encoding$Shape), encoding$Status)
  } else {
    colors <- c(Scored = "#1F78B4", Review = "#B65E16")
    shapes <- c(Scored = 16, Review = 17)
  }
  if (identical(type, "interval")) {
    display_order <- suppressWarnings(as.numeric(
      df$DisplayOrder %||% seq_len(nrow(df))
    ))
    df$PersonFactor <- factor(
      as.character(df$Person),
      levels = as.character(df$Person[order(display_order)])
    )
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        y = .data$PersonFactor, colour = .data$ReviewStatus,
        shape = .data$ReviewStatus
      )
    ) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .data$Lower, xend = .data$Upper,
          yend = .data$PersonFactor
        ),
        linewidth = 0.65
      ) +
      ggplot2::geom_point(ggplot2::aes(x = .data$Estimate), size = 2.2) +
      ggplot2::scale_colour_manual(values = colors) +
      ggplot2::scale_shape_manual(values = shapes) +
      .mfrmr_gg_theme() +
      ggplot2::theme(
        legend.position = "bottom", legend.title = ggplot2::element_blank()
      )
    p <- .mfrmr_gg_add_references(p, payload$reference_lines)
    return(.mfrmr_gg_labs(
      p, payload, x = "Posterior EAP (logits)", y = NULL,
      fallback = "Portable calibration score intervals"
    ))
  }
  if (identical(type, "precision")) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$Observations, y = .data$SD,
        colour = .data$ReviewStatus, shape = .data$ReviewStatus
      )
    ) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::scale_colour_manual(values = colors) +
      ggplot2::scale_shape_manual(values = shapes) +
      .mfrmr_gg_theme() +
      ggplot2::theme(
        legend.position = "bottom", legend.title = ggplot2::element_blank()
      )
    labels <- df[df$ReviewFlag & nzchar(as.character(df$Person)), , drop = FALSE]
    if (isTRUE(payload$settings$label_review %||% TRUE) && nrow(labels) > 0L) {
      if (!"LabelX" %in% names(labels)) labels$LabelX <- labels$Observations
      if (!"LabelY" %in% names(labels)) labels$LabelY <- labels$SD
      p <- p + ggplot2::geom_segment(
        data = labels,
        ggplot2::aes(
          x = .data$Observations, y = .data$SD,
          xend = .data$LabelX, yend = .data$LabelY
        ),
        inherit.aes = FALSE, colour = unname(colors["Review"]),
        linewidth = 0.3
      ) + ggplot2::geom_text(
        data = labels,
        ggplot2::aes(
          x = .data$LabelX, y = .data$LabelY, label = .data$Person
        ),
        inherit.aes = FALSE, colour = unname(colors["Review"]),
        vjust = -0.5, show.legend = FALSE, check_overlap = FALSE
      )
    }
    return(.mfrmr_gg_labs(
      p, payload, x = "Valid response rows", y = "Posterior SD (logits)",
      fallback = "Conditional score precision"
    ))
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$Estimate, y = .data$QuadratureEdgeMass,
      colour = .data$ReviewStatus, shape = .data$ReviewStatus
    )
  ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::scale_shape_manual(values = shapes) +
    .mfrmr_gg_theme() +
    ggplot2::theme(
      legend.position = "bottom", legend.title = ggplot2::element_blank()
    )
  labels <- df[df$ReviewFlag & nzchar(as.character(df$Person)), , drop = FALSE]
  if (isTRUE(payload$settings$label_review %||% TRUE) && nrow(labels) > 0L) {
    if (!"LabelPosition" %in% names(labels)) labels$LabelPosition <- 3L
    if (!"LabelX" %in% names(labels)) labels$LabelX <- labels$Estimate
    labels$LabelHjust <- ifelse(
      labels$LabelPosition == 2L, 1,
      ifelse(labels$LabelPosition == 4L, 0, 0.5)
    )
    p <- p + ggplot2::geom_text(
      data = labels,
      ggplot2::aes(
        x = .data$LabelX, label = .data$Person,
        hjust = .data$LabelHjust
      ),
      vjust = -0.7, show.legend = FALSE, check_overlap = FALSE
    )
  }
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  .mfrmr_gg_labs(
    p, payload, x = "Posterior EAP (logits)",
    y = "Posterior mass at outer quadrature nodes",
    fallback = "Quadrature edge-mass review"
  )
}

.mfrmr_gg_bias <- function(x) {
  payload <- x$data
  ranked <- identical(payload$plot, "ranked")
  df <- if (ranked) payload$ranked_table else payload$scatter_data
  df <- df[is.finite(df$BiasSize), , drop = FALSE]
  if (ranked) df <- utils::head(df, payload$plot_settings$TopN %||% 40L)
  pal <- resolve_palette(payload$palette, c(normal = "#2B6CB0", flag = "#B2182B"))
  df$Colour <- ifelse(df$Flag %in% TRUE, pal[["flag"]], pal[["normal"]])
  df$Position <- seq_len(nrow(df))
  p <- if (ranked) ggplot2::ggplot(df, ggplot2::aes(x = .data$BiasSize, y = .data$Position)) else
    ggplot2::ggplot(df, ggplot2::aes(x = .data$ObsExpAverage, y = .data$BiasSize))
  if (isTRUE(payload$plot_settings$CIRequested) && all(c("CI_Lower", "CI_Upper", "CI_Complete") %in% names(df))) {
    intervals <- df[is.finite(df$CI_Lower) & is.finite(df$CI_Upper), , drop = FALSE]
    intervals$LineType <- ifelse(intervals$CI_Complete %in% TRUE, "solid", "dashed")
    mapping <- if (ranked) ggplot2::aes(x = .data$CI_Lower, xend = .data$CI_Upper,
      y = .data$Position, yend = .data$Position, colour = .data$Colour, linetype = .data$LineType) else
      ggplot2::aes(x = .data$ObsExpAverage, xend = .data$ObsExpAverage,
        y = .data$CI_Lower, yend = .data$CI_Upper, colour = .data$Colour, linetype = .data$LineType)
    p <- p + ggplot2::geom_segment(data = intervals, mapping = mapping, inherit.aes = FALSE) +
      ggplot2::scale_linetype_identity()
  }
  p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$Colour), na.rm = TRUE) +
    ggplot2::scale_colour_identity() + .mfrmr_gg_theme()
  if (ranked) p <- p + ggplot2::scale_y_continuous(breaks = df$Position, labels = df$Pair)
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  .mfrmr_gg_labs(p, payload, x = if (ranked) "Bias size (logits)" else "Observed - expected average",
    y = if (ranked) NULL else "Bias size (logits)")
}

.mfrmr_gg_halo <- function(x) {
  payload <- x$data
  if (identical(payload$plot, "matrix")) return(.mfrmr_gg_matrix(payload$matrix, payload))
  df <- payload$table
  pal <- resolve_palette(payload$palette, c(halo = "#B65E16", non_halo = "#2B6CB0",
    warning = "#B2182B", review = "#B65E16", ok = "#2B6CB0", neutral = "grey60"))
  if (identical(payload$plot, "edge_distribution")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$EdgeType, y = .data$AbsEstimate, fill = .data$EdgeType)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1), alpha = 0.35) +
      ggplot2::scale_fill_manual(values = pal[c("halo", "non_halo")], guide = "none") +
      ggplot2::scale_x_discrete(limits = c("halo", "non_halo"), labels = c("Same rater", "Other raters"))
    xlab <- NULL; ylab <- "Absolute correlation"
  } else {
    df$Position <- seq_len(nrow(df))
    df$Value <- df[[payload$metric]]
    df$Colour <- unname(pal[ifelse(df$ReviewStatus %in% c("warning", "review", "ok"), df$ReviewStatus, "neutral")])
    available <- is.finite(df$Value)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Position, y = .data$Value)) +
      ggplot2::geom_col(data = df[available, , drop = FALSE], ggplot2::aes(fill = .data$Colour)) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_x_continuous(breaks = df$Position, labels = df$Rater)
    if (any(!available)) p <- p + ggplot2::geom_text(data = df[!available, , drop = FALSE],
      ggplot2::aes(x = .data$Position), y = 0, label = "Unavailable", vjust = -0.5,
      inherit.aes = FALSE, size = 3)
    xlab <- payload$settings$rater_facet %||% "Rater"; ylab <- payload$metric
  }
  p <- p + .mfrmr_gg_theme()
  .mfrmr_gg_labs(p, payload, x = xlab, y = ylab)
}

.mfrmr_gg_categories <- function(x) {
  payload <- x$data
  counts <- identical(payload$plot, "counts")
  pal <- resolve_palette(payload$palette, c(counts = "#9ecae1", expected = "#08519c",
    step_line = "#238b45"))
  if (counts) {
    df <- as.data.frame(payload$category_table)
    if (!nrow(df) || !all(c("Category", "Count") %in% names(df)))
      stop("Saved category counts are required for this view.", call. = FALSE)
    df$Position <- seq_len(nrow(df))
    labels <- as.character(df$Category)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Position, y = .data$Count)) +
      ggplot2::geom_col(data = df[is.finite(df$Count), , drop = FALSE], fill = pal[["counts"]])
    if ("ExpectedCount" %in% names(df)) {
      df$ExpectedCount[!is.finite(df$ExpectedCount)] <- NA_real_
      p$data <- df
      p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$ExpectedCount), colour = pal[["expected"]], na.rm = TRUE) +
        ggplot2::geom_point(ggplot2::aes(y = .data$ExpectedCount), shape = 21, fill = "white",
          colour = pal[["expected"]], na.rm = TRUE)
    }
    missing <- !is.finite(df$Count)
    xlab <- "Category (model score)"; ylab <- "Count"
  } else {
    df <- calc_step_order(payload$threshold_table)
    if (!nrow(df)) stop("Saved step estimates are required for this view.", call. = FALSE)
    df$Estimate[!is.finite(df$Estimate)] <- NA_real_
    df$Position <- seq_len(nrow(df))
    labels <- if (all(df$StepFacet == "Common")) as.character(df$Step) else
      paste(df$StepFacet, df$Step, sep = ": ")
    ends <- which(!is.na(df$Ordered))
    segments <- data.frame(X = ends - 1L, Xend = ends,
      Y = df$Estimate[ends - 1L], Yend = df$Estimate[ends])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Position, y = .data$Estimate)) +
      ggplot2::geom_segment(data = segments,
        ggplot2::aes(x = .data$X, xend = .data$Xend, y = .data$Y, yend = .data$Yend),
        inherit.aes = FALSE, colour = pal[["step_line"]]) +
      ggplot2::geom_point(colour = pal[["step_line"]], na.rm = TRUE)
    missing <- !is.finite(df$Estimate)
    xlab <- "Step within owner level"; ylab <- "Stored step estimate"
  }
  if (any(missing)) p <- p + ggplot2::geom_text(data = df[missing, , drop = FALSE],
    ggplot2::aes(x = .data$Position), y = 0, label = "Unavailable",
    inherit.aes = FALSE, vjust = -0.5, size = 3)
  p <- p + ggplot2::scale_x_continuous(breaks = df$Position, labels = labels) +
    .mfrmr_gg_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = payload$label_angle %||% 45, hjust = 1))
  .mfrmr_gg_labs(p, payload, x = xlab, y = ylab)
}

.mfrmr_gg_linking <- function(x) {
  payload <- x$data
  style <- resolve_plot_preset(payload$preset %||% "standard")
  view <- payload$plot %||% x$name
  if (identical(view, "heatmap")) {
    df <- as.data.frame(as.table(payload$matrix), stringsAsFactors = FALSE)
    names(df) <- c("Element", "Wave", "Drift")
    df$Element <- factor(df$Element, levels = rev(rownames(payload$matrix)))
    extent <- if (any(is.finite(df$Drift))) max(abs(df$Drift), na.rm = TRUE) else 1
    if (extent == 0) extent <- 1
    colours <- payload$legend$value
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Wave, y = .data$Element,
      fill = .data$Drift)) + ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = colours[1L], mid = "white",
        high = colours[2L], limits = c(-extent, extent), na.value = "grey90",
        name = "Drift (logits)") + .mfrmr_gg_theme()
    p <- .mfrmr_gg_labs(p, payload, x = "Wave", y = "Facet level")
  } else if (identical(view, "chain")) {
    df <- as.data.frame(payload$table)
    df$Position <- seq_len(nrow(df))
    df$NextOffset <- c(tail(df$Cumulative_Offset, -1L), NA_real_)
    edges <- is.finite(df$Cumulative_Offset) & is.finite(df$NextOffset)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Position, y = .data$Cumulative_Offset)) +
      ggplot2::geom_segment(data = df[edges, , drop = FALSE],
        ggplot2::aes(xend = .data$Position + 1, yend = .data$NextOffset),
        colour = style$accent_primary) +
      ggplot2::geom_point(na.rm = TRUE, colour = style$accent_primary) +
      ggplot2::scale_x_continuous(breaks = df$Position, labels = df$Wave) +
      .mfrmr_gg_theme()
    p <- .mfrmr_gg_add_references(p, payload$reference_lines)
    p <- .mfrmr_gg_labs(p, payload, x = "Wave", y = "Cumulative offset (logits)")
  } else {
    forest <- identical(x$name, "anchor_drift_forest")
    common <- identical(x$name, "equating_chain_common_anchors")
    df <- as.data.frame(if (forest || common) payload$data else payload$table)
    if (!forest && !common) df <- df[order(-abs(df$Drift), na.last = TRUE), , drop = FALSE]
    df$Position <- seq_len(nrow(df))
    df$Value <- if (forest) df$Estimate else if (common) df$N_Common else df$Drift
    labels <- if (forest) df$Label else if (common) df$Pair else
      paste(paste(df$Facet, df$Level, sep = ":"), df$Wave, sep = " @ ")
    df$Colour <- if (forest || common) style$accent_primary else
      ifelse(df$Flag %in% TRUE, style$warn, style$accent_primary)
    shown <- is.finite(df$Value)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Value, y = .data$Position))
    if (forest) {
      intervals <- is.finite(df$CI_Lower) & is.finite(df$CI_Upper)
      p <- p + ggplot2::geom_segment(data = df[intervals, , drop = FALSE],
        ggplot2::aes(x = .data$CI_Lower, xend = .data$CI_Upper, yend = .data$Position),
        colour = style$accent_primary)
    }
    if (common) p <- p + ggplot2::geom_segment(data = df[shown, , drop = FALSE],
      ggplot2::aes(x = 0, xend = .data$Value, yend = .data$Position),
      linewidth = 4, colour = style$accent_primary)
    p <- p + ggplot2::geom_point(data = df[shown, , drop = FALSE],
      ggplot2::aes(colour = .data$Colour), size = 2) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_y_reverse(breaks = df$Position, labels = labels) +
      .mfrmr_gg_theme()
    p <- .mfrmr_gg_add_references(p, payload$reference_lines)
    p <- .mfrmr_gg_labs(p, payload, x = if (forest) "Estimate (logits)" else if (common)
      "N common" else "Drift (logits)", y = NULL)
  }
  attr(p, "mfrmr_plot_data") <- x
  p
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

.mfrmr_gg_residual_pca <- function(x) {
  .require_mfrmr_ggplot2()
  payload <- x$data
  df <- payload$data
  view <- payload$plot
  required <- switch(view,
    scree = c("Component", "Eigenvalue"),
    parallel_scree = c("Component", "Eigenvalue", "ParallelCutoff", "ParallelMean"),
    parallel_excess = c("Component", "ExcessOverParallelCutoff", "ExceedsParallelCutoff"),
    loadings = c("Variable", "Loading"),
    stop("Unsupported residual-PCA view. Recreate the plot from its result.", call. = FALSE))
  if (!is.data.frame(df) || !nrow(df) || !all(required %in% names(df))) {
    stop("Residual-PCA conversion requires the saved values for the selected view.", call. = FALSE)
  }
  legend <- payload$legend
  if (!all(c("label", "value") %in% names(legend)) || nrow(legend) < 2L) {
    stop("Recreate the residual-PCA plot to retain its display encoding.", call. = FALSE)
  }
  if (view == "parallel_scree") {
    labels <- legend$label[1:3]
    if (length(labels) != 3L || anyNA(labels)) stop("The reference legend is incomplete.", call. = FALSE)
    long <- do.call(rbind, lapply(seq_along(labels), function(i) {
      data.frame(Component = df$Component,
        Value = df[[c("Eigenvalue", "ParallelCutoff", "ParallelMean")[i]]],
        Series = factor(labels[i], levels = labels))
    }))
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$Component, y = .data$Value,
      colour = .data$Series, linetype = .data$Series, shape = .data$Series)) +
      ggplot2::geom_line() + ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::scale_colour_manual(values = stats::setNames(legend$value[1:3], labels)) +
      ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "dashed", "dotted"), labels)) +
      ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17, NA), labels)) +
      ggplot2::labs(colour = NULL, linetype = NULL, shape = NULL)
    guide <- ggplot2::guide_legend(ncol = 2, byrow = TRUE)
    p <- p + ggplot2::guides(colour = guide, linetype = guide, shape = guide)
    xlab <- "Component"; ylab <- "Eigenvalue"
  } else if (view == "scree") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Component, y = .data$Eigenvalue)) +
      ggplot2::geom_line(colour = legend$value[1]) +
      ggplot2::geom_point(colour = legend$value[1])
    xlab <- "Component"; ylab <- "Eigenvalue"
  } else {
    positive <- if (view == "loadings") df$Loading >= 0 else df$ExceedsParallelCutoff
    df$Direction <- factor(ifelse(positive, legend$label[1], legend$label[2]),
      levels = legend$label[1:2])
    if (view == "loadings") {
      df$Variable <- factor(df$Variable, levels = unique(df$Variable))
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Loading, y = .data$Variable,
        fill = .data$Direction)) + ggplot2::geom_col(orientation = "y")
      xlab <- "Loading"; ylab <- NULL
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Component,
        y = .data$ExcessOverParallelCutoff, fill = .data$Direction)) + ggplot2::geom_col()
      xlab <- "Component"; ylab <- "Eigenvalue minus cutoff"
    }
    p <- p + ggplot2::scale_fill_manual(values = stats::setNames(
      legend$value[1:2], legend$label[1:2]), drop = FALSE) + ggplot2::labs(fill = NULL)
  }
  if (view != "loadings") p <- p + ggplot2::scale_x_continuous(breaks = df$Component)
  p <- p + .mfrmr_gg_theme() + ggplot2::theme(legend.position = "bottom")
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  p <- .mfrmr_gg_labs(p, payload, x = xlab, y = ylab)
  p <- p + ggplot2::labs(alt = "Exploratory residual structure, not a dimensionality decision.")
  attr(p, "mfrmr_plot_data") <- x
  p
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

.mfrmr_gg_multivariate_gstudy <- function(payload) {
  df <- payload$table
  style <- resolve_plot_preset(payload$preset %||% "standard")
  if (identical(payload$plot, "covariance")) {
    df$Row <- factor(df$Row, levels = rev(payload$score))
    df$Column <- factor(df$Column, levels = payload$score)
    df$Panel <- factor(df$Panel, levels = payload$panels$Panel, labels = payload$panels$Label)
    colors <- grDevices::colorRampPalette(payload$heat_colors)(101)
    index <- round((df$Estimate / payload$fill_limits[2L] + 1) * 50) + 1L
    df$Ink <- ifelse(colMeans(grDevices::col2rgb(colors[index])) < 140, "white", "black")
    df$CellLabel <- trimws(formatC(df$Estimate, format = "g", digits = 3))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Column, y = .data$Row)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$Estimate), colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = .data$CellLabel, colour = .data$Ink), size = 3) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_fill_gradientn(colours = colors, limits = payload$fill_limits,
        breaks = c(payload$fill_limits[1L], 0, payload$fill_limits[2L]),
        labels = function(x) trimws(formatC(x, format = "g", digits = 3)), name = "Variance / covariance") +
      ggplot2::scale_x_discrete(labels = function(x) vapply(x, function(s)
        paste(strwrap(s, 14), collapse = "\n"), character(1)), expand = c(0, 0)) +
      ggplot2::scale_y_discrete(labels = function(x) vapply(x, function(s)
        paste(strwrap(s, 14), collapse = "\n"), character(1)), expand = c(0, 0)) +
      ggplot2::facet_wrap(~Panel, ncol = payload$ncol, labeller = ggplot2::label_wrap_gen(28)) +
      ggplot2::coord_fixed() + .mfrmr_gg_theme() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")
    return(.mfrmr_gg_labs(p, payload, x = NULL, y = NULL))
  }
  sources <- unique(df$Source)
  df$Position <- length(sources) + 1L - match(df$Source, sources)
  faceted <- !is.null(payload$panels) && nrow(payload$panels) > 1L
  if (faceted) df$Panel <- factor(df$Panel, levels = payload$panels$Panel, labels = payload$panels$Label)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Estimate, y = .data$Position)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60")
  if (!is.null(payload$intervals)) {
    p <- p + ggplot2::geom_errorbar(data = df[df$IntervalAvailable, , drop = FALSE],
      ggplot2::aes(xmin = .data$Lower, xmax = .data$Upper),
      orientation = "y", width = .18, colour = style$foreground) +
      ggplot2::geom_point(data = df[!df$IntervalAvailable, , drop = FALSE],
        shape = 0, size = 3.75, colour = style$foreground)
  } else p <- p + ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = .data$Estimate, yend = .data$Position), colour = style$foreground)
  p <- p +
    ggplot2::geom_point(ggplot2::aes(shape = .data$Status), size = 2.5, colour = style$foreground) +
    ggplot2::scale_shape_manual(values = c("Negative estimate" = 1, "Nonnegative estimate" = 16), guide = "none") +
    ggplot2::scale_y_continuous(breaks = rev(seq_along(sources)), labels = payload$labels,
      limits = c(.5, length(sources) + .5), expand = ggplot2::expansion(mult = 0)) +
    .mfrmr_gg_theme()
  if (faceted) {
    p <- p + ggplot2::facet_wrap(~Panel, ncol = payload$ncol, labeller = ggplot2::label_wrap_gen(28))
  }
  .mfrmr_gg_labs(p, payload, x = payload$x_label, y = NULL)
}

.mfrmr_gg_multivariate_d_study <- function(payload) {
  if (identical(payload$plot, "error_components")) {
    df <- payload$series
    df$SourceLabel <- factor(df$SourceLabel, levels = rev(payload$source_labels))
    df$Plan <- factor(df$Plan, levels = payload$legend$label)
    df$Panel <- factor(payload$panels$Label[match(df$Panel, payload$panels$Panel)],
      levels = payload$panels$Label)
    dodge <- ggplot2::position_dodge(width = .6)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$SourceLabel, y = .data$Contribution,
      colour = .data$Plan, shape = .data$Plan, group = .data$Plan)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey75") +
      ggplot2::geom_point(position = dodge, size = 2.2) +
      ggplot2::facet_wrap(~Panel, ncol = payload$ncol, scales = "free_y") +
      ggplot2::scale_y_continuous(limits = range(c(0, df$Contribution))) +
      ggplot2::coord_flip() +
      ggplot2::scale_colour_manual(values = stats::setNames(payload$legend$value, payload$legend$label)) +
      ggplot2::scale_shape_manual(values = stats::setNames(payload$point_shapes, payload$legend$label)) +
      ggplot2::labs(colour = NULL, shape = NULL) +
      .mfrmr_gg_theme() + ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1),
        shape = ggplot2::guide_legend(ncol = 1))
    return(.mfrmr_gg_labs(p, payload, y = payload$x_label))
  }
  df <- payload$series
  faceted <- !is.null(payload$panel_by)
  df$Panel <- if (faceted) factor(payload$panels$Label[match(df$Panel, payload$panels$Panel)],
    levels = payload$panels$Label) else
    factor(payload$panel_labels[match(df$Metric, payload$metric)], levels = payload$panel_labels)
  groups <- payload$legend$label
  df$Group <- factor(df$Group, levels = groups)
  df$Value[df$Status != "Available" | !is.finite(df$Value)] <- NA_real_
  coefficients <- identical(payload$plot, "coefficients")
  y_lim <- .mfrm_mvds_y_limits(df, coefficients, payload$reference_lines$value)
  ticks <- sort(unique(df$X))
  if (length(ticks) > 15L) ticks <- pretty(range(ticks))
  ticks <- ticks[ticks >= min(df$X) & ticks <= max(df$X) & ticks == floor(ticks)]
  # Keep internal NA rows so lines break at unavailable scenarios. A series
  # with fewer than two available points has no line to draw.
  line_rows <- stats::ave(as.integer(is.finite(df$Value)), df$Panel, df$Group, FUN = sum) >= 2L
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$X, y = .data$Value,
    colour = .data$Group, linetype = .data$Group, shape = .data$Group, group = .data$Group))
  if (!is.null(payload$intervals)) {
    bands <- df[df$IntervalAvailable, , drop = FALSE]
    single <- stats::ave(bands$X, bands$IntervalRun, FUN = function(x) length(unique(x))) == 1L
    p <- p + ggplot2::geom_ribbon(data = bands[!single, , drop = FALSE],
      ggplot2::aes(x = .data$X, ymin = .data$Lower, ymax = .data$Upper,
        fill = .data$Group, colour = .data$Group, linetype = .data$Group,
        group = .data$IntervalRun),
      inherit.aes = FALSE, alpha = .14, linewidth = .35, outline.type = "both", show.legend = FALSE) +
      ggplot2::geom_linerange(data = bands[single, , drop = FALSE],
        ggplot2::aes(x = .data$X, ymin = .data$Lower, ymax = .data$Upper,
          colour = .data$Group, linetype = .data$Group),
        inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = stats::setNames(payload$legend$value, groups), guide = "none")
  }
  p <- p + ggplot2::geom_line(data = df[line_rows, , drop = FALSE], linewidth = 0.7, na.rm = TRUE) +
    ggplot2::geom_point(data = if (is.null(payload$intervals)) df else
      df[df$IntervalAvailable, , drop = FALSE], size = 2.2, na.rm = TRUE) +
    ggplot2::scale_x_continuous(limits = range(df$X), breaks = ticks) +
    ggplot2::scale_y_continuous(limits = y_lim,
      breaks = if (coefficients && identical(y_lim, c(0, 1))) seq(0, 1, 0.2) else ggplot2::waiver()) +
    ggplot2::scale_colour_manual(values = stats::setNames(payload$legend$value, groups), drop = FALSE, name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(payload$line_types, groups), drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = stats::setNames(payload$point_shapes, groups), drop = FALSE, name = NULL) +
    .mfrmr_gg_theme() +
    ggplot2::theme(legend.position = if (payload$show_legend %||% (length(payload$group_var) > 0L)) "bottom" else "none")
  if (!is.null(payload$intervals)) for (j in seq_along(groups)) {
    p <- p + ggplot2::geom_point(data = df[!df$IntervalAvailable & df$Group == groups[j], , drop = FALSE],
      shape = c(1, 2, 0, 5, 1, 2)[(j-1L) %% 6L+1L], size = 2.2, na.rm = TRUE, show.legend = FALSE)
  }
  if (isTRUE(payload$show_legend) && length(payload$score) > 1L) {
    guide <- ggplot2::guide_legend(ncol = 2, byrow = TRUE)
    p <- p + ggplot2::guides(colour = guide, linetype = guide, shape = guide)
  }
  if (faceted || length(payload$metric) > 1L) p <- p +
    ggplot2::facet_wrap(~Panel, ncol = payload$ncol %||% 1L, drop = FALSE)
  p <- .mfrmr_gg_add_references(p, payload$reference_lines)
  for (panel in levels(df$Panel)) {
    s <- df[df$Panel == panel, , drop = FALSE]
    if (all(is.na(s$Value)) && !any(s$IntervalAvailable %||% FALSE)) {
      note <- data.frame(Panel = s$Panel[1L], X = mean(range(df$X)), Y = mean(y_lim),
        Label = .mfrm_mvds_unavailable_message(s$Status, compact = faceted))
      p <- p + ggplot2::geom_text(data = note,
        ggplot2::aes(x = .data$X, y = .data$Y, label = .data$Label),
        inherit.aes = FALSE, size = 3.2)
    }
  }
  # Earlier payloads lack a caption; component notes are already merged by labs.
  if (is.null(payload$caption)) {
    payload$caption <- "Points are requested scenarios; lines are guides."
    if (nrow(payload$unavailable)) payload$caption <- paste(payload$caption,
      sprintf("Unavailable estimates in %d scenario(s); inspect the D-study table for reasons.",
        length(unique(payload$unavailable$Scenario))))
  }
  .mfrmr_gg_labs(p, payload, x = payload$x_label %||% paste("Number of", tolower(payload$x_var)),
    y = payload$y_label %||% if (coefficients) "Dependability (higher is better)" else "SEM in score units (lower is better)")
}

# Earlier saved plot payloads are upgraded only for rendering, without refitting.
.mfrmr_gg_extended_payload <- function(payload, name) {
  if (!is.null(payload$display$style)) return(.mfrmr_gg_extended_estimates(payload))
  scores <- name %in% c("testlet_scores", "random_rater_scores", "person_scores")
  raters <- name %in% c("random_rater_severity", "random_rater_interval_bootstrap")
  bootstrap <- name == "random_rater_interval_bootstrap"
  labels <- payload$labels %||% if (scores) payload$table$Person else if (raters) payload$table$Rater else payload$table$Level
  title <- payload$title %||% if (scores) "Conditional Person scores" else if (bootstrap) "Observed-rater prediction intervals" else if (raters) "Shared-rater severity" else paste("Fixed", payload$settings$facet, "effects")
  xlab <- payload$xlab %||% if (scores) "Person ability (logits)" else if (raters) "Rater severity (logits)" else "Facet severity (logits)"
  note <- payload$caption %||% if (scores) sprintf(
    "%.0f%% conditional intervals | calibration fixed | %sopen: prior only", 100 * payload$settings$level,
    if (name == "random_rater_scores") "rater Laplace | " else ""
  ) else if (bootstrap) sprintf("Nominal %.0f%% pointwise | %d planned refits | arrows: unbounded%s",
    100 * payload$settings$plotted_level, payload$settings$nsim,
    if (isTRUE(payload$settings$comparison)) "\nSolid: bootstrap; dashed: ordinary normal." else "") else if (raters) {
    if (isTRUE(payload$checks$EstimatedVarianceBoundary) || isTRUE(payload$checks$EstimatedPersonVarianceBoundary))
      "Estimated variance is zero; regular intervals unavailable" else if (!isTRUE(payload$checks$NumericalReady) || !isTRUE(payload$checks$InformationPositive))
      "Numerical review required; intervals unavailable" else
      "Saved 95% normal approximation, with calibration uncertainty to first order; nominal coverage is not established"
  } else "Approximate 95% intervals | missing whiskers: unavailable"
  updated <- mfrm_extended_estimate_plot(payload$table, labels, name, title, xlab, note,
    FALSE, payload$settings, show_title = payload$display$show_title %||% FALSE,
    show_notes = payload$display$show_notes %||% FALSE)
  .mfrmr_gg_extended_estimates(updated$data)
}

#' Convert draw-free mfrmr plot data to ggplot2
#'
#' `as_ggplot()` is an optional renderer for an `mfrm_plot_data` object or an
#' object whose `plot()` method supports `draw = FALSE`. Base graphics remain
#' the default; the returned `ggplot` can be restyled or composed downstream.
#'
#' Bubble conversion uses the saved radii, facet colours, facet order and
#' reference lines. It preserves relative circle radii, rather than replacing
#' them with equal-sized points. Base circles use plot units; ggplot point
#' sizes use physical units, so absolute sizes need not match across devices. Payloads lacking matching saved radii
#' or facet colours must be recreated with [plot(x, type = "bubble")][plot.mfrm_fit].
#'
#' Dedicated conversions are provided for Wright maps, theta-to-expected-score
#' pathways, fit-statistic-to-measure pathways, category characteristic curves,
#' bubble charts, DIF/DFF summaries and heatmaps, portable-calibration
#' score review plots, testlet calibration, shared-rater severity and bootstrap
#' intervals, conditional Person scores, and
#' multivariate D-study comparisons. Matched ordinary/extended facet comparisons
#' preserve centered effects, paired/difference views, excluded rows, display
#' controls and alternative text; see [plot.mfrm_extended_comparison()].
#' Portable-calibration score views preserve their selected axes, stored
#' endpoints, review symbols and unplotted dispositions with either the default
#' or `component = "data"`. Use [plot_data()] for other components; conversion
#' does not rescore Persons or change their interval level.
#' Agreement conversions retain the selected exact/correlation/difference
#' metric, saved pair ordering and review thresholds, model-expected overlay,
#' and missing comparisons. Agreement heatmaps use the saved facet labels and
#' palette with a fixed metric scale. Use `plot_data()` for other tables in
#' an agreement payload; automatic conversion does not guess a different
#' statistic from the first numeric column.
#' Linking conversions preserve the selected drift, wave estimates and stored
#' interval endpoints, cumulative offsets, or common-element counts. Missing
#' offsets do not become zero or connect a line across an unavailable link.
#' Drift heatmaps use their saved residual matrix. Use [plot_data()] for other
#' tables; conversion does not recompute offsets or their uncertainty.
#' Threshold-sensitivity tile and curve payloads also preserve their selected
#' view, display controls and text alternatives; see
#' [plot.mfrm_screening_sensitivity()].
#' Response residual displays preserve their posterior or conditional plug-in
#' probability definition and lack of reference cutoffs. Paired views retain
#' available Infit when Outfit is undefined; see [plot.mfrm_response_diagnostics()].
#' Testlet conversions retain unavailable
#' rows, prior-only symbols and the conditional-interval note; they do not
#' estimate diagnostics or add calibration uncertainty. Extended-model interval,
#' precision and empirical-distribution views preserve display settings, data
#' exclusions and text alternatives; see [plot.mfrm_testlet_scores()]. Bootstrap
#' conversion retains infinite endpoints as arrows and ordinary intervals as
#' a dashed comparison. D-study conversions
#' preserve G/Phi or SEM panels, fixed-count groups, score units, and unavailable
#' estimates. They retain supplied D-study coefficient ribbons and interval
#' gaps. Converting a saved plot object does not refit models or calculate
#' intervals. G-study component conversions retain raw negative points,
#' component intervals and square outlines for unavailable intervals.
#' For a G-study or D-study input, `intervals = TRUE` is forwarded to its
#' plot method and computes intervals before conversion; also supply
#' `assumption = "normal"`.
#' For multivariate D-study plots, an explicit `component = "series"` keeps
#' this dedicated conversion; use [plot_data()] for other tables.
#' Difference-interval plots from [mfrm_multivariate_d_compare()] use their
#' base `plot()` method or [plot_data()]; automatic conversion is not supported.
#' Severity profiles preserve their saved normal bounds,
#' ordered levels, guide bands, fixed points, unavailable states and precision
#' notes with the default or `component = "data"`. The word "profile" does not
#' imply a likelihood-profile interval. Conversion does not recalculate bounds
#' or qualify exploratory JML or fallback/regularized uncertainty for inference.
#' Residual-PCA scree, reference, excess and loadings plots have dedicated
#' conversions with the default or `component = "data"`. They preserve saved
#' eigenvalues, reference cutoffs, signed loadings and display encodings; method
#' identifiers are not inferred as plotted variables. No PCA is recomputed.
#' External-feature PCA scree, scores and loadings views have dedicated
#' conversions. They preserve selected axes, retained-component symbols,
#' saved group colours and shapes, labels and transformation metadata.
#' Scores retain equal axis units; loading coefficients are not correlations.
#' The default and `component = "table"` use the same complete PCA view.
#' Use `ggplot2::labs(title = NULL, subtitle = NULL)` to remove headings from
#' the returned ggplot. No PCA or clustering is refitted.
#' External-feature dendrogram conversion retains the stored merges, heights,
#' leaf order and group boxes, including tied-height cuts. The default and
#' `component = "tree"` use the same complete view. Labels follow the saved
#' setting; no tree is refitted or partition selected. Source metadata remain
#' available with [plot_data()]. Heights are not significance or branch support.
#' Imputation co-membership heatmaps have dedicated conversion with the default
#' or `component = "matrix"`. Values retain their all-imputation denominator,
#' selected ID order and fixed zero-to-one scale. Unavailable cells use grey
#' fill and crosses; zero is a valid fraction. The saved legend, labels and
#' imputation count are reused, without recomputing or pooling partitions.
#' Silhouette conversion preserves widths on a fixed minus-one-to-one scale,
#' their saved order and overall-mean line. Numeric profiles preserve
#' original-unit means, medians and counts, using separate symbols with slight
#' vertical offsets. Categorical profiles preserve the category order, unused
#' levels, group counts and within-group proportions on a zero-to-one scale.
#' The default and `component = "table"` keep the full selected view;
#' categorical profiles also accept `component = "matrix"`. These summaries
#' do not estimate uncertainty, rater quality or sampling stability.
#' For main-effects `mfrm_d_study` results, use the base `plot()` method or
#' [plot_data()] for custom graphics; automatic conversion is refused because
#' generic column selection does not preserve those design comparisons.
#' Fit-measure `measure_ci` plots retain each facet/level row and its saved
#' normal-band endpoints, including exploratory JML bands. Open diamonds mark
#' fixed values and crosses mark estimates without intervals. Conversion does
#' not recompute bounds or change their inferential qualification. The default
#' and `component = "table"` preserve the complete view and its returned notes.
#' Pooled fixed-facet MI intervals have dedicated conversion with the default
#' or `component = "table"`. Saved t-interval endpoints, degrees of freedom,
#' contrasts and complete-data information cautions are retained, with no
#' repooling or interval recalculation. Open diamonds mark fixed targets without
#' intervals, crosses mark missing intervals, and arrows mark infinite bounds.
#' Arrow tips are plotting limits, not replacement finite confidence limits.
#' Selecting `component` does not enable unsupported conversions: it cannot
#' preserve D-study differences through generic column selection. Use [plot_data()] to extract the table
#' and specify the axes, intervals and grouping explicitly in custom graphics.
#' Residual Q-Q plots require the base renderer or explicitly specified custom
#' quantile axes; automatic conversion is unavailable, even with `component`.
#' Residual-correlation heatmaps retain the selected facet on both
#' axes, a fixed -1 to 1 colour scale, and grey unavailable cells; their pair
#' table and notes retain overlap counts and interpretation qualifications.
#' Saved rating-scale counts and thresholds, and category-structure counts,
#' retain the selected values, missing outcomes, expected-count overlay and
#' saved palette. Use the default or the matching `category_table` /
#' `threshold_table` component. Threshold figures show point estimates; stored
#' SEs and CI-use flags remain in the payload, not as curve or gap intervals.
#' Category-curve bundles (including their overview), category-structure
#' boundary views and threshold ladders require their base renderer or custom
#' graphics from [plot_data()]; automatic conversion is unavailable even with
#' `component`. Generic column selection would mix distinct curve targets or
#' lose the selected view. Fitted CCC plots from
#' [plot(fit, type = "ccc")][plot.mfrm_fit] retain their dedicated conversion.
#'
#' Bias scatter/ranked views preserve their axes, saved review rules and
#' conditional interval bounds, including dashed search-limited bounds.
#' Bias heatmaps, histograms and facet profiles require base/custom graphics;
#' the generic converter cannot combine their distinct layouts safely.
#' Halo distributions, selected per-rater summaries and signed retained-edge
#' matrices preserve the selected metric, missing comparisons and caveats.
#' Halo network conversion is unavailable. Save the complete analysis result
#' and recreate its plot to obtain current interpretation notes.
#'
#' Shrinkage plots from [plot(fit, type = "shrinkage")][plot.mfrm_fit]
#' preserve both estimates for each facet level, their adjustment direction
#' and saved descriptive bands. Use the default or `component = "data"`.
#' Conversion does not turn plug-in bands into calibrated confidence intervals.
#' Recreate older plots from a fit with a refreshed shrinkage report first.
#'
#' Other draw-free payloads use a conservative tabular fallback; inspect
#' [plot_data_components()] when automatic inference is not
#' appropriate.
#' Titles, subtitles and captions are omitted by default and retained as
#' metadata. Use `plot_data(p, "notes")` to retrieve descriptions and cautions.
#' Explicitly enabled text is wrapped at 72 columns; use `ggplot2::labs()`
#' for custom labels or narrower exports.
#' Fit plots created with `show_title = FALSE` or `show_notes = FALSE` retain
#' those settings on conversion. To set them when converting a fit directly,
#' pass the flags through `...`. Interpretation and display notes remain
#' available in `attr(plot, "mfrmr_notes")`, when present in the source payload.
#' Wright, pathway, and CCC conversions share the fit-family series palette,
#' honour monochrome presets and supplied palette overrides, and preserve
#' line-type or point-shape distinctions. CCC categories retain their source
#' order. When `slope_aes = "colour"`, slopes instead use a continuous viridis
#' scale (a grey gradient in monochrome); categories still have line types.
#'
#' Fixed-facet RSM/PCM intervals retain the selected covariance method, level,
#' status symbols and optional ordinary comparison. Titles, subtitles, captions,
#' reference lines and legends can be changed through the source plot arguments.
#' Methods use both line types and vertical offsets, including in monochrome.
#'
#' @param x An `mfrm_plot_data` object, or an mfrmr object with a draw-free
#'   plot method.
#'   Corrected-JML slope, location and step plots retain their point or
#'   cumulative-distribution view, labels, reference and interpretation note;
#'   conversion does not add confidence intervals.
#' @param type Optional plot type passed to `plot()` for a non-plot-data input.
#'   A saved plot-data input already selects a view and rejects `type`.
#' @seealso [mfrmr_output_guide()] with `scope = "plots"` for selected
#'   purpose-based routes, conversion status and alternatives.
#' @param component Optional tabular payload component to convert.
#' @param ... Arguments passed to the draw-free plot method. CCC conversion
#'   additionally accepts `slope_aes`, `facet_by`, and `show_overlay`.
#'   A saved `mfrm_plot_data` object already contains its view and display
#'   settings: it rejects extra arguments (except those three CCC controls
#'   for a complete CCC view). To change the level or source plot settings,
#'   create a new payload from the fitted model or statistical result. For
#'   appearance changes after conversion, use ggplot tools such as
#'   `ggplot2::labs(title = NULL)`.
#'
#' @section Interactive views:
#' The [visual gallery](https://ryuya-dot-com.github.io/mfrmr/articles/visual-gallery.html)
#' shows how to pass a supported ggplot to `plotly::ggplotly()` for browser-based
#' hovering, zooming and legend filtering. Plotly is an optional downstream
#' renderer. Keep the original mfrmr plot and its [plot_data()] notes alongside
#' the widget; conversion does not recalculate estimates or intervals. Check
#' labels, custom geoms and interval layers after conversion. The gallery includes
#' static figures, numeric tables and an HTML export recipe.
#'
#' @return A `ggplot2` plot object, with a `mfrmr_notes` attribute when the
#'   source plot payload contains a `notes` table.
#'   Paired Wright/CCC payloads from [plot(reference, comparison = comparison)][plot.mfrm_fit] retain their selected
#'   comparison or difference view, group selection and monochrome panel policy.
#' @examples
#' # Long-running workflow; run explicitly when needed.
#' # example("as_ggplot", run.dontrun = TRUE)
#' \dontrun{
#' toy <- load_mfrmr_data("example_core")
#' # A balanced slice retains every Rater and Criterion while running quickly.
#' toy <- toy[toy$Person %in% unique(toy$Person)[1:12], , drop = FALSE]
#' fit <- fit_mfrm(toy, "Person",
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
#' @name as_ggplot.mfrm_signal_detection
#' @aliases as_ggplot.mfrm_signal_detection
#' @usage \method{as_ggplot}{mfrm_signal_detection}(x, type = NULL, component = NULL, ...)
#' @rawNamespace S3method(as_ggplot,mfrm_signal_detection,as_ggplot.mfrm_design_evaluation)
NULL

#' @rdname as_ggplot
#' @name as_ggplot.mfrm_signal_detection_plot_data
#' @aliases as_ggplot.mfrm_signal_detection_plot_data
#' @usage \method{as_ggplot}{mfrm_signal_detection_plot_data}(x, type = NULL, component = NULL, ...)
#' @rawNamespace S3method(as_ggplot,mfrm_signal_detection_plot_data,as_ggplot.mfrm_design_evaluation_plot_data)
NULL

#' @rdname as_ggplot
#' @export
as_ggplot.mfrm_plot_data <- function(x, type = NULL, component = NULL, ...) {
  if (identical(x$name, "network_ergm") ||
      (identical(x$name, "rater_network") && any(x$data$metric %in% c("SWS", "PUS")))) {
    stop("Automatic network-fit conversion is unavailable; use plot() or plot_data() to retain the selected metric, intervals and unavailable states in custom graphics.", call. = FALSE)
  }
  if (x$name %in% c("subset_connectivity", "rater_network", "halo_network") &&
      identical(x$data$plot, "network")) {
    stop("Automatic network conversion is unavailable; use plot() with layout/graphical controls, or plot_data() for node, edge and coordinate tables.", call. = FALSE)
  }
  if (identical(x$name, "residual_network")) {
    stop("Use plot() on the saved residual network for network/dendrogram views; as_ggplot() supports its type = 'matrix' view only.", call. = FALSE)
  }
  if (!is.null(type)) stop("A saved plot already selects a view. Recreate it from the fit or result to select another view.", call. = FALSE)
  if (x$name %in% c("anchor_drift", "anchor_drift_forest", "equating_chain_common_anchors")) {
    rlang::check_dots_empty()
    selected <- if (x$name != "anchor_drift") "data" else
      if (identical(x$data$plot, "heatmap")) "matrix" else "table"
    if (!is.null(component) && !identical(component, selected)) stop(
      "Use the selected linking view, or plot_data() for other components.", call. = FALSE)
    .require_mfrmr_ggplot2()
    return(.mfrmr_gg_linking(x))
  }
  if (identical(x$name, "shrinkage")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "data")) stop(
      "Use component = 'data' for the shrinkage comparison, or plot_data() for other components.", call. = FALSE)
    .require_mfrmr_ggplot2()
    return(.mfrmr_gg_shrinkage(x))
  }
  if (x$name %in% c("table13_bias", "halo_network")) {
    rlang::check_dots_empty()
    allowed <- if (x$name == "table13_bias") c(scatter = "scatter_data", ranked = "ranked_table") else
      c(edge_distribution = "table", halo_summary = "table", matrix = "matrix")
    if (!x$data$plot %in% names(allowed)) stop(
      "Automatic ggplot conversion of this rater-effect view is unavailable. Use plot() on the saved result, or plot_data() for custom graphics.", call. = FALSE)
    if (!is.null(component) && !identical(component, unname(allowed[x$data$plot]))) stop(
      "Use the selected rater-effect view, or plot_data() for other components.", call. = FALSE)
    .require_mfrmr_ggplot2()
    return(if (x$name == "table13_bias") .mfrmr_gg_bias(x) else .mfrmr_gg_halo(x))
  }
  if (x$name %in% c("rating_scale", "category_structure", "category_curves", "threshold_ladder")) {
    rlang::check_dots_empty()
    supported <- (x$name %in% c("rating_scale", "category_structure") && identical(x$data$plot, "counts")) ||
      (x$name == "rating_scale" && identical(x$data$plot, "thresholds"))
    if (!supported) stop(
      "Automatic ggplot conversion of this category view is unavailable. Use plot() on the original fit/result for base graphics, or plot_data() for explicitly grouped custom graphics.", call. = FALSE)
    expected <- if (x$data$plot == "counts") "category_table" else "threshold_table"
    if (!is.null(component) && !identical(component, expected)) stop(
      "Use the selected category view, or plot_data() for other components.", call. = FALSE)
    .require_mfrmr_ggplot2()
    return(.mfrmr_gg_categories(x))
  }
  if (x$name %in% c("interrater", "rater_agreement_heatmap")) {
    rlang::check_dots_empty()
    .require_mfrmr_ggplot2()
    expected <- if (x$name == "interrater") "table" else "matrix"
    if (!is.null(component) && !identical(component, expected))
      stop("Use the selected agreement view, or plot_data() for other components.", call. = FALSE)
    return(if (x$name == "interrater") .mfrmr_gg_interrater(x) else .mfrmr_gg_agreement_heatmap(x))
  }
  if (identical(x$name, "rater_severity_profile")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "data")) stop(
      "Use component = 'data' for the severity profile, or plot_data() for other components.", call. = FALSE)
    .require_mfrmr_ggplot2()
    return(.mfrmr_gg_severity(x))
  }
  if (identical(x$name, "residual_qq")) {
    stop("Automatic ggplot conversion of residual Q-Q plots is unavailable. Use plot(fit, type = \"residual_qq\") or plot_data() to specify both quantile axes explicitly.", call. = FALSE)
  }
  if (identical(x$name, "residual_pca")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "data")) {
      stop("Use component = 'data' for the selected residual-PCA view, or plot_data() for other components.", call. = FALSE)
    }
    return(.mfrmr_gg_residual_pca(x))
  }
  if (identical(x$name, "gpcm_slope_profile")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table"))
      stop("Use component = 'table' for the complete profile view, or plot_data() for other tables.", call. = FALSE)
    return(mfrm_gg_gpcm_profile(x))
  }
  if (identical(x$name, "population_intervals")) {
    rlang::check_dots_empty()
    if (!is.null(component)) stop("Use plot_data() to select population tables; as_ggplot() preserves the selected interval view.", call. = FALSE)
    return(.mfrm_gg_population_intervals(x))
  }
  if (identical(x$name, "facet_interval_methods")) {
    if (!is.null(component)) stop("Use plot_data() to select a fixed-facet table; as_ggplot() preserves the complete interval view.", call. = FALSE)
    rlang::check_dots_empty()
    return(mfrm_gg_facet_intervals(x))
  }
  if (identical(x$name, "posterior_intervals")) {
    rlang::check_dots_empty()
    if (!is.null(component)) stop("Use plot_data() to select posterior tables; as_ggplot() preserves the selected interval view.", call. = FALSE)
    return(.mfrm_gg_posterior_intervals(x))
  }
  if (identical(x$name, "dstudy_intervals")) {
    rlang::check_dots_empty()
    if (!is.null(component)) stop("Use plot_data() to select tables; as_ggplot() preserves all interval rows.", call. = FALSE)
    return(.mfrm_gg_dstudy_intervals(x))
  }
  if (x$name %in% c("feature_pca_scree", "feature_pca_scores", "feature_pca_loadings")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) {
      stop("Use component = 'table' for the complete PCA view, or plot_data() for other components.", call. = FALSE)
    }
    return(.mfrmr_gg_feature_pca(x))
  }
  if (identical(x$name, "cluster_dendrogram")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "tree")) {
      stop("Use component = 'tree' for the complete dendrogram, or plot_data() for other components.", call. = FALSE)
    }
    return(.mfrmr_gg_dendrogram(x))
  }
  if (identical(x$name, "cluster_co_membership")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "matrix")) {
      stop("Use component = 'matrix' for the complete co-membership view, or plot_data() for other components.", call. = FALSE)
    }
    return(.mfrmr_gg_co_membership(x))
  }
  if (identical(x$name, "pooled_facet_intervals")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) {
      stop("Use component = 'table' for the complete pooled interval view, or plot_data() for other components.", call. = FALSE)
    }
    return(.mfrmr_gg_pooled(x))
  }
  if (identical(x$name, "screening_performance")) {
    stop("Use plot() for screening performance, or plot_data() for custom graphics; automatic ggplot conversion is not available.", call. = FALSE)
  }
  if (identical(x$name, "multivariate_d_comparison")) {
    stop("Use plot() for D-study difference intervals, or plot_data() for custom graphics; automatic ggplot conversion is not available.", call. = FALSE)
  }
  if (x$name %in% c("cluster_silhouette", "cluster_profile")) {
    rlang::check_dots_empty()
    allowed <- if (identical(x$name, "cluster_profile") && !is.null(x$data$matrix)) c("table", "matrix") else "table"
    if (!is.null(component) && !component %in% allowed) {
      stop("Use the default component for the complete cluster view, or plot_data() for other components.", call. = FALSE)
    }
    return(.mfrmr_gg_cluster_summary(x))
  }
  if (identical(x$name, "d_study")) {
    stop("Automatic ggplot conversion is not available for mfrm_d_study plots. ",
      "Use plot(x) to preserve the selected design comparisons, or plot_data() for custom graphics.", call. = FALSE)
  }
  .require_mfrmr_ggplot2()
  payload <- x$data %||% list()
  dots <- list(...)
  if (identical(x$name, "design_recommendation")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table"))
      stop("Use component = 'table' for criterion outcomes; use plot_data() for the other tables.", call. = FALSE)
    tab <- payload$table
    tab$Design <- factor(tab$DesignRow, levels = payload$designs$DesignRow,
      labels = paste0(payload$designs$DesignRow, ifelse(payload$designs$Selected, " *", "")))
    tab$Check <- factor(tab$Check, levels = rev(unique(tab$Check)))
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = .data$Design, y = .data$Check, fill = .data$Status)) +
      ggplot2::geom_tile(colour = "#666666") +
      ggplot2::scale_fill_manual(values = c(pass = "#2A9D8F", fail = "#D55E00",
        unavailable = "#999999", not_requested = "#FFFFFF"), drop = FALSE) +
      .mfrmr_gg_theme() +
      ggplot2::labs(x = "Ranked design row (* selected)", y = NULL, fill = "Criterion")
    attr(p, "mfrmr_plot_data") <- x
    attr(p, "mfrmr_notes") <- payload$notes
    return(p)
  }
  if (identical(x$name, "bias_collection")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) stop(
      "Use component = 'table' for the complete overview, or plot_data() for other components.", call. = FALSE)
    return(.mfrmr_gg_bias_collection(payload))
  }
  if (identical(x$name, "fit_measures") && identical(payload$plot, "measure_ci")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) stop(
      "Use component = 'table' for measure bands, or plot_data() for other components.", call. = FALSE)
    return(.mfrmr_gg_fit_measure_bands(x))
  }
  if (identical(x$name,"extended_model_map")) {
    rlang::check_dots_empty()
    if(!is.null(component) && !identical(component,"table")) stop("Use plot_data() for other map components.",call.=FALSE)
    return(.mfrmr_gg_model_map(payload))
  }
  if (identical(x$name, "response_diagnostics")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) stop("Use component = 'table', or plot_data() for other diagnostic components.", call. = FALSE)
    return(.mfrmr_gg_response_diagnostics(payload))
  }
  if (identical(x$name, "screening_sensitivity")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) stop("Use component = 'table', or plot_data() for other screening components.", call. = FALSE)
    return(.mfrmr_gg_screening_sensitivity(payload))
  }
  if (identical(x$name, "extended_model_comparison")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) stop("Use component = 'table', or plot_data() for other comparison components.", call. = FALSE)
    return(.mfrmr_gg_extended_comparison(payload))
  }
  if (x$name %in% c("testlet_calibration", "testlet_scores", "random_rater_scores", "person_scores", "random_rater_severity", "random_rater_interval_bootstrap", "jml_adjustment")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) stop(
      "Extended-model graphics use component = 'table'; use plot_data() for other components.", call. = FALSE)
    return(.mfrmr_gg_extended_payload(payload, x$name))
  }
  if (identical(x$name, "multivariate_gstudy")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "table")) stop(
      "G-study conversion uses component = 'table'; use plot_data() for other components.", call. = FALSE)
    return(.mfrmr_gg_multivariate_gstudy(payload))
  }
  if (identical(x$name, "multivariate_d_study")) {
    rlang::check_dots_empty()
    if (!is.null(component) && !identical(component, "series")) stop(
      "D-study conversion uses component = 'series'; use plot_data() for tables and other components.", call. = FALSE)
    return(.mfrmr_gg_multivariate_d_study(payload))
  }
  conversion_args <- if (is.null(component) && x$name %in% c(
    "category_characteristic_curves", "category_characteristic_curves_overlay"
  )) c("slope_aes", "facet_by", "show_overlay") else character()
  if (length(dots) && (is.null(names(dots)) || anyDuplicated(names(dots)) ||
      any(!names(dots) %in% conversion_args))) {
    stop("Unsupported arguments for a saved plot. Recreate the plot from its fit or result to change settings; use ggplot2::labs() for converted labels. Only complete CCC views accept slope_aes, facet_by and show_overlay.", call. = FALSE)
  }
  if (is.null(component) && identical(x$name, "fair_average") && !is.null(payload$plot_data)) {
    df <- payload$plot_data
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$X, y = .data$Y))
    ci <- df[is.finite(df$Lower) & is.finite(df$Upper), , drop = FALSE]
    if (payload$plot == "difference") {
      p <- p + ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$X, yend = .data$Y), colour = "grey65") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::scale_y_continuous(breaks = df$Y, labels = df$Label)
    } else if (payload$plot == "scatter") {
      limits <- range(df$X, df$Y, ci$Lower, ci$Upper, finite = TRUE)
      p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::coord_equal(xlim = limits, ylim = limits)
    }
    if (nrow(ci)) {
      p <- p + if (payload$plot == "measure") {
        ggplot2::geom_segment(data = ci, ggplot2::aes(y = .data$Lower, yend = .data$Upper, xend = .data$X, colour = .data$Facet))
      } else ggplot2::geom_segment(data = ci, ggplot2::aes(x = .data$Lower, xend = .data$Upper, yend = .data$Y, colour = .data$Facet))
    }
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$Facet, shape = .data$Facet)) +
      ggplot2::scale_colour_manual(values = stats::setNames(payload$encoding$Color, payload$encoding$Facet)) +
      ggplot2::scale_shape_manual(values = stats::setNames(payload$encoding$Shape, payload$encoding$Facet)) +
      .mfrmr_gg_theme()
    if (nrow(payload$encoding) == 1L || payload$plot == "difference") p <- p + ggplot2::theme(legend.position = "none")
    # Long reference/uncertainty explanations remain available outside the figure.
    payload$subtitle <- if (!is.null(payload$ci_note)) {
      if (identical(payload$plot, "difference")) "Approximate fair-score intervals; observed means held fixed."
      else "Approximate fair-score intervals omit some estimation uncertainty."
    } else NULL
    if (identical(payload$interpretation_status, "review_only")) {
      payload$subtitle <- paste(c("REVIEW ONLY", payload$subtitle), collapse = " - ")
    }
    return(.mfrmr_gg_labs(p, payload, x = payload$xlab, y = payload$ylab))
  }
  if (is.null(component) && identical(x$name, "paired_model_comparison")) {
    return(.mfrmr_gg_comparison(payload))
  }
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
  calibration_score_type <- sub(
    "^calibration_score_", "", as.character(x$name %||% "")
  )
  if (calibration_score_type %in% c("interval", "precision", "edge_mass")) {
    if (!is.null(component) && !identical(component, "data")) stop(
      "Use component = 'data' for the selected portable-score view, or plot_data() for other components.", call. = FALSE)
    return(.mfrmr_gg_calibration_score(payload, calibration_score_type))
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
